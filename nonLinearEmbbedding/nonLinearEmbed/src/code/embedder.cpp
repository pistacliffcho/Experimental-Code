

class Node{
public:
  std::vector<double> coord;
  double eta;
};

class Embedder{
public:
  int k;
  // Rank of embedding
  int n;
  // Node of nodes
  double h;
  // Step size for numeric derivatives
  std::vector<Node> nodes;
  // Vector of all node info. Contains coordinates and eta
  double calcDiv(int i, int j);
  // Base function for calculating divergence/distance
  double calcEdgeProb(int i, int j);
  // Calculates edge probability 
  double calcLLKcont(int i, int j, bool hasEdge);
  // Calculates the log-likelihood contribution for i, j
  void calcDervs(int i, int j, bool hasEdge);
  // Calculates derivatives with respect to node i ONLY

  void calcDervs(int i, int j, bool hasEdge, double w);
  // Calculates derivatives with respect to node i ONLY
  
  void sgdUpdate(int i, int j, double alpha, double w, bool hasEdge);
  // Does a single SGD update *only on node i*
  
  void adamUpdate(int i, int j, double alpha, double w, bool hasEdge);
  // Does a single adam update *only on node i*
  
  std::vector<double> derv;
  // Vector for storing derivatives 
  
  std::vector<Adam*> adamList;
  
  Embedder(Rcpp::NumericMatrix coords, 
           Rcpp::NumericVector eta, 
           double h);
  ~Embedder();
};

// For now, we will just look at Euclidean Distance 
double Embedder::calcDiv(int i, int j){
  double ans = 0;
  Node* n_i = &nodes[i];
  Node* n_j = &nodes[j];
  double xi, xj;
  for(int ii = 0; ii < k; ii++){
    xi = n_i->coord[ii];
    xj = n_j->coord[ii];
    ans += std::abs(xi-xj);
  }
  return(ans);
}

double expit(double x){
  double ans = exp(x) / (1.0 + exp(x));
  return(ans);
}
double Embedder::calcEdgeProb(int i, int j){
  double div = calcDiv(i,j);
  double eta_sum = nodes[i].eta + nodes[j].eta;
  double ans = exp(-div - exp(-eta_sum));
  return(ans);
}

double Embedder::calcLLKcont(int i, int j, bool hasEdge){
  double prob = calcEdgeProb(i,j);
  double ans;
  if(hasEdge){ ans = log(prob); }
  else { ans = log(1.0-prob);}  
  return(ans);
}

// For now, just calculating numeric derivatives. Not fast!
void Embedder::calcDervs(int i, int j, bool hasEdge){
  derv.resize(k + 1);
  Node* n_i = &nodes[i];
  double llk_h, llk_l;
  for(int ii = 0; ii < k; ii++){
    n_i->coord[ii] += h;
    llk_h = calcLLKcont(i,j, hasEdge);
    n_i->coord[ii] -= 2.0 * h;
    llk_l = calcLLKcont(i,j, hasEdge);
    derv[ii] = (llk_h - llk_l) / (2.0 * h);
    n_i->coord[ii] += h;
  }  
  n_i->eta += h;
  llk_h = calcLLKcont(i,j, hasEdge);
  n_i->eta -= 2.0 * h;
  llk_l = calcLLKcont(i,j, hasEdge);
  n_i->eta += h;
  derv[k] = (llk_h - llk_l) / (2.0 * h);
}

void Embedder::calcDervs(int i, int j, bool hasEdge, double w){
	calcDervs(i,j, hasEdge);
	for(int ii = 0; ii <(k+1); ii++) derv[ii] *= w;
}

// Returns the min(absolute(x, abs_rng)) * sign
double ab_thresh(double x, double abs_rng){
  if(x < 0){
    if(x > abs_rng) return(x);
    return(-abs_rng);
  }
  if(x < abs_rng) return(x);
  return(abs_rng);
}

double sign(double x){
  if(x > 0) return(1.0);
  return(-1.0);
}

void Embedder::sgdUpdate(int i, int j, double alpha, 
                         double w, bool hasEdge){
  calcDervs(i, j, hasEdge);
  Node* ni = &nodes[i];
  double update;
  for(int ii = 0; ii < k; ii++){
    update = alpha * w * derv[ii];
    ni->coord[ii] += update;
  }
  update = alpha * w * derv[k];
  ni->eta += update;
}

void Embedder::adamUpdate(int i, int j, double alpha, 
                         double w, bool hasEdge){
  calcDervs(i, j, hasEdge, w);
  Node* ni = &nodes[i];
  Adam* this_adam = adamList[i];
  this_adam->calc_update(alpha, &derv[0]);
  std::vector<double> delta = this_adam->getDelta();
  for(int ii = 0; ii < k; ii++){ ni->coord[ii] -= delta[ii]; }
  // Subtracting delta because we are using LLK, not negative LLK
  ni->eta -= delta[k];
}



Embedder::Embedder(Rcpp::NumericMatrix R_coords, 
                   Rcpp::NumericVector R_eta, 
                   double R_h){
  n = R_eta.size();
  k = R_coords.size()/n;
  int numAdamPars = k + 1;
  nodes.resize(n);
  adamList.resize(n);
  for(int i = 0; i < n; i++){
    nodes[i].eta = R_eta[i];
    nodes[i].coord.resize(k);
    for(int j = 0; j < k; j++){
      nodes[i].coord[j] = R_coords(i,j);
    }
    adamList[i] = new Adam(numAdamPars);
  }
  h = R_h;
}

Embedder::~Embedder(){
	int nAdams = adamList.size();
	for(int i = 0; i < nAdams; i++){
		delete adamList[i];
	}
}

// [[Rcpp::export]]
double estLLK(Rcpp::NumericVector is, 
              Rcpp::NumericVector js,
              Rcpp::LogicalVector hasEdges, 
              Rcpp::NumericVector ws, 
              Rcpp::NumericMatrix coord,
              Rcpp::NumericVector etas){
  double h = 0;
  unsigned int n_updates = is.size();
  if(n_updates != js.size()){ Rcpp::stop("n_updates != js.size"); }
  if(n_updates != hasEdges.size()){ Rcpp::stop("n_updates != hasEdges.size"); }
  if(n_updates != ws.size()){ Rcpp::stop("n_updates != ws.size");}

  Embedder updater(coord, etas, h);
  double ans = 0;
  for(int i = 0; i < n_updates; i++){
    bool c_hasEdge = hasEdges[i] == TRUE;
    ans += ws[i] * updater.calcLLKcont(is[i] - 1, js[i] - 1, c_hasEdge);
  }
  ans = ans / n_updates;
  return(ans);
}

// [[Rcpp::export]]
List sgd_updates(Rcpp::NumericVector is, 
                 Rcpp::NumericVector js,
                 Rcpp::LogicalVector hasEdges,
                 Rcpp::NumericVector ws, 
                 Rcpp::NumericVector alphas,
                 Rcpp::NumericMatrix coord, 
                 Rcpp::NumericVector etas, 
            double h){
  unsigned int n_updates = is.size();
  if(n_updates != js.size()){ Rcpp::stop("n_updates != js.size"); }
  if(n_updates != hasEdges.size()){ Rcpp::stop("n_updates != hasEdges.size"); }
  if(n_updates != ws.size()){ Rcpp::stop("n_updates != ws.size");}
  if(n_updates != alphas.size()){Rcpp::stop("n_updates != alpha.size");}
  
  Embedder updater(coord, etas, h);
  for(int i = 0; i < n_updates; i++){
    bool c_hasEdge = hasEdges[i] == TRUE;
    updater.sgdUpdate(is[i] - 1, js[i] - 1, alphas[i], 
                      ws[i], c_hasEdge);
    updater.sgdUpdate(js[i] - 1, is[i] - 1, alphas[i], 
                      ws[i], c_hasEdge);
    
  }
  int n = etas.size();
  int k = coord.size()/n;
  
  Rcpp::NumericMatrix coord_out(n,k);
  Rcpp::NumericVector etas_out(n);
  for(int i = 0; i < n; i++){
    etas_out[i] = updater.nodes[i].eta;
    for(int j = 0; j < k; j++){
      coord_out(i,j) = updater.nodes[i].coord[j];
    }
  }
  Rcpp::List ans = Rcpp::List::create(coord_out, etas_out);
  return(ans);
}


// [[Rcpp::export]]
List adam_updates(Rcpp::NumericVector is, 
                 Rcpp::NumericVector js,
                 Rcpp::LogicalVector hasEdges,
                 Rcpp::NumericVector ws, 
                 Rcpp::NumericVector alphas,
                 Rcpp::NumericMatrix coord, 
                 Rcpp::NumericVector etas, 
            double h){
  unsigned int n_updates = is.size();
  if(n_updates != js.size()){ Rcpp::stop("n_updates != js.size"); }
  if(n_updates != hasEdges.size()){ Rcpp::stop("n_updates != hasEdges.size"); }
  if(n_updates != ws.size()){ Rcpp::stop("n_updates != ws.size");}
  if(n_updates != alphas.size()){Rcpp::stop("n_updates != alpha.size");}
  
  Embedder updater(coord, etas, h);
  for(int i = 0; i < n_updates; i++){
    bool c_hasEdge = hasEdges[i] == TRUE;
    updater.adamUpdate(is[i] - 1, js[i] - 1, alphas[i], 
                      ws[i], c_hasEdge);
    updater.adamUpdate(js[i] - 1, is[i] - 1, alphas[i], 
                      ws[i], c_hasEdge);
    
  }
  int n = etas.size();
  int k = coord.size()/n;
  
  Rcpp::NumericMatrix coord_out(n,k);
  Rcpp::NumericVector etas_out(n);
  for(int i = 0; i < n; i++){
    etas_out[i] = updater.nodes[i].eta;
    for(int j = 0; j < k; j++){
      coord_out(i,j) = updater.nodes[i].coord[j];
    }
  }
  Rcpp::List ans = Rcpp::List::create(coord_out, etas_out);
  return(ans);
}
