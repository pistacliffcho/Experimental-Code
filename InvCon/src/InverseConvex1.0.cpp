#include <Rcpp.h>
#include <math.h>
#include <stdio.h>
using namespace Rcpp;



// [[Rcpp::export]]

double ic_llk(NumericVector x, NumericVector betas, NumericVector weight) /* Likelihood Function*/
	{
	int k = x.size();
	double n = 0;
	double main_comp = 0;
	for(int i = 0; i < k; i++)
		{
		main_comp = main_comp - log(betas[i]) * weight[i];
		n = n + weight[i];
		}
	double numer = 0;
	double db = 0;
	double dx = 0;
	double mass = 0;
	for(int i = 0; i < (k-1); i++)
		{
		db = betas[i+1] - betas[i];
		dx = x[i+1] - x[i];
		numer = log(betas[i+1]) - log(betas[i]);
		if(db == 0)
			mass = mass + 1/betas[i] * dx;
		if(db != 0)
			mass = mass + numer * dx/db;
		}
	double lk = main_comp - n * log(mass);
	return(lk);
	}

// [[Rcpp::export]]

double ic_mass(NumericVector x, NumericVector betas) /*Computes Regularizing Constant*/
	{
	int k = x.size();
	double numer = 0;
	double db = 0;
	double dx = 0;
	double mass = 0;
	for(int i = 0; i < (k-1); i++)
		{
		db = betas[i+1] - betas[i];
		dx = x[i+1] - x[i];
		numer = log(betas[i+1]) - log(betas[i]);
		if(db == 0)
			mass = mass + 1/betas[i] * dx;
		if(db != 0)
			mass = mass + numer * dx/db;
		}
	return(mass);
	}
	
// [[Rcpp::export]]	

void move_active(int ind, double h, IntegerVector act_inds, NumericVector x, NumericVector betas) /* Moves Heights in ActiveSet Parameterization */
	{
	int a_k = act_inds.size();
	int k = betas.size();
	betas[ind-1] = betas[ind-1] + h;
	int ind_l = -1;
	int ind_r = -1;
	double slp;
	
	if(ind-1 < 0)
		{
		printf("Warning: ind-1 < 0 in move_active!\n");
		return;
		}
	
	if(ind > k)
		{
		printf("Warning: ind > k in move_active!\n");
		return;
		}
	
	if(a_k == 2 & ind == act_inds[0])
		{
		slp = (betas[act_inds[1]-1] - betas[act_inds[0] -1 ]) / (x[act_inds[1]-1] -x[act_inds[0]-1]);
		for(int i = act_inds[1]; i >= act_inds[0]; i--)
			betas[i-1] = betas[act_inds[1] - 1] + slp * (x[i-1] - x[act_inds[1]-1]);
		return;
		}
	
	if(a_k == 2 & ind == act_inds[1])
		{
		slp = (betas[act_inds[1]-1] - betas[act_inds[0] -1 ]) / (x[act_inds[1]-1] -x[act_inds[0]-1]);
		for(int i = act_inds[0]; i <= act_inds[1]; i ++)		
			betas[i-1] = betas[act_inds[0] - 1] + slp * (x[i-1] - x[act_inds[0]-1]);
		return;
		}
	
	
	if(ind > act_inds[0])
		{
		for(int i = 1; i < a_k; i++)
			{
			if(act_inds[i]	>= ind )
				{
				ind_l = act_inds[i-1];
				break;
				}
			}
		if(ind_l < ind - 1 & ind_l > -1)
			{
			slp = (betas[ind-1] - betas[ind_l-1])/(x[ind-1] - x[ind_l-1]);
			
			for(int i = ind_l; i < ind-1; i++) 
				{
				betas[i] = betas[ind_l-1] + slp * (x[i] - x[ind_l-1]);
				}
			}
		}
	if(ind < act_inds[a_k-1])
		{
		for(int i = a_k - 2; i >= 0; i--)
			{
			if(act_inds[i] <= ind)
				{
				ind_r = act_inds[i+1];
				break;
				}
			}
	 	if(ind_r > ind + 1 & ind_r > -1)
	 		{
	 		slp = (betas[ind_r-1] - betas[ind-1])/(x[ind_r-1] - x[ind-1]);
	 		for(int i = ind; i < ind_r-1; i++)
	 			{
	 			betas[i] = betas[ind-1] + slp * (x[i] - x[ind-1]);
	 			}
	 		}
		}	
	}
	
double min(double a, double b)
	{
	if(a <= b)
		return(a);
	return(b);
	}	
	
// [[Rcpp::export]]		
	
void ic_act_ders(int ind, NumericVector d_vec, IntegerVector act_inds, NumericVector x, NumericVector betas, NumericVector weight, double h = 0.0001)  /*Computes Numerical Derivatives*/
	{
	int k = x.size();
	if(ind < 1)
		{
		printf("warning: ind < 1! ind = %d\n", ind);
		return;
		}
	if(ind > k)
		{
		printf("warning: ind > k-1! k = %d, ind = %d\n",k, ind );
		d_vec[0] = 0;
		d_vec[1] = -1;
		return;
		}
	double h_use = min(h, betas[ind-1]/2); 
	double lk_0 = ic_llk(x, betas, weight);
	move_active(ind, h_use, act_inds, x, betas);
	double lk_h = ic_llk(x, betas, weight);
	move_active(ind, -2*h_use, act_inds, x, betas);
	double lk_l = ic_llk(x, betas, weight);
	move_active(ind, h_use, act_inds, x, betas);

	d_vec[0] = (lk_h - lk_l)/(2*h_use);
	d_vec[1] = (lk_h + lk_l - 2* lk_0)/(h_use*h_use);
	}

// [[Rcpp::export]]
		
void ic_active(LogicalVector active, NumericVector x, NumericVector betas, double slack = 0.00000001) /*Finds Active Points*/
	{
	double slp_old;
	double slp_new= -INFINITY;
	int k = x.size();
	for(int i = 0; i < k-1; i++)
		{
		if(betas[i] == -INFINITY & betas[i+1] == -INFINITY)
			{
			active[i] = FALSE;
			slp_new = -INFINITY;
			continue;
			}
		slp_old = slp_new;
		slp_new = (betas[i+1] - betas[i])/(x[i+1] - x[i]);
		active[i] = slp_new  > slp_old + slack;
		}
	active[k-1] = betas[k-1] < INFINITY; 
	}	

// [[Rcpp::export]]		

void ic_aerr_vec(NumericVector e_vec, IntegerVector act_inds, NumericVector x, NumericVector betas, NumericVector weight) /*Calculates the analytic derivatives of individual points */
	{
	int k = x.size();
	double mass = ic_mass(x, betas);
	double n = 0;
	double t1 = 0;
	double t2 = 0;
	double dx = 0;
	double db = 0;
	for(int i = 0; i < k; i++)
		n = n + weight[i];
	for(int i = 1; i < k-1; i++)
		{
		dx = x[i+1] - x[i];
		db = (betas[i+1] - betas[i]);
		if(db != 0)
			t1 = dx/(db * db) * (-betas[i+1]/betas[i] + 1 + log(betas[i+1]) - log(betas[i]) );
		if(db == 0)
			t1 = -1/(betas[i]*betas[i]) * dx * 0.5;
		dx = x[i] - x[i-1];
		db = (betas[i] - betas[i-1]);
		if(db == 0)
			t2 = -1/(betas[i] * betas[i]) * dx * 0.5;
		if(db != 0)
			t2 = dx/(db*db) * (1 - betas[i-1]/betas[i] - log(betas[i]) + log(betas[i-1]) );
		e_vec[i] = -1/betas[i] * weight[i] - n/mass * (t1 + t2);
		}
		
	dx = x[1] - x[0];
	db = (betas[1] - betas[0]);
	if(db != 0)
		t1 = dx/(db * db) * (-betas[1]/betas[0] + 1 + log(betas[1]) - log(betas[0]) );
	if(db == 0)
		t1 = -1/(betas[0]*betas[0]) * dx * 0.5;
	e_vec[0] = -1/betas[0] * weight[0] - n/mass * t1;
	
	dx = x[k-1] - x[k-2];
	db = (betas[k-1] - betas[k-2]);
	if(db == 0)
		t2 =- dx * 1/(betas[k-1]*betas[k-1]) * 0.5 ;
	if(db != 0)
		t2 = dx/(db*db) * (1 - betas[k-2]/betas[k-1] - log(betas[k-1]) + log(betas[k-2]) );
	e_vec[k-1] = -1/betas[k-1] * weight[k-1] - n/mass * t2;
	}

// [[Rcpp::export]]	

void act_d_vec(NumericVector x, NumericVector d_vec, IntegerVector act_inds, NumericVector a_d_vec) /*Transforms derivatives in basic set to active set derivatives*/
{
	int k = x.size();
	int a_k = act_inds.size();
	for(int i = 0; i < a_k; i++)
		act_inds[i]--;
	int cur_l = act_inds[0];
	int cur_r = act_inds[1];
	int act_cnt = 1; 
	for(int i = 0; i < k; i++)
		a_d_vec[i] = 0;
	double l_sum = 0;
	double r_sum = 0;
	bool New_Act = TRUE;
	

	for(int i = cur_l+1; i < cur_r; i++)
		r_sum = r_sum + d_vec[i] * (x[cur_r] - x[i])/(x[cur_r] - x[cur_l]);
	a_d_vec[cur_l] = r_sum + d_vec[cur_l];
	for(int i = act_inds[0] + 1; i < act_inds[a_k-1]; i++)
		{

		l_sum = (l_sum + d_vec[i-1]) * (x[i-1] - x[cur_l])/(x[i]-x[cur_l]) ;
		r_sum = r_sum  * (x[cur_r] - x[i-1])/(x[cur_r] - x[i]) - d_vec[i];



		if(New_Act == TRUE)
			{
			New_Act = FALSE;
			cur_l = act_inds[act_cnt-1];
			l_sum = 0;
			if(cur_r == i)
				r_sum = 0;


			} 
		if(i == cur_r)
			{
			act_cnt++;
			if(act_cnt == a_k)
				printf("Warning: act_cnt == a_k in a_d_vec!\n");
			cur_r = act_inds[act_cnt];
			r_sum = 0;
			

			
			for(int j = i + 1; j < cur_r; j++)
				{
					r_sum = r_sum + d_vec[j] * (x[cur_r] - x[j])/(x[cur_r] - x[i]);
				}	
			New_Act = TRUE;
			}
			a_d_vec[i] = l_sum + r_sum + d_vec[i];	
		}
	int k_last = act_inds[a_k-1];
	l_sum = (l_sum + d_vec[k_last-1]) * (x[k_last-1] - x[cur_l])/(x[k_last]-x[cur_l]);
	a_d_vec[k_last] = l_sum + d_vec[k_last];
	
	
	for(int i = 0; i < a_k; i++)
		act_inds[i]++;

}

	
// [[Rcpp::export]]			
	
double p_ic(double x, NumericVector fit_x, NumericVector fit_dens)  /* calculates estimated cdf from fit */
	{
	double p = 0;
	double denom = 0;
	double db = 0;
	double dx = 0;
	double mass = 0;
	int it = 0;
	int k = fit_x.size();
	if(x <= fit_x[0])
		return(0);
	if(x >= fit_x[k-1])
		return(1);
	
	double raw_p = 0;
	double tot_p = 0;	
		
	for(int i = 0; i < k-1; i++)
		{
		db = 1/fit_dens[i+1] - 1/fit_dens[i];
		dx = fit_x[i+1] - fit_x[i];
		denom = log(1/fit_dens[i+1]) - log(1/fit_dens[i]);
		
		if(x < fit_x[i+1] & x >= fit_x[i])
			{
			it = i;
			raw_p = tot_p;
			}
		
		if(db == 0)
			tot_p = tot_p + 1/fit_dens[i] * dx ;
		if(db != 0)
			tot_p = tot_p + denom * dx/db;
		}
	double b_new = 1/fit_dens[it] + (x - fit_x[it])  * (1/fit_dens[it + 1] - 1/fit_dens[it]) / (fit_x[it+1] - fit_x[it]);
	db = b_new - 1/fit_dens[it];
	dx = x - fit_x[it];
	denom = log(b_new) - log(1/fit_dens[it]);
	if(db == 0)
		raw_p = raw_p + 1/fit_dens[it] * dx ;
	if(db != 0)
		raw_p = raw_p + denom * dx/db;
	return(raw_p/tot_p);
	}
	
		
// [[Rcpp::export]]

double d_ic(double x, NumericVector fit_x, NumericVector fit_dens)    /* calculates estimated density */
	{
		
	int k = fit_x.size();
	double numer = 0;
	double db = 0;
	double dx = 0;
	double mass = 0;
	for(int i = 0; i < (k-1); i++)
		{
		db = 1/fit_dens[i+1] - 1/fit_dens[i];
		dx = fit_x[i+1] - fit_x[i];
		numer = log(1/fit_dens[i+1]) - log(1/fit_dens[i]);
		if(db == 0)
			mass = mass + 1/fit_dens[i] * dx;
		if(db != 0)
			mass = mass + numer * dx/db;
		}
	
		
	if(x < fit_x[0])
		return(0);
	if(x > fit_x[k-1])
		return(0);
	dx = 0;
	double b = 0;
	double g = 0;
	for(int i = 0; i < (k-1); i++)
		{
		if(x <= fit_x[i+1])
			{
			dx = fit_x[i+1] - fit_x[i];
			db = 1/fit_dens[i+1] - 1/fit_dens[i];
			g = 1/fit_dens[i] + db/dx * (x - fit_x[i]);
			return(1/(g*mass) );
			}
		}
	}		
		
	
