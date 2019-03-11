mutable struct AdamInfo
    b1::Real
    b2::Real
    b1t::Real
    b2t::Real
    eps::Real

    k::Int

    mt::Array{Real, 1}
    vt::Array{Real, 1}
    delta::Array{Real, 1}

    AdamInfo(nPars::Int) = new(0.9, 0.999, 1.0, 1.0, 0.00001,
                                    nPars, zeros(nPars), zeros(nPars),
                                    zeros(nPars))

end

updateDelta! = function(grad, alpha::Real, ai::AdamInfo)
    # Updating b1t and b2t due to new iteration
    ai.b1t *= ai.b1
    ai.b2t *= ai.b2

    # Updating mt and vt
    for i in 1:ai.k
        this_grad = grad[i]
        ai.mt[i] = ai.b1 * ai.mt[i] + (1.0 - ai.b1) * this_grad
        ai.vt[i] = ai.b2 * ai.vt[i] + (1.0 - ai.b2) * this_grad * this_grad

        hat_mt = ai.mt[i] / (1.0 - b1t)
        hat_vt = ai.vt[i] / (1.0 - b2t)

        ai.delta[i] = - alpha * hat_mt / (sqrt(hat_vt) + ai.eps)
    end
end

getDelta = function(ai::AdamInfo)
    return ai.delta
end
