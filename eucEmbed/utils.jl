randVec = function(;lower::Real = -1.0, upper::Real = 1.0, n::Int = 1)
    unif_rands = rand(n)
    rng = upper - lower
    ans = lower .+ unif_rands .* rng
    return(ans)
end

expit = function(x::Real)
    ans = exp(x) / (1.0 + exp(x))
    return(ans)
end
