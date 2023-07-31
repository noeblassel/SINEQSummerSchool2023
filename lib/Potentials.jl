module Potentials

    export mueller_brown, grad_mueller_brown, LJClusterInteraction2D,lj_energy,lj_grad,lj_energy_threaded,lj_grad_threaded,entropic_switch,grad_entropic_switch
    using Threads, Parameters

### Mueller-Brown potential


    const A_mb = [-200,-100,-170,15]
    const a_mb = [-1,-1,-6.5,0.7]
    const b_mb = [0,0,11,0.6]
    const c_mb = [-10,-10,-6.5,0.7]
    const x0_mb = [1,0,-0.5,-1]
    const y0_mb = [0,0.5,1.5,1]


mueller_brown(x,y) = sum(@. A_mb*exp(a_mb*(x-x0_mb)^2 + b*(x-x0_mb)*(y-y0_mb)+c*(y-y0_mb)^2))
mueller_brown(X) = mueller_brown(X...)

function grad_mueller_brown(x,y)
    v = @. A_mb*exp(a_mb*(x-x0_mb)^2 + b_mb*(x-x0_mb)*(y-y0_mb)+c_mb*(y-y0_mb)^2)
    return [sum(@. (2a_mb*(x-x0_mb)+b_mb*(y-y0_mb))*v), sum(@. (2c_mb*(y-y0_mb)+b_mb*(x-x0_mb))*v) ]
end

grad_mueller_brown(X) = grad_mueller_brown(X...)

Base.@kwdef struct MuellerBrownProperties
    minima = []
    saddles = []
    normal_hyperplanes = []
end

### Lennard-Jones Cluster

Base.@kwdef struct LJClusterInteraction2D{N}
    σ=1.0
    σ6=σ^6
    σ12=σ6^2
    ε=1.0
    α=1.0 # sharpness of harmonic confining potential
end

function lj_energy(X,inter::LJClusterInteraction2D{N}) where {N}
    V = 0.0
    for i=1:N-1
        for j=i+1:N
            inv_r6=inv(sum(abs2,X[:,i]-X[:,j]))^3
            V += (inter.σ12*inv_r6-inter.σ6)*inv_r6
        end
    end
    return 4inter.ε*V+inter.α*sum(abs2,X)/2 # add confining potential
end

function lj_grad(X,inter::LJClusterInteraction2D{N}) where {N}
    F = zeros(2,N)
    for i=1:N-1
        r = zeros(2)
        f = zeros(2)
        for j=i+1:N
            r = X[:,i]-X[:,j]
            inv_r2 = inv(sum(abs2,r))
            inv_r4 = inv_r2^2
            inv_r8 = inv_r4^2
            
            f = (6inter.σ6 - 12inter.σ12*inv_r2*inv_r4)*inv_r8*r

            F[:,i] .+= f
            F[:,j] .-= f
        end
    end

    return 4inter.ε*F + inter.α*X #add confining potential
end

function lj_energy_threaded(X,inter::LJClusterInteraction2D{N}) where {N}
    V_threaded = zeros(nthreads())
    @threads for i=1:N-1
        for j=2:N
            inv_r6=inv(sum(abs2,X[:,i]-X[:,j]))^3
            V_threaded[threadid()] += (inter.σ12*inv_r6^2-inter.σ6*inv_r6)
        end
    end
    return 4inter.ε*sum(V_threaded)+inter.α*sum(abs2,X)/2 # add confining potential
end

function lj_grad_threaded(X,inter::LJClusterInteraction2D{N}) where {N}
    F_threaded = fill(zeros(2,N),nthreads())
    @threads for i=1:N-1
        r = zeros(2)
        f = zeros(2)
        for j=i+1:N
            r = X[:,i]-X[:,j]
            inv_r2 = inv(sum(abs2,r))
            inv_r4 = inv_r2^2
            inv_r8 = inv_r4^2
            
            f = (6inter.σ6 - 12inter.σ12*inv_r2*inv_r4)*inv_r8*r

            F_threaded[threadid()][:,i] .+= f
            F_threaded[threadid()][:,j] .-= f
        end
    end

    return 4inter.ε*sum(F_threaded) + inter.α*X #add confining potential
end

### Entropic switch potential

    const es_saddles = [-0.6172723078764598 0.6172723078764598 0.0; 1.1027345175080963 1.1027345175080963 -0.19999999999972246]
    const es_pos_eigvecs = [0.6080988038706289 -0.6080988038706289 -1.0; 0.793861351075306 0.793861351075306 0.0]
    const es_minima = [-1.0480549928242202 0.0 1.0480549928242202; -0.042093666306677734 1.5370820044494633 -0.042093666306677734]
    

    function entropic_switch(x, y)
        tmp1 = x^2
        tmp2 = (y - 1 / 3)^2
        return 3 * exp(-tmp1) * (exp(-tmp2) - exp(-(y - 5 / 3)^2)) - 5 * exp(-y^2) * (exp(-(x - 1)^2) + exp(-(x + 1)^2)) + 0.2 * tmp1^2 + 0.2 * tmp2^2
    end

    function entropic_switch(q)
        return entropic_switch(q...)
    end

    function grad_entropic_switch(x, y)

        tmp1 = exp(4*x)
        tmp2 = exp(-x^2 - 2*x - y^2 - 1)
        tmp3 = exp(-x^2)
        tmp4 = exp(-(y-1/3)^2)
        tmp5 = exp(-(y-5/3)^2)

        dx = 0.8*x^3 + 10*(tmp1*(x - 1) + x + 1)*tmp2 - 6*tmp3*x*(tmp4 - tmp5)

        dy = 10*(tmp1 + 1)*y*tmp2 + 3*tmp3*(2*tmp5*(y - 5/3) - 2*tmp4*(y - 1/3)) + 0.8*(y - 1/3)^3

        return [dx, dy]
    end

    function grad_entropic_switch(q)
        return grad_entropic_switch(q...)
    end




end