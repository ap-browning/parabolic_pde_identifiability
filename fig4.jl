## FIGURE 4

using Distributions
using Plots, StatsPlots
using Optim
using QuadGK
using Random

include("defaults.jl")

# True parameter values
c = 1.0
d = 0.05

# Eigenfunctions and eigenvalues
v(n,x) = sin(n * π * x)
λ(n;c=c,d=d) = c - d * n^2 * π^2

# Identifiable function
lrfunc = ld -> log(λ(1) + exp(ld) * π^2)

# Number of terms in expansion(s)
N = 8

# Model solution
function solve_model(c,d,C)
    (x,t) -> sum(C[n] * exp(λ(n;c,d) * t) * v(n,x) for n = 1:N)
end

# Data grid
X = range(0.0,1.0,11)[2:end-1]
T = range(0.0,2.0,11)

# Noise parameters (known)
η = 10.0
σ = 0.2
dnoise = MvNormal(σ^2 * [exp(-η * abs(xᵢ - xⱼ)) for xᵢ in X, xⱼ in X])

# Initial condition (Gaussian)
function u₀(x,ω)
    f_unscaled(x) = pdf(Normal(0.5,ω),x) - pdf(Normal(0.5,ω),0.0)
    f_unscaled(x) / f_unscaled(0.5)
end

# Initial condition (coefficients)
function init_cond_coef(ω,N=N)
    f = x -> u₀(x,ω)
    [2quadgk(x -> f(x) * v(n,x),0.0,1.0)[1] for n = 1:N]
end

# Profile likelihood surface for each ω
lc_grid = range(-1.0,2.0,100)
ld_grid = range(-6.0,0.0,101)
F = Array{Matrix}(undef,3)

# Look through each ω
ω = [0.1,0.2,0.3]
for (i,ωᵢ) = enumerate(ω)

    # Parameter set
    C = init_cond_coef(ωᵢ)

    # Generate data
    ldata = [solve_model(c,d,C)(x,t) for x in X, t in T] + 
                rand(MersenneTwister(),dnoise,length(T))

    # Full likelihood
    function loglike(lp)
        c,d = exp.(lp[1:2])
        C = lp[3:end]
        u = solve_model(c,d,C)
        M = [u(x,t) for x in X, t in T]
        loglikelihood(dnoise, M - ldata)
    end

    # Profile likelihood
    function logproflike(lc,ld)
        func(C) = -loglike([lc;ld;C])
        res = optimize(func, C, LBFGS(); autodiff = :forward)
        -res.minimum
    end

    # Maximise over a grid
    F[i] = [logproflike(lc,ld) for lc = lc_grid, ld = ld_grid]
    F[i] = F[i] .- maximum(F[i])

end


## Figure 4

row1 = [plot(title="ω = $ωᵢ") for ωᵢ in ω]
row2 = [plot(title="ω = $ωᵢ") for ωᵢ in ω]

# Options
    # Colorscheme    
    cs = range(RGBA(1,1,1),RGBA(0,0,1),15)
    # Lowest likelihood value to plot    
    thresh = -6.0
    # 95% CI threshold
    ci_thresh = quantile(Chisq(2),0.95) / 2

# Create plots
for (i,ωᵢ) = enumerate(ω)
    # Row 1
    C = init_cond_coef(ωᵢ)
    plot!(row1[i],x -> u₀(x,ωᵢ),xlim=(0.0,1.0),lw=2.0,c=:black,label="IC")
    plot!(row1[i],x -> sum(C[n] * v(n,x) for n = 1:N),c=:red,lw=2.0,ls=:dash,label="IC match")
    plot!(row1[i],x -> C[1] * v(1,x),c=:blue,lw=2.0,label="v₁",xlabel="x",ylabel="u(x,0)")

    # Row 2
        # Likelihood surface
        contourf!(row2[i],ld_grid,lc_grid,F[i],c=cs,lw=0.0,clim=(thresh,0.0),fα=0.8)
        # 95% confidence interval
        contour!(row2[i],ld_grid,lc_grid,F[i],levels=[-ci_thresh],lw=2.0,c=:black)
        plot!([],[],c=:black,lw=2.0,label="95% CI")
        # A
        plot!(row2[i],lrfunc,lw=2.0,c=:black,ls=:dash,label="A")
        # True value
        scatter!([log(d)],[log(c)],c=:orange,m=:diamond,label="True",xlabel="log d",ylabel="log c")
end

# Save figure
fig4 = plot(row1...,row2...,layout=grid(2,3),size=(800,450),colorbar=false)
add_plot_labels!(fig4)
savefig(fig4,"fig4.svg")

# # Save colorbar
# fig4 = plot(row1...,row2...,layout=grid(2,3),size=(800,450))
# add_plot_labels!(fig4)
# savefig(fig4,"fig4_cbar.svg")

