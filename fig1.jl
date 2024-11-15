## FIGURE 1

using DifferentialEquations
using Plots
using LinearAlgebra

include("defaults.jl")

## Figure 1(a): ODE model

    # Two parameter sets
    a,b,c,d = 1.0,-1.0,1.0,-1.0
    B = [a b; c d]
    A₁ = 1.0*I(2)
    A₂ = A₁ + B

    # Initial condition
    X₀ = [-b/a; 1]
    F₁ = [t -> (exp(A₁ * t) * X₀)[i] for i = 1:2]
    F₂ = [t -> (exp(A₂ * t) * X₀)[i] for i = 1:2]

    fig1a = plot(F₁[1],xlim=(0.0,2.5),c=:red,lw=2.0,label="x₁(t)")
    plot!(fig1a,F₂[1],xlim=(0.0,2.5),c=:black,ls=:dash,lw=2.0)
    plot!(fig1a,xlim=(0.0,2.0),widen=true,ylim=(0.0,8.0))
    plot!(fig1a,xlabel="t",ylabel="x(t)")

## Figure 1(b): PDE model

    # Two parameter sets
    c₁,d₁ = 1.0,0.05
    c₂ = 2.0
    d₂ = d₁ + (c₂ - c₁) / π^2

    # Eigenvalues, eigenfunctions
    v(n,x) = sin((2n - 1) * π * x)
    λ(n;r=r,d=d) = r - d * (2n - 1)^2 * π^2

    # Initial condition
    N = 4; C = zeros(N); C[1] = 1.0
    f₀(x) = sum(C[n] * v(n,x) for n = 1:length(C))

    # Solution
    pde_sol(r,d,x,t) = sum(C[n] * exp(λ(n;r,d) * t) * v(n,x) for n = 1:N)

    # Plot
    t = range(0.0,2.0,6)
    F₁ = [x -> pde_sol(c₁,d₁,x,tᵢ) for tᵢ in t]
    F₂ = [x -> pde_sol(c₂,d₂,x,tᵢ) for tᵢ in t]

    fig1b = plot(F₁,xlim=(0.0,1.0),xwiden=true,palette=reverse(palette(:Blues_7))[2:end],lw=2.0,label=t')
    plot!(fig1b,F₂,xlim=(0.0,1.0),c=:black,ls=:dash,lw=2.0,label=["p2" fill("",1,length(t)-1)])
    plot!(fig1b,
        xlabel="x", ylabel="u(x,t)",yticks=0:0.5:3.0,ylim=(0.0,3.0),ywiden=true
    )

## Figure 1
fig1 = plot(fig1a,fig1b,size=(600,200))
savefig("$(@__DIR__)/fig1.svg")