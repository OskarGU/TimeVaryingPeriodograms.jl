# Data generating processes used to compare methods





""" 
    simTvARMA(ϕ, θ, μ, σ) 

Simulate a time-varying ARMA process 

- ϕ is a T×p matrix with time evolutions of the p AR coefficients  
- θ is a T×q matrix with time evolutions of the q MA coefficients  
- σ is a T×1 vector with time evolutions of the standard deviation of the noise
- μ is a T×1 vector with time evolutions of the mean of the process

# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> plot(y)
```
""" 
function simTvARMA(ϕ, θ, μ, σ)
    T, p = size(ϕ)
    q = 1
    if length(size(θ)) > 1
      q = size(θ)[2]
    end 
    lagMax = max(p,q)
    if length(σ) > 1
        e = [rand(Normal(0, σ[t])) for t ∈ 1:T+q]
    else 
        e = rand(Normal(0, σ), T+lagMax)
    end
  
    y = zeros(T+p)
    y[1:p] = rand(Normal(μ[1], σ[1]), p)
    if length(μ) > 1
        μ = [μ[1]*ones(p);μ]
    else
        μ = μ*ones(T+p)
    end

    for t ∈ p+1:T+p
        y[t] = μ[t-p] + ϕ[t-p,:]'*reverse(y[t-p:t-1] .- μ[t-p]) + 
            vcat(1, θ[t-p, :])'*reverse(e[t-q:t])
    end
    return y[end-T+1:end]
end


""" 
    wikleTV6AR(T, σ, A, μ)

Simulate spectrogram and time series from a time-varying autoregressive model of order 6.

Length of time series is `T`, innovation standard deviation is `σ`, and the AR coefficients are compute from the three-element vector `A`. The mean of the time series is `μ`.

# Examples
```julia-repl
julia> T = 1024; σ = 2; A = [1.1, 1.12, 1.1]; μ = 0;
julia> y, spdMat, coeffϕ = wikleTV6AR(T, σ, A, μ)
julia> heatmap(log.(spdMat'), c = :viridis, clim = (-3,1.5), 
    size = (1000,400), title="True Spectrogram");
julia> plot(y, label = nothing, title = "Simulated time series from TvAR(6) model", 
        xlabel = "time", ylabel = "y")
```
""" 
function wikleTV6AR(T, σ, A, μ)

    A₁, A₂, A₃ = A
    θ₁ = zeros(T)
    θ₂ = zeros(T)
    θ₃ = zeros(T)

    for t ∈ 1:T
        θ₁[t] = 0.05 .+ (0.1/(T-1)) * t
        θ₂[t] = 0.25
        θ₃[t] = 0.45 .- (0.1/(T-1)) * t
    end

    B₁ = A₁ .* exp.(2π .* im .* θ₁)
    B₂ = A₂ .* exp.(2π .* im .* θ₂)
    B₃ = A₃ .* exp.(2π .* im .* θ₃)

    coeffϕ = zeros(T,6);
    for t ∈ 1:T
        pol = fromroots([B₁[t], conj(B₁[t]),
                    B₂[t], conj(B₂[t]),
                    B₃[t], conj(B₃[t])]).coeffs
        coeffϕ[t,:] = -pol[2:end]/pol[1]
    end

    # True log spectrogram
    ω = 0:0.01:π;
    spdMat = zeros(T, length(ω))
    for t ∈ 1:T
        spdMat[t,:] .= SpecDensARMA.(ω, Ref(coeffϕ[t,:]), [0.0], σ^2)
    end

    y = simTvARMA(coeffϕ, zeros(T), μ, σ)

    return y, spdMat, coeffϕ

end


""" 
    DahlhausAR2(T, σ, μ)

Simulate spectrogram and time series from a time-varying autoregressive model of order 2 from Dahlhaus handbook chapter..

Length of time series is `T`, innovation standard deviation is `σ`. The mean of the time series is `μ`(default is 0).

# Examples
```julia-repl
julia> T = 1024; σ = 2; μ = 0;
julia> y, spdMat, coeffϕ = DahlhausAR2(T, σ, μ)
julia> heatmap(log.(spdMat'), c = :viridis, clim = (-3,1.5), 
    size = (1000,400), title="True Spectrogram");
julia> plot(y, label = nothing, title = "Simulated time series from TvAR(2) model", 
        xlabel = "time", ylabel = "y")
```
""" 
function DahlhausAR2(T, σ, μ=0)
    ϕ₁ = -1.8*cos.(1.5.-cos.(4π*range(0,1,length=T)));
    ϕ₂ = -0.9;
    ϕ = hcat(ϕ₁, ones(T)*ϕ₂);
    # True log spectrogram
    ω = 0:0.001:π;
    spdMat = zeros(T, length(ω))
    for t ∈ 1:T
        spdMat[t,:] .= SpecDensARMA.(ω, Ref(ϕ[t,:]), [0.0], σ^2)
    end
    y = simTvARMA(ϕ, zeros(T), μ, σ)
    return y, spdMat, ϕ
end