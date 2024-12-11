"""
```
    tvPeriodogram(y, N, S, nFreq = nothing, taper = nothing; taperArgs...)
```
Computes the block-periodogram (Dahlhaus, 2012) of the time series `y` using a window size `N`, and step-size `S` on the Fourier frequencies. NFreq can be used to evaluate other frequencies and taper can be used to apply a taper to the data, see package DSP for details. 

Returns a spectrogram of the data. 
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia>  tvPeriodogram(y, 25, 15)
```
""" 
function tvPeriodogram(y::AbstractVector, N, S, nFreq = nothing, taper = nothing; taperArgs...)
    if isnothing(taper) 
      win_normalized = nothing
    else
      win = taper(Int64(N); taperArgs...);
      win_normalized = win/sqrt(mean(win.^2))
    end
    if isnothing(nFreq) 
      numFrequencies = N
    else
      numFrequencies = nextfastfft(2*nFreq)
    end
    noverlap = N - S
    specG = spectrogram(y, N, noverlap ; onesided=true,  nfft=numFrequencies, fs=1, window = win_normalized);
    
   # spectrogram uses frequencies (-0.5;0.5) then they truncate the spectral density to (0,0.5)
   # and double it. So we need to revert the doubling and strech it to (0,π) ==> 2*2π=4π. 
   return specG.power ./ 4π , 2π*specG.freq
end



""" 
```
  mvPeriodogram(y, m, startFreq) 
```
Computes the moving local periodogram of the time series `y` using a window size `m`, starting at frequency `startFreq`. 

Returns the local periodogram data MI_t data and a vector λ with the frequencies where MI_t is computed at frequency λ[1 + (t + startFreq - 2) mod m)].
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI, λ = mvPeriodogram(y, 10, 2)
```
""" 
function mvPeriodogram(y, m, startFreq)
  modulus(t) = 1 + mod(t + startFreq - 2, m) 
  T = length(y) - 2*m
  normConst = 1/(2π*(2m+1))
  λ = (2*(1:m))/(2m+1)
  MI = zeros(T)
  ν = 0:2m
  for t in 1:T
    MI[t] = abs(sum(y[t:(t+2m)].*exp.(-im*π*ν*λ[modulus(t)])))^2
  end
  MI = normConst*MI
  return MI, λ
end



""" 
```  
mvTapPeriodogram(y, m, α=1, startFreq=1) 
```
Same as mvPeriodogram but with a Tukey taper applied to the data. Use α=1 for Hanning taper.

Returns the local periodogram data MI_t data and a vector λ with the frequencies where MI_t is computed at frequency λ[1 + (t + startFreq - 2) mod m)].
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI, λ = mvTapPeriodogramDiffStartFreq(y, 10)
```
""" 
function mvTapPeriodogram(y, m, α=1, startFreq=1) # Add argument for α instead?
  modulus(t) = 1 + mod(t + startFreq - 2, m) 
  T = length(y) - 2*m
  tap = tukey(2m+1,α)
  normalizedTaper = tap ./ sqrt.(mean(tap.^2))
  normConst = 1/(2π*(2m+1))
  λ = (2*(1:m))/(2m+1)
  MI = zeros(T)
  ν = 0:2m
  for t in 1:T
    yTap = normalizedTaper .* y[t:(t+2m)]
    #yTap = rescaleVar(yTap, y[t:(t+2m)]) 
    MI[t] = abs(sum(yTap .* exp.(-im*π*ν*λ[modulus(t)])))^2
  end
  MI = normConst*MI
  return MI, λ
end



""" 
```
  mvPrewhitePeriodogram(y, m, p=nothing, startFreq) 
```
Same as mvPeriodogram but with a Tukey taper applied to the data. Use α=1 for Hanning taper.

Returns the local periodogram data MI_t data and a vector λ with the frequencies where MI_t is computed at frequency λ[1 + (t + startFreq - 2) mod m)].
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI, λ = mvTapPeriodogramDiffStartFreq(y, 10)
```
""" 
function mvPrewhitePeriodogram(y, m, p=nothing, startFreq=1)
  modulus(t) = 1 + mod(t + startFreq - 2, m) 
  T = length(y) - 2*m
  normConst = 1/(2π*(2m+1))
  λ = (2*(1:m))/(2m+1)
  MI = zeros(T)
  ν = 0:2m
  for t in 1:T
    x, ϕ, σ = prewhite(y[t:(t+2m)], p)
    # Need to be careful with the transfer function here. What freq interval is used, double or single side.
    MI[t] = abs(sum(x .* exp.(-im*π*ν*λ[modulus(t)])))^2 .* transFunkAR(λ[modulus(t)]*π, ϕ)*2π
  end
  MI = normConst*MI
  return MI, λ
end



""" 
```
  mvPrewhiteTapPeriodogram(y, m, p, startFreq) 
```
Same as mvPrewhitePeriodogram but with prewhitening. p is the lag length used in the AR-model, if left blanc the AR order is estimated on each segment.

Returns the local periodogram data MI_t data and a vector λ with the frequencies where MI_t is computed at frequency λ[1 + (t + startFreq - 2) mod m)].
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI, λ = mvPreWPeriodogram(y, 10, 2, 1)
```
""" 
function mvPrewhiteTapPeriodogram(y, m, p=nothing, startFreq=1)
  modulus(t) = 1 + mod(t + startFreq - 2, m) 
  T = length(y) - 2*m
  α = 1 # note, α=1 is same as hanning, α=0 is rectangular.
  tap = tukey(2m+1,α)
  normalizedTaper = tap ./ sqrt.(mean(tap.^2))
  normConst = 1/(2π*(2m+1))
  λ = (2*(1:m))/(2m+1)
  MI = zeros(T)
  ν = 0:2m
  for t in 1:T
    x, ϕ, σ = prewhite(y[t:(t+2m)], p)
    xTap = normalizedTaper .* x
    # Need to be careful with the transfer function here. What freq interval is used, double or single side.
    MI[t] = abs(sum(xTap .* exp.(-im*π*ν*λ[modulus(t)])))^2 .* transFunkAR(λ[modulus(t)]*π, ϕ)*2π
  end
  MI = normConst*MI
  return MI, λ
end



""" 
```
  predictiveDFT(y, ω; p=nothing, ϕ=nothing)
```
Predictive DFT from Rao and Yang (2021), see r-package "cspec" on CRAN. Default is fitting the best AR model and compute the predictive DFT based on the predictions att frequencies ´ω´. ´p´ will fix the AR-order and ´ϕ´ gives the AR-coefficients.


Returns the discrete Fourier coefficients associated with the predictions.
# Examples
```julia-repl
julia> rand(Normal(500))
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI, λ = mvPreWPeriodogram(y, 10, 2, 1)
```
""" 
function predictiveDFT(y, ω; p=nothing, ϕ=nothing)
  T=length(y)
  m = length(ω)

  if isnothing(ϕ)
      if isnothing(p)
        p = whichLag(lagSelectARBurg(y, 2, 6), "HQ")
      end
      ϕ = -lpc(y, p, LPCBurg())[1]
  end
  
  if !isnothing(ϕ) p = length(ϕ) end
  ϕfunc = zeros(Complex,m)
  if (length(ϕ) == 0) 
      pred = 0
  end
  if length(ϕ) > 0
      for w ∈ 1:m 
          temp = exp.(-im * (1:p) * ω[w])
          ϕfunc[w] = 1 - sum(ϕ .* temp)
      end
      JL = zeros(Complex, m)
      JR = zeros(Complex, m)
      for w ∈ 1:m 
          tempL = zeros(Complex, p)
          tempR = zeros(Complex, p)
          for ell ∈ 1:p
              tempL[ell] = sum(ϕ[ell:p] .* exp.(-im*(0:(p - ell)).*ω[w]))          
              tempR[ell] = sum(ϕ[ell:p] .* exp.(im*(1:(p - ell + 1)).*ω[w]))
          end
          JL[w] = sum(y[1:p] .* tempL)
          JR[w] = sum(reverse(y)[1:p] .* tempR)
      end
      JL = (JL ./ ϕfunc) ./ sqrt(T)
      JR = (JR ./ conj(ϕfunc)) ./ sqrt(T)
      pred = JL + JR
  end
  return pred
end



""" 
```
    BoundCorrectedMvPeriodogram(y, m, p, startFreq) 
```
Same as mvPeriodogram but with boundary correction in Rao and Yang (2021). ´p´ is the lag length used in the AR-model, if left blanc the AR order is estimated on each segment.

Returns the local periodogram data MI_t data and a vector λ with the frequencies where MI_t is computed at frequency λ[1 + (t + startFreq - 2) mod m)].
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI, λ = BoundCorrectedMvPeriodogram(y, 10, 2, 1)
```
""" 
function BoundCorrectedMvPeriodogram(y, m, p=nothing, startFreq=1)
  modulus(t) = 1 + mod(t + startFreq - 2, m) 
  T = length(y) - 2*m
  n = 2m+1
  normConst = 1/(2π*n)
  λ = (2*(1:m))/n
  MI = zeros(T)
  ν = 0:2m
  for t in 1:T
      #boundCorrectDFT = sum(y[t:(t+2m)].*exp.(-im*π*ν*λ[modulus(t)])) + sqrt(n)*predictiveDFT(y[t:(t+2m)], π*λ[modulus(t)], p)[1]
      compDFT = completeDFTSingle(y[t:(t+2m)], π*λ[modulus(t)],; p, ϕ = nothing)*sqrt(n)
      regularDFT = SingleDft2(y[t:(t+2m)], π*λ[modulus(t)])*sqrt(n)
      MI[t] = real(compDFT*conj(regularDFT))
  end
  MI = normConst*MI
  return max.(MI,0.1^8), λ
end


""" 
  tapBoundCorrectedMvPeriodogram(y, m, p, startFreq) 

Same as ´BoundCorrectedMvPeriodogram´ but with tapering using a Hanning window.

Returns the local periodogram data MI_t data and a vector λ with the frequencies where MI_t is computed at frequency λ[1 + (t + startFreq - 2) mod m)].
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI, λ = tapBoundCorrectedMvPeriodogram(y, 10, 2, 1)
```
""" 
function tapBoundCorrectedMvPeriodogram(y, m, p=nothing, startFreq=1)
    modulus(t) = 1 + mod(t + startFreq - 2, m) 
    T = length(y) - 2*m
    n = 2m+1
    α = 1 # note, α=1 is same as hanning, α=0 is rectangular.
    tap = tukey(n,α)
  
    normalizedTaper = tap ./ sqrt.(mean(tap.^2))
    normConst = 1/(2π*n)
    λ = (2*(1:m))/n
    MI = zeros(T)
    ν = 0:2m
    for t in 1:T
        yTap = normalizedTaper .* y[t:(t+2m)]
        yTap = rescaleVar(yTap, y[t:(t+2m)]) 
        #boundCorrectDFT = sum(y[t:(t+2m)].*exp.(-im*π*ν*λ[modulus(t)])) + sqrt(n)*predictiveDFT(y[t:(t+2m)], π*λ[modulus(t)], p)[1]
        compDFT = completeDFTSingle(y[t:(t+2m)], π*λ[modulus(t)],; p, ϕ = nothing)*sqrt(n)
        regularDFT = SingleDft2(yTap, π*λ[modulus(t)])*sqrt(n)
        MI[t] = real(compDFT*conj(regularDFT))
    end
    MI = normConst*MI
    return max.(MI,0.1^6), λ
end


#""" 
#  fft2(y)
#
#Inverse Fourier transform following Rao and Yang (2021), see r-package "cspec" on CRAN. I find it a bit strange that they use the inverse Fourier transform here.
#
#Returns the inverse Fourier transform of the data.
#""" 
function fft2(y)
  n=length(y)
  newy = vcat(y[end], y[1:end-1])
  dftTemp = DSP.ifft(newy) .* n
  dft = vcat(dftTemp[2:end], dftTemp[1])/sqrt(n)
  return dft
end


""" 
   SingleDft2(y, fr)

Computing the DFT on a single frequency, `fr`. This is needed in the boundary correction for the dynamic Whittle likellihood..

Returns the discrete Fourier transform of the data.
"""
function SingleDft2(y, fr)
  n=length(y)
  dftmat = exp.(-im .* Vector(fr*(1:n)))
  dft = (dftmat' * y) / sqrt(n)
  return dft
end


""" 
```
   completeDFTSingle(y, fr; p = nothing, ϕ = nothing)
```
Computing the complete DFT on a single frequency, `fr`. This is needed in the boundary correction for the dynamic Whittle likellihood.

Returns the discrete Fourier transform of the data.
"""
function completeDFTSingle(y, fr; p = nothing, ϕ = nothing)
  SingleDft2(y, fr)[1] + predictiveDFT(y, fr; p, ϕ)[1]
end



""" 
```
   completePeriodogram(y; fr=nothing, p = nothing, ϕ = nothing, α=nothing, threshold=0.1^8)
```
Computing the complete periodogram on frequencies, `fr` (Fourier frquencies if `nothing`). This is needed in the boundary correction for the dynamic Whittle likellihood.

Returns the complete periodogram of the data given the model.
"""
function completePeriodogram(y; fr=nothing, p = nothing, ϕ = nothing, α=nothing, threshold=0.1^8)
  n=length(y)
  if isnothing(fr)
      if iseven(n)
          fr = (0:n÷2)/(n÷2)*π 
      else
          fr = (0:(n-1)/2)/(n/2) .* π
      end
  end
  if !isnothing(α)
      tap = tukey(n,α)
      normalizedTaper = tap ./ sqrt.(mean(tap.^2))
      yTap = normalizedTaper .* y
      yTap = y
      Ri =SingleDft2.(Ref(yTap), fr)
  else
      Ri = SingleDft2.(Ref(y), fr)
  end
  Le = completeDFTSingle.(Ref(y), fr; p, ϕ)
  return max.(real(Le .* conj(Ri))/2π, threshold)
end


""" 
```
    tapBoundCorrectTvBlockPeriodogram(y,N,S,p=nothing,α=nothing,thresH=0.1^8)
```
Same as tvPeriodogram but with boundary correction from Rao and Yang (2021). ´p´ is the lag length used in the AR-model, if left blanc the AR order is estimated on each segment.

Returns the local periodogram data MI_t at the Fourier frequencies.
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI = tapBoundCorrectTvBlockPeriodogram(y,25,10,p=nothing,α=nothing,thresH=0.1^8)
```
""" 
function tapBoundCorrectTvBlockPeriodogram(y,N,S,p=nothing,α=nothing,thresH=0.1^8)
  M = (length(y) - N)÷S
  tFirst = N÷2+1
  if isodd(N)
      MI = zeros(M+1, N÷2+1)
      for t in 1:M+1
        MI[t,:] = completePeriodogram(y[(t-1)*S+1:tFirst + (t-1)*S + N÷2]; fr=nothing, p, ϕ = nothing, α)
      end      
  else   
      MI = zeros(M+1, N÷2+1)
      for t in 1:M+1
          MI[t,:] = completePeriodogram(y[(t-1)*S+1:tFirst + (t-1)*S + N÷2 - 1]; fr=nothing, p, ϕ = nothing, α)[1:N÷2+1]
      end
  end 

  return MI
end



""" 
```
    preWhiteTvBlockPeriodogram(y,N,S,p=nothing)
```
Same as tvPeriodogram but with pre-whitening. ´p´ is the lag length used in the AR-model, if left blanc the AR order is estimated on each segment.

Returns the local periodogram data MI_t at the Fourier frequencies.
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI = preWhiteTvBlockPeriodogram(y,25,10,p=nothing)
```
""" 
function preWhiteTvBlockPeriodogram(y,N,S,p=nothing)
  M = (length(y) - N)÷S
  tFirst = N÷2+1
  MI = zeros(M+1, N÷2+1)
  if isodd(N)
      for t in 1:M+1
        MI[t,:] = whitePeriodogram(y[(t-1)*S+1:tFirst + (t-1)*S + N÷2], p)
      end      
  else   
      for t in 1:M+1
          MI[t,:] = whitePeriodogram(y[(t-1)*S+1:tFirst + (t-1)*S + N÷2 - 1], p)
      end
  end 

  return MI
end



""" 
```
    reScaleTapTvBlockPeriodogram(y,N,S)
```
Same as tvPeriodogram but with rescaling of the tapered data to match the variance of the non-scaled data.
Returns the local periodogram data MI_t at the Fourier frequencies.
# Examples
```julia-repl
julia> T = 500; ϕ = [sin.(2π*(1:T)/T) -0.1*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
julia> y = simTvARMA(ϕ, θ, μ, σ);
julia> MI = reScaleTapTvBlockPeriodogram(y,25,10)
```
""" 
function reScaleTapTvBlockPeriodogram(y,N,S,α=1)
  M = (length(y) - N)÷S
  tFirst = N÷2+1
  MI = zeros(M+1, N÷2+1)
  tap = tukey(N,α)
  normalizedTaper = tap ./ sqrt.(mean(tap.^2))

  if isodd(N)
      MI = zeros(M+1, N÷2)
      for t in 1:M+1
        ySeg = y[(t-1)*S+1:tFirst + (t-1)*S + N÷2]
        yTap = normalizedTaper .* ySeg
        yTap = rescaleVar(yTap, ySeg) 
        MI[t,:] = whitePeriodogram(yTap, p)
      end      
  else   
      MI = zeros(M+1, N÷2+1)
      for t in 1:M+1
        ySeg = y[(t-1)*S+1:tFirst + (t-1)*S + N÷2 - 1]
        yTap = normalizedTaper .* ySeg
        yTap = rescaleVar(yTap, ySeg) 
          MI[t,:] = whitePeriodogram(yTap, p)
      end
  end 

  return MI
end



###############################
# Functions to compute the exact likelihood for AR-processes, used to select lag length in each segment.
#""" 
#```
#    arAutocov(ϕ,σ²)
#```
#compute the theoretical p first autocovariances of an AR-process with parameters `ϕ` and `σ²`.
## Examples
#```julia-repl
#julia> arAutocov([1.25; -0.5], 4)
#```
#""" 
function arAutocov(ϕ,σ²)
  p = length(ϕ)
  F = vcat(ϕ', I(p)[1:end-1,:])
  σ²*inv(I(p^2) - kron(F,F))[1:p,1]
end


function logLikFirstPobs(y,c,ϕ,σ²)
  p=length(ϕ)
  Γ = Matrix(SymmetricToeplitz(arAutocov(ϕ,σ²)))
  μ = ones(p)*c/(1-sum(ϕ))
  logpdf(MvNormal(μ,Γ), y[1:p])
end

function exactARlogLik(y, c, ϕ, σ²)
  p = length(ϕ)
  T = length(y)
  # Initialize with the log likelihood for the initial observations.
  logLik = logLikFirstPobs(y,c,ϕ,σ²)
  for t in p+1:T
      μ = c + ϕ'*reverse(y[t-p:t-1]) 
      logLik += logpdf(Normal(μ, sqrt(σ²)), y[t])
  end
  return logLik
end

# Checking whether we are in the stationary region.
function is_stationary(ϕ)
  poly = Polynomial([1, -ϕ...])  
  roots_poly = roots(poly)  
  return all(abs.(roots_poly) .> 1), roots_poly  
end  

function miniFun(param, p=3, y=x)
  c = param[1]; ϕ=param[2:p+1]; σ²=exp(param[p+2])  
  if !is_stationary(ϕ)
      return 10^10#Inf
  end 
  -exactARlogLik(y, c, ϕ, σ²)
end 

BIC(k,T,logLik) = k*log(T) - 2*logLik
HQ(k,T,logLik) = 2k*log(log(T)) - 2logLik
AIC(k,logLik) = 2k -2*logLik

function lagSelectAR(y, pmin, pmax)
  BICp = []
  AICp = []
  HQp = []
  for i in pmin:pmax
      optFun(θ) = miniFun(θ, i, y)
      # Maximize the log-likelihood function
      initial_guess = vcat(0.5, 0.1*rand(i), 1.0)
      logLik = -optimize(optFun, initial_guess, LBFGS()).minimum
      push!(AICp, AIC(2+i, logLik))
      push!(BICp, BIC(2+i, length(y), logLik))
      push!(HQp,   HQ(2+i, length(y), logLik))
  end
  df = DataFrame(hcat(Vector(pmin:pmax), AICp,BICp,HQp), :auto)
  rename!(df, ["lags", "AIC", "BIC", "HQ"])
  return df
end

# Burgs method is orders of magnitude faster, and should be good for short time series.
function lagSelectARBurg(y, pmin, pmax) 
  BICp = []
  AICp = []
  HQp = []
  for i in pmin:pmax
      Burg = lpc(y, i, LPCBurg())
      ϕ = -Burg[1]
      σ² = Burg[2]
      logLik = exactARlogLik(y, 0, ϕ, σ²)
      push!(AICp, AIC(2+i, logLik))
      push!(BICp, BIC(2+i, length(y), logLik))
      push!(HQp,   HQ(2+i, length(y), logLik))
  end
  df = DataFrame(hcat(Vector(pmin:pmax), AICp,BICp,HQp), :auto)
  rename!(df, ["lags", "AIC", "BIC", "HQ"])
  return df
end

function whichLag(lSARIn, IC)
  lSARIn."lags"[argmin(lSARIn[!, Symbol(IC)])] 
end

# SKA VI TVINGA INTERCEPTET ATT VARA NOLL?
function ARCoeffExactMLE(y, p)
    optFun(θ) = miniFun(θ, p, y)
    # Maximize the log-likelihood function
    initial_guess = vcat(0.5, 0.1*rand(p), 1.0)
    return optimize(optFun, initial_guess, LBFGS()).minimizer[2:p+1]
end

function ARCoeffAndConstExactMLE(y, p)
  optFun(θ) = miniFun(θ, p, y)
  # Maximize the log-likelihood function
  initial_guess = vcat(0.5, 0.1*rand(p), 1.0)
  res = optimize(optFun, initial_guess, LBFGS()).minimizer[1:p+1]
  return res[1], res[2:p+1] 
end


# Prewhitening correction
function prewhite(y, p=nothing)

  if isnothing(p)
    p = whichLag(lagSelectARBurg(y, 2, 6), "HQ")
  end

  X = hcat([y[p-pp+1:end-pp] for pp in 1:p]...)
  #ϕHat = inv(X'X)*X'y[p+1:end]  
  #c, ϕHat = ARCoeffAndConstExactMLE(y, p)
  #res = y[p+1:end] .- c - X*ϕHat

  #ϕHat = ARCoeffExactMLE(y, p)
  ϕHat = -lpc(y, p, LPCBurg())[1]
  res = y[p+1:end] - X*ϕHat
  # compute the first p residuals.
  σHat = std(res)
  res = vcat(rand(Normal(0, σHat), p), res) 
  return res, ϕHat, σHat
end

function rawPeriodogram(y)
  T=length(y)
  f = fft(y)
  p = abs.(f[1:T÷2+1]).^2
  # Here I do normalization into (-0.5, 0.5) but only take the positive freq.
  # This makes it comparable to the DSP.jl periodograms but it might not be relevant 
  # in this project...
  normP = 2*p ./ T
  normP /=4π
  #normP[end] = normP[end]*2
  #normP[1] = normP[1]*2
  normP
end



""" 
    whitePeriodogram(y, p)

Returns the periodogram from the prewhitening approach using AR-order p.

""" 
function whitePeriodogram(y, p)
  T = length(y)
  y_white, ϕHat, σhat = prewhite(y, p)
  pWhite = rawPeriodogram(y_white)
  ω = ((0:T÷2)/(T/2))*π
  spdWhite = SpecDensARMA.(ω, Ref(ϕHat), [0.0], 1)*2π
  pWhite.*spdWhite
end


######################
# Preperiodograms.
""" 
    prePeriodogram(y, ω) 

ω is a selected vector of frequencies. Returns a (length(ω) x T-1) matrix containing the prePeriodogram.

""" 
function prePeriodogram(y, ω)
    T = length(y)
    #j = 1:Int(round((T-0.95)/2));
    j=1:length(ω)
    #ω = (2*π .*j) ./ T;
#    preI = complex(zeros(Int((T)/2), T-1));
    preI = complex(zeros(length(ω), T-1));

    # Note, one needs to start Julia with multiple threads to make use of the multithreading.
    Threads.@threads for t in ProgressBar(2:(T-1)) 
      #for j in 1:Int((T)/2)
      for j in 1:length(ω)
    
        tempFreq = 0.0
        if iseven(t)
          for k in Int(max(-((t-1)*2 - 1), -((T-t-1)*2 - 1))):Int(min(((t-1)*2 - 1), ((T-t-1)*2 - 1)))
            if iseven(k)
              tempFreq += (1/2π)*y[Int(t + abs(k)/2)]*y[Int(t - abs(k)/2)]*exp(-im*k*ω[j])  
            end
    
            if !iseven(k)
              tempFreq += (1/2π)*(y[Int(t + 0.5 + abs(k)/2)]*y[Int(t + 0.5 - abs(k)/2)] + 
                                   y[Int(t - 0.5 + abs(k)/2)]*y[Int(t - 0.5 - abs(k)/2)])*0.5*exp(-im*k*ω[j])  
            end
          end
        end
    
        # Does this even differ from when t is even? Remove?
        if !iseven(t)
          for k in Int(max(-((t-1)*2 - 1), -((T-t-1)*2 - 1))):Int(min(((t-1)*2 - 1),((T-t-1)*2 - 1)))
            if iseven(k)
              tempFreq += (1/2π)*y[Int(t + abs(k)/2)]*y[Int(t - abs(k)/2)]*exp(-im*k*ω[j])  
            end
    
            if !iseven(k)
              tempFreq += (1/2π)*(y[Int(t + 0.5 + abs(k)/2)]*y[Int(t + 0.5 - abs(k)/2)] + 
                                   y[Int(t - 0.5 + abs(k)/2)]*y[Int(t - 0.5 - abs(k)/2)])*0.5*exp(-im*k*ω[j])  
            end
          end
        end
    
        preI[j, t] = tempFreq
    
      end
    end
    return(preI)                                                                                
end



# Weighted moving average - to be used on preperiodogram to smooth over time.
# x - time series to be smoothed.
# vector of weights.
function moveAvg(x, w)
  len = length(w)
  y = copy(x)
  for i ∈ len÷2+1 : length(x)-len÷2
    y[i] = sum(x[i-len÷2:i+len÷2] .* w)
  end
  return y
end


#########################
# Everitt et al. version of the pre-periodogram for a single frequency.

function EverittIω(y, ω, ρ, ν)
  T = length(y);
  μ = zeros(T); μ[1] = y[1];
  A = zeros(T)im; A[1] = y[1]*y[2];
  B = zeros(T)im; B[1] = y[1]*exp(-im*ω)
  Σ = zeros(T);
  yBar = zeros(T); 
  for t ∈ 1:T-1
      μ[t+1] = (1 - ρ)*μ[t] + ρ*y[t+1];
      yBar[t+1] = y[t+1] - μ[t+1]
      A[t+1] = (1 - ρ)A[t] + ρ*yBar[t+1]*B[t]
      B[t+1] = ν*exp(-im*ω)*(B[t] + yBar[t+1])
      Σ[t+1] = (1 - ρ)*Σ[t] + ρ*yBar[t+1]^2
  end
  return (1/2π)*(Σ + 2A)  
end


# Everitt et al. version of the pre-periodogram.
# y - input time series, ωgrid - grid of frequencies to be evaluated, ρ - Parameter that controls smoothing ove time,
# ν - parameter that controlls exponential taper.
""" 
    EverittI(y, ωgrid, ρ, ν)

Everitt et al. version of the pre-periodogram. 
- `y` - input time series, 
- `ωgrid` - grid of frequencies to be evaluated,  
- `ρ` - Parameter that controls smoothing ove time,
- `ν` - parameter that controlls exponential taper.

""" 
function EverittI(y, ωgrid, ρ, ν)
  numω = length(ωgrid)
  andI = zeros(length(y), numω)im
  for i ∈ 1:numω   
      andI[:, i] .= EverittIω(y, ωgrid[i], ρ, ν)
  end
  return andI
end


""" 
InterpolateSpectrogram(blockSpec, N, S)

Takes a segmented locally stationary spectral density (as we get from our inference), blockSpec, with the number of 
observations used in each segment, N, and the time distance between the segments centers, S. 
Outputs an interpolated spectrogram with the same number of frequencies, but with T time-points.  
# Examples
```julia-repl
julia> ;
```
"""
# Interpolate spectrogram from block-Whittle estimates
function InterpolateSpectrogram(blockSpec, N, S)
    numFreq, M = size(blockSpec)
    T = S*(M-1) + N
    specInterpolate = zeros(T, numFreq)
    for i ∈ 1:M
        if i == 1
            specInterpolate[1:N÷2 + S÷2, :] .= blockSpec'[1, :]'
        else
            specInterpolate[(i-1)*S + N÷2 - S÷2:(i-1)*S + N÷2 + S÷2, :] .= blockSpec'[i, :]'
        end
      
    end
    specInterpolate[(M-1)*S + N÷2 + S÷2 + 1:end, :] .= blockSpec'[M, :]'
    return specInterpolate'
end

# Moved it here from spectools to reduce on dependency on the clusters.
""" 
    SpecDensARMA(ω, ϕ, θ, σ²) 

Compute spectral density for the univariate ARMA model. 

- ω is a radial frequency
- ϕ is a vector of AR coefficients
- θ is a vector of MA coefficients
- σ² is the noise variance

# Examples
The spectral density for an AR(1) process with unit noise variance is
```doctests 
julia> SpecDensARMA(0.5, 0.9, 0, 1)
0.6909224383713601
```
""" 
function SpecDensARMA(ω, ϕ, θ, σ²)
    ARpoly =  Polynomial([1;-ϕ], :z)
    MApoly =  Polynomial([1;θ], :z) 
    specDens = (σ²/(2π))*(abs(MApoly(exp(-im*ω)))^2/abs(ARpoly(exp(-im*ω)))^2)
	return specDens
end 



function SpecDensARMAFast(ω::AbstractVector, ϕ, θ, σ²)
  ARpoly = Polynomial([1;-ϕ], :z)
  MApoly = Polynomial([1;θ], :z) 
  specDens = (σ²/(2π))*(abs.(MApoly.(exp.(-im*ω))).^2 ./ abs.(ARpoly.(exp.(-im*ω))).^2)
  return specDens   
end

transFunkAR(ω, ϕ) = SpecDensARMA(ω, ϕ,[0], 1) 

# Helpfunction to rescale the variance of a vector x to the variance of a vector y.
function rescaleVar(x, y)
  return x .* std(y) ./ std(x)
end