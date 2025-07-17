
figFolder = "c:\\Users\\oskarg\\Box Sync\\PostDoc\\My Devfolder\\TimeVaryingPeriodograms\\docs\\src\\figs\\"


using TimeVaryingPeriodograms, Plots, DSP, ProgressBars

gr(grid = true, lw = 1, legendfontsize=11,
    xtickfontsize=9, ytickfontsize=9, xguidefontsize=10, yguidefontsize=10,
    titlefontsize = 14, markersize = 5, markerstrokecolor = :auto, c=:viridis
)


T = 1000; ϕ = [sin.(2π*(1:2:2T)/T).+0.4 -0.8*ones(T)]; θ = zeros(T); σ = 1; μ = 2;
y = simTvARMA(ϕ, θ, 0, σ);
py = plot(y, label="simulated data", c=:black, xlab=time)
pϕ = plot(ϕ[:,1], label="AR parameters", c=:black, xlab="time", ylim=(-0.99,1.9))
plot!(ϕ[:,2], c=:black,label=false)

plot(py,pϕ, size = (1000,350), bottom_margin=4Plots.mm) 
#savefig(figFolder*"dataCoeff.png")
 
ω=0:0.01:π
spectrogram=zeros(T,length(ω))
[spectrogram[t,:] = SpecDensARMAFast(Vector(0:0.01:π), ϕ[t,:], θ[t], 1.) for t in 1:1000]
heatmap(1:T, ω, log.(spectrogram'), c=:viridis, clim=(-4,2), 
       xlab=time,ylab="Frequency", title="log spectrogram", size=(600,400))
#savefig(figFolder*"spectrogram.png")


m=20 
MI, λ = mvPeriodogram(y, m, 1)
λr = repeat(λ*π, (T-2m)÷m + 1)
scatterRaw = scatter(m+1:T-m, λr, marker_z=log.(MI), label =false, 
        title = "raw", ylab="frequency"
)



MIBC, λ = BoundCorrectedMvPeriodogram(y, m, 2, 1)
scatterBC = scatter(m+1:T-m, λr, marker_z=log.(MIBC), label =false, 
        title = "boundary corrected", xlab="time"
)



MIPreW, λ = mvPrewhitePeriodogram(y, m, 2, 1)
scatterPreW = scatter(m+1:T-m, λr, marker_z=log.(MIPreW), label =false, 
        title = "prewhitened", xlab="time", ylab="frequency"
)


MITap, λ = mvTapPeriodogram(y, m, 1, 1) 
scatterTap = scatter(m+1:T-m, λr, marker_z=log.(MITap),  label =false, 
        title = "tapered"
)


#savefig(figFolder*"mvPeriodograms.png")

plot(scatterRaw,scatterTap,scatterPreW,scatterBC, layout = (2,2), clim=(-4,2),
     size = (1000,600), left_margin=4Plots.mm
)

N=2m+1
S=20
BW = tvPeriodogram(y, N, S)
tapBW = tvPeriodogram(y, N, S, nothing, hanning)
preWBW =  preWhiteTvBlockPeriodogram(y,N,S,2)
BCBW =  tapBoundCorrectTvBlockPeriodogram(y,N,S,2,nothing)




hBW = heatmap(m+1:S:T-m, λ*π,  log.(BW[1][2:end,:]), title="raw", ylab="frequency")
hTap = heatmap(m+1:S:T-m, λ*π,  log.(tapBW[1][2:end,:]), title = "tapered")
hPreW = heatmap(m+1:S:T-m, λ*π,  log.(preWBW'[2:end,:]),title="prewhitened", xlab="time", ylab="frequency")
hBC = heatmap(m+1:S:T-m, λ*π,  log.(BCBW'[2:end,:]), title = "boundary corrected", xlab="time")
plot(hBW,hTap,hPreW,hBC, layout = (2,2), size = (1000,600), 
    left_margin=4Plots.mm, background_color=:darkgrey, clim=(-4,2)
)
#savefig(figFolder*"blockPeriodograms.png")


preP = prePeriodogram(y,λ*π)

hPreP = heatmap(1:T-1,λ*π, real(preP'), clim=(-1,2), title="original")


everittPeriodogram = EverittI(y, λ*π , 0.4, 0.92)
hEverit = heatmap(1:T,λ*π, real(everittPeriodogram)', clim=(-0.1,0.4), title="Everitt et al.") 

plot(hPreP,hEverit, layout = (1,2), size = (1000,300), left_margin=4Plots.mm, 
     bottom_margin=4Plots.mm, background_color=:darkgrey, 
     xlab="time", ylab="frequency", clim=(-0.5,3)
)


#savefig(figFolder*"prePeriodograms.png")


