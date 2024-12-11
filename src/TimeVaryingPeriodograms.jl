module TimeVaryingPeriodograms


using Distributions, DSP, Polynomials, FFTW, ProgressBars, LinearAlgebra, ToeplitzMatrices,DataFrames
# Write your package code here.
include("DSSDGPs.jl")
export wikleTV6AR, DahlhausAR2, simTvARMA

include("DSSInfer.jl")
export mvPeriodogram, mvTapPeriodogram, mvPrewhitePeriodogram, mvPrewhiteTapPeriodogram, BoundCorrectedMvPeriodogram, tapBoundCorrectedMvPeriodogram,
tvPeriodogram, prePeriodogram, predictiveDFT, mvPeriodogram,  completeTapPeriodogram, tapBoundCorrectTvBlockPeriodogram,whitePeriodogram,
preWhiteTvBlockPeriodogram, prePeriodogram, EverittI, SpecDensARMA,SpecDensARMAFast,
rawPeriodogram, transFunkAR, BoundCorrectTvBlockPeriodogram,
completeDFTSingle, SingleDft2, whitePeriodogram, EverittI, completePeriodogram,reScaleTapTvBlockPeriodogram


end
