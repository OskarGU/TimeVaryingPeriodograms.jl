# Periodograms
Time varying periodograms that can be used as data input in dynamic models that uses spectral likelihoods.

## Moving periodograms
Based on the periodogram in Häfner and Kirch (2017) - "Moving Fourier Analysis for Locally Stationary Processes with the Bootstrap in View". Consider adjustments in form of tapering, pre-whitening and boundary correction.

```@docs
mvPeriodogram
mvTapPeriodogram
mvPrewhitePeriodogram
mvPrewhiteTapPeriodogram
BoundCorrectedMvPeriodogram
tapBoundCorrectedMvPeriodogram
```

## Block periodograms
Based on the periodogram used in block-Whittle likelihood, see e.g. Dahlhaus (2012) - "Locally Stationary Processes" using the same adjustments as above.
```@docs
tvPeriodogram
preWhiteTvBlockPeriodogram
tapBoundCorrectTvBlockPeriodogram
reScaleTapTvBlockPeriodogram
```


## Pre-periodograms
Preperiodogram based on Neumann and Von Sachs (1997) and Everitt et al. (2013).
```@docs
prePeriodogram
EverittI
```

# Time invariant periodograms
Time invariant periodograms that are used as parts of the above. Can be interesting for comparison. 
```@docs
whitePeriodogram
completePeriodogram
```