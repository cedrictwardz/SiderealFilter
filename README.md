# SiderealFilter

sidereal.py is a short python code that can be used to apply a sidereal filter (SF) to 30 sec position time series according the method described in Twardzik et al. (2019).  

## Instructions for use:

### Summary: 

The code uses 30 sec kinematic position time series (XXXX.k.tdp). It takes the positions before a given large earthquake to build the sidereal filter. From experience, 6 days prior the earthquake are needed to build a robust filter. The SF is built by optimally shifting and stacking the days before the earthquake using cross-correlations. By doing that, the aim is to enhance the signature of multipaths and overcome the issue of using a constant constellation throughout the processing. Once built, the SF slides along the time series and is applied in order to remove the multipaths over the full time series. 

### Input variables: 

It is possible to act like there is no earthquake by setting [eq] to some time in the future.

### Input files: 

Note that it is assumed that the positions are given with respect to the first epoch in the time series. Thus, the first epoch is always zero.

## Note:
[1] This code has not been optimised and I do not guarantee that it is bug free when using your own time series with their own particularities. Use it at you own risk.

[2] This code has been developed for 30 sec position time series, so I do not guarantee that it provides meaningful results for shorter or longer sampling interval.

[3] Please cite Twardzik et al. (2019) when the code has been used for publication (doi:10.1038/s41598-019-39038-z).

[4] The code does not handle foreshocks and aftershocks. It acts as if there is only one earthquake.

[5] The code uses standard libraries with the exception of the xcorr function which is taken from obspy (doi:10.5281/zenodo.1040770).
