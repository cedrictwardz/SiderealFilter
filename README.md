# SiderealFilter

sidereal.py is a short python code that can be used to apply a sidereal filter (SF) to 30 sec position time series according the method described in Twardzik et al. (2019).  

## Instructions for use:

### Summary: 

The code uses 30 sec kinematic position time series (XXXX.k.tdp). It uses the days before a given large earthquake to build the sidereal filter. From experience, 6 days prior the earthquake are needed to build a robust filter. The SF is built by optimally shifting and stacking the days before the earthquake using cross-correlations. By doing that, the aim is to enhance the signature of multipaths and to overcome the issue of using a constant constellation throughout the processing. Once built, the SF is slide along the time series and applied in order to remove the multipaths over the full time series. 

### Input variables: 

The variables [eq] and [evla,evlo,evdp] at the beginning of the script need to be changed manually. It is possible to act like there is no earthquake by setting [eq] to some time in the future, but this has not been tested. [evla,evlo,evdp] are not used.

There are several variables that are hard-coded in the script. In particular, there are safeguards that are used to check if there is a indeed a strong enough sidereal signature (line 137-139). If the criteria are not passed, the filter won't be applied.

### Input files: 

[1] "stations.in" gives the name and location of the receivers. The file is not really useful but it is there anyway and needed. It is formatted the following: 

LATITUDE LONGITUDE NAME

[2] "XXXX.k.tdp" are the position time series. Here is how the file should be formated:

YYYY mm dd HH MM SS.SSS GPS_TIME(S) NORTH(m) EAST(m) UP(m) 

Note that it is assumed that the positions are given with respect to the first epoch in the time series. Thus, the first epoch is always zero.

## Note:
[1] This code has not been optimised and I do not guarantee that it is bug free when using your own time series with their own particularities. Use it at you own risk.

[2] This code has been developed for 30 sec position time series, so I do not guarantee that it provides meaningful results for shorter or longer sampling interval.

[3] Please cite Twardzik et al. (2019) when the code has been used for publication (doi:10.1038/s41598-019-39038-z).

[4] The code does not handle foreshocks and aftershocks. It acts as if there is only one earthquake.

[5] The code uses standard libraries with the exception of the xcorr function which is taken from obspy (doi:10.5281/zenodo.1040770).
