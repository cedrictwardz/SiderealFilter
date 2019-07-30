# SiderealFilter

sidereal.py is a short python code that can be used to apply a sidereal filter (SF) to 30 sec position time series according the method described in Twardzik et al. (2019).  

## Instructions for use:

### Summary: 

The code uses 30 sec kinematic position time series (XXXX.k.tdp). It uses the days before a given large earthquake to build the sidereal filter. From experience, 6 days prior the earthquake are needed to build a robust filter. The SF is built by shifting and stacking the days before the earthquake with the aim that it will enhance the signature of multipaths. Once done, the SF is slide along the time series in order to remove the multipaths over the full time series. 

### Input variables: 

The variables [eq] and [evla,evlo,evdp] at the beginning of the script need to be changed manually. It is possible to act like there is no earthquake by setting [eq] to some time in the future, but this has not been tested.

There are several variables that are hard-coded in the script. In particular, there are safeguards that are used to check if there is a indeed a sidereal signature (line 137-139).

### Input files: 

[1] "stations.in" gives the name and location of the receivers. It is not really useful but it is there anyway. The file is formatted the following: 

LATITUDE LONGITUDE NAME

[2] "XXXX.k.tdp" are the position time series. Here is how the file should be formated:

YYYY mm dd HH MM SS.SSS GPS_TIME(S) NORTH(m) EAST(m) UP(m) 

Note that it is assumed that the positions are given with respect to the first epoch in the time series

## Note:
[1] This code has not been optimised and I do not guarantee that it is bug free. Use it at you own risk.

[2] This code has been developed for 30 sec position time series, so I do not guarantee that it provides meaningfull results from shorter or longer sampling interval.

[3] Please cite Twardzik et al. (2019) when the code has been used for publication (doi:10.1038/s41598-019-39038-z).

[4] The code do not handle potential foreshocks and aftershokcs. It acts as if there is only one earthquake.

[5] The code uses standard library with the exception of the xcorr function which is taken from obspy (doi:10.5281/zenodo.1040770)
