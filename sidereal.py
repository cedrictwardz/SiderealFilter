#!/usr/bin/env python3
  
import numpy as np
import scipy.stats as stats
import scipy.signal as signal
import datetime as datetime
import matplotlib.pyplot as plt
import obspy.signal.cross_correlation as xcorr

def filter(date, time, north, east, up, eq):

    # Find the sampling interval
    dt = np.diff(time).min()

    # Remove the mean
    north -= north[date < eq].mean()
    east  -= east [date < eq].mean()
    up    -= up   [date < eq].mean()

    # Compute a linear trend
    nslope, nintercept, _, _, _ = stats.linregress(time[date<eq],north[date<eq])
    eslope, eintercept, _, _, _ = stats.linregress(time[date<eq],east [date<eq])
    uslope, uintercept, _, _, _ = stats.linregress(time[date<eq],up   [date<eq])

    # Remove the linear trend
    north = north - (nslope*time + nintercept)
    east  = east  - (eslope*time + eintercept)
    up    = up    - (uslope*time + uintercept)

    # Create the full time series
    full_time        = np.arange(time[0], time[-1]+dt, dt)
    full_date        = np.array([date[0] + datetime.timedelta(seconds=i*dt) for i in range(full_time.size)])
    full_north       = np.full(full_time.size, np.nan)
    full_east        = np.full(full_time.size, np.nan)
    full_up          = np.full(full_time.size, np.nan)
    ij, i, j         = np.intersect1d(time, full_time, assume_unique=True, return_indices=True)
    full_north[j]    = north[i]
    full_east [j]    = east [i]
    full_up   [j]    = up   [i]

    # Remove the coseismic offsets
    ib             = np.where(full_date <= eq + datetime.timedelta(seconds=30.0))[0]
    ia             = np.where(full_date >= eq + datetime.timedelta(minutes= 5.5))[0]
    offset_n       = full_north[ia[0]] - full_north[ib[-1]]
    offset_e       = full_east [ia[0]] - full_east [ib[-1]]
    offset_u       = full_up   [ia[0]] - full_up   [ib[-1]]
    full_north[ia] = full_north[ia] - offset_n
    full_east [ia] = full_east [ia] - offset_e
    full_up   [ia] = full_up   [ia] - offset_u
    
    # Remove the gaps in the time series using simple linear interpolation
    nans             = np.isnan(full_north)
    full_north[nans] = np.interp(full_time[nans], full_time[~nans], full_north[~nans])
    full_east [nans] = np.interp(full_time[nans], full_time[~nans], full_east [~nans])
    full_up   [nans] = np.interp(full_time[nans], full_time[~nans], full_up   [~nans])

    # First guess on the max. time shift allowed
    max_shift = np.int(full_time[-1]/86400.0 * 236.0 / dt + 1.5001)

    # Extract the first day
    ntrace  = [1, 1, 1]
    t0, t1  = full_date[0], full_date[0]+datetime.timedelta(days=1)
    index   = (t0 <= full_date) & (full_date < t1)
    n_stack = full_north[index]
    e_stack = full_east [index]
    u_stack = full_up   [index]

    # Start the stack at zero
    t_stack  = np.arange(n_stack.size) * dt
    n_stack -= n_stack[0]
    e_stack -= e_stack[0]
    u_stack -= u_stack[0]
    
    # Stack the days before the earthquake
    while True:

        # Set some logicals
        do_not_n = False
        do_not_e = False
        do_not_u = False

        # Extract the following day
        t0 = t1
        t1 = t0 + datetime.timedelta(days=1)

        # Select the relevant date plus some overlap
        index = (t0-datetime.timedelta(hours=1.0) <= full_date) & (full_date < t1+datetime.timedelta(hours=1.0))

        # Ensure that we are not getting data from after the mainshock
        if np.any(full_date[index] >= eq): break

        # Select the relevant data
        n  = full_north[index]
        e  = full_east [index]
        u  = full_up   [index]

        # Cross-correlation between the current day and the stack
        n_cc = xcorr.correlate(n_stack/ntrace[0], n, n_stack.size)
        e_cc = xcorr.correlate(e_stack/ntrace[1], e, e_stack.size)
        u_cc = xcorr.correlate(u_stack/ntrace[2], u, u_stack.size)

        # Noise level of the cross-correlograms
        n_noise = n_cc.std()
        e_noise = e_cc.std()
        u_noise = u_cc.std()

        # Find the optimal shift to maximize the cross-correlation
        n_shift, n_val = xcorr.xcorr_max(n_cc, abs_max=False)
        e_shift, e_val = xcorr.xcorr_max(e_cc, abs_max=False)
        u_shift, u_val = xcorr.xcorr_max(u_cc, abs_max=False)

        # Check if the trace should be use for stacking
        if (n_shift == 0) or (n_val < 4.0/3.0*n_noise) or (np.abs(n_shift) > max_shift): do_not_n = True
        if (e_shift == 0) or (e_val < 4.0/3.0*e_noise) or (np.abs(e_shift) > max_shift): do_not_e = True
        if (u_shift == 0) or (u_val < 4.0/3.0*u_noise) or (np.abs(u_shift) > max_shift): do_not_u = True

        # Count the number of trace in the stack
        ntrace[0] += (not do_not_n)
        ntrace[1] += (not do_not_e)
        ntrace[2] += (not do_not_u)

        # Extract the relevant portion of the trace
        i0 = np.int(3600.0/dt + 0.5001)
        if not do_not_n: n_stack += (n[i0-n_shift:i0-n_shift+n_stack.size] - n[i0-n_shift])
        if not do_not_e: e_stack += (e[i0-e_shift:i0-e_shift+e_stack.size] - e[i0-e_shift])
        if not do_not_u: u_stack += (u[i0-u_shift:i0-u_shift+u_stack.size] - u[i0-u_shift])
          
    # Remove the mean
    n_stack -= n_stack.mean()
    e_stack -= e_stack.mean()
    u_stack -= u_stack.mean()

    # Calculate a linear trend
    n_slope, n_intercept, _, _, _ = stats.linregress(t_stack, n_stack)
    e_slope, e_intercept, _, _, _ = stats.linregress(t_stack, e_stack)
    u_slope, u_intercept, _, _, _ = stats.linregress(t_stack, u_stack)

    # Remove the linear trend
    n_stack = n_stack - (n_slope*t_stack + n_intercept)
    e_stack = e_stack - (e_slope*t_stack + e_intercept)
    u_stack = u_stack - (u_slope*t_stack + u_intercept)

    # Average all the traces
    n_stack /= ntrace[0]
    e_stack /= ntrace[1]
    u_stack /= ntrace[2]

    # Safeguard in case not enough traces have been used to build the stack
    n_apply = True
    e_apply = True
    u_apply = True
    if ntrace[0] <= 2: print("Not many traces have been used to create the sidereal filter (N)"); n_apply = False
    if ntrace[1] <= 2: print("Not many traces have been used to create the sidereal filter (E)"); e_apply = False
    if ntrace[2] <= 2: print("Not many traces have been used to create the sidereal filter (U)"); u_apply = False

    # Initialize the sidereal filter
    t_filter = full_time.copy()
    n_filter = np.full(full_time.size, np.nan)
    e_filter = np.full(full_time.size, np.nan)
    u_filter = np.full(full_time.size, np.nan)
    ndays    = 0

    # Build the sidereal filter
    while True:
      
        # Extract the relevant dates
        t0    = full_date[0] + datetime.timedelta(days=ndays+0)
        t1    = full_date[0] + datetime.timedelta(days=ndays+1)
        index = (t0 <= full_date) & (full_date < t1)
        t     = full_time [index]
        n     = full_north[index]
        e     = full_east [index]
        u     = full_up   [index]

        # Remove the mean
        n -= n.mean()
        e -= e.mean()
        u -= u.mean()

        # Calculate a linear trend
        n_slope, n_intercept, _, _, _ = stats.linregress(t, n)
        e_slope, e_intercept, _, _, _ = stats.linregress(t, e)
        u_slope, u_intercept, _, _, _ = stats.linregress(t, u)

        # Remove the linear trend
        n = n - (n_slope*t + n_intercept)
        e = e - (e_slope*t + e_intercept)
        u = u - (u_slope*t + u_intercept)

        # Compute the cross-correlation between the given day and the stack
        n_cc = xcorr.correlate(n, n_stack[:n.size], n.size)
        e_cc = xcorr.correlate(e, e_stack[:e.size], e.size)
        u_cc = xcorr.correlate(u, u_stack[:u.size], u.size)

        # Set cross-correlation to zero outside a given range to ensure a reasonable time shift
        n_cc[:n.size-max_shift ] = 0.0
        e_cc[:e.size-max_shift ] = 0.0
        u_cc[:u.size-max_shift ] = 0.0
        n_cc[ n.size+max_shift:] = 0.0
        e_cc[ e.size+max_shift:] = 0.0
        u_cc[ u.size+max_shift:] = 0.0

        # Find the shift that maximize the cross-correlation
        n_shift, _ = xcorr.xcorr_max(n_cc, abs_max=False)
        e_shift, _ = xcorr.xcorr_max(e_cc, abs_max=False)
        u_shift, _ = xcorr.xcorr_max(u_cc, abs_max=False)

        # Create the sidereal filter
        t           = np.arange(e.size)*dt + ndays*86400.0
        ij, i, j    = np.intersect1d(t_filter, t, assume_unique=True, return_indices=True)
        n_filter[i] = (np.roll(n_stack, n_shift) * signal.tukey(n_stack.size, alpha=0.05))[j]
        e_filter[i] = (np.roll(e_stack, e_shift) * signal.tukey(e_stack.size, alpha=0.05))[j]
        u_filter[i] = (np.roll(u_stack, u_shift) * signal.tukey(u_stack.size, alpha=0.05))[j]

        # Exit strategy
        ndays = ndays + 1
        if t1 >= full_date[-1]: break
          
    # Create the full time series
    full_time        = np.arange(time[0], time[-1]+dt, dt)
    full_date        = np.array([date[0] + datetime.timedelta(seconds=i*dt) for i in range(full_time.size)])
    full_north       = np.full(full_time.size, np.nan)
    full_east        = np.full(full_time.size, np.nan)
    full_up          = np.full(full_time.size, np.nan)
    ij, i, j         = np.intersect1d(time, full_time, assume_unique=True, return_indices=True)
    full_north[j]    = north[i]
    full_east [j]    = east [i]
    full_up   [j]    = up   [i]

    # Compute the RMS before sidereal filter
    n_rms0 = north[date < eq].std()
    e_rms0 = east [date < eq].std()
    u_rms0 = up   [date < eq].std()

    # Apply the sidereal filter
    full_up    = full_up    - u_filter
    full_east  = full_east  - e_filter
    full_north = full_north - n_filter
    if u_apply: up   [i] = full_up   [j]
    if e_apply: east [i] = full_east [j]
    if n_apply: north[i] = full_north[j]

    # Compute the RMS after sidereal filter
    n_rms1 = north[date < eq].std()
    e_rms1 = east [date < eq].std()
    u_rms1 = up   [date < eq].std()

    # Print a friendly message
    print('RMS before and after S.F. (north) = %7.2f/%7.2f (%7.2f)' % (n_rms0,n_rms1,(n_rms0-n_rms1)/n_rms0*100.0))
    print('RMS before and after S.F. (east ) = %7.2f/%7.2f (%7.2f)' % (e_rms0,e_rms1,(e_rms0-e_rms1)/e_rms0*100.0))
    print('RMS before and after S.F. (up   ) = %7.2f/%7.2f (%7.2f)' % (u_rms0,u_rms1,(u_rms0-u_rms1)/u_rms0*100.0))

    # All done
    return north, east, up
