#!/usr/bin/env python3

import numpy as np
import scipy.stats as stats
import scipy.signal as signal
import datetime as datetime
import matplotlib.pyplot as plt
import obspy.signal.cross_correlation as xcorr

def filter(date, time, north, east, up, eq, showme=False):
    """
    @author: Cedric Twardzik
    @contact: cedric.twardz(at)gmail.com
    @inputs: date   [type datetime]: dates of the positions
             time   [type float]: times of the positions (s)
             north  [type float]: position in the north component
             east   [type float]: position in the east component
             up     [type float]: position in the vertical component
             eq     [type datetime]: date of the mainshock 
             showme [type logical]: show the stacks and the filter 
    """

    # Find the sampling interval
    dt = time[1] - time[0]

    # Remove the coseismic offset
    ib         = date < eq
    ia         = date > eq
    offset_n   = north[ia][0] - north[ib][-1]
    offset_e   = east [ia][0] - east [ib][-1]
    offset_u   = up   [ia][0] - up   [ib][-1]
    north[ia] -= offset_n
    east [ia] -= offset_e
    up   [ia] -= offset_u

    # Create the full time series
    full_time     = np.arange(time[0], time[-1]+dt, dt)
    full_date     = np.array([date[0] + datetime.timedelta(seconds=i*dt) for i in range(full_time.size)])
    full_north    = np.full(full_time.size, np.nan)
    full_east     = np.full(full_time.size, np.nan)
    full_up       = np.full(full_time.size, np.nan)
    ij, i, j      = np.intersect1d(time, full_time, assume_unique=True, return_indices=True)
    full_north[j] = north[i]
    full_east [j] = east [i]
    full_up   [j] = up   [i]

    # Sidereal day (Choi et al. 2004)
    sidereal_shift = datetime.datetime(2000, 1, 2, 0, 0, 0) - datetime.datetime(2000, 1, 1, 23, 56, 4, 0)
    sidereal_shift = np.int(sidereal_shift.total_seconds()/dt + 0.5001)

    # Maximum sidereal shift allowed
    ndays           = (full_date[-1] - full_date[0]).days
    sidereal_shift *= (ndays*4)

    # Remove the gaps in the time series using simple linear interpolation
    nans             = np.isnan(full_north)
    full_north[nans] = np.interp(full_time[nans], full_time[~nans], full_north[~nans])
    full_east [nans] = np.interp(full_time[nans], full_time[~nans], full_east [~nans])
    full_up   [nans] = np.interp(full_time[nans], full_time[~nans], full_up   [~nans])

    # Extract the first day
    tstart       = full_date[0]
    tstop        = full_date[0] + datetime.timedelta(days=1)
    index        = (tstart <= full_date) & (full_date < tstop)
    stack_north  = full_north[index] * full_north[index].std()
    stack_east   = full_east [index] * full_east [index].std()
    stack_up     = full_up   [index] * full_up   [index].std()
    stack_time   = full_time [index]
    stack_size   = stack_time.size
    ntrace_north = full_north[index].std()
    ntrace_east  = full_east [index].std()
    ntrace_up    = full_up   [index].std()
    
    # Show the stacks
    if showme:
        plt.close()
        fig, ax = plt.subplots(3, 1, sharex='col')
        ax[0].plot(stack_time, stack_north/ntrace_north, 'k-', lw=1.0, alpha=0.5)
        ax[1].plot(stack_time, stack_east /ntrace_east , 'k-', lw=1.0, alpha=0.5)
        ax[2].plot(stack_time, stack_up   /ntrace_up   , 'k-', lw=1.0, alpha=0.5)

    # Initialize the day counter
    iday = 1

    # Stack the days before the earthquake
    while True:

        # Extract the following day
        tstart = tstop
        tstop  = tstart + datetime.timedelta(days=1)
        index  = (tstart-datetime.timedelta(hours=1.0) <= full_date) & (full_date < tstop+datetime.timedelta(hours=1.0))
        i0     = np.int(3600.0/dt + 0.5001)

        # Ensure that we are not getting any data from after the mainchock
        if np.any(full_date[index] >= eq): break

        # Select the relevant data
        n = full_north[index]
        e = full_east [index]
        u = full_up   [index]

        # Cross-correlation between the current day and the stack (position)
        #cc_north = xcorr.correlate(stack_north/ntrace_north, n, stack_size)
        #cc_east  = xcorr.correlate(stack_east /ntrace_east , e, stack_size)
        #cc_up    = xcorr.correlate(stack_up   /ntrace_up   , u, stack_size)

        # Cross-correlation between the current day and the stack (velocity)
        cc_north = xcorr.correlate(np.diff(stack_north/ntrace_north), np.diff(n), stack_size-1)
        cc_east  = xcorr.correlate(np.diff(stack_east /ntrace_east ), np.diff(e), stack_size-1)
        cc_up    = xcorr.correlate(np.diff(stack_up   /ntrace_up   ), np.diff(u), stack_size-1)

        # Find the optimal shift to maximize the cross-correlation
        shift_north, ccmax_north = xcorr.xcorr_max(cc_north, abs_max=False)
        shift_east , ccmax_east  = xcorr.xcorr_max(cc_east , abs_max=False)
        shift_up   , ccmax_up    = xcorr.xcorr_max(cc_up   , abs_max=False)

        # Show the stacks
        if showme:
            ax[0].plot(stack_time, n[i0-shift_north:i0+stack_size-shift_north], 'k-', lw=1.0, alpha=0.5)
            ax[1].plot(stack_time, e[i0-shift_east :i0+stack_size-shift_east ], 'k-', lw=1.0, alpha=0.5)
            ax[2].plot(stack_time, u[i0-shift_up   :i0+stack_size-shift_up   ], 'k-', lw=1.0, alpha=0.5)

        # Add the trace to the stack
        stack_north += n[i0-shift_north:i0+stack_size-shift_north] * n[i0-shift_north:i0+stack_size-shift_north].std()
        stack_east  += e[i0-shift_east :i0+stack_size-shift_east ] * e[i0-shift_east :i0+stack_size-shift_east ].std()
        stack_up    += u[i0-shift_up   :i0+stack_size-shift_up   ] * u[i0-shift_up   :i0+stack_size-shift_up   ].std()

        # Update the normalization factor
        ntrace_north += n[i0-shift_north:i0+stack_size-shift_north].std()
        ntrace_east  += e[i0-shift_east :i0+stack_size-shift_east ].std()
        ntrace_up    += u[i0-shift_up   :i0+stack_size-shift_up   ].std()

        # Update the day counter
        iday += 1

    # Normalize the final stack
    stack_north /= ntrace_north
    stack_east  /= ntrace_east
    stack_up    /= ntrace_up
    
    # Show the final stack
    if showme:
        ax[0].plot(stack_time, stack_north, 'r-', lw=2.0, alpha=0.9)
        ax[1].plot(stack_time, stack_east , 'r-', lw=2.0, alpha=0.9)
        ax[2].plot(stack_time, stack_up   , 'r-', lw=2.0, alpha=0.9)
        plt.show()

    # Initialize the sidereal filter
    filter_date  = full_date.copy()
    filter_time  = full_time.copy()
    filter_north = np.full(filter_time.size, np.nan)
    filter_east  = np.full(filter_time.size, np.nan)
    filter_up    = np.full(filter_time.size, np.nan)

    # Initialize the day counter
    iday = 0

    # Build the sidereal filter
    while True:

        # Extract the relevant dates
        tstart = full_date[0] + datetime.timedelta(days=iday)
        tstop  = full_date[0] + datetime.timedelta(days=iday+1)
        index  = (tstart <= full_date) & (full_date < tstop)

        # Ensure we have a full date
        if index.sum() != 2880: break

        # Remove the log-trend before cross-correlation
        n = full_north[index]
        e = full_east [index]
        u = full_up   [index]

        # Cross-correlation between the current day and the stack (position)
        #cc_north = xcorr.correlate(stack_north, n, stack_size)
        #cc_east  = xcorr.correlate(stack_east , e, stack_size)
        #cc_up    = xcorr.correlate(stack_up   , u, stack_size)

        # Cross-correlation between the current day and the stack (velocity)
        cc_north = xcorr.correlate(np.diff(stack_north), np.diff(n), stack_size-1)
        cc_east  = xcorr.correlate(np.diff(stack_east ), np.diff(e), stack_size-1)
        cc_up    = xcorr.correlate(np.diff(stack_up   ), np.diff(u), stack_size-1)

        # Find the optimal shift to maximize the cross-correlation
        shift_north, ccmax_north = xcorr.xcorr_max(cc_north, abs_max=False)
        shift_east , ccmax_east  = xcorr.xcorr_max(cc_east , abs_max=False)
        shift_up   , ccmax_up    = xcorr.xcorr_max(cc_up   , abs_max=False)

        # Insert the stack
        filter_north[index] = np.roll(stack_north * signal.tukey(stack_size, 0.05), -shift_north)
        filter_east [index] = np.roll(stack_east  * signal.tukey(stack_size, 0.05), -shift_east )
        filter_up   [index] = np.roll(stack_up    * signal.tukey(stack_size, 0.05), -shift_up   )

        # Update the day counter
        iday += 1

    # Remove the mean of the filter
    filter_north -= np.nanmean(filter_north)
    filter_east  -= np.nanmean(filter_east )
    filter_up    -= np.nanmean(filter_up   )
    
    # Remove the sidereal filter
    full_north -= filter_north
    full_east  -= filter_east
    full_up    -= filter_up

    # Get the relevant part of the filter
    ij, i, j = np.intersect1d(time, full_time, assume_unique=True, return_indices=True)
    north[i] = full_north[j]
    east [i] = full_east [j]
    up   [i] = full_up   [j]

    # Add the coseismic offset
    north[ia] += offset_n
    east [ia] += offset_e
    up   [ia] += offset_u

    # All done
    return north, east, up
