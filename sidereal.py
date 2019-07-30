#!/usr/bin/env python2
  
import numpy                          as np
import glob                           as glob
import datetime                       as dt
import scipy.stats                    as stats
import scipy.signal                   as signal
import matplotlib.pyplot              as plt
import obspy.signal.cross_correlation as xcorr

# --------- #
# Mainshock #
# --------- #
eq             = dt.datetime(2015,9,16,22,54,32,860000)
evla,evlo,evdp = -31.573,-71.674,22.4
stalist        = open('stations.in','r').readlines()

# ------------------------- #
# Loop over the time series #
# ------------------------- #
for kfile in glob.glob('*.k.tdp'):
    # -------------------- #
    # Load the time series #
    # -------------------- #
    ts    = np.loadtxt(kfile)
    nt    = ts.shape[0]
    delta = ts[1,6] - ts[0,6]
    date  = np.asarray([dt.datetime(np.int(ts[i,0]),np.int(ts[i,1]),np.int(ts[i,2]),np.int(ts[i,3]),np.int(ts[i,4]),np.int(ts[i,5])) for i in range(nt)])
    time  = ts[:,6] - ts[0,6]
    north = ts[:,7] * 1000.0
    east  = ts[:,8] * 1000.0
    up    = ts[:,9] * 1000.0
    stnm  = kfile.split('/')[-1].split('.')[0]
    # --------------------------------------- #
    # Extract the coordinates of the receiver #
    # --------------------------------------- #
    stla,stlo = [[np.float(stalist[i].split()[0]),np.float(stalist[i].split()[1])] for i in range(len(stalist)) if stnm in stalist[i].split()[2]][0]
    # --------------- #
    # Remove the mean #
    # --------------- #
    north = north - np.mean(north[date<eq])
    east  = east  - np.mean(east [date<eq])
    up    = up    - np.mean(up   [date<eq])
    # ------------------------ #
    # Calculate a linear trend #
    # ------------------------ #
    nslope,nintercept,r_value,p_value,std_err = stats.linregress(time[date<eq],north[date<eq])
    eslope,eintercept,r_value,p_value,std_err = stats.linregress(time[date<eq],east [date<eq])
    uslope,uintercept,r_value,p_value,std_err = stats.linregress(time[date<eq],up   [date<eq])
    # --------------------- #
    # Remove a linear trend #
    # --------------------- #
    north = north - (nslope*time + nintercept)
    east  = east  - (eslope*time + eintercept)
    up    = up    - (uslope*time + uintercept)
    # ----------------------------------------- #
    # Find the last point before the earthquake #
    # ----------------------------------------- #
    ib    = np.where(date <= eq+dt.timedelta(seconds=30.0))[0]
    nlast = north[ib[-1]]
    elast = east [ib[-1]]
    ulast = up   [ib[-1]]
    # -------------------------------- #
    # Remove gaps from the time series #
    # -------------------------------- #
    fulltime         = np.arange(time[0],time[-1]+delta,delta)
    fulldate         = np.asarray( [date[0] + dt.timedelta(seconds=i*delta) for i in range(fulltime.size)] )
    fullnorth        = np.ones(fulltime.size) * nlast
    fulleast         = np.ones(fulltime.size) * elast
    fullup           = np.ones(fulltime.size) * ulast
    index            = np.where(np.isin(fulltime,time))[0]
    fullnorth[index] = north.copy()
    fulleast [index] = east .copy()
    fullup   [index] = up   .copy()
    # ------------------------------------------------------ #
    # Remove the co-seismic offset from the full time series #
    # ------------------------------------------------------ #
    ib            = np.where(fulldate <= eq+dt.timedelta(seconds=30.0))[0]
    ia            = np.where(fulldate >= eq+dt.timedelta(minutes= 5.5))[0]
    offset_north  = fullnorth[ia[0]] - fullnorth[ib[-1]]
    offset_east   = fulleast [ia[0]] - fulleast [ib[-1]]
    offset_up     = fullup   [ia[0]] - fullup   [ib[-1]]
    fullnorth[ia] = fullnorth[ia] - offset_north
    fulleast [ia] = fulleast [ia] - offset_east
    fullup   [ia] = fullup   [ia] - offset_up
    # ------------------------------- #
    # Guess the maximum shift allowed #
    # ------------------------------- #
    max_shift = int(fulltime[-1]/86400.0) * 236.0
    max_shift = int(max_shift/delta+1.5001)
    # --------------------- #
    # Extract the first day #
    # --------------------- #
    ntrace = [1,1,1]
    t0,t1  = fulldate[0],fulldate[0]+dt.timedelta(days=1)
    nstack = fullnorth[(fulldate >= t0) & (fulldate < t1)]
    estack = fulleast [(fulldate >= t0) & (fulldate < t1)]
    ustack = fullup   [(fulldate >= t0) & (fulldate < t1)]
    # ----------------------- #
    # Start the stack at zero #
    # ----------------------- #
    nstack = nstack - nstack[0]
    estack = estack - estack[0]
    ustack = ustack - ustack[0]
    # ------------------------------------- #
    # Stack the days before the earthquakes #
    # ------------------------------------- #
    while True:
        # -------- #
        # Logicals #
        # -------- #
        donot_north = False
        donot_east  = False
        donot_up    = False
        # -------------------- #
        # Extract the next day #
        # -------------------- #
        t0 = t1
        t1 = t0+dt.timedelta(days=1)
        # *** ENSURE WE ARE NOT GETTING DATA FROM AFTER THE EARTHQUAKE *** #
        if t1+dt.timedelta(hours=1) >= eq or t0-dt.timedelta(hours=1) >= eq: break
        n  = fullnorth[(fulldate >= t0-dt.timedelta(hours=1)) & (fulldate < t1+dt.timedelta(hours=1))]
        e  = fulleast [(fulldate >= t0-dt.timedelta(hours=1)) & (fulldate < t1+dt.timedelta(hours=1))]
        u  = fullup   [(fulldate >= t0-dt.timedelta(hours=1)) & (fulldate < t1+dt.timedelta(hours=1))]
        # ------------------------------------------------------------------------------------------ #
        # Cross-correlation between the current day and the stack and maximize the cross-correlation #
        # ------------------------------------------------------------------------------------------ #
        nmax              = nstack.size
        ncc,ecc,ucc       = xcorr.correlate(nstack/ntrace[0],n,nmax),xcorr.correlate(estack/ntrace[1],e,nmax),xcorr.correlate(ustack/ntrace[2],u,nmax)
        nnoiz,enoiz,unoiz = np.std(ncc),np.std(ecc),np.std(ucc)
        nshift,n_val      = xcorr.xcorr_max(ncc,abs_max=False)
        eshift,e_val      = xcorr.xcorr_max(ecc,abs_max=False)
        ushift,u_val      = xcorr.xcorr_max(ucc,abs_max=False)
        # ------------------------------------- #
        # Ensure that the shifts are reasonable #
        # ------------------------------------- #
        if (nshift == 0) or (n_val < 4.0/3.0*nnoiz) or (np.abs(nshift) > max_shift): donot_north = True
        if (eshift == 0) or (e_val < 4.0/3.0*enoiz) or (np.abs(nshift) > max_shift): donot_east  = True
        if (ushift == 0) or (u_val < 4.0/3.0*enoiz) or (np.abs(nshift) > max_shift): donot_up    = True
        # ---------------------------- #
        # Extract the relevant portion #
        # ---------------------------- #
        nmax  = nstack.size
        n0    = int(3600.0/delta+0.5001)
        n,e,u = n[n0-nshift:n0-nshift+nmax],e[n0-eshift:n0-eshift+nmax],u[n0-eshift:n0-eshift+nmax]
        # -------------------------------- #
        # Add the current day to the stack #
        # -------------------------------- #
        if not donot_north: nstack += n - n[0] ; ntrace[0] += 1
        if not donot_east : estack += e - e[0] ; ntrace[1] += 1
        if not donot_up   : ustack += u - u[0] ; ntrace[2] += 1
    # -------------------- #
    # Normalize the stacks #
    # -------------------- #
    nstack,estack,ustack = nstack/ntrace[0],estack/ntrace[1],ustack/ntrace[2]
    ####################### SIDEREAL FILTER TEMPLATE ######################
    
    # ---------------------------------------- #
    # Safeguard if not enough traces were used #
    # ---------------------------------------- #
    napply,eapply,uapply = True,True,True
    if ntrace[0] <= 2: print "Not many traces have been used to create the sidereal filter (N)"; napply = False
    if ntrace[1] <= 2: print "Not many traces have been used to create the sidereal filter (E)"; eapply = False
    if ntrace[2] <= 2: print "Not many traces have been used to create the sidereal filter (U)"; uapply = False
    
    ####################### BUILD THE SIDEREAL FILTER #######################
    # ------------------------------ #
    # Initialize the sidereal filter #
    # ------------------------------ #
    tfilt = fulltime.copy()
    nfilt = np.zeros(fullnorth.size)
    efilt = np.zeros(fulleast .size)
    ufilt = np.zeros(fullup   .size)
    ndays = 0
    while True:
        # ---------------- #
        # Extract each day #
        # ---------------- #
        t0 = fulldate[0] + dt.timedelta(days=ndays+0)
        t1 = fulldate[0] + dt.timedelta(days=ndays+1)
        if (t1 > fulldate[-1]): t1 = fulldate[-1]
        n = fullnorth[(fulldate >= t0) & (fulldate < t1)]
        e = fulleast [(fulldate >= t0) & (fulldate < t1)]
        u = fullup   [(fulldate >= t0) & (fulldate < t1)]
        t = np.arange(n.size) * delta
        # ------------------------ #
        # Calculate a linear trend #
        # ------------------------ #
        nslope,nintercept,r_value,p_value,std_err = stats.linregress(t,n)
        eslope,eintercept,r_value,p_value,std_err = stats.linregress(t,e)
        uslope,uintercept,r_value,p_value,std_err = stats.linregress(t,u)
        # --------------------- #
        # Remove a linear trend #
        # --------------------- #
        n = n - (nslope*t + nintercept)
        e = e - (eslope*t + eintercept)
        u = u - (uslope*t + uintercept)
        # ----------------------------- #
        # Compute the cross-correlation #
        # ----------------------------- #
        nmax        = e.size
        ncc,ecc,ucc = xcorr.correlate(n,nstack[:nmax],nmax),xcorr.correlate(e,estack[:nmax],nmax),xcorr.correlate(u,ustack[:nmax],nmax)
        # ----------------------------------- #
        # Ensure that the shift is reasonable #
        # ----------------------------------- #
        ncc[:nmax-1/2-max_shift ] = 0.0
        ncc[ nmax-1/2+max_shift:] = 0.0
        ecc[:nmax-1/2-max_shift ] = 0.0
        ecc[ nmax-1/2+max_shift:] = 0.0
        ucc[:nmax-1/2-max_shift ] = 0.0
        ucc[ nmax-1/2+max_shift:] = 0.0
        # -------------- #
        # Find the shift #
        # -------------- #
        nshift,val = xcorr.xcorr_max(ncc,abs_max=False)
        eshift,val = xcorr.xcorr_max(ecc,abs_max=False)
        ushift,val = xcorr.xcorr_max(ucc,abs_max=False)
        # -------------------------- #
        # Create the sidereal filter #
        # -------------------------- #
        t          = np.arange(e.size)*delta + ndays*86400.0
        idx        = np.where(np.isin(tfilt,t))[0]
        nfilt[idx] = (np.roll(nstack,nshift) * signal.tukey(nstack.size,alpha=0.05))[:idx.size]
        efilt[idx] = (np.roll(estack,eshift) * signal.tukey(estack.size,alpha=0.05))[:idx.size]
        ufilt[idx] = (np.roll(ustack,ushift) * signal.tukey(ustack.size,alpha=0.05))[:idx.size]
        # ------------- #
        # Exit strategy #
        # ------------- #
        if t1 >= fulldate[-1]: break
        ndays = ndays + 1
    ####################### BUILD THE SIDEREAL FILTER #######################

    # -------------------------------- #
    # Remove gaps from the time series #
    # -------------------------------- #
    fulltime         = np.arange(time[0],time[-1]+delta,delta)
    fulldate         = np.asarray( [date[0] + dt.timedelta(seconds=i*delta) for i in range(fulltime.size)] )
    fullnorth        = np.ones(fulltime.size) * np.float('nan')
    fulleast         = np.ones(fulltime.size) * np.float('nan')
    fullup           = np.ones(fulltime.size) * np.float('nan')
    index            = np.where(np.isin(fulltime,time))[0]
    fullnorth[index] = north
    fulleast [index] = east
    fullup   [index] = up
    # -------------------------------------- #
    # Compute the RMS before sidereal filter #
    # -------------------------------------- #
    nrms0 = np.std(north[(date<eq)])
    erms0 = np.std(east [(date<eq)])
    urms0 = np.std(up   [(date<eq)])
    # ------------------------- #
    # Apply the sidereal filter #
    # ------------------------- #
    fullup    = fullup    - ufilt
    fulleast  = fulleast  - efilt
    fullnorth = fullnorth - nfilt
    idx       = np.where(np.isin(tfilt,time))[0]
    if uapply: up    = fullup   [idx]
    if eapply: east  = fulleast [idx]
    if napply: north = fullnorth[idx]
    # ------------------------------------- #
    # Compute the RMS after sidereal filter #
    # ------------------------------------- #
    nrms1 = np.std(north[(date<eq)])
    erms1 = np.std(east [(date<eq)])
    urms1 = np.std(up   [(date<eq)])
    # --------------------------------------- #
    # Print the effect of the sidereal filter #
    # --------------------------------------- #
    print 'STATION %s' % stnm
    print 'RMS before and after S.F. (north) = %7.2f/%7.2f (%7.2f)' % (nrms0,nrms1,(nrms0-nrms1)/nrms0*100.0)
    print 'RMS before and after S.F. (east ) = %7.2f/%7.2f (%7.2f)' % (erms0,erms1,(erms0-erms1)/erms0*100.0)
    print 'RMS before and after S.F. (up   ) = %7.2f/%7.2f (%7.2f)' % (urms0,urms1,(urms0-urms1)/urms0*100.0)
