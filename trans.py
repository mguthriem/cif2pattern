def SimTrans3_tof(fit_params, tof, y, e):
    '''
    %SimTrans calculates transmission spectrum from two crystals
    %   lam - array containing wavelengths to calc over
    %   hkl - contains all Nref hkl's that calculation is performed for
    %   bgd - array containing coefficients of polynomial for background
    %   sf - overall scale factor
    %   pktype - 1 = gauss; 2 = lorentz; ...
    %   UB1 - UB matrix for first crystal
    %   setang1 - setting angles for first crystal (deviations from ideal UB
    %               location).
    %   pkpars1 - position(lambda), position(d-spacing), position(ttheta), width, intensity for each Nref reflections
    %   UB2,setang2,pkpars2 - as above, for second crystal
    %
    % M. Guthrie 21st Jan 2014
    %
    % calculate background profile
    % determine number of coeffs bgd
    %
    % version 2 constrains intensities for equivalent dips to be the same
    % M. Guthrie 3 April 2014
    %
    % M. Guthrie 7 April 2014, realised I was making an (obvious) mistake by
    % adding the transmissions from the two diamonds. Clearly, they should be
    % multiplied. I've implemented the change...will see what difference it
    % makes.
    %
    % M. Guthrie 9 April 2014, introduced possibility to refine L2 and also a
    % scale factor for calculated dip wavelengths (to account for diamond
    % compressibility).
    %
    %M. Guthrie 12 March 2020: (version a) want to upgrade simulated spectrum with developments since
    % 2014. This will include: 
    % 1) including the stress tensor to allow peak positions to shift asymmetrically according to (the massive) diamond strain 
    % 2) allowing more background coefficiets for more complex baxkground shape
    % 3) attempt to allow refinement of "difa", currently passed as via a global
    % In doing this, it makes sense to re-number all of the parameters, as it's extremely messy at the moment
    % The parameters that depend on the number of observed reflections should be moved to the end
    % of the array so that the numbering of other static parameters doesn't change with every successive
    % edit of the software. I will just leave sufficient space to allow for future tweaking...
    %
    % M. Guthrie 16 March 2020 on suggestion of Thomas Holm Rod, have switched to using lmfit
    % this looks much better for setting up as parameters are specified by a library
    % and not some random location in an array. Also, setting up of constraints and 
    % boundary conditions looks much more intuitive.
    %
    %
    % M. Guthrie 19 March 2020modified fitting so that it occurs in TOF space, this is more correct as
    % that is the measurement parameter. To allow easy flipping of lambda and TOF
    % I defined the functions lam2tof and tof2lam. 
    % note in previous versions, I used the term difa in analogy with the gsas term
    % however, I realise the resulting expression lam = A*T + difa*T^2 is confusing as
    % gsas operates with the inverse of this function (and with d instead of lambda)
    % long story short: I'm keeping difa, but it is now a different quantity...
    % 
    % M. Guthrie 22 Dec 2020 (version a): added a psuedovoigt peak shape adjusted gaussian
    % peak definitions for consistency.
    %
    % M. Guthrie 23 Dec 2020: added tick marks for calculated peak positions
    %
    % M. Guthrie 07 Jan 2021: fixed bug in pkmult setting during fitting
    % adjusted use of pktype to allow pure gauss, pure lorentz or psuedovoigt. Found another bug in peak definitions where
    % pkmult1 was used as multiplier for diamond 2 dips.
    % Added file output containing dip: hkl, lambda
    %
    % M. Guthrie 26 Jan 2021: added Bk2BkExpConvPV peak type. this was quite a bit of hacking as I had to create additional peak par terms
    % (I ended up increasing to 10 although I only need 6, for future proofing). Also, am using a bit of a clunky way to evaluate the 
    % peaks (using a mantid ws), currently the only way I can figure out to do it. 
    %
    % M Guthrie 5 Feb 2021: major change to peak profile definitions, now using GSAS TOF profile 3 as reference
    % increased number of possible parameter to 31 and renamed everything for more consistency.
    '''
    
    
    
    global hkl1, hkl2
    global UB1, pkcalcint1
    global UB2, pkcalcint2
    global pktype,nbgd,bgdtype  #pktype = 1 is Guass, = 2 is Lorenzt, 3 = psuedovoigt and mix par refined.
    global lam, shftlam
    global L1
    global ttot,t1,t2,bgdprof
    global fxsamediam
    global neqv1, eqvlab1, neqv2, eqvlab2
    global function_verbose
    global runNumber,directory,insTag
    
    ioioName = insTag + str(runNumber) + '_monitors_ioio' #this workspace is used to provide instrument info for Bk2BkExpConvPV peak shape

    nref1 = hkl1.shape[0]  # % number of reflections to integrate over
    nref2 = hkl2.shape[0]  # % number of reflections to integrate over
    # % returns array with same dim as input labelling equivs
    eqvlab1, neqv1 = findeqvs(hkl1)
    eqvlab2, neqv2 = findeqvs(hkl2)

    bgd = np.zeros(24) # up to 24 term chebychev background used 
    pkmult1 = np.zeros(neqv1)
    pkmult2 = np.zeros(neqv2)  
    sf = fit_params['sf'].value                 
    relsf =fit_params['relsf'].value            
    setang1 = np.zeros(3)
    setang1[0] = fit_params['setang1_alp'].value        #setting angles for diamond 1
    setang1[1] = fit_params['setang1_bet'].value
    setang1[2] = fit_params['setang1_gam'].value
    setang2 = np.zeros(3)
    setang2[0] = fit_params['setang2_alp'].value        #setting angles for diamond 2
    setang2[1] = fit_params['setang2_bet'].value
    setang2[2] = fit_params['setang2_gam'].value
    epsilon = np.zeros(shape=(3,3),dtype=float)
    epsilon[0,0] = fit_params['Estr11'].value
    epsilon[1,1] = fit_params['Estr22'].value
    epsilon[2,2] = fit_params['Estr33'].value
    
    epsilon[0,1] = fit_params['Estr12'].value
    epsilon[0,2] = fit_params['Estr13'].value
    epsilon[1,2] = fit_params['Estr23'].value

    alp = fit_params['pkCoef00'].value   
    bet0 = fit_params['pkCoef01'].value 
    bet1 = fit_params['pkCoef02'].value  
    sig0 = fit_params['pkCoef03'].value
    sig1 = fit_params['pkCoef04'].value
    sig2 = fit_params['pkCoef05'].value 
    gam0 = fit_params['pkCoef06'].value
    gam1 = fit_params['pkCoef07'].value
    gam2 = fit_params['pkCoef08'].value
    eta = fit_params['pkCoef21'].value
    regwid1 = fit_params['pkCoef26'].value
    regwid2 = fit_params['pkCoef27'].value
    
    L2 = fit_params['L2'].value              #secondary flight path (to monitor)
    difa = fit_params['difa'].value            #quadratic TOF to lambda term
    for i in range(24):
        bgd[i] = fit_params['bgd'+str(i+1).zfill(3)]
    for i in range(neqv1):
        pkmult1[i] = fit_params['pkmult1_'+str(i+1).zfill(3)]
    for i in range(neqv2):
        pkmult2[i] = fit_params['pkmult2_'+str(i+1).zfill(3)]
    if fxsamediam == 1:
        if neqv1==neqv2:
            pkmult2 = pkmult1*relsf
        elif neqv2>neqv1:
            pkmult2[0:neqv1] = pkmult1[0:neqv1]*relsf
        elif neqv1>neqv2:
            pkmult2[0:neqv2]=pkmult1[0:neqv2]*relsf

    # number of data points to calculate over fit over
    npt = tof.shape[0]
    
    # calculate information for peaks for crystal 1 using hkl,UB1, setang,
    # pkpos    
    a, b, c, q = pkposcalc_strain(hkl1, UB1, setang1, epsilon)
    diamTOF1 = lam2tof(a,(L1+L2),difa) #get time of flight for each hkl
    lam1 = a
    srtlam1 = np.sort(lam1)[::-1] #array of lambdas sorted with lambda descending
    isrtlam1 = np.argsort(lam1)[::-1] #corresponding sorted indices
    srthkl1 = hkl1[isrtlam1] #and sorted hkl
    srtd1 = b[isrtlam1] #sorted d-spacing
    srtQ1 = q[isrtlam1,:] #sorted Q coords
    srtTT1 = c[isrtlam1] #sorted TTheta
    srtdiamTOF1 = diamTOF1[isrtlam1] #sorted TOFs for diamond 1
    srteqvlab1 = eqvlab1[isrtlam1] #sorted equivalent lables
    # calculate information for peaks for crystal 2 using hkl,UB1, setang,
    # pkpos
    a, b, c, q = pkposcalc_strain(hkl2, UB2, setang2, epsilon)
    diamTOF2 = lam2tof(a,L1+L2,difa)
    lam2 = a
    srtlam2 = np.sort(lam2)[::-1] #useful to have file output sorted with lambda descending
    isrtlam2 = np.argsort(lam2)[::-1] #corresponding sorted indices
    srthkl2 = hkl2[isrtlam2] #and sorted hkl       
    srtd2 = b[isrtlam2] #sorted d-spacing
    srtQ2 = q[isrtlam2,:] #sorted Q coords
    srtTT2 = c[isrtlam2]
    srtdiamTOF2 = diamTOF2[isrtlam2] #sorted TOFs for diamond 2
    srteqvlab2 = eqvlab2[isrtlam2] #sorted equivalent lables
            
    yvals = np.divide(lam1,lam1)*1.03;        
    ticks1 = CreateWorkspace(OutputWorkspace='ticks1',DataX=lam1, DataY=yvals, NSpec=1,UnitX='Wavelength')        

    yvals = np.divide(lam2,lam2)*1.03;        
    ticks2 = CreateWorkspace(OutputWorkspace='ticks2',DataX=lam2, DataY=yvals, NSpec=1,UnitX='Wavelength')        
    
    # used to add custom ticks if necessary
    #mystery = np.array([2.245,2.274,2.318,2.500])
    #yvals = np.divide(mystery,mystery)*1.04;
    #ticks3 = CreateWorkspace(OutputWorkspace='ticks3',DataX=mystery, DataY=yvals, NSpec=1,UnitX='Wavelength')        

    if bgdtype ==1: #Chebychev background
#        print('Chebychev background being fitted')
        print('nBgd= ',nbgd)
        print('BGD: ',bgd[0:nbgd])
        bgdprof = np.polynomial.chebyshev.chebval(tof,bgd[0:nbgd])
        
    tofsclfact = 2.00e-3 

#####################################################################################
# calculate dip spectrum here and write peak info to file
#####################################################################################

    peakInfoFname1 = directory + insTag + str(runNumber) + '_TransDipInfo1.log'    
    peakInfoFname2 = directory + insTag + str(runNumber) + '_TransDipInfo2.log'  
    
    f1 = open(peakInfoFname1,'w')
    f2 = open(peakInfoFname2,'w')
    
    regionTOFmin = 9620.0
    regionTOFmax = 10970.0
    regionwidscale = np.array([regwid1,regwid2])# first affect diam 1, second diam 2
    regionintscale = np.array([1.0,1.0])#

    # calculate peaks for crystal 1

    peakInfoFname1 = directory + insTag + str(runNumber) + '_TransDipInfo1.log'    
    peakInfoFname2 = directory + insTag + str(runNumber) + '_TransDipInfo2.log'  
    
    f1 = open(peakInfoFname1,'w')
    f1.write('   h   k   l     lam      Qx      Qy      Qz  ttheta     clk    mult eqv\n')
    f2 = open(peakInfoFname2,'w')
    f2.write('   h   k   l     lam      Qx      Qy      Qz  ttheta     clk    mult eqv\n')

    # calculate peaks for diamond 1
    dips1 = np.zeros(npt)  # initialise array containing dip profile
    for i in range(nref1):

        TOF = srtdiamTOF1[i]
        scl = TOF*tofsclfact
        Sigma2 = sig0**2+sig1**2*TOF**2+sig2**2*TOF**4  #gaussian variance
        if TOF>regionTOFmin and TOF<regionTOFmax:
            Sigma2 = Sigma2*regionwidscale[0]
            scl = scl*regionintscale[0]

        Alpha = alp/TOF
        Beta = bet0 + bet1/TOF**4
        Gamma = gam0+gam1*TOF+gam2*TOF**2
        FWHM = 2*np.sqrt(2*np.log(2))*np.sqrt(Sigma2)
       

        Amp = 1e4*scl*pkmult1[int(srteqvlab1[i])]
        
        if pktype == 1: # Gauss
            func = Gaussian()
            Height = Amp/(np.sqrt(2*np.pi*Sigma2))
            pin = [Height,TOF,np.sqrt(Sigma2)]
        elif pktype == 2: # Lorentz
            func = Lorentzian()
            pin = [Amp,TOF,FWHM]
        elif pktype == 3: # psuedovoigt
            func = PseudoVoigt()
            pin = [eta,Amp,TOF,FWHM]
        elif pktype == 4: # Bk2BkExpConvPV            
            func = Bk2BkExpConvPV()#name = Bk2BkExpConvPV,X0=-0,Intensity=0,Alpha=1,Beta=1,Sigma2=1,Gamma=0
            pin = [TOF,Amp,Alpha,Beta,Sigma2,Gamma]

        ws = mtd['tmp_peak_ws']
        [func.function.setParameter(ii,p) for ii,p in enumerate(pin)]
        ws_eval = func(ws)
        testyin = ws_eval.readY(1)     
        dips1= dips1 - testyin
        

        if srtQ1[i,0]<=0: #in mantid this mean rhs of vertical centre, looking along beam
            clock = 360.0- np.degrees( np.arccos(srtQ1[i,1]/np.linalg.norm(srtQ1[i,:])))
        else: #in mantid this mean lhs of vertical centre, looking along beam
            clock = np.degrees( np.arccos(srtQ1[i,1]/np.linalg.norm(srtQ1[i,:]))) 
        f1.write('%4u%4u%4u%8.4f%8.4f%8.4f%8.4f%8.1f%8.1f%8.4f%4u\n'%(srthkl1[i,0],srthkl1[i,1],\
        srthkl1[i,2],srtlam1[i],srtQ1[i,0],srtQ1[i,1],srtQ1[i,2],srtTT1[i],clock,\
        pkmult1[int(srteqvlab1[i])],int(srteqvlab1[i])+1))
            
    # calculate peaks for diamond 2
    dips2 = np.zeros(npt)  # initialise array containing profile
    for i in range(nref2):
        TOF = srtdiamTOF2[i]
        scl = TOF*tofsclfact
        Sigma2 = sig0**2+sig1**2*TOF**2+sig2**2*TOF**4  #gaussian variance
        if TOF>regionTOFmin and TOF<regionTOFmax:
            Sigma2 = Sigma2*regionwidscale[1]
            scl = scl*regionintscale[1]
        Alpha = alp/TOF
        Beta = bet0 + bet1/TOF**4
        Gamma = gam0+gam1*TOF+gam2*TOF**2
        FWHM = 2*np.sqrt(2*np.log(2))*np.sqrt(Sigma2)
       

        Amp = 1e4*scl*pkmult2[int(srteqvlab2[i])]
        
        if pktype == 1: # Gauss
            func = Gaussian()
            Height = Amp/(np.sqrt(2*np.pi*Sigma2))
            pin = [Height,TOF,np.sqrt(Sigma2)]     
        elif pktype == 2: # Lorentz
            func = Lorentzian()
            pin = [Amp,TOF,FWHM]
        elif pktype == 3: # psuedovoigt
            func = PseudoVoigt()
            pin = [eta,Amp,TOF,FWHM]
        elif pktype == 4: # Bk2BkExpConvPV            
            func = Bk2BkExpConvPV()#name = Bk2BkExpConvPV,X0=-0,Intensity=0,Alpha=1,Beta=1,Sigma2=1,Gamma=0
            pin = [TOF,Amp,Alpha,Beta,Sigma2,Gamma]

        ws = mtd['tmp_peak_ws']
        [func.function.setParameter(ii,p) for ii,p in enumerate(pin)]
        ws_eval = func(ws)
        testyin = ws_eval.readY(1)     
        dips2= dips2 - testyin
                
        if srtQ2[i,0]<=0: #in mantid this mean rhs of vertical centre, looking along beam
            clock = 360.0 - np.degrees( np.arccos(srtQ2[i,1]/np.linalg.norm(srtQ2[i,:])))
        else: #in mantid this mean lhs of vertical centre, looking along beam
            clock = np.degrees( np.arccos(srtQ2[i,1]/np.linalg.norm(srtQ2[i,:]))) 
        f2.write('%4u%4u%4u%8.4f%8.4f%8.4f%8.4f%8.1f%8.1f%8.4f%4u\n'%(srthkl2[i,0],srthkl2[i,1],\
        srthkl2[i,2],srtlam2[i],srtQ2[i,0],srtQ2[i,1],srtQ2[i,2],srtTT2[i],clock,\
        pkmult2[int(srteqvlab2[i])],int(srteqvlab2[i])+1))
  
        
    f1.close()
    f2.close()

    # calculate final profile
    t1 = 1+sf*dips1 #individual transmission 
    t2 = 1+sf*dips2 #individual transmission
    ttot = t2*t1*bgdprof #note order of multiplication doesn't have an affect (i.e. either diamond could be upstream of the other)
    # explanation:
    #Beam from left to right
    # I_0 ------> [diam1]----->I_1----->[diam2]-----I_2 [d/s monitor]
    #
    # so I_1 = t_1*I_0, and
    # I_2 = t_2*I_1 = t_2*t_1*I_0
    #
    # following Guthrie et al J Appl Cryst.( 2017) 50 https://doi.org/10.1107/S1600576716018185
    # I_0 taken to equal bgdprof and t_1 and t_2 = 1.0 where there is no dip.
    
    if y.size == 0:
        return ttot
    if e.size == 0:
        return ttot - y
    resid = np.divide((ttot-y),e)
    return resid #lmfit requies array to be returned