! Thomas Meissner 
! 01/12/2026
! combines C-W band and L-band surface emission RTM
! log(frequeny) interpolation between L-C (VH) and L-X (S3/S4)
!
! References
!
! [MW2004]       Meissner, T. and F. Wentz, The complex dielectric constant of pure and sea water from microwave satellite observations, 
!                IEEE TGRS, 2004, 42(9), 1836 – 1849, doi:10.1109/TGRS.2004.831888.
!
! [MW2012]	     Meissner, T. and F. Wentz, The emissivity of the ocean surface between 6 and 90 GHz 
!                over a large range of wind speeds and Earth incidence angles, 
!                IEEE TGRS, 2012, 50(8), 3004 – 3026, doi: 10.1109/TGRS.2011.2179662.   
!
! [MWR2014]  	 Meissner, T., F. Wentz, F. and L. Ricciardulli, The emission and scattering of L-band microwave radiation 
!                from rough ocean surfaces and wind speed measurements from Aquarius, 
!                J. Geophys. Res. Oceans, 2014, 119, doi:10.1002/2014JC009837.
!
!

module RSS_SURFACE_EMISSION_RTM_ALLFREQ_MODULE
save 
public

!public :: find_surface_tb_allfreq, dielectric_meissner_wentz, fdem0_meissner_wentz, fd_emiss_allfreq, fd_scatterm_all, fd_tcos_eff  


! external files     

! [MW 2012] Table 2
character(len=200), parameter     :: file_coeffs_wind_isotropic = '\\THOMAS\C\RSS_RTM\external_files\finetune_emiss_wind.dat' 

! [MW 2012] Tables 3 + 4
character(len=200), parameter     :: file_coeffs_wind_direction = '\\THOMAS\C\RSS_RTM\external_files\fit_emiss_phir_wind.dat' 

! [MW 2012] Tables 8, 9, 10, 11, 12
character(len=200), parameter     :: file_coeffs_sctterm        = '\\THOMAS\C\RSS_RTM\external_files\mk_scatterm_table_all.dat' 

! [MW 2012] eq. (15) for fast computation
character(len=200), parameter     :: file_em0_ref_freq_sst      = '\\THOMAS\C\RSS_RTM\external_files\em0_ref_freq_sst.dat' 


character(len=200), parameter     :: emiss_coeff_harm_file_AQ   = '\\Thomas\C\AQUARIUS\surface_roughness\tables\V_may_2013\coeffs\deW_harm_coeffs_V9A_MI.dat'

character(len=200), parameter	  :: emiss_coeff_harm_file_SMAP = '\\THOMAS\C\SMAP\SSS_algo_V3\roughness_model\tables\dew_phi_VH34_harmonic_tab.dat'


! Aquarius GMF parameters
integer(4), parameter                       :: n_rad=3
real(4), parameter                          :: sst_ref0=20.0, sss_ref0=35.0, freq_aq=1.413, wcasp=11.0, wextrapol=17.0, teff0=290.
real(4), parameter, dimension(n_rad)        :: tht_ref_AQ = (/29.36, 38.44, 46.29/)
real(4), parameter                          :: tht_ref_SMAP=40.0


contains

!
! Routines:

! 1.  find_surface_tb_ALLFREQ        Master routine. Caclulates all components of the ocean surface emissivity. Output user selected.

! 2.  dielectric_meissner_wentz      Dielectric model of sea and pure water [MW 2004], with minor updates in [MW 2012].

! 3.  fdem0_meissner_wentz           Caclulates emissivity of specular surface [MW 2004]

! 4.  fd_emiss                       Calculates emissivity of specular surface (v/h), isotropic wind induced emissivity (v/h) and wind direction signal (v/h/S3/S4) 
!                                    [MW 2012], sections IV + VI.
!
! 5.  fd_scatterm_all                Calculates correction for downwelling scattered atmospheric radiation [MW 2012], section V.
!
! 6.  fd_tcos_eff                    Calculates effective cold space temperature taking into account the deviation between Rayleigh-Jeans and Planck law
!                                    as function of frequency [MWD 2011], Appendix D.
!
! 7.  get_emiss_wind                 Calculates isotropic wind induced emissivity at tabulated frequency values and reference EIA of 55.2. 
!                                    [MW 2012], section IV, Table 2.
!                               
! 8.  get_aharm_phir                 Calulates harmonic ocefficients of wind direction signal at tabulated frequency values and reference EIA of 55.2.
!                                    [MW 1012], section VI, Tables 3 + 4.
!
! 9.  get_aharm_phir_nad             Calculates harmonic ocefficients of wind direction signal at tabulated frequency values at nadir [MW 2012, section VI, eq. (26)]
!
! 10. get_sst_fac                    Precomputed value of E0(SST)/E0(SST_ref=20C) for faster computation [MW 2012], eq. (15) 
!
! 11. fd_xmea_win                    Provides wind speed polynomials for [MW 2012], eqs. (14) + (25) 



! Master interface routine
subroutine find_surface_tb_ALLFREQ ( freq,tht,surtep,sal,ssws,phir,tran,tbdw,tc,        &
                                     e0,ewind,omega,edirstokes,eharm,tbscat,tbsurf)

!   input 
!   name   descritpion                                            type      dimension     unit         valid physical range   Default value              
!   freq   frequency                                              real(4)   scalar        GHz          [1.0, 100.0]     

!   tht    earth incidence angle                                  real(4)   scalar        deg          [0,65]

!   surtep sea surface temperature                                real(4)   scalar        Kelvin       [271.15,313.15]                                   

!   sal    sea surface salinity                                   real(4)   scalar        ppt          [0,45]                optional    35.0 

!   ssws   sea surface wind speed                                 real(4)   scalar        m/s          [0,100]               optional                                           

!   phir   sea surface wind direction relative to azimuthal look  real(4)   scalar        deg          [0,360[               optional
!   upwind=0, downwind=180 deg

!   tran   atmospheric trnasmittance                              real(4)   scalar                     [0,1]                 optional

!   tbdw   downwelling atmospheric brightness temperature         real(4)   scalar        Kelvin       >=0                   optional

!   tc     cold space temperature                                 real(4)   scalar        Kelvin       >=0                   optional   value from fd_tcos_eff  
!
!   output
!   name    descritpion                                           type      dimension     unit         physical range              

!   e0      specular sea surface emissivity  (1)=v, (2)=h         real(4)   vector(2)                  [0,1]                 optional 

!   ewind   isotropic wind induced sea surface emissivity         real(4)   vector(2)                  [0,1]                 optional
!           (1)=v, (2)=h


!   omega   omega term (correction for scattered downwelling      real(4)   vector(2)                  [0,1]                 optional
!           reflection to include reflection from 
!           non specular directions)
!           calculated from geometric optics model using 
!           updated slope distribution
!           see amsr atbd  



!   edirstokes wind directional emissivity signal                 real(4)   vector(4)                  [0,1]                 optional
!           modified stokes vector
!           (1)=v, (2)=h, (3)=s3, (4)=s4         


!   eharm harmonic coefficients of wind direction signal          real(4)   array(2,4)                 [0,1]                 optional  
!           modified stokes vector
!           index 1: 1=1st 2=2nd harmonic                          
!           index 2: (1)=v, (2)=h, (3)=s3, (4)=s4         
             

!   tbscat  scattered downweling radiation                        real(4)   vector(4)     Kelvin       >=0                   optional
!           modified stokes vector
!           (1)=v, (2)=h, (3)=s3, (4)=s4  

!   tbsurf  total brightness temperature from sea surface         real(4)   vector(4)     Kelvin       >=0                   optional        
!           emitted and reflected
!           modified stokes vector
!           (1)=v, (2)=h, (3)=s3, (4)=s4  


!   the following polarization basis vector convention for modified stokes parameters is used:
!   h = (k cross n) / abs(k cross n) (k: propagation direction from E to S/S, 
!   n: local mean sea surface normal upward from earth)                                        
!   v = h cross k
!   p = (v +  h)/sqrt(2)
!   m = (v -  h)/sqrt(2)
!   r = (v + ih)/sqrt(2) (ieee convention)
!   l = (v - ih)/sqrt(2) (ieee convention)
!   s3 = p - m  
!   s4 = l - r

implicit none
save

    real(4), intent(in)                               ::    freq
    real(4), intent(in)                               ::    tht
    real(4), intent(in)                               ::    surtep
    real(4), intent(in), optional                     ::    sal
    real(4), intent(in), optional                     ::    ssws
    real(4), intent(in), optional                     ::    phir
    real(4), intent(in), optional                     ::    tran
    real(4), intent(in), optional                     ::    tbdw
    real(4), intent(in), optional                     ::    tc

    real(4), intent(out), dimension(2),   optional	  ::    e0 
    real(4), intent(out), dimension(2),   optional	  ::    ewind 
    real(4), intent(out), dimension(2),   optional	  ::    omega 
    real(4), intent(out), dimension(4),   optional	  ::    edirstokes
    real(4), intent(out), dimension(2,4), optional	  ::    eharm
    real(4), intent(out), dimension(4),   optional	  ::    tbscat 
    real(4), intent(out), dimension(4),   optional	  ::    tbsurf 


    real(4)			                                  ::    xtc, xomegabar, sst, xsal, xssws, xphir
    real(4), dimension(2)                             ::    xe0, xewind, xscat, xomega
    real(4)				                              ::    costht, path, opacty
    real(4), dimension(4)                             ::    xestokes, xtbscat, xetot, xrtot
    real(4), dimension(2,4)                           ::    xeharm
    real(4)                                           ::    xfreq
    
    if (freq<1 .or. freq>100.) then
    write(*,*) freq,' freq OOB in RSS_SURFACE_EMISSION_RTM_ALLFREQ. 1<=freq<=100  required'
    stop
    endif
    
    xfreq=freq     
    sst=surtep-273.15

    if (present(ssws)) then
      xssws=ssws
    else
      xssws=0.0
    endif

    if (present(phir)) then
      xphir=phir
    else
      xphir=-999.0
    endif


    if (present(sal)) then
      xsal=sal
    else
      xsal=35.0
    endif 

    ! effective cold space temperature including deviation from Rayleigh-Jeans law
    if (present(tc)) then
      xtc=tc
    else
      call fd_tcos_eff(freq, xtc)
    endif


    ! bug fix TM 06/29/2026: include sss 
    call fd_emiss_allfreq(freq=xfreq,tht=tht,sst=sst,sal=xsal,wind=xssws,phir=xphir, &
                  emiss_0=xe0, emiss_wind=xewind, emiss_phi=xestokes, eharm=xeharm)


    if (present(e0))           e0=xe0
    if (present(ewind))        ewind=xewind
    if (present(edirstokes))   edirstokes=xestokes
    
    if (present(eharm))        eharm=xeharm

    if (present(tbscat) .or. present(tbsurf) .or. present(omega)) then

        if ( .not.(present(tran)))   stop  ' need tran for computing tbscat. pgm stopped.'
        if ( .not.(present(tbdw)))   stop  ' need tbdw for computing tbscat. pgm stopped.'

        costht=cosd(tht)
        path=1.00035/sqrt(costht*costht+7.001225e-4)   
        ! (1+hratio)/sqrt(costht**2+hratio*(2+hratio)), hratio=.00035 [MW 2012] eq. (3)
        
        opacty=-alog(tran)/path
        call fd_scatterm_all(freq,tht,xssws,opacty, xscat)
        xomega(1:2) = xscat(1:2) /(tbdw + tran*xtc - xtc) ! [MW 2012] eq. (21)

        xetot(1:4) = (/xe0(1),xe0(2),0.0,0.0/) + (/xewind(1),xewind(2),0.0,0.0/) + xestokes(1:4)
        ! limit emissivity to [0,1]
        where(xetot(1:2)>=1.0) xetot(1:2)=1.0
        where(xetot(1:2)<=0.0) xetot(1:2)=0.0
      
        xrtot(1:2) = 1.0-xetot(1:2)
        xrtot(3:4) =    -xetot(3:4)
        xomegabar= (xomega(1)*xrtot(1) + xomega(2)*xrtot(2))/(xrtot(1)+xrtot(2))
        xtbscat(1:2) = ((1.0+xomega(1:2))*(tbdw+tran*xtc) - xomega(1:2)*xtc) * xrtot(1:2) 
        xtbscat(3:4) = ((1.0+xomegabar ) *(tbdw+tran*xtc) - xomegabar  *xtc) * xrtot(3:4)   
      
        if (present(tbscat)) tbscat = xtbscat
        if (present(tbsurf)) tbsurf = xetot*surtep + xtbscat
        if (present(omega))  omega=xomega

    endif
    
return
end subroutine find_surface_tb_ALLFREQ 


subroutine fd_emiss_ALLFREQ(freq,tht,sst,wind,phir,sal,     emiss_0,emiss_wind,emiss_phi,eharm)
!  Calculates emissivity of specular surface (v/h), isotropic wind induced emissivity (v/h) and wind direction signal (v/h/S3/S4) 
!  [MW 2012, sections IV + VI.]. [MWR 2014]
!  1 GHz <= freq <= 100 GHz
!  Frequency interpolation

implicit none

     real(4), intent(in)                    :: freq,tht,sst
     real(4), intent(in), optional          :: wind
     real(4), optional, intent(in)          :: sal
     real(4), optional, intent(in)          :: phir

     
     real(4), optional, intent(out)         :: emiss_0(2),emiss_wind(2),emiss_phi(4),eharm(2,4)

     real(4)                                :: em0(2), xemiss_wind(2), xwind, xsal, xphir
     real(4)                                :: xeharm(2,4), xeharm1(2,4), xeharm2(2,4)
     real(4)                                :: yemiss_wind(2), yeharm(2,4) 
     real(4)                                :: zemiss_wind(2), zeharm(2,4), zemiss_phi(4)  


     real(4), parameter                     :: freq0=freq_aq, freq1=6.8, freq2=10.7  ! interpolation points
     
     real(4)                                :: xfreq, xlog, xlog0, xlog1, xlog2, brief, qtht
 
     if (freq<1 .or. freq>100.) then
     write(*,*) freq,' freq OOB in RSS_SURFACE_EMISSION_RTM_ALLFREQ. 1<=freq<=100  required'
     stop
     endif
     
     xfreq=freq
          
     if (present(wind)) then
        xwind=wind
     else
        xwind=0.0
     endif

     if (present(phir)) then
        xphir=phir
     else
        xphir=-999.0
     endif

     if (present(sal)) then
        xsal=sal
     else
        xsal=35.0
     endif  

     qtht=tht
	 if(qtht.gt.65) qtht=65.  !qtht is just used for extrpolation for tht>thtref and i limit it to 60 deg

     ! emiss_0 (MW) is valid for all frequencies
     if (present(emiss_0)) then
        call fdem0_meissner_wentz(freq=freq,tht=tht,sst=sst,salinity=xsal, em0=em0) 
        emiss_0=em0
     endif

     if (xfreq<freq_aq) xfreq=freq_aq

     if (xfreq>freq2) then
     ! X band and higher
     call fd_emiss(freq=xfreq,tht=tht,sst=sst,wind=xwind,phir=xphir,sal=xsal,     emiss_wind=zemiss_wind,eharm=zeharm)     
     endif
     
     if (xfreq<=freq1) then !interpolate VH between L and C band
     ! C band 
     call fd_emiss(freq=freq1,tht=tht,sst=sst,wind=xwind,phir=xphir,sal=xsal,     emiss_wind=xemiss_wind,eharm=xeharm1)
     xeharm(1:2,1) = xeharm1(1:2,1) !V  @ C-band
     xeharm(1:2,2) = xeharm1(1:2,2) !H  @ C-band      
     
     ! L band
     call fd_emiss_LBAND (wspd=xwind, sst=sst, tht=tht,    emiss_wind=yemiss_wind, eharm=yeharm)
        
     ! interpolation between L and C for V and H 
     xlog0 = alog(freq0)
     xlog1 = alog(freq1)
     xlog  = alog(xfreq)
     
     brief = (xlog-xlog0)/(xlog1-xlog0)
     zemiss_wind   = yemiss_wind*(1.0-brief)  + xemiss_wind*brief
     zeharm(1:2,1) = yeharm(1:2,1)*(1.0-brief)+ xeharm(1:2,1)*brief !V
     zeharm(1:2,2) = yeharm(1:2,2)*(1.0-brief)+ xeharm(1:2,2)*brief !H
     endif  
     

     if (xfreq<=freq2) then  !  interpolate S3/S4 between L and X band
     ! X band
     call fd_emiss(freq=freq2,tht=tht,sst=sst,wind=xwind,phir=xphir,sal=xsal,                            eharm=xeharm2)
     xeharm(1:2,3) = xeharm2(1:2,3) !S3 @ X-band
     xeharm(1:2,4) = xeharm2(1:2,4) !S4 @ X-band
     
     ! L band
     call fd_emiss_LBAND (wspd=xwind, sst=sst, tht=tht,    emiss_wind=yemiss_wind, eharm=yeharm)
     
     ! interpolation between L and X for S3 and S4. 
     ! no interpolation for V and H.
     xlog0 = alog(freq0)
     xlog2 = alog(freq2)
     xlog  = alog(xfreq)
     
     brief = (xlog-xlog0)/(xlog2-xlog0)    
     zeharm(1:2,3) = yeharm(1:2,3)*(1.0-brief)+ xeharm(1:2,3)*brief !S3
     zeharm(1:2,4) = yeharm(1:2,4)*(1.0-brief)+ xeharm(1:2,4)*brief !S4
     endif 
     
     
     if (xfreq>=freq1 .and. xfreq<=freq2) then  !V/H
     call fd_emiss(freq=xfreq,tht=tht,sst=sst,wind=xwind,phir=xphir,sal=xsal,     emiss_wind=xemiss_wind,eharm=xeharm)     
     zemiss_wind = xemiss_wind
     zeharm(1:2,1:2) = xeharm(1:2,1:2)
     endif     
     
     if (xphir>-998.) then
     zemiss_phi(1) = zeharm(1,1)*cosd(xphir) + zeharm(2,1)*cosd(2.0*xphir)
     zemiss_phi(2) = zeharm(1,2)*cosd(xphir) + zeharm(2,2)*cosd(2.0*xphir)
     zemiss_phi(3) = zeharm(1,3)*sind(xphir) + zeharm(2,3)*sind(2.0*xphir)
     zemiss_phi(4) = zeharm(1,4)*sind(xphir) + zeharm(2,4)*sind(2.0*xphir)
     else
     zemiss_phi=0.0
     endif
     
     if (present(emiss_wind)) emiss_wind=zemiss_wind
     if (present(emiss_phi))  emiss_phi =zemiss_phi
     if (present(eharm))      eharm     =zeharm     
     

return
end subroutine fd_emiss_ALLFREQ

 
subroutine fd_emiss(freq,tht,sst,wind,phir,sal,     emiss_0,emiss_wind,emiss_phi,emiss_tot,eharm)
!  Calculates emissivity of specular surface (v/h), isotropic wind induced emissivity (v/h) and wind direction signal (v/h/S3/S4) 
!  Frequency >= 6.8 GHz
!  [MW 2012, sections IV + VI.]
implicit none

     real(4), intent(in)                    :: freq,tht,sst
     real(4), intent(in), optional          :: wind
     real(4), optional, intent(in)          :: sal
     real(4), optional, intent(in)          :: phir

     
     real(4), optional, intent(out)         :: emiss_0(2),emiss_tot(2),emiss_wind(2),emiss_phi(4),eharm(2,4)

     real(4)                                :: xemiss_tot(2),xemiss(2),xemiss_phi(4),xsal 

     integer(4), parameter                  :: nstoke=4


     integer(4)                             :: ifreq1,ifreq2,ipol,istoke,iharm
     real(4)                                :: em0(2)
     real(4)                                :: qtht,wt,emiss1(2),emiss2(2),enad,h1,h2
     real(4)                                :: xphir
     real(4), save                          :: cos1phi,cos2phi,sin1phi,sin2phi
     real(4)                                :: aharm1(2,nstoke),aharm2(2,nstoke),aharm(2,nstoke),amp1,amp2,amp,anad(2,nstoke)

     real(4), save                          :: phirsv = 1.e30
     real(4), parameter                     :: thtref = 55.2	 
     real(4), dimension(2), parameter       :: xexp =(/4.,1.5/) ![MW 2012] section IV C.
     real(4), parameter, dimension(2,nstoke):: xexp_phir=reshape((/ 2.,2., 1.,4., 1.,4., 2.,2./), (/2,4/)) ![MW 2012] Table 5
     real(4), parameter, dimension(6)       :: freq0 =(/ 6.8,  10.7,  18.7,  23.8,  37.0, 85.5/)  
     !now referenced to windsat and ssmi

     if (freq<6.8 .or. freq>100.) then
     write(*,*) freq, ' emiss routine only set up for Windsat and AMSR freq.'
     stop
     endif
     
     
     if (present(phir)) then
        xphir=phir
     else
        xphir=-999.0
     endif

     if (present(sal)) then
        xsal=sal
     else
        xsal=35.0
     endif  

     qtht=tht
	 if(qtht.gt.65) qtht=65.  !qtht is just used for extrpolation for tht>thtref and i limit it to 60 deg

     if (present(emiss_0) .or. present(emiss_tot)) then
        call fdem0_meissner_wentz(freq=freq,tht=tht,sst=sst,salinity=xsal, em0=em0) 
     endif
      

    ifreq1=1
    if(freq.gt.freq0(2)) ifreq1=2
    if(freq.gt.freq0(3)) ifreq1=3  !between 18.7 and 37 ghz
    if(freq.gt.freq0(5)) ifreq1=5
 
    if(ifreq1.ne.3) then
        ifreq2=ifreq1+1
    else
        ifreq2=ifreq1+2
    endif
 
    wt=(freq-freq0(ifreq1))/(freq0(ifreq2)-freq0(ifreq1))
    if(freq.gt.freq0(ifreq2)) wt=1  !only occurs for freq>85.5

!   isotropic wind-induced emissivity
    if (present(emiss_wind) .or. present(emiss_tot)) then
        if(.not. (present(wind))) stop ' need wind speed for computing wind induced emissivity. pgm stopped.'
        call get_emiss_wind(ifreq1,sst,wind, emiss1)
        call get_emiss_wind(ifreq2,sst,wind, emiss2)
        xemiss=(1-wt)*emiss1 + wt*emiss2 
        !emiss is thtref value interpolated to input freq

        enad=0.5*(xemiss(1)+xemiss(2))

        do ipol=1,2
            if(tht.le.thtref) then
                xemiss(ipol)=enad         + (xemiss(ipol)-enad)*( tht/thtref)**xexp(ipol)
            else
                xemiss(ipol)=xemiss(ipol) + (xemiss(ipol)-enad)*(qtht-thtref)*xexp(ipol)/thtref  
            endif
        enddo  !ipol

    endif ! present isotropic wind induced emissivity

!   find emiss_phi
    if (present(emiss_phi) .or. present(emiss_tot) .or. present(eharm)) then
        if(.not.(present(wind))) stop ' need wind speed for computing wind induced emissivity. pgm stopped.'

        if(.not.present(eharm) .and.(xphir.lt.-998. .or. wind.le.3) ) then !-999. default for doing no correction
            xemiss_phi=0.0
        else  !find emiss_phi
            call get_aharm_phir(ifreq1,sst,wind, aharm1)  !aharm in terms of true stokes
            call get_aharm_phir(ifreq2,sst,wind, aharm2)  !aharm in terms of true stokes
            aharm=(1-wt)*aharm1 + wt*aharm2  !aharm is thtref value interpolated to input freq
            !     get nadir harmonic
            call get_aharm_phir_nad(ifreq1,freq0(ifreq1),sst,wind, amp1) 
            call get_aharm_phir_nad(ifreq2,freq0(ifreq2),sst,wind, amp2)  
            amp=(1-wt)*amp1 + wt*amp2
            anad=0  !most elements are zero
            anad(2,2)=  amp
            anad(2,3)= -amp
            do istoke=1,nstoke
            do iharm=1,2
                if(tht.le.thtref) then
                    aharm(iharm,istoke)= anad(iharm,istoke) +(aharm(iharm,istoke)-anad(iharm,istoke))* (tht/thtref)**xexp_phir(iharm,istoke)
                else
                    aharm(iharm,istoke)=aharm(iharm,istoke) +(aharm(iharm,istoke)-anad(iharm,istoke))*(qtht-thtref)*xexp_phir(iharm,istoke)/thtref  
                endif
            enddo  !iharm
            enddo  !istoke
            !     convert back from true stokes to v and h
            do iharm=1,2
                h1=aharm(iharm,1) + 0.5*aharm(iharm,2) !(v+h)/2 + (v-h)/2=v
                h2=aharm(iharm,1) - 0.5*aharm(iharm,2) !(v+h)/2 - (v-h)/2=h
                aharm(iharm,1)=h1
                aharm(iharm,2)=h2
            enddo

            !if(abs(xphir-phirsv).gt.0.01) then
                phirsv=xphir
                cos1phi=cosd(  xphir)
                cos2phi=cosd(2*xphir)
                sin1phi=sind(  xphir)
                sin2phi=sind(2*xphir)
            !endif

            xemiss_phi(1:2)=aharm(1,1:2)*cos1phi + aharm(2,1:2)*cos2phi
            xemiss_phi(3:4)=aharm(1,3:4)*sin1phi + aharm(2,3:4)*sin2phi
            xemiss_phi(4) = - xemiss_phi(4) ! IEEE convention

        endif ! (xphir.lt.-998. .or. wind.le.3)
        
    endif  ! emiss_phi

    if (present(emiss_tot))  then
        xemiss_tot = em0 + xemiss  + xemiss_phi(1:2)
        emiss_tot  = xemiss_tot
    endif
 
    if (present(emiss_0))    emiss_0   =em0
    if (present(emiss_wind)) emiss_wind=xemiss
    if (present(emiss_phi))  emiss_phi =xemiss_phi

    if (present(eharm)) then
        eharm=aharm
        eharm(:,4)=-aharm(:,4) ! IEEE
        if (wind<=3.0) eharm=0.0
    endif

return
end subroutine fd_emiss



subroutine get_emiss_wind(ifreq,sst,wind,   emiss)
! [MW 2012] section IV
implicit none

    integer(4), intent(in)                ::  ifreq
    real(4), intent(in)                   ::  sst,wind
    real(4), dimension(2), intent(out)    ::  emiss

    real(4), save                         ::  acoef(5,2,6)
    real(4)                               ::  sst_fac(2)
    real(8)                               ::  xmea(5)

    integer(4), save :: istart=1

    if(istart.eq.1) then
      istart=0
      open(unit=3,file=file_coeffs_wind_isotropic,form='binary',status='old',action='read')
      read(3) acoef
      close(3)
    endif

    if(ifreq.eq.4) stop 'ifreq oob in get_emiss_wind, pgm stopped'

    call  fd_xmea_win(wind, xmea)
    emiss(1)=dot_product(acoef(:,1,ifreq),xmea) 
    emiss(2)=dot_product(acoef(:,2,ifreq),xmea) 

    call get_sst_fac(ifreq,1,sst, sst_fac(1))
    call get_sst_fac(ifreq,2,sst, sst_fac(2))

    emiss=emiss*sst_fac


return
end subroutine get_emiss_wind




subroutine get_aharm_phir(ifreq,sst,wind,   aharm)
! [MW 2012] section VI
implicit none

    integer(4), parameter                           :: nstoke=4
    
    integer(4), intent(in)                          :: ifreq 
    real(4), intent(in)                             :: sst, wind
    
    real(4), dimension(2,nstoke), intent(out)       :: aharm

    integer(4)                                      :: istoke,iharm
    
    real(4)                                         :: h1,h2
    real(4), save                                   :: bcoef(5,2,nstoke,6)
    real(4)                                         :: sst_fac(nstoke)
    real(8)                                         :: xmea(5)

    integer(4), save  :: istart=1
 
    if(istart.eq.1) then
        istart=0
        open(unit=3,file=file_coeffs_wind_direction,form='binary',status='old',action='read')
        read(3) bcoef
        close(3) 
    endif

    if(ifreq.eq.4) stop 'ifreq oob in get_aharm_phir, pgm stopped'
    
    call  fd_xmea_win(wind, xmea)

    call get_sst_fac(ifreq,1,sst,   sst_fac(1))
    call get_sst_fac(ifreq,2,sst,   sst_fac(2))
    sst_fac(3:4)=0.5*(sst_fac(1)+sst_fac(2))

    do istoke=1,nstoke
        aharm(1,istoke)=sst_fac(istoke)*dot_product(xmea,bcoef(:,1,istoke,ifreq))
        aharm(2,istoke)=sst_fac(istoke)*dot_product(xmea,bcoef(:,2,istoke,ifreq))
    enddo  !istoke

!     convert to true stokes paramters,ie (v+h)/2 and v-h, rather than v and h in order to do tht adjustment
    do iharm=1,2
        h1=0.5*(aharm(iharm,1)+aharm(iharm,2))
        h2=     aharm(iharm,1)-aharm(iharm,2)
        aharm(iharm,1)=h1
        aharm(iharm,2)=h2
    enddo
  
return
end subroutine get_aharm_phir



subroutine get_aharm_phir_nad(ifreq,freq,sst,wind,      amp)
! [MW 2012], section VI] eq. (26).
implicit none

    integer(4), intent(in)    :: ifreq
    real(4), intent(in)       :: freq,sst,wind
    real(4), intent(out)      :: amp
    
    real(4)                   :: amp_10_nad,ywind,qfreq
    real(4)                   :: sst_fac

    qfreq=freq
    if(qfreq.gt.37) qfreq=37

    if(freq.lt.3) then
        amp_10_nad=.2/290.
    else
        amp_10_nad=2*(1. - 0.9*alog10(30./qfreq))/290.
    endif
 
    ywind=wind
    if(wind.lt. 0) ywind= 0
    if(wind.gt.15) ywind=15

    amp=amp_10_nad*ywind*(ywind - ywind**2/22.5)/55.5556
    call get_sst_fac(ifreq,0,sst, sst_fac)
    amp=amp*sst_fac

return
end subroutine get_aharm_phir_nad



subroutine get_sst_fac(ifreq,ipol,sst,      sst_fac)
! [MW 2012], eq. (15)
!ipol=0 denotes nadir value, sst_fac=em0(sst)/em0(sst=20)
!rcoef values for the ratio of nadir em0(sst)/em0(sst=20) were computed offline 
!so that they are available for fast computation

implicit none

    integer(4), intent(in)  :: ifreq,ipol
    real(4), intent(in)     :: sst
    
    real(4), intent(out)    :: sst_fac

    integer(4), save        :: istart=1
    real(4)                 :: xmea(3)
    real(4), save           :: rcoef(3,0:2,6)


    if(istart.eq.1) then
        istart=0
        open(unit=3,file=file_em0_ref_freq_sst,status='old',form='binary',action='read')
        read(3) rcoef
        close(3)
    endif

    xmea(1)= sst-20
    xmea(2)=xmea(1)*xmea(1)
    xmea(3)=xmea(1)*xmea(2)
    sst_fac=1 + dot_product(rcoef(:,ipol,ifreq),xmea)

return
end subroutine get_sst_fac


subroutine fd_xmea_win(wind,        xmea)
! Provides wind speed polynomials for [MW 2012, eqs. (14) + (25)] 
implicit none
 
    real(4), intent(in)                     :: wind
    real(8), dimension(5), intent(out)      :: xmea
    
    real(4)                                 :: x,dif
    real(4), parameter                      :: wcut =20.
    
    x=wind
    if(x.lt.0) x=0
    
    xmea(1)=x
    if(x.le.wcut) then
        xmea(2)=xmea(1)*x
        xmea(3)=xmea(2)*x
        xmea(4)=xmea(3)*x
        xmea(5)=xmea(4)*x
    else
        dif=x-wcut
        xmea(2)=2*dif*wcut       + wcut**2
        xmea(3)=3*dif*wcut**2    + wcut**3
        xmea(4)=4*dif*wcut**3    + wcut**4
        xmea(5)=5*dif*wcut**4    + wcut**5
    endif
 
return
end subroutine fd_xmea_win



subroutine fd_scatterm_all(freq,tht,wind,opacty,    xscat)
! [MW 2012], section V.
implicit none

    real(4), intent(in)                             :: freq,tht,wind,opacty
    real(4), dimension(2), intent(out)              :: xscat
    
    real(4)                                         :: xlog_freq,xscat1(2),xscat2(2)
    real(4)                                         :: a1,a2,b1,b2,c1,c2,brief,d1,d2
    real(4), save                                   :: scatterm(91,50,26,13,2)

    integer(4), save                                :: istart=1
    integer(4)                                      :: i1,i2,j1,j2,k1,k2,l1,l2


    if(istart.eq.1) then
        istart=0
        open(unit=3,file=file_coeffs_sctterm,status='old',form='binary',action='read')
        read(3) scatterm
        close(3)
    endif

    ! check inputs

    if(freq.lt.1 .or. freq.gt.200) stop 'freq oob in fd_scatterm, pgm stopped'
    if(tht .lt.0 .or.  tht.gt. 90) stop 'tht  oob in fd_scatterm, pgm stopped'
    if(wind.lt.0 .or. wind.gt.100) then    
    write(*,*) wind
    stop 'wind oob in fd_scatterm, pgm stopped'
    endif
    if(opacty.lt.0)                stop 'opacty oob in fd_scatterm, pgm stopped'

    xlog_freq=alog10(freq)
    
    ! multi-linear interpolation from table values
    
    brief=tht
    if(brief.gt.89.99) brief=89.99
    i1=1+brief
    i2=i1+1
    a1=i1-brief
    a2=1.-a1
    
    brief=wind
    if(brief.gt.24.99) brief=24.99
    j1=1+brief
    j2=j1+1
    b1=j1-brief
    b2=1-b1
    
    brief=xlog_freq/0.2
    if(brief.gt.11.99) brief=11.99
    k1=1+brief
    k2=k1+1
    c1=k1-brief
    c2=1-c1
    
    brief=opacty/0.025
    if(brief.gt.48.99) brief=48.99
    l1=1+brief
    l2=l1+1
    d1=l1-brief
    d2=1-d1
    
    xscat1= &
    a1*b1*(c1*scatterm(i1,l1,j1,k1,:)+c2*scatterm(i1,l1,j1,k2,:))+ &
    a1*b2*(c1*scatterm(i1,l1,j2,k1,:)+c2*scatterm(i1,l1,j2,k2,:))+ &
    a2*b1*(c1*scatterm(i2,l1,j1,k1,:)+c2*scatterm(i2,l1,j1,k2,:))+ &
    a2*b2*(c1*scatterm(i2,l1,j2,k1,:)+c2*scatterm(i2,l1,j2,k2,:))

    xscat2= &
    a1*b1*(c1*scatterm(i1,l2,j1,k1,:)+c2*scatterm(i1,l2,j1,k2,:))+ &
    a1*b2*(c1*scatterm(i1,l2,j2,k1,:)+c2*scatterm(i1,l2,j2,k2,:))+ &
    a2*b1*(c1*scatterm(i2,l2,j1,k1,:)+c2*scatterm(i2,l2,j1,k2,:))+ &
    a2*b2*(c1*scatterm(i2,l2,j2,k1,:)+c2*scatterm(i2,l2,j2,k2,:))
    
    xscat=d1*xscat1 + d2*xscat2
 
return
end subroutine fd_scatterm_all



subroutine fd_emiss_LBAND (wspd, sst, tht,    emiss_wind, eharm)
! wind roughness emissivity
! Aquarius V5.0 for V/H. 
! SMAP V3.0 for S3/S4
! interpolate/extrapolate EIA.


implicit none

real(4), intent(in)                                 ::  wspd ! [m/s]
real(4), intent(in)                                 ::  sst  ! SST [Celsius]
real(4), intent(in)                                 ::  tht  ! Earth Incidence Angle [deg]

real(4), intent(out), dimension(2), optional        ::  emiss_wind 
! isotropic part of wind induced surface emissivity [0,1]

real(4), intent(out), dimension(2,4), optional      ::  eharm
! directional signal [0,1]

real(4), dimension(0:2,4)                           ::  xeharm

real(4), dimension(0:2,4)                           ::  aharm_L_0, aharm_L_1, aharm_L_2, aharm_L_3, aharm_L_x
real(4), dimension(0:2,4)                           ::  daharm


real(4), dimension(0:2,4)                           ::  yharm_0, yharm_1
real(4), dimension(4)                               ::  xem0, yem0
real(4)                                             ::  brief
real(4), dimension(4)                               ::  yy

real(4), dimension(0:n_rad)                         ::  thtfix

integer(4)                                          ::  iharm, ipol, irad_0, irad_1

    ! SST form factor
    ! assume geometric optics as we do it for higher frequencies
    call fdem0_meissner_wentz(freq_aq,tht,sst,     sss_ref0, xem0(1:2)) 
    call fdem0_meissner_wentz(freq_aq,tht,sst_ref0,sss_ref0, yem0(1:2)) 
    xem0(3) = (xem0(1)+xem0(2))/2.0 
    xem0(4) = (xem0(1)+xem0(2))/2.0     
    yem0(3) = (yem0(1)+yem0(2))/2.0 
    yem0(4) = (yem0(1)+yem0(2))/2.0
         

    ! Aquarius V5.0 for V/H
    thtfix(1:3) = tht_ref_AQ(1:3)
    thtfix(0)   = 0.0  
    call fd_AQ_emiss_harmonics(1,wspd,		    aharm_L_1, daharm) 
    call fd_AQ_emiss_harmonics(2,wspd,		    aharm_L_2, daharm) 
    call fd_AQ_emiss_harmonics(3,wspd,		    aharm_L_3, daharm) 
    
    ! nadir
    aharm_L_0(0,1) = (aharm_L_1(0,1)+aharm_L_1(0,2))/2.0  
    !I assume (V+H)/2 of horn 1 at nadir
    aharm_L_0(0,2) = (aharm_L_1(0,1)+aharm_L_1(0,2))/2.0  
    !I assume (V+H)/2 of horn 1 at nadir
    aharm_L_0(1,1) = (aharm_L_1(1,1)+aharm_L_1(1,2))/2.0  
    !I assume (V+H)/2 of horn 1 at nadir
    aharm_L_0(1,2) = (aharm_L_1(1,1)+aharm_L_1(1,2))/2.0  
    !I assume (V+H)/2 of horn 1 at nadir    
    aharm_L_0(2,1) = (aharm_L_1(2,1)-aharm_L_1(2,2))/2.0  
    !This guarantess that A2(S1) = 0 at nadir, as required by reflection symmetry of Maxwell eqs. (Yueh)  
    aharm_L_0(2,2) = (aharm_L_1(2,2)-aharm_L_1(2,1))/2.0  
    !This guarantess that A2(S1) = 0 at nadir, as required by reflection symmetry of Maxwell eqs. (Yueh)     
    ! note that A2(V) and A2(H) have opposite signs.

    do iharm=0,2
    do ipol=1,2
    if (tht<tht_ref_AQ(1)) then
        irad_0 = 0
        irad_1 = 1
        yharm_1(iharm,ipol) = aharm_L_1(iharm,ipol) ! *290K
        yharm_0(iharm,ipol) = aharm_L_0(iharm,ipol) ! *290K
    else if (tht>=tht_ref_AQ(1) .and. tht<tht_ref_AQ(2)) then
        irad_0 = 1
        irad_1 = 2   
        yharm_1(iharm,ipol) = aharm_L_2(iharm,ipol) ! *290K
        yharm_0(iharm,ipol) = aharm_L_1(iharm,ipol) ! *290K       
    else 
        irad_0 = 2
        irad_1 = 3    
        yharm_1(iharm,ipol) = aharm_L_3(iharm,ipol) ! *290K
        yharm_0(iharm,ipol) = aharm_L_2(iharm,ipol) ! *290K                    
    endif    
         
    brief = (tht-thtfix(irad_0))/(thtfix(irad_1)-thtfix(irad_0)) 
    ! linear interpolation/extrapolation to/form Aquarius EIA

    yy(ipol) = yharm_0(iharm,ipol)*(1.0-brief) + yharm_1(iharm,ipol)*brief          
    xeharm(iharm,ipol) = yy(ipol)*(xem0(ipol)/yem0(ipol))/teff0    
    enddo !ipol   
    enddo !iharm


    ! SMAP for S3/S4
    call fd_SMAP_harmonics(wspd,		aharm_L_x) ! EIA=40 
    
    ! I assume that there is no S4 1st harmonic at any EIA (even though the SMAP harmoniscs show one).
    aharm_L_0(1,4) = 0.0
    aharm_L_x(1,4) = 0.0
    
    ! A2(S4) = 0 at nadir, as required by reflection symmetry of Maxwell eqs. (Yueh)     
    aharm_L_0(2,4) = 0.0
    
    ! I assume that A1(S3) at nadir is the same as for SMAP
    aharm_L_0(1,3) = aharm_L_x(1,3)
    
    ! A2(S3) = -A2(S2) (S2=V-H) at nadir, as required by reflection symmetry of Maxwell eqs. (Yueh) 
    aharm_L_0(2,3) = aharm_L_0(2,2)-aharm_L_0(2,1)

    do iharm=1,2
    do ipol=3,4
           
    brief = (tht-tht_ref_SMAP)/(0.0-tht_ref_SMAP) 
    ! linear interpolation/extrapolation to/form SMAP EIA

    yy(ipol) = aharm_L_x(iharm,ipol)*(1.0-brief) + aharm_L_0(iharm,ipol)*brief  
    xeharm(iharm,ipol) = yy(ipol)*(xem0(ipol)/yem0(ipol))/teff0    
    enddo !ipol   
    enddo !iharm
   
    
    if (present(emiss_wind)) then
        emiss_wind(1:2)=xeharm(0,1:2)
    endif

    if (present(eharm)) then
        do iharm=1,2
        eharm(iharm,1:4)=xeharm(iharm,1:4)
        enddo
    endif


return
end subroutine fd_emiss_LBAND



! Aquarius V5.0
subroutine fd_AQ_emiss_harmonics(irad,wspd,		aharm, daharm) 
implicit none

	integer(4), intent(in)					            :: irad
	real(4), intent(in)						            :: wspd
	
    real(4), dimension(0:2,2), intent(out)	            ::  aharm(0:2,2)  !1=V, 2=H
    real(4), dimension(0:2,2), intent(out), optional   	:: daharm(0:2,2)  !1=V, 2=H
 
    integer(4), parameter                               :: n_rad=3, npoly=5


	real(4), dimension(2)				                ::  A0,  A1,  A2  !1=V, 2=H
	real(4), dimension(2)				                :: dA0, dA1, dA2  !1=V, 2=H
	
	
	real(4)									            :: ww
	integer(4)								            :: ipol, iharm
	real(4)									            :: fval, dval

	real(4)         						            :: w0, w1, w2 ! linear extrapolation/cutoff points
	
    
	integer(4), save                                    :: istart=1	
    real(8), dimension(0:2,2,n_rad,npoly), save         :: acoef      ! harmonic coefficients for radiometer wind direction signal
    real(8), dimension(0:2,2,n_rad), save               :: wspd_max_a ! high wind speed for radiometer wind speed signal

	
	if (istart==1) then
	    istart=0
	    open(unit=3,file=emiss_coeff_harm_file_AQ, form='binary', action='read', status='old')
	    read(3) acoef
	    read(3) wspd_max_a
	    ! overwrite
	    wspd_max_a = wextrapol
	    close(3)
	endif 
	

	A0=0.0
	A1=0.0
	A2=0.0

	do ipol=1,2 

	! A0
	iharm=0
	w0 = wspd_max_a(iharm,ipol,irad)
	ww = wspd
	if (wspd >= w0) ww=w0 ! extrapolation at w0 
	fval = &
	ww*acoef(iharm,ipol,irad,1)  +  (ww**2)*acoef(iharm,ipol,irad,2) +       (ww**3)*acoef(iharm,ipol,irad,3) +       (ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5) 		   
	dval = &
	   acoef(iharm,ipol,irad,1)  + (2.0*ww)*acoef(iharm,ipol,irad,2) + (3.0*(ww**2))*acoef(iharm,ipol,irad,3) + (4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5) 
	
	if (wspd<=w0) then
		A0(ipol) = fval
	else
		A0(ipol) = fval + dval*(wspd-w0)
	endif
    
    dA0(ipol) = dval
   
    
	! A1
	iharm=1
	w1 = wspd_max_a(iharm,ipol,irad)
	ww = wspd
	if (wspd >= w1) ww=w1 ! cutoff at w1 
	fval = ww*acoef(iharm,ipol,irad,1) + (ww**2) *acoef(iharm,ipol,irad,2) +       (ww**3)*acoef(iharm,ipol,irad,3)      +       (ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5)   
	dval =    acoef(iharm,ipol,irad,1) + (2.0*ww)*acoef(iharm,ipol,irad,2) + (3.0*(ww**2))*acoef(iharm,ipol,irad,3)      + (4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5)  
	
	 A1(ipol)  = fval
	dA1(ipol)  = dval

	! A2
	iharm=2
	w2 = wspd_max_a(iharm,ipol,irad)
	ww = wspd
	if (wspd >= w2) ww=w2 ! cutoff at w2 
	fval = ww*acoef(iharm,ipol,irad,1) + (ww**2)      *acoef(iharm,ipol,irad,2) +        (ww**3)*acoef(iharm,ipol,irad,3) +       (ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5)  
	dval =    acoef(iharm,ipol,irad,1) + (2.0*ww)     *acoef(iharm,ipol,irad,2) +  (3.0*(ww**2))*acoef(iharm,ipol,irad,3) + (4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5) 

	A2(ipol)  = fval
	dA2(ipol) = dval

	enddo !ipol


	do ipol=1,2
		aharm(0,ipol) = A0(ipol)
		aharm(1,ipol) = A1(ipol)
		aharm(2,ipol) = A2(ipol)
	enddo

	if (present(daharm)) then
	do ipol=1,2
		daharm(0,ipol) = dA0(ipol)
		daharm(1,ipol) = dA1(ipol)
		daharm(2,ipol) = dA2(ipol)
	enddo
	endif


return
end subroutine fd_AQ_emiss_harmonics


subroutine fd_SMAP_harmonics(wspd,		aharm) 
implicit none

	real(4), intent(in)						            :: wspd
	
    real(4), dimension(0:2,4), intent(out)	            :: aharm  !1=V, 2=H, 3=S3, 4=S4
    ! the 0th harmonic is set ot 0
 
    integer(4), parameter                               :: npoly=5


	real(4), dimension(4)				                :: A0,  A1,  A2  !1=V, 2=H, 3=S3, 4=S4
	
	
	real(4)									            :: ww
	integer(4)								            :: ipol, iharm
	real(4)									            :: fval

	real(4), parameter  					            :: wcut=24.5 ! cutoff for 1st and 2nd harm 

	
	integer(4), save                                    :: istart=1	
    real(8), dimension(0:2,4,npoly), save               :: acoef      
    ! harmonic coefficients for radiometer wind direction signal
    ! V/H/S3/S4
	
	if (istart==1) then
	    istart=0
	    open(unit=3,file=emiss_coeff_harm_file_SMAP, form='binary', action='read', status='old')
	    read(3) acoef
	    close(3)
	endif 
	
	A0=0.0
	A1=0.0
	A2=0.0

	do ipol=1,4 

	! A0 set to 0
	! The A0 is taken form the Aquarius V5 GMF
	A0(ipol) = 0.0
    
	! A1
	iharm=1
	ww = wspd
	if (wspd >= wcut) ww=wcut ! cutoff at wcut 
	fval = ww*acoef(iharm,ipol,1) + (ww**2) *acoef(iharm,ipol,2) +       (ww**3)*acoef(iharm,ipol,3) +   (ww**4)*acoef(iharm,ipol,4) +   (ww**5)*acoef(iharm,ipol,5)   
    A1(ipol)  = fval

	! A2
	iharm=2
    ww = wspd
	if (wspd >= wcut) ww=wcut ! cutoff at wcut 
	fval = ww*acoef(iharm,ipol,1) + (ww**2) *acoef(iharm,ipol,2) +       (ww**3)*acoef(iharm,ipol,3) +   (ww**4)*acoef(iharm,ipol,4) +   (ww**5)*acoef(iharm,ipol,5)  
	A2(ipol)  = fval

	enddo !ipol

	do ipol=1,4
		aharm(0,ipol) = A0(ipol)
		aharm(1,ipol) = A1(ipol)
		aharm(2,ipol) = A2(ipol)
	enddo

return
end subroutine fd_SMAP_harmonics





subroutine fdem0_meissner_wentz(freq,tht,sst,salinity,      em0)
! Compute specular emissivity using Frsenel equations and MW dielectric model.
implicit none

    real(4), intent(in)                         :: freq,tht,sst,salinity

    real(4), dimension(2), intent(out)          :: em0
    
    real(4), parameter                          :: f0=17.97510
 

    real(4)                                     :: costht,sinsqtht
    real(4)                                     :: e0s,e1s,e2s,n1s,n2s,sig
    
    complex(4)                                  :: permit,esqrt,rh,rv
    complex(4), parameter                       :: j=(0.,1.)
    
 
    call dielectric_meissner_wentz(sst,salinity,  e0s,e1s,e2s,n1s,n2s,sig)

    costht=cosd(tht)
    sinsqtht=1.-costht*costht


!   debye law (2 relaxation wavelengths)
    permit = (e0s - e1s)/(1.0 - j*(freq/n1s)) + (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s + j*sig*f0/freq
    permit = conjg(permit)
    
    esqrt=csqrt(permit-sinsqtht)
    rh=(costht-esqrt)/(costht+esqrt)
    rv=(permit*costht-esqrt)/(permit*costht+esqrt)
    em0(1)  =1.-rv*conjg(rv)
    em0(2)  =1.-rh*conjg(rh)
 
return
end subroutine fdem0_meissner_wentz

 


subroutine dielectric_meissner_wentz(sst_in,s,   e0s,e1s,e2s,n1s,n2s,sig)
!
!     complex dielectric constant: eps
!     [MW 2004, MW 2012].
!     
!     Changes from [MW 2012]:
!     1. Typo (sign) in the printed version of coefficient d3 in Table 7. Its value should be -0.35594E-06.
!     2. Changed SST behavior of coefficient b2 from:
!     b2 = 1.0 + s*(z(10) + z(11)*sst) to
!     b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30)) 
!
!!
!     input:
!     name   parameter  unit  range
!     sst      sst        [c]   -25 c to 40 c for pure water
!                               -2  c to 34 c for saline water
!     s      salinity   [ppt]  0 to 40
!
!     output:
!     eps    complex dielectric constant
!            negative imaginary part to be consistent with wentz1 convention
!

implicit none


    real(4), intent(in)  :: sst_in,s
    real(4), intent(out) :: e0s,e1s,e2s,n1s,n2s,sig
 
    real(4), dimension(11), parameter :: &
      x=(/ 5.7230e+00, 2.2379e-02, -7.1237e-04, 5.0478e+00, -7.0315e-02, 6.0059e-04, 3.6143e+00, &
           2.8841e-02, 1.3652e-01,  1.4825e-03, 2.4166e-04 /)
    
    real(4), dimension(13), parameter :: &
      z=(/ -3.56417e-03,  4.74868e-06,  1.15574e-05,  2.39357e-03, -3.13530e-05, &
            2.52477e-07, -6.28908e-03,  1.76032e-04, -9.22144e-05, -1.99723e-02, &
            1.81176e-04, -2.04265e-03,  1.57883e-04  /)  ! 2004

    real(4), dimension(3), parameter :: a0coef=(/ -0.33330E-02,  4.74868e-06,  0.0e+00/)
    real(4), dimension(5), parameter :: b1coef=(/0.23232E-02, -0.79208E-04, 0.36764E-05, -0.35594E-06, 0.89795E-08/)
 
    real(4) :: e0,e1,e2,n1,n2
    real(4) :: a0,a1,a2,b1,b2
    real(4) :: sig35,r15,rtr15,alpha0,alpha1

    real(4) :: sst,sst2,sst3,sst4,s2
    
    sst=sst_in
    if(sst.lt.-30.16) sst=-30.16  !protects against n1 and n2 going zero for very cold water
    
    sst2=sst*sst
    sst3=sst2*sst
    sst4=sst3*sst

    s2=s*s
 
    !     pure water
    e0    = (3.70886e4 - 8.2168e1*sst)/(4.21854e2 + sst) ! stogryn et al.
    e1    = x(1) + x(2)*sst + x(3)*sst2
    n1    = (45.00 + sst)/(x(4) + x(5)*sst + x(6)*sst2)
    e2    = x(7) + x(8)*sst
    n2    = (45.00 + sst)/(x(9) + x(10)*sst + x(11)*sst2)
    
    !     saline water
    !     conductivity [s/m] taken from stogryn et al.
    sig35 = 2.903602 + 8.60700e-2*sst + 4.738817e-4*sst2 - 2.9910e-6*sst3 + 4.3047e-9*sst4
    r15   = s*(37.5109+5.45216*s+1.4409e-2*s2)/(1004.75+182.283*s+s2)
    alpha0 = (6.9431+3.2841*s-9.9486e-2*s2)/(84.850+69.024*s+s2)
    alpha1 = 49.843 - 0.2276*s + 0.198e-2*s2
    rtr15 = 1.0 + (sst-15.0)*alpha0/(alpha1+sst)
    
    sig = sig35*r15*rtr15
    
    !    permittivity
    a0 = exp(a0coef(1)*s + a0coef(2)*s2 + a0coef(3)*s*sst)  
    e0s = a0*e0
    
    if(sst.le.30) then
        b1 = 1.0 + s*(b1coef(1) + b1coef(2)*sst + b1coef(3)*sst2 + b1coef(4)*sst3 + b1coef(5)*sst4)
    else
        b1 = 1.0 + s*(9.1873715e-04 + 1.5012396e-04*(sst-30))
    endif
      
    n1s = n1*b1
    
    a1  = exp(z(7)*s + z(8)*s2 + z(9)*s*sst)
    e1s = e1*a1

    b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30))
    n2s = n2*b2
    
    a2 = 1.0  + s*(z(12) + z(13)*sst)
    e2s = e2*a2
    
return
end subroutine  dielectric_meissner_wentz

 

subroutine fd_tcos_eff(freq,    tcos_eff)
!     Calculates effective cold space temperature taking into account 
!     the deviation between Rayleigh-Jeans and Planck law
!     as function of frequency [MWD 2011], Appendix D.
!   
!     for these routine the term b is the flux 
!    (2*h*f**3/(c**2*(dexp(h*f/(k*t))-1))) divided by  2*k*f**2/c**2
implicit none


    real(4), intent(in)         :: freq
    real(4), intent(out)        :: tcos_eff


    real(8), parameter          :: tcos=2.73
    real(8), parameter          :: teff=63.  !selected to provide optimum fit over 60-300 k range
    real(8), parameter          :: h=6.6260755d-34
    real(8), parameter          :: k= 1.380658d-23
    real(8), parameter          :: a=h/k
    real(8)                     :: x,b1,b2

    x=a*freq*1.d9
    b1=x/(dexp(x/tcos)-1)
    b2=x/(dexp(x/teff)-1)
    tcos_eff=b1-b2+teff

return
end subroutine fd_tcos_eff







end module RSS_SURFACE_EMISSION_RTM_ALLFREQ_MODULE
