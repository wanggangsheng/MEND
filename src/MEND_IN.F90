SUBROUTINE MENDIN(sPAR_SCE,sINI)
!c$debug
!c
!c   THIS SUBROUTINE READS 
!c    (1) MEND MODEL INITILIZATION
!c    (2) INPUT DATA/CALIBRATION DATA
!c    (3) INPUT VARIABLES FOR MODEL OPTIMIZATION
!c    (4) MODEL PARAMETERS 
!c  
!c  AUTHOR: GANGSHENG WANG
!C  Environmental Sciences Division
!C  Oak Ridge National Laboratory
!C  Oak Ridge, TN 37831-6301
!C  March, 2013 
!C  Updated: May 5, 2015
!C  Updated: Apr 6, 2017, use 'MEND_namelist.nlm' (created by Junyi Liang) to replace 'MEND.ini' & 'MENDcn.ini'

    USE MOD_OPT_TYPE
    USE MOD_MEND_TYPE
    USE MOD_MEND,   only: fSWC2SWP,fSWP2SWC,sINP_Read,subMEND_INI_Read,sOUT_ALL_tscale!!sOUT_OPT  !!function
    USE MOD_String, only: StrCompress
    USE MOD_USRFS,  only: nDaysbwDates,nMonsbwDates,nDaysofMon
    USE MOD_USRFS,  only: sDate_After,sInt2Str
    USE MOD_USRFS,  only: Array_Normalize !!function 

    IMPLICIT NONE
    !!ARGUMENTS:
    TYPE(sSCE_PAR), intent(inout):: sPAR_SCE
    TYPE(sMEND_INI),intent(inout):: sINI

    !!LOCAL VARIABLES:
    INTEGER iFin  !input file unit

    real pcenta
    integer, allocatable:: iPar(:), iPar_opt(:)
    character*10 pcntrl,deflt,usrsp
    character*4 reduc,initl,ysflg,noflg
    integer i, j, k, lp, ierror, iwarn
    integer eof !!end of file
    INTEGER is_total  !!=1: need to convert to hourly rate; =0: directly assign value

    data deflt/' DEFAULT  '/
    data usrsp/'USER SPEC.'/
    data ysflg/'YES '/
    data noflg/'NO  '/


    INTEGER             :: iRead
    real(8)             :: rRead, rRead2
    character(len=50)   :: sRead,sRead1
    character(len=200)  :: sRead2
    
    integer             :: ndays,nmons  
    CHARACTER(LEN=2)    :: str2
    character(len=8)    :: sDate 
    character(len=10)   :: sUnits,ststep
    character(len=200)  :: sfilename_full,sFile_inp,sFile_out
    
    character(len=20)   :: Name_POOL(const_nPOOL), Name_MNPOOL(const_nPOOL_MN)
    character(len=20)   :: Name_FLUX(const_nFLUX), Name_MNFLUX(const_nFLUX_MN)
    character(len=20)   :: Name_PAR(const_nPAR), Name_RATE(const_nRATE)
  
    !!------------------------------------------------------------------------
    !! namelist variables: beg
    integer             :: ierr
    character(len=50)   :: msg
    character(len=20)   :: Pname(const_nPAR0)
    real(8)             :: Pinitial(const_nPAR0),Plow(const_nPAR0),Phigh(const_nPAR0)
    integer             :: Pcal(const_nPAR0)
    CHARACTER(len=10)   :: sSite
    CHARACTER(len=5)    :: sBIOME
    CHARACTER(len=5)    :: sSOM
    REAL(8)             :: Altitude, GPPref
    INTEGER             :: iGPPscaler
    integer             :: iModel
    integer             :: iSA_range(2)
    CHARACTER(len=3)    :: ssMEND
    CHARACTER(len=20)   :: Dir_Input, Dir_Output
    character(len=8)    :: ssDate_beg_all, ssDate_end_all, ssDate_beg_sim, ssDate_end_sim

    character(len=20)   :: sfilename_ST(20), sfilename_SM(20), sfilename_type1(20), sfilename_pH(20)
    integer             :: ifdata_ST, nfile_ST
    character(len=10)   :: sUnits_ST
    character(len=10)   :: step_ST

    real(8)             :: vg_SWCres,vg_SWCsat,vg_alpha,vg_n, Ksat, Lambda!!van-Genuchten equation
    integer             :: ifdata_SM, nfile_SM
    character(len=10)   :: sUnits_SM
    character(len=10)   :: step_SM
    integer             :: ifdata_type1, nfile_type1
    character(len=10)   :: sUnits_type1
    character(len=10)   :: step_type1
    real(8)             :: Fraction_type1(3), Input_type2(3), Input_type3(3)
    character(len=8)    :: ssDate_beg_inp2, ssDate_end_inp2
    REAL(8)             :: sSIN_C12_C14, sSIN_Multiplier
    integer             :: ifdata_pH, nfile_pH
    character(len=10)   :: sUnits_pH
    character(len=10)   :: step_pH
    real(8)             :: spH_constant
    CHARACTER(len=20)   :: sSOIL_INI_file !soil initialization file
    integer             :: siScenario
    real(8)             :: sSTP_delta
    real(8)             :: sSWC_logis(3)
    real(8)             :: sSIN_logis (4)
!    integer nPar !, nVAR
    CHARACTER(len=20)   :: Cali_var_title(10)
    integer             :: Cali_varid(const_nVAR0)
    CHARACTER(len=20)   :: Cali_VAR(const_nVAR0)
    CHARACTER(len=20)   :: Cali_Units(const_nVAR0)
    integer             :: Cali_Calibrate(const_nVAR0)
    integer             :: Cali_tstep(const_nVAR0)
    CHARACTER(len=20)   :: Cali_obs_file(const_nVAR0)
    integer             :: Cali_obs_file_column(const_nVAR0)
    CHARACTER(len=20)   :: Cali_OBJ(const_nVAR0)
    integer             :: Cali_OBJ_Weight(const_nVAR0)
    real(8)             :: Cali_OBJ_Tolerance(const_nVAR0)
    integer             :: SCE_Parameters(5)
    integer             :: SCE_control_Parameters(6)
    integer             :: siKinetics, siHR, siTmp_Func, siENZN_Allocate
    real(8)             :: CN_ratio_input(3)
    integer             :: ifdata_NH4, nfile_NH4
    character(len=10)   :: sUnits_NH4
    character(len=10)   :: step_NH4
    integer             :: ifdata_NO3, nfile_NO3
    character(len=10)   :: sUnits_NO3
    character(len=10)   :: step_NO3
    character(len=20)   :: sfilename_NH4(1), sfilename_NO3(1)
    real(8)             :: ST_constant, SM_constant, Input_type1_constant, NH4_constant, NO3_constant

    namelist /mend_config/ sSite, sBIOME, sSOM, iModel, ssMEND, iSA_range, Altitude,GPPref, iGPPscaler, &
                       Dir_Input, Dir_Output, &
                       ssDate_beg_all, ssDate_end_all, ssDate_beg_sim, ssDate_end_sim, &
                       ifdata_ST, sUnits_ST, step_ST, nfile_ST, sfilename_ST, ST_constant, &
                       vg_SWCres, vg_SWCsat, vg_alpha, vg_n, Ksat,Lambda,&
                       ifdata_SM, sUnits_SM, step_SM, nfile_SM, sfilename_SM, SM_constant, &
                       ifdata_type1, sUnits_type1, step_type1, nfile_type1, sfilename_type1, Input_type1_constant, &
                       Fraction_type1, Input_type2, Input_type3, ssDate_beg_inp2, ssDate_end_inp2, &
                       sSIN_C12_C14, sSIN_Multiplier, &
                       ifdata_pH, sUnits_pH, step_pH, nfile_pH, sfilename_pH, spH_constant, &
                       sSOIL_INI_file, siScenario, sSTP_delta, sSWC_logis, sSIN_logis, &
                       Cali_var_title, Cali_varid, Cali_VAR, Cali_Units, Cali_Calibrate, &
                       Cali_tstep, Cali_obs_file, Cali_obs_file_column, Cali_OBJ, Cali_OBJ_Weight, Cali_OBJ_Tolerance, &
                       SCE_Parameters, SCE_control_Parameters, &
                       siKinetics, siHR, siTmp_Func, siENZN_Allocate, &
                       Pname, Pinitial, Plow, Phigh, Pcal, CN_ratio_input, &
                       ifdata_NH4, sUnits_NH4, step_NH4, nfile_NH4, sfilename_NH4, NH4_constant, &
                       ifdata_NO3, sUnits_NO3, step_NO3, nfile_NO3, sfilename_NO3, NO3_constant
    !! nVar/nPar has been defined in 'MOD_MEND_TYPE.F90': const_nVAR0, const_nPAR0
    !! namelist variables: end
    !!------------------------------------------------------------------------
      
    !!VARIABLES NAMES FOR OUTPUT FILES: BEG
    DATA Name_POOL    /"TM","TOM","SOM","POMO","POMH","POMT","MOMT","MOM","QOM",&
                       "DOM","DOMS","DOMT","MB","MBA","MBD",&
                       "ENZ","ENZD","ENZP1","ENZP2","ENZM",&
                       "ENZNmt","ENZNfix","ENZNH4","ENZNO3","ENZNO2","ENZNO","ENZN2O","TM_err"/
                       
    DATA Name_MNPOOL  /"CO2","Nmine","Nmine_Free","Nmine_Solid","NH4tot","NH4ads","NH4", &
                       "NOx","NO3+NO2","NO3","NO2","NO","N2O","N2", "NGas"/
                       
    DATA Name_FLUX    /"TOTout","TOTinp","POMinp1","POMinp2","DOMinp",&
                       "POMdec1","POMdec2","POMdec2DOM1","POMdec2DOM2","POMdec2MOM1","POMdec2MOM2",&
                       "MOMdec","MOM2DOM","QOM2DOM","DOM2QOM","DOM2QOMnet","QOM2MOM","MOM2QOM","MOM2QOMnet",&
                       "DOM2MBA","MBA_mortality", "MBA2EP1","MBA2EP2","MBA2EM","MBA_PM",&
                       "EP2DOM1","EP2DOM2","EM2DOM","MBA2DOM","MBA2POM1","MBA2POM2","MBA2MBD","MBD2MBA",&
                       "MBA2ENZNm1","MBA2ENZNm2","MBA2ENZNm3","MBA2ENZNm4","MBA2ENZNm5","MBA2ENZNm6",&
                       "ENZNm2DOM1","ENZNm2DOM2","ENZNm2DOM3","ENZNm2DOM4","ENZNm2DOM5","ENZNm2DOM6"/
                       
    DATA Name_MNFLUX  /"CO2_soil","CO2_root","CO2_gmo","CO2_gm","CO2_growth","CO2_maintn","CO2_ovflow", &
                       "CO2_maintn_MBA","CO2_maintn_MBD","CO2_ovflow_MBA","CO2_ovflow_MBD", &
                       "Nmine_dep","NH4_dep","NO3_dep", &
                       "Nmn_net","Nmn","Nmn_MBA","Nmn_MBD","Nim","Nim_NH4","Nim_NO3", &
                       "NFix","Nitrif","Denit_NO3","Denit_NO2","Denit_NO","Denit_N2O", &
                       "NO_efflux","N2O_efflux","N2_efflux", &
                       "Nim_VG","Nim_NH4_VG","Nim_NO3_VG",&
                       "NO32_Leaching","NO3_Leaching","NO2_Leaching"/
                       
    DATA Name_PAR     /"fRa","fINP","VP1","VP2","VM","KP1","KP2","KM","Qmax","Kba","Kdes","Kads","Kp2u","Ku2p", &
                       "rEP1","rEP2","rEM","pEP","pEM","fD","gD",&
                       "Vg","alpha","Vm","KD","Yg","Q10","gamma","rMORT",&
                       "beta","Vm_dorm","VmA2D","VmD2A","SWP_A2D","tau","SWP_D2A","wdorm", &
                       "VNup_MB","VNup_VG","rNleach","bNup_VG","KsNH4_MB","KsNO3_MB","YgN", &
                       "Qmax_NH4","Kba_NH4","KsNH4_VG","KsNO3_VG", &
                       "fpENZN","VNif","VNit","VDenit1", "VDenit2","VDenit3","VDenit4",&
                       "KsNif","KsNit","KsDenit1", "KsDenit2","KsDenit3","KsDenit4"/
                       
    DATA Name_RATE    /"kSOM","kPOM1","kPOM2","kPOM","kMOM","kQOM","kDOM","kMBA","kMBA_in","kMBD",&
                       "kMBD_in","kMB","kMB_in","phi","Active_Fraction","CUE","NUE", &
                       "rNH4","fDOM","fNO3_Leaching", "GPPscaler_VG",&
                       "EA_POM1","EA_POM2","EA_MOM",&
                       "EA_NFix","EA_Nit","EA_Denit_NO3","EA_Denit_NO2","EA_Denit_NO","EA_Denit_N2O",&
                       "kNFix","kNit","kDenit_NO3","kDenit_NO2","kDenit_NO","kDenit_N2O",&
                       "TM_err","TMbeg","TMend","TOTinp","TOTout", &
                       "STP", "SWC", "SWP", "pH"/
    !!VARIABLES NAMES FOR OUTPUT FILES: END

    write (*,*) '>>ENTER SUBROUTINE <MEND_INI>'
    sINI%Name_POOL = Name_POOL
    sINI%Name_MNPOOL = Name_MNPOOL
    sINI%Name_FLUX = Name_FLUX
    sINI%Name_MNFLUX = Name_MNFLUX
    sINI%Name_PAR = Name_PAR
    sINI%Name_RATE = Name_RATE
      
!c  INITIALIZE I/O VARIABLES
!      iFin = 10
!      open(unit=iFin,file='MEND.ini',status='old')

    ierror = 0
    iwarn = 0
      
!!--------------------------------------------------------------------------------------
      !Read the namelist
    !ierr = 0
!        print*, 'ierr = ', ierr
    open (10,file="MEND_namelist.nml",status='OLD',recl=80,delim='APOSTROPHE')

    read(10,nml=mend_config, iostat=ierr,iomsg=msg)
    if(ierr/=0) then
        print*,"Namelist Error",ierr
        print*, msg
        stop
    end if
!!--------------------------------------------------------------------------------------
      
    !Define sINI values
    sINI%SITE = sSite
    sINI%altitude = Altitude
    sINI%BIOME = sBIOME
    sINI%SOM = sSOM
    sINI%iModel = iModel
    sINI%iSA_range = iSA_range
    
    
    if(StrCompress(ssMEND).eq.'C') then
          sINI%Carbon_only = .true.
    else
          sINI%Carbon_only = .false.
    end if
    sINI%dirinp = StrCompress(Dir_Input)
    sINI%dirout = StrCompress(Dir_Output)

!    CALL system('mkdir '//sINI%dirout)
    sRead2 = "mkdir "//trim(sINI%dirout)
    CALL SYSTEM(trim(sRead2))

    sINI%sDate_beg_all = ssDate_beg_all
    sINI%sDate_end_all = ssDate_end_all
    sINI%sDate_beg_sim = ssDate_beg_sim
    sINI%sDate_end_sim = ssDate_end_sim

    sINI%SIN_frac(1:3) = Fraction_type1
    sINI%SIN_other(1,1:3) = Input_type2
    sINI%SIN_other(2,1:3) = Input_type3
    sINI%sDate_beg_inp2 = ssDate_beg_inp2
    sINI%sDate_end_inp2 = ssDate_end_inp2

    sINI%SIN_C12_C14 = sSIN_C12_C14
    sINI%SIN_Multiplier = sSIN_Multiplier
    sINI%SIN_frac(3) = 1.d0 - sINI%SIN_frac(1) - sINI%SIN_frac(2) !!ensure total = 100%

    sINI%SOIL_INI_file = sSOIL_INI_file
    sINI%iScenario = siScenario
    sINI%STP_delta = sSTP_delta
    sINI%SWC_logis(1:3) = sSWC_logis
    sINI%SIN_logis(1:4) = sSIN_logis

    ALLOCATE(sINI%VARopt(const_nVAR0))
    ALLOCATE(sINI%VARstep(const_nVAR0))
    ALLOCATE(sINI%VARfile(const_nVAR0))
    ALLOCATE(sINI%VARcol(const_nVAR0))
    ALLOCATE(sINI%VARobj(const_nVAR0))
    ALLOCATE(sINI%VARobjw(const_nVAR0))
    
    sINI%VARopt = Cali_Calibrate
    sINI%VARstep = Cali_tstep 
    sINI%VARfile = Cali_obs_file
    sINI%VARcol = Cali_obs_file_column
    sINI%VARobj = Cali_OBJ
    sINI%VARobjw = Cali_OBJ_Weight
    sINI%VARobj_tol = Cali_OBJ_Tolerance
    
    sPAR_SCE%nRun = SCE_Parameters(1)
    sPAR_SCE%maxn = SCE_Parameters(2)
    sPAR_SCE%kstop = SCE_Parameters(3)
    sPAR_SCE%ngs = SCE_Parameters(4)
    sPAR_SCE%ideflt = SCE_Parameters(5)
    
    sPAR_SCE%pcento = 0.0001 !!SCE_Parameters(3)

    !c  IF ideflt IS EQUAL TO 1, READ THE SCE CONTROL PARAMETERS

    if (sPAR_SCE%ideflt .eq. 1) Then
        sPAR_SCE%npg = SCE_control_Parameters(1)
        sPAR_SCE%nps = SCE_control_Parameters(2)
        sPAR_SCE%nspl = SCE_control_Parameters(3)
        sPAR_SCE%mings = SCE_control_Parameters(4)
        sPAR_SCE%iniflg = SCE_control_Parameters(5)
        sPAR_SCE%iprint = SCE_control_Parameters(6)
        pcntrl = usrsp
    else
        pcntrl = deflt
    end if
!      write (*,*)sRead
!      write (*,810)sPAR_SCE%npg,sPAR_SCE%nps,sPAR_SCE%nspl,sPAR_SCE%mings,sPAR_SCE%iniflg,sPAR_SCE%iprint
  

    sINI%iKinetics = siKinetics
    sINI%iHR = siHR
    sINI%iTmp_Func = siTmp_Func
    sINI%iENZN_Allocate = siENZN_Allocate
    sINI%Ksat = Ksat
    sINI%Lambda = Lambda

    sPAR_SCE%nPar = const_nPAR0
    sINI%nPar = sPAR_SCE%nPar

    ALLOCATE(iPar_opt(sPAR_SCE%nPar))
    ALLOCATE(sPAR_SCE%parName(sPAR_SCE%nPar))
    ALLOCATE(sPAR_SCE%a(sPAR_SCE%nPar))
    ALLOCATE(sPAR_SCE%bl(sPAR_SCE%nPar))
    ALLOCATE(sPAR_SCE%bu(sPAR_SCE%nPar))
    sPAR_SCE%parName    = Pname
    sPAR_SCE%a          = Pinitial
    sPAR_SCE%bl         = Plow
    sPAR_SCE%bu         = Phigh
    iPar_opt            = Pcal 
    
    sINI%Name_PAR0 = sPAR_SCE%parName
        
    ndays = nDaysbwDates(sINI % sDate_beg_all, sINI % sDate_beg_sim)
    if (ndays .lt. 1) then
        print*, "ERROR: Date_beg_sim = ", sINI % sDate_beg_sim, " < Date_beg_all = ", sINI % sDate_beg_all
        stop
    end if

    if (sINI % iModel .eq. 1) then !!optimization
        ndays = nDaysbwDates(sINI % sDate_end_all, sINI % sDate_end_sim)
        if (ndays .gt. 0) then
            sINI % sDate_end_sim = sINI % sDate_end_all
        else
            write(*, *) "Ignore the WARNING, Simulation Period is OK (within the Input-Data Period)"
        end if
    end if
    
    
    
!!==========================================================================

      
      !!define file units
      sPAR_SCE%iFout_all = 11  !!OPT_ALL
      sPAR_SCE%iFout_end = 12  !!OPT_END
      sPAR_SCE%iFout_ini = 13  !!OPT_INI
      
      sINI%iFout_all = sPAR_SCE%iFout_all
      sINI%iFout_end = sPAR_SCE%iFout_end
      sINI%iFout_ini = sPAR_SCE%iFout_ini
      
      sINI%iFout_mcmc = 15
      
      sINI%iFout_SIM_obs    = 21
      sINI%iFout_SIM_day    = 22
      sINI%iFout_SIM_mon    = 23
      sINI%iFout_VAR_hour   = 24
      sINI%iFout_FLX_hour   = 25
      sINI%iFout_RATE_hour  = 26
      sINI%iFout_PAR_hour   = 27
      sINI%iFout_ITW_hour   = 28
      
      !!3 output files for model optimization
      !!--------------------------------------------------------------------------
      sINI%dirout = trim(sINI%dirout)//"/"//trim(sINI%SITE)//"_"
      sfilename_full = trim(sINI%dirout)//"OPT_all.out"
      open(unit=sPAR_SCE%iFout_all,file=sfilename_full,status='unknown')
      sfilename_full = trim(sINI%dirout)//"OPT_end.out"
      open(unit=sPAR_SCE%iFout_end,file=sfilename_full,status='unknown')
      sfilename_full = trim(sINI%dirout)//"OPT_ini.out"
      open(unit=sPAR_SCE%iFout_ini,file=sfilename_full,status='unknown')
      
      sfilename_full = trim(sINI%dirout)//"MCMC.out"
      open(unit=sINI%iFout_mcmc,file=sfilename_full,status='unknown')
      !!--------------------------------------------------------------------------
      
      !!output files for response variables
      !!--------------------------------------------------------------------------
      sfilename_full = trim(sINI%dirout)//'SIM_obs.out'
      open(unit = sINI%iFout_SIM_obs, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_SIM_obs,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_SIM_obs,'(2a5,a15,3a20)')"ID","VAR","Date","OBS_avg","SIM_avg","SIM_sd"
      
      sfilename_full = trim(sINI%dirout)//'SIM_day.out'
      open(unit = sINI%iFout_SIM_day, file = sfilename_full, status = 'unknown')  !!daily simulation output; see MOD_MEND::subMEND_RUN() 
      write(sINI%iFout_SIM_day,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_SIM_day,'(a10,a50)')"Day","AVG & STDDEV for VAR[1:n]"
      
      sfilename_full = trim(sINI%dirout)//'SIM_mon.out'
      open(unit = sINI%iFout_SIM_mon, file = sfilename_full, status = 'unknown')  !!daily simulation output; see MOD_MEND::subMEND_RUN() 
      write(sINI%iFout_SIM_mon,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_SIM_mon,'(a10,a50)')"Mon","AVG & STDDEV for VAR[1:n]"
      
      sfilename_full = trim(sINI%dirout)//"VAR_hour.out"
      open(unit = sINI%iFout_VAR_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_VAR_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_VAR_hour,'(a10,200a20)')"Hour",(trim(Name_POOL(i))//"_C",i=1,const_nPOOL), &
                                                (trim(Name_MNPOOL(i)),i=1,const_nPOOL_MN), &
                                                (trim(Name_POOL(i))//"_N",i=1,const_nPOOL), &
                                                (trim(Name_POOL(i))//"_CN",i=1,const_nPOOL)

      sfilename_full = trim(sINI%dirout)//"FLX_hour.out"
      open(unit = sINI%iFout_FLX_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_FLX_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_FLX_hour,'(a10,200a20)')"Hour",(trim(Name_FLUX(i))//"_C",i=1,const_nFLUX), &
                                                (trim(Name_MNFLUX(i)),i=1,const_nFLUX_MN), &
                                                (trim(Name_FLUX(i))//"_N",i=1,const_nFLUX)

      sfilename_full = trim(sINI%dirout)//"RATE_hour.out"
      open(unit = sINI%iFout_rate_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_RATE_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_RATE_hour,'(a10,100a20)')"Hour",(trim(Name_RATE(i)),i=1,const_nRATE)
      
      sfilename_full = trim(sINI%dirout)//"PAR_hour.out"
      open(unit = sINI%iFout_par_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_PAR_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_PAR_hour,'(a10,100a20)')"Hour",(trim(Name_PAR(i)),i=1,const_nPAR)
      
      sfilename_full = trim(sINI%dirout)//trim("ITW_hour.dat")
      open(unit = sINI%iFout_ITW_hour, file = sfilename_full, status = "unknown")
      write(sINI%iFout_ITW_hour,*)"Data_Period = ",sINI%sDate_beg_all, " -- ",sINI%sDate_end_all
      write(sINI%iFout_ITW_hour,'(a10,5a20)')"Hour","SIN_mg/cm3/h","STP_oC","SWC","SWP_MPa","pH"
      
      !!--------------------------------------------------------------------------
      
      write(sPAR_SCE%iFout_ini,700)
  700 format(10x,'SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION',&
     &       /,10x,46(1h=))
      
      nmons = nMonsbwDates(sINI%sDate_beg_all,sINI%sDate_end_all)
      ndays = nDaysbwDates(sINI%sDate_beg_all,sINI%sDate_end_all)
      sINI%nHour = ndays*24
      sINI%nHour_sim = 24*nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_end_sim)
      ALLOCATE(sINI%STP(sINI%nHour))
      ALLOCATE(sINI%SWC(sINI%nHour))
      ALLOCATE(sINI%SWP(sINI%nHour))
      ALLOCATE(sINI%SpH(sINI%nHour))
      ALLOCATE(sINI%SIN(sINI%nHour)) !!SOC input
        
      !!INPUT INFO:
      !!Soil temperature data
      if(ifdata_ST.eq.1) then  !!read data
          sUnits = StrCompress(sUnits_ST)
          ststep = StrCompress(step_ST)  !!TODO NEXT: convert data with time-step (ststep=monthly,daily) to hourly
          is_total = 0
          CALL sINP_Read(nfile_ST,sfilename_ST,sINI%dirinp,ststep,is_total,nmons,sINI%nHour,sINI%STP)
      else  !!constant temperature
          sINI%STP = ST_constant   
      end if
      
      !!Soil water info: Hourly or Daily
      
      sINI % porosity = vg_SWCsat      
      sINI % waterRetention_vg(1) = vg_SWCres
      sINI % waterRetention_vg(2) = vg_SWCsat
      sINI % waterRetention_vg(3) = vg_alpha
      sINI % waterRetention_vg(4) = vg_n
      
      sINI % SWCFC = fSWP2SWC(const_SWPFC,"MPa",vg_SWCres,vg_SWCsat,vg_alpha,vg_n)

      !!--------------------------------------------------
      sfilename_full = trim(sINI%dirout)//"water_retention_curve.out"
      open(unit = 101, file = sfilename_full, status = 'unknown')
      write(101,'(2A20)')"SWC","SWP_MPa"
      iRead = 800
      do k = 1,iRead
          rRead = 0.001 + (k - 1)*0.8/iRead
          rRead2 = fSWC2SWP(rRead,vg_SWCres,vg_SWCsat,vg_alpha,vg_n,const_SWPmin)
          write(101, '(f20.3,e20.6)')rRead, rRead2
      end do
      close(101)
      !!--------------------------------------------------
      
      !! soil moisture data
      if(ifdata_SM.eq.1) then  !!read hourly data         
          sUnits = StrCompress(sUnits_SM)
          ststep = StrCompress(step_SM)
          is_total = 0
          CALL sINP_Read(nfile_SM,sfilename_SM,sINI%dirinp,ststep,is_total,nmons,sINI%nHour,sINI%SWC)
          
          if (trim(sUnits).eq."percent") then  !!%, percent
            sINI%SWC = sINI%SWC/1.d2
          end if
          
          do k = 1, sINI%nHour
              sINI%SWP(k) = fSWC2SWP(sINI%SWC(k),vg_SWCres,vg_SWCsat,vg_alpha,vg_n,const_SWPmin)
          end do
          
      else  !!constant SM
          sINI%SWC = SM_constant 
          sINI%SWP = fSWC2SWP(SM_constant,vg_SWCres,vg_SWCsat,vg_alpha,vg_n,const_SWPmin)
      end if
      
      !!External input, e.g., litter fall, fertilizer: Daily or Monthly
      
      if(ifdata_type1.eq.1) then  !!nfile = 1 if ststep='monthly'
          sUnits = StrCompress(sUnits_type1)
          ststep = StrCompress(step_type1)  !!convert data with time-step (ststep=monthly,daily) to hourly
          is_total = 1  !!usually litter_input is the total amount during a period
          CALL sINP_Read(nfile_type1,sfilename_type1,sINI%dirinp,ststep,is_total,nmons,sINI%nHour,sINI%SIN)
      
      else !!constant litter input [mgC/cm3/h]
          sINI%SIN = Input_type1_constant
      end if
      
      !!soil pH
      if(ifdata_pH.eq.1) then  !!nfile = 1 if ststep='monthly'
          sUnits = StrCompress(sUnits_pH)
          ststep = StrCompress(step_pH)  !!convert data with time-step (ststep=monthly,daily) to hourly
          is_total = 0  
          CALL sINP_Read(nfile_pH,sfilename_pH,sINI%dirinp,ststep,is_total,nmons,sINI%nHour,sINI%SpH)     
      else !!constant pH
          sINI%SpH = spH_constant
      end if
      
    
      sfilename_full = trim(sINI%dirinp)//trim(sINI%SOIL_INI_file)
      call subMEND_INI_Read(sINI%dINI,sfilename_full)
      sINI % soilDepth    = sINI%dINI(1)  !![cm], soil depth
!      write(*,*)"dINI=",sINI%dINI
      
      sINI%iGPPscaler = iGPPscaler
      sINI%GPPref = GPPref
      
      write(sPAR_SCE%iFout_ini,*)
      write(sPAR_SCE%iFout_ini,*)"sINI%iModel = ", sINI%iModel
      write(sPAR_SCE%iFout_ini,*)"Carbon (C) or Carbon-Nitrogen (CN) Model = ", ssMEND
      write(sPAR_SCE%iFout_ini,*)"sINI%iKinetics = ", sINI%iKinetics
      write(sPAR_SCE%iFout_ini,*)"sINI%iHR       = ", sINI%iHR
      write(sPAR_SCE%iFout_ini,*)"sINI%iTmp_Func = ", sINI%iTmp_Func
      write(sPAR_SCE%iFout_ini,*)"sINI%iENZN_Allocate = ", sINI%iENZN_Allocate
      write(sPAR_SCE%iFout_ini,*)"sINI%iGPPscaler = ", sINI%iGPPscaler
      write(sPAR_SCE%iFout_ini,*)"sINI%GPPref = ", sINI%GPPref," mg C/cm2/d"
      write(sPAR_SCE%iFout_ini,*)
      write(sPAR_SCE%iFout_ini,'(A50,2f10.4)')"SWC_Saturation, SWC_Field_Capacity =", sINI%porosity, sINI%SWCFC
      write(sPAR_SCE%iFout_ini,'(A50,f10.4,A5)')"Saturated Hydraulic Conductivity (Ksat) =", sINI%Ksat," cm/h"
      write(sPAR_SCE%iFout_ini,'(A50,f10.4)')"Log-Slope for Hydraulic Conductivity (Lambda) =", sINI%Lambda
      write(sPAR_SCE%iFout_ini,*)
      
      sINI%GPPref = GPPref / sINI % soilDepth / 24.d0  !!covert mgC/cm2/d to mgC/cm3/h
      !!write hourly input data into 1 file
      sINI%SIN       = sINI%SIN      /sINI % soilDepth  !!covert mgC/cm2/[T] to mgC/cm3/[T]
      sINI%SIN_other = sINI%SIN_other/sINI % soilDepth
      do k = 1,ndays
          CALL sDate_After(k,sINI%sDate_beg_all,sDate)
          do j=1,24 !!sINI%nHour
              CALL sInt2Str(j,2,str2)
              i = (k - 1)*24 + j
              write(sINI%iFout_ITW_hour,'(A10,5f20.6)')sDate//str2,sINI%SIN(i),sINI%STP(i),sINI%SWC(i),sINI%SWP(i),sINI%SpH(i)
          end do
      end do
      close(sINI%iFout_ITW_hour)
      
      !!compute daily/monthly/yearly STP,SWP,SIN-----------------------------------------------------------BEG
      !    print*,">>>Inputs, Temperature, Water Content & Potential:"
      sFile_inp = trim(sINI%dirout)//"ITW_hour.dat"
      sFile_out = trim(sINI%dirout)//"ITW_day.dat"
!        sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sINI%sDate_beg_all, sINI%sDate_beg_end,tstep,flag_avg)
      CALL sOUT_ALL_tscale(sFile_inp,sFile_out,2,5, sINI%sDate_beg_all, sINI%sDate_end_all,1,1)
      CALL system('gzip -f '//sFile_inp)
      
      sFile_inp = trim(sINI%dirout)//"ITW_day.dat"
      sFile_out = trim(sINI%dirout)//"ITW_mon.dat"
      CALL sOUT_ALL_tscale(sFile_inp,sFile_out,2,5, sINI%sDate_beg_all, sINI%sDate_end_all,2,1)
      
      sFile_inp = trim(sINI%dirout)//"ITW_day.dat"
      sFile_out = trim(sINI%dirout)//"ITW_year.dat"
      CALL sOUT_ALL_tscale(sFile_inp,sFile_out,2,5, sINI%sDate_beg_all, sINI%sDate_end_all,4,1)
      !!compute daily/monthly/yearly STP,SWP,SIN-----------------------------------------------------------END  
      
      !!Calibration Variables & Data

      write(sPAR_SCE%iFout_ini,*)
      write(sPAR_SCE%iFout_ini,*)"MODEL CALIBRATION VARIABLES & DATA-FILE (see MEND_namelist.nml):"
      write(sPAR_SCE%iFout_ini,'(5a10,2a20,3a10)')Cali_var_title
      write(sPAR_SCE%iFout_ini,'(120(1h-))')

      sINI%nVARopt = 0
      do i=1,const_nVAR0
          
          if(sINI%VARopt(i).gt.0) then
              sINI%nVARopt = sINI%nVARopt + 1
              write(sPAR_SCE%iFout_ini,'(I5,5x,2a10,2(I5,5x),a20,I5,10x,a10,2f10.2)') &
                                      Cali_varid(i),StrCompress(Cali_VAR(i)),StrCompress(Cali_Units(i)), &
                                      sINI%VARopt(i),sINI%VARstep(i),StrCompress(sINI%VARfile(i)),&
                                      sINI%VARcol(i),StrCompress(sINI%VARobj(i)),sINI%VARobjw(i),sINI%VARobj_tol(i)
              sINI%Name_VARopt(sINI%nVARopt) = Cali_VAR(i)
!              sINI%Name_VARopt = ADJUSTR(sINI%Name_VARopt)
!!1 CO2 mgC-cm3-h   1 	1   HR.obs 	2   NSEC	10
          end if
      end do
      
      if(sINI%nVARopt.lt.1) then
          write(*,*)'No Variables Available for Model Optimization!!!'
      end if
         
      ALLOCATE(sINI%VARopt_int(sINI%nVARopt, 3))
      ALLOCATE(sINI%rOBJ(sINI%nVARopt))
      ALLOCATE(sINI%rOBJw(sINI%nVARopt))
      sINI%VARopt_int = 0  !!initialization

      j=0
      do i=1,const_nVAR0
          if(sINI%VARopt(i).gt.0) then
              j=j+1
              sINI%VARopt_int(j,1) = i  !!index of par_opt
          end if
      end do

      do i=1,sINI%nVARopt
          j = sINI%VARopt_int(i,1)
          sINI%rOBJw(i) = sINI%VARobjw(j)
          sINI%VARopt_int(i,3) = sINI%VARstep(j)  !!time-step
          sfilename_full = trim(sINI%dirinp)//trim(sINI%VARfile(j))
          open(101,file=sfilename_full,status='old')
          read(101,*)sRead,iRead    !!# of observations
          sINI%VARopt_int(i,2) = iRead
          close(101)
!          print*, sINI%VARopt_int(i,:)
      end do
      
      sINI%nOBS_tot = sum(sINI%VARopt_int(1:sINI%nVARopt,2))
      ALLOCATE(sINI%dOBS_opt(sINI%nOBS_tot,3)) !!date,obs,iVARopt
      ALLOCATE(sINI%dSIM_opt(sINI%nOBS_tot,3)) !!date,sim,sim_sd
      sINI%dOBS_opt = const_FillValue
      sINI%dSIM_opt = const_FillValue
      CALL Array_Normalize(sINI%nVARopt,sINI%rOBJw,const_FillValue)    
      
      !!-------------------------------------------------
 
!  800 format(2i5,f10.4,3i5)
!  810 format(6i5)

      sPAR_SCE%nOpt = 0
      do j = 1, sPAR_SCE%nPar
        sPAR_SCE%parName(j) = ADJUSTR(sPAR_SCE%parName(j))
!        print*, sPAR_SCE%parName(j)
!        print*, iPar(j), sPAR_SCE%parName(j), sPAR_SCE%a(j), sPAR_SCE%bl(j), sPAR_SCE%bu(j), iPar_opt(j)
        
        if(iPar_opt(j).gt.0) then
            sPAR_SCE%nOpt = sPAR_SCE%nOpt + 1
        end if
      end do
  830 format(3f10.3)
        
      !!---------------------------------------------------------
!      IF(.NOT.sINI%Carbon_only) THEN
          ALLOCATE(sINI%SIN_NH4(sINI%nHour))
          ALLOCATE(sINI%SIN_NO3(sINI%nHour))
          
          sINI%CN_LITT_avg = CN_ratio_input(1)
          sINI%CN_WOOD_avg = CN_ratio_input(2)
          sINI%CN_ROOT_avg = CN_ratio_input(3)
               
          !!NH4 input [mg N/cm2/mon]: Monthly
          
          if(ifdata_NH4.eq.1) then  !!nfile = 1 if ststep='monthly'
              sUnits = StrCompress(sUnits_NH4)
              ststep = StrCompress(step_NH4)  !!convert data with time-step (ststep=monthly,daily) to hourly
              is_total = 1  !!usually litter_input is the total amount during a period
              CALL sINP_Read(nfile_NH4,sfilename_NH4,sINI%dirinp,ststep,is_total,nmons,sINI%nHour,sINI%SIN_NH4)

          else !!constant input [mgN/cm2/h]
              sINI%SIN_NH4 = NH4_constant
          end if
          
          !!NO3 input [mg N/cm2/mon]: Monthly
          
          if(ifdata_NO3.eq.1) then  !!nfile = 1 if ststep='monthly'
              sUnits = StrCompress(sUnits_NO3)
              ststep = StrCompress(step_NO3)  !!convert data with time-step (ststep=monthly,daily) to hourly
              is_total = 1  !!usually litter_input is the total amount during a period
              CALL sINP_Read(nfile_NO3,sfilename_NO3,sINI%dirinp,ststep,is_total,nmons,sINI%nHour,sINI%SIN_NO3)

          else !!constant input [mgN/cm2/h]
              sINI%SIN_NO3 = NO3_constant
          end if
          
      IF(.NOT.sINI%Carbon_only) THEN    
          !!write hourly input data into 1 file
          sfilename_full = trim(sINI%dirout)//trim("Ndep_hour.dat")
          iRead = 29
          open(unit = 29, file = sfilename_full, status = "unknown")
          write(iRead,*)"Data_Period = ",sINI%sDate_beg_all, " -- ",sINI%sDate_end_all
          write(iRead,'(a10,5a20)')"Hour","SIN_NH4_mg/cm2/h","SIN_NO3_mg/cm2/h"
          do k = 1,ndays
              CALL sDate_After(k,sINI%sDate_beg_all,sDate)
              do j=1,24 !!sINI%nHour
                  CALL sInt2Str(j,2,str2)
                  i = (k - 1)*24 + j
                  write(iRead,'(A10,5f20.10)')sDate//str2,sINI%SIN_NH4(i),sINI%SIN_NO3(i)
              end do
          end do
          close(iRead)
          
          sINI%SIN_NH4 = sINI%SIN_NH4/sINI%soilDepth  !!convert [mg N/cm2/h] to [mg N/cm3/h]
          sINI%SIN_NO3 = sINI%SIN_NO3/sINI%soilDepth  !!convert [mg N/cm2/h] to [mg N/cm3/h]
!          close(iFin)
      END IF !!(.NOT.sINI%Carbon_only)
      !!---------------------------------------------------------
      
      ALLOCATE(sPAR_SCE%iOpt(sPAR_SCE%nOpt))
      i = 0
      do j = 1, sPAR_SCE%nPar
          if(iPar_opt(j).gt.0) then
              i = i + 1
              sPAR_SCE%iOpt(i) = j
          end if
      end do
!      print*, sPAR_SCE%iOpt

!c  IF ideflt IS EQUAL TO 0, SET THE SCE CONTROL PARAMETERS TO
!c  THE DEFAULT VALUES
      if (sPAR_SCE%ideflt .eq. 0) then
        sPAR_SCE%npg = 2*sPAR_SCE%nOpt + 1
        sPAR_SCE%nps = sPAR_SCE%nOpt + 1
        sPAR_SCE%nspl = sPAR_SCE%npg
        sPAR_SCE%mings = sPAR_SCE%ngs
        sPAR_SCE%iniflg = 0
        sPAR_SCE%iprint = 0
      end if


!c  CHECK IF THE SCE CONTROL PARAMETERS ARE VALID
      if (sPAR_SCE%ngs .lt. 1 .or. sPAR_SCE%ngs .ge. 1320) then
        write(sPAR_SCE%iFout_ini,900) sPAR_SCE%ngs
  900   format(//,1x,'**ERROR** NUMBER OF COMPLEXES IN INITIAL ',&
     &         ' POPULATION ',i5,' IS NOT A VALID CHOICE')
        ierror = ierror + 1
      end if

      if (sPAR_SCE%kstop .lt. 0 .or. sPAR_SCE%kstop .ge. 20) then
        write(sPAR_SCE%iFout_ini,901) sPAR_SCE%kstop
  901   format(//,1x,'**WARNING** THE NUMBER OF SHUFFLING LOOPS IN',&
     &  ' WHICH THE CRITERION VALUE MUST CHANGE ',/,13x,'SHOULD BE',&
     &  ' GREATER THAN 0 AND LESS THAN 10.  ','kstop = ',i2,&
     &  ' WAS SPECIFIED.'/,13x,'BUT kstop = 5 WILL BE USED INSTEAD.')
        iwarn = iwarn + 1
        sPAR_SCE%kstop=5
      end if

      if (sPAR_SCE%mings .lt. 1 .or. sPAR_SCE%mings .gt. sPAR_SCE%ngs) then
        write(sPAR_SCE%iFout_ini,902) sPAR_SCE%mings
  902   format(//,1x,'**WARNING** THE MINIMUM NUMBER OF COMPLEXES ',&
     &         i2,' IS NOT A VALID CHOICE. SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%mings = sPAR_SCE%ngs
      end if

      if (sPAR_SCE%npg .lt. 2 .or. sPAR_SCE%npg .gt. 1320/max(sPAR_SCE%ngs,1)) then
        write(sPAR_SCE%iFout_ini,903) sPAR_SCE%npg
  903   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A COMPLEX ',&
     &         I4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%npg = 2*sPAR_SCE%nOpt+1
      end if

      if (sPAR_SCE%nps.lt.2 .or. sPAR_SCE%nps.gt.sPAR_SCE%npg .or. sPAR_SCE%nps.gt.50) then
        write(sPAR_SCE%iFout_ini,904) sPAR_SCE%nps
  904   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A SUB-',&
     &  'COMPLEX ',i4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%nps = sPAR_SCE%nOpt + 1
      end if

      if (sPAR_SCE%nspl .lt. 1) then
        write(sPAR_SCE%iFout_ini,905) sPAR_SCE%nspl
  905   format(//,1x,'**WARNING** THE NUMBER OF EVOLUTION STEPS ',&
     &         'TAKEN IN EACH COMPLEX BEFORE SHUFFLING ',I4,/,13x,&
     &         'IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%nspl = sPAR_SCE%npg
      end if

!c  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPULATION
      sPAR_SCE%npt = sPAR_SCE%ngs * sPAR_SCE%npg

      if (sPAR_SCE%npt .gt. 1320) then
        write(sPAR_SCE%iFout_ini,906) sPAR_SCE%npt
  906   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN INITIAL ',&
     &         'POPULATION ',i5,' EXCEED THE POPULATION LIMIT,',/,13x,&
     &         'SET NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS')
        iwarn = iwarn + 1
        sPAR_SCE%ngs = 2
        sPAR_SCE%npg = 2*sPAR_SCE%nOpt + 1
        sPAR_SCE%nps = sPAR_SCE%nOpt + 1
        sPAR_SCE%nspl = sPAR_SCE%npg
      end if

!c  PRINT OUT THE TOTAL NUMBER OF ERROR AND WARNING MESSAGES
      if (ierror .ge. 1) write(sPAR_SCE%iFout_ini,907) ierror
  907 format(//,1x,'*** TOTAL NUMBER OF ERROR MESSAGES IS ',i2)

      if (iwarn .ge. 1) write(sPAR_SCE%iFout_ini,908) iwarn
  908 format(//,1x,'*** TOTAL NUMBER OF WARNING MESSAGES IS ',i2)

      if (sPAR_SCE%mings .lt. sPAR_SCE%ngs) then
        reduc = ysflg
      else
        reduc = noflg
      end if

      if (sPAR_SCE%iniflg .ne. 0) then
        initl = ysflg
      else
        initl = noflg
      end if


!c  PRINT SHUFFLED COMPLEX EVOLUTION OPTIMIZATION OPTIONS
  104 write(sPAR_SCE%iFout_ini,910)
  910 format(//,2x,'SCE CONTROL',5x,'MAX TRIALS',5x,&
     &'REQUIRED IMPROVEMENT',5x,'RANDOM',/,3x,'PARAMETER',8x,&
     &'ALLOWED',6x,'PERCENT',4x,'NO. LOOPS',6x,'SEED',/,&
     &2x,11(1h-),5x,10(1H-),5x,7(1h-),4x,9(1h-),5x,6(1h-))

      pcenta=sPAR_SCE%pcento*100.
      write(sPAR_SCE%iFout_ini,912) pcntrl,sPAR_SCE%maxn,pcenta,sPAR_SCE%kstop,sPAR_SCE%iseed
  912 format(3x,a10,7x,i5,7x,f6.2,9x,i2,9x,i5)
      write(sPAR_SCE%iFout_ini,914) sPAR_SCE%ngs,sPAR_SCE%npg,sPAR_SCE%npt,sPAR_SCE%nps,sPAR_SCE%nspl
  914 format(//,18x,'SCE ALGORITHM CONTROL PARAMETERS',/,18x,32(1H=),&
     &//,2x,'NUMBER OF',5x,'POINTS PER',5x,'POINTS IN',6x,'POINTS PER',&
     &4x,'EVOL. STEPS',/,2x,'COMPLEXES',6X,'COMPLEX',6x,'INI. POPUL.',&
     &5x,'SUB-COMPLX',4x,'PER COMPLEX',/,2x,9(1h-),5x,10(1h-),4x,&
     &11(1h-),5x,10(1h-),4x,11(1h-),5x,/,2x,5(i5,10x))
      write(sPAR_SCE%iFout_ini,915) reduc,sPAR_SCE%mings,initl
  915 format(//,15x,'COMPLX NO.',5x,'MIN COMPLEX',5x,'INI. POINT',/,&
     &15x,'REDUCTION',6x,'NO. ALLOWED',6x,'INCLUDED',/,&
     &15x,10(1h-),5x,11(1h-),5x,10(1h-),/,18x,a4,6x,i8,13x,a4)
      write(sPAR_SCE%iFout_ini,916)
  916 format(//,8x,'INITIAL PARAMETER VALUES AND PARAMETER BOUNDS',/,&
     &       8x,45(1h=),//,2x,'ID',2x,'PARAMETER',6x,'OPT-Y/N',2x,'INITIAL VALUE',5x,&
     &       'LOWER BOUND',5x,'UPPER BOUND',/,15(1h-),3x,10(1h-),2x,&
     &       13(1h-),5x,11(1h-),5x,11(1h-))
     
      do 920 i = 1, sPAR_SCE%nPar
        write(sPAR_SCE%iFout_ini,918) i, sPAR_SCE%parName(i),iPar_opt(i),sPAR_SCE%a(i),sPAR_SCE%bl(i),sPAR_SCE%bu(i)
!        write(*,918) sPAR_SCE%parName(i),sPAR_SCE%a(i),sPAR_SCE%bl(i),sPAR_SCE%bu(i),iPar_opt(i)
  920 continue
  918   format(I4, a10,4x,I10,(3x,f12.6),2(4x,f12.6))
!  918   format(a10,4x,(3x,f12.6),2(4x,f12.6),5x,I10)
  
      if (ierror .ge. 1) then
      write(sPAR_SCE%iFout_ini,922)
  922 format(//,'*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE',&
     &       ' OF INPUT DATA ERROR ***')
      stop
      end if
      
      
      write (*,*) '>>>EXIT SUBROUTINE <MEND_INI>'

      return
END !!SUBROUTINE MENDIN

