MODULE MOD_USRFS
    !!A Collection of User-Defined Functions & Subroutines
    !!Author: Wang, Gangsheng
    !!ORNL
    !!wangg@ornl.gov
    !!updated: May 5, 2015
    !------------------------------------------------------------------
    !! [1] OBJECTIVE FUNCTION
    !! [2] ARRAY FUNCTION
    !! [3] DATETIME FUNCTION
    !------------------------------------------------------------------
    IMPLICIT NONE
    REAL(8), PARAMETER:: dLOW  = 1.0D-30
    REAL(8), PARAMETER:: dHIGH = 1.0D+30
    
    ! ----------
    ! Visibility
    ! ----------
    ! Everything private by default
    PRIVATE
    ! Public procedures
    PUBLIC :: fWAVG
    PUBLIC :: WOBJ_INI
    PUBLIC :: fNSE
    PUBLIC :: f1NSE
    PUBLIC :: fMARE
    PUBLIC :: fMARE_tolerance
    PUBLIC :: f1NSE_norm
    PUBLIC :: fMARE_norm
    PUBLIC :: fMAP
    PUBLIC :: fSSE  !! sum of squared errors
    PUBLIC :: fSSE_Normalize  !! SSE normalized by variance of observations
    PUBLIC :: fMSE_Normalize  !! Mean Squared Errors
    PUBLIC :: fNRMSE
    PUBLIC :: fSUM
    PUBLIC :: fAVG
    PUBLIC :: f1RAVG !! =1-fAVG_sim/fAVG_obs
    PUBLIC :: f1RAVG_ratio  !!dabs(ratio-fAVG_sim/fAVG_obs)
    PUBLIC :: fSUM2 !! sum of array elements between iBeg and iEnd
    PUBLIC :: fAVG2 !! avg of array elements between iBeg and iEnd
    PUBLIC :: fLxy
    PUBLIC :: fCORR !!correlation coefficient
    PUBLIC :: f1CORR !! = 1-fCORR
    PUBLIC :: fSTDDEV !!standard deviation
    PUBLIC :: fVARIANCE !!=fSTDDEV
    PUBLIC :: fSIGN    !! = 1 for the same sign (positive or negative), =0 for opposite sign

    PUBLIC :: fMAXarray !!max element in an array
    PUBLIC :: fMINarray
    PUBLIC :: fRANGEarray
    PUBLIC :: Array_Normalize
    PUBLIC :: selectINT
    PUBLIC :: sort
    PUBLIC :: sort1
    PUBLIC :: indexx

    PUBLIC :: nDaysUPMon
    PUBLIC :: nDaysofMon
    PUBLIC :: nDaysofYears
    PUBLIC :: nDaysofYear
    PUBLIC :: nMonsbwDates
    PUBLIC :: nYearsbwDates
    PUBLIC :: nMonths
    PUBLIC :: nDaysbwDates
    PUBLIC :: iJulianDay
    PUBLIC :: sDate_After
    PUBLIC :: sYM_After
    PUBLIC :: sDate2YMD
    PUBLIC :: sYMD2Date
    PUBLIC :: sInt2Str
    PUBLIC :: Sec2HMS
    PUBLIC :: iRandSeedGen
    
    PUBLIC :: ran1
    PUBLIC :: gasdev0
    PUBLIC :: gasdev        !!Normal distribution random number
    
!    PUBLIC :: fPa

CONTAINS

    !!=================================================================
    !! [1] OBJECTIVE FUNCTION: BEGIN
    !------------------------------------------------------------------

    REAL(8) FUNCTION fWAVG(n, ww, xx, dFill)
        !!compute weighting average
        !!ARGUMENTS:
        INTEGER n
        REAL(8), DIMENSION(n) :: ww, xx
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) sum1      
        INTEGER j, k
        INTEGER jj(n)
        REAL(8), DIMENSION(:), ALLOCATABLE :: ww2
        k = 0
        do j = 1, n
            if (xx(j) .ne. dFill) then !!has a value different than dFill (e.g., -9999)
                k = k + 1
                jj(k) = j
            end if
        end do

        if (k < 1) then
            fWAVG = dFill
        else
            ALLOCATE(ww2(k))
            do j = 1, k
                ww2(j) = ww(jj(j))
            end do
            CALL Array_Normalize(k, ww2, dFill)
            fWAVG = 0
            do j = 1, k
                fWAVG = fWAVG + ww2(j) * xx(jj(j))
            end do
        end if
        return
    END !!FUNCTION
    !------------------------------------------------------------------
    SUBROUTINE WOBJ_INI(nobj, wobj, nobj0, wobj0)
        !!assign value to wobj based on wobj0
        !!ARGUMENTS:
        INTEGER nobj, nobj0
        REAL(8), DIMENSION(nobj) :: wobj
        REAL(8), DIMENSION(nobj0) :: wobj0
        
        if (nobj .lt. 1) then
            return
        elseif (nobj .le. nobj0) then
            wobj(1:nobj) = wobj0(1:nobj)
        else
            wobj(1:nobj0) = wobj0(1:nobj0)
            wobj(nobj0 + 1:nobj) = wobj0(nobj0)
        end if
        return
    END !!SUBROUTINE !!WOBJ_INI
    !------------------------------------------------------------------
    REAL(8) function fNSE(nArray, Dataobs, Datasim, dFill)
        !Nash-Sutcliffe Efficiency Coefficient; Determination Coefficient
        !data excluded from calculation if Dataobs(i) = const_FillValue 
        !    REAL(8), parameter:: const_FillValue = -999d0
        !!ARGUMENTS:
        INTEGER, intent(in) :: nArray
        REAL(8) dFill
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        
        !!LOCAL VARIABLES:
        REAL(8) AVGobs, sum11, sum12
        INTEGER i, j, m
        INTEGER iData(nArray)

        m = 0
        sum11 = 0
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
                sum11 = sum11 + Dataobs(i)
            end if
        end do
        if (m .lt. 2) then
            fNSE = dFill
            return
        end if
        AVGobs = sum11/DBLE(m) !SUM(Dataobs)/nArray

        sum11 = 0
        sum12 = 0
        do i = 1, m
            j = iData(i)
            sum11 = sum11 + (Datasim(j) - Dataobs(j))*(Datasim(j) - Dataobs(j))
            sum12 = sum12 + (Dataobs(j) - AVGobs)*(Dataobs(j) - AVGobs)
        end do

        fNSE = 1 - sum11/sum12
        return
    END !!FUNCTION  
    !------------------------------------------------------------------
    REAL(8) function f1NSE(nArray, Dataobs, Datasim, dFill)
        !f1NSE = 1.0 - fNSE
        !NSE:Nash-Sutcliffe Efficiency Coefficient; Determination Coefficient
        !data excluded from calculation if Dataobs(i) = const_FillValue 
        !    REAL(8), parameter:: const_FillValue = -999d0
        !!ARGUMENTS:
        INTEGER, intent(in) :: nArray
        REAL(8) dFill
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        
        !!LOCAL VARIABLES:
        REAL(8) AVGobs, sum11, sum12
        INTEGER i, j, m
        INTEGER iData(nArray)

        m = 0
        sum11 = 0
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
                sum11 = sum11 + Dataobs(i)
            end if
        end do
        if (m .lt. 2) then
            f1NSE = dFill
            return
        end if
        AVGobs = sum11/REAL(m) !SUM(Dataobs)/nArray

        sum11 = 0
        sum12 = 0
        do i = 1, m
            j = iData(i)
            sum11 = sum11 + (Datasim(j) - Dataobs(j))*(Datasim(j) - Dataobs(j))
            sum12 = sum12 + (Dataobs(j) - AVGobs)*(Dataobs(j) - AVGobs)
        end do

        f1NSE = sum11/sum12
        return
    END !!FUNCTION 
    !------------------------------------------------------------------
    REAL(8) function f1NSE_norm(nArray, Dataobs, Datasim, dFill)
        !f1NSE = 1.0 - fNSE
        !Normalize Dataobs and Datasim, respectively, then calculate f1NSE
        !NSE:Nash-Sutcliffe Efficiency Coefficient; Determination Coefficient
        !data excluded from calculation if Dataobs(i) = const_FillValue 
        !    REAL(8), parameter:: const_FillValue = -999d0
        !!ARGUMENTS:
        INTEGER, intent(in) :: nArray
        REAL(8) dFill
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        
        
        !!LOCAL VARIABLES:
        REAL(8) AVGobs, sum11, sum12
        REAL(8), ALLOCATABLE :: obs(:), sim(:)
        REAL(8) MINobs, MAXobs, MINsim, MAXsim

        INTEGER i, j, m
        INTEGER iData(nArray)

        m = 0
        
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
            end if
        end do
        
        if (m .lt. 2) then
            f1NSE_norm = dFill
            return
        end if
        
        ALLOCATE(obs(m))
        ALLOCATE(sim(m))
        obs = Dataobs(iData)
        sim = Datasim(iData)
        
        MINobs = minval(obs)
        MAXobs = maxval(obs)
        MINsim = minval(sim)
        MAXsim = maxval(sim)
        
        !! Normalize data
!        obs = (obs-MINobs)/(MAXobs - MINobs)
!        sim = (sim-MINsim)/(MAXsim - MINsim)
        obs = obs/MAXobs
        sim = sim/MAXsim        
        
        AVGobs = fAVG(m, obs, dFill)

        sum11 = 0
        sum12 = 0
        do j = 1, m
            sum11 = sum11 + (sim(j) - obs(j))*(sim(j) - obs(j))
            sum12 = sum12 + (obs(j) - AVGobs)*(obs(j) - AVGobs)
        end do

        f1NSE_norm = sum11/sum12
        return
    END !!FUNCTION
    
    !------------------------------------------------------------------
    REAL(8) function fMARE_norm(nArray, Dataobs, Datasim, dFill)
        !!Mean Absolute Relative Error
        !!Normalize Dataobs and Datasim, respectively, then calculate MARE
        !	implicit REAL*8 (a-h,o-z)
        !	dimension Dataobs(1000),Datasim(1000)
        !!ARGUMENTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) sum1
        INTEGER i, j, m, iData(nArray)
        
        REAL(8), ALLOCATABLE :: obs(:), sim(:)
        REAL(8) MINobs, MAXobs, MINsim, MAXsim

        m = 0
        
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
            end if
        end do
        
        if (m .lt. 1) then
            fMARE_norm = dFill
            return
        end if
        
        ALLOCATE(obs(m))
        ALLOCATE(sim(m))
        obs = Dataobs(iData)
        sim = Datasim(iData)
        
        MINobs = minval(obs)
        MAXobs = maxval(obs)
        MINsim = minval(sim)
        MAXsim = maxval(sim)
        
        !! Normalize data
!        obs = (obs-MINobs)/(MAXobs - MINobs)
!        sim = (sim-MINsim)/(MAXsim - MINsim)
        obs = obs/MAXobs
        sim = sim/MAXsim 
        
        sum1 = 0
        do j = 1, m
            sum1 = sum1 + DABS(sim(j) - obs(j))/max(DABS(obs(j)), dLOW)
        end do
        fMARE_norm = sum1/REAL(m)
        return
    END !!FUNCTION
    
    !------------------------------------------------------------------
    REAL(8) function fMARE(nArray, Dataobs, Datasim, dFill)
        !!Mean Absolute Relative Error
        !	implicit REAL*8 (a-h,o-z)
        !	dimension Dataobs(1000),Datasim(1000)
        !!ARGUMENTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) sum1
        INTEGER i, j, m, iData(nArray)

        m = 0
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
            end if
        end do
        if (m .lt. 1) then
            fMARE = dFill
            return
        end if
        sum1 = 0
        do i = 1, m
            j = iData(i)
            sum1 = sum1 + DABS(Datasim(j) - Dataobs(j))/max(DABS(Dataobs(j)), dLOW)
        end do
        fMARE = sum1/REAL(m)
        return
    END !!FUNCTION
    
    REAL(8) function fMARE_tolerance(nArray, Dataobs, Datasim, tolerance, dFill)
        !!Mean Absolute Relative Error < tolerance
        !	implicit REAL*8 (a-h,o-z)
        !	dimension Dataobs(1000),Datasim(1000)
        !!ARGUMENTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        REAL(8) :: tolerance, dFill
        
        !!LOCAL VARIABLES:
        REAL(8) fRE
        
        fRE = fMARE(nArray, Dataobs, Datasim, dFill)
        
        if (fRE .eq. dFill) then
            fMARE_tolerance = dFill
            return
        end if
        
        if(DABS(fRE).le.DABS(tolerance)) then
            fMARE_tolerance = 0.0D0
        else
            fMARE_tolerance = DABS(fRE)
        end if
        
        return
    END !!FUNCTION
    
    !------------------------------------------------------------------
    REAL(8) function fSIGN(nArray, Dataobs, Datasim, dSign, dFill)
    !! = 1 for the same sign (positive or negative), =0 for opposite sign
    !! dSign = +1 means ySim should be > yBase, , define by Cali_OBJ_Tolerance in namelist
    !! dSign = -1 means ySim should be < yBase
        !!ARGUMEMNTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        REAL(8) :: dSign, dFill

        !!LOCAL VARIABLES:
        REAL(8) AVGobs, AVGsim, dtemp

        AVGobs = fAVG(nArray, Dataobs, dFill)
        AVGsim = fAVG(nArray, Datasim, dFill)

        dtemp = (AVGsim - AVGobs)*dSign
        if(dtemp.gt.0) then
            fSign = 0  !! satisfy the objective
        else
            fSign = 1  !! fail to satisfy the objective
        end if

        return

    END
    !------------------------------------------------------------------
    REAL(8) function fMAP(nArray, Dataobs, Datasim, dFillValue)
        !Mean Absolute Percent Error (MAP), best MAP = 0, [Kothamasu et al., 2004]
        !data excluded from calculation if Dataobs(i) = const_FillValue 
        !    REAL(8), parameter:: const_FillValue = -999d0
        !!ARGUMENTS:
        INTEGER, intent(in) :: nArray
        REAL(8), intent(in) :: Dataobs(nArray), Datasim(nArray), dFillValue
        
        !!LOCAL VARIABLES:
        REAL(8) sum11
        INTEGER i, m, iData(nArray)

        m = 0
        do i = 1, nArray
            if (Dataobs(i) .ne. dFillValue.and.Datasim(i) .ne. dFillValue) then
                m = m + 1
                iData(m) = i
                !            sum11 = sum11 + Dataobs(i)
            end if
        end do
        if (m .lt. 1) then
            fMAP = dFillValue
            return
        end if
        !    AVGobs = sum11/m !SUM(Dataobs)/nArray

        sum11 = 0

        do i = 1, m
            sum11 = sum11 + abs((Datasim(iData(i)) - Dataobs(idata(i)))/Dataobs(iData(i)))
        end do

        fMAP = sum11/m
        return
    END !!FUNCTION 

    !------------------------------------------------------------------
    REAL(8) function fSSE(nArray, Dataobs, Datasim, ifNorm, dFill)
        !!Sum of Squared Errors
        !!ifNorm = 1: Normalize Dataobs and Datasim by their mean value, respectively, 
        !	implicit REAL*8 (a-h,o-z)
        !	dimension Dataobs(1000),Datasim(1000)
        !!ARGUMENTS:
        INTEGER nArray,ifNorm
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) sum1
        INTEGER i, j, m, iData(nArray)
        
        REAL(8), ALLOCATABLE :: obs(:), sim(:)
        REAL(8) AVGobs, AVGsim

        m = 0
        
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
            end if
        end do
        
        if (m .lt. 1) then
            fSSE = dFill
            return
        end if
        
        ALLOCATE(obs(m))
        ALLOCATE(sim(m))
        obs = Dataobs(iData)
        sim = Datasim(iData)
        
        if(ifNorm.eq.1) then
            AVGobs = fAVG(m,obs,dFill)
            AVGsim = fAVG(m,sim,dFill)
            obs = obs/max(1.0d-20,AVGobs)
            sim = sim/max(1.0d-20,AVGsim)
        end if
  
        !! Normalize data
!        obs = (obs-MINobs)/(MAXobs - MINobs)
!        sim = (sim-MINsim)/(MAXsim - MINsim)
!        obs = obs/MAXobs
!        sim = sim/MAXsim 
        
        sum1 = 0
        do j = 1, m
            sum1 = sum1 + (sim(j) - obs(j))*(sim(j) - obs(j))
        end do
        fSSE = sum1
        return
    END !!FUNCTION 
    
        !------------------------------------------------------------------
    REAL(8) function fSSE_Normalize(nArray, Dataobs, Datasim, ifNorm, dFill)
        !!SSE/[2*VAR)
        !!SSE: Sum of Squared Errors 
        !!ifNorm = 1: Normalize Dataobs and Datasim by their mean value, respectively, 
        !	implicit REAL*8 (a-h,o-z)
        !	dimension Dataobs(1000),Datasim(1000)
        !!ARGUMENTS:
        INTEGER nArray,ifNorm
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) sum1, var1
        INTEGER i, j, m, iData(nArray)
        
        REAL(8), ALLOCATABLE :: obs(:), sim(:)
        REAL(8) AVGobs, AVGsim

        m = 0
        
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
            end if
        end do
        
        if (m .lt. 1) then
            fSSE_Normalize = dFill
            return
        end if
        
        ALLOCATE(obs(m))
        ALLOCATE(sim(m))
        obs = Dataobs(iData)
        sim = Datasim(iData)
        
        if(ifNorm.eq.1) then
            AVGobs = fAVG(m,obs,dFill)
            AVGsim = fAVG(m,sim,dFill)
            obs = obs/max(1.0d-20,AVGobs)
            sim = sim/max(1.0d-20,AVGsim)
        end if
        
        AVGobs = fAVG(m,obs,dFill)
        AVGsim = fAVG(m,sim,dFill)
        
        var1 = fVARIANCE(m,obs,dFill)
        if(var1.lt.1.0d-20) var1 = AVGobs*AVGobs 
        
        sum1 = 0.d0
        do j = 1, m
            sum1 = sum1 + (sim(j) - obs(j))*(sim(j) - obs(j))
        end do
!        write(*,*)'SSE=',m,sum1,AVGobs,var1
        fSSE_Normalize = sum1/(2.d0*var1)
        return
    END !!FUNCTION 
    
         !------------------------------------------------------------------
    REAL(8) function fMSE_Normalize(nArray, Dataobs, Datasim, ifNorm, dFill)
        !!1/m*SSE/[2*VAR)
        !!MSE: Mean Squared Errors 
        !!ifNorm = 1: Normalize Dataobs and Datasim by their mean value, respectively, 
        !	implicit REAL*8 (a-h,o-z)
        !	dimension Dataobs(1000),Datasim(1000)
        !!ARGUMENTS:
        INTEGER nArray,ifNorm
        REAL(8), DIMENSION(nArray) :: Dataobs, Datasim
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) sum1, var1
        INTEGER i, j, m, iData(nArray)
        
        REAL(8), ALLOCATABLE :: obs(:), sim(:)
        REAL(8) AVGobs, AVGsim

        m = 0
        
        do i = 1, nArray
            if (Dataobs(i) .ne. dFill.and.Datasim(i) .ne. dFill) then
                m = m + 1
                iData(m) = i
            end if
        end do
        
        if (m .lt. 1) then
            fMSE_Normalize = dFill
            return
        end if
        
        ALLOCATE(obs(m))
        ALLOCATE(sim(m))
        obs = Dataobs(iData)
        sim = Datasim(iData)
        
        if(ifNorm.eq.1) then
            AVGobs = fAVG(m,obs,dFill)
            AVGsim = fAVG(m,sim,dFill)
            obs = obs/max(1.0d-20,AVGobs)
            sim = sim/max(1.0d-20,AVGsim)
        end if
        
        var1 = fVARIANCE(m,obs,dFill)
        
        sum1 = 0
        do j = 1, m
            sum1 = sum1 + (sim(j) - obs(j))*(sim(j) - obs(j))
        end do
        fMSE_Normalize = sum1/2*var1/dble(m)
        return
    END !!FUNCTION 
    !------------------------------------------------------------------
    REAL(8) FUNCTION fNRMSE(nn, dobs, dobs_sd, dsim, dFill)
        !!Normalized Root-Mean-Square Error
        !!ARGUMENTS:      
        INTEGER nn
        REAL(8) dFill
        REAL(8), DIMENSION(nn) :: dobs, dobs_sd, dsim
        
        !!LOCAL VARIABLES:
        INTEGER j,k
        REAL(8) sum1, dobs_range
        REAL(8), PARAMETER :: dtol = 1.0d-6
        

        dobs_range = fRANGEarray(nn, dobs, dFill)

        sum1 = 0
        k = 0
        DO j = 1, nn
            if (dobs_sd(j) .eq. dFill.or.dobs_sd(j) .le. dtol) then
                dobs_sd(j) = dobs_range
            end if

            if (dobs(j) .ne. dFill.and.dsim(j).ne.dFill) then
                if (dobs_sd(j) .ne. dFill) then
                    sum1 = sum1 + &
                    (dsim(j) - dobs(j))*(dsim(j) - dobs(j)) &
                    /(dobs_sd(j) * dobs_sd(j))
                    k = k + 1
                end if
            end if
            if (k .gt. 0) then
                fNRMSE = sqrt(sum1/k)
            else
                fNRMSE = dFill
            end if
        end do
        return
    END !!FUNCTION

    !------------------------------------------------------------------
    REAL(8) function fSUM(nArray, xx, dFill)
        !!compute SUM of all elements in an Array. 
        !        implicit REAL(8)*8 (a-h,o-z)
        !!ARGUMENTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: xx
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        INTEGER jua
        
        fSUM = 0.0
        do jua = 1, nArray
            if (xx(jua) .ne. dFill) then
                fSUM = fSUM + xx(jua) !!min(dHIGH,max(dLOW,xx(jua)))
            end if
        end do
        return
    END !!FUNCTION
    !------------------------------------------------------------------
    REAL(8) function fAVG(nArray, xx, dFill)
        !!compute SUM of all elements in an Array. 
        !        implicit REAL*8 (a-h,o-z)
        !!ARGUMENTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: xx
        REAL(8) dFill,dx
        
        !!LOCAL VARIABLES:
        INTEGER jua,k
        
        fAVG = 0.0
        k = 0
        do jua = 1, nArray
            if (xx(jua) .ne. dFill) then
                k = k + 1
                if(abs(xx(jua)).ge.dHIGH.or.abs(xx(jua)).le.dLOW) then  !!wgs [12/12/2018]: avoid big/small values
                    dx = 0.d0 
                else
                    dx = xx(jua)
                end if
                fAVG = fAVG + dx !!min(dHIGH,max(dLOW,xx(jua)))
            end if
        end do
        if (k < 1) then
            fAVG = dFill
        else
            fAVG = fAVG/k
        end if
        return
    END
    !------------------------------------------------------------------
    REAL(8) function f1RAVG(nn, dobs, dsim, dFill)
        !!f1RAVG = abs(1.0 - fAVG_sim/fAVG_obs) 
        !        implicit REAL*8 (a-h,o-z)
        !!ARGUMENTS:
        INTEGER nn
        REAL(8), DIMENSION(nn) :: dobs, dsim
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) fAVG_obs, fAVG_sim

        fAVG_obs = fAVG(nn, dobs, dFill)
        fAVG_sim = fAVG(nn, dsim, dFill)
        if (fAVG_obs .eq. 0.or.fAVG_obs .eq. dFill.or.fAVG_sim .eq. dFill) then
            f1RAVG = dFill
            return
        else
            f1RAVG = dabs(1.0d0 - fAVG_sim/fAVG_obs)
            return
        end if
    END
    !------------------------------------------------------------------
    REAL(8) function f1RAVG_ratio(nn, dobs, dsim, ratio,dFill)
        !!f1RAVG = abs(1.0 - fAVG_sim/fAVG_obs) 
        !        implicit REAL*8 (a-h,o-z)
        !!ARGUMENTS:
        INTEGER nn
        REAL(8), DIMENSION(nn) :: dobs, dsim
        REAL(8) dFill, ratio
        
        !!LOCAL VARIABLES:
        REAL(8) fAVG_obs, fAVG_sim

        fAVG_obs = fAVG(nn, dobs, dFill)
        fAVG_sim = fAVG(nn, dsim, dFill)
        if (fAVG_obs .eq. 0.or.fAVG_obs .eq. dFill.or.fAVG_sim .eq. dFill) then
            f1RAVG_ratio = dFill
            return
        else
            f1RAVG_ratio = dabs(ratio - fAVG_sim/fAVG_obs)
            return
        end if
    END
    !------------------------------------------------------------------	
    REAL(8) function fSUM2(nArray, xx, iBegin, iEnd, dFill)
        !!compute AVERAGE of selected elements (iBegin to iEnd) in an Array. 
        !        implicit REAL*8 (a-h,o-z)
        !!ARGUMENTS:
        INTEGER nArray, iBegin, iEnd
        REAL(8), DIMENSION(nArray) :: xx
        REAL(8) dFill
        
        !1LOCAL VARIABLES:
        INTEGER jua

        fSUM2 = 0.0
        do jua = iBegin, iEnd
            if (xx(jua) .ne. dFill) then
                fSUM2 = fSUM2 + xx(jua) !!min(dHIGH,max(dLOW,xx(jua)))
            end if
        end do
        return
    END
    !------------------------------------------------------------------	
    REAL(8) function fAVG2(nArray, xx, iBegin, iEnd, dFill)
        !!compute AVERAGE of selected elements (iBegin to iEnd) in an Array. 
        !        implicit REAL*8 (a-h,o-z)
        !!ARGUMENTS:
        INTEGER nArray, iBegin, iEnd
        REAL(8), DIMENSION(nArray) :: xx
        REAL(8) dFill, dx
        
        !!LOCAL VARIABLES:
        INTEGER jua, k
        
        fAVG2 = 0.0
        k = 0
        do jua = iBegin, iEnd
            if (xx(jua) .ne. dFill.and.(.not.(isnan(xx(jua))))) then
                k = k + 1
                if(abs(xx(jua)).ge.dHIGH.or.abs(xx(jua)).le.dLOW) then  !!wgs [12/12/2018]: avoid big/small values
                    dx = 0.d0
                else
                    dx = xx(jua)
                end if
                fAVG2 = fAVG2 + dx !!min(dHIGH,max(dLOW,xx(jua)))
            end if
        end do
        if (k .gt. 0) then
            fAVG2 = fAVG2/k !!(iEnd - iBegin + 1)
        else
            fAVG2 = dFill
        end if
        return
    END
    !------------------------------------------------------------------

    REAL(8) function fLxy(nArray, xx, yy, dFill)
        !!ARGUMENTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: xx, yy
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8), DIMENSION(nArray) :: xx1, yy1
        INTEGER jua, k   
        REAL(8) sumxy, sumx, sumy
        !        REAL(8) fSUM2 !!FUNCTION

        sumxy = 0.0
        k = 0
        do jua = 1, nArray
            if (xx(jua) .ne. dFill.and.yy(jua) .ne. dFill) then
                k = k + 1
                xx1(k) = xx(jua)
                yy1(k) = yy(jua)
                sumxy = sumxy + xx(jua) * yy(jua)
            end if
        end do
        if (k < 1) then
            fLxy = dFill
        else
            sumx = fSUM2(nArray, xx1, 1, k, dFill)!!fSUM(nArray,xx,dFill)
            sumy = fSUM2(nArray, yy1, 1, k, dFill)!!fSUM(nArray,yy,dFill)
            fLxy = sumxy - sumx * sumy/dble(k) !!FuncLxy = sumxy - sumx*sumy/nArray
            
!            write(*,*)"xx=",xx1
!            write(*,*)"yy=",yy1
!            write(*,*)"fLxy:",sumxy, sumx, sumy,sumx * sumy/dble(k),fLxy
        end if
        return
    END
    !------------------------------------------------------------------
    REAL(8) function fCORR(nArray, xx, yy, dFill)
        !!CORRELATION COEFFICIENT
        !!ARGUMENTS:
        INTEGER nArray
        REAL(8), DIMENSION(nArray) :: xx, yy
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) aLxx, aLyy, aLxy

        aLxx = fLxy(nArray, xx, xx, dFill)
        aLyy = fLxy(nArray, yy, yy, dFill)
        aLxy = fLxy(nArray, xx, yy, dFill)
        if (aLxx .eq. dFill.or.aLyy .eq. dFill) then
            fCORR = dFill
        else
            fCORR = aLxy/SQRT(aLxx * aLyy)
        end if
        return
    END
    !------------------------------------------------------------------
    REAL(8) function f1CORR(nArray, xx0, yy0, iLog,dFill)
        !!1-CORRELATION COEFFICIENT
        !!ARGUMENTS:
        !! iLog = 1: log10 transformation
        INTEGER nArray, iLog
        REAL(8), DIMENSION(nArray) :: xx0, yy0
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8), DIMENSION(nArray) :: xx,yy
        REAL(8) aLxx, aLyy, aLxy
        
        

        aLxx = fLxy(nArray, xx0, xx0, dFill)
        aLyy = fLxy(nArray, yy0, yy0, dFill)
!        aLxy = fLxy(nArray, xx, yy, dFill)
        if (aLxx .eq. dFill.or.aLyy .eq. dFill) then
            f1CORR = dFill
        else
            if(iLog.eq.0) then
                xx = xx0
                yy = yy0
            else
                xx = DLOG10(max(dLOW,xx0))
                yy = DLOG10(max(dLOW,yy0))
            end if
            
            
            aLxx = fLxy(nArray, xx, xx, dFill)
            aLyy = fLxy(nArray, yy, yy, dFill)
            aLxy = fLxy(nArray, xx, yy, dFill)
            
            f1CORR = 1.D0 - aLxy/SQRT(aLxx * aLyy)
            
!            print*,xx0
!            print*,xx
!            print*,yy0
!            print*,yy
!            print*,aLxx,aLyy,aLxy,aLxy/SQRT(aLxx * aLyy)
        end if
        return
    END
    !------------------------------------------------------------------      
    REAL(8) function fVARIANCE(n, xx, dFill)
        !!ARGUMENTS:
        INTEGER n
        REAL(8), DIMENSION(n) :: xx
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        INTEGER j, k
        REAL(8) avg, avg2
        avg = fAVG(n, xx, dFill)
        avg2 = 0d0
        k = 0
        do j = 1, n
            if (xx(j) .ne. dFill) then
                k = k + 1
                avg2 = avg2 + xx(j) * xx(j)
            end if
        end do
        if (k < 1) then
            fVARIANCE = dFill
        else
            fVARIANCE = avg2/REAL(k) - avg * avg
        end if
        
        if(dabs(fVARIANCE).le.dLOW .or. fVARIANCE.lt.0.D0) then
            fVARIANCE = 0.0D0
        end if
        
        return
    END !!function fVARIANCE     
    !------------------------------------------------------------------      
    REAL(8) function fSTDDEV(n, xx, dFill)
        !!standard deviation
        !!ARGUMENTS:
        INTEGER n
        REAL(8), DIMENSION(n) :: xx
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) variance
        variance = fVARIANCE(n, xx, dFill)
                
        if (variance .eq. dFill) then
            fSTDDEV = dFill
        else
            fSTDDEV = dsqrt(variance)
        end if
        return
    END !!function fSTDDEV
    !------------------------------------------------------------------      
    REAL(8) function fSTDERR(n, xx, dFill)
        !!standard error
        !!ARGUMENTS:
        INTEGER n
        REAL(8), DIMENSION(n) :: xx
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        REAL(8) sd
        sd = fSTDDEV(n, xx, dFill)
        if (sd .eq. dFill) then
            fSTDERR = dFill
        else
            fSTDERR = sd/REAL(n)
        end if
        return
    END !!function fSTDERR
    !------------------------------------------------------------------
    !! [1] OBJECTIVE FUNCTION: END
    !!=================================================================

    !!=================================================================
    !! [2] ARRAY FUNCTION: BEGIN
    !------------------------------------------------------------------

    REAL(8) FUNCTION fMAXarray(nn, dobs, dFill)
        !!ARGUMENTS:
        REAL(8) dFill
        INTEGER nn
        REAL(8), DIMENSION(nn) :: dobs
        
        !!LOCAL VARIABLES:
        INTEGER j
        fMAXarray = dFill
        do j = 1, nn
            if (dobs(j) .ne. dFill.and.dobs(j) .gt. fMAXarray) then
                fMAXarray = dobs(j)
            end if
        end do
        return
    END !!FUNCTION
    !------------------------------------------------------------------
    REAL(8) FUNCTION fMINarray(nn, dobs, dFill)
        !!ARGUMENTS:
        REAL(8) dFill
        INTEGER nn
        REAL(8), DIMENSION(nn) :: dobs
        
        !!LOCAL VARIABLES:
        INTEGER j
        
        fMINarray = dHIGH
        do j = 1, nn
            if (dobs(j) .ne. dFill.and.dobs(j) .lt. fMINarray) then
                fMINarray = dobs(j)
            end if
        end do
        if (fMINarray .eq. dHIGH) then
            fMINarray = dFill
        end if
        return
    END !!FUNCTION
    !------------------------------------------------------------------
    REAL(8) FUNCTION fRANGEarray(nn, dobs, dFill)
        !!Normalized Root-Mean-Square Error
        !!ARGUMENTS:
        REAL(8) dFill
        INTEGER nn
        REAL(8), DIMENSION(nn) :: dobs
        
        !!LOCAL VARIABLES
        REAL(8) dobs_max, dobs_min, dobs_range

        dobs_max = fMAXarray(nn, dobs, dFill)
        dobs_min = fMINarray(nn, dobs, dFill)
        if (dobs_max .gt. dobs_min)then
            fRANGEarray = dobs_max - dobs_min
        elseif (dobs_max .ne. dFill) then
            fRANGEarray = dobs_max
        else !!dobs_max = dobs_min = dFill
            fRANGEarray = dFill
        end if
        return
    END !!FUNCTION
    !------------------------------------------------------------------      
    SUBROUTINE Array_Normalize(n, xx, dFill)
        !!Normalize an Array to make the sum = 1.0
        !!ARGUMENTS:
        INTEGER n
        REAL(8), DIMENSION(n) :: xx
        REAL(8) dFill
        
        !!LOCAL VARIABLES:
        INTEGER j
        REAL(8) sum1
        !            REAL(8) fSUM  !!function
        sum1 = fSUM(n, xx, dFill)
        do j = 1, n
            if (xx(j) .ne. dFill) then
                xx(j) = xx(j)/sum1
            else
                xx(j) = 0
            end if
        end do
        return
    END !!subroutine Array_Normalize
    
   
    !------------------------------------------------------------------
    !!    subroutine: randomly select m distinct INTEGERs from the range [1..npg]      
    SUBROUTINE selectINT(npg, m, lcs)
        !!ARGUMENTS:
        INTEGER, intent(in)  :: npg, m
        INTEGER, intent(out) :: lcs(m)
        
        !!LOCAL VARIABLES:
        INTEGER k, k1, lpos
        
        REAL(8) randx
        if (m .ge. npg) then
            do k = 1, npg
                lcs(k) = k
            end do
        else !!m < npg
            randx = rand() !rand = ran1(iseed1)
            lcs(1) = 1 + int(npg + 0.5 - sqrt((npg + .5)**2 - &
            npg * (npg + 1) * randx))
            do k = 2, m
                60 randx = rand() !rand = ran1(iseed1)
                lpos = 1 + int(npg + 0.5 - sqrt((npg + .5)**2 - &
                npg * (npg + 1) * randx))
                do k1 = 1, k - 1
                    if (lpos .eq. lcs(k1)) go to 60
                end do
                lcs(k) = lpos
            end do
        end if
        !      write(*,*)'subroutine lcs = ',lcs
        return
    END !!subroutine selectINT

  !!------------------------------------------------------------------
    SUBROUTINE sort(n,m,rb,ra)
    !c$debug      
    !c
    !c
    !c  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
    !c  BY W.H. PRESS ET AL., pp. 233-234
    !c
    !c  LIST OF VARIABLES
    !c     ra(.) = array to be sorted
    !c     rb(.,.) = arrays ordered corresponding to rearrangement of ra(.)
    !c     wk(.,.), iwk(.) = local varibles

        !!ARGUMENTS:
          INTEGER n,m
!          REAL(8) ra(n),rb(2000,m)
          REAL(8) ra(n),rb(:,:)

          !!LOCAL VARIABLES
!          REAL(8) wk(2000,m)
!          REAL(8),ALLOCATABLE :: wk(:,:)
          REAL(8) wk(n,m)
          INTEGER iwk(n)
          INTEGER i,j
    !      dimension ra(n),rb(2000,m),iwk(n)
    !      REAL(8), allocatable:: wk(:,:)
!          allocate(wk(n, m))
          CALL indexx(n, ra, iwk)
          do i = 1, n
            wk(i,1) = ra(i)
          end do

          do i = 1, n
            ra(i) = wk(iwk(i),1)
          end do

          do i = 1, n
            do j = 1, m
                wk(i,j) = rb(i,j)
            end do
          end do

          do i = 1, n
            do j = 1, m
                rb(i,j) = wk(iwk(i),j)
            end do
          end do
!          DEALLOCATE(wk)
    !c  END OF SUBROUTINE SORT
          return
    END
      
!!------------------------------------------------------------------
SUBROUTINE sort1(n,ira)
!c$debug      
!c
!c
!c  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!c  BY W.H. PRESS ET AL., pp. 231
!c
!c  LIST OF VARIABLES
!c     ira(.) = INTEGER array to be sorted
    !!ARGUMENTS:
      INTEGER n
      INTEGER:: ira(n)

    !!LOCAL VARIABLES:  
      INTEGER l,ir,i,j
!      dimension ra(n)

!      INTEGER ra, rra
      INTEGER irra

      l = (n / 2) + 1
      ir = n
   10 continue
      if (l .gt. 1) then
      l = l - 1
      irra = ira(l)
      else
      irra = ira(ir)
      ira(ir) = ira(1)
      ir = ir - 1
      if (ir .eq. 1) then
      ira(1) = irra
      return
      end if
      end if
      i = l
      j = l + l
   20 if (j .le. ir) then
      if (j .lt. ir) then
      if (ira(j) .lt. ira(j + 1)) j = j + 1
      end if
      if (irra .lt. ira(j)) then
      ira(i) = ira(j)
      i = j
      j = j + j
      else
      j = ir + 1
      end if
      goto 20
      end if
      ira(i) = irra
      goto 10

!c  END OF SUBROUTINE SORT1
END

!!------------------------------------------------------------------
SUBROUTINE indexx(n,arr,indx)
    !!ARGUMENTS:
    INTEGER n
    INTEGER indx(n)
    REAL(8) arr(n)

    !!LOCAL VARIABLES:
    INTEGER, PARAMETER :: M = 7
    INTEGER, PARAMETER :: NSTACK=50
!        Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that arr(indx(j))
!        is in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr are not changed.
    INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
    REAL(8) a
    do 11 j=1,n
        indx(j)=j
 11 continue !! end do !!11
    jstack=0
    l=1
    ir=n
  1 if(ir-l.lt.M)then
    do 13 j=l+1,ir
    indxt=indx(j)
    a=arr(indxt)
    do 12 i=j-1,l,-1
        if(arr(indx(i)).le.a)goto 2
        indx(i+1)=indx(i)
 12 continue !! end do !12
    i=l-1
  2 indx(i+1)=indxt
 13 continue !!  end do !!13
    if(jstack.eq.0) return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
    else
    k=(l+ir)/2
    itemp=indx(k)
    indx(k)=indx(l+1)
    indx(l+1)=itemp
    if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
    end if
    if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
    end if
    if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
    end if
    i=l+1
    j=ir
    indxt=indx(l+1)
    a=arr(indxt)
  3 continue
    i=i+1
    if(arr(indx(i)).lt.a)goto 3
  4 continue
    j=j-1
    if(arr(indx(j)).gt.a)goto 4
    if(j.lt.i)goto 5
    itemp=indx(i)
    indx(i)=indx(j)
    indx(j)=itemp
    goto 3
  5 indx(l+1)=indx(j)
    indx(j)=indxt
    jstack=jstack+2
    if(jstack.gt.NSTACK) print*, "NSTACK too small in indexx"

    if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
    else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
    end if
    end if
    goto 1
END 

    !!-------------------------------------------------------------------        
    !! [2] ARRAY FUNCTION: END
    !!===================================================================

    !!===================================================================     
    !! [3] DATETIME FUNCTION: BEGIN
    !!-------------------------------------------------------------------
    INTEGER FUNCTION nDaysUPMon(iYear, iMon)
        !!Julian day of the last day of a month (iMon) in a year (iYear)
        !!ARGUMENTS:
        INTEGER iYear, iMon
        
        !!LOCAL VARIABLES:
        INTEGER i, nDOM1(12), nDOM2(12)
        !            INTEGER nDaysofYear !!FUNCTION
        DATA nDOM1 /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
        Data nDOM2 /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

        SELECT CASE(iMON)
        CASE(1)
            nDaysUPMon = nDOM1(1)
        CASE (2:12)
            nDaysUPMon = 0
            do i = 1, iMON
                nDaysUPMon = nDaysUPMon + nDOM1(i)
            end do
            if (nDaysofYear(iYear) == 366) nDaysUPMon = nDaysUPMon + 1
        END SELECT

        return
    END !!FUNCTION nDaysUPMon

    !------------------------------------------------------------------

    INTEGER FUNCTION nDaysofMon(iYear, iMon)
        !!# of days in a specific (iYear, iMon)
        !!ARGUMENTS:
        INTEGER iYear, iMon
        
        !!LOCAL VARIABLES:
        INTEGER i, nDOM1(12), nDOM2(12)
        DATA nDOM1 /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
        Data nDOM2 /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

        do i = 1, 12
            nDaysofMon = nDOM1(iMon)
        end do

        if (nDaysofYear(iYear) == 366.and.iMon == 2) then
            nDaysofMon = nDOM2(iMon)
        end if

        return
    END !!FUNCTION nDaysofMon
    !------------------------------------------------------------------
    INTEGER FUNCTION nDaysofYears(iyear_beg, iyear_end)
        !!ARGUMENTS
        INTEGER iyear_beg, iyear_end
        
        !!LOCAL VARIABLES:
        INTEGER i

        nDaysofYears = -1
        if (iyear_end .lt. iyear_beg) then
            write(*, *) 'iyear_end < iyear_beg...'
        else if (iyear_end .eq. iyear_beg) then
            nDaysofYears = nDaysofYear(iyear_beg)
        else
            nDaysofYears = 0
            do i = iyear_beg, iyear_end
                nDaysofYears = nDaysofYears + nDaysofYear(i)
            end do
        end if
        return
    END !!function nDaysofYears
    !------------------------------------------------------------------
    Integer FUNCTION nDaysofYear(iyear)
        !!# of days in a year
        !c	calculate the days of the year
        !!ARGUMENTS:
        INTEGER iyear
        
        if (mod(iyear, 100) .eq. 0) then
            if (mod(iyear, 400) .eq. 0) then
                ndaysofyear = 366
            else
                ndaysofyear = 365
            end if
        else
            if (mod(iyear, 4) .eq. 0) then
            ndaysofyear = 366
        else
            ndaysofyear = 365
            end if
        end if
    END !!FUNCTION ndaysofyear
    !------------------------------------------------------------------
    Integer FUNCTION nMonsbwDates(sDate_beg, sDate_end)
        !!# OF MONTHs between two Year_Month (boundaries included)
        !!ARGUMENTS:
        CHARACTER(len = 8) sDate_beg, sDate_end
        
        !!LOCAL VARIABLES:
        INTEGER iYear_beg, iYear_end, iMon_beg, iMon_end, iDay_beg, iDay_end

        CALL sDate2YMD(sDate_beg, iYear_beg, iMon_beg, iDay_beg)
        CALL sDate2YMD(sDate_end, iYear_end, iMon_end, iDay_end)

        nMonsbwDates = (iYear_end - iYear_beg - 1) * 12 + &
        (12 - iMon_beg + 1) + iMon_end
        return
    END !!FUNCTION nMonsbwDates
    
    !------------------------------------------------------------------
    Integer FUNCTION nYearsbwDates(sDate_beg, sDate_end)
        !!# OF YEARs between two Year_Month (boundaries included)
        !!ARGUMENTS:
        CHARACTER(len = 8) sDate_beg, sDate_end
        
        !!LOCAL VARIABLES:
        INTEGER iYear_beg, iYear_end, iMon_beg, iMon_end, iDay_beg, iDay_end

        CALL sDate2YMD(sDate_beg, iYear_beg, iMon_beg, iDay_beg)
        CALL sDate2YMD(sDate_end, iYear_end, iMon_end, iDay_end)

        nYearsbwDates = iYear_end - iYear_beg + 1
        return
    END !!FUNCTION nMonsbwDates
    !------------------------------------------------------------------
    Integer FUNCTION nMonths(iYearMon_begin, iYearMon_end)
        !!# OF MONTHs between two Year_Month (boundaries included)
        !!ARGUMENTS:
        INTEGER iYearMon_begin, iYearMon_end
        !!LOCAL VARIABLES:
        INTEGER iYear_begin, iYear_end, iMon_begin, iMon_end
        
        iYear_begin = iYearMon_begin/100
        iYear_end = iYearMon_end/100
        iMon_begin = iYearMon_begin - iYear_begin * 100
        iMon_end = iYearMon_end - iYear_end * 100
        nMonths = (iYear_end - iYear_begin - 1) * 12 + &
        (12 - iMon_begin + 1) + iMon_end
        return
    END !!FUNCTION nMonths
    
    !------------------------------------------------------------------
    Integer FUNCTION nDaysbwDates(sDate_beg, sDate_end)
        !!# OF days between two dates (boundaries included)
        !!ARGUMENTS:
        CHARACTER(len = 8) sDate_beg, sDate_end
        !!LOCAL VARIABLES:
        INTEGER iYear_beg, iYear_end, iMon_beg, iMon_end, iDay_beg, iDay_end
        INTEGER jday_beg, jday_end

        nDaysbwDates = -1
        CALL sDate2YMD(sDate_beg, iYear_beg, iMon_beg, iDay_beg)
        CALL sDate2YMD(sDate_end, iYear_end, iMon_end, iDay_end)
        If (iYear_end .lt. iYear_beg) then
!            Write (*, *) "WARNING: Date_beg=", sDate_beg, "> Date_end=",sDate_end, " in 'nDaysbwDates'"
            Go to 200
        else if (iYear_end .eq. iYear_beg) then
            if (iMon_end .lt. iMon_beg) then
!            Write (*, *) "WARNING: Date_beg=", sDate_beg, "> Date_end=",sDate_end, " in 'nDaysbwDates'"
            Go to 200
        else if (iMon_end .eq. iMon_beg) then
            if (iDay_end .lt. iDay_beg) then
!            Write (*, *) "WARNING: Date_beg=", sDate_beg, "> Date_end=",sDate_end, " in 'nDaysbwDates'"
            Go to 200
            end if
            end if
        end if

        jday_beg = iJulianDay(iYear_beg, iMon_beg, iDay_beg)
        jday_end = iJulianDay(iYear_end, iMon_end, iDay_end)

        if (iYear_end .eq. iYear_beg) then
            nDaysbwDates = jday_end - jday_beg + 1
        else if (iYear_end .eq. (iYear_beg + 1)) then
            nDaysbwDates = (nDaysofYear(iYear_beg) - jday_beg + 1) + jday_end
        else
            nDaysbwDates = (nDaysofYear(iYear_beg) - jday_beg + 1) + jday_end
            nDaysbwDates = nDaysbwDates + nDaysofYears(iYear_beg + 1, iYear_end - 1)
        end if
    200 Return
    END !!FUNCTION nDaysbwDates

    !------------------------------------------------------------------
    INTEGER FUNCTION iJulianDay(iYear, iMonth, iDay)
        !!ARGUMENTS:
        INTEGER iYear, iMonth, iDay
        
        !!LOCAL VARIABLES:
        INTEGER nday
!        INTEGER incmon(13), leapinc(13)

        if (iMonth .eq. 1) then
            iJulianDay = iDay
        else
            nday = nDaysUPMon(iYear, iMonth - 1)
            iJulianDay = nday + iDay
        end if
        Return
    End !!FUNCTION iJulianDay

    !!-----------------------------------------------------------------
    SUBROUTINE sDate_After(nday,sDate_beg,sDate_end)
        !!return the DATE (sDate_end) that is nday after sDate_beg
        !!e.g., sDate_beg = 7/9/2015, sDate_After(2,sDate_beg,sDate_end) = 7/10/2015
        !!sDate_beg: YYYYMMDD, beginning date

        !!ARGUMENTS:
        INTEGER, intent(in) :: nday
        CHARACTER(len=8), intent(in)    :: sDate_beg
        CHARACTER(len=8), intent(out)   :: sDate_end

        !!LOCAL VARIABLES:
        INTEGER iyr0, imo0, ida0, iyr, imo, ida
        INTEGER nda_in_ym
        INTEGER i,j !!,nda,nmo, ibeg, iend,nbe

        CALL sDate2YMD(sDate_beg,iyr0,imo0,ida0)

        iyr = iyr0
        imo = imo0
        ida = ida0  - 1
        do j=1,nday
            nda_in_ym = nDaysofMon(iyr,imo)
            ida = ida + 1
            if(ida.le.nda_in_ym) then
                ida = ida
            else
                ida = 1
                imo = imo + 1
                if(imo.le.12) then
                    iyr = iyr
                else
                    iyr = iyr+1
                    imo = 1
                end if
            end if
        end do
        CALL sYMD2Date(iyr,imo,ida,sDate_end)
        RETURN
    END !!subroutine sOUT_OPT
     !!-----------------------------------------------------------------
    SUBROUTINE sYM_After(nmon,sYM_beg,sYM_end)
        !!return the yyyymm (sYM_end) that is nmon after sYM_beg
        !!e.g., sYM_beg = 9/2015, sYM_After(2,sYM_beg,sYM_end) = 10/2015
        !!sYM_beg: YYYYMM, beginning year-month

        !!ARGUMENTS:
        INTEGER, intent(in) :: nmon
        CHARACTER(len=6), intent(in)    :: sYM_beg
        CHARACTER(len=6), intent(out)   :: sYM_end

        !!LOCAL VARIABLES:
        CHARACTER(len=8) sDate_beg
        INTEGER iyr0, imo0, ida0, iyr, imo
        INTEGER i,j 
        
        sDate_beg = sYM_beg//"01"
!        write(sDate,'(A6,I1,I1)')sYM_beg,0,iDay 
        CALL sDate2YMD(sDate_beg,iyr0,imo0,ida0)

        iyr = iyr0
        imo = imo0 - 1
        do j=1,nmon
            imo = imo + 1
            if(imo.gt.12) then
                iyr = iyr+1
                imo = 1
            end if
        end do
        if(imo.lt.10) then
            write(sYM_end,'(I4,I1,I1)')iyr,0,imo 
        else
            write(sYM_end,'(I4,I2)')iyr,imo
        end if
        RETURN
    END !!subroutine sOUT_OPT
    !------------------------------------------------------------------
    Subroutine sDate2YMD(sDate, iYear, iMon, iDay)
        !!ARGUMENTS:
        CHARACTER(len = 8) sDate
        INTEGER iYear, iMon, iDay
        
        !!LOCAL VARIABLES:
        character(len = 4) sRead

        sRead = sDate(1:4)
        read(sRead, *) iYear
        sRead = sDate(5:6)
        read(sRead, *) iMon
        sRead = sDate(7:8)
        read(sRead, *) iDay
        return
    END !!subroutine sDate2YMD
    !------------------------------------------------------------------
    Subroutine sYMD2Date(iYear, iMon, iDay, sDate)
        !!ARGUMENTS:
        CHARACTER(len = 8), intent(out):: sDate
!        CHARACTER(len = 1), intent(out):: sFormat
        INTEGER,            intent(in) :: iYear, iMon, iDay
        
        !!LOCAL VARIABLES:
        character(len = 6) sYM
        if(iMon.lt.10) then
            write(sYM,'(I4,I1,I1)')iYear,0,iMon 
        else
            write(sYM,'(I4,I2)')iYear,iMon
        end if
        
        if(iDay.lt.10) then
            write(sDate,'(A6,I1,I1)')sYM,0,iDay 
        else
            write(sDate,'(A6,I2)')sYM,iDay
        end if
        return
    END !!subroutine sDate2YMD
    !------------------------------------------------------------------
    Subroutine sInt2Str(inum_in,nLen,str_out)
        !!ARGUMENTS:
        CHARACTER(*), intent(out):: str_out
!        CHARACTER(len = 1), intent(out):: sFormat
        INTEGER,            intent(in) :: inum_in,nLen
        
        !!LOCAL VARIABLES:
!        character(len = 6) sYM
        INTEGER ndigit,nzero,i
        CHARACTER(len=100) sFormat
        
        str_out = "0"
        ndigit = int(DLOG10(dble(inum_in)))
        ndigit = ndigit + 1
        if (ndigit>=nLen) then
!            write(*,*)"Error: # of digits of INPUT_integer >= Length of OUT_string!"
            write(sFormat,*)"(I",ndigit,")" 
            write(str_out,sFormat)inum_in
        else
            nzero = nLen - ndigit
            write(sFormat,*)"(A",nzero, ", I",ndigit,")"
            do i = 1, nzero - 1
                str_out = str_out//"0"
            end do
            write(str_out,sFormat)str_out,inum_in
        end if
        return
    END !!subroutine sDate2YMD
    !------------------------------------------------------------------
    SUBROUTINE Sec2HMS(tSec, tHMS)
        !!convert Seconds to Hours-Minutes-Seconds
        !!ARGUMENTS:
        INTEGER tSec
        INTEGER, dimension(3) :: tHMS
        
        tHMS(1) = tSec/3600 !!hours
        tHMS(2) = (tSec - 3600 * tHMS(1))/60
        tHMS(3) = tSec - 3600 * tHMS(1) - 60 * tHMS(2)
        return
    END !!subroutine
    

    !------------------------------------------------------------------     
    !! [3] DATETIME FUNCTION: END
    !!=================================================================

         
    !!=================================================================
    !! [4] RANDOM NUMBER FUNCTION: BEG
    !------------------------------------------------------------------
    
    !!------------------------------------------------------------------ 
    SUBROUTINE iRandSeedGen(iSeed)
        !        VALUE(1):	The year 
        !        VALUE(2):	The month 
        !        VALUE(3):	The day of the month 
        !        VALUE(4):	Time difference with UTC in minutes 
        !        VALUE(5):	The hour of the day 
        !        VALUE(6):	The minutes of the hour 
        !        VALUE(7):	The seconds of the minute 
        !        VALUE(8):	The milliseconds of the second 
        !!ARGUMENTS:
        INTEGER, intent(inout) :: iSeed
        
        !!LOCAL VARIABLES:
        Integer time(8)
        
        CALL DATE_AND_TIME(values = time) ! Get the current time
        iSeed = time(4) * (360000 * time(5) + 6000 * time(6) + 100 * time(7) + time(8))
        RETURN
    END
    !------------------------------------------------------------------
    
    REAL(8) FUNCTION gasdev(idum)
    !!MODIFIED BY WGS, REPLACING gasdev0
    !!HAS BEEN TESTED BY WGS
    
    !c  GENERATE A STANDARD NORMAL DISTRIBUTION RANDOM NUMBER ~ NORM(0,1)
    !c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
    !      implicit REAL*8 (a-h,o-z)
    !      INTEGER,data:: iset / 0 /
        !!ARGUMETNS:
          INTEGER idum

        !!LOCAL VARIABLES:
          REAL(8) v1,v2,r,fac,gset
          INTEGER iset 
          DATA iset / 0 /
          
          if (iset .eq. 0) then
        1 v1 = (2. * rand()) - 1.
          v2 = (2. * rand()) - 1.
          r = (v1 ** 2) + (v2 ** 2)
          if (r .ge. 1.) goto 1
          fac = sqrt(- ((2. * log(r)) / r))
          gset = v1 * fac
          gasdev = v2 * fac
          iset = 1
          else
          gasdev = gset
          iset = 0
          end if

    !c  END OF SUBROUTINE GASDEV
          return
    END !!gasdev
    
    !!------------------------------------------------------------------
    REAL(8) FUNCTION ran1(idum)
!$debug
!c
!c
!c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!      implicit real*8 (a-h,o-z)
      INTEGER idum
      INTEGER j, iff,ix1,ix2,ix3
      
      common /rancom/ ix1,ix2,ix3
      data iff /0/
      REAL(8) r(97)
      INTEGER, PARAMETER :: m1 = 259200, ia1 = 7141, ic1 = 54773 
      INTEGER, PARAMETER :: m2 = 134456, ia2 = 8121, ic2 = 28411
      INTEGER, PARAMETER :: m3 = 243000, ia3 = 4561, ic3 = 51349
      REAL(8), PARAMETER :: rm1 = 3.8580247e-6, rm2 = 7.4373773e-6
      
      
      if ((idum .lt. 0) .or. (iff .eq. 0)) then
      iff = 1
      ix1 = mod(ic1 - idum,m1)
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = mod(ix1,m2)
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix3 = mod(ix1,m3)
      do 11 j = 1, 97
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = mod((ia2 * ix2) + ic2,m2)
      r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
   11 continue
      idum = 1
      end if
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = mod((ia2 * ix2) + ic2,m2)
      ix3 = mod((ia3 * ix3) + ic3,m3)
      j = 1 + ((97 * ix3) / m3)
!      if ((j .gt. 97) .or. (j .lt. 1)) pause
      if ((j .gt. 97) .or. (j .lt. 1)) j=97
      ran1 = r(j)
      r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1

!c  END OF SUBROUTINE RAN1
    RETURN
    END

!!------------------------------------------------------------------
    REAL(8) FUNCTION gasdev0(idum)
!$debug
!c
!c
!c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!      implicit real*8 (a-h,o-z)
      INTEGER idum
      INTEGER iset
      REAL(8) v1,v2,r,fac,gset
      data iset / 0 /
      if (iset .eq. 0) then
    1     v1 = (2. * ran1(idum)) - 1.
          v2 = (2. * ran1(idum)) - 1.
          r = (v1 ** 2) + (v2 ** 2)
          if (r .ge. 1.) goto 1
          fac = sqrt(- ((2. * log(r)) / r))
          gset = v1 * fac
          gasdev0 = v2 * fac
          iset = 1
      else
          gasdev0 = gset
          iset = 0
      end if

!c  END OF SUBROUTINE GASDEV
      return
      END


    !------------------------------------------------------------------     
    !! [4] RANDOM NUMBER FUNCTION: END
    !!=================================================================


END MODULE MOD_USRFS


!!===========================================================
!        SUBROUTINE Write_filecio(dir_swatio,nbyr1,iyear1,iprint1)
!            IMPLICIT NONE
!            CHARACTER*20 dir_swatio
!            INTEGER nbyr1,iyear1,iprint1
!            INTEGER j,nline,eof,idata
!!            INTEGER,DIMENSION(:),ALLOCATABLE ::idata
!            CHARACTER*150,DIMENSION(:),ALLOCATABLE ::sdata
!            write(*,*)">>>UPDATE PARs in 'file.cio': ", &
!            "NBYR(Line 8),IYR(Line 9),IPRINT(Line 59)"
!            open(unit=101, &
!            file=trim(dir_swatio)//'file.cio',status='old')
!            nline=0
!            read(101,*,iostat=eof)
!            do while(.not.(eof<0))
!               nline = nline+1 
!               read(101,*,iostat=eof)
!            enddo
!            close(101)
!!            ALLOCATE(idata(nline))
!            ALLOCATE(sdata(nline))
!            
!            open(unit=101,file=trim(dir_swatio)//'file.cio',status='old')
!            do j=1,7
!                read(101,'(a150)')sdata(j)
!            end do
!            read(101,'(I16,a150)')idata,sdata(8)
!            read(101,'(I16,a150)')idata,sdata(9)
!            do j=10,58
!                read(101,'(a150)')sdata(j)
!            end do
!            read(101,'(I16,a150)')idata,sdata(59)
!            do j=60,nline
!                read(101,'(a150)')sdata(j)
!            end do
!            close(101)
!            
!            open(unit=101,file=trim(dir_swatio)//'file.cio')
!            do j=1,7
!                write(101,'(a150)')sdata(j)
!            end do
!            write(101,'(I16,a150)')nbyr1,sdata(8)
!            write(101,'(I16,a150)')iyear1,sdata(9)
!            do j=10,58
!                write(101,'(a150)')sdata(j)
!            end do
!            write(101,'(I16,a150)')iprint1,sdata(59)
!            do j=60,nline
!                write(101,'(a150)')sdata(j)
!            end do
!            close(101)
!            RETURN
!        END SUBROUTINE
!!===========================================================
!        SUBROUTINE SGLB_INI_VAR()
!            USE parm
!            sGLB%dir_swatio = ''                          !!dir for swat in & out files
!            sGLB%itype_opt = 0                               !!type of data for calibration; 0-subbasin with flow,1-floodgate,2-Reservoir with volume, 3-subbasin with SWC !!original: iRESFLG
!            sGLB%nsub_opt = 0                                !!# of subbasins included in calibration !!noptsub,icodesub,icoderes,ndaysim,idaysim
!            sGLB%icode_sub = 0                               !!code of the subbasin with flow gage/reservoir, etc.
!            sGLB%icode_res = 0                               !!code of the reservoir in SWAT, ATTENTION, different from the code of corresponding subbasin
!            sGLB%nyear_warmup = 0                            !!# of years for SWAT warm-up
!            sGLB%nday_sim = 0                                !!# of days for simulation, particularly, for calibration, counting through SWAT simulations (should be= nday_sim0, but computed by different ways)
!            sGLB%iday_sim = 0                                !!the ith day in nday_sim
!            sGLB%imon_sim = 0                                !!the ith day in nday_sim
!            sGLB%iyear0 = 0                                 !!Beginning year of simulation
!            sGLB%nyear_sim = 0                              !!# of years simulated, excluding 'nyear_warmup'
!            sGLB%nday_sim0 = 0                              !!# of days simulated, computed by iyear0,nyear_sim, nyear_warmup
!            sGLB%nmon_sim = 0                             !!# of months simulated,e xcluding warmup
!            sGLB%iSWAT = 0                                  !!the ith SWAT simulation run
!            sGLB%iModel = 0                                 !!0: run SWAT; >0: SCEUA optimization
!            sGLB%rOBJ0 = 0                                     !!reference OBJ, e.g., observed runoff coefficient (RC)
!            sGLB%wOBJ(:) = 0
!            return
!        END SUBROUTINE
!!===========================================================
!        SUBROUTINE SGLB_INI_ARRAY()
!            USE parm
!            sGLB%ivar_par(:) = 0 !!variation method: 1-absolute value; 2-relative value, multiplied by a factor
!            sGLB%iopt_par(:) = 0 !!optimize or not: 1-optimize, 0-keep original value
!            sGLB%isub_opt(:) = 0 !!subbasin IDs to be calibrated 
!            sGLB%xua(:) = 0         !!parameter values 
!            sGLB%subsurlag(:) = 4.0   !!surlag for each subbasin 
!            sGLB%par_min(:) = 0    !!min and max of absolute parameter values
!            sGLB%par_max(:) = 0     !!min and max of absolute parameter values
!            sGLB%sub_area(:) = 0   !!subbasin area (km2)
!            sGLB%iGateDam(:) = 0 !!0-SubBasin,1-FloodGate,2-Res/Dam
!            sGLB%cwus(:) = 0        !!wateruse
!            sGLB%crwus(:) = 0       !!wateruse coef
!!            sGLB%data_sim(:,:) = 0 !!simulated time series, may include multiple variables, e.g., multiple subbasins within a HUC8 
!!            sGLB%data_obs(:,:) = 0 !!observed time series for calibration
!            return
!        END SUBROUTINE
!!===========================================================

!        SUBROUTINE SGLB_WRITE(sGLB1)
!          USE parm
!          TYPE(sSCE_GLOBAL) sGLB1
!          PARAMETER(n = 5)
!          write(*,*)"itype_opt = ",sGLB1%itype_opt                               !!type of data for calibration; 0-subbasin with flow,1-floodgate,2-Reservoir with volume, 3-subbasin with SWC !!original: iRESFLG
!          write(*,*)"nsub_opt = ", sGLB1%nsub_opt                                !!# of subbasins included in calibration !!noptsub,icodesub,icoderes,ndaysim,idaysim
!          write(*,*)"icode_sub = ", sGLB1%icode_sub                               !!code of the subbasin with flow gage/reservoir, etc.
!          write(*,*)"icode_res = ", sGLB1%icode_res                               !!code of the reservoir in SWAT, ATTENTION, different from the code of corresponding subbasin
!          write(*,*)"nyear_warmup  = ", sGLB1%nyear_warmup                            !!# of years for SWAT warm-up
!          write(*,*)"nday_sim = ", sGLB1%nday_sim                                !!# of days for simulation, particularly, for calibration, counting through SWAT simulations (should be= nday_sim0, but computed by different ways)
!          write(*,*)"iday_sim = ", sGLB1%iday_sim                                !!the ith day in nday_sim
!          write(*,*)"imon_sim = ", sGLB1%imon_sim                                !!the ith day in nday_sim
!          write(*,*)"iyear0 = ", sGLB1%iyear0                                 !!Beginning year of simulation
!          write(*,*)"nyear_sim = ", sGLB1%nyear_sim                              !!# of years simulated, excluding 'nyear_warmup'
!          write(*,*)"nday_sim0 = ", sGLB1%nday_sim0                              !!# of days simulated, computed by iyear0,nyear_sim, nyear_warmup
!          write(*,*)"nmon_sim = ", sGLB1%nmon_sim                             !!# of months simulated,e xcluding warmup
!          write(*,*)"iSWAT = ", sGLB1%iSWAT                                  !!the ith SWAT simulation run
!          write(*,*)"rOBJ0 = ", sGLB1%rOBJ0                                     !!reference OBJ, e.g., observed runoff coefficient (RC)
!          write(*,*)"wOBJ = ",sGLB1%wOBJ
!          write(*,*)"ivar_par = ", sGLB1%ivar_par !!variation method: 1-absolute value; 2-relative value, multiplied by a factor
!          write(*,*)"iopt_par = ", sGLB1%iopt_par !!optimize or not: 1-optimize, 0-keep original value
!          write(*,*)"isub_opt = ", sGLB1%isub_opt !!subbasin IDs to be calibrated 
!          write(*,*)"iGateDam = ", sGLB1%iGateDam(1:n) !!0-SubBasin,1-FloodGate,2-Res/Dam
!          write(*,*)"xua = ", sGLB1%xua         !!parameter values 
!          write(*,*)"subsurlag = ", sGLB1%subsurlag(1:n)   !!surlag for each subbasin 
!          write(*,*)"cwus = ", sGLB1%cwus(1:n)        !!wateruse
!          write(*,*)"crwus = ", sGLB1%crwus(1:n)       !!wateruse coef
!          write(*,*)"par_min = ", sGLB1%par_min     !!min and max of absolute parameter values
!          write(*,*)"par_max = ", sGLB1%par_max     !!min and max of absolute parameter values
!          write(*,*)"sub_area = ", sGLB1%sub_area(1:n)   !!subbasin area (km2)
!          write(*,*)"data_sim = ", sGLB1%data_sim(1,1:n) !!simulated time series, may include multiple variables, e.g., multiple subbasins within a HUC8 
!          write(*,*)"data_obs = ", sGLB1%data_obs(1,1:n) !!observed time series for calibration
!        END
!!===========================================================
!        Subroutine Read_Sub_Area
!     &        (datafile,sub_area)
!        IMPLICIT NONE
!        character*100 datafile
!        INTEGER nSub,i,itemp
!        REAL,DIMENSION(:), ALLOCATABLE :: sub_area
!        open(unit=101,file = datafile,status='old')
!        read(101,*)  !!"No. of Subbasins"
!        read(101,*)nSub
!        ALLOCATE(sub_area(nSub))
!        read(101,*)  !!head
!        do i = 1, nSub
!            read(101,*)itemp, sub_area(i)
!        end do
!        close(101)
!        RETURN
!        END
!=======================================================================
!      REAL function fHMLE(ndata, qobs,q)
!!!Heteroscedastic Maximum Likelihood Error
!!      implicit REAL*8 (a-h,o-z)
!!!      common /fnblk/ rlamda, ad
!!!      common /block1/ ndata, ns, iobj
!!!      common /block2/ p(1000), qobs(1000), s0(2)
!!!      common /block3/ q(1000), r(1000), b(1000), s(1000)
!      IMPLICIT NONE
!      REAL rlamda, ad, a, ex, rd,rn,w
!      INTEGER ndata,i,ict,isign,lcount
!      REAL q(ndata), qobs(ndata)
!      REAL ra(2)
!      REAL, PARAMETER::eps=5.d-02, del=5.d-02
!      INTEGER iflag
!      DATA iflag /0/
!
!!c  COMPUTE THE MEAN OF LOGARITHM OF OBSERVED FLOWS
!      if (iflag .eq. 0) then
!        ad = 0.d0
!        do 10 i = 1, ndata
!          ad = ad + dlog(dble(qobs(i)))
!   10   continue
!        ad = ad / dble(ndata)
!        rlamda = 1.d0
!        iflag = 1
!      end if
!
!!c  ESTIMATE THE LAMDA VALUE
!      lcount = 0
!      ict = 1
!      ra(1) = 0.d0
!      ra(2) = 0.d0
!   25 continue
!      lcount = lcount + 1
!      if(lcount .gt. 40) then
!        write(*,*) 'LAMDA ITERATION GO OVER 40', rlamda, ra(1), ra(2)
!        go to 50
!      end if
!      rd = 0.d0
!      rn = 0.d0
!      do 30 i = 1, ndata
!        a = dlog(dble(qobs(i))) / ad
!        w = qobs(i)**(2*(rlamda-1.d0))
!        rd = rd + w*(qobs(i) - q(i))**2
!        rn = rn + w*(qobs(i) - q(i))**2 * a
!   30 continue
!      ra(ict) = rn / rd - 1.d0
!      if (dabs(dble(ra(ict))) .le. eps) go to 50
!      isign = -1
!      if (ra(ict) .lt. 0.d0) isign = 1
!      rlamda = rlamda + isign * del
!      if (ict .eq. 2) go to 35
!      ict = 2
!      go to 25
!
!   35 continue
!      if (ra(1)*ra(2) .lt. 0.d0) go to 40
!      ra(1) = ra(2)
!      go to 25
!
!   40 continue
!      rlamda = rlamda - isign * del / 2.d0
!
!!c  COMPUTE HMLE
!   50 continue
!      fHMLE = 0.d0
!      ex = 2. * (rlamda - 1.)
!      do 60 i = 1, ndata
!        fHMLE = fHMLE + qobs(i)**ex * (qobs(i) -q(i))**2
!   60 continue
!      fHMLE = fHMLE / dble(ndata)
!      fHMLE = fHMLE / dexp(dble(ex * a))
!
!	if(fHMLE.gt.1.0E20) fHMLE = 1.0E20
!      return
!      end
!
!!!=======================================================================
!      REAL function fHMLE1(ndata, qobs,q)
!!      implicit REAL*8 (a-h,o-z)
!!!      common /fnblk/ rlamda, ad
!!!      common /block1/ ndata, ns, iobj
!!!      common /block2/ p(1000), qobs(1000), s0(2)
!!!      common /block3/ q(1000), r(1000), b(1000), s(1000)
!      IMPLICIT NONE
!      REAL rlamda, ad, ex, a
!      INTEGER ndata,i
!      REAL q(ndata), qobs(ndata)
!      REAL ra(2)
!      REAL, PARAMETER:: eps=5.d-02, del=5.d-02
!      INTEGER iflag
!      data iflag /0/
!
!!c  COMPUTE THE MEAN OF LOGARITHM OF OBSERVED FLOWS
!      if (iflag .eq. 0) then
!        ad = 0.d0
!        do 10 i = 1, ndata
!          ad = ad + dlog(dble(qobs(i)))
!   10   continue
!        ad = ad / dble(ndata)
!        rlamda = 0
!!!        iflag = 1
!      end if
!
!!c  COMPUTE HMLE
!   50 continue
!      fHMLE1 = 0.d0
!      ex = 2. * (rlamda - 1.)
!      do 60 i = 1, ndata
!        fHMLE1 = fHMLE1 + qobs(i)**ex * (qobs(i) -q(i))**2
!   60 continue
!      fHMLE1 = fHMLE1 / dble(ndata)
!      fHMLE1 = fHMLE1 / dexp(dble(ex * a))  !!PROBLEM: how to compute 'a'???
!
!	if(fHMLE1.gt.1.0E20) fHMLE1 = 1.0E20
!      return
!      end

!!========================================================
!!===================================================================
!      SUBROUTINE par_new(par,isub_cur,x3,ifopt,ivar3,nsub3,isub3)
!!!update SWAT parameter (single value)
!!!par: old value
!!!x3: new value
!!!isub_cur: curren subbasin ID in SWAT loop
!!!ifopt: whether the parameter to be optimized (1) or not (0)
!!!ivar3: variation method for parameter values: 1-absolute value; 2-relative value (old value multiplied by x3)
!!!nsub3: # of subbasins whose parameters to be calibrated (changed)
!!!isub3(nsub3): array to calibrated subbasins
!!      USE parm
!      IMPLICIT NONE
!      INTEGER jua,isub_cur,ifopt,ivar3,nsub3 !!loop;current subbasin in SWAT simulation, if optimization, variation method, #of calibrated subbasins
!      INTEGER isub3(nsub3) !!array of calibrated subbasins
!      REAL par,par_old,x3
!      par_old = par
!      if(ifopt.gt.0) then
!          do jua = 1, nsub3
!              if(isub3(jua).eq.isub_cur) then
!                  if(ivar3.eq.2) then  !relative value
!                    par = par_old*x3  !relative value, multiplied by a factor
!                  else
!                    par = x3 !absolute value
!                  end if
!              end if
!          end do
!      end if
!      return
!      END
!!!===================================================================
!      SUBROUTINE par_bsn_new(par,x3,ifopt,ivar3)
!!!update SWAT parameter (single value)
!!!par: old value
!!!x3: new value
!!!isub_cur: curren subbasin ID in SWAT loop
!!!ifopt: whether the parameter to be optimized (1) or not (0)
!!!ivar3: variation method for parameter values: 1-absolute value; 2-relative value (old value multiplied by x3)
!!!nsub3: # of subbasins whose parameters to be calibrated (changed)
!!!isub3(nsub3): array to calibrated subbasins
!!      USE parm
!      IMPLICIT NONE
!      INTEGER ifopt,ivar3  !!jua,isub_cur,nsub3 !!loop;current subbasin in SWAT simulation, if optimization, variation method, #of calibrated subbasins
!!!      INTEGER isub3(nsub3) !!array of calibrated subbasins
!      REAL par,par_old,x3
!      par_old = par
!      if(ifopt.gt.0) then
!!!          do jua = 1, nsub3
!!!              if(isub3(jua).eq.isub_cur) then
!                  if(ivar3.eq.2) then  !relative value
!                    par = par_old*x3  !relative value, multiplied by a factor
!                  else
!                    par = x3 !absolute value
!                  end if
!!!              end if
!!!          end do
!      end if
!      return
!      END
!!!===================================================================
!      SUBROUTINE par_res_new(par,ires_cur,x3,ifopt,ivar3,ires)
!!!update SWAT reservoir parameter (single value)
!!!par: old value
!!!x3: new value
!!!ires_cur: curren reservoir ID in SWAT loop
!!!ifopt: whether the parameter to be optimized (1) or not (0)
!!!ivar3: variation method for parameter values: 1-absolute value; 2-relative value (old value multiplied by x3)
!!!ires: calibrated reservoir #
!!      USE parm
!      IMPLICIT NONE
!      INTEGER ires_cur,ifopt,ivar3,ires !!loop;current subbasin in SWAT simulation, if optimization, variation method, #of calibrated subbasins
!!!      INTEGER isub3(nsub3) !!array of calibrated subbasins
!      REAL par,par_old,x3
!      par_old = par
!      if(ifopt.gt.0) then
!!!          do jua = 1, nsub3
!              if(ires.eq.ires_cur) then
!                  if(ivar3.eq.2) then  !relative value
!                    par = par_old*x3  !relative value, multiplied by a factor
!                  else
!                    par = x3 !absolute value
!                  end if
!              end if
!!!          end do
!      end if
!      return
!      END
!!!===================================================================
!      SUBROUTINE pars_new(pars,nly,isub_cur,x3,ifopt,ivar3,nsub3,isub3)
!!!update SWAT parameter (multiple values, e.g.,SOL_AWC, SOL_K, the same parameter with different values at multiple soil layers)
!!!pars: old values (e.g., at multiple soil layers)
!!!x3: new value (single value, e.g., only for 1st soil layer), parameter values for other layers will be changed proportionally
!!!isub_cur: curren subbasin ID in SWAT loop
!!!ifopt: whether the parameter to be optimized (1) or not (0)
!!!ivar3: variation method for parameter values: 1-absolute value; 2-relative value (old value multiplied by x3)
!!!nsub3: # of subbasins whose parameters to be calibrated (changed)
!!!isub3(nsub3): array to calibrated subbasins
!      IMPLICIT NONE
!      INTEGER nly !!# of soil layers
!      REAL,dimension(nly):: pars,pars_old
!      INTEGER j,jua,isub_cur,ifopt,ivar3,nsub3 !!loop;current subbasin in SWAT simulation, if optimization, variation method, #of calibrated subbasins
!      INTEGER isub3(nsub3) !!array of calibrated subbasins
!      REAL x3
!      pars_old = pars
!      if(ifopt.gt.0) then
!          do jua = 1, nsub3
!              if(isub3(jua).eq.isub_cur) then
!                  if(ivar3.eq.2) then  !relative value
!                    do j = 1, nly
!                        pars(j) = pars_old(j)*x3  !relative value, multiplied by a factor
!                    end do
!                  else  !absolute value
!                    pars(1) = x3  !absolute value
!                    do j = 2, nly
!                        pars(j) = x3*pars_old(j)/max(1E-6,pars_old(1))  !absolute value
!                    end do
!                  end if
!              end if
!          end do
!      end if
!      RETURN 
!      END

!!==============================================================================
!Subroutine Read_Month_Runoff_HUC8                               &
!             (datafile,sHUC8,iYearMon_begin,iYearMon_end,ndata,runoff)
!        IMPLICIT NONE
!    !!"mv01d_row_data.txt", 
!    !!# of columns: 1357 = 1st column + 1356 (113 hydrological years*12 months); Oct-1900~Sept-2013  
!    !!# of rows: 2111 = 1st row + 2110 HUC8
!        INTEGER,PARAMETER::nRow = 2110,nCol = 1356
!        INTEGER,PARAMETER::iYearMon1 = 190010, iYearMon2 = 201309
!        character*100 datafile 
!        character*8 sHUC8,sHUC8_read
!        INTEGER iYearMon_begin,iYearMon_end,ndata
!        INTEGER nMonths !!function CALL
!        
!        INTEGER i,iCol_begin,iCol_end
!        REAL runoff(ndata)
!        REAL runoff_all(nCol)
!        INTEGER ndata0
!        
!        iCol_begin = nMonths(iYearMon1,iYearMon_begin)
!        iCol_end = nMonths(iYearMon1,iYearMon_end)
!        ndata0 = iCol_end - iCol_begin + 1
!        if(ndata0.ne.ndata) then
!            write(*,*)"ndata (No. of data) should be ",ndata0
!            return
!        end if
!        open(unit=101,file = datafile,status='old')
!        read(101,*)  !!head
!        do i = 1,nRow
!            read(101,*)sHUC8_read,runoff_all(1:nCol)
!            if(sHUC8_read.eq.sHUC8) goto 200
!        end do
!200     runoff(1:ndata) = runoff_all(iCol_begin:iCol_end)
!        close(101)
!        return
!        END
!!!===========================================================
!        Subroutine Read_Month_NPS  &
!             (datafile,iYearMon_begin,iYearMon_end,mCOL,nMonth,dOBS)
!        IMPLICIT NONE
!    !!Read Monthly Nutrient
!    !!e.g. "NPS_03609750.dat", 
!    !!# of columns: 1357 = 1st column + 1356 (113 hydrological years*12 months); Oct-1900~Sept-2013  
!    !!# of rows: 2111 = 1st row + 2110 HUC8
!!        INTEGER,PARAMETER::nRow = 2110,nCol = 1356
!!        INTEGER,PARAMETER::iYearMon1 = 190010, iYearMon2 = 201309
!!        INTEGER, PARAMETER::mNPS = 5  !!# of NPS for NPS data:TN,TP,Sediment,NO2+NO3,NO3
!!        INTEGER mNPS
!        INTEGER mCOL !!= 1+2*mNPS  !!Mean Value + SD1st; 1st COLumn: Flow_cms,
!        character*100 datafile 
!        character*8 sHUC8,sHUC8_read
!        INTEGER iYearMon1,iYearMon_begin,iYearMon_end,nMonth
!        INTEGER nMonths !!function CALL
!        
!        INTEGER i,iRow_begin,iRow_end
!        REAL,DIMENSION(mCOL,nMonth):: dOBS
!!        REAL runoff_all(nCol)
!        INTEGER eof,ivoid,nMon0
!        
!!        mCOL = 1+2*mNPS !!Mean Value + SD1st; 1st COLumn: Flow_cms,
!        
!        open(unit=101,file = datafile,status='old')
!        read(101,*)  !!head
!        read(101,*)iYearMon1
!        read(101,*) !!head: Flow_cms	TN_ton	TP_ton 	SiO2_ton NO2+3_ton	NO3_ton	
!        
!        iRow_begin = nMonths(iYearMon1,iYearMon_begin)
!        iRow_end = nMonths(iYearMon1,iYearMon_end)
!        nMon0 = iRow_end - iRow_begin + 1
!        if(nMon0.ne.nMonth) then
!            write(*,*)"nMonth (No. of Months) should be ",nMon0
!            return
!        end if
!        
!        do i = 1,iRow_begin-1
!            read(101,*)  !!skip data
!        end do
!        dOBS(1:mCOL,1:nMonth) = -9999
!        do i = 1,nMonth
!            read(101,*,iostat=eof)ivoid,dOBS(1:mCOL,i)
!            if(eof<0) exit
!!            write(*,*)ivoid,dOBS(1:mCOL,i)
!        end do
!        close(101)
!        return
!        END
!!!===========================================================
!        SUBROUTINE Read_Daily_Obs &
!            (file_obs,iyr0,nyr_spinup,mobs,nday,dobs,mobs_read,dFill)
!            IMPLICIT NONE
!            character*100 file_obs
!            INTEGER iyr0, nyr_spinup  !!year_begin for SWAT simulation,# of years for spinup
!            INTEGER mobs,nday,mobs_read !!# of COL & ROW in dobs, mobs_read is the # of COL to be read
!            REAL dFill  !!-9999
!            REAL, DIMENSION(mobs,nday)::dobs
!            INTEGER nDaysofYear  !!FUNCTION
!            
!            INTEGER iFobs,iyr0_obs,j,k,eof
!            character*20 svoid
!
!            iFobs = 1001
!            dobs(:,:) = dFill
!            open(unit=iFobs,file=file_obs,status='unknown')
!            read(iFobs,*)svoid   !!yyyymmdd, start date
!            svoid = trim(svoid)
!            svoid = svoid(1:4)
!            read(svoid,*)iyr0_obs  !!the beginning year of the observations
!            do j=iyr0_obs,(iyr0+nyr_spinup - 1)  !!skip the data for warmup
!                do k =1,nDaysofYear(j)
!                    read(iFobs,*)
!                end do
!            end do
!            
!            do j = 1, nday  !!nbyr: # of simulation years; ATTENTION: please append lines of ANY data to meet the total # of nbyr*366 as not every year has 366 days
!                read(iFobs,*,iostat=eof)dobs(1:mobs_read,j)!!q: streamflow
!                if(eof<0) exit
!            end do
!            close(iFobs)
!            return
!        END SUBROUTINE
!!==============================================================================