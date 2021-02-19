MODULE MOD_OPT
! File:   MOD_OPT.F90
! Author: GANGSHENG WANG @ ORNL
! Updated: May 5, 2015
! Created on February 26, 2013, 11:02 AM    
USE MOD_OPT_TYPE
USE MOD_MEND
USE MOD_USRFS, ONLY: sort,sort1,indexx,gasdev,selectINT

    IMPLICIT NONE
    PRIVATE :: cce
    PRIVATE :: getpnt
    PRIVATE :: parstt
    PRIVATE :: comp
    PRIVATE :: chkcst
    
    PUBLIC :: SCEUA
!    PUBLIC :: gasdev
    
!    PUBLIC :: sort
!    PUBLIC :: sort1
!    PUBLIC :: indexx
    
    
    CONTAINS
    
    
SUBROUTINE SCEUA(sPAR_SCE, sPAR, sINI, sOUT)

!    (a,bl,bu,nopt,maxn,kstop,pcento,iseed,ngs,npg,nps,nspl,mings,iniflg,iprint)
!c$debug
!c
!c  MODIFIED by GANGSHENG WANG
!C  Environmental Sciences Division
!C  Oak Ridge National Laboratory
!C  Oak Ridge, TN 37831-6301
!C  March, 2013
!C  Updated: April 28, 2015

!c  SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
!c     -- Version 2.1
!c
!c  by QINGYUN DUAN
!c  DEPARTMENT OF HYDROLOGY & WATER RESOURCES
!c  UNIVERSITY OF ARIZONA, TUCSON, AZ 85721
!c  (602) 621-9360, email: duan@hwr.arizona.edu
!c
!c  WRITTEN IN OCTOBER 1990.
!c  REVISED IN AUGUST 1991
!c  REVISED IN APRIL 1992
!c
!c  STATEMENT BY AUTHOR:
!c  --------------------
!c
!c     This general purpose global optimization program is developed at
!c     the Department of Hydrology & Water Resources of the University
!c     of Arizona.  Further information regarding the SCE-UA method can
!c     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
!c     at the address and phone number listed above.  We request all
!c     users of this program make proper reference to the paper entitled
!c     'Effective and Efficient Global Optimization for Conceptual
!c     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
!c     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.
!c
!c
!c  LIST OF INPUT ARGUEMENT VARIABLES
!c
!c     a(.) = initial parameter set
!c     bl(.) = lower bound on parameters
!c     bu(.) = upper bound on parameters
!c     nopt = number of parameters to be optimized
!c
!c
!c  LIST OF SCE ALGORITHMIC CONTROL PARAMETERS:
!c
!c     ngs = number of complexes in the initial population
!c     npg = number of points in each complex
!c     npt = total number of points in initial population (npt=ngs*npg)
!c     nps = number of points in a sub-complex
!c     nspl = number of evolution steps allowed for each complex before
!c         complex shuffling
!c     mings = minimum number of complexes required, if the number of
!c         complexes is allowed to reduce as the optimization proceeds
!c     iseed = initial random seed
!c     iniflg = flag on whether to include the initial point in population
!c         = 0, not included
!c         = 1, included
!c     iprint = flag for controlling print-out after each shuffling loop
!c         = 0, print information on the best point of the population
!c         = 1, print information on every point of the population
!c
!c
!c  CONVERGENCE CHECK PARAMETERS
!c
!c     maxn = max no. of trials allowed before optimization is terminated
!c     kstop = number of shuffling loops in which the criterion value must
!c         chang by the given percentage before optimization is terminated
!c     pcento = percentage by which the criterion value must change in
!c         given number of shuffling loops
!c     ipcnvg = flag indicating whether parameter convergence is reached
!c         (i.e., check if gnrng is less than 0.001)
!c         = 0, parameter convergence not satisfied
!c         = 1, parameter convergence satisfied
!c
!c
!c  LIST OF LOCAL VARIABLES
!c     x(.,.) = coordinates of points in the population
!c     xf(.) = function values of x(.,.)
!c     xx(.) = coordinates of a single point in x
!c     cx(.,.) = coordinates of points in a complex
!c     cf(.) = function values of cx(.,.)
!c     s(.,.) = coordinates of points in the current simplex
!c     sf(.) = function values of s(.,.)
!c     bestx(.) = best point at current shuffling loop
!c     bestf = function value of bestx(.)
!c     worstx(.) = worst point at current shuffling loop
!c     worstf = function value of worstx(.)
!c     xnstd(.) = standard deviation of parameters in the population
!c     gnrng = normalized geometric mean of parameter ranges
!c     lcs(.) = indices locating position of s(.,.) in x(.,.)
!c     bound(.) = bound on ith variable being optimized
!c     ngs1 = number of complexes in current population
!c     ngs2 = number of complexes in last population
!c     iseed1 = current random seed
!c     criter(.) = vector containing the best criterion values of the last
!c         10 shuffling loops
!c
!      implicit none!REAL*8 (a-h,o-z)

!c  ARRAYS FROM THE INPUT DATA
      !!ARGUMENTS:
      TYPE(sSCE_PAR) ,intent(inout):: sPAR_SCE
      TYPE(sMEND_PAR),intent(inout):: sPAR
      TYPE(sMEND_INI),intent(inout):: sINI
      TYPE(sMEND_OUT),intent(inout):: sOUT

      !!LOCAL VARIABLES:
      INTEGER i, j, k, k1, k2, l, lpos, iCALL, ibound 
      INTEGER igs, ngs1, ngs2, nloop, loop, ipcnvg, iseed1
      INTEGER nPar, nOpt, maxn, kstop, iseed
      INTEGER ngs, npg, nps, nspl, mings, npt, npt1, iniflg, iprint
      REAL(8) fobj, fa, bestf, worstf
      REAL(8) denomi, timeou, gnrng, randx, pcento
      
      INTEGER iOpt(sPAR_SCE%nOpt)
      REAL(8) a(sPAR_SCE%nPar),bl(sPAR_SCE%nPar),bu(sPAR_SCE%nPar)
      REAL(8) xx(sPAR_SCE%nPar),bestx(sPAR_SCE%nPar),worstx(sPAR_SCE%nPar)
      REAL(8) xnstd(sPAR_SCE%nPar),bound(sPAR_SCE%nPar),criter(20),unit(sPAR_SCE%nPar)
      
!      REAL(8) x(2000,sPAR_SCE%nPar),xf(2000)
!      REAL(8) cx(2000,sPAR_SCE%nPar),cf(2000)
!      REAL(8) s(50,sPAR_SCE%nPar),sf(50)  
      REAL(8) x(sPAR_SCE%npt,sPAR_SCE%nPar),xf(sPAR_SCE%npt)  
      REAL(8) cx(sPAR_SCE%npg,sPAR_SCE%nPar),cf(sPAR_SCE%npg) 
      REAL(8) s(sPAR_SCE%nps,sPAR_SCE%nPar),sf(sPAR_SCE%nps)  
      INTEGER lcs(sPAR_SCE%nps) !!lcs(50)                     
      
      CHARACTER*256 format510,format520,format521, format610,format630,format650,format660
     
      write (*,*) '>>ENTER SUBROUTINE <SCEUA>'
  
      a     = sPAR_SCE%a
      bl    = sPAR_SCE%bl
      bu    = sPAR_SCE%bu
      nPar  = sPAR_SCE%nPar
      nOpt  = sPAR_SCE%nOpt
      iOpt  = sPAR_SCE%iOpt
      maxn  = sPAR_SCE%maxn
      kstop = sPAR_SCE%kstop
      iseed = sPAR_SCE%iseed
      ngs   = sPAR_SCE%ngs
      npg   = sPAR_SCE%npg
      npt   = sPAR_SCE%npt
      nps   = sPAR_SCE%nps
      nspl  = sPAR_SCE%nspl
      pcento = sPAR_SCE%pcento
      mings  = sPAR_SCE%mings
      iniflg = sPAR_SCE%iniflg
      iprint = sPAR_SCE%iprint
      
      
      
!      print*, sPAR_SCE%parName
            
      write(format510,*)"(/,' CRITERION',",nPar,"(a15),/1x,",&  !(6x,a4)
     &                   (nPar+1)*10,"(1h-))"
      write(format520,*)"(f10.4,",nPar,"f15.8)"

!      write(format521,*)"(",nPar,"f15.8,","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
      write(format521,*)"(",nPar,"(f15.8,','),","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
      write(format610,*)"(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,",&
     &       "'BESTF',3x,'WORSTF',3x,'PAR-RNG',1x,",nPar,"(a15))"
      write(format630,*)"(i5,1x,i5,3x,i5,3g10.3,",nPar,"(f15.8))"
!      write(format660,*)"(15x,g10.3,20x,",nPar,"(f15.8))"  
      write(format660,*)"(g10.3,",nPar,"(f15.8))"  
      write(format650,*)"(/a10,",nPar,"(a15))"

!c  INITIALIZE VARIABLES
      
      nloop = 0
      loop = 0
      igs = 0

!c  INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER

      iseed1 = iseed !-abs(iseed)

!c  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
!      npt = ngs * npg
      ngs1 = ngs
      npt1 = npt

      write(sPAR_SCE%iFout_ini,400)
      write (*,*) ' *  Evolution Loop # ',nloop
      
      write(sPAR_SCE%iFout_all,format650)'fOBJ',sPAR_SCE%parName

!c  COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
      do j = 1, nPar
        bound(j) = bu(j) - bl(j)
        unit(j) = 1.0
        xx(j) = a(j)              !wgs
      end do

      !!------------------------------------------------------------------------
        !c Check if the Initial Value of Parameter is within the defined range
        l = 0
        do j = 1, nPar
            CALL chkcst(1,a(j),bl(j),bu(j),ibound)
            if (ibound .ge. 1) then
                l = l +1
                write(*,*)"Constrains are Violated for Initial Value of Parameter: ",&
     &                      j,sPAR_SCE%parName(j),a(j), bl(j), bu(j)
            end if
        end do
        if (l.ge.1) STOP
        !!------------------------------------------------------------------------   
!c  COMPUTE THE FUNCTION VALUE OF THE INITIAL POINT
!      fa = functn(nopt,a)

      fa = fMEND_OBJ(a, sPAR, sINI, sOUT)

!c  PRINT THE INITIAL POINT AND ITS CRITERION VALUE
      write(sPAR_SCE%iFout_ini,500)
      write(sPAR_SCE%iFout_ini,format510) sPAR_SCE%parName!(xname(j),j=1,nPar)
      write(sPAR_SCE%iFout_ini,format520) fa,a !(a(j),j=1,nPar)
      
!c  GENERATE AN INITIAL SET OF npt1 POINTS IN THE PARAMETER SPACE
!c  IF iniflg IS EQUAL TO 1, SET x(1,.) TO INITIAL POINT a(.)
      
      if (iniflg .eq. 1) then
        do j = 1, nPar
          x(1,j) = a(j)
        end do
        xf(1) = fa

!c  ELSE, GENERATE A POINT RANDOMLY AND SET IT EQUAL TO x(1,.)
      else
        CALL getpnt(nPar,nOpt,iOpt,1,iseed1,xx,bl,bu,unit,bl)

        do j=1, nPar
          x(1,j) = xx(j)
        end do
        xf(1) = fMEND_OBJ(xx, sPAR, sINI, sOUT) !functn(nopt,xx)
      end if
      iCALL = 1
      if (iCALL .ge. maxn) go to 9000

!c  GENERATE npt1-1 RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
!c  SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES
      do i = 2, npt1
        CALL getpnt(nPar,nOpt,iOpt,1,iseed1,xx,bl,bu,unit,bl)
        do j = 1, nPar
          x(i,j) = xx(j)
        end do
        xf(i) = fMEND_OBJ(xx, sPAR, sINI, sOUT) !functn(nopt,xx)
        iCALL = iCALL + 1
        if (iCALL .ge. maxn) then
          npt1 = i
          go to 45
        end if
      end do

!c  ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE
!      print*, npt1
   45 CALL sort(npt1,nPar,x,xf)

!c  RECORD THE BEST AND WORST POINTS
      do j = 1, nPar
        bestx(j) = x(1,j)
        worstx(j) = x(npt1,j)
      end do
      bestf = xf(1)
      worstf = xf(npt1)

!c  COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
      CALL parstt(nPar,nOpt,iOpt,npt1,x,xnstd,bound,gnrng,ipcnvg)

!c  PRINT THE RESULTS FOR THE INITIAL POPULATION
!      write(sPAR_SCE%iFout_ini,600)
      write(sPAR_SCE%iFout_ini,format610) sPAR_SCE%parName 
      write(sPAR_SCE%iFout_ini,format630) nloop,iCALL,ngs1,bestf,worstf,gnrng,&
     &               (bestx(j),j=1,nPar)

      if (iprint .eq. 1) then
!        write(sPAR_SCE%iFout_all,650) nloop
        do i = 1, npt1
          write(sPAR_SCE%iFout_all,format660) xf(i),x(i,1:nPar)!!(x(i,j),j=1,nPar)
        end do

      end if

      if (iCALL .ge. maxn) go to 9000
      if (ipcnvg .eq. 1) go to 9200

!c  BEGIN THE MAIN LOOP ----------------
 1000 continue
      nloop = nloop + 1

      write (*,*) ' *  Evolution Loop # ',nloop

!c  BEGIN LOOP ON COMPLEXES
      do 3000 igs = 1, ngs1

!c  ASSIGN POINTS INTO COMPLEXES
      do k1 = 1, npg
        k2 = (k1-1) * ngs1 + igs
        do j = 1, nPar
          cx(k1,j) = x(k2,j)
        end do
        cf(k1) = xf(k2)
      end do

!c  BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
      do 2000 loop = 1, nspl

!c  CHOOSE A SUB-COMPLEX (nps points) ACCORDING TO A LINEAR
!c  PROBABILITY DISTRIBUTION
!!wgs[7/18/2015]: replaced by a subroutine: BEG          
!      if (nps .eq. npg) then
!        do k = 1, nps
!          lcs(k) = k
!        end do
!        go to 85
!      end if
!
!
!      randx = rand() !ran1(iseed1)
!      lcs(1) = 1 + dint(npg + 0.5 - dsqrt( (npg+.5)**2 - &
!     &         npg * (npg+1) * randx ))
!      do k = 2, nps
!   60 randx = rand()
!
!     
!        lpos = 1 + dint(npg + 0.5 - dsqrt((npg+.5)**2 - &
!     &         npg * (npg+1) * randx ))
!        do k1 = 1, k-1
!          if (lpos .eq. lcs(k1)) go to 60
!        end do
!        lcs(k) = lpos
!      end do
!!wgs[7/18/2015]: replaced by a subroutine: END 
          !!above 'lcs' generation is replaced by a subroutine below
          CALL selectINT(npg, nps, lcs)
!c  ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
      CALL sort1(nps,lcs)

!c  CREATE THE SUB-COMPLEX ARRAYS
   85 do k = 1, nps
        do j = 1, nPar
          s(k,j) = cx(lcs(k),j)
        end do
        sf(k) = cf(lcs(k))
      end do
!      write(*,'(/,100I10)')lcs(1:nps)
!      write(*,'(100F10.4)')sf(1:nps)
!c  USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
      CALL cce(sPAR_SCE, sPAR, sINI, sOUT, s,sf,xnstd,iCALL,iseed1)
!      CALL cce(nPar,nOpt,iOpt,nps,s,sf,bl,bu,xnstd,iCALL,maxn,iseed1)

!c  IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
!c  INTO THE COMPLEX
      do k = 1, nps
        do j = 1, nPar
          cx(lcs(k),j) = s(k,j)
        end do
        cf(lcs(k)) = sf(k)
      end do

!c  SORT THE POINTS
      CALL sort(npg,nPar,cx,cf)

!c  IF MAXIMUM NUMBER OF RUNS EXCEEDED, BREAK OUT OF THE LOOP
      if (iCALL .ge. maxn) go to 2222

!c  END OF INNER LOOP ------------
 2000 continue
 2222 continue

!c  REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)
      do k1 = 1, npg
        k2 = (k1-1) * ngs1 + igs
        do j = 1, nPar
          x(k2,j) = cx(k1,j)
        end do
        xf(k2) = cf(k1)
      end do
      if (iCALL .ge. maxn) go to 3333

!c  END LOOP ON COMPLEXES
 3000 continue

!c  RE-SORT THE POINTS
! print*, npt1
 3333 CALL sort(npt1,nPar,x,xf)

!c  RECORD THE BEST AND WORST POINTS
      do j = 1, nPar
        bestx(j) = x(1,j)
        worstx(j) = x(npt1,j)
      end do
      bestf = xf(1)
      worstf = xf(npt1)

!c  TEST THE POPULATION FOR PARAMETER CONVERGENCE
      CALL parstt(nPar,nOpt,iOpt,npt1,x,xnstd,bound,gnrng,ipcnvg)

!c  PRINT THE RESULTS FOR CURRENT POPULATION
!      if (mod(nloop,5) .ne. 0) go to 501
!      write(sPAR_SCE%iFout_all,format610) sPAR_SCE%parName!(xname(j),j=1,nopt)
!  501 continue
      write(sPAR_SCE%iFout_ini,format630) nloop,iCALL,ngs1,bestf,worstf,gnrng,&
     &               (bestx(j),j=1,nPar)

      if (iprint .eq. 1) then
!        write(sPAR_SCE%iFout_all,650) nloop
        do i = 1, npt1
          write(sPAR_SCE%iFout_all,format660) xf(i),x(i,1:nPar)!!(x(i,j),j=1,nPar)
        end do
      end if

!c  TEST IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED
      if (iCALL .ge. maxn) go to 9000

!c  COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
      criter(20) = bestf
      if (nloop .lt. (kstop+1)) go to 132
      denomi = dabs(criter(20-kstop) + criter(20)) / 2.
      timeou = dabs(criter(20-kstop) - criter(20)) / denomi
      if (timeou .lt. pcento) go to 9100
  132 continue
      do l = 1, 19
        criter(l) = criter(l+1)
      end do

!c  IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
      if (ipcnvg .eq. 1) go to 9200

!c  NONE OF THE STOPPING CRITERIA IS SATISFIED, CONTINUE SEARCH

!c  CHECK FOR COMPLEX NUMBER REDUCTION
      if (ngs1 .gt. mings) then
        ngs2 = ngs1
        ngs1 = ngs1 - 1
        npt1 = ngs1 * npg
        !!wgs[7/17/2015]: correction: change 'nOpt' to 'npt1'
!        CALL comp(nPar,nOpt,ngs1,ngs2,npg,x,xf,cx,cf)
        CALL comp(nPar,npt1,ngs1,ngs2,npg,x,xf,cx,cf)
! subroutine comp(n,npt,ngs1,ngs2,npg,a,af,b,bf)
      end if

!c  END OF MAIN LOOP -----------
      go to 1000

!c  SEARCH TERMINATED
 9000 continue
      write(sPAR_SCE%iFout_ini,800) maxn,loop,igs,nloop
      go to 9999
 9100 continue
      write(sPAR_SCE%iFout_ini,810) pcento*100.,kstop
      go to 9999
 9200 write(sPAR_SCE%iFout_ini,820) gnrng*100.
 9999 continue

!c  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
      write(sPAR_SCE%iFout_ini,830)
      write(sPAR_SCE%iFout_ini,format510) sPAR_SCE%parName
      write(sPAR_SCE%iFout_ini,format520) bestf,bestx  
      
      fObj = fMEND_OBJ(bestx, sPAR, sINI, sOUT)
      write(sPAR_SCE%iFout_end,format521) bestx,fobj,sINI%rOBJ
      write(*,format521) bestx,fobj,sINI%rOBJ 
      
      sPAR_SCE%bestObj = bestf
      sPAR_SCE%bestPar = bestx
      
      write (*,*) '>>EXIT SUBROUTINE <SCEUA>'
!c  END OF SUBROUTINE SCEUA
      return
  400 format(//,2x,50(1h=),/,2x,'ENTER THE SHUFFLED COMPLEX EVOLUTION',&
     &       ' GLOBAL SEARCH',/,2x,50(1h=))
  500 format(//,'*** PRINT THE INITIAL POINT AND ITS CRITERION ',&
     &       'VALUE ***')
!  510 format(/,' CRITERION',12(6x,a4),/1x,60(1h-))
!  520 format(g10.3,12f10.3)
  530 format(10x,12(6x,a4))
  540 format(10x,12f10.3)
  600 format(//,1x,'*** PRINT THE RESULTS OF THE SCE SEARCH ***')
!  610 format(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,'BEST F',3x,
!     &       'WORST F',3x,'PAR RNG',1x,8(6x,a4))
  620 format(49x,8(6x,a4))
!  630 format(i5,1x,i5,3x,i5,3g10.3,8(f10.3))
  640 format(49x,8(f10.3))
  650 format(/,1x,'POPULATION AT LOOP ',i3,/,1x,22(1h-))
!  660 format(15x,g10.3,20x,8(f10.3))
  800 format(//,1x,'*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE',&
     &       ' LIMIT ON THE MAXIMUM',/,5x,'NUMBER OF TRIALS ',i5,&
     &       ' EXCEEDED.  SEARCH WAS STOPPED AT',/,5x,'SUB-COMPLEX ',&
     &       i3,' OF COMPLEX ',i3,' IN SHUFFLING LOOP ',i3,' ***')
  810 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION',&
     &       ' VALUE HAS NOT CHANGED ',/,5x,f5.2,' PERCENT IN',i3,&
     &       ' SHUFFLING LOOPS ***')
  820 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION',&
     &       ' HAS CONVERGED INTO ',/,4x,f5.2,' PERCENT OF THE',&
     &       ' FEASIBLE SPACE ***')
  830 format(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS',&
     &       ' CRITERION VALUE ***')
     
END !!SUBROUTINE SCEUA

      
      
!!------------------------------------------------------------------
SUBROUTINE cce(sPAR_SCE, sPAR, sINI, sOUT, s,sf,xnstd,iCALL,iseed)
!(nPar,nOpt,iOpt,nps,s,sf,bl,bu,xnstd,iCALL,maxn,iseed)
!USE MOD_MEND
!USE MOD_OPT
!c$debug
!c
!c  ALGORITHM GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
!c
!c  SUB-COMPLEX VARIABLES
!      implicit REAL*8 (a-h,o-z)

!      dimension s(50,16),sf(50),bu(16),bl(16),xnstd(16)

!c  LIST OF LOCAL VARIABLES
!c    sb(.) = the best point of the simplex
!c    sw(.) = the worst point of the simplex
!c    w2(.) = the second worst point of the simplex
!c    fw = function value of the worst point
!c    ce(.) = the centroid of the simplex excluding wo
!c    snew(.) = new point generated from the simplex
!c    iviol = flag indicating if constraints are violated
!c          = 1 , yes
!c          = 0 , no
    !!ARGUMENTS:
      TYPE(sSCE_PAR) ,intent(inout):: sPAR_SCE
      TYPE(sMEND_PAR),intent(inout):: sPAR
      TYPE(sMEND_INI),intent(inout):: sINI
      TYPE(sMEND_OUT),intent(inout):: sOUT

!      parameter (c1=0.8,c2=0.4)
      INTEGER iCALL, iseed
!      REAL(8):: s(50,sPAR_SCE%nPar),sf(50),xnstd(sPAR_SCE%nPar)
      REAL(8):: s(:,:),sf(:),xnstd(:)
      
    !!LOCAL VARIABLES:  
      INTEGER :: iOpt(sPAR_SCE%nOpt)
      REAL(8) :: bu(sPAR_SCE%nPar),bl(sPAR_SCE%nPar)
      REAL(8) sw(sPAR_SCE%nPar),sb(sPAR_SCE%nPar),ce(sPAR_SCE%nPar),snew(sPAR_SCE%nPar)
      INTEGER nPar,nOpt,nps,maxn,n,m
      REAL(8) alpha,beta,fw,fnew
      INTEGER i,j,k,ibound

!c  EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
      nPar = sPAR_SCE%nPar
      nOpt = sPAR_SCE%nOpt
      iOpt = sPAR_SCE%iOpt
      nps = sPAR_SCE%nps
      maxn = sPAR_SCE%maxn
      bl = sPAR_SCE%bl
      bu = sPAR_SCE%bu
      
      n = nps
      m = nOpt
      alpha = 1.0d0
      beta = 0.5d0
      
      do i = 1, nPar
          snew(i) = s(1, i)
          sb(i) = s(1, i)
      end do
!c  IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
!c  COMPUTE THE CENTROID ce OF THE REMAINING POINTS
!c  COMPUTE step, THE VECTOR BETWEEN wo AND ce
!c  IDENTIFY THE WORST FUNCTION VALUE fw
      do k = 1, m
        j = iOpt(k)
        sb(j) = s(1,j)
        sw(j) = s(n,j)
        ce(j) = 0.0
        do i = 1, n-1
          ce(j) = ce(j) + s(i,j)
        end do
        ce(j) = ce(j)/dble(n-1)
      end do
      fw = sf(n)

!c  COMPUTE THE NEW POINT snew

!c  FIRST TRY A REFLECTION STEP
      do k = 1, m
        j = iOpt(k)
        snew(j) = ce(j) + alpha * (ce(j) - sw(j))
      end do

!c  CHECK IF snew SATISFIES ALL CONSTRAINTS
      CALL chkcst(nPar,snew,bl,bu,ibound)


!c  snew IS OUTSIDE THE BOUND,
!c  CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
!c  A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
!c  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
      if (ibound .ge. 1) CALL getpnt(nPar,nOpt,iOpt,2,iseed,snew,bl,bu,xnstd,sb)

!c  COMPUTE THE FUNCTION VALUE AT snew
      fnew = fMEND_OBJ(snew, sPAR, sINI, sOUT) !functn(nopt,snew)
      iCALL = iCALL + 1
!c
!c  COMPARE fnew WITH THE WORST FUNCTION VALUE fw
!c
!c  fnew IS LESS THAN fw, ACCEPT THE NEW POINT snew AND RETURN
      if (fnew .le. fw) go to 2000
      if (iCALL .ge. maxn) go to 3000

!c  fnew IS GREATER THAN fw, SO TRY A CONTRACTION STEP
      do k = 1, m
          j = iOpt(k)
        snew(j) = ce(j) - beta * (ce(j) - sw(j))
      end do

!c  COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
      
      fnew = fMEND_OBJ(snew, sPAR, sINI, sOUT) !functn(nopt,snew)
      iCALL = iCALL + 1
!c  COMPARE fnew TO THE WORST VALUE fw
!c  IF fnew IS LESS THAN OR EQUAL TO fw, THEN ACCEPT THE POINT AND RETURN
      if (fnew .le. fw) go to 2000
      if (iCALL .ge. maxn) go to 3000


!c  IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
!c  ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
!c  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
 1000 CALL getpnt(nPar,nOpt,iOpt,2,iseed,snew,bl,bu,xnstd,sb)

!c  COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
      fnew = fMEND_OBJ(snew, sPAR, sINI, sOUT) !functn(nopt,snew)
      iCALL = iCALL + 1


!c  REPLACE THE WORST POINT BY THE NEW POINT
 2000 continue
      do j = 1, nPar
        s(n,j) = snew(j)
      end do
      sf(n) = fnew
 3000 continue

!c  END OF SUBROUTINE CCE
      return
END !!SUBROUTINE CCE

!!------------------------------------------------------------------
SUBROUTINE getpnt(nPar,nOpt,iOpt,idist,iseed,x,bl,bu,std,xi)
!c$debug      
!c
!c     This subroutine generates a new point within feasible region
!c
!c     x(.) = new point
!c     xi(.) = focal point
!c     bl(.) = lower bound
!c     bu(.) = upper bound
!c     std(.) = standard deviation of probability distribution
!c     idist = probability flag
!c           = 1 - uniform distribution
!c           = 2 - Gaussian distribution
    
    !!ARGUMENTS:
      INTEGER nPar,nOpt,idist,iseed
      INTEGER iOpt(nOpt)
      REAL(8) x(nPar),bl(nPar),bu(nPar),std(nPar),xi(nPar)
      
    !!LOCAL VARIABLES:  
      INTEGER k,j,ibound
      REAL(8) randx

    1 do k=1, nOpt
        j = iOpt(k)  
    2   if (idist .eq. 1) randx = rand() !ran1(iseed)
        if (idist .eq. 2) then
            randx = gasdev(iseed)
!            print*,">>>Gaussian=",randx,iseed
        end if
        x(j) = xi(j) + std(j) * randx * (bu(j) - bl(j))

!c     Check explicit constraints
        
        CALL chkcst(1,x(j),bl(j),bu(j),ibound)
        if (ibound .ge. 1) go to 2
      end do

!c     Check implicit constraints
      
      CALL chkcst(nPar,x,bl,bu,ibound)
      if (ibound .ge. 1) go to 1

      return
END 

!!------------------------------------------------------------------
SUBROUTINE parstt(nPar,nOpt,iOpt,npt,x,xnstd,bound,gnrng,ipcnvg)
!c$debug      
!c
!c  SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
    !!ARGUMENTS:
      INTEGER nPar,nOpt,npt,ipcnvg
      INTEGER:: iOpt(nOpt)
      REAL(8):: gnrng
!      REAL(8):: x(2000,nPar),xnstd(nPar),bound(nPar)
      REAL(8):: x(:,:),xnstd(nPar),bound(nPar)
      
      !!LOCAL VARIABLES:
      INTEGER i,j,k
      REAL(8), parameter:: delta = 1.0d-20,peps=1.0d-3
      REAL(8):: xmean(nPar),xmax(nPar),xmin(nPar)
      REAL(8) gsum,xsum1,xsum2

!c  COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
      gsum = 0.d0
      do j = 1, nOpt
        k = iOpt(j)
        xmax(k) = -1.0d+20
        xmin(k) = 1.0d+20
        xsum1 = 0.d0
        xsum2 = 0.d0
        do i = 1, npt
          xmax(k) = dmax1(x(i,k), xmax(k))
          xmin(k) = dmin1(x(i,k), xmin(k))
          xsum1 = xsum1 + x(i,k)
          xsum2 = xsum2 + x(i,k)*x(i,k)
        end do
        xmean(k) = xsum1 / dble(npt)
        xnstd(k) = (xsum2 / dble(npt) - xmean(k)*xmean(k))
        if (xnstd(k) .le. delta) xnstd(k) = delta
        xnstd(k) = dsqrt(xnstd(k))
        xnstd(k) = xnstd(k) / bound(k)
        gsum = gsum + dlog( delta + (xmax(k)-xmin(k))/bound(k) )
      end do
      gnrng = dexp(gsum/dble(nOpt))

!c  CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
      ipcnvg = 0
      if (gnrng .le. peps) then
        ipcnvg = 1
      end if

!c  END OF SUBROUTINE PARSTT
      return
END

      
      
!!------------------------------------------------------------------
SUBROUTINE comp(n,npt,ngs1,ngs2,npg,a,af,b,bf)
!c$debug      
!c
!c
!c  THIS SUBROUTINE REDUCE INPUT MATRIX a(n,ngs2*npg) TO MATRIX
!c  b(n,ngs1*npg) AND VECTOR af(ngs2*npg) TO VECTOR bf(ngs1*npg)
!      implicit REAL*8 (a-h,o-z)
    !!ARGUMENTS:
!      REAL(8):: a(2000,n),af(2000),b(2000,n),bf(2000)
      REAL(8):: a(:,:),af(:),b(:,:),bf(:)
      INTEGER n,npt,ngs1,ngs2,npg
    
      !!LOCAL VARIABLES:
      INTEGER igs,ipg,k1,k2,j,i
!      dimension a(2000,16),af(2000),b(2000,16),bf(2000)
      do igs=1, ngs1
        do ipg=1, npg
          k1=(ipg-1)*ngs2 + igs
          k2=(ipg-1)*ngs1 + igs
          do i=1, n
            b(k2,i) = a(k1,i)
          end do
          bf(k2) = af(k1)
        end do
      end do

      do j=1, npt
        do i=1, n
          a(j,i) = b(j,i)
        end do
        af(j) = bf(j)
      end do

!c  END OF SUBROUTINE COMP
      return
END
      
!!------------------------------------------------------------------
SUBROUTINE chkcst(nPar,x,bl,bu,ibound)
!c
!c     This subroutine check if the trial point satisfies all
!c     constraints.
!c
!c     ibound - violation indicator
!c            = -1 initial value
!c            = 0  no violation
!c            = 1  violation
!c     nopt = number of optimizing variables
!c     ii = the ii'th variable of the arrays x, bl, and bu
!c
    !!ARGUMENTS:
      INTEGER nPar,ibound
      REAL(8) x(nPar),bl(nPar),bu(nPar)
      
      !!LOCAL VARIABLES:
      INTEGER ii

      ibound = -1

!c     Check if explicit constraints are violated

      do ii=1, nPar
        if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) then
!            write(*,*)"Constrains are violated for Parameter: ",ii
            go to 10
        end if
      end do
      if (nPar .eq. 1) go to 9

!c     Check if implicit constraints are violated
!c     (no implicit constraints for this function)
!c
!c     No constraints are violated
!c      
    9 ibound = 0
      return

!c     At least one of the constraints are violated
     
   10 ibound = 1
   
!      do k = 1, nPar
!          print*, k, x(k), bl(k), bu(k)
!      end do
!      print*, ibound
    
      return
  END


END MODULE MOD_OPT     
!!c==============================================================
!      real*8 function ran1(idum)
!!c$debug
!!c
!!c
!!c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!      implicit real*8 (a-h,o-z)
!      dimension r(97)
!      common /rancom/ ix1,ix2,ix3
!      parameter (m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 =
!     &3.8580247e-6)
!      parameter (m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 =
!     &7.4373773e-6)
!      parameter (m3 = 243000, ia3 = 4561, ic3 = 51349)
!      data iff / 0 /
!      if ((idum .lt. 0) .or. (iff .eq. 0)) then
!      iff = 1
!      ix1 = mod(ic1 - idum,m1)
!      ix1 = mod((ia1 * ix1) + ic1,m1)
!      ix2 = mod(ix1,m2)
!      ix1 = mod((ia1 * ix1) + ic1,m1)
!      ix3 = mod(ix1,m3)
!      do 11 j = 1, 97
!      ix1 = mod((ia1 * ix1) + ic1,m1)
!      ix2 = mod((ia2 * ix2) + ic2,m2)
!      r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
!   11 continue
!      idum = 1
!      end if
!      ix1 = mod((ia1 * ix1) + ic1,m1)
!      ix2 = mod((ia2 * ix2) + ic2,m2)
!      ix3 = mod((ia3 * ix3) + ic3,m3)
!      j = 1 + ((97 * ix3) / m3)
!      if ((j .gt. 97) .or. (j .lt. 1)) pause
!!      if ((j .gt. 97) .or. (j .lt. 1)) j = 97
!      ran1 = r(j)
!      r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
!
!!c  END OF SUBROUTINE RAN1
!      return
!      end
!!----------------------------------------------------------------------      
!      subroutine indexx(n, arrin, indx)
!!c$debug      
!!c
!!c
!!c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!         
!      implicit real*8 (a-h,o-z)
!      dimension arrin(n), indx(n)
!      print*,"indexx = ", n
!      do 11 j = 1, n
!        indx(j) = j
!   11 continue
!      l = (n / 2) + 1
!      ir = n
!   10 continue
!      if (l .gt. 1) then
!        l = l - 1
!        indxt = indx(l)
!        q = arrin(indxt)
!      else
!        indxt = indx(ir)
!        q = arrin(indxt)
!        indx(ir) = indx(1)
!        ir = ir - 1
!        if (ir .eq. 1) then
!            indx(1) = indxt
!            return
!        end if
!      end if
!      i = l
!      j = l + l
!   20 if (j .le. ir) then
!        if (j .lt. ir) then
!            if (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
!        end if
!        if (q .lt. arrin(indx(j))) then
!            indx(i) = indx(j)
!            i = j
!            j = j + j
!        else
!            j = ir + 1
!        end if
!        goto 20
!      end if
!      indx(i) = indxt
!      goto 10
!
!!c  END OF SUBROUTINE INDEXX
!      end

!!------------------------------------------------------------------
!SUBROUTINE sort(n,m,rb,ra)
!!c$debug      
!!c
!!c
!!c  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!!c  BY W.H. PRESS ET AL., pp. 233-234
!!c
!!c  LIST OF VARIABLES
!!c     ra(.) = array to be sorted
!!c     rb(.,.) = arrays ordered corresponding to rearrangement of ra(.)
!!c     wk(.,.), iwk(.) = local varibles
!
!    !!ARGUMENTS:
!      INTEGER n,m
!      REAL(8) ra(2000),rb(2000,m)
!      
!      !!LOCAL VARIABLES
!      REAL(8) wk(2000,m)
!      INTEGER iwk(2000)
!      INTEGER i,j
!!      dimension ra(n),rb(2000,m),iwk(n)
!!      REAL(8), allocatable:: wk(:,:)
!!      allocate(wk(n, m))
!      CALL indexx(n, ra, iwk)
!      do 11 i = 1, n
!      wk(i,1) = ra(i)
!   11 continue
!      do 12 i = 1, n
!      ra(i) = wk(iwk(i),1)
!   12 continue
!      do 14 j = 1, m
!      do 13 i = 1, n
!      wk(i,j) = rb(i,j)
!   13 continue
!   14 continue
!      do 16 j = 1, m
!      do 15 i = 1, n
!      rb(i,j) = wk(iwk(i),j)
!   15 continue
!   16 continue
!
!!c  END OF SUBROUTINE SORT
!      return
!END
!      
!!!------------------------------------------------------------------
!SUBROUTINE sort1(n,ira)
!!c$debug      
!!c
!!c
!!c  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!!c  BY W.H. PRESS ET AL., pp. 231
!!c
!!c  LIST OF VARIABLES
!!c     ira(.) = INTEGER array to be sorted
!    !!ARGUMENTS:
!      INTEGER n
!      INTEGER:: ira(n)
!      
!    !!LOCAL VARIABLES:  
!      INTEGER l,ir,i,j
!!      dimension ra(n)
!
!!      INTEGER ra, rra
!      INTEGER irra
!
!      l = (n / 2) + 1
!      ir = n
!   10 continue
!      if (l .gt. 1) then
!      l = l - 1
!      irra = ira(l)
!      else
!      irra = ira(ir)
!      ira(ir) = ira(1)
!      ir = ir - 1
!      if (ir .eq. 1) then
!      ira(1) = irra
!      return
!      end if
!      end if
!      i = l
!      j = l + l
!   20 if (j .le. ir) then
!      if (j .lt. ir) then
!      if (ira(j) .lt. ira(j + 1)) j = j + 1
!      end if
!      if (irra .lt. ira(j)) then
!      ira(i) = ira(j)
!      i = j
!      j = j + j
!      else
!      j = ir + 1
!      end if
!      goto 20
!      end if
!      ira(i) = irra
!      goto 10
!
!!c  END OF SUBROUTINE SORT1
!END
!      
!!!------------------------------------------------------------------
!SUBROUTINE indexx(n,arr,indx)
!    !!ARGUMENTS:
!    INTEGER n
!    INTEGER indx(n)
!    REAL(8) arr(n)
!    
!    !!LOCAL VARIABLES:
!    INTEGER, PARAMETER :: M = 7
!    INTEGER, PARAMETER :: NSTACK=50
!!        Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that arr(indx(j))
!!        is in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr are not changed.
!    INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
!    REAL(8) a
!    do 11 j=1,n
!        indx(j)=j
! 11 continue !! end do !!11
!    jstack=0
!    l=1
!    ir=n
!  1 if(ir-l.lt.M)then
!    do 13 j=l+1,ir
!    indxt=indx(j)
!    a=arr(indxt)
!    do 12 i=j-1,l,-1
!        if(arr(indx(i)).le.a)goto 2
!        indx(i+1)=indx(i)
! 12 continue !! end do !12
!    i=l-1
!  2 indx(i+1)=indxt
! 13 continue !!  end do !!13
!    if(jstack.eq.0) return
!    ir=istack(jstack)
!    l=istack(jstack-1)
!    jstack=jstack-2
!    else
!    k=(l+ir)/2
!    itemp=indx(k)
!    indx(k)=indx(l+1)
!    indx(l+1)=itemp
!    if(arr(indx(l)).gt.arr(indx(ir)))then
!        itemp=indx(l)
!        indx(l)=indx(ir)
!        indx(ir)=itemp
!    end if
!    if(arr(indx(l+1)).gt.arr(indx(ir)))then
!        itemp=indx(l+1)
!        indx(l+1)=indx(ir)
!        indx(ir)=itemp
!    end if
!    if(arr(indx(l)).gt.arr(indx(l+1)))then
!        itemp=indx(l)
!        indx(l)=indx(l+1)
!        indx(l+1)=itemp
!    end if
!    i=l+1
!    j=ir
!    indxt=indx(l+1)
!    a=arr(indxt)
!  3 continue
!    i=i+1
!    if(arr(indx(i)).lt.a)goto 3
!  4 continue
!    j=j-1
!    if(arr(indx(j)).gt.a)goto 4
!    if(j.lt.i)goto 5
!    itemp=indx(i)
!    indx(i)=indx(j)
!    indx(j)=itemp
!    goto 3
!  5 indx(l+1)=indx(j)
!    indx(j)=indxt
!    jstack=jstack+2
!    if(jstack.gt.NSTACK) print*, "NSTACK too small in indexx"
!
!    if(ir-i+1.ge.j-l)then
!        istack(jstack)=ir
!        istack(jstack-1)=i
!        ir=j-1
!    else
!        istack(jstack)=j-1
!        istack(jstack-1)=l
!        l=i
!    end if
!    end if
!    goto 1
!END 
!

