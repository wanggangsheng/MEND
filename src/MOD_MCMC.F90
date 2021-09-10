MODULE MOD_MCMC
! File:   MOD_MCMC.F90
! Author: GANGSHENG WANG @ OU
! Updated: December 28, 2018
! Created on February 26, 2013, 11:02 AM    
USE MOD_OPT_TYPE
USE MOD_MEND
!USE MOD_USRFS, ONLY: sort,sort1,indexx,gasdev,selectINT

    IMPLICIT NONE
!    PRIVATE :: cce
    PRIVATE :: getpnt
    PRIVATE :: parstt
!    PRIVATE :: comp
    PRIVATE :: chkcst
    
    PUBLIC :: MCMC

    
    
    CONTAINS
    
    
SUBROUTINE MCMC(sPAR_SCE, sPAR, sINI, sOUT)

!    (a,bl,bu,nopt,maxn,kstop,pcento,iseed,ngs,npg,nps,nspl,mings,iniflg,iprint)
!c$debug
!c
!c  by GANGSHENG WANG
!C  University of Oklahoma
!C  December, 2018



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
!c     maxn = max no. of trials allowed 
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
      INTEGER i, j, iseed1, kUpdate  !!j, k, k1, k2, l, lpos, iCALL,
!      INTEGER igs, ngs1, ngs2, nloop, loop, ipcnvg, kstop
      INTEGER nPar, nOpt, maxn, iseed, iniflg, iprint
!      INTEGER ngs, npg, nps, nspl, mings, npt, npt1
      REAL(8) alpha, randx, acceptance_rate
!      REAL(8) fobj, fa, bestf, worstf
!      REAL(8) denomi, timeou, gnrng, randx, pcento
      
      INTEGER iOpt(sPAR_SCE%nOpt)
      REAL(8) a(sPAR_SCE%nPar),bl(sPAR_SCE%nPar),bu(sPAR_SCE%nPar)
      REAL(8) xx(sPAR_SCE%nPar),bestx(sPAR_SCE%nPar), xx_new(sPAR_SCE%nPar)  
      REAL(8) unit(sPAR_SCE%nPar)  !!xnstd(sPAR_SCE%nPar),bound(sPAR_SCE%nPar),criter(20),
            
      REAL(8) SSE_old, SSE_new, SSE_diff,SSE_min
!      
      CHARACTER*256 format510,format520,format521, format530, format610,format630,format660
     
      write (*,*) '>>ENTER SUBROUTINE <MCMC>'
  
      a = sPAR_SCE%a
      bl = sPAR_SCE%bl
      bu = sPAR_SCE%bu
      nPar = sPAR_SCE%nPar
      nOpt = sPAR_SCE%nOpt
      iOpt = sPAR_SCE%iOpt
      maxn = sPAR_SCE%maxn
      
      iseed = sPAR_SCE%iseed
      
!      kstop = sPAR_SCE%kstop
!      pcento = sPAR_SCE%pcento
!      ngs = sPAR_SCE%ngs
!      npg = sPAR_SCE%npg
!      npt = sPAR_SCE%npt
!      nps = sPAR_SCE%nps
!      nspl = sPAR_SCE%nspl
!      mings  = sPAR_SCE%mings
      iniflg  = sPAR_SCE%iniflg
      iprint = sPAR_SCE%iprint
      
      
!      print*, sPAR_SCE%parName
            
      write(format510,*)"(",(nPar+4),"(a15)",")"
      write(format520,*)"(f10.4,",nPar,"f15.8)"
      write(format530,*)"(a20,",nPar,"g10.3)"

!      write(format521,*)"(",nPar,"f15.8,","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
      write(format521,*)"(",nPar,"(f15.8,','),","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
      write(format610,*)"(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,",&
     &       "'BESTF',3x,'WORSTF',3x,'PAR-RNG',1x,",nPar,"(a15))"
      write(format630,*)"(i5,1x,i5,3x,i5,3g10.3,",nPar,"(f15.8))"
      write(format660,*)"(4g15.4,",nPar,"(f15.8))"    

!c  INITIALIZE VARIABLES
      
!      nloop = 0
!      loop = 0
!      igs = 0

!c  INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER
      
      write(sINI%iFout_mcmc,format510)"alpha","randx", "SSE_old", "SSE_new",sPAR_SCE%parName 

      iseed1 = iseed !-abs(iseed)
      unit = 1.0D-1
!      SSE_old = 1.0D+6
      xx = a
      
      if (iniflg .eq. 0) then
          CALL getpnt(nPar,nOpt,iOpt,1,iseed1,xx,bl,bu,unit,bl)
      end if
!      write(*,format530)"xx0 = ",xx
      
      SSE_old = fMEND_OBJ(xx, sPAR, sINI, sOUT)
      
      kUpdate = 0
      SSE_min = SSE_old
      do i = 1,maxn
!          write(*,*)">>>MCMC I = ", i
          xx_new = xx
          CALL getpnt(nPar,nOpt,iOpt,1,iseed1,xx_new,bl,bu,unit,xx)
!          write(*,format530)"xx_new = ",xx_new
          SSE_new = fMEND_OBJ(xx_new, sPAR, sINI, sOUT)
          SSE_diff = SSE_new - SSE_old
          
!          write(*,*)"SSE= ",SSE_old,SSE_new,SSE_diff
!          write(*,*)"exp= ",dexp(-SSE_diff)
          
          alpha = min(1.d0,dexp(-SSE_diff))
          randx = rand()
          if(alpha > randx) then
              kUpdate = kUpdate + 1
              if(SSE_new < SSE_min) then
                  SSE_min = SSE_new
                  sPAR_SCE%bestPar = xx_new
              end if
              
              write(sINI%iFout_mcmc,format660) alpha,randx, SSE_old, SSE_new, (xx_new(j),j=1,nPar)
              
              SSE_old = SSE_new
!              bestx = xx_new
              xx = xx_new
              write(*,'(A10,I6,A15,I6,A4,I6)')">>>MCMC: ",maxn, ": k-Update = ", kUpdate, " of ",i
   
          end if
          
      end do  !!i = 1,maxn

!c  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTIO

!c  PRINT THE RESULTS FOR THE INITIAL POPULATION


!
!      fObj = fMEND_OBJ(bestx, sPAR, sINI, sOUT)
!      write(sPAR_SCE%iFout2,format521) bestx,fobj,sINI%rOBJ
!      write(*,format521) bestx,fobj,sINI%rOBJ 
!      
!      sPAR_SCE%bestObj = bestf
!      sPAR_SCE%bestPar = bestx
      acceptance_rate = kUpdate*1.d0/maxn*100
      write (*,'(A25,f6.2,A1)')"MCMC Acceptance Rate = ",acceptance_rate,"%"
      write (*,*) '>>EXIT SUBROUTINE <MCMC>'
!c  END OF SUBROUTINE SCEUA
      return

END !!SUBROUTINE SCEUA



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
      CHARACTER*256 format530
      
      write(format530,*)"(a10,",nPar,"g10.3)"

    1 do k=1, nOpt
        j = iOpt(k)  
    2   if (idist .eq. 1) randx = rand() !ran1(iseed)
        if (idist .eq. 2) then
            randx = gasdev(iseed)
!            print*,">>>Gaussian=",randx,iseed
        end if
        x(j) = xi(j) + std(j) * (randx - 0.5) * (bu(j) - bl(j))

!c     Check explicit constraints
        
        CALL chkcst(1,x(j),bl(j),bu(j),ibound)
!        write(*,*)"getpnt: ",k,j,x(j),bl(j),bu(j),ibound
        if (ibound .ge. 1) go to 2
      end do

!c     Check implicit constraints
      
      CALL chkcst(nPar,x,bl,bu,ibound)
      if (ibound .ge. 1) then
          go to 1
      end if

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
        if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) go to 10
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


END MODULE MOD_MCMC  
!!c==============================================================
