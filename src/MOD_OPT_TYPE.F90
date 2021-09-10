MODULE MOD_OPT_TYPE       
! File:   STRUCT_OPT.F90
! Author: GANGSHENG WANG @ ORNL
! Updated: May 5, 2015
! Created on March 5, 2013, 4:06 PM
!

!-----------------------------------------------------------------------------!   
!SCE-UA Algorithm
    TYPE sSCE_PAR
        INTEGER nPar                                    !# of parameters
        INTEGER nOpt                                    !# of optimized parameters
        INTEGER nRun                                    !# of SCEUA runs, each with a unique random-seed
        INTEGER maxn                                    !max no. of trials allowed before optimization is terminated
        INTEGER kstop                                   !number of shuffling loops in which the criterion value must change by the given percentage before optimization is terminated
        INTEGER ngs                                     !number of complexes in the initial population
        INTEGER npt                                     !total number of points in initial population (npt=ngs*npg)
        INTEGER npg                                     !number of points in each complex
        INTEGER nps                                     !number of points in a sub-complex
        INTEGER nspl                                    !number of evolution steps allowed for each complex before complex shuffling
        INTEGER mings                                   !minimum number of complexes required, if the number of complexes is allowed to reduce as the optimization proceeds
        INTEGER ideflt                                  !dIF ideflt IS EQUAL TO 0, SET THE SCE CONTROL PARAMETERS TO THE DEFAULT VALUES
        INTEGER iniflg                                  !flag on whether to include the initial point in population; = 0, not included;  = 1, included
        INTEGER iprint                                  !flag for controlling print-out after each shuffling loop;= 0, print information on the best point of the population; = 1, print information on every point of the population
        INTEGER iseed                                   !random seed
        INTEGER iFout_all                                  !output file1: print the results of SCEUA
        INTEGER iFout_end                                  !output file2: summarize optimal parameter values from nRun
        INTEGER iFout_ini                                  !output file3: other text output
        REAL(8) pcento                                  !percentage by which the criterion value must change in given number of shuffling loops
        REAL(8) bestObj                                 !best objective function value
        INTEGER, ALLOCATABLE:: iOpt(:)                  !index of optimized parameters
        REAL(8), ALLOCATABLE:: a(:)                     !initial values
        REAL(8), ALLOCATABLE:: bl(:)                    !lower bound values
        REAL(8), ALLOCATABLE:: bu(:)                    !upper bound values
        REAL(8), ALLOCATABLE:: bestPar(:)               !best parameter values (nPar)
        CHARACTER(LEN = 10), ALLOCATABLE:: parName(:)   !parameter names
    END TYPE sSCE_PAR
!-----------------------------------------------------------------------------!   
END MODULE MOD_OPT_TYPE  

