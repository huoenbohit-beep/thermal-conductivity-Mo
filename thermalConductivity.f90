!huoenbo
!Mo



    program thermalConductivity
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      compute the thermal conductivity according to the method in 
!*      "An estimate for thermal diffusivity in highly irradiated tungsten using MolecularDynamics simulation"
!*      D.R. Mason et al (2021)
!*
!*      The most important input is the potential energies, which should be in a simple text file
!*      line 1 :        number_of_atoms                                                                                         
!*      line 2+:        e1 e2 e3 ...
!*

        use Lib_CommandLineArguments            !   used to read command line arguments and produce a simple -h option
        use iso_fortran_env
        implicit none



    !---    input parameters - note that these initial values are for tungsten, in order to help interpret the paper.
        character(len=256)              ::      filename = ""               !   input file containing energies
        
        real(kind=real64)               ::      b0           = 2.724d0     !   nearest neighbour separation                     ( A )
        real(kind=real64)               ::      e0           = -6.82d0      !   energy per atom perfect lattice, inc thermal expansion, exc 3/2 kT ( eV )
        real(kind=real64)               ::      Delta        = 0.029d0      !   scaling of Maxwell-Boltzmann energy spread with temperature sigma^2 = Delta kB T (eV)
        
        real(kind=real64)               ::      S_0          = 2.5d0       !   impurity scattering rate r_imp = S0 |E| ( 1/fs/eV )
        real(kind=real64)               ::      S_1          = 1.9d-4     !   e-ph scattering rate r_eph = S1 T       ( 1/fs/K )
        real(kind=real64)               ::      S_2          = 1.706d-7     !   e-e scattering rate ree = S2 T^2        ( 1/fs/K^2 )           
        real(kind=real64)               ::      ceonT        = 1.49d-8    !   electronic heat capacity per atom       ( ev/K^2 )
        real(kind=real64)               ::      vF           = 8.72d0       !   Fermi velocity                          ( A/fs )
        
real(kind=real64) :: T = 10.0d0
        logical                         ::      MDsim        = .true.       !   true if we expect thermal variation 3/2 kT, false if CG data local minimum 
        real(kind=real64)               ::      EMIN = -0.5d0 , EMAX = 2.0d0    !   maximum and minimum excess energy per atom eV
        integer                         ::      NBIN = 700                      !   number of bins        
        logical                         ::      opHist       = .true.      !   output the full histogram?        
        
        type(CommandLineArguments)      ::      cla
        
        
        
    !---    physical constants
        real(kind=real64),parameter     ::      KB = 0.00008617d0           !   Boltzmann's constant (eV/K)    
        real(kind=real64),parameter     ::      EVFSATOWM = 1.602d6         !   eV/fs 1/A 1/K in W/m/K
        real(kind=real128)              ::      sqrt2onPI = sqrt( 0.5_real128/atan(1.0_real128) )
        real(kind=real64),parameter     ::      PI = 3.141592653590d0       !    
        real(kind=real64),parameter     ::      L = 2.44d-8                 !   Lorentz constant 
        
    !---    properties read from input file
        integer                         ::      nAtoms                      !   number of atoms    
        real(kind=real64),dimension(:),allocatable  ::      xe              !   array storing potential energy of atoms
        
        
        
    !---    derived constants
        real(kind=real64)               ::      kappa_const                 !   thermal conductivity constant 1/3 (ce/T) T vF^2/Omega0
        real(kind=real64)               ::      deltaE                      !   histogram energy window 
        real(kind=real64)               ::      sigma                       !   Maxwell-Boltzmann broadening
        
        
    !---    dummy variables
        logical                         ::      ok 
        integer                         ::      ii,ioerr,bb
        real(kind=real64)               ::      ee,kk,rr,e2Bar
        
        
    !---    output variables
        real(kind=real64)                   ::      nDefects        !   statistically determined count of atoms with athermal pot eng
        real(kind=real64)                   ::      rBar            !   mean scattering rate
        real(kind=real64)                   ::      eBar            !   mean excess potential energy
        real(kind=real64)                   ::      kappa           !   thermal conductivity x volume per atom   
        real(kind=real64)                   ::      alpha           !   thermal diffusivity
        real(kind=real64)                   ::      rho             !   electrical resisitivity
        real(kind=real64)                   ::      eVar            !   variance of potential energy
        integer,dimension(:),allocatable    ::      hist            !   histogram of atom potential energy counts
        
        
        
        
    !---    read in command line arguments
        cla = CommandLineArguments_ctor(20)
        call setProgramDescription( cla, "thermalConductivity.exe - computes thermal conductivity from potential energies &             
                                          \n     Note that kappa needs to be divided by volume per atom (in A^3) to give answer in standard units       &
                                          \n     This keeps thermal diffusivity independent of simulation volume." )
        call setProgramVersion( cla, "1.0.0" )
        call setCategories(cla,(/ "required    ","potential   ","simulation  ","scattering  ","histogram   " /))
        
        call get( cla,"f"    ,filename     ,LIB_CLA_REQUIRED,"      input filename                                                                 ",1 )  
        call get( cla,"Delta",Delta        ,LIB_CLA_OPTIONAL,"scaling of Maxwell-Boltzmann energy spread with temperature sigma^2 = Delta kB T (eV)",2 )  
        call get( cla,"b0"   ,b0           ,LIB_CLA_OPTIONAL,"   nearest neighbour spacing               ( A )                                     ",2 )  
        call get( cla,"e0"   ,e0           ,LIB_CLA_OPTIONAL,"   energy per atom perfect lattice, inc thermal expansion, exc 3/2 kT ( eV )         ",2 )  
        call get( cla,"T"    ,T            ,LIB_CLA_OPTIONAL,"    temperature                             (K)                                      ",3 )  
        call get( cla,"MD"   ,MDsim        ,LIB_CLA_OPTIONAL,"   true if we expect thermal variation 3/2 kT, false if CG data local minimum        ",3 )  
        call get( cla,"S_0"  ,S_0          ,LIB_CLA_OPTIONAL,"  impurity scattering rate r_imp = S0 |E| ( 1/fs/eV )                                ",4 )  
        !call get( cla,"S_1"  ,S_1          ,LIB_CLA_OPTIONAL,"  e-ph scattering rate r_eph = S1 T       ( 1/fs/K )                                 ",4 )  
        !call get( cla,"S_2"  ,S_2          ,LIB_CLA_OPTIONAL,"  e-e scattering rate ree = S2 T^2        ( 1/fs/K^2 )                               ",4 )  
        call get( cla,"ceonT",ceonT        ,LIB_CLA_OPTIONAL,"electronic heat capacity constant       ( ev/K^2 )                                   ",4 )  
        call get( cla,"vF"   ,vF           ,LIB_CLA_OPTIONAL,"   Fermi velocity                          ( A/fs )                                  ",4 )  
        call get( cla,"EMIN" ,EMIN         ,LIB_CLA_OPTIONAL," minimum excess energy per atom          ( eV )                                      ",5 )  
        call get( cla,"EMAX" ,EMAX         ,LIB_CLA_OPTIONAL," maximum excess energy per atom          ( eV )                                      ",5 )  
        call get( cla,"NBIN" ,NBIN         ,LIB_CLA_OPTIONAL," number of histogram bins                                                            ",5 )  
        call get( cla,"hist" ,opHist       ,LIB_CLA_OPTIONAL," output the full histogram?                                                          ",5 )  
        
        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
        
        
        
        
    !---    compute derived constants from input variables
        kappa_const = (ceonT/3)*T*vF*vF 
        sigma = sqrt( Delta*KB*T )
        deltaE = (EMAX - EMIN)/NBIN
        
        
        
    !---    output a friendly message about the input parameters
        print *,""
        print *,"thermalConductivity.exe info - input parameters"
        print *,"file  : ",filename    
        print *,"Delta : ",Delta       
        print *,"b0    : ",b0          
        print *,"e0    : ",e0          
        print *,"T     : ",T           
        print *,"MD    : ",MDsim       
        print *,"S_0   : ",S_0         
        print *,"S_1   : ",S_1         
        print *,"S_2   : ",S_2         
        print *,"ceonT : ",ceonT       
        print *,"vF    : ",vF          
        print *,"EMIN  : ",EMIN        
        print *,"EMAX  : ",EMAX        
        print *,"NBIN  : ",NBIN      
        print *,""  
        
        
        

                print *,"[number of atoms]"
                print *,"e1 e2 e3 ... eN"  
                stop
            end if                
        close(unit=400)               
        
    !   remove expected potential energy per atom, in order to find _excess_ energy only   
        xe = xe - e0            
        
        
    !---    construct histogram 
        allocate(hist(0:NBIN))
        hist = 0                   
        eBar = 0.0d0
        e2Bar = 0.0d0
        do ii = 1,nAtoms
            ee = xe(ii)                             !   energy of atom i                                
            bb = nint( (ee - EMIN)/deltaE )         !   find correct bin for energy
            bb = max( 0, min( NBIN,bb ) )           !   limit bin range
            hist(bb) = hist(bb) + 1
            eBar = eBar + ee
            e2Bar = e2Bar + ee*ee
        end do
        eBar = eBar / nAtoms
        e2Bar = e2Bar / nAtoms
        print *,"thermalConductivity.exe info - number of atoms in MIN bin (",EMIN,") = ",hist(0)
        print *,"thermalConductivity.exe info - number of atoms in MAX bin (",EMAX,") = ",hist(NBIN)
        
        
        
        
    !---    compute number of defects and average scattering rate
        print *,""
        if (opHist) write(*,fmt='(5a16)')   "energy","<n>","n","n_defect","r"
        do bb = 0,NBIN
            ee = EMIN + deltaE*bb
            
            
            call computeDefectsFromHist( hist(bb),nAtoms,ee,deltaE,T, sigma, b0,vF,S_0,S_1,S_2 , kk,rr )
            
            nDefects = nDefects + kk
            rBar = rBar + rr
            
            if (opHist) then
                rho = BroadenedMaxwellBoltzmann( ee,T,sigma )*deltaE*nAtoms 
                write(*,fmt='(2f16.6,i16,2f16.6)') ee,rho,hist(bb),kk,rr
            end if
            
        end do
        rBar = rBar/nAtoms
        if (opHist) print *,""
        
    !---    compute thermal conductivity x volume per atom - divide by Omega0 to get "correct" units
        kappa = kappa_const / rBar
        
    !---    compute thermal diffusivity
        alpha = kappa /(3*kB)
        
    !---    compute excess potential energy variance
        eVar = max(0.0d0, e2Bar - eBar*eBar)            !   max function to avoid rounding errors.
        
    !---    compute resisitivity
        rho = L * T / kappa
        
        
    !---    output answer
        print *,""
        write(*,fmt='(3(a,i16))')   "atom count             ",nAtoms
        
        write(*,fmt='(3(a,g16.6))') "mean excess energy     ",eBar ,"                   +/- ",sqrt( eVar )," eV"
        if (MDsim) &            
        write(*,fmt='(3(a,g16.6))') "thermal energy 3/2 kT  ",1.5d0*kB*T," eV,           delta = ",eBar-1.5d0*kB*T," eV"
        
        write(*,fmt='(3(a,g16.6))') "excess energy variance ",eVar ," eV^2 "
        write(*,fmt='(3(a,g16.6))') "mean scattering rate   ",rBar ," 1/fs"
        !write(*,fmt='(3(a,g16.6))') "therm. cond. x Omega0  ",kappa," eV/fs A^2 1/K = ",kappa*EVFSATOWM," Omega0 W/m/K "
        write(*,fmt='(3(a,g16.6))') "therm. cond. x Omega0  ",kappa," eV/fs 1/A 1/K (A^3) = ",kappa*EVFSATOWM," W/m/K (A^3)"
        
        
        
        write(*,fmt='(3(a,g16.6))') "resistivity / Omega0   ",rho/EVFSATOWM," Ohm m (1/A^3) "
        
        
        write(*,fmt='(3(a,g16.6))') "thermal diffusivity    ",alpha," A^2/fs              = ",alpha*1.0d-5   ," m^2/s"
        write(*,fmt='(3(a,g16.6))') "number of defects      ",nDefects,"                     = ",nDefects*100.0/nAtoms   ,"%"
        
        
        
        
        
        print *,""
        print *,"done"
        print *,""
        
        
    contains
!---^^^^^^^^





        pure real(kind=real64) function r_i( b0,vF,e,s_0,s_1,s_2,T )    
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the athermal atom scattering rate
            real(kind=real64),intent(in)            ::      b0              !   nearest neighbour separation
            real(kind=real64),intent(in)            ::      vF              !   fermi velocity
            real(kind=real64),intent(in)            ::      S_0,S_1,S_2     !   scattering rate constants
            real(kind=real64),intent(in)            ::      T               !   temperature
            real(kind=real64),intent(in)            ::      e               !   excess potential energy
        
        
            r_i = b0*( S_0*abs(e) + T*(S_1 + T*S_2) ) + vF
            r_i = vF*( S_0*abs(e) + T*(S_1 + T*S_2) ) / r_i            
            
            return
        end function r_i
            
        
        pure real(kind=real64) function r_theta( b0,vF,s_1,s_2,T )    
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the thermal atom scattering rate
            real(kind=real64),intent(in)            ::      b0              !   nearest neighbour separation
            real(kind=real64),intent(in)            ::      vF              !   fermi velocity
            real(kind=real64),intent(in)            ::      S_1,S_2         !   scattering rate constants
            real(kind=real64),intent(in)            ::      T               !   temperature
        
        
            r_theta = b0*( T*(S_1 + T*S_2) ) + vF
            r_theta = vF*( T*(S_1 + T*S_2) ) / r_theta    
            
            return
        end function r_theta
        

            
        subroutine computeDefectsFromHist( n,nAtoms,E,dE,T, sigma, b0,vF,S_0,S_1,S_2 , kbar,r )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given that there are nAtoms in the system, 
    !*      and that n have energy between E and E+dE
    !*      find the expected thermal count from the broadened M-B distribution
    !*      and the thermal and defect scattering rates.
    !*      compute the expected scattering rate and the expected defect count
            integer,intent(in)                  ::      n,nAtoms                        !   atom counts
            real(kind=real64),intent(in)        ::      E,dE                            !   energy window          
            real(kind=real64),intent(in)        ::      T                               !   temperature
            real(kind=real64),intent(in)        ::      sigma, b0,vF,S_0,S_1,S_2        !   constants 
            real(kind=real64),intent(out)       ::      kbar                            !   expected defect count
            real(kind=real64),intent(out)       ::      r                               !   expected rate
            
            real(kind=real64)           ::      nBar  
            real(kind=real64)           ::      rtheta,ri
             

            
            real(kind=real128)              ::      sqrt2s  
            
            if (T>0.0d0) then
                beta = 1/(kB*T)
            else
                pp = 0
                if (e==0) pp = 1.0d0
            end if
             
            dd = e/sigma
            ee = e - 2*beta*sigma*sigma
            sqrt2s = sqrt( 2.0_real128 )*sigma
            
            pp = exp( -dd*dd/2 )*sqrt2onPI*sigma*ee                                                       &
               + exp( -beta*(ee+e) )*(sigma*sigma + ee*ee)*(1 + erf( ee/sqrt2s ))  
            
            
            pp = pp * 2*beta*beta*beta                       
            
            BroadenedMaxwellBoltzmann = max(pp,0.0d0)        !   possible for rounding error to give < 0 for very low T    
                                
            return
        end function BroadenedMaxwellBoltzmann        
        
        
 
            
    end program thermalConductivity
        
        
        
            
        
        
                                                                             

