c     RELATIVISTIC UNITS (c=\hbar=m=1) are used everywhere in the program
c     except the files, which contain potentials. Since these files
c     initially come from Ilya, they contain potentials in
c     atomic units, mutiplied by R. (i.e. pure Coulomb potential
c     must be constant there)

      include 'inc.par'
      common /momentum_projection/ amu,amj_max
      common /r_nuc/ r01,r02
      common /TargProj_prop/ Proj_mass,Targ_mass,Proj_vel,
     &b_ImpactParam
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /ferm/ cf1,af1,cf2,af2
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_mod/ nuc_model
      common /Barycentres/ RadiusOne,RadiusTwo
      common /common_dkb/ dkb
      common /staterangeinput/ state_range_input
      common /step_progress/ xi_stepslower,xi_stepsupper,ii_xi
      common /recoil/ Recoil_on_off
      common /neg_cont/ Manual_ncont_states
      real*8, dimension(:),allocatable:: eigval,eigval_e,eigval_o,
     &aaeigval,d_number_states_mj,d_number_states_mj_even,
     &d_number_states_mj_odd
      real*8, dimension(:,:),allocatable:: e,dmat,
     &eigval_mj,eigval_e_mj,eigval_o_mj,projMatMultipole
      real*8, dimension(:,:,:),allocatable:: wave,wave_new_at_dip,
     &wave_new,wave_new_prev,wave_new_even,wave_new_even_prev,
     &wave_new_at_dip_even,wave_new_odd,wave_new_odd_prev,
     &wave_new_at_dip_odd,alt_dmat
      real*8, dimension(:,:,:,:),allocatable:: vmat,wave_new_mj,
     &wave_new_even_mj,wave_new_odd_mj,wave_even_odd_combined
      real*8, dimension(:,:,:,:,:),allocatable::dvdRmatdkb1
      real*8, dimension(:,:,:,:,:,:),allocatable::dvdRmatdkb2
      integer, dimension(:,:),allocatable:: number_states,
     &number_states_perm
      complex*16, dimension(:,:),allocatable:: bb,dd,pp,
     &aaeigvecr,mm,mmeven,mmodd,bb_mjj,bb_mjj_even,bb_mjj_odd
      complex*16, dimension(:), allocatable :: ddmatnorm,coeff,
     &coefffornorm!,coefffornorm_prev
      real*8 i_xi,i_xi_next,i_xi_prev, energy_lowest_bound
      integer en_eval1,en_eval1j,xi_stepslower,xi_stepsupper,
     &c_ivalues(8)
      complex*16 summe
      character*2 How_many_kappas
      character*3 en_eval1_fileoutput,en_eval1j_fileoutput
      character*4 Charge_ofnuc_1,Charge_ofnuc_2,How_fast
      character*8 c_date
      character(len=:), allocatable:: status_q
      character*10 c_time,c_zone
      character*13 inp_dir
      common /input_directory/ inp_dir
      logical Manual_Coeff_Input,dkb,Sudden_approx,
     &state_range_input,Plots_wanted,Recoil_on_off,
     &init_recoil,Manual_ncont_states,mkdirs,dipping_wave_allocated,
     &unfreeze_basis
      inp_dir = 'input_output/'
      call factor()

      dipping_wave_allocated = .false.
      unfreeze_basis = .false.
      energy_lowest_bound = 1.d0
      mkdirs=makedirqq('Coeffs')
      mkdirs=makedirqq('Mat_norms')

      open(1,file=inp_dir//'inp.inp',status='old')
      read(1,*) z_nuc1,r_sq1
      read(1,*) z_nuc2,r_sq2
      read(1,*) amu,amj_max
      read(1,*) nuc_model
      read(1,*) xi_stepslower,xi_stepsupper
      read(1,*) xi_range
      read(1,*) Proj_mass
      read(1,*) Targ_mass
      read(1,*) Proj_vel
      read(1,*) b_ImpactParam
      read(1,*) state_range_input
      read(1,*) Manual_ncont_states
      read(1,*) dkb
      read(1,*) Manual_Coeff_Input
      read(1,*) Sudden_approx
      read(1,*) Plots_wanted
      read(1,*) Recoil_on_off
      close(1)

      init_recoil=Recoil_on_off

      if(z_nuc2.gt. z_nuc1) then
         open(1,file=inp_dir//'inp.inp',status='old')
         read(1,*) z_nuc2,r_sq2
         read(1,*) z_nuc1,r_sq1
         close(1)
      endif
      if(Proj_mass.gt. Targ_mass) then
         open(1,file=inp_dir//'inp.inp',status='old')
         read(1,*)
         read(1,*)
         read(1,*)
         read(1,*)
         read(1,*)
         read(1,*)
         read(1,*) Targ_mass
         read(1,*) Proj_mass
         close(1)
      endif

      nvmat=2
      if(dkb)nvmat=9

      if(nuc_model.eq.2)then
         r_sq1=r_sq1*dsqrt(5.d0/3.d0)
         r_sq2=r_sq2*dsqrt(5.d0/3.d0)
      endif

cc      write(*,*) 'Z1=',z_nuc1,'Z2=',z_nuc2
cc      write(*,*) 'M=',amu
c      select case(nuc_model)
c      case(0)
cc         write(*,*) 'point-like nucleus'
c      case(1)
cc         write(*,*) 'homogeneously charged shell'
c      case(2)
cc         write(*,*) 'homogeneously charged sphere'
c      case(3)
cc         write(*,*) 'Fermi nuclear charge distribution'
c      end select
      open(1,file=inp_dir//'data.inp',status='old')
      read(1,*) rmin1,rmax1 ! size of the box
      read(1,*) up_energy
      read(1,*) wavenumb,amplitude
      close(1)

      aNat_unit_Elec_field=1.32329d18/dsqrt(c)
      permittivity= 8.8541878176d-12
      speed_of_light=299792458d0
      Oscillation_time=2.d0*pi/wavenumb
      E_amp=dsqrt(2.d0*amplitude/permittivity/speed_of_light)
      if(wavenumb.eq.0.d0)then
         amplitude=1.d0
      else
         amplitude=E_amp/aNat_unit_Elec_field/wavenumb
      endif

      open(1,file=inp_dir//'nm.inp',status='old')
      read(1,*) nm
      close(1)
cc      write(*,*) 'Number of splines',nm
      open(1,file=inp_dir//'nkap.inp',status='old')
      read(1,*) nkap
      close(1)

cc      write(*,*) 'different kappas',nkap

c      allocate(ang_sigx_plus(-nkap:nkap,0:2*nkap,-nkap:nkap))
c      allocate(ang_sigx_minus(-nkap:nkap,0:2*nkap,-nkap:nkap))
c      call ang_sigxx(ang_sigx,nkap,0.5d0)
c      do k1=-nkap,nkap
c         do k2=-nkap,nkap
c         if(k1*k2.ne.0)then
c            do L=0,2*nkap
c            if(ang_sigx(k1,L,k2).ne.0.d0)then
c            write(*,*)'METHOD 1',ang_sigx(k1,L,k2)
c            write(*,*)k1,L,k2
c            endif
c            if(Coeff_maker2(0.5d0,1.d0*k1,1.d0*k2,l).ne.0.d0)then
c            write(*,*)'METHOD 2',Coeff_maker2(0.5d0,1.d0*k1,1.d0*k2,l)
c            write(*,*)k1,L,k2
c            endif
c            enddo
c            pause
c         endif
c         enddo
c      enddo

c      do k1=-4,4
c         do k2=-4,4
c            do l=0,3
c            if(Coeff_maker(amu,k1*1.d0,k2*1.d0,L).ne.0.d0)then
c            write(*,*)k1,k2,l,Coeff_maker(amu,k1*1.d0,k2*1.d0,L)
c            endif
c            enddo
c         enddo
c      enddo
c      stop

      call D01BAZ(-1.d0, 1.d0, 4, W4, t4, IFAIL)
      call D01BAZ(-1.d0, 1.d0, 8, W8, t8, IFAIL)
      call D01BAZ(-1.d0, 1.d0, 16, W16, t16, IFAIL)
      call D01BAZ(-1.d0, 1.d0, 32, W32, t32, IFAIL)
      call D01BAZ(-1.d0, 1.d0, 64, W64, t64, IFAIL)
      call D01BAZ(-1.d0, 1.d0, 6, W6, t6, IFAIL)

c      if(Nuc_model.eq.3)then
      af1=0.5233875553d0
      dsh=7.d0/3.d0*pi**4*af1**4-5.d0/3.d0*r_sq1**2*pi**2*af1**2
      bsh=5.d-1*(10.d0/3.d0*pi**2*af1**2-5.d0/3.d0*r_sq1**2)

      if(dsqrt(bsh**2-dsh)-bsh.lt.0.d0)then
         cf1=1.d0
         write(*,*) 'ERROR 1 IN CF'
         pause
      else
         cf1=dsqrt(dsqrt(bsh**2-dsh)-bsh)
      endif
      af2=0.5233875553d0
      dsh=7.d0/3.d0*pi**4*af2**4-5.d0/3.d0*r_sq2**2*pi**2*af2**2
      bsh=5.d-1*(10.d0/3.d0*pi**2*af2**2-5.d0/3.d0*r_sq2**2)

      if(dsqrt(bsh**2-dsh)-bsh.lt.0.d0) then
         cf2=1.d0
         write(*,*) 'ERROR 2 IN CF'
         pause
      else
         cf2=dsqrt(dsqrt(bsh**2-dsh)-bsh)
      endif

      r01=r_sq1*0.0025896063d0
      r02=r_sq2*0.0025896063d0

c     call test_potentials(nkap)
c     For the case of homogeneously charged sphere we have to
c     multiply root-mean-square radius by sqrt(5.d0/3.d0)

      if(Nuc_model.gt.1) then
         r0=r0*dsqrt(5.d0/3.d0)
      endif

      nu=nm+ns+2
      az1=z_nuc1/c
      az2=z_nuc2/c
      i_xi_start=(xi_stepslower+1)*xi_range*1.d0/xi_stepsupper
      call date_and_time(c_date,c_time,c_zone, c_ivalues)
      rmin=rmin1
      rmax=rmax1
      call NewInterNuclDist(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,Targ_mass,
     &Proj_vel,b_ImpactParam,0.d0,dist_min)
      gamma_inf=1.d0/(1.d0-Proj_vel**2.d0)**0.5d0

      open(333,file='Mat_norms/D_Mat_Norm.dat')
      allocate(number_states_perm(2,-nkap:nkap))

      do ii_xi=xi_stepslower,xi_stepsupper,2

c      if(ii_xi.ge.xi_stepslower+2)rmin=(rmin1/RadiusOne)*
c     & (2.d0*Starting_distance*Proj_mass/(Proj_mass+Targ_mass))

c    i_xi=(ii_xi+1)*xi_range*1.d0/
c     &max(xi_stepsupper,abs(xi_stepslower))
c    i_xi_prev=i_xi - 2.d0*xi_range/
c     &max(xi_stepsupper,abs(xi_stepslower))
c    i_xi_next=i_xi + 2.d0*xi_range/
c     &max(xi_stepsupper,abs(xi_stepslower))

        i_xi=derf((ii_xi+1)*dsqrt(pi)/
     &  max(xi_stepsupper,abs(xi_stepslower)))*xi_range
        i_xi_prev=derf((ii_xi-1)*dsqrt(pi)/
     &  max(xi_stepsupper,abs(xi_stepslower)))*xi_range
        i_xi_next=derf((ii_xi+3)*dsqrt(pi)/
     &  max(xi_stepsupper,abs(xi_stepslower)))*xi_range

c      xi_startP=(ii_xi)*xi_range*1.d0/
c     &max(xi_stepsupper,abs(xi_stepslower))
c      xi_endP=(ii_xi+2)*xi_range*1.d0/
c     &max(xi_stepsupper,abs(xi_stepslower))

        step_faktor=0.5d0*(i_xi_next-i_xi_prev)

        xi_startP=derf((ii_xi)*dsqrt(pi)/
     &  max(xi_stepsupper,abs(xi_stepslower)))*xi_range
        xi_endP=derf((ii_xi+2)*dsqrt(pi)/
     &  max(xi_stepsupper,abs(xi_stepslower)))*xi_range

c call Axis_angle_Theta(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,Targ_mass,
c     &Proj_vel,0.11394267720000001d0,0.d0,axis_angle_theta1)
c call DofThetawrtXi(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,Targ_mass,
c     &Proj_vel,0.11394267720000001d0,0.d0,dThetadXi)

        if(.not.sudden_approx) then
          call NewInterNuclDist(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*i_xi,distance)
          call NewInterNuclDist(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*i_xi_prev,
     &    distance_prev)
          call NewInterNuclDist(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*i_xi_next,
     &    distance_next)
          call Axis_angle_Theta(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*i_xi,
     &    axis_angle_theta1)
          call NewInterNuclDist(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*xi_startP,
     &    distance_startP)
          call NewInterNuclDist(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*xi_endP,distance_endP)
        else
          distance=dabs(proj_vel*i_xi)
          distance_prev=dabs(proj_vel*i_xi_prev)
          distance_next=dabs(proj_vel*i_xi_next)
        endif

        delta_distance=distance_endP-distance_startP

        if(sudden_approx) then
          distance=distance/2.d0+rmin*1.d-5
          distance_prev=distance_prev/2.d0+rmin*1.d-5
          distance_next=distance_next/2.d0+rmin*1.d-5
        else
          distance=distance/2.d0
          distance_startP=distance_startP/2.d0
        endif

        if(ii_xi.eq. xi_stepslower)then
          Starting_distance=distance
        endif

        R_std_dev=dabs(1.d0/(max(az1,az2)))
        Translation_factor=1.d0
        !dexp(-(2.d0*distance)**2/2.d0/R_std_dev**2)
        !0.5d0*derfc(R_std_dev*
c     &(dabs(2.d0*distance)-(dabs(2.d0*Starting_distance)-
c     &min(1.d0/R_std_dev,Starting_distance))))
         write(*,*)'TRANSLATION FACTOR',Translation_factor

        RadiusOne=2.d0*distance*Proj_mass/(Proj_mass+Targ_mass)
        RadiusTwo=2.d0*distance*Targ_mass/(Proj_mass+Targ_mass)

        write(*,*)'RMAX', rmax
        write(*,*)'RMIN' ,rmin

        write(*,*)'RADIUS ONE', RadiusOne
        write(*,*)'RADIUS TWO', RadiusTwo
        write(*,*)'INTERNUCLEAR DISTANCE',distance*2.d0/0.0025896d0

        if(.not.sudden_approx) then
          call DofRwrtXi(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,Targ_mass,
     &    Proj_vel,b_ImpactParam,1.d0*i_xi,dRdXi)
          call DofTwrtXi(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,Targ_mass,
     &    Proj_vel,b_ImpactParam,1.d0*i_xi,dTdXi)
          call DofThetawrtXi(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*i_xi,dthetadt)
        else
          dRdXi=proj_vel
          dTdXi=1.d0
        endif

        if(.not.sudden_approx)then
          call New_Elapsed_time(1.d0*z_nuc1,1.d0*z_nuc2,Proj_mass,
     &    Targ_mass,Proj_vel,b_ImpactParam,1.d0*xi_endP,elapsed_time)
        else
          elapsed_time=i_xi
        endif

        write(*,*) 'ELAPSED TIME', elapsed_time
        write(*,*) 'DRDXI',dRdXi
        write(*,*) 'DTDXI',dTdXi
        write(*,*) 'VEL',dRdXi/dTdXi
        write(*,*) 'XI',i_xi
        write(*,*) 'DELTA XI', step_faktor

c      The constants for the Fermi charge distribution. Does not
c      work for small Z

c      if(Nuc_model.eq.3)then
c      af1=0.5233875553d0
c      dsh=7.d0/3.d0*pi**4*af1**4-5.d0/3.d0*r_sq1**2*pi**2*af1**2
c      bsh=5.d-1*(10.d0/3.d0*pi**2*af1**2-5.d0/3.d0*r_sq1**2)

c      if(dsqrt(bsh**2-dsh)-bsh.lt.0.d0)then
c         cf1=1.d0
c         write(*,*) 'ERROR 1 IN CF'
c         pause
c      else
c         cf1=dsqrt(dsqrt(bsh**2-dsh)-bsh)
c      endif
c      af2=0.5233875553d0
c      dsh=7.d0/3.d0*pi**4*af2**4-5.d0/3.d0*r_sq2**2*pi**2*af2**2
c      bsh=5.d-1*(10.d0/3.d0*pi**2*af2**2-5.d0/3.d0*r_sq2**2)

c      if(dsqrt(bsh**2-dsh)-bsh.lt.0.d0)then
c         cf2=1.d0
c         write(*,*) 'ERROR 2 IN CF'
c         pause
c      else
c         cf2=dsqrt(dsqrt(bsh**2-dsh)-bsh)
c      endif
c      endif

c      r01=r_sq1*0.0025896063d0
c      r02=r_sq2*0.0025896063d0

c      call test_potentials(nkap)
c      For the case of homogeneously charged sphere we have to
c      multiply root-mean-square radius by sqrt(5.d0/3.d0)

c      if(Nuc_model.gt.1) r0=r0*dsqrt(5.d0/3.d0)

c      NOTE The allocation of vmat has been expanded. vmat is not only
c     the v matrix from Johnson's paper, but also the nvmat'th vmat is
c      the dv/dr terms required for the CC matrix.

        allocate(e(2*nm,-nkap:nkap))
        allocate(vmat(nm,nm,0:2*nkap,nvmat))
        allocate(dmat(2*nm,2*nm))
        allocate(dvdRmatdkb1(nm,nm,-nkap:nkap,0:2*nkap,2))
        allocate(dvdRmatdkb2(nm,nm,-nkap:nkap,-nkap:nkap,0:2*nkap,2))
        allocate(number_states(2,-nkap:nkap))
        allocate(wave(2*nm,2*nm,-nkap:nkap))
        allocate(alt_dmat(2*nm,2*nm,-nkap:nkap))

        write(*,*) 'SIZE OF MATRIX',4*(nm**2)*(2*nkap+2)*8/1000000,'MB'
        write(*,*) 'ENTERING B_SPLINE_CALCULATION'
c         Does anything in here depend on xi? Only up_energy changes.
        call b_spline_calculation_no_laser
     &    (nstates,nsto,nste,nm,nu,nkap,number_states,rmin,rmax,
     &    wave,vmat,nvmat,e,up_energy,dmat,
     &    dvdRmatdkb1,dvdRmatdkb2,alt_dmat)
        if(ii_xi.eq.xi_stepslower)then
          number_states_perm=number_states
          nstates_perm=nstates
          nste_perm=nste
          nsto_perm=nsto
        endif
        number_states=number_states_perm
        nstates=nstates_perm
        nste=nste_perm
        nsto=nsto_perm

        call sign_of_states_monopole(number_states,
     &  wave,nstates,nm,nkap,ii_xi,xi_stepslower,rmin,rmax,'a')
        if(z_nuc1.ne.z_nuc2)then
          allocate(eigval(nstates))
          allocate(wave_new(nstates,2*nm,-nkap:nkap))
          i_even_odd_normal = 0
          call buildMultipoleBasis(nm,nkap,nstates,number_states,wave,
     &    vmat,nvmat,e, ii_xi, xi_stepslower,rmin,rmax,eigval,
     &    wave_new,i_even_odd_normal)
        else
          allocate(eigval_o(nsto))
          allocate(eigval_e(nste))
          allocate(wave_new_even(nste,2*nm,-nkap:nkap))
          allocate(wave_new_odd(nsto,2*nm,-nkap:nkap))
          i_even_odd_normal = 1
          call buildMultipoleBasis(nm,nkap,nste,number_states,wave,
     &    vmat,nvmat,e, ii_xi, xi_stepslower,rmin,rmax,eigval_e,
     &    wave_new_even,i_even_odd_normal)
          i_even_odd_normal = 2
          call buildMultipoleBasis(nm,nkap,nsto,number_states,wave,
     &    vmat,nvmat,e, ii_xi, xi_stepslower,rmin,rmax,eigval_o,
     &    wave_new_odd,i_even_odd_normal)
        endif
        d_amuOrig=amu
        d_amjmaxOrig=amj_max
        call DiracAngularJ(dble(nkap),d_mjMax)
        n_jstates=1
        if(b_ImpactParam .gt. 0.d0)then
          n_jstates=nint(d_mjMax+0.5d0)
          amu=d_mjMax
          amj_max=amu
        endif
        if(z_nuc1.ne.z_nuc2)then
          allocate(eigval_mj(nstates,-n_jstates:n_jstates))
          eigval_mj=0.d0
          allocate(wave_new_mj(nstates,2*nm,-nkap:nkap,
     &    -n_jstates:n_jstates))
          wave_new_mj=0.d0
          if(b_ImpactParam .gt. 0.d0)then
            do n_jstate=n_jstates,1,-1
              amu=amu-1.d0
              amj_max=amj_max-1.d0
              eigval_mj(:,n_jstate)=eigval
              eigval_mj(:,-n_jstate)=eigval
              wave_new_mj(:,:,:,n_jstate)=wave_new
              wave_new_mj(:,:,:,-n_jstate)=wave_new
            enddo
            nsts=2*n_jstates*nstates
            allocate(bb_mjj(nsts,nsts))
            allocate(d_number_states_mj(2*n_jstates*nstates))
            call Rotating_INaxis(nstates,nm,nkap,dthetadt,
     &      nsts,bb_mjj,dmat,alt_dmat,wave_new,d_mjMax,
     &      d_number_states_mj,n_jstates)
          endif
        else
          allocate(eigval_e_mj(nste,-n_jstates:n_jstates))
          eigval_e_mj=0.d0
          allocate(wave_new_even_mj(nste,2*nm,-nkap:nkap,
     &    -n_jstates:n_jstates))
          wave_new_even_mj=0.d0
          allocate(eigval_o_mj(nsto,-n_jstates:n_jstates))
          eigval_o_mj=0.d0
          allocate(wave_new_odd_mj(nste,2*nm,-nkap:nkap,
     &    -n_jstates:n_jstates))
          wave_new_odd_mj=0.d0
          if(b_ImpactParam .gt. 0.d0)then
            do n_jstate=n_jstates,1,-1
              amu=amu-1.d0
              amj_max=amj_max-1.d0
              eigval_e_mj(:,n_jstate)=eigval_e
              eigval_e_mj(:,-n_jstate)=eigval_e
              eigval_o_mj(:,n_jstate)=eigval_o
              eigval_o_mj(:,-n_jstate)=eigval_o
              wave_new_even_mj(:,:,:,n_jstate)=wave_new_even
              wave_new_even_mj(:,:,:,-n_jstate)=wave_new_even
              wave_new_odd_mj(:,:,:,n_jstate)=wave_new_odd
              wave_new_odd_mj(:,:,:,-n_jstate)=wave_new_odd
            enddo
            nsts=2*n_jstates*nste
            allocate(bb_mjj_even(nsts,nsts))
            allocate(d_number_states_mj_even(2*n_jstates*nste))
            call Rotating_INaxis_even(nste,nm,nkap,dthetadt,
     &      nsts,bb_mjj_even,dmat,alt_dmat,wave_new_even,d_mjMax,
     &      d_number_states_mj_even,n_jstates)
            nsts=2*n_jstates*nsto
            allocate(bb_mjj_odd(nsts,nsts))
            allocate(d_number_states_mj_odd(2*n_jstates*nsto))
            call Rotating_INaxis_odd(nsto,nm,nkap,dthetadt,
     &      nsts,bb_mjj_odd,dmat,alt_dmat,wave_new_odd,d_mjMax,
     &      d_number_states_mj_odd,n_jstates)
          endif
        endif

        if(b_ImpactParam .gt. 0.d0)then
          amu=d_amuOrig
          amj_max=d_amjmaxOrig
        endif
        deallocate(number_states)
        deallocate(wave)
        deallocate(e)

        write(*,*) 'ENTER MAIN MATRIX'
        if(z_nuc1.ne.z_nuc2)then
          if (b_ImpactParam .gt. 0.d0)then
            ! Expand the size of eigval and wavenew
            deallocate(eigval)
            deallocate(wave_new)
            allocate(eigval(2*n_jstates*nstates))
            allocate(wave_new(2*n_jstates*nstates,2*nm,-nkap:nkap))
C           This function redefines nstates=2*n_jstates*nstates
            call redefineEigvalWaveNew(n_jstates,eigval,eigval_mj,
     &      wave_new,wave_new_mj,nstates,nm,nkap)
            deallocate(eigval_mj)
            deallocate(wave_new_mj)
          else
            allocate(d_number_states_mj(nstates))
            d_number_states_mj = amu
          endif
          allocate(mm(nstates,nstates))
          CALL MatrixMMGenerator_neqZ(dTdXi,dRdXi,eigval,nkap,vmat,
     &    wave_new,nstates,nm,nvmat,mm,dvdRmatdkb1,dvdRmatdkb2,
     &    d_number_states_mj)
          deallocate(d_number_states_mj)
          if(b_ImpactParam .gt. 0.d0)then
            mm = mm + bb_mjj
            deallocate(bb_mjj)
          endif
        else
          if (b_ImpactParam .gt. 0.d0)then
            deallocate(eigval_e)
            deallocate(wave_new_even)
            allocate(eigval_e(2*n_jstates*nste))
            allocate(wave_new_even(2*n_jstates*nste,2*nm,-nkap:nkap))
C           This function redefines nste=2*n_jstates*nste
            call redefineEigvalWaveNew(n_jstates,eigval_e,eigval_e_mj,
     &      wave_new_even,wave_new_even_mj,nste,nm,nkap)
            deallocate(eigval_e_mj)
            deallocate(wave_new_even_mj)
            deallocate(eigval_o)
            deallocate(wave_new_odd)
            allocate(eigval_o(2*n_jstates*nsto))
            allocate(wave_new_odd(2*n_jstates*nsto,2*nm,-nkap:nkap))
C           This function redefines nsto=2*n_jstates*nsto
            call redefineEigvalWaveNew(n_jstates,eigval_o,eigval_o_mj,
     &      wave_new_odd,wave_new_odd_mj,nsto,nm,nkap)
            deallocate(eigval_o_mj)
            deallocate(wave_new_odd_mj)
          else
            allocate(d_number_states_mj_even(nste))
            allocate(d_number_states_mj_odd(nsto))
            d_number_states_mj_even = amu
            d_number_states_mj_odd = amu
          endif
          allocate(mmeven(nste,nste))
          CALL MatrixMMGenerator_eqZeven(dTdXi,dRdXi,eigval_e,nkap,
     &    vmat,wave_new_even,nste,nm,nvmat,mmeven,dvdRmatdkb1,
     &    dvdRmatdkb2,d_number_states_mj_even)
          deallocate(d_number_states_mj_even)
          if(b_ImpactParam .gt. 0.d0)then
            mmeven = mmeven + bb_mjj_even
            deallocate(bb_mjj_even)
          endif
          allocate(mmodd(nsto,nsto))
          CALL MatrixMMGenerator_eqZodd(dTdXi,dRdXi,eigval_o,nkap,
     &    vmat,wave_new_odd,nsto,nm,nvmat,mmodd,dvdRmatdkb1,
     &    dvdRmatdkb2,d_number_states_mj_odd)
          deallocate(d_number_states_mj_odd)
          if(b_ImpactParam .gt. 0.d0)then
            mmodd = mmodd + bb_mjj_odd
            deallocate(bb_mjj_odd)
          endif
        endif
        write(*,*) 'EXIT MAIN MATRIX'

        deallocate(dvdRmatdkb1)
        deallocate(dvdRmatdkb2)
        deallocate(vmat)
        if(z_nuc1.eq.z_nuc2)then
          nstates=nste+nsto
        endif

        deallocate(dmat)

        if(ii_xi.eq.xi_stepslower) then
          allocate(coeff(nstates))
   !         allocate(coeff_prev(nstates))
          allocate(coefffornorm(nstates))
   !         allocate(coefffornorm_prev(nstates))
        endif

        write(*,*) 'CREATING B'
        allocate(bb(nstates,nstates))
        bb=0.d0
        if(z_nuc1.eq.z_nuc2)then
          bb(1:nste,1:nste)=mmeven
          deallocate(mmeven)
          bb(1+nste:,1+nste:)=mmodd
          deallocate(mmodd)
        else
          bb=mm
          deallocate(mm)
        endif

!      do i=1,nstates
!         do j=i+1,nstates
!         if(bb(i,j).ne.-bb(j,i))then
!         write(*,*)i,j
!         write(*,*)bb(i,j),-bb(j,i)
!         endif
!         enddo
!      write(*,*)bb(i,i),i
!      enddo
!      pause

        write(*,*) 'ALLOCATING TIME DEPENDENT MATRICES'
        allocate(aaeigval(nstates))
        allocate(aaeigvecr(nstates,nstates))
        aaeigval=0.d0
        aaeigvecr=0.d0
        allocate(pp(nstates,nstates))
        allocate(dd(nstates,nstates))
        allocate(ddmatnorm(nstates))
        write(*,*) 'DIAGONALISING B',2*8*(nstates)**2/
     &  10**6,'MB'
        call c_diagonal(nstates,bb,aaeigval,aaeigvecr)
        deallocate(bb)
        pp=0.d0
        do i=1,nstates
          do j=1,nstates
            pp(i,j)=cdexp(aaeigval(i)*(0,-1)*step_faktor)*
     &              dconjg(aaeigvecr(j,i))
          enddo
        enddo
        deallocate(aaeigval)
        dd=0.d0
        do i=1,nstates
          do j=1,nstates
            do k=1,nstates
              dd(i,j)=dd(i,j)+aaeigvecr(i,k)*pp(k,j)
            enddo
          enddo
        enddo

        deallocate(pp)
        deallocate(aaeigvecr)

        do i=1,nstates
          summe=0.d0
          do j=1,nstates
            summe=summe+dd(i,j)
          enddo
          ddmatnorm(i)=summe
        enddo
        ddmatnorm_1=maxval(cdabs(ddmatnorm))

        write(333,*)i_xi,'MATRIX NORM',ddmatnorm_1
        deallocate(ddmatnorm)

        if(ii_xi.eq.xi_stepslower)then
          coeff=0.d0
!      coeff_prev=0.d0
!      coefffornornm_prev=0.d0
          write(Charge_ofnuc_1,173)int(z_nuc1)
  173     format(I4.4)
          write(Charge_ofnuc_2,178)int(z_nuc2)
  178     format(I4.4)
          write(How_many_kappas,194)nkap
  194     format(I2.2)
          write(How_fast,195)int(Proj_vel*1000.d0)
  195     format(I4.4)
          if(z_nuc1.ne.z_nuc2)then
            i=1
            do while(eigval(i).le.-1.d0)
              i=i+1
            enddo
            lowest_bound=i
!         coefffornorm(lowest_bound)=1.d0
            coeff(lowest_bound)=1.d0
          else
            i=1
            do while(eigval_e(i).le.-1.d0)
              i=i+1
            enddo
            lowest_bound_e=i
!         coefffornorm(lowest_bound_e)=1.d0/dsqrt(2.d0)
!         coefffornorm_prev(lowest_bound_e)=1.d0/dsqrt(2.d0)
            coeff(lowest_bound_e)=1.d0/dsqrt(2.d0)
            i=1
            do while(eigval_o(i).le.-1.d0)
              i=i+1
            enddo
            lowest_bound_o=nste+i
!       coefffornorm(lowest_bound_o)=1.d0/dsqrt(2.d0)
!         coefffornorm_prev(lowest_bound_o)=1.d0/dsqrt(2.d0)
            coeff(lowest_bound_o)=1.d0/dsqrt(2.d0)
          endif
!      coeff=coefffornorm
!      coeff_prev=coeff
        endif

        if(z_nuc1.ne.z_nuc2)then
          energy_lowest_bound = eigval(lowest_bound)
          if ((energy_lowest_bound .lt. -1.d0) .and.
     &        .not.(dipping_wave_allocated)) then
            allocate(wave_new_at_dip(nstates,2*nm,-nkap:nkap))
            dipping_wave_allocated = .true.
            wave_new_at_dip = wave_new_prev
          endif
          if ((energy_lowest_bound .gt. -1.d0) .and.
     &        unfreeze_basis) then
            deallocate(wave_new_at_dip)
            dipping_wave_allocated = .false.
          endif
          if (ii_xi.gt.xi_stepslower)then
            deallocate(wave_new_prev)
          endif
          allocate(wave_new_prev(nstates,2*nm,-nkap:nkap))
          wave_new_prev = wave_new
        else
          energy_lowest_bound = eigval_e(lowest_bound_e)
          if (eigval_o(lowest_bound_o - nste) < energy_lowest_bound)then
            energy_lowest_bound = eigval_o(lowest_bound_o - nste)
          endif
          if ((energy_lowest_bound .lt. -1.d0) .and.
     &        .not.(dipping_wave_allocated)) then
            allocate(wave_new_at_dip_even(nste,2*nm,-nkap:nkap))
            dipping_wave_allocated = .true.
            wave_new_at_dip_even = wave_new_even_prev
            allocate(wave_new_at_dip_odd(nsto,2*nm,-nkap:nkap))
            wave_new_at_dip_odd = wave_new_odd_prev
          endif
          if ((energy_lowest_bound .gt. -1.d0) .and.
     &        unfreeze_basis) then
            deallocate(wave_new_at_dip_even)
            deallocate(wave_new_at_dip_odd)
            dipping_wave_allocated = .false.
          endif
          if (ii_xi.gt.xi_stepslower)then
            deallocate(wave_new_even_prev)
            deallocate(wave_new_odd_prev)
          endif
          allocate(wave_new_even_prev(nste,2*nm,-nkap:nkap))
          allocate(wave_new_odd_prev(nsto,2*nm,-nkap:nkap))
          wave_new_even_prev = wave_new_even
          wave_new_odd_prev = wave_new_odd
        endif

        if (energy_lowest_bound .gt. -1.d0 .and.
     &  dipping_wave_allocated)then
          unfreeze_basis = .true.
        endif
        if ((energy_lowest_bound .lt. -1.d0) .or. unfreeze_basis)then
          allocate(projMatMultipole(nstates,nstates))
          if(z_nuc1.eq.z_nuc2)then
            call projection_matrix_frozen_basis(nstates,nm,
     &      nkap,alt_dmat,wave_new,wave_new_at_dip,projMatMultipole)
          else
            allocate(wave_even_odd_combined(2,nstates,2*nm,-nkap:nkap))
            wave_even_odd_combined(1,1:nste,:,:) = wave_new_even
            wave_even_odd_combined(1,1+nste:,:,:) = wave_new_odd
            wave_even_odd_combined(2,1:nste,:,:) = wave_new_at_dip_even
            wave_even_odd_combined(2,1+nste:,:,:) = wave_new_at_dip_odd
            call projection_matrix_frozen_basis(nstates,nm,
     &      nkap,alt_dmat,wave_even_odd_combined(1,:,:,:),
     &      wave_even_odd_combined(2,:,:,:),projMatMultipole)
            deallocate(wave_even_odd_combined)
          endif
        endif
        deallocate(alt_dmat)
        if(z_nuc1.eq.z_nuc2)then
          deallocate(wave_new_even)
          deallocate(wave_new_odd)
        else
          deallocate(wave_new)
        endif

C       Project forward to the moving basis
        if ((energy_lowest_bound .lt. -1.d0)
     &  .or. unfreeze_basis)then
          coefffornorm=0.d0
          do i=1,nstates
            summe=0.d0
            do k=1,nstates
              summe=summe+projMatMultipole(i,k)*coeff(k)
            enddo
            coefffornorm(i)=summe
          enddo
          deallocate(projMatMultipole)
          coeff=coefffornorm
        endif
        if (unfreeze_basis .and. .not.dipping_wave_allocated)then
          unfreeze_basis = .false.
        endif
CC******THE FINAL STEP!!!!!!!!!!!!!!!!
        coefffornorm=0.d0
        do i=1,nstates
          summe=0.d0
          do k=1,nstates
            summe=summe+dd(i,k)*coeff(k)
          enddo
          coefffornorm(i)=summe
        enddo
        coeff=coefffornorm
C       Project back to the frozen basis
        if (energy_lowest_bound .lt. -1.d0)then
          coefffornorm=0.d0
          do i=1,nstates
            summe=0.d0
            do k=1,nstates
              summe=summe+projMatMultipole(i,k)*coeff(k)
            enddo
            coefffornorm(i)=summe
          enddo
          deallocate(projMatMultipole)
          coeff=coefffornorm
        endif
c*******END OF THE FINAL STEP!!!!!!!!!!!
        if(z_nuc1.eq.z_nuc2)then
          write(*,*)cdabs(coeff(lowest_bound_e))**2,
     &    cdabs(coeff(lowest_bound_o))**2
        else
          write(*,*)cdabs(coeff(lowest_bound))**2
        endif

        iii_xi=ii_xi
        if(ii_xi.eq. 0)iii_xi=1
        if(z_nuc1.eq.z_nuc2)then
          do i=1,min(nste,nsto)
            en_eval1=2*i-1
            write(en_eval1_fileoutput, 74)en_eval1
   74       format(I3.3)
            en_eval1j=2*i
            write(en_eval1j_fileoutput, 75)en_eval1j
   75       format(I3.3)
            if(ii_xi.eq.xi_stepslower)then
              allocate(character(len=7)::status_q)
              status_q='replace'
            else
              allocate(character(len=3)::status_q)
              status_q='old'
            endif
            open(50513+en_eval1,
     &      file='Coeffs/Coeff_Prob_Raw_'//en_eval1_fileoutput//
     &      '_charges_'//Charge_ofnuc_1//'_'//Charge_ofnuc_2//
     &      '_nkap_'//How_many_kappas//'_vel_'//How_fast//
     &      '_Date_'//c_date//'.dat',status=status_q,position=
     &      'append')
            open(50513+en_eval1j,
     &      file='Coeffs/Coeff_Prob_Raw_'//en_eval1j_fileoutput//
     &      '_charges_'//Charge_ofnuc_1//'_'//Charge_ofnuc_2//
     &      '_nkap_'//How_many_kappas//'_vel_'//How_fast//
     &      '_Date_'//c_date//'.dat',status=status_q,position=
     &      'append')
            open(2256+en_eval1,
     &      file='Coeffs/Coeff_Prob_En_'//en_eval1_fileoutput//
     &      '_charges_'//Charge_ofnuc_1//'_'//Charge_ofnuc_2//
     &      '_nkap_'//How_many_kappas//'_vel_'//How_fast//
     &      '_Date_'//c_date//'.dat',status=status_q,position=
     &      'append')
            open(2256+en_eval1j,
     &      file='Coeffs/Coeff_Prob_En_'//en_eval1j_fileoutput//
     &      '_charges_'//Charge_ofnuc_1//'_'//Charge_ofnuc_2//
     &      '_nkap_'//How_many_kappas//'_vel_'//How_fast//
     &      '_Date_'//c_date//'.dat',status=status_q,position=
     &      'append')

            if(ii_xi.gt.xi_stepslower)then
              write(50513+en_eval1,'(4f19.11)')(iii_xi/abs(iii_xi))*
     &        2.d0*distance_startP/0.0025896063d0,
     &        dble(coeff(i)),aimag(coeff(i)),eigval_e(i)
              write(50513+en_eval1j,'(4f19.11)')(iii_xi/abs(iii_xi))*
     &        2.d0*distance_startP/0.0025896063d0,
     &        dble(coeff(i+nste)),aimag(coeff(i+nste)),eigval_o(i)
              write(2256+en_eval1,'(3f19.11)')(iii_xi/abs(iii_xi))*
     &        2.d0*distance_startP/0.0025896063d0,
     &        cdabs(coeff(i))**2,eigval_e(i)
              write(2256+en_eval1j,'(3f19.11)') (iii_xi/abs(iii_xi))*
     &        2.d0*distance_startP/0.0025896063d0,
     &        cdabs(coeff(i+nste))**2,eigval_o(i)
            endif
            close(50513+en_eval1)
            close(50513+en_eval1j)
            close(2256+en_eval1)
            close(2256+en_eval1j)
            deallocate(status_q)
          enddo
        else
          do i=1,nstates
            en_eval1=i
            write(en_eval1_fileoutput, 76)en_eval1
   76       format(I3.3)
            if(ii_xi.eq.xi_stepslower)then
              allocate(character(len=7)::status_q)
              status_q='replace'
            else
              allocate(character(len=3)::status_q)
              status_q='old'
            endif
            open(50513+en_eval1,
     &      file='Coeffs/Coeff_Prob_Raw_'//en_eval1_fileoutput//
     &      '_charges_'//Charge_ofnuc_1//'_'//Charge_ofnuc_2//
     &      '_nkap_'//How_many_kappas//'_vel_'//How_fast//
     &      '_Date_'//c_date//'.dat',status=status_q,position=
     &      'append')
            open(2256+en_eval1,
     &      file='Coeffs/Coeff_Prob_En_'//en_eval1_fileoutput//
     &      '_charges_'//Charge_ofnuc_1//'_'//Charge_ofnuc_2//
     &      '_nkap_'//How_many_kappas//'_vel_'//How_fast//
     &      '_Date_'//c_date//'.dat',status=status_q,position=
     &      'append')

            if(ii_xi.gt.xi_stepslower)then
              write(50513+i,'(4f19.11)') (iii_xi/abs(iii_xi))*2.d0*
     &        distance_startP/0.0025896063d0,
     &        dble(coeff(i)),aimag(coeff(i)),eigval(i)
              write(2256+i,'(3f19.11)') (iii_xi/abs(iii_xi))*2.d0*
     &        distance_startP/0.0025896063d0,
     &        cdabs(coeff(i))**2,eigval(i)
            endif
            close(50513+en_eval1)
            close(2256+en_eval1)
            deallocate(status_q)
          enddo
        endif

        if(ii_xi.eq.xi_stepslower)then
          open(127,file='Coeffs/Ionisation_Probability_charges_'
     &    //Charge_ofnuc_1//'_'//Charge_ofnuc_2//
     &    '_nkap_'//How_many_kappas//'_vel_'//How_fast//
     &    '_Date_'//c_date//'.dat')
        endif

        Prob_ionisation=0.d0
        Prob_electron_creation=0.d0
        if(z_nuc1.ne.z_nuc2)then
          do i=1,nstates
            if(eigval(i).gt. 1.d0)then
              Prob_ionisation=Prob_ionisation+cdabs(coeff(i))**2
            endif
            if(eigval(i).lt.-1.d0)then
              Prob_electron_creation=Prob_electron_creation+
     &        cdabs(coeff(i))**2
            endif
          enddo
          deallocate(eigval)
        elseif(z_nuc1.eq.z_nuc2)then
          do i=1,nste
            if(eigval_e(i).gt. 1.d0)then
              Prob_ionisation=Prob_ionisation+cdabs(coeff(i))**2
            endif
            if(eigval_e(i).lt.-1.d0)then
              Prob_electron_creation=Prob_electron_creation+
     &        cdabs(coeff(i))**2
            endif
          enddo
          deallocate(eigval_e)
          do i=1,nsto
            if(eigval_o(i).gt. 1.d0)then
              Prob_ionisation=Prob_ionisation+
     &        cdabs(coeff(i+nste))**2
            endif
            if(eigval_o(i).lt.-1.d0)then
              Prob_electron_creation=Prob_electron_creation+
     &        cdabs(coeff(i+nste))**2
            endif
          enddo
          deallocate(eigval_o)
        endif

        if (ii_xi.ne. 0)then
          write(127,*)(ii_xi/abs(ii_xi))*
     &    distance_endP/0.0025896063d0,Prob_ionisation,
     &    Prob_electron_creation
        elseif(ii_xi.eq. 0)then
          write(127,*)
     &    distance_endP/0.0025896063d0,Prob_ionisation,
     &    Prob_electron_creation
        endif

        vectornorm_1=0.d0
        do i=1,nstates
          vectornorm_1=vectornorm_1+cdabs(coeff(i))**2
        enddo
        write(333,*)i_xi,'MULTIPOLE NORM',vectornorm_1
c     Is renormalisation recommended for pair production?
        coeff=coeff/dsqrt(vectornorm_1)

        if(ii_xi.eq.0)then
          iii_xi=1
        else
          iii_xi=ii_xi
        endif

        write(*,*)'I-N DISTANCE',distance_endP/0.0025896063
        write(*,*)'IONISATION PROBABILITY',Prob_ionisation
        write(*,*)'E CREATION PROBABILITY',Prob_electron_creation
!      pause
        if(ii_xi.eq.xi_stepslower)then
          up_energy=up_energy*1.25d0
        endif

        deallocate(dd)

      enddo
      stop
      end

      subroutine draw_sigma(nm,nkap,nu,tknot,wave)
      include 'inc.par'
      real*8 ro(ns),ro1(ns),wave(2*nm,-nkap:nkap),tknot(nu)
      common /momentum_projection/ amu,amj_max
      real*8, dimension(:,:,:),allocatable:: plot

      allocate(plot(5,16*(nu-2*ns+1),0:nugl))

      plot=0.d0

      do kk=-nkap,nkap
        if(kk.ne.0)then
          num=0
          if(kk.gt.0)then
            l1=kk
            l2=kk-1
          else
            l1=-kk-1
            l2=-kk
          endif
          al1=dble(l1)
          al2=dble(l2)
          aj=iabs(kk)-5.d-1

          do is=nu-ns,ns,-1
            do ip=8,1,-1
              xx=tknot(is)+(tknot(is+1)-tknot(is))*dble(ip)/8.d0
              call dsplines(xx,ro,ro1,is,tknot,nu)

              g=0.d0
              f=0.d0

              if(is.eq.ns)then
                do i=2,ns
                  g=g+wave(i-1,kk)*ro(i)
                  f=f+wave(nm+i-1,kk)*ro(i)
                enddo
              elseif(is.eq.nu-ns)then
                do i=1,ns-1
                  g=g+wave(i-1+is-ns,kk)*ro(i)
                  f=f+wave(nm+i-1+is-ns,kk)*ro(i)
                enddo
              else
                do i=1,ns
                  g=g+wave(i-1+is-ns,kk)*ro(i)
                  f=f+wave(nm+i-1+is-ns,kk)*ro(i)
                enddo
              endif
              g=g/xx
              f=f/xx
              num=num+1
              do j=0,nugl
                ang=pi-j*pi/nugl
                dc=dcos(ang)
                ds=dsin(ang)
                gfact1=clebsh(al1,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &          spharm(l1,idint(amu-5.d-1),dc)
                gfact2=clebsh(al1,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &          spharm(l1,idint(amu+5.d-1),dc)
                ffact1=clebsh(al2,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &          spharm(l2,idint(amu-5.d-1),dc)
                ffact2=clebsh(al2,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &          spharm(l2,idint(amu+5.d-1),dc)

                if((xx.ne.plot(1,num,j)).and.(kk.ne.-nkap))then
                  write(*,*) 'something goes wrong with num'
                  write(*,*) xx,plot(1,num,j)
                  write(*,*) num,j
                  pause
                endif

                plot(1,num,j)=xx
                plot(2,num,j)=plot(2,num,j)+g*gfact1
                plot(3,num,j)=plot(3,num,j)+g*gfact2
                plot(4,num,j)=plot(4,num,j)+f*ffact1
                plot(5,num,j)=plot(5,num,j)+f*ffact2
              enddo
            enddo
          enddo
        endif
      enddo

      do i=1,num
        write(1,1) -plot(1,i,0),(plot(k,i,0),k=2,5),
     &  plot(2,i,0)**2+plot(3,i,0)**2+plot(4,i,0)**2+
     &  plot(5,i,0)**2
      enddo
      do i=num,1,-1
        write(1,1) plot(1,i,nugl),(plot(k,i,nugl),k=2,5),
     &  plot(2,i,nugl)**2+plot(3,i,nugl)**2+plot(4,i,nugl)**2+
     &  plot(5,i,nugl)**2
      enddo
      do i=1,num
        r=plot(1,i,0)
        do j=0,nugl
          ang=pi-j*pi/nugl
          dc=dcos(ang)
          ds=dsin(ang)
          x=r*dc
          y=r*ds
          if((dabs(x).lt.3.d0).and.(dabs(y).lt.3.d0))then
            dens=plot(2,i,j)**2+plot(3,i,j)**2+plot(4,i,j)**2+
     &      plot(5,i,j)**2
            write(11,*) x,y,dens
            write(12,*) x,y,dlog(dens)
            if((j.ne.0).and.(j.ne.nugl))then
              ang=pi+j*pi/nugl
              dc=dcos(ang)
              ds=dsin(ang)
              x=r*dc
              y=r*ds
              write(11,*) x,y,dens
              write(12,*) x,y,dlog(dens)
            endif
          endif
        enddo
      enddo

      deallocate(plot)
 1    format(e14.6,2x,e14.6,2x,e14.6,2x,e14.6,2x,e14.6,2x,e14.6)
      return
      end

      subroutine draw_sigma_raw(nm,nkap,nu,tknot,wave)
      include 'inc.par'
      real*8 ro(ns),ro1(ns),wave(2*nm,-nkap:nkap),tknot(nu)
      common /momentum_projection/ amu,amj_max
      real*8, dimension(:,:,:),allocatable:: plot

      allocate(plot(5,16*(nu-2*ns+1),0:nugl))

      plot=0.d0

      do kk=-nkap,nkap
        if(kk.ne.0)then
          num=0
          if(kk.gt.0)then
            l1=kk
            l2=kk-1
          else
            l1=-kk-1
            l2=-kk
          endif
          al1=dble(l1)
          al2=dble(l2)
          aj=iabs(kk)-5.d-1
          do is=nu-ns,ns,-1
            do ip=8,1,-1
              xx=tknot(is)+(tknot(is+1)-tknot(is))*dble(ip)/8.d0
              call dsplines(xx,ro,ro1,is,tknot,nu)

              g=0.d0
              f=0.d0

              if(is.eq.ns)then
                do i=2,ns
                  g=g+wave(i-1,kk)*ro(i)
                  f=f+wave(nm+i-1,kk)*ro(i)
                enddo
              elseif(is.eq.nu-ns)then
                do i=1,ns-1
                  g=g+wave(i-1+is-ns,kk)*ro(i)
                  f=f+wave(nm+i-1+is-ns,kk)*ro(i)
                enddo
              else
                do i=1,ns
                  g=g+wave(i-1+is-ns,kk)*ro(i)
                  f=f+wave(nm+i-1+is-ns,kk)*ro(i)
                enddo
              endif
              g=g/xx
              f=f/xx
              num=num+1
              do j=0,nugl
                ang=pi-j*pi/nugl
                dc=dcos(ang)
                ds=dsin(ang)
                gfact1=clebsh(al1,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &          spharm(l1,idint(amu-5.d-1),dc)
                gfact2=clebsh(al1,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &          spharm(l1,idint(amu+5.d-1),dc)
                ffact1=clebsh(al2,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &          spharm(l2,idint(amu-5.d-1),dc)
                ffact2=clebsh(al2,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &          spharm(l2,idint(amu+5.d-1),dc)

                if((xx.ne.plot(1,num,j)).and.(kk.ne.-nkap))then
                  write(*,*) 'something goes wrong with num'
                  write(*,*) xx,plot(1,num,j)
                  write(*,*) num,j
                  pause
                endif

                plot(1,num,j)=xx
                plot(2,num,j)=plot(2,num,j)+g*gfact1
                plot(3,num,j)=plot(3,num,j)+g*gfact2
                plot(4,num,j)=plot(4,num,j)+f*ffact1
                plot(5,num,j)=plot(5,num,j)+f*ffact2
              enddo
            enddo
          enddo
        endif
      enddo

c      do i=1,num
c         write(1,1) -plot(1,i,0),(plot(k,i,0),k=2,5),
c     &        plot(2,i,0)**2+plot(3,i,0)**2+plot(4,i,0)**2+
c     &        plot(5,i,0)**2
c      enddo
c      do i=num,1,-1
c         write(1,1) plot(1,i,nugl),(plot(k,i,nugl),k=2,5),
c     &        plot(2,i,nugl)**2+plot(3,i,nugl)**2+plot(4,i,nugl)**2+
c     &        plot(5,i,nugl)**2
c      enddo
      do i=1,num
        r=plot(1,i,0)
        do j=0,nugl
          ang=pi-j*pi/nugl
          dc=dcos(ang)
          ds=dsin(ang)
          x=r*dc
          y=r*ds
          if((dabs(x).lt. 4.d0).and.(dabs(y).lt. 4.d0))then
c               dens=plot(2,i,j)**2+plot(3,i,j)**2+plot(4,i,j)**2+
c     &              plot(5,i,j)**2
            write(511,11) x,y,plot(2,i,j),plot(3,i,j),plot(4,i,j),
     &      plot(5,i,j)
c               write(12,*) x,y,dlog(dens)
            if((j.ne.0).and.(j.ne.nugl))then
              ang=pi+j*pi/nugl
              dc=dcos(ang)
              ds=dsin(ang)
              x=r*dc
              y=r*ds
              write(511,11) x,y,plot(2,i,j),plot(3,i,j),plot(4,i,j),
     &        plot(5,i,j)
              write(512,11) x,y,dlog(plot(2,i,j)),dlog(plot(3,i,j)),
     &        dlog(plot(4,i,j)),dlog(plot(5,i,j))
            endif
          endif
        enddo
      enddo

      deallocate(plot)
 11   format(6(e12.5,1x))
 1    format(e14.6,2x,e14.6,2x,e14.6,2x,e14.6,2x,e14.6,2x,e14.6)
      return
      end

      subroutine draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      include 'inc.par'
      real*8 ro(ns),ro1(ns),wave(2*nm,2*nm,-nkap:nkap),tknot(nu)
      common /momentum_projection/ amu,amj_max

      if(kk.gt.0)then
        l1=kk
        l2=kk-1
      else
        l1=-kk-1
        l2=-kk
      endif
      al1=dble(l1)
      al2=dble(l2)
      aj=iabs(kk)-5.d-1

      gfact1=clebsh(al1,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &spharm(l1,idint(amu-5.d-1),-1.d0)
      gfact2=clebsh(al1,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &spharm(l1,idint(amu+5.d-1),-1.d0)
      ffact1=clebsh(al2,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &spharm(l2,idint(amu-5.d-1),-1.d0)
      ffact2=clebsh(al2,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &spharm(l2,idint(amu+5.d-1),-1.d0)
      do is=nu-ns,ns,-1
        do ip=ns,1,-1
          xx=tknot(is)+(tknot(is+1)-tknot(is))*dble(ip)/dble(ns)
          call dsplines(xx,ro,ro1,is,tknot,nu)

          g=0.d0
          f=0.d0

          if(is.eq.ns)then
            do i=2,ns
              g=g+wave(i-1,jj,kk)*ro(i)
              f=f+wave(nm+i-1,jj,kk)*ro(i)
            enddo
          elseif(is.eq.nu-ns)then
            do i=1,ns-1
              g=g+wave(i-1+is-ns,jj,kk)*ro(i)
              f=f+wave(nm+i-1+is-ns,jj,kk)*ro(i)
            enddo
          else
            do i=1,ns
              g=g+wave(i-1+is-ns,jj,kk)*ro(i)
              f=f+wave(nm+i-1+is-ns,jj,kk)*ro(i)
            enddo
          endif
          g=g/xx
          f=f/xx
          write(1,1) -xx,g*gfact1,g*gfact2,f*ffact1,f*ffact2,
     &    g**2*(gfact1**2+gfact2**2)+
     &    f**2*(ffact1**2+ffact2**2)
         enddo
      enddo

      gfact1=clebsh(al1,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &spharm(l1,idint(amu-5.d-1),1.d0)
      gfact2=clebsh(al1,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &spharm(l1,idint(amu+5.d-1),1.d0)
      ffact1=clebsh(al2,amu-5.d-1,5.d-1,5.d-1,aj,amu)*
     &spharm(l2,idint(amu-5.d-1),1.d0)
      ffact2=clebsh(al2,amu+5.d-1,5.d-1,-5.d-1,aj,amu)*
     &spharm(l2,idint(amu+5.d-1),1.d0)
      do is=ns,nu-ns
        do ip=1,ns
          xx=tknot(is)+(tknot(is+1)-tknot(is))*dble(ip)/dble(ns)
          call dsplines(xx,ro,ro1,is,tknot,nu)

          g=0.d0
          f=0.d0

          if(is.eq.ns)then
            do i=2,ns
              g=g+wave(i-1,jj,kk)*ro(i)
              f=f+wave(nm+i-1,jj,kk)*ro(i)
            enddo
          elseif(is.eq.nu-ns)then
            do i=1,ns-1
              g=g+wave(i-1+is-ns,jj,kk)*ro(i)
              f=f+wave(nm+i-1+is-ns,jj,kk)*ro(i)
            enddo
          else
            do i=1,ns
              g=g+wave(i-1+is-ns,jj,kk)*ro(i)
              f=f+wave(nm+i-1+is-ns,jj,kk)*ro(i)
            enddo
          endif
          g=g/xx
          f=f/xx
          write(1,1) xx,g*gfact1,g*gfact2,f*ffact1,f*ffact2,
     &    g**2*(gfact1**2+gfact2**2)+
     &    f**2*(ffact1**2+ffact2**2)
        enddo
      enddo
 1    format(e14.6,2x,e14.6,2x,e14.6,2x,e14.6,2x,e14.6,2x,e14.6)
      return
      end

      subroutine r_diagonal_simple(nm,a,eignum,eigvect)
      implicit real*8(a-h,o-z)
      real*8 a(nm,nm),eigvect(nm,nm)
      real*8 VL,VU,ABSTOL,eignum(nm)
      integer, dimension(:),allocatable:: ISUPPZ,IWORK
      character*1 JOBZ,RANGE,UPLO
      real*8, dimension(:),allocatable:: WORK
c$$$      real*8, dimension(:,:),allocatable:: atmp1,WORK(26*nm)
c$$$
c$$$      allocate(atmp1(nm,nm))
c$$$      atmp1=a

      allocate(ISUPPZ(2*nm))
      LWORK=26*nm
      allocate(WORK(LWORK))
      LIWORK=10*nm
      allocate(IWORK(LIWORK))
      JOBZ='V'
      RANGE='A'
      UPLO='U'
      VL=-1.d3
      VU=1.d3
      IL=1
      IU=nm
      ABSTOL=1.d-16

      LWORK=-1
      LIWORK=-1
      call DSYEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   IWORK, LIWORK, INFO )

      LWORK=idint(WORK(1))
      LIWORK=IWORK(1)
      deallocate(IWORK)
      deallocate(WORK)
      allocate(WORK(LWORK))
      allocate(IWORK(LIWORK))

      call DSYEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   IWORK, LIWORK, INFO )

      deallocate(ISUPPZ)
      deallocate(WORK)
      deallocate(IWORK)
cc      write(*,*) 'EXITING R_DIAG'
      return
      end


      subroutine r_diagonal(choice,nm,a,eignum,eigvect)
      implicit real*8(a-h,o-z)
      real*8 a(nm,nm),eigvect(nm,nm)
      real*8 VL,VU,ABSTOL,eignum(nm)
      integer, dimension(:),allocatable:: ISUPPZ,IWORK
      character*1 JOBZ,RANGE,UPLO
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8, dimension(:),allocatable:: WORK
      character choice
c$$$      real*8, dimension(:,:),allocatable:: atmp1,WORK(26*nm)
c$$$
c$$$      allocate(atmp1(nm,nm))
c$$$      atmp1=a

      allocate(ISUPPZ(2*nm))
      LWORK=26*nm
      allocate(WORK(LWORK))
      LIWORK=10*nm
      allocate(IWORK(LIWORK))
      JOBZ='V'
      RANGE='A'
      UPLO='U'
      VL=-1.d3
      VU=1.d3
      IL=1
      IU=nm
      ABSTOL=1.d-16

      LWORK=-1
      LIWORK=-1
      call DSYEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   IWORK, LIWORK, INFO )

      LWORK=idint(WORK(1))
      LIWORK=IWORK(1)
      deallocate(IWORK)
      deallocate(WORK)
      allocate(WORK(LWORK))
      allocate(IWORK(LIWORK))

      call DSYEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   IWORK, LIWORK, INFO )

c      write(*,*) M,'found vectors'

cc      write(*,*) 'ground state energy is'
      jj=1
      do while (dabs(eignum(jj)).gt.1.d0)
        jj=jj+1
      enddo
cc      write(*,*) eignum(jj)

      select case(choice)
      case('a')
        open(1,file='energies.dat')
      case('e')
        open(1,file='energies_even.dat')
      case('o')
        open(1,file='energies_odd.dat')
      end select
      write(1,*) 'NUMERICAL  VALUE ',eignum(jj)
      if((z_nuc1.ne.0.d0).and.(z_nuc2.eq.0.d0))then
        write(1,*) 'ANALYTICAL VALUE ',dsqrt(1.d0-az1**2)
        write(1,*) 'RELATIVE ACCURACY',
     &  dabs((dsqrt(1.d0-az1**2)-eignum(jj))/eignum(jj))
      elseif((z_nuc2.ne.0.d0).and.(z_nuc1.eq.0.d0))then
        write(1,*) 'ANALYTICAL VALUE',dsqrt(1.d0-az2**2)
        write(1,*) 'RELATIVE ACCURACY',
     &  dabs((dsqrt(1.d0-az2**2)-eignum(jj))/eignum(jj))
      endif

      write(1,*) eignum
      close(1)

c$$$      write(*,*) 'Testing eigenvectors'
c$$$      do i=1,nm
c$$$         write(*,*) i,'-th vector'
c$$$         do j=1,nm
c$$$            sum=0.d0
c$$$            do k=1,nm
c$$$               sum=sum+atmp1(j,k)*eigvect(k,i)
c$$$            enddo
c$$$            write(*,*) sum,eignum(i)*eigvect(j,i),
c$$$     &           dabs((sum-eignum(i)*eigvect(j,i))/sum)
c$$$         enddo
c$$$         pause
c$$$      enddo

      deallocate(ISUPPZ)
      deallocate(WORK)
      deallocate(IWORK)
cc      write(*,*) 'EXITING R_DIAG'
      return
      end

      subroutine c_diagonal_zheevx(nm,a,eignum,eigvect)
      implicit complex*16(a-h,o-z)
      complex*16 a(nm,nm),eigvect(nm,nm)
      real*8 VL,VU,ABSTOL,eignum(nm),aa
      integer INFO
      integer, dimension(:),allocatable:: IWORK,IFAIL
      character*1 JOBZ,RANGE,UPLO
      complex*16, dimension(:),allocatable:: WORK
      real*8, dimension(:),allocatable:: RWORK

      LWORK=2*nm
      allocate(WORK(LWORK))
      LRWORK=7*nm
      allocate(RWORK(LRWORK))
      LIWORK=5*nm
      allocate(IWORK(LIWORK))
      allocate(IFAIL(nm))
      JOBZ='V'
      RANGE='A'
      UPLO='U'
      VL=-1.d3
      VU=1.d3
      IL=1
      IU=nm
      ABSTOL=1.d-16

      LWORK=-1

      call ZHEEVX( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm,
     &                   WORK, LWORK,
     &                   RWORK, IWORK, IFAIL, INFO )

      aa=WORK(1)
      LWORK=idint(aa)
      deallocate(WORK)
      allocate(WORK(LWORK))
c      write(*,*) M,'found vectors'
      call ZHEEVX( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm,
     &                   WORK, LWORK,
     &                   RWORK, IWORK, IFAIL, INFO )

      deallocate(WORK)
      deallocate(RWORK)
      deallocate(IWORK)
      deallocate(IFAIL)

      return
      end

      subroutine c_diagonal(nm,a,eignum,eigvect)
      implicit complex*16(a-h,o-z)
      complex*16 a(nm,nm),eigvect(nm,nm)
      real*8 VL,VU,ABSTOL,eignum(nm),aa
      integer, dimension(:),allocatable:: ISUPPZ,IWORK
      character*1 JOBZ,RANGE,UPLO
      complex*16, dimension(:),allocatable:: WORK
      real*8, dimension(:),allocatable:: RWORK

      LWORK=2*nm
      allocate(WORK(LWORK))
      LRWORK=24*nm
      allocate(RWORK(LRWORK))
      allocate(ISUPPZ(2*nm))
      LIWORK=10*nm
      allocate(IWORK(LIWORK))
      JOBZ='V'
      RANGE='A'
      UPLO='U'
      VL=-1.d3
      VU=1.d3
      IL=1
      IU=nm
      ABSTOL=1.d-16

      LWORK=-1
      LRWORK=-1
      LIWORK=-1

      call ZHEEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   RWORK, LRWORK, IWORK, LIWORK, INFO )

      aa=WORK(1)
      LWORK=idint(aa)
      LRWORK=idint(RWORK(1))
      LIWORK=IWORK(1)
      deallocate(WORK)
      deallocate(RWORK)
      deallocate(IWORK)
      allocate(WORK(LWORK))
      allocate(RWORK(LRWORK))
      allocate(IWORK(LIWORK))
c      write(*,*) M,'found vectors'
      call ZHEEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   RWORK, LRWORK, IWORK, LIWORK, INFO )

      deallocate(WORK)
      deallocate(RWORK)
      deallocate(IWORK)
      deallocate(ISUPPZ)

      return
      end

      subroutine c_diagonal_lower(nm,a,eignum,eigvect)
      implicit complex*16(a-h,o-z)
      complex*16 a(nm,nm),eigvect(nm,nm)
      real*8 VL,VU,ABSTOL,eignum(nm),aa
      integer, dimension(:),allocatable:: ISUPPZ,IWORK
      character*1 JOBZ,RANGE,UPLO
      complex*16, dimension(:),allocatable:: WORK
      real*8, dimension(:),allocatable:: RWORK

      LWORK=2*nm
      allocate(WORK(LWORK))
      LRWORK=24*nm
      allocate(RWORK(LRWORK))
      allocate(ISUPPZ(2*nm))
      LIWORK=10*nm
      allocate(IWORK(LIWORK))
      JOBZ='V'
      RANGE='A'
      UPLO='L'
      VL=-1.d3
      VU=1.d3
      IL=1
      IU=nm
      ABSTOL=1.d-16

      LWORK=-1
      LRWORK=-1
      LIWORK=-1

      call ZHEEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   RWORK, LRWORK, IWORK, LIWORK, INFO )

      aa=WORK(1)
      LWORK=idint(aa)
      LRWORK=idint(RWORK(1))
      LIWORK=IWORK(1)
      deallocate(WORK)
      deallocate(RWORK)
      deallocate(IWORK)
      allocate(WORK(LWORK))
      allocate(RWORK(LRWORK))
      allocate(IWORK(LIWORK))
c      write(*,*) M,'found vectors'
      call ZHEEVR( JOBZ, RANGE, UPLO, nm, A, nm, VL, VU, IL, IU,
     &                   ABSTOL, M, eignum,eigvect, nm, ISUPPZ,
     &                   WORK, LWORK,
     &                   RWORK, LRWORK, IWORK, LIWORK, INFO )

      deallocate(WORK)
      deallocate(RWORK)
      deallocate(IWORK)
      deallocate(ISUPPZ)

      return
      end

      subroutine plot_functions(nstates,nm,nkap,num_st,nu,
     &wave,e,tknot,eigval,eigvec)
      include 'inc.par'
      real*8 wave(2*nm,2*nm,-nkap:nkap),e(2*nm,-nkap:nkap),tknot(nu),
     &eigval(nstates),eigvec(nstates,nstates)
      integer num_st(-nkap:nkap,2*nm),dist_as_integer
      common /momentum_projection/ amu,amj_max
      common /r_nuc/ r01,r02
      common /dist/distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /nuc_mod/ nuc_model
      real*8, dimension(:,:),allocatable:: wcf
      character*5 dist
      allocate (wcf(2*nm,-nkap:nkap))
      wcf=0.d0

      dist_as_integer=int(2.d0*distance/0.0025896063d0)
      write(dist,173)dist_as_integer
  173 format(I5.5)

      kk=-1
      open(1,file='plot_rad_ground_'//dist//'.dat')
      jj=1
      do while(e(jj,kk).lt.-1.d0)
        jj=jj+1
      enddo

      write(1,*) '# energy',e(jj,kk),'number',jj

      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)

      close(1)

      open(1,file='plot_rad_2s_'//dist//'.dat')
      jj=jj+1
      write(1,*) '# energy',e(jj,kk),'number',jj

      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)

      kk=1
      open(1,file='plot_rad_2p1ground_'//dist//'.dat')
      jj=1
      do while(e(jj,kk).lt.-1.d0)
        jj=jj+1
      enddo
      write(1,*) '# energy',e(jj,kk),'number',jj
      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)

      open(1,file='plot_rad_2p1_'//dist//'.dat')
      jj=jj+1
      write(1,*) '# energy',e(jj,kk),'number',jj
      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)

      kk=-2
      open(1,file='plot_rad_2p3ground_'//dist//'.dat')
      jj=1
      do while(e(jj,kk).lt.-1.d0)
        jj=jj+1
      enddo
      write(1,*) '# energy',e(jj,kk),'number',jj
      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)

      open(11,file='dens_1sigma_'//dist//'.dat')
      open(12,file='dens_log_1sigma_'//dist//'.dat')
      open(511,file='dens_1sigma_'//dist//'_raw.dat')
      open(512,file='dens_log_1sigma_'//dist//'_raw.dat')
      open(1,file='plot_1sigma_'//dist//'.dat')
      open(2,file='states_1sigma_'//dist//'.dat')
      write(2,*) '#coe, N,kap,ene'

      jj=1
      do while(eigval(jj).lt.-1.d0)
        jj=jj+1
      enddo

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='states_1sigma_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)
      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_2sigma_'//dist//'.dat')
      open(12,file='dens_log_2sigma_'//dist//'.dat')
      open(511,file='dens_2sigma_'//dist//'_raw.dat')
      open(512,file='dens_log_2sigma_'//dist//'_raw.dat')
      open(2,file='wave_2sigma_'//dist//'.dat')
      write(2,*) '#coe, N,kap,parity,ene'
      open(1,file='plot_2sigma_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_2sigma_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)

      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_3sigma_'//dist//'.dat')
      open(12,file='dens_log_3sigma_'//dist//'.dat')
      open(511,file='dens_3sigma_'//dist//'_raw.dat')
      open(512,file='dens_log_3sigma_'//dist//'_raw.dat')
      open(2,file='states_3sigma_'//dist//'.dat')
      write(2,*) '#coe, N,kap,ene'
      open(1,file='plot_3sigma_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_3sigma_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)

      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_4sigma_'//dist//'.dat')
      open(12,file='dens_log_4sigma_'//dist//'.dat')
      open(511,file='dens_4sigma_'//dist//'_raw.dat')
      open(512,file='dens_log_4sigma_'//dist//'_raw.dat')
      open(2,file='states_4sigma_'//dist//'.dat')
      write(2,*) '#coe, N,kap,ene'
      open(1,file='plot_4sigma_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_4sigma_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0
      close(1)
      close(2)
      close(11)
      close(12)
      close(511)
      close(512)
      deallocate(wcf)
      return
      end

      subroutine write_output(ns,nstates,nm,nkap,num_st,nu,
     &wave,e,tknot,eigval,eigvec)
      real*8 wave(2*nm,2*nm,-nkap:nkap),e(2*nm,-nkap:nkap),tknot(nu),
     &eigval(nstates),eigvec(nstates,nstates)
      integer num_st(-nkap:nkap,2*nm)
      common /momentum_projection/ amu,amj_max
      common /r_nuc/ r01,r02
      common /dist/distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /nuc_mod/ nuc_model

      open(1,file='wave_functions_monopol.out',form='unformatted')
      write(1) wave
      close(1)
      open(1,file='wave_functions_CI.out',form='unformatted')
      write(1) eigvec
      close(1)
      open(1,file='energies_monopol.out',form='unformatted')
      write(1) e
      close(1)
      open(1,file='B_spline_knots.out',form='unformatted')
      write(1) tknot
      close(1)
      open(1,file='energies_CI.out',form='unformatted')
      write(1) eigval
      close(1)
      open(1,file='numbers_monopol.out',form='unformatted')
      write(1) num_st
      close(1)

      open(1,file='parameters.out')
      write(1,*) ns,nstates,nm,nkap,nu,amu,r01,r02,distance,
     &z_nuc1,az1,z_nuc2,az2,nuc_model
      close(1)
      return
      end

      subroutine plot_functions_even(nstates,nm,nkap,num_st,nu,
     &wave,e,tknot,eigval,eigvec)
      include 'inc.par'
      real*8 wave(2*nm,2*nm,-nkap:nkap),e(2*nm,-nkap:nkap),tknot(nu),
     &eigval(nstates),eigvec(nstates,nstates)
      integer num_st(-nkap:nkap,2*nm),dist_as_integer
      common /momentum_projection/ amu,amj_max
      common /r_nuc/ r01,r02
      common /dist/distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /nuc_mod/ nuc_model
      real*8, dimension(:,:),allocatable:: wcf
      character*5 dist
      allocate (wcf(2*nm,-nkap:nkap))
      wcf=0.d0

      dist_as_integer=int(2.d0*distance/0.0025896063d0)
      write(dist,173)dist_as_integer
  173 format(I5.5)

      kk=-1
      open(1,file='plot_rad_ground_'//dist//'.dat')
      jj=1
      do while(e(jj,kk).lt.-1.d0)
        jj=jj+1
      enddo

      write(1,*) '# energy',e(jj,kk),'number',jj

      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)

      close(1)

      open(1,file='plot_rad_2s_'//dist//'.dat')
      jj=jj+1
      write(1,*) '# energy',e(jj,kk),'number',jj

      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)

      kk=1
      open(1,file='plot_rad_2p1ground_'//dist//'.dat')
      jj=1
      do while(e(jj,kk).lt.-1.d0)
        jj=jj+1
      enddo
      write(1,*) '# energy',e(jj,kk),'number',jj
      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)

      open(1,file='plot_rad_2p1_'//dist//'.dat')
      jj=jj+1
      write(1,*) '# energy',e(jj,kk),'number',jj
      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)

      kk=-2
      open(1,file='plot_rad_2p3ground_'//dist//'.dat')
      jj=1
      do while(e(jj,kk).lt.-1.d0)
        jj=jj+1
      enddo
      write(1,*) '# energy',e(jj,kk),'number',jj
      call draw_radial(nm,nkap,nu,tknot,wave,jj,kk)
      close(1)


      open(1,file='plot_1sigma_g_'//dist//'.dat')
      open(2,file='states_1sigma_g_'//dist//'.dat')
      write(2,*) '#coe, N,kap,ene'

      jj=1
      do while(eigval(jj).lt.-1.d0)
        jj=jj+1
      enddo

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(11,file='dens_1sigma_g_'//dist//'.dat')
      open(12,file='dens_log_1sigma_g_'//dist//'.dat')
      open(511,file='dens_1sigma_g_'//dist//'_raw.dat')
      open(512,file='dens_log_1sigma_g_'//dist//'_raw.dat')
      open(3,file='states_1sigma_g_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)
      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_2sigma_g_'//dist//'.dat')
      open(12,file='dens_log_2sigma_g_'//dist//'.dat')
      open(511,file='dens_2sigma_g_'//dist//'_raw.dat')
      open(512,file='dens_log_2sigma_g_'//dist//'_raw.dat')
      open(2,file='wave_2sigma_g_'//dist//'.dat')
      write(2,*) '#coe, N,kap,parity,ene'
      open(1,file='plot_2sigma_g_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_2sigma_g_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)

      close(2)

      close(511)
      close(512)
      close(11)
      close(12)
      open(11,file='dens_3sigma_g_'//dist//'.dat')
      open(511,file='dens_3sigma_g_'//dist//'_raw.dat')
      open(2,file='states_3sigma_g_'//dist//'.dat')
      open(12,file='dens_log_3sigma_g_'//dist//'.dat')
      open(512,file='dens_log_3sigma_g_'//dist//'_raw.dat')
      write(2,*) '#coe, N,kap,ene'
      open(1,file='plot_3sigma_g_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_3sigma_g_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)

      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_4sigma_g_'//dist//'.dat')
      open(511,file='dens_4sigma_g_'//dist//'_raw.dat')
      open(2,file='states_4sigma_g_'//dist//'.dat')
      open(12,file='dens_log_4sigma_g_'//dist//'.dat')
      open(512,file='dens_log_4sigma_g_'//dist//'_raw.dat')
      write(2,*) '#coe, N,kap,ene'
      open(1,file='plot_4sigma_g_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_4sigma_g_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0
      close(1)
      close(2)
      close(511)
      close(512)
      close(11)
      close(12)
      deallocate(wcf)
      return
      end

      subroutine write_output_even(ns,nstates,nm,nkap,num_st,nu,
     &wave,e,tknot,eigval,eigvec)
      real*8 wave(2*nm,2*nm,-nkap:nkap),e(2*nm,-nkap:nkap),tknot(nu),
     &eigval(nstates),eigvec(nstates,nstates)
      integer num_st(-nkap:nkap,2*nm)
      common /momentum_projection/ amu,amj_max
      common /r_nuc/ r01,r02
      common /dist/distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /nuc_mod/ nuc_model

      open(1,file='wave_functions_monopol.out',form='unformatted')
      write(1) wave
      close(1)
      open(1,file='wave_functions_CI_g.out',form='unformatted')
      write(1) eigvec
      close(1)
      open(1,file='energies_monopol.out',form='unformatted')
      write(1) e
      close(1)
      open(1,file='B_spline_knots.out',form='unformatted')
      write(1) tknot
      close(1)
      open(1,file='energies_CI_g.out',form='unformatted')
      write(1) eigval
      close(1)
      open(1,file='numbers_monopol_g.out',form='unformatted')
      write(1) num_st
      close(1)

      open(1,file='parameters.out')
      write(1,*) ns,nstates,nm,nkap,nu,amu,r01,r02,distance,
     &z_nuc1,az1,z_nuc2,az2,nuc_model
      close(1)
      return
      end

      subroutine plot_functions_odd(nstates,nm,nkap,num_st,nu,
     &wave,e,tknot,eigval,eigvec)
      include 'inc.par'
      real*8 wave(2*nm,2*nm,-nkap:nkap),e(2*nm,-nkap:nkap),tknot(nu),
     &eigval(nstates),eigvec(nstates,nstates)
      integer num_st(-nkap:nkap,2*nm),dist_as_integer
      common /momentum_projection/ amu,amj_max
      common /r_nuc/ r01,r02
      common /dist/distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /nuc_mod/ nuc_model
      real*8, dimension(:,:),allocatable:: wcf
      character*5 dist
      allocate (wcf(2*nm,-nkap:nkap))
      wcf=0.d0

      dist_as_integer=int(2.d0*distance/0.0025896063d0)
      write(dist,173)dist_as_integer
  173 format(I5.5)

      open(1,file='plot_1sigma_u_'//dist//'.dat')
      open(2,file='states_1sigma_u_'//dist//'.dat')
      write(2,*) '#coe, N,kap,ene'

      jj=1
      do while(eigval(jj).lt.-1.d0)
        jj=jj+1
      enddo

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(11,file='dens_1sigma_u_'//dist//'.dat')
      open(12,file='dens_log_1sigma_u_'//dist//'.dat')
      open(511,file='dens_1sigma_u_'//dist//'_raw.dat')
      open(512,file='dens_log_1sigma_u_'//dist//'_raw.dat')
      open(3,file='states_1sigma_u_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)
      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_2sigma_u_'//dist//'.dat')
      open(12,file='dens_log_2sigma_u_'//dist//'.dat')
      open(511,file='dens_2sigma_u_'//dist//'_raw.dat')
      open(512,file='dens_log_2sigma_u_'//dist//'_raw.dat')
      open(2,file='wave_2sigma_u_'//dist//'.dat')
      write(2,*) '#coe, N,kap,parity,ene'
      open(1,file='plot_2sigma_u_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_2sigma_u_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)

      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_3sigma_u_'//dist//'.dat')
      open(12,file='dens_log_3sigma_u_'//dist//'.dat')
      open(511,file='dens_3sigma_u_'//dist//'_raw.dat')
      open(512,file='dens_log_3sigma_u_'//dist//'_raw.dat')
      open(2,file='states_3sigma_u_'//dist//'.dat')
      write(2,*) '#coe, N,kap,ene'
      open(1,file='plot_3sigma_u_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_3sigma_u_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0

      close(1)

      close(2)

      close(11)
      close(12)
      close(511)
      close(512)
      open(11,file='dens_4sigma_u_'//dist//'.dat')
      open(12,file='dens_log_4sigma_u_'//dist//'.dat')
      open(511,file='dens_4sigma_u_'//dist//'_raw.dat')
      open(512,file='dens_log_4sigma_u_'//dist//'_raw.dat')
      open(2,file='states_4sigma_u_'//dist//'.dat')
      write(2,*) '#coe, N,kap,ene'
      open(1,file='plot_4sigma_u_'//dist//'.dat')
      jj=jj+1

      write(1,*) '# energy',eigval(jj),'number',jj

      do i=1,2*nm
        do j=-nkap,nkap
          if(num_st(j,i).ne.0.d0)then
            n=num_st(j,i)
            do l=1,2*nm
              wcf(l,j)=wcf(l,j)+eigvec(n,jj)*wave(l,i,j)
            enddo
            if(dabs(eigvec(n,jj)).gt.1.d-12)then
              if(((j.gt.0).and.((j/2)*2.ne.j)).or.
     &        ((j.lt.0).and.((j/2)*2.eq.j)))then
                write(2,*) eigvec(n,jj),i,j,'odd ',e(i,j)
              else
                write(2,*) eigvec(n,jj),i,j,'even',e(i,j)
              endif
            endif
          endif
        enddo
      enddo

      open(3,file='wave_4sigma_u_'//dist//'.dat')
      write(3,*) nm,nkap,nu,tknot,wcf,eigval(jj)
      close(3)
      call draw_sigma(nm,nkap,nu,tknot,wcf)
      call draw_sigma_raw(nm,nkap,nu,tknot,wcf)
      wcf=0.d0
      close(1)
      close(2)
      close(11)
      close(12)
      close(511)
      close(512)
      deallocate(wcf)
      return
      end

      subroutine write_output_odd(nstates,nm,nkap,num_st,
     &eigval,eigvec)
      real*8 eigval(nstates),eigvec(nstates,nstates)
      integer num_st(-nkap:nkap,2*nm)
      common /momentum_projection/ amu,amj_max
      common /r_nuc/ r01,r02
      common /dist/distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /nuc_mod/ nuc_model

      open(1,file='wave_functions_CI_u.out',form='unformatted')
      write(1) eigvec
      close(1)
      open(1,file='energies_CI_u.out',form='unformatted')
      write(1) eigval
      close(1)
      open(1,file='numbers_monopol_u.out',form='unformatted')
      write(1) num_st
      close(1)

      return
      end

c Added in order to get the states better defined 110810.

      function State_orgnzr(i,kk,max_kappa)
      integer State_orgnzr,i,kk,max_kappa
      if (kk .ne. 0) then
        State_orgnzr=Invtd_Kppa_arry(kk)+
     &  Invtd_Kppa_arry(max_kappa)*(i-1)
      endif
      return
      end

      function Invtd_Kppa_arry(kk)
      integer Invtd_Kppa_arry,kk
      if (kk .ne. 0) then
        if (kk .gt. 0) then
          Invtd_Kppa_arry=2*kk
        elseif (kk .lt. 0) then
          Invtd_Kppa_arry=2*abs(kk)-1
        endif
      endif
      return
      end
