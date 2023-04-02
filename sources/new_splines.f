!     This subroutine calculates the wave functions and energies

      subroutine b_spline_calculation(nstates,nsto,nste,
     &nm,nu,nkap,number_states,rmin,rmax,wave1,vmat,nvmat,
     &e,up_energy,dmat1,dvdRmat_dkb1,dvdRmat_dkb2,
     &alternate_dmat,v2sbj1a,
     &v2sbj2a,v2sbj3a,v2sbj4a,wavenumb)
      include 'inc.par'
      real*8, dimension(:),allocatable:: u,sbj,cc,t
      real*8, dimension(:,:),allocatable:: amat,wave,dmat,b,b_2,dd,d3
      real*8, dimension(:,:,:),allocatable:: div1,div2,v,
     &v_2,dvd,dv_1,dvb,v2sbj1,v2
      real*8, dimension(:,:,:,:),allocatable:: v2sbj2,v2sbj3,
     &v_22b,v_22d
      real*8, dimension(:,:,:,:,:),allocatable:: v2sbj4,v_22smalla,
     &v_22smallb
      real*8 ro(ns),ro1(ns),ro2a(ns)
      real*8 wave1(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),e(2*nm,-nkap:nkap),
     &dmat1(2*nm,2*nm),dvdRmat_dkb1(nm,nm,-nkap:nkap,0:2*nkap,2),
     &dvdRmat_dkb2(nm,nm,-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &alternate_dmat(2*nm,2*nm,-nkap:nkap),v2sbj1a(nm,nm,0:2*nkap),
     &v2sbj2a(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj3a(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj4a(-nkap:nkap,-nkap:nkap,nm,nm,0:2*nkap)
      common /weights/ w4n(4),t4n(4),w8(8),t8(8),w16(16),t16(16),
     &w32(32),t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 ttt(nuz),www(nuz),dv_dR(0:2*nkap),el(2*nm)
      common /dist/distance,Starting_Distance
      integer number_states(2,-nkap:nkap)
      integer One_s_state_location,xi_stepslower,xi_stepsupper
      logical dkb,Manual_ncont_states
      common /neg_cont/ Manual_ncont_states
      common /common_dkb/ dkb
      common /Barycentres/ RadiusOne,RadiusTwo
      common /step_progress/ xi_stepslower,xi_stepsupper,ii_xi
      common /TargProj_prop/ Proj_mass,Targ_mass,Proj_vel,
     &b_ImpactParam

      pi4=1.d0/(2.d0)

      charge_radius=Targ_mass**(1.d0/3.d0)*1.2d0*0.0025896063d0
      Starting_R1=2.d0*Starting_Distance*Proj_mass/(Proj_mass+Targ_mass)
      Starting_R2=2.d0*Starting_Distance*Targ_mass/(Proj_mass+Targ_mass)

      nul=nu
      alternate_dmat=0.d0
      dmat1=0.d0
      dvdRmat_dkb1=0.d0
      dvdRmat_dkb2=0.d0
      e=0.d0

c b2,v2, are added for the
c calculation of the dV/dR matrix

      allocate(amat(2*nm,2*nm))
      allocate(dmat(2*nm,2*nm))
      allocate(cc(2*nm))
      allocate(wave(2*nm,2*nm))
      allocate(b(nm+2,nm+2))
c      allocate(b2(nm+2,nm+2,0:2*nkap))
      allocate(div2(nm+2,nm+2,-nkap:nkap))
      allocate(div1(nm+2,nm+2,-nkap:nkap))
      allocate(v(nm+2,nm+2,0:2*nkap))
      allocate(v2(nm+2,nm+2,0:2*nkap))
      allocate(v2sbj1(nm+2,nm+2,0:2*nkap))

      allocate(sbj(0:2*nkap))
      allocate(u(0:2*nkap))
      if(dkb)then
        allocate(b_2(nm+2,nm+2))
        allocate(dd(nm+2,nm+2))
        allocate(d3(nm+2,nm+2))
        allocate(v_2(nm+2,nm+2,0:2*nkap))
        allocate(v_22b(-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(v_22d(-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(v_22smalla(-nkap:nkap,-nkap:nkap,nm+2,
     &  nm+2,0:2*nkap))
        allocate(v_22smallb(-nkap:nkap,-nkap:nkap,nm+2,
     &  nm+2,0:2*nkap))
        allocate(dvd(nm+2,nm+2,0:2*nkap))
        allocate(dv_1(nm+2,nm+2,0:2*nkap))
        allocate(dvb(nm+2,nm+2,0:2*nkap))
        allocate(v2sbj2(-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(v2sbj3(-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(v2sbj4(-nkap:nkap,-nkap:nkap,nm+2,nm+2,0:2*nkap))
      endif

      call getIntegralWeights(ttt, www)

!      select case (nuz)
!      case(4)
!        do i=1,4
!          ttt(i)=t4n(i)
!          www(i)=w4n(i)
!        enddo
!      case(8)
!        do i=1,8
!          ttt(i)=t8(i)
!          www(i)=w8(i)
!        enddo
!      case(6)
!        do i=1,6
!          ttt(i)=t6(i)
!          www(i)=w6(i)
!        enddo
!      case(16)
!        do i=1,16
!          ttt(i)=t16(i)
!          www(i)=w16(i)
!        enddo
!      case(32)
!        do i=1,32
!          ttt(i)=t32(i)
!          www(i)=w32(i)
!        enddo
!      case(64)
!        do i=1,64
!          ttt(i)=t64(i)
!          www(i)=w64(i)
!        enddo
!      end select

      v=0.d0
      b=0.d0
      div1=0.d0
      div2=0.d0
      amat=0.d0
      dmat=0.d0
      cc=0.d0
      vmat=0.d0
      if(dkb) then
        b_2=0.d0
        dd=0.d0
        d3=0.d0
        v_2=0.d0
c        v_22=0.d0, v_22 is set to zero elsewhere.
        dvd=0.d0
        dv_1=0.d0
        dvb=0.d0
      endif
      allocate(t(nu))
      call spline_arranger(nu,rmin,rmax,t)

      if(z_nuc1.eq.z_nuc2)then
        istepforll=2
      else
        istepforll=1
      endif
c     The calculation of the matrix
      do inx1=ns,nm+2
        ac=(t(inx1)+t(inx1+1))/2.d0
        bc=(t(inx1+1)-t(inx1))/2.d0

        do j=nuz,1,-1
          xx=ac+bc*ttt(j)
          bbc=bc*www(j)

          if(dkb)then
            call ddsplines(xx,ro,ro1,ro2a,inx1,t,nu)
          else
            call dsplines(xx,ro,ro1,inx1,t,nu)
          endif
          if(z_nuc1.eq.z_nuc2)then
            call all_potential(nkap,xx,u)
          elseif(z_nuc1.ne.z_nuc2 )then
            call all_potential_baryonic(nkap,xx,u)
          endif

          do i=inx1-ns+1,inx1
            do k=inx1-ns+1,inx1
              b(i,k)=b(i,k)+ro(i-inx1+ns)*ro(k-inx1+ns)*bbc
              if(dkb)then
                b_2(i,k)=b_2(i,k)+
     &          ro(i-inx1+ns)*ro(k-inx1+ns)*bbc/xx/xx
                dd(i,k)=dd(i,k)+
     &          ro1(i-inx1+ns)*ro1(k-inx1+ns)*bbc
                d3(i,k)=d3(i,k)+
     &          ro1(i-inx1+ns)*ro2a(k-inx1+ns)*bbc
              endif
              do ll=0,2*nkap,istepforll
                v(i,k,ll)=v(i,k,ll)+
     &          ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*u(ll)

                if(dkb)then
                  v_2(i,k,ll)=v_2(i,k,ll)+
     &            ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*u(ll)/xx/xx
                  dvd(i,k,ll)=dvd(i,k,ll)+
     &            ro1(i-inx1+ns)*ro1(k-inx1+ns)*bbc*u(ll)
                  dv_1(i,k,ll)=dv_1(i,k,ll)+
     &            (ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &            ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &            bbc*u(ll)/xx
                  if((i.gt.1).and.(k.gt.1).and.(i-1.le.nm).and.
     &            (k-1.le.nm))then
                    vmat(i-1,k-1,ll,6)=vmat(i-1,k-1,ll,6)+
     &              (-ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &              bbc*u(ll)/xx/4.d0
                    vmat(i-1,k-1,ll,7)=vmat(i-1,k-1,ll,7)+
     &              (-ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &              bbc*u(ll)/2.d0
                    vmat(i-1,k-1,ll,8)=vmat(i-1,k-1,ll,8)+
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              bbc*u(ll)/xx/2.d0

                  endif
                  dvb(i,k,ll)=dvb(i,k,ll)+
     &            (ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &            ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &            bbc*u(ll)

                endif
              enddo
              do kk=-nkap,nkap
                if(kk.ne.0)then
                  if(dkb)then
                    div1(i,k,kk)=div1(i,k,kk)+(
     &              ro(i-inx1+ns)*ro2a(k-inx1+ns)*kk/xx-
     &              ro1(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*(kk-1)/xx/xx-
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*kk*(kk-1)/xx/xx/xx
     &              )*bbc
                    div2(i,k,kk)=div2(i,k,kk)+(
     &              ro(i-inx1+ns)*ro2a(k-inx1+ns)*kk/xx+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*(kk+1)/xx/xx-
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*kk*(kk+1)/xx/xx/xx
     &              )*bbc
                  else
                    div1(i,k,kk)=div1(i,k,kk)+(
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*kk/xx-
     &              ro(i-inx1+ns)*ro1(k-inx1+ns))*bbc
                    div2(i,k,kk)=div2(i,k,kk)+(
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*kk/xx+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*bbc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      enddo


      if(dkb) then
        do i=1,nm
          do k=1,nm
            l=i+1
            m=k+1
            do ll=0,2*nkap,istepforll
              vmat(i,k,ll,1)=v(l,m,ll)
              vmat(i,k,ll,2)=dvd(l,m,ll)/4.d0
              vmat(i,k,ll,3)=v_2(l,m,ll)/4.d0
              vmat(i,k,ll,4)=dv_1(l,m,ll)/4.d0
              vmat(i,k,ll,5)=dvb(l,m,ll)/2.d0
            enddo
          enddo
        enddo
      else
        do i=1,nm
          do k=1,nm
            l=i+1
            m=k+1
            amat(i,k)=b(l,m)+v(l,m,0)*pi4
            amat(i+nm,k+nm)=-b(l,m)+v(l,m,0)*pi4
            dmat(i,k)=b(l,m)
            dmat(i+nm,k+nm)=b(l,m)
            do ll=0,2*nkap,istepforll
              vmat(i,k,ll,1)=v(l,m,ll)
            enddo
          enddo
        enddo
      endif

c****** MATRIX OF dV/dR
c This part is the calculation of my so called dV/dR matrix element,
c similar to the above, yet somewhat altered.


c      b2=0.d0
      v2=0.d0
      v2sbj1=0.d0

      if(dkb) then
        v_22b=0.d0
        v_22d=0.d0
        v_22smalla=0.d0
        v_22smallb=0.d0
        v2sbj2=0.d0
        v2sbj3=0.d0
        v2sbj4=0.d0
      endif

      R=2.d0*distance
      lstep=1
      if(z_nuc1.eq.z_nuc2)lstep=2

      do inx1=ns,nm+2
        ac=(t(inx1)+t(inx1+1))/2.d0
        bc=(t(inx1+1)-t(inx1))/2.d0

        do j=nuz,1,-1
          xx=ac+bc*ttt(j)
          bbc=bc*www(j)

          call spherical_bessel_j(wavenumb*xx,sbj,nkap)

          if(z_nuc1.eq.z_nuc2)then
            call dvdr_potential(nkap,xx,dv_dR)
          else
            call dvdr_baryonic(nkap,xx,dv_dR)
          endif
          call dsplines(xx,ro,ro1,inx1,t,nu)

          do i=inx1-ns+1,inx1
            do k=inx1-ns+1,inx1
              do ll=0,2*nkap,lstep
                v2(i,k,ll)=v2(i,k,ll)+ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*
     &          dv_dR(ll)
              enddo
              do ll=0,2*nkap
                v2sbj1(i,k,ll)=v2sbj1(i,k,ll)+ro(i-inx1+ns)*
     &          ro(k-inx1+ns)*bbc*sbj(ll)
              enddo
            enddo
          enddo
          if(dkb) then
            do i=inx1-ns+1,inx1
              do k=inx1-ns+1,inx1
                do k1=-nkap,nkap
                  k2=k1
                  if(k1.ne. 0) then
                    do ll=0,2*nkap,lstep
                      v_22b(k1,i,k,ll)=v_22b(k1,i,k,ll)+
     &                ro(k-inx1+ns)*(ro1(i-inx1+ns)-
     &                k1*ro(i-inx1+ns)/xx)*bbc*dv_dR(ll)/2.d0

                      v_22d(k1,i,k,ll)=v_22d(k1,i,k,ll)+
     &                ro(k-inx1+ns)*(ro1(i-inx1+ns)+
     &                k1*ro(i-inx1+ns)/xx)*bbc*dv_dR(ll)/2.d0
                    enddo
                  endif
                  do k22=-nkap,nkap
                    if(k1*k22.ne. 0) then
                      do ll=0,2*nkap,lstep
                        v_22smalla(k1,k22,i,k,ll)=
     &                  v_22smalla(k1,k22,i,k,ll)+
     &                  (ro1(i-inx1+ns)*ro1(k-inx1+ns)+
     &                  k1*k22*ro(i-inx1+ns)*ro(k-inx1+ns)/xx/xx-
     &                  k1*(ro(i-inx1+ns)*ro1(k-inx1+ns))/xx-
     &                  k22*(ro1(i-inx1+ns)*ro(k-inx1+ns))/xx)*
     &                  bbc*dv_dR(ll)/4.d0

                        v_22smallb(k1,k22,i,k,ll)=
     &                  v_22smallb(k1,k22,i,k,ll)+
     &                  (ro1(i-inx1+ns)*ro1(k-inx1+ns)+
     &                  k1*k22*ro(i-inx1+ns)*ro(k-inx1+ns)/xx/xx+
     &                  k1*(ro(i-inx1+ns)*ro1(k-inx1+ns))/xx+
     &                  k22*(ro1(i-inx1+ns)*ro(k-inx1+ns))/xx)*
     &                  bbc*dv_dR(ll)/4.d0
                      enddo
                    endif
                  enddo
                enddo
                do ll=0,2*nkap
                  do k1=-nkap,nkap
                    k2=k1
                    if(k1.ne. 0) then
                      v2sbj2(k1,i,k,ll)=v2sbj2(k1,i,k,ll)+
     &                ro(k-inx1+ns)*(
     &                ro1(i-inx1+ns)+
     &                k1*ro(i-inx1+ns)/xx)*
     &                bbc*sbj(ll)/2.d0

                      v2sbj3(k1,i,k,ll)=v2sbj3(k1,i,k,ll)+
     &                ro(k-inx1+ns)*(
     &                ro1(i-inx1+ns)-
     &                k1*ro(i-inx1+ns)/xx)*
     &                bbc*sbj(ll)/2.d0
                    endif
                    do k22=-nkap,nkap
                      if(k1*k22.ne. 0) then
                        v2sbj4(k1,k22,i,k,ll)=
     &                  v2sbj4(k1,k22,i,k,ll)+
     &                  (ro1(i-inx1+ns)*ro1(k-inx1+ns)-
     &                  k1*k22*ro(i-inx1+ns)*ro(k-inx1+ns)/xx/xx+
     &                  k1*(ro(i-inx1+ns)*ro1(k-inx1+ns))/xx-
     &                  k22*(ro1(i-inx1+ns)*ro(k-inx1+ns))/xx)*
     &                  bbc*sbj(ll)/4.d0
                      endif
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
      deallocate(t)


c******* VERY IMPORTANT  ********
c the matrix elements of dV/dR are stored in vmat at positions starting from
c element(nm+1,nm+1) this was done to because I didn't want to add any extra
c terms into the subroutine b_spline_calculation
      do i=1,nm
        do k=1,nm
          l=i+1
          m=k+1
          if(dkb) then
            do kk1=-nkap,nkap
              if(kk1.ne. 0) then
                do ll=0,2*nkap,lstep
                  dvdRmat_dkb1(i,k,kk1,ll,1)=v_22b(kk1,l,m,ll)
                  dvdRmat_dkb1(i,k,kk1,ll,2)=v_22d(kk1,l,m,ll)
                enddo
                do ll=0,2*nkap
                  v2sbj2a(kk1,i,k,ll)=v2sbj2(kk1,l,m,ll)
                  v2sbj3a(kk1,i,k,ll)=v2sbj3(kk1,l,m,ll)
                enddo
              endif
              do kk2=-nkap,nkap
                if(kk2*kk1.ne. 0) then
                  do ll=0,2*nkap,lstep
                    dvdRmat_dkb2(i,k,kk1,kk2,ll,1)=
     &              v_22smalla(kk1,kk2,l,m,ll)
                    dvdRmat_dkb2(i,k,kk1,kk2,ll,2)=
     &              v_22smallb(kk1,kk2,l,m,ll)
                  enddo
                  do ll=0,2*nkap
                    v2sbj4a(kk1,kk2,i,k,ll)=v2sbj4(kk1,kk2,l,m,ll)
                  enddo
                endif
              enddo
            enddo
          endif
          do ll=0,2*nkap,lstep
            vmat(i,k,ll,nvmat)=v2(l,m,ll)
          enddo
          do ll=0,2*nkap
            v2sbj1a(i,k,ll)=v2sbj1(l,m,ll)
          enddo
        enddo
      enddo

c****** END OF MATRIX dV/dR

CCC      write(*,*) 'DEALLOCATING V'
      deallocate(v)
      deallocate(v2)
      deallocate(v2sbj1)
      nstates=0
      nsto=0
      nste=0
      number_states=0

      do kk=-nkap,nkap
        if(kk.ne.0)then
          do i=1,nm
            do k=1,nm
              l=i+1
              m=k+1

              if(dkb) then
                amat(i,k)=b(l,m)+vmat(i,k,0,1)*pi4+7.5d-1*dd(l,m)+
     &          7.5d-1*kk*(kk+1)*b_2(l,m)+vmat(i,k,0,2)*pi4+
     &          kk*kk*vmat(i,k,0,3)*pi4+kk*vmat(i,k,0,4)*pi4

                amat(i+nm,k+nm)=-b(l,m)+vmat(i,k,0,1)*pi4-
     &          7.5d-1*dd(l,m)-
     &          7.5d-1*kk*(kk-1)*b_2(l,m)+vmat(i,k,0,2)*pi4+
     &          kk*kk*vmat(i,k,0,3)*pi4-kk*vmat(i,k,0,4)*pi4

                amat(i,k+nm)=vmat(i,k,0,5)*pi4+
     &          (d3(l,m)-d3(m,l))/8.d0+
     &          (div1(l,m,kk)+div2(m,l,kk))/8.d0

                amat(k+nm,i)=amat(i,k+nm)

                dmat(i,k)=b(l,m)+dd(l,m)/4.d0+kk*b_2(l,m)/
     &          4.d0*(kk+1)
                dmat(i+nm,k+nm)=b(l,m)+dd(l,m)/4.d0+
     &          kk*b_2(l,m)/4.d0*(kk-1)

c       alternate_dmat(i,k,kk)=dmat(i,k)
c       alternate_dmat(i+nm,k+nm,kk)=dmat(i+nm,k+nm)
              else
                amat(i,k+nm)=(div1(l,m,kk)+div2(l,m,kk))/2.d0

                amat(k+nm,i)=amat(i,k+nm)
              endif
            enddo
          enddo


c DMAT should now be a visible result of calling b_splines
c Use this for checking orthogonality, comment out when not in use.

          do i=1,2*nm
            do k=1,2*nm
              dmat1(i,k)=dmat(i,k)
              alternate_dmat(i,k,kk)=dmat(i,k)
            enddo
          enddo
c End of DMAT

c     Solving of the generalized problem for eigenvectors
c I am guessing that wave is vector v in eq 16 of Johnson

          call solve_equation(nm,dmat,amat,cc,wave)

          if(cc(nm+1).gt. -1.d0) then
            do i=1,2*nm
              e(i,kk)=cc(i)
              do j=1,2*nm
                wave1(j,i,kk)=wave(j,i)
              enddo
            enddo
          elseif(cc(nm+1).lt. -1.d0) then
            nrange=nm+1
            el=0.d0
            do j=1,nm+1
              if (cc(j).gt. -up_energy .and. cc(j).lt. -1.d0) then
                do l=1,nm/2+1
                  do m=1,nm/2+1
                    el(j)=el(j)+
     &              (wave(l,j)*wave(m,j)*dmat1(l,m)+wave(l+nm,j)*
     &              wave(m+nm,j)*dmat1(l+nm,m+nm))
                  enddo
                enddo
              endif
            enddo
            One_s_state_location=maxloc(el, dim=1, mask= el .lt. 1.d0)
            do i=1,2*nm
              if (i .lt. One_s_state_location) then
                e(i,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i,kk)=wave(j,i)
                enddo
              elseif(i .eq. One_s_state_location) then
                e(nm+1,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,nm+1,kk)=wave(j,i)
                enddo
              elseif(i .gt. One_s_state_location .and.
     &        i .le. nm+1) then
                e(i-1,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i-1,kk)=wave(j,i)
                enddo
              else
                e(i,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i,kk)=wave(j,i)
                enddo
              endif
            enddo
          endif
          jjj=1
          if(Manual_ncont_states)then
             do while(jjj.le. 2*nm .and. dabs(cc(jjj)).gt.up_energy)
              jjj=jjj+1
             enddo
             number_states(1,kk)=jjj
           else
             number_states(1,kk)=nm+1
             jjj=nm+1
          endif

          do while(jjj.le.2*nm .and. cc(jjj).le.up_energy)
            if(dabs(0.5d0).le.(iabs(kk)-5.d-1))then
              nstates=nstates+1
              if(((kk.gt.0).and.((kk/2)*2.eq.kk))
     &        .or.((kk.lt.0).and.((kk/2)*2.ne.kk)))then
                nste=nste+1
              else
                nsto=nsto+1
              endif
              number_states(2,kk)=jjj
            endif
            jjj=jjj+1
          enddo
        endif
      enddo

      deallocate(amat)
      deallocate(dmat)
      deallocate(cc)
      deallocate(wave)
      deallocate(b)
      deallocate(div2)
      deallocate(div1)
      deallocate(u)
      deallocate(sbj)
      if(dkb)then
        deallocate(b_2)
        deallocate(dd)
        deallocate(d3)
        deallocate(v_2)
        deallocate(v_22b)
        deallocate(v_22d)
        deallocate(v_22smalla)
        deallocate(v_22smallb)
        deallocate(v2sbj2)
        deallocate(v2sbj3)
        deallocate(v2sbj4)
        deallocate(dvd)
        deallocate(dv_1)
        deallocate(dvb)
      endif
      return
      end

      subroutine store_new_wave(nm,nstates,nkap,num_st,wave,eigvec,
     &wave_new)
      include 'inc.par'
      real*8 wave(2*nm,2*nm,-nkap:nkap)
      real*8 eigvec(nstates,nstates),
     &wave_new(nstates,2*nm,-nkap:nkap)
      integer num_st(-nkap:nkap,2*nm)

      wave_new=0.d0
      do num=1,nstates
        do kap=-nkap,nkap
          do mon=1,2*nm
            num_mon=num_st(kap,mon)
            if(num_mon.ne.0)then
              do nbs=1,2*nm
                wave_new(num,nbs,kap)=wave_new(num,nbs,kap)+
     &          eigvec(num_mon,num)*wave(nbs,mon,kap)
              enddo
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine store_new_wave_mj(nm,nstates,nkap,
     &wave_new,wave_new_mj)
      include 'inc.par'
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wave_new_mj(nstates,2*nm,-nkap:nkap,-nkap:nkap-1)

      wave_new_mj=0.d0
      do num=1,nstates
        do kap=-nkap,nkap
          if(kap.ne. 0)then
            do m_prj=-abs(kap),abs(kap)-1
              do nbs=1,2*nm
                wave_new_mj(num,nbs,kap,m_prj)=wave_new(num,nbs,kap)
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine store_new_wave_old(nm,nstates,nkap,num_st,wave,
     &eigvec,wave_new,rmin,rmax,ii_xi,xi_stepslower)
      include 'inc.par'
      real*8 wave(2*nm,2*nm,-nkap:nkap),eigvec(nstates,nstates),
     &wave_new(nstates,2*nm,-nkap:nkap)
      real*8, dimension(:),allocatable:: t
      real*8 ro(ns),ro1(ns),ro2(ns),WFValue_node1,
     &sampl_point,signn(nstates),
     &gradient(nstates)
      integer num_st(-nkap:nkap,2*nm)
      integer xi_stepslower
      logical dkb
      real*8, dimension(:), allocatable:: Line_value
      common /TargProj_prop/ Proj_mass,Targ_mass,Proj_vel,
     &b_ImpactParam
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /Barycentres/ RadiusOne,RadiusTwo
      common /dist/ distance,Starting_Distance
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max

c      Starting_R1=2.d0*Starting_Distance*Proj_mass/(Proj_mass+Targ_mass)
c      Starting_R2=2.d0*Starting_Distance*Targ_mass/(Proj_mass+Targ_mass)

      if(ii_xi.gt.xi_stepslower)then
        open(22,file='D_of_wf_at_small_r.dat')
        allocate(line_value(2))
        do i=1,nstates
c      read(22,'(A)',end=300)Line_readin
          read(22,*)Line_value(1),Line_value(2)
          signn(i)=Line_value(1)
          gradient(i)=Line_value(2)
        enddo
c 300  continue
        deallocate(line_value)

        close(22)
      endif

      nu=nm+ns+2
      allocate(t(nu))

      call spline_arranger(nu,rmin,rmax,t)

      wave_new=0.d0
      do num=1,nstates
        do kap=-nkap,nkap
          do mon=1,2*nm
            num_mon=num_st(kap,mon)
            if(num_mon.ne.0)then
              do nbs=1,2*nm
                wave_new(num,nbs,kap)=wave_new(num,nbs,kap)+
     &          eigvec(num_mon,num)*wave(nbs,mon,kap)
              enddo
            endif
          enddo
        enddo
      enddo

      sampl_point=1.01d0*rmin

      i2=1
      do while(t(i2).lt. sampl_point)
        i2=i2+1
      enddo
      if (t(i2).gt. sampl_point) then
        i2=i2-1
      endif
ccc i2=i2-1
      call ddsplines(sampl_point,ro,ro1,ro2,i2,t,nu)
c      Kapp=-nint(dabs(amu)+0.5d0)
c      write(*,*) 'nearest spline',i2, t(i2),sampl_point,i2-ns+1,Kapp
c      pause
      open(22,file='D_of_wf_at_small_r.dat')
      do i=1,nstates
        WFValue_node1=0.d0      !Value of G
        WFValue_node2=0.d0      !Value of DG
        WFValue_node3=0.d0      !Value of F
        WFValue_node4=0.d0      !Value of DF
c     MAYBE IT'S BETTER TO SAMPLE THE DERIVATIVE OF THE FUNCTION
c     NEAR 0.
        do kapp=-nkap,nkap
          if(kapp.ne.0)then
            do i_spl=max(2,i2-ns+1),min(i2,nm+1)
              if(dkb)then
                WFValue_node2=WFValue_node2+
     &          wave_new(i,i_spl-1,Kapp)*ro1(i_spl-i2+ns)+
     &          wave_new(i,i_spl+nm-1,Kapp)*0.5d0*(
     &          ro2(i_spl-i2+ns)-Kapp*ro1(i_spl-i2+ns)/sampl_point
     &          +Kapp*ro(i_spl-i2+ns)/sampl_point/sampl_point)
                WFValue_node1=WFValue_node1+
     &          wave_new(i,i_spl-1,Kapp)*ro(i_spl-i2+ns)+
     &          wave_new(i,i_spl+nm-1,Kapp)*0.5d0*(
     &          ro1(i_spl-i2+ns)-Kapp*ro(i_spl-i2+ns)/sampl_point)
                WFValue_node3=WFValue_node3+
     &          wave_new(i,i_spl+nm-1,Kapp)*ro(i_spl-i2+ns)+
     &          wave_new(i,i_spl-1,Kapp)*0.5d0*(
     &          ro1(i_spl-i2+ns)+Kapp*ro(i_spl-i2+ns)/sampl_point)
                WFValue_node4=WFValue_node4+
     &          wave_new(i,i_spl-1+nm,Kapp)*ro1(i_spl-i2+ns)+
     &          wave_new(i,i_spl-1,Kapp)*0.5d0*(
     &          ro2(i_spl-i2+ns)+Kapp*ro1(i_spl-i2+ns)/sampl_point-
     &          Kapp*ro(i_spl-i2+ns)/sampl_point/sampl_point)
              else
                WFValue_node2=WFValue_node2+
     &          wave_new(i,i_spl-1,Kapp)*ro1(i_spl-i2+ns)
                WFValue_node1=WFValue_node1+
     &          wave_new(i,i_spl-1,Kapp)*ro(i_spl-i2+ns)
                WFValue_node3=WFValue_node3+
     &          wave_new(i,i_spl-1+nm,Kapp)*ro(i_spl-i2+ns)
                WFValue_node4=WFValue_node4+
     &          wave_new(i,i_spl-1+nm,Kapp)*ro1(i_spl-i2+ns)
              endif
            enddo
          endif
        enddo
        DofFonG=WFValue_node4/WFValue_node1-
     &  WFValue_node3*WFValue_node2/WFValue_node1**2
c      write(*,*)i,DofFonG,WFValue_node1,WFValue_node2

        if(ii_xi.gt.xi_stepslower+3)then
          prediction=gradient(i)+signn(i)
          write(*,*)'PREDICTION',prediction,'ACTUAL',WFValue_node2,
     &    'PREVIOUS',signn(i)
          variance=dabs(WFValue_node2-prediction)
c            if(WFValue_node2*signn(i).lt. 0.d0)then
          if(variance.gt. dabs(gradient(i)*ii_xi))then
            do kap=-nkap,nkap
              if(kap.ne. 0) then
                do j=1,2*nm
                  wave_new(i,j,kap)=-wave_new(i,j,kap)
                enddo
              endif
            enddo
          endif
          WFValue_node2=-WFValue_node2
        endif

        if(ii_xi.eq.xi_stepslower)then
          delta_value=0.d0
        else
          delta_value=(WFValue_node2-signn(i))
        endif

        write(22,*)WFValue_node2,delta_value
        write(*,*)WFValue_node2,delta_value

      enddo
      close(22)
c      write(*,*) 'nearest spline',i2, t(i2),sampl_point,i2-ns+1,Kapp
c      pause
      deallocate(t)
c      pause

      return
      end

      subroutine ddsplines(x,ro,ro1,ro2,in,t1,nul)
      include 'inc.par'
      real*8 ro(ns),t1(nul),ro1(ns),ro2(ns)

      ro(ns)=1.d0
      ro1(ns)=0.d0
      ro2(ns)=0.d0
      do i=2,ns
         js=1

         ro2(ns-i+js)=-2.d0*ro1(ns-i+js+1)/(t1(in+js)-t1(in+js+1-i))+
     *        ro2(ns-i+js+1)*(t1(in+js)-x)/(t1(in+js)-t1(in+js+1-i))

         ro1(ns-i+js)=ro(ns-i+js+1)*(-1.d0)/(t1(in+js)-t1(in+js+1-i))+
     *        ro1(ns-i+js+1)*(t1(in+js)-x)/(t1(in+js)-t1(in+js+1-i))

         ro(ns-i+js)=ro(ns-i+js+1)*(t1(in+js)-x)/(t1(in+js)-
     *        t1(in+js+1-i))

         do js=2,i-1
            r1=2.d0*ro1(ns-i+js)/(t1(in+js-1)-t1(in+js-i))
     *           +ro2(ns-i+js)*(x-t1(in+js-i))/(t1(in+js-1)-t1(in+js-i))
            ro2(ns-i+js)=-2.d0*ro1(ns-i+js+1)/
     *           (t1(in+js)-t1(in+js+1-i))
     *           +ro2(ns-i+js+1)*(t1(in+js)-x)/
     *           (t1(in+js)-t1(in+js+1-i))
            ro2(ns-i+js)=ro2(ns-i+js)+r1

            r1=ro(ns-i+js)/(t1(in+js-1)-t1(in+js-i))
     *           +ro1(ns-i+js)*(x-t1(in+js-i))/(t1(in+js-1)-t1(in+js-i))
            ro1(ns-i+js)=-ro(ns-i+js+1)/(t1(in+js)-t1(in+js+1-i))
     *           +ro1(ns-i+js+1)*(t1(in+js)-x)/(t1(in+js)-t1(in+js+1-i))
            ro1(ns-i+js)=ro1(ns-i+js)+r1

            r1=ro(ns-i+js)*(x-t1(in+js-i))/(t1(in+js-1)-t1(in+js-i))
            ro(ns-i+js)=ro(ns-i+js+1)*(t1(in+js)-x)/
     *           (t1(in+js)-t1(in+js+1-i))
            ro(ns-i+js)=ro(ns-i+js)+r1
         enddo
         ro2(ns)=2.d0*ro1(ns)/(t1(in+i-1)-t1(in))
     *        +ro2(ns)*(x-t1(in))/(t1(in+i-1)-t1(in))

         ro1(ns)=ro(ns)/(t1(in+i-1)-t1(in))
     *        +ro1(ns)*(x-t1(in))/(t1(in+i-1)-t1(in))

         ro(ns)=ro(ns)*(x-t1(in))/(t1(in+i-1)-t1(in))
      enddo

      return
      end

      subroutine dsplines(x,ro,ro1,in,t1,nul)
      include 'inc.par'
      real*8 ro(ns),t1(nul),ro1(ns)

      ro(ns)=1.d0
      ro1(ns)=0.d0
      do i=2,ns
         js=1

         ro1(ns-i+js)=ro(ns-i+js+1)*(-1.d0)/(t1(in+js)-t1(in+js+1-i))+
     *        ro1(ns-i+js+1)*(t1(in+js)-x)/(t1(in+js)-t1(in+js+1-i))

         ro(ns-i+js)=ro(ns-i+js+1)*(t1(in+js)-x)/(t1(in+js)-
     *        t1(in+js+1-i))

         do js=2,i-1
            r1=ro(ns-i+js)/(t1(in+js-1)-t1(in+js-i))
     *           +ro1(ns-i+js)*(x-t1(in+js-i))/(t1(in+js-1)-t1(in+js-i))
            ro1(ns-i+js)=-ro(ns-i+js+1)/(t1(in+js)-t1(in+js+1-i))
     *           +ro1(ns-i+js+1)*(t1(in+js)-x)/(t1(in+js)-t1(in+js+1-i))
            ro1(ns-i+js)=ro1(ns-i+js)+r1

            r1=ro(ns-i+js)*(x-t1(in+js-i))/(t1(in+js-1)-t1(in+js-i))
            ro(ns-i+js)=ro(ns-i+js+1)*(t1(in+js)-x)/
     *           (t1(in+js)-t1(in+js+1-i))
            ro(ns-i+js)=ro(ns-i+js)+r1
         enddo
         ro1(ns)=ro(ns)/(t1(in+i-1)-t1(in))
     *        +ro1(ns)*(x-t1(in))/(t1(in+i-1)-t1(in))

         ro(ns)=ro(ns)*(x-t1(in))/(t1(in+i-1)-t1(in))
      enddo
      return
      end

c This secondary dsplines subroutine is to be used for calculating matrix
c elements of < n | dV/dR | k > as no such matrix element is already in
c use with the program Anton originally designed. The reason for this is
c the existence of different factors in the integral should r<R or r>R.

      subroutine ddsplines_for_dVdRMatElems(x,ro,ro1,ro2,in,t1,nul,R,
     & nkap,dv_dR)
      include 'inc.par'
      integer nkap
      real*8 ro(ns),t1(nul),ro1(ns),ro2(ns),
     &     dv_dR(0:2*nkap)
c      real*8, allocatable, dimension(:,:) :: ro2,ro12
c      open(1,file='nkap.inp',status='old')
c      read(1,*) nkap
c      close(1)
c      allocate(ro2(ns,2*nkap))
c      allocate(ro12(ns,2*nkap))
      call ddsplines(x,ro,ro1,ro2,in,t1,nul)
      do L=0,2*nkap
c$$$        if (t1(in+i).lt.R) then
c$$$        ro2(i,L)=(-L-1)*ro(i)/R
c        ro12(i,L)=(-L-1)*ro(i)/R
c$$$        elseif (t1(in+i+ns).gt.R) then
c$$$        ro2(i,L)=L*ro(i)/R
c        ro12(i,L)=L*ro1(i)/R
CCC**** Why this extra factor of (2L+1)/2 you ask? I think Anton has
CCC**** intermixed parts
CCC**** of the radial part of the potential with the angular parts, his
CCC**** radial part seems to contain an extra 2 out the front, which
CCC**** really belongs to the angular part, as far as I can see, hence I am
CCC**** neutralising that here.
         if (x.lt. R*0.5d0) then
            dv_dR(L)=(dble(-L-1))*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=(-L-1)*ro(i)/R
        else
           dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=L*ro(i)/R
        endif
      enddo
      return
      end

      subroutine dsplines_for_dVdRMatElems(x,ro,ro1,in,t1,nul,R,nkap
     &    ,dv_dR)
      include 'inc.par'
      integer nkap
      real*8 ro(ns),t1(nul),ro1(ns),
     &     dv_dR(0:2*nkap)
c      real*8, allocatable, dimension(:,:) :: ro2,ro12
c      open(1,file='nkap.inp',status='old')
c      read(1,*) nkap
c      close(1)
c      allocate(ro2(ns,2*nkap))
c      allocate(ro12(ns,2*nkap))
      call dsplines(x,ro,ro1,in,t1,nul)
      do L=0,2*nkap
         if (x.lt. R*0.5d0) then
            dv_dR(L)=(dble(-L-1))*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=(-L-1)*ro(i)/R
        else
           dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=L*ro(i)/R
        endif
      enddo
      return
      end

ccc******************** dV/dR Matrix elements for the baryonic expansion

      subroutine ddsplines_for_dVdRMatElems_Baryonic(
     &x,ro,ro1,ro2,in,t1,nul,R,nkap,dv_dR)
      include 'inc.par'
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /dist/distance,Starting_Distance
      common /Barycentres/ RadiusOne,RadiusTwo
      integer nkap
      real*8 u(0:2*nkap),ro(ns),t1(nul),ro1(ns),ro2(ns),
     &dv_dR(0:2*nkap)

      call ddsplines(x,ro,ro1,ro2,in,t1,nul)
      if(x.ge. RadiusOne .and. x.lt. RadiusTwo) then
      call all_potential_baryonic(nkap,x,u)
      endif

      do L=0,2*nkap
        if (x.lt. RadiusOne) then
          dv_dR(L)=(dble(-L-1))*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=(-L-1)*ro(i)/R
        elseif(x.ge. RadiusTwo) then
          dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=L*ro(i)/R
        else
ccc dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
          dv_dR(L)=-az1*dble(L)*dexp(dble(L)*dlog(RadiusOne))/
     &    dexp(dble(L+1)*dlog(x))/R-
     &    (-1)**L*az2*dble(-L-1)*dexp(dble(L)*dlog(x))/
     &    dexp(dble(L+1)*dlog(RadiusTwo))/R
          dv_dR(L)=dv_dR(L)/u(L)

c dv_dR(L)=-az1*dble(L)*RadiusOne**L/x**(L+1)/R-
c     &     (-1)**L*az2*dble(-L-1)*x**L/RadiusTwo**(L+1)/R
        endif
      enddo
      return
      end



      subroutine dsplines_for_dVdRMatElems_Baryonic(
     &x,ro,ro1,in,t1,nul,R,nkap,dv_dR)
      include 'inc.par'
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /dist/distance,Starting_Distance
      common /Barycentres/ RadiusOne,RadiusTwo
      integer nkap
      real*8 u(0:2*nkap),ro(ns),t1(nul),ro1(ns),!ro12(ns,0:2*nkap),
     &dv_dR(0:2*nkap)


c      real*8, allocatable, dimension(:,:) :: ro2,ro12
c      open(1,file='nkap.inp',status='old')
c      read(1,*) nkap
c      close(1)
c      allocate(ro2(ns,2*nkap))
c      allocate(ro12(ns,2*nkap))

      call dsplines(x,ro,ro1,in,t1,nul)

      if(x.ge. RadiusOne .and. x.lt. RadiusTwo) then
      call all_potential_baryonic(nkap,x,u)
      endif

      do L=0,2*nkap
        if (x.lt. RadiusOne) then
          dv_dR(L)=(dble(-L-1))*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=(-L-1)*ro(i)/R
        elseif(x.ge. RadiusTwo) then
          dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=L*ro(i)/R
        else
ccc dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
          dv_dR(L)=-az1*dble(L)*dexp(dble(L)*dlog(RadiusOne))/
     &    dexp(dble(L+1)*dlog(x))/R-
     &    (-1)**L*az2*dble(-L-1)*dexp(dble(L)*dlog(x))/
     &    dexp(dble(L+1)*dlog(RadiusTwo))/R

c dv_dR(L)=-az1*dble(L)*RadiusOne**L/x**(L+1)/R-
c     &     (-1)**L*az2*dble(-L-1)*x**L/RadiusTwo**(L+1)/R
          dv_dR(L)=dv_dR(L)/u(L)
        endif
      enddo

c      do L=0,2*nkap
c$$$        if (t1(in+i).lt.R) then
c$$$        ro2(i,L)=(-L-1)*ro(i)/R
c        ro12(i,L)=(-L-1)*ro(i)/R
c$$$        elseif (t1(in+i+ns).gt.R) then
c$$$        ro2(i,L)=L*ro(i)/R
c        ro12(i,L)=L*ro1(i)/R
CCC**** Why this extra factor of (2L+1)/2 you ask? I think Anton has
CCC**** intermixed parts
CCC**** of the radial part of the potential with the angular parts, his
CCC**** radial part seems to contain an extra 2 out the front, which
CCC**** really belongs to the angular part, as far as I can see, hence I am
CCC**** neutralising that here.
c         if (x.lt. RadiusOne) then
c       dv_dR(L)=(dble(-L-1))*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=(-L-1)*ro(i)/R
c      elseif(x.ge. RadiusTwo) then
c       dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
c        ro12(i,L)=L*ro(i)/R
c      else
ccc dv_dR(L)=dble(L)*dble(2*L+1)/(R*2.d0)
c dv_dR(L)=-az1*dble(L)*RadiusOne**(L-1)*
c     & (z_nuc2/(z_nuc1+z_nuc2))/x**(L+1)-
c     & (-1)**L*az2*x**L*dble((-L-1))*(z_nuc1/(z_nuc1+z_nuc2))/
c     &  RadiusTwo**(L+2)
c dv_dR(L)=dv_dR(L)/u(L)
c      endif
c      enddo
      return
      end


      SUBROUTINE D01BAZ(A, B, NPTS, WEIGHT, ABSCIS, IFAIL)
      DOUBLE PRECISION A, B
      INTEGER IFAIL, NPTS
      DOUBLE PRECISION ABSCIS(NPTS), WEIGHT(NPTS)
      DOUBLE PRECISION HFRNGE, PNTMID
      INTEGER I, IIJJ, N, NL, NN, NPTSA
      DOUBLE PRECISION ABST(136), WTST(136)
      INTEGER NSTOR(16)
      DATA WTST( 1),WTST( 2),WTST( 3),WTST( 4),WTST( 5)/+
     *0.200000000000000000000000000000D   1,+
     *0.100000000000000000000000000000D   1,+
     *0.555555555555555555555555555555D   0,+
     *0.888888888888888888888888888888D   0,+
     *0.347854845137453857373063949221D   0/
      DATA WTST( 6),WTST( 7),WTST(08),WTST(09),WTST(10)/+
     *0.652145154862546142626936050778D   0,+
     *0.236926885056189087514264040719D   0,+
     *0.478628670499366468041291514835D   0,+
     *0.568888888888888888888888888888D   0,+
     *0.171324492379170345040296142172D   0/
      DATA WTST(11),WTST(12),WTST(13),WTST(14),WTST(15)/+
     *0.360761573048138607569833513837D   0,+
     *0.467913934572691047389870343989D   0,+
     *0.101228536290376259152531354309D   0,+
     *0.222381034453374470544355994426D   0,+
     *0.313706645877887287337962201986D   0/
      DATA WTST(16),WTST(17),WTST(18),WTST(19),WTST(20)/+
     *0.362683783378361982965150449277D   0,+
     *0.666713443086881375935688098933D  -1,+
     *0.149451349150580593145776339657D   0,+
     *0.219086362515982043995534934228D   0,+
     *0.269266719309996355091226921569D   0/
      DATA WTST(21),WTST(22),WTST(23),WTST(24),WTST(25)/+
     *0.295524224714752870173892994651D   0,+
     *0.471753363865118271946159614850D  -1,+
     *0.106939325995318430960254718193D   0,+
     *0.160078328543346226334652529543D   0,+
     *0.203167426723065921749064455809D   0/
      DATA WTST(26),WTST(27),WTST(28),WTST(29),WTST(30)/+
     *0.233492536538354808760849898924D   0,+
     *0.249147045813402785000562436042D   0,+
     *0.351194603317518630318328761381D  -1,+
     *0.801580871597602098056332770628D  -1,+
     *0.121518570687903184689414809072D   0/
      DATA WTST(31),WTST(32),WTST(33),WTST(34),WTST(35)/+
     *0.157203167158193534569601938623D   0,+
     *0.185538397477937813741716590125D   0,+
     *0.205198463721295603965924065661D   0,+
     *0.215263853463157790195876443316D   0,+
     *0.271524594117540948517805724560D  -1/
      DATA WTST(36),WTST(37),WTST(38),WTST(39),WTST(40)/+
     *0.622535239386478928628438369943D  -1,+
     *0.951585116824927848099251076022D  -1,+
     *0.124628971255533872052476282192D   0,+
     *0.149595988816576732081501730547D   0,+
     *0.169156519395002538189312079030D   0/
      DATA WTST(41),WTST(42),WTST(43),WTST(44),WTST(45)/+
     *0.182603415044923588866763667969D   0,+
     *0.189450610455068496285396723208D   0,+
     *0.176140071391521183118619623518D  -1,+
     *0.406014298003869413310399522749D  -1,+
     *0.626720483341090635695065351870D  -1/
      DATA WTST(46),WTST(47),WTST(48),WTST(49),WTST(50)/+
     *0.832767415767047487247581432220D  -1,+
     *0.101930119817240435036750135480D   0,+
     *0.118194531961518417312377377711D   0,+
     *0.131688638449176626898494499748D   0,+
     *0.142096109318382051329298325067D   0/
      DATA WTST(51),WTST(52),WTST(53),WTST(54),WTST(55)/+
     *0.149172986472603746787828737001D   0,+
     *0.152753387130725850698084331955D   0,+
     *0.123412297999871995468056670700D  -1,+
     *0.285313886289336631813078159518D  -1,+
     *0.442774388174198061686027482113D  -1/
      DATA WTST(56),WTST(57),WTST(58),WTST(59),WTST(60)/+
     *0.592985849154367807463677585001D  -1,+
     *0.733464814110803057340336152531D  -1,+
     *0.861901615319532759171852029837D  -1,+
     *0.976186521041138882698806644642D  -1,+
     *0.107444270115965634782577342446D   0/
      DATA WTST(61),WTST(62),WTST(63),WTST(64),WTST(65)/+
     *0.115505668053725601353344483906D   0,+
     *0.121670472927803391204463153476D   0,+
     *0.125837456346828296121375382511D   0,+
     *0.127938195346752156974056165224D   0,+
     *0.701861000947009660040706373885D  -2/
      DATA WTST(66),WTST(67),WTST(68),WTST(69),WTST(70)/+
     *0.162743947309056706051705622063D  -1,+
     *0.253920653092620594557525897892D  -1,+
     *0.342738629130214331026877322523D  -1,+
     *0.428358980222266806568786466061D  -1,+
     *0.509980592623761761961632446895D  -1/
      DATA WTST( 71),WTST( 72),WTST( 73),WTST( 74),WTST( 75)/+
     *0.586840934785355471452836373001D  -1,+
     *0.658222227763618468376500637069D  -1,+
     *0.723457941088485062253993564784D  -1,+
     *0.781938957870703064717409188283D  -1,+
     *0.833119242269467552221990746043D  -1/
      DATA WTST( 76),WTST( 77),WTST( 78),WTST( 79),WTST( 80)/+
     *0.876520930044038111427714627518D  -1,+
     *0.911738786957638847128685771116D  -1,+
     *0.938443990808045656391802376681D  -1,+
     *0.956387200792748594190820022041D  -1,+
     *0.965400885147278005667648300635D  -1/
      DATA WTST( 81),WTST( 82),WTST( 83),WTST( 84),WTST( 85)/+
     *0.315334605230583863267731154389D  -2,+
     *0.732755390127626210238397962178D  -2,+
     *0.114772345792345394895926676090D  -1,+
     *0.155793157229438487281769558344D  -1,+
     *0.196161604573555278144607196522D  -1/
      DATA WTST( 86),WTST( 87),WTST( 88),WTST( 89),WTST( 90)/+
     *0.235707608393243791405193013784D  -1,+
     *0.274265097083569482000738362625D  -1,+
     *0.311672278327980889020657568463D  -1,+
     *0.347772225647704388925485859638D  -1,+
     *0.382413510658307063172172565237D  -1/
      DATA WTST( 91),WTST( 92),WTST( 93),WTST( 94),WTST( 95)/+
     *0.415450829434647492140588223610D  -1,+
     *0.446745608566942804194485871258D  -1,+
     *0.476166584924904748259066234789D  -1,+
     *0.503590355538544749578076190878D  -1,+
     *0.528901894851936670955050562646D  -1/
      DATA WTST( 96),WTST( 97),WTST( 98),WTST( 99),WTST(100)/+
     *0.551995036999841628682034951916D  -1,+
     *0.572772921004032157051502346847D  -1,+
     *0.591148396983956357464748174335D  -1,+
     *0.607044391658938800529692320278D  -1,+
     *0.620394231598926639041977841375D  -1/
      DATA WTST(101),WTST(102),WTST(103),WTST(104),WTST(105)/+
     *0.631141922862540256571260227502D  -1,+
     *0.639242385846481866239062018255D  -1,+
     *0.644661644359500822065041936577D  -1,+
     *0.647376968126839225030249387365D  -1,+
     *0.178328072169643294729607914497D  -2/
      DATA WTST(106),WTST(107),WTST(108),WTST(109),WTST(110)/+
     *0.414703326056246763528753572855D  -2,+
     *0.650445796897836285611736039998D  -2,+
     *0.884675982636394772303091465973D  -2,+
     *0.111681394601311288185904930192D  -1,+
     *0.134630478967186425980607666859D  -1/
      DATA WTST(111),WTST(112),WTST(113),WTST(114),WTST(115)/+
     *0.157260304760247193219659952975D  -1,+
     *0.179517157756973430850453020011D  -1,+
     *0.201348231535302093723403167285D  -1,+
     *0.222701738083832541592983303841D  -1,+
     *0.243527025687108733381775504090D  -1/
      DATA WTST(116),WTST(117),WTST(118),WTST(119),WTST(120)/+
     *0.263774697150546586716917926252D  -1,+
     *0.283396726142594832275113052002D  -1,+
     *0.302346570724024788679740598195D  -1,+
     *0.320579283548515535854675043478D  -1,+
     *0.338051618371416093915654821107D  -1/
      DATA WTST(121),WTST(122),WTST(123),WTST(124),WTST(125)/+
     *0.354722132568823838106931467152D  -1,+
     *0.370551285402400460404151018095D  -1,+
     *0.385501531786156291289624969468D  -1,+
     *0.399537411327203413866569261283D  -1,+
     *0.412625632426235286101562974736D  -1/
      DATA WTST(126),WTST(127),WTST(128),WTST(129),WTST(130)/+
     *0.424735151236535890073397679088D  -1,+
     *0.435837245293234533768278609737D  -1,+
     *0.445905581637565630601347100309D  -1,+
     *0.454916279274181444797709969712D  -1,+
     *0.462847965813144172959532492322D  -1/
      DATA WTST(131),WTST(132),WTST(133),WTST(134),WTST(135)/+
     *0.469681828162100173253262857545D  -1,+
     *0.475401657148303086622822069442D  -1,+
     *0.479993885964583077281261798713D  -1,+
     *0.483447622348029571697695271580D  -1,+
     *0.485754674415034269347990667839D  -1/
      DATA WTST(136)/+0.486909570091397203833653907347D  -1/
      DATA ABST( 1),ABST( 2),ABST( 3),ABST( 4),ABST( 5)/+
     *0.000000000000000000000000000000D   0,+
     *0.577350269189625764509148780501D   0,+
     *0.774596669241483377035853079956D   0,+
     *0.000000000000000000000000000000D   0,+
     *0.861136311594052575223946488892D   0/
      DATA ABST( 6),ABST( 7),ABST( 8),ABST( 9),ABST(10)/+
     *0.339981043584856264802665759103D   0,+
     *0.906179845938663992797626878299D   0,+
     *0.538469310105683091036314420700D   0,+
     *0.000000000000000000000000000000D   0,+
     *0.932469514203152027812301554493D   0/
      DATA ABST(11),ABST(12),ABST(13),ABST(14),ABST(15)/+
     *0.661209386466264513661399595019D   0,+
     *0.238619186083196908630501721680D   0,+
     *0.960289856497536231683560868569D   0,+
     *0.796666477413626739591553936475D   0,+
     *0.525532409916328985817739049189D   0/
      DATA ABST(16),ABST(17),ABST(18),ABST(19),ABST(20)/+
     *0.183434642495649804939476142360D   0,+
     *0.973906528517171720077964012084D   0,+
     *0.865063366688984510732096688423D   0,+
     *0.679409568299024406234327365114D   0,+
     *0.433395394129247190799265943165D   0/
      DATA ABST(21),ABST(22),ABST(23),ABST(24),ABST(25)/+
     *0.148874338981631210884826001129D   0,+
     *0.981560634246719250690549090149D   0,+
     *0.904117256370474856678465866119D   0,+
     *0.769902674194304687036893833212D   0,+
     *0.587317954286617447296702418940D   0/
      DATA ABST(26),ABST(27),ABST(28),ABST(29),ABST(30)/+
     *0.367831498998180193752691536643D   0,+
     *0.125233408511468915472441369463D   0,+
     *0.986283808696812338841597266704D   0,+
     *0.928434883663573517336391139377D   0,+
     *0.827201315069764993189794742650D   0/
      DATA ABST(31),ABST(32),ABST(33),ABST(34),ABST(35)/+
     *0.687292904811685470148019803019D   0,+
     *0.515248636358154091965290718551D   0,+
     *0.319112368927889760435671824168D   0,+
     *0.108054948707343662066244650219D   0,+
     *0.989400934991649932596154173450D   0/
      DATA ABST(36),ABST(37),ABST(38),ABST(39),ABST(40)/+
     *0.944575023073232576077988415534D   0,+
     *0.865631202387831743880467897712D   0,+
     *0.755404408355003033895101194847D   0,+
     *0.617876244402643748446671764048D   0,+
     *0.458016777657227386342419442983D   0/
      DATA ABST(41),ABST(42),ABST(43),ABST(44),ABST(45)/+
     *0.281603550779258913230460501460D   0,+
     *0.950125098376374401853193354249D  -1,+
     *0.993128599185094924786122388471D   0,+
     *0.963971927277913791267666131197D   0,+
     *0.912234428251325905867752441203D   0/
      DATA ABST(46),ABST(47),ABST(48),ABST(49),ABST(50)/+
     *0.839116971822218823394529061701D   0,+
     *0.746331906460150792614305070355D   0,+
     *0.636053680726515025452836696226D   0,+
     *0.510867001950827098004364050955D   0,+
     *0.373706088715419560672548177024D   0/
      DATA ABST(51),ABST(52),ABST(53),ABST(54),ABST(55)/+
     *0.227785851141645078080496195368D   0,+
     *0.765265211334973337546404093988D  -1,+
     *0.995187219997021360179997409700D   0,+
     *0.974728555971309498198391993008D   0,+
     *0.938274552002732758523649001708D   0/
      DATA ABST(56),ABST(57),ABST(58),ABST(59),ABST(60)/+
     *0.886415527004401034213154341982D   0,+
     *0.820001985973902921953949872669D   0,+
     *0.740124191578554364243828103099D   0,+
     *0.648093651936975569252495786910D   0,+
     *0.545421471388839535658375617218D   0/
      DATA ABST(61),ABST(62),ABST(63),ABST(64),ABST(65)/+
     *0.433793507626045138487084231913D   0,+
     *0.315042679696163374386793291319D   0,+
     *0.191118867473616309158639820757D   0,+
     *0.640568928626056260850430826247D  -1,+
     *0.997263861849481563544981128665D   0/
      DATA ABST(66),ABST(67),ABST(68),ABST(69),ABST(70)/+
     *0.985611511545268335400175044630D   0,+
     *0.964762255587506430773811928118D   0,+
     *0.934906075937739689170919134835D   0,+
     *0.896321155766052123965307243719D   0,+
     *0.849367613732569970133693004967D   0/
      DATA ABST( 71),ABST( 72),ABST( 73),ABST( 74),ABST( 75)/+
     *0.794483795967942406963097298970D   0,+
     *0.732182118740289680387426665091D   0,+
     *0.663044266930215200975115168663D   0,+
     *0.587715757240762329040745476401D   0,+
     *0.506899908932229390023747474377D   0/
      DATA ABST( 76),ABST( 77),ABST( 78),ABST( 79),ABST( 80)/+
     *0.421351276130635345364119436172D   0,+
     *0.331868602282127649779916805730D   0,+
     *0.239287362252137074544603209165D   0,+
     *0.144471961582796493485186373598D   0,+
     *0.483076656877383162348125704405D  -1/
      DATA ABST( 81),ABST( 82),ABST( 83),ABST( 84),ABST( 85)/+
     *0.998771007252426118600541491563D   0,+
     *0.993530172266350757547928750849D   0,+
     *0.984124583722826857744583600026D   0,+
     *0.970591592546247250461411983800D   0,+
     *0.952987703160430860722960666025D   0/
      DATA ABST( 86),ABST( 87),ABST( 88),ABST( 89),ABST( 90)/+
     *0.931386690706554333114174380101D   0,+
     *0.905879136715569672822074835671D   0,+
     *0.876572020274247885905693554805D   0,+
     *0.843588261624393530711089844519D   0,+
     *0.807066204029442627082553043024D   0/
      DATA ABST( 91),ABST( 92),ABST( 93),ABST( 94),ABST( 95)/+
     *0.767159032515740339253855437522D   0,+
     *0.724034130923814654674482233493D   0,+
     *0.677872379632663905211851280675D   0,+
     *0.628867396776513623995164933069D   0,+
     *0.577224726083972703817809238540D   0/
      DATA ABST( 96),ABST( 97),ABST( 98),ABST( 99),ABST(100)/+
     *0.523160974722233033678225869137D   0,+
     *0.466902904750958404544928861650D   0,+
     *0.408686481990716729916225495814D   0,+
     *0.348755886292160738159817937270D   0,+
     *0.287362487355455576735886461316D   0/
      DATA ABST(101),ABST(102),ABST(103),ABST(104),ABST(105)/+
     *0.224763790394689061224865440174D   0,+
     *0.161222356068891718056437390783D   0,+
     *0.970046992094626989300539558536D  -1,+
     *0.323801709628693620333222431521D  -1,+
     *0.999305041735772139456905624345D   0/
      DATA ABST(106),ABST(107),ABST(108),ABST(109),ABST(110)/+
     *0.996340116771955279346924500676D   0,+
     *0.991013371476744320739382383443D   0,+
     *0.983336253884625956931299302156D   0,+
     *0.973326827789910963741853507352D   0,+
     *0.961008799652053718918614121897D   0/
      DATA ABST(111),ABST(112),ABST(113),ABST(114),ABST(115)/+
     *0.946411374858402816062481491347D   0,+
     *0.929569172131939575821490154559D   0,+
     *0.910522137078502805756380668008D   0,+
     *0.889315445995114105853404038272D   0,+
     *0.865999398154092819760783385070D   0/
      DATA ABST(116),ABST(117),ABST(118),ABST(119),ABST(120)/+
     *0.840629296252580362751691544695D   0,+
     *0.813265315122797559741923338086D   0,+
     *0.783972358943341407610220525213D   0,+
     *0.752819907260531896611863774885D   0,+
     *0.719881850171610826848940217831D   0/
      DATA ABST(121),ABST(122),ABST(123),ABST(124),ABST(125)/+
     *0.685236313054233242563558371031D   0,+
     *0.648965471254657339857761231993D   0,+
     *0.611155355172393250248852971018D   0,+
     *0.571895646202634034283878116659D   0,+
     *0.531279464019894545658013903544D   0/
      DATA ABST(126),ABST(127),ABST(128),ABST(129),ABST(130)/+
     *0.489403145707052957478526307021D   0,+
     *0.446366017253464087984947714758D   0,+
     *0.402270157963991603695766771260D   0,+
     *0.357220158337668115950442615046D   0,+
     *0.311322871990210956157512698560D   0/
      DATA ABST(131),ABST(132),ABST(133),ABST(134),ABST(135)/+
     *0.264687162208767416373964172510D   0,+
     *0.217423643740007084149648748988D   0,+
     *0.169644420423992818037313629748D   0,+
     *0.121462819296120554470376463492D   0,+
     *0.729931217877990394495429419403D  -1/
      DATA ABST(136)/+0.243502926634244325089558428537D  -1/
      DATA NSTOR(1), NSTOR(2), NSTOR(3), NSTOR(4) /1,2,3,4/
      DATA NSTOR(5), NSTOR(6), NSTOR(7), NSTOR(8) /5,6,8,10/
      DATA NSTOR(9), NSTOR(10), NSTOR(11), NSTOR(12) /12,14,16,20/
      DATA NSTOR(13), NSTOR(14), NSTOR(15), NSTOR(16) /24,32,48,64/
      DO 20 I=1,NPTS
         WEIGHT(I) = 0.0D0
         ABSCIS(I) = 0.0D0
   20 CONTINUE
      N = 0
      NPTSA = 0
      IFAIL = 0
      DO 60 I=1,16
         IF (NPTS.LT.NSTOR(I)) GO TO 80
         N = N + (NPTSA+1)/2
         NPTSA = NSTOR(I)
         IF (NPTS.EQ.NSTOR(I)) GO TO 100
   60 CONTINUE
   80 IFAIL = 1
  100 HFRNGE = 0.5D0*(B-A)
      PNTMID = 0.5D0*(A+B)
      NL = NPTSA/2
      IF (NL.LT.1) GO TO 140
      DO 120 NN=1,NL
         N = N + 1
         IIJJ = NPTSA + 1 - NN
         ABSCIS(NN) = HFRNGE*ABST(N) + PNTMID
         WEIGHT(NN) = HFRNGE*WTST(N)
         ABSCIS(IIJJ) = -HFRNGE*ABST(N) + PNTMID
         WEIGHT(IIJJ) = HFRNGE*WTST(N)
  120 CONTINUE
  140 IF (NPTSA.LE.(NL+NL)) GO TO 160
      N = N + 1
      ABSCIS(NL+1) = HFRNGE*ABST(N) + PNTMID
      WEIGHT(NL+1) = HFRNGE*WTST(N)
  160 RETURN
      END

      subroutine solve_equation(nm,b,a,energy,wave)
      include 'inc.par'
      real*8 b(2*nm,2*nm),a(2*nm,2*nm),energy(2*nm),wave(2*nm,2*nm)
      real*8, dimension(:,:),allocatable:: vspom4,vsp4,vsp5
      real*8, dimension(:),allocatable:: vsp_blas2
      integer, dimension(:),allocatable:: IFAIL,vsp_blas3
      character*1 JOBZ,RANGE,UPLO,RR
      character*13 inp_dir
      common /blas_inp/ IL,IU,RR
      common /input_directory/ inp_dir


      allocate(vspom4(2*nm,2*nm))
      vspom4=0.d0
      allocate(vsp4(2*nm,2*nm))
      vsp4=0.d0
      allocate(vsp5(2*nm,2*nm))
      vsp5=0.d0



      allocate(IFAIL(2*nm))
      allocate(vsp_blas3(10*nm))
      allocate(vsp_blas2(16*nm))

      open(1,file=inp_dir//'blas.inp')
      read(1,*) JOBZ
      read(1,*) RANGE
      read(1,*) UPLO
      read(1,*) VL
      read(1,*) VU
      read(1,*) IL
      read(1,*) IU
      read(1,*) ABSTOL
      close(1)

      RR=RANGE
      wave=0.d0
      energy=0.d0


      do i=1,2*nm
         do j=1,2*nm
            vspom4(i,j)=a(i,j)
            vsp4(i,j)=b(i,j)
         enddo
      enddo

      NUM=0

      lwork1=-1
      call DSYGVX(1,JOBZ,RANGE,UPLO,2*nm,vspom4,2*nm,vsp4,2*nm,
     $     VL,VU,2*nm/2+1-IL,2*nm/2+1+IU,ABSTOL,NUM,energy,
     $     vsp5,2*nm,vsp_blas2,
     $     lwork1, vsp_blas3,IFAIL,INFO )

c vsp5 is the eigenvector associated with an energy eigenvalue.
c column 1 of vsp5 for example is the e-vector of the 1st e-val.
      lwork1=idint(vsp_blas2(1))

      deallocate(vsp_blas2)
      allocate(vsp_blas2(lwork1))


      call DSYGVX(1,JOBZ,RANGE,UPLO,2*nm,vspom4,2*nm,vsp4,2*nm,
     $     VL,VU,2*nm/2+1-IL,2*nm/2+1+IU,ABSTOL,NUM,energy,
     $     vsp5,2*nm,vsp_blas2,
     $     lwork1, vsp_blas3,IFAIL,INFO )

      if(info.ne.0) then
cc         write(*,*) INFO-2*nm,2*nm, 'ERROR IF NOT EQUAL TO ZERO!'
cc         write(*,*) b(INFO-2*nm,INFO-2*nm)
      else
c         write(*,*) 'NORMAL EXIT FROM DSYGVX',INFO
      endif
      do i=1,2*nm
         do j=1,NUM
            wave(i,j)=vsp5(i,j)
         enddo
      enddo

      deallocate(vspom4)
      deallocate(vsp4)
      deallocate(vsp5)
      deallocate(vsp_blas2)
      deallocate(IFAIL)
      deallocate(vsp_blas3)
      return
      end


      subroutine form_matrix(nm,nkap,nstates,number_states,num_st,
     &amat,wave,vmat,nvmat,e,ang,Translation_factor)
      include 'inc.par'
      integer number_states(2,-nkap:nkap)
      real*8 amat(nstates,nstates)
      real*8 wave(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),e(2*nm,-nkap:nkap),
     &ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      integer num_st(-nkap:nkap,2*nm)
      real*8, dimension(:,:,:),allocatable:: vmatr

c      write(*,*) 'COMPOSING CI MATRIX'

      num_st=0

      amat=0.d0

      n1=0
      do k1=-nkap,-1
        do j1=number_states(1,k1),number_states(2,k1)
          n1=n1+1
          amat(n1,n1)=e(j1,k1)
          num_st(k1,j1)=n1
        enddo
      enddo

      do k1=1,nkap
        do j1=number_states(1,k1),number_states(2,k1)
c            if(j1.ne.nm+1)then
          n1=n1+1
          amat(n1,n1)=e(j1,k1)
          num_st(k1,j1)=n1
        enddo
      enddo


      do k1=-nkap,nkap
        do k2=-nkap,nkap
c            write(*,*) 'COMPOSING KAPPAS',k1,k2
          nnn2=number_states(1,k2)
          nn2=number_states(2,k2)-number_states(1,k2)+1
          allocate(vmatr(nm,nn2,2))
          if(k1*k2.ne.0)then
            call stor_matr(k1,k2,
     &      number_states(1,k2),nn2,
     &      vmatr,nm,nkap,vmat,nvmat,wave,ang,1,1)
            do j1=number_states(1,k1),number_states(2,k1)
              do j2=number_states(1,k2),number_states(2,k2)
                n1=num_st(k1,j1)
                n2=num_st(k2,j2)
                if(n1*n2.ne.0)then
                  do i=1,nm
                    amat(n1,n2)=amat(n1,n2)+
     &              wave(i,j1,k1)*vmatr(i,j2-nnn2+1,1)
                    amat(n1,n2)=amat(n1,n2)+
     &              wave(i+nm,j1,k1)*vmatr(i,j2-nnn2+1,2)
                  enddo
                  if(n1.ne.n2)then
                    amat(n1,n2)=amat(n1,n2)*Translation_factor
!                        amat(n2,n1)=amat(n1,n2)
                  endif
c$$$                        write(*,*) 'STATES NUMBER',j1,j2
c$$$                        write(*,*) 'STATES KAPPAS',k1,k2
c$$$                        write(*,*) 'STATES ENERGIES',e(j1,k1),e(j2,k2)
c$$$                        write(*,*) 'MATRIX ELEMENT',n1,n2,amat(n1,n2)
c$$$   pause
c$$$                        if(amat(n1,n2).ne.0.d0) pause
                endif
              enddo
            enddo
          endif
          deallocate(vmatr)
        enddo
      enddo

c      write(*,*) 'MATRIX COMPOSED'

      return
      end

      subroutine stor_matr
     &  (k1,k2,n2m,n2g,v_matr,nm,nkap,vmat,nvmat,wave,ang,lstrt,lstp)
      include 'inc.par'
      real*8 vmat(nm,nm,0:2*nkap,nvmat),wave(2*nm,2*nm,-nkap:nkap),
     &ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      real*8 v_matr(nm,n2g,2)
      logical dkb
      common /common_dkb/dkb

      v_matr=0.d0

      if(dkb)then
         do i=1,nm
            do n2=n2m,n2g+n2m-1
               do j=max(1,i-ns+1),min(i+ns-1,nm)
                  do l=lstrt,2*nkap,lstp
                     if(ang(k1,k2,l).ne.0.d0)then
                        v_matr(i,n2-n2m+1,1)=v_matr(i,n2-n2m+1,1)+(
     &                       wave(j,n2,k2)*vmat(i,j,l,1)+
     &                       wave(j+nm,n2,k2)*(
     &                       (vmat(i,j,l,5)-vmat(i,j,l,7))/2.d0-
     &                       k2*vmat(i,j,l,8))
     &                       )*ang(k1,k2,l)
                        v_matr(i,n2-n2m+1,2)=v_matr(i,n2-n2m+1,2)+(
     &                       wave(j,n2,k2)*(
     &                       (vmat(i,j,l,5)+vmat(i,j,l,7))/2.d0-
     &                       k1*vmat(i,j,l,8)
     &                       )+
     &                       wave(j+nm,n2,k2)*(
     &                       vmat(i,j,l,2)+k1*k2*vmat(i,j,l,3)-
     &                       k1*(vmat(i,j,l,4)-vmat(i,j,l,6))/2.d0-
     &                       k2*(vmat(i,j,l,4)+vmat(i,j,l,6))/2.d0
     &                       )
     &                       )*ang(k1,k2,l)
                     endif
                     if(ang(-k1,-k2,l).ne.0.d0)then
                        v_matr(i,n2-n2m+1,2)=v_matr(i,n2-n2m+1,2)+(
     &                       wave(j+nm,n2,k2)*vmat(i,j,l,1)+
     &                       wave(j,n2,k2)*(
     &                       (vmat(i,j,l,5)-vmat(i,j,l,7))/2.d0+
     &                       k2*vmat(i,j,l,8))
     &                       )*ang(-k1,-k2,l)
                        v_matr(i,n2-n2m+1,1)=v_matr(i,n2-n2m+1,1)+(
     &                       wave(j+nm,n2,k2)*(
     &                       (vmat(i,j,l,5)+vmat(i,j,l,7))/2.d0+
     &                       k1*vmat(i,j,l,8)
     &                       )+
     &                       wave(j,n2,k2)*(
     &                       vmat(i,j,l,2)+k1*k2*vmat(i,j,l,3)+
     &                       k1*(vmat(i,j,l,4)-vmat(i,j,l,6))/2.d0+
     &                       k2*(vmat(i,j,l,4)+vmat(i,j,l,6))/2.d0
     &                       )
     &                       )*ang(-k1,-k2,l)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else
         do i=1,nm
            do n2=n2m,n2g+n2m-1
               do j=max(1,i-ns+1),min(i+ns-1,nm)
                  do l=lstrt,2*nkap,lstp
                     if(ang(k1,k2,l).ne.0.d0)then
                        v_matr(i,n2-n2m+1,1)=v_matr(i,n2-n2m+1,1)+
     &                  wave(j,n2,k2)*vmat(i,j,l,1)*ang(k1,k2,l)
                     endif
                     if(ang(-k1,-k2,l).ne.0.d0)then
                        v_matr(i,n2-n2m+1,2)=v_matr(i,n2-n2m+1,2)+
     &                  wave(j+nm,n2,k2)*vmat(i,j,l,1)*ang(-k1,-k2,l)
                     endif
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end

      subroutine form_matrix_even(nm,nkap,nstates,number_states,num_st,
     &amat,wave,vmat,nvmat,e,ang,Translation_factor)
      include 'inc.par'
      integer number_states(2,-nkap:nkap)
      real*8 amat(nstates,nstates)
      real*8 wave(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),e(2*nm,-nkap:nkap),
     &ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      integer num_st(-nkap:nkap,2*nm)
      real*8, dimension(:,:,:),allocatable:: vmatr

c      write(*,*) 'COMPOSING EVEN CI MATRIX'

      num_st=0

      amat=0.d0

      n1=0
      kfact=(-1)**nkap
      do k=nkap,1,-1
        k1=k*kfact
        kfact=-kfact
        do j1=number_states(1,k1),number_states(2,k1)
          n1=n1+1
          amat(n1,n1)=e(j1,k1)
          num_st(k1,j1)=n1
        enddo
      enddo

      kfact1=(-1)**nkap
      do k11=nkap,1,-1
        k1=k11*kfact1
        kfact1=-kfact1
        kfact2=(-1)**nkap
        do k22=nkap,1,-1
        k2=k22*kfact2
        kfact2=-kfact2
        nnn2=number_states(1,k2)
        nn2=number_states(2,k2)-number_states(1,k2)+1
        allocate(vmatr(nm,nn2,2))
        call stor_matr(k1,k2,nnn2,nn2,
     &  vmatr,nm,nkap,vmat,nvmat,wave,ang,2,2)
        do j1=number_states(1,k1),number_states(2,k1)
          do j2=number_states(1,k2),number_states(2,k2)
            n1=num_st(k1,j1)
            n2=num_st(k2,j2)
            if(n1*n2.ne.0)then
               do i=1,nm
                 amat(n1,n2)=amat(n1,n2)+
     &           wave(i,j1,k1)*vmatr(i,j2-nnn2+1,1)
                 amat(n1,n2)=amat(n1,n2)+
     &           wave(i+nm,j1,k1)*vmatr(i,j2-nnn2+1,2)
               enddo
               if(n1.ne.n2)then
                 amat(n1,n2)=amat(n1,n2)*Translation_factor
!                           amat(n2,n1)=amat(n1,n2)
               endif
             endif
           enddo
         enddo
         deallocate(vmatr)
        enddo
      enddo
c      write(*,*) 'MATRIX COMPOSED'

      return
      end

      subroutine
     &form_matrix_odd(nm,nkap,nstates,number_states,num_st,
     &amat,wave,vmat,nvmat,e,ang,Translation_factor)
      include 'inc.par'
      integer number_states(2,-nkap:nkap)
      real*8 amat(nstates,nstates)
      real*8 wave(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),e(2*nm,-nkap:nkap),
     &ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      integer num_st(-nkap:nkap,2*nm)
      real*8, dimension(:,:,:),allocatable:: vmatr

c      write(*,*) 'COMPOSING ODD CI MATRIX'

      num_st=0

      amat=0.d0

      n1=0
      kfact=(-1)**(nkap+1)
      do k=nkap,1,-1
      k1=k*kfact
      kfact=-kfact
        do j1=number_states(1,k1),number_states(2,k1)
          n1=n1+1
          amat(n1,n1)=e(j1,k1)
          num_st(k1,j1)=n1
        enddo
      enddo

      kfact1=(-1)**(nkap+1)
      do k11=nkap,1,-1
        k1=k11*kfact1
        kfact1=-kfact1
        kfact2=(-1)**(nkap+1)
        do k22=nkap,1,-1
          k2=k22*kfact2
          kfact2=-kfact2
          nnn2=number_states(1,k2)
          nn2=number_states(2,k2)-number_states(1,k2)+1
          allocate(vmatr(nm,nn2,2))
          call stor_matr(k1,k2,nnn2,nn2,
     &    vmatr,nm,nkap,vmat,nvmat,wave,ang,2,2)
            do j1=number_states(1,k1),number_states(2,k1)
              do j2=number_states(1,k2),number_states(2,k2)
                n1=num_st(k1,j1)
                n2=num_st(k2,j2)
                if(n1*n2.ne.0)then
                  do i=1,nm
                    amat(n1,n2)=amat(n1,n2)+
     &              wave(i,j1,k1)*vmatr(i,j2-nnn2+1,1)
                    amat(n1,n2)=amat(n1,n2)+
     &              wave(i+nm,j1,k1)*vmatr(i,j2-nnn2+1,2)
                  enddo
                if(n1.ne.n2)then
                  amat(n1,n2)=amat(n1,n2)*Translation_factor
!                           amat(n2,n1)=amat(n1,n2)
                endif
              endif
            enddo
          enddo
          deallocate(vmatr)
        enddo
      enddo
c      write(*,*) 'MATRIX COMPOSED'

      return
      end

      subroutine
     &form_matrix_2(nm,nkap,nstates,number_states,num_st,
     &amat,wave,vmat,nvmat,ang,nmj,Translation_factor,dtdxi)
      include 'inc.par'
      integer number_states(2,-nkap:nkap)
      real*8 amat(nstates,nstates,nmj)
      real*8 wave(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      integer num_st(-nkap:nkap,2*nm)
      real*8, dimension(:,:,:),allocatable:: vmatr

c      write(*,*) 'COMPOSING CI MATRIX'

      num_st=0

      amat=0.d0

      n1=0
      do k1=-nkap,nkap
        if(k1.ne.0)then
          do j1=number_states(1,k1),number_states(2,k1)
            n1=n1+1
c               if(abs(k1).ge.mj)then
c               amat(n1,n1,mj)=e(j1,k1)
c               endif
            num_st(k1,j1)=n1
          enddo
        endif
      enddo

      do mj=nmj,1,-1
        do k1=-nkap,nkap
          do k2=-nkap,nkap
c            write(*,*) 'COMPOSING KAPPAS',k1,k2
            nnn2=number_states(1,k2)
            nn2=number_states(2,k2)-number_states(1,k2)+1
            if(k1*k2.ne.0.and.abs(k1).ge.mj.and.abs(k2).ge.mj)then
              allocate(vmatr(nm,nn2,2))
              call stor_matr(k1,k2,number_states(1,k2),nn2,
     &        vmatr,nm,nkap,vmat,nvmat,wave,ang,1,1)
              do j1=number_states(1,k1),number_states(2,k1)
                do j2=number_states(1,k2),number_states(2,k2)
                  n1=num_st(k1,j1)
                  n2=num_st(k2,j2)
                  if(n1*n2.ne.0)then
                    do i=1,nm
                      amat(n1,n2,mj)=amat(n1,n2,mj)+
     &                wave(i,j1,k1)*vmatr(i,j2-nnn2+1,1)
                      amat(n1,n2,mj)=amat(n1,n2,mj)+
     &                wave(i+nm,j1,k1)*vmatr(i,j2-nnn2+1,2)
                    enddo
                    amat(n1,n2,mj)=amat(n1,n2,mj)*Translation_factor
                    amat(n2,n1,mj)=amat(n1,n2,mj)
c$$$                        write(*,*) 'STATES NUMBER',j1,j2
c$$$                        write(*,*) 'STATES KAPPAS',k1,k2
c$$$                        write(*,*) 'STATES ENERGIES',e(j1,k1),e(j2,k2)
c$$$                        write(*,*) 'MATRIX ELEMENT',n1,n2,amat(n1,n2)
c$$$   pause
c$$$                        if(amat(n1,n2).ne.0.d0) pause
                  endif
                enddo
              enddo
              deallocate(vmatr)
            endif
          enddo
        enddo
      enddo
      amat=amat*dtdxi

c      write(*,*) 'MATRIX COMPOSED'

      return
      end


      subroutine form_matrix_even_2(nm,nkap,nstates,number_states,
     &num_st,amat,wave,vmat,nvmat,ang,nmj,Translation_factor,
     &dtdxi)
      include 'inc.par'
      integer number_states(2,-nkap:nkap)
      real*8 amat(nstates,nstates,nmj)
      real*8 wave(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),
     &ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      integer num_st(-nkap:nkap,2*nm)
      real*8, dimension(:,:,:),allocatable:: vmatr

c      write(*,*) 'COMPOSING EVEN CI MATRIX'

      num_st=0

      amat=0.d0

      n1=0
      kfact=(-1)**nkap
      do k=nkap,1,-1
      k1=k*kfact
        do j1=number_states(1,k1),number_states(2,k1)
          n1=n1+1
c               if(abs(k1).ge.mj)then
c               amat(n1,n1,mj)=e(j1,k1)
c               endif
          num_st(k1,j1)=n1
        enddo
        kfact=-kfact
      enddo

c      do k1=nmin,-1,2
c         do k2=k1,nkap
      do mj=nmj,1,-1
        kfact1=(-1)**nkap
        do k=nkap,1,-1
          k1=k*kfact1
          kfact1=-kfact1
          kfact2=(-1)**nkap
            do kk=nkap,1,-1
            k2=kk*kfact2
            kfact2=-kfact2
            if(abs(k1).ge.mj.and.abs(k2).ge.mj)then
              nnn2=number_states(1,k2)
              nn2=number_states(2,k2)-number_states(1,k2)+1
              allocate(vmatr(nm,nn2,2))
              call stor_matr(k1,k2,number_states(1,k2),nn2,
     &        vmatr,nm,nkap,vmat,nvmat,wave,ang,2,2)
              do j1=number_states(1,k1),number_states(2,k1)
                do j2=number_states(1,k2),number_states(2,k2)
                  n1=num_st(k1,j1)
                  n2=num_st(k2,j2)
                  if(n1*n2.ne.0)then
                    do i=1,nm
                      amat(n1,n2,mj)=amat(n1,n2,mj)+
     &                wave(i,j1,k1)*vmatr(i,j2-nnn2+1,1)
                      amat(n1,n2,mj)=amat(n1,n2,mj)+
     &                wave(i+nm,j1,k1)*vmatr(i,j2-nnn2+1,2)
                    enddo
                    amat(n1,n2,mj)=amat(n1,n2,mj)*Translation_factor
!                        amat(n2,n1,mj)=amat(n1,n2,mj)

c                     write(*,*) 'STATES NUMBER',j1,j2
c                     write(*,*) 'STATES KAPPAS',k1,k2
c                     write(*,*) 'STATES ENERGIES',e(j1,k1),e(j2,k2)
c                     write(*,*) 'MATRIX ELEMENT',n1,n2,amat(n1,n2,mj)
c   pause
c                        if(amat(n1,n2,mj).ne.0.d0) pause
                  endif
                enddo
              enddo
              deallocate(vmatr)
            endif
          enddo
        enddo
      enddo
      amat=amat*dtdxi

      return
      end

      subroutine form_matrix_odd_2(nm,nkap,nstates,number_states,num_st,
     &amat,wave,vmat,nvmat,ang,nmj,Translation_factor,dtdxi)
      include 'inc.par'
      integer number_states(2,-nkap:nkap)
      real*8 amat(nstates,nstates,nmj)
      real*8 wave(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),
     &ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      integer num_st(-nkap:nkap,2*nm)
      real*8, dimension(:,:,:),allocatable:: vmatr

c      write(*,*) 'COMPOSING ODD CI MATRIX'

      num_st=0

      amat=0.d0

      n1=0
      kfact=(-1)**(nkap+1)
      do k=nkap,1,-1
        k1=k*kfact
        kfact=-kfact
        do j1=number_states(1,k1),number_states(2,k1)
          n1=n1+1
c               if(abs(k1).ge.mj)then
c               amat(n1,n1,mj)=e(j1,k1)
c               endif
          num_st(k1,j1)=n1
        enddo
      enddo

      do mj=nmj,1,-1
        kfact1=(-1)**(nkap+1)
        do k=nkap,1,-1
          k1=kfact1*k
          kfact1=-kfact1
          kfact2=(-1)**(nkap+1)
          do kk=nkap,1,-1
            k2=kk*kfact2
            kfact2=-kfact2
            if(abs(k1).ge.mj.and.abs(k2).ge.mj)then
              nnn2=number_states(1,k2)
              nn2=number_states(2,k2)-number_states(1,k2)+1
              allocate(vmatr(nm,nn2,2))
              call stor_matr(k1,k2,number_states(1,k2),nn2,
     &        vmatr,nm,nkap,vmat,nvmat,wave,ang,2,2)
              do j1=number_states(1,k1),number_states(2,k1)
                do j2=number_states(1,k2),number_states(2,k2)
                  n1=num_st(k1,j1)
                  n2=num_st(k2,j2)
                  if(n1*n2.ne.0)then
                    do i=1,nm
                      amat(n1,n2,mj)=amat(n1,n2,mj)+
     &                wave(i,j1,k1)*vmatr(i,j2-nnn2+1,1)
                      amat(n1,n2,mj)=amat(n1,n2,mj)+
     &                wave(i+nm,j1,k1)*vmatr(i,j2-nnn2+1,2)
                    enddo
                    amat(n1,n2,mj)=amat(n1,n2,mj)*Translation_factor
!                        amat(n2,n1,mj)=amat(n1,n2,mj)
                  endif
                enddo
              enddo
              deallocate(vmatr)
            endif
          enddo
        enddo
      enddo
      amat=amat*dtdxi

      return
      end

      subroutine b_spline_basis_only(nstates,nsto,nste,nm,nu,nkap,
     &number_states,rmin,rmax,wave1,vmat,nvmat,e,up_energy,dmat1,
     &alternate_dmat,distance)
      include 'inc.par'
      real*8, dimension(:),allocatable:: u,cc,t
      real*8, dimension(:,:),allocatable:: amat,wave,dmat,b,b_2,dd,d3
      real*8, dimension(:,:,:),allocatable:: div1,div2,v,v_2,dvd,
     &dv_1,dvb
      real*8 ro(ns),ro1(ns),ro2a(ns)
      real*8 wave1(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),e(2*nm,-nkap:nkap),
     &dmat1(2*nm,2*nm),alternate_dmat(2*nm,2*nm,-nkap:nkap),
     &ttt(nuz),www(nuz),el(2*nm)
      integer number_states(2,-nkap:nkap),One_s_state_location
      logical dkb,Manual_ncont_states
      common /weights/ w4n(4),t4n(4),w8(8),t8(8),w16(16),t16(16),
     &w32(32),t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /neg_cont/ Manual_ncont_states
      common /common_dkb/ dkb
      common /Barycentres/ RadiusOne,RadiusTwo
      common /TargProj_prop/ Proj_mass,Targ_mass,Proj_vel,
     &b_ImpactParam

      pi4=1.d0/(2.d0)
      charge_radius=Targ_mass**(1.d0/3.d0)*1.2d0*0.0025896063d0
      nul=nu
      alternate_dmat=0.d0
      dmat1=0.d0
      e=0.d0

      allocate(amat(2*nm,2*nm))
      allocate(dmat(2*nm,2*nm))
      allocate(cc(2*nm))
      allocate(wave(2*nm,2*nm))
      allocate(b(nm+2,nm+2))
      allocate(div2(nm+2,nm+2,-nkap:nkap))
      allocate(div1(nm+2,nm+2,-nkap:nkap))
      allocate(v(nm+2,nm+2,0:2*nkap))
      allocate(u(0:2*nkap))
      if(dkb)then
        allocate(b_2(nm+2,nm+2))
        allocate(dd(nm+2,nm+2))
        allocate(d3(nm+2,nm+2))
        allocate(v_2(nm+2,nm+2,0:2*nkap))
        allocate(dvd(nm+2,nm+2,0:2*nkap))
        allocate(dv_1(nm+2,nm+2,0:2*nkap))
        allocate(dvb(nm+2,nm+2,0:2*nkap))
      endif

      call getIntegralWeights(ttt, www)

!      select case (nuz)
!      case(4)
!        do i=1,4
!          ttt(i)=t4n(i)
!          www(i)=w4n(i)
!        enddo
!      case(8)
!        do i=1,8
!          ttt(i)=t8(i)
!          www(i)=w8(i)
!        enddo
!      case(6)
!        do i=1,6
!          ttt(i)=t6(i)
!          www(i)=w6(i)
!        enddo
!      case(16)
!        do i=1,16
!          ttt(i)=t16(i)
!          www(i)=w16(i)
!        enddo
!      case(32)
!        do i=1,32
!          ttt(i)=t32(i)
!          www(i)=w32(i)
!        enddo
!      case(64)
!        do i=1,64
!          ttt(i)=t64(i)
!          www(i)=w64(i)
!        enddo
!      end select

      v=0.d0
      b=0.d0
      div1=0.d0
      div2=0.d0
      amat=0.d0
      dmat=0.d0
      cc=0.d0
      vmat=0.d0
      if(dkb) then
        b_2=0.d0
        dd=0.d0
        d3=0.d0
        v_2=0.d0
        dvd=0.d0
        dv_1=0.d0
        dvb=0.d0
      endif
      allocate(t(nu))
      call spline_arranger(nu,rmin,rmax,t)

      istepforll=1
      if(z_nuc1.eq.z_nuc2)istepforll=2

      do inx1=ns,nm+2
        ac=(t(inx1)+t(inx1+1))/2.d0
        bc=(t(inx1+1)-t(inx1))/2.d0
        do j=nuz,1,-1
          xx=ac+bc*ttt(j)
          bbc=bc*www(j)
          if(dkb)then
            call ddsplines(xx,ro,ro1,ro2a,inx1,t,nu)
          else
            call dsplines(xx,ro,ro1,inx1,t,nu)
          endif
          if(z_nuc1.eq.z_nuc2)then
            call all_potential_dist(nkap,xx,u,distance)
          elseif(z_nuc1.ne.z_nuc2)then
            call all_potential_baryonic_dist(nkap,xx,u,distance)
          endif

          do i=inx1-ns+1,inx1
            do k=inx1-ns+1,inx1

              b(i,k)=b(i,k)+ro(i-inx1+ns)*ro(k-inx1+ns)*bbc
              if(dkb)then
                b_2(i,k)=b_2(i,k)+
     &          ro(i-inx1+ns)*ro(k-inx1+ns)*bbc/xx/xx
                dd(i,k)=dd(i,k)+
     &          ro1(i-inx1+ns)*ro1(k-inx1+ns)*bbc
                d3(i,k)=d3(i,k)+
     &          ro1(i-inx1+ns)*ro2a(k-inx1+ns)*bbc
              endif
              do ll=0,2*nkap,istepforll
                v(i,k,ll)=v(i,k,ll)+
     &          ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*u(ll)
                if(dkb)then
                  v_2(i,k,ll)=v_2(i,k,ll)+
     &            ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*u(ll)/xx/xx
                  dvd(i,k,ll)=dvd(i,k,ll)+
     &            ro1(i-inx1+ns)*ro1(k-inx1+ns)*bbc*u(ll)
                  dv_1(i,k,ll)=dv_1(i,k,ll)+
     &            (ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &            ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &            bbc*u(ll)/xx
                  if((i.gt.1).and.(k.gt.1).and.(i-1.le.nm).and.
     &            (k-1.le.nm))then
                    vmat(i-1,k-1,ll,6)=vmat(i-1,k-1,ll,6)+
     &              (-ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &              bbc*u(ll)/xx/4.d0
                    vmat(i-1,k-1,ll,7)=vmat(i-1,k-1,ll,7)+
     &              (-ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &              bbc*u(ll)/2.d0
                    vmat(i-1,k-1,ll,8)=vmat(i-1,k-1,ll,8)+
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              bbc*u(ll)/xx/2.d0

                  endif
                  dvb(i,k,ll)=dvb(i,k,ll)+
     &            (ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &            ro1(i-inx1+ns)*ro(k-inx1+ns))
     &            *bbc*u(ll)

                endif
              enddo
              do kk=-nkap,nkap
                if(kk.ne.0)then
                  if(dkb)then
                    div1(i,k,kk)=div1(i,k,kk)+(
     &              ro(i-inx1+ns)*ro2a(k-inx1+ns)*kk/xx-
     &              ro1(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*(kk-1)/xx/xx-
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*kk*(kk-1)/xx/xx/xx
     &              )*bbc
                    div2(i,k,kk)=div2(i,k,kk)+(
     &              ro(i-inx1+ns)*ro2a(k-inx1+ns)*kk/xx+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*(kk+1)/xx/xx-
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*kk*(kk+1)/xx/xx/xx
     &              )*bbc
                  else
                    div1(i,k,kk)=div1(i,k,kk)+(
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*kk/xx-
     &              ro(i-inx1+ns)*ro1(k-inx1+ns))*bbc
                    div2(i,k,kk)=div2(i,k,kk)+(
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*kk/xx+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*bbc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(t)
      if(dkb) then
        do i=1,nm
          do k=1,nm
            l=i+1
            m=k+1
            do ll=0,2*nkap,istepforll
              vmat(i,k,ll,1)=v(l,m,ll)
              vmat(i,k,ll,2)=dvd(l,m,ll)/4.d0
              vmat(i,k,ll,3)=v_2(l,m,ll)/4.d0
              vmat(i,k,ll,4)=dv_1(l,m,ll)/4.d0
              vmat(i,k,ll,5)=dvb(l,m,ll)/2.d0
            enddo
          enddo
        enddo
      else
        do i=1,nm
          do k=1,nm
            l=i+1
            m=k+1
            amat(i,k)=b(l,m)+v(l,m,0)*pi4

            amat(i+nm,k+nm)=-b(l,m)+v(l,m,0)*pi4

            dmat(i,k)=b(l,m)
            dmat(i+nm,k+nm)=b(l,m)
            do ll=0,2*nkap,istepforll
              vmat(i,k,ll,1)=v(l,m,ll)
            enddo
          enddo
        enddo
      endif
      deallocate(v)
      nstates=0
      nsto=0
      nste=0
      number_states=0
      do kk=-nkap,nkap
        if(kk.ne.0)then
          do i=1,nm
            do k=1,nm
              l=i+1
              m=k+1
              if(dkb) then
                amat(i,k)=b(l,m)+vmat(i,k,0,1)*pi4+7.5d-1*dd(l,m)+
     &          7.5d-1*kk*(kk+1)*b_2(l,m)+vmat(i,k,0,2)*pi4
     &          +kk*kk*vmat(i,k,0,3)*pi4+kk*vmat(i,k,0,4)*pi4

                amat(i+nm,k+nm)=-b(l,m)+vmat(i,k,0,1)*pi4
     &          -7.5d-1*dd(l,m)-
     &          7.5d-1*kk*(kk-1)*b_2(l,m)+vmat(i,k,0,2)*pi4
     &          +kk*kk*vmat(i,k,0,3)*pi4-kk*vmat(i,k,0,4)*pi4

                amat(i,k+nm)=vmat(i,k,0,5)*pi4+
     &          (d3(l,m)-d3(m,l))/8.d0+
     &          (div1(l,m,kk)+div2(m,l,kk))/8.d0

                amat(k+nm,i)=amat(i,k+nm)

                dmat(i,k)=b(l,m)+dd(l,m)/4.d0+kk*b_2(l,m)/
     &          4.d0*(kk+1)
                dmat(i+nm,k+nm)=b(l,m)+dd(l,m)/4.d0+
     &          kk*b_2(l,m)/4.d0*(kk-1)
              else
                amat(i,k+nm)=(div1(l,m,kk)+div2(l,m,kk))/2.d0
                amat(k+nm,i)=amat(i,k+nm)
              endif
            enddo
          enddo
          do i=1,2*nm
            do k=1,2*nm
              dmat1(i,k)=dmat(i,k)
              alternate_dmat(i,k,kk)=dmat(i,k)
            enddo
          enddo
          call solve_equation(nm,dmat,amat,cc,wave)
          if(cc(nm+1).gt. -1.d0) then
            do i=1,2*nm
              e(i,kk)=cc(i)
              do j=1,2*nm
                wave1(j,i,kk)=wave(j,i)
              enddo
            enddo
          elseif(cc(nm+1).lt. -1.d0) then
            nrange=nm+1
            el=0.d0
            do j=1,nm+1
              if (cc(j).gt. -2.d0 .and. cc(j).lt. -1.d0)then
                do l=1,nm/2+1
                  do m=1,nm/2+1
                    el(j)=el(j)+
     &              (wave(l,j)*wave(m,j)*dmat1(l,m)+wave(l+nm,j)*
     &              wave(m+nm,j)*dmat1(l+nm,m+nm))
                  enddo
                enddo
              endif
            enddo
            One_s_state_location=maxloc(el, dim=1, mask= el .lt. 1.d0)
            do i=1,2*nm
              if (i .lt. One_s_state_location) then
                e(i,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i,kk)=wave(j,i)
                enddo
              elseif(i .eq. One_s_state_location) then
                e(nm+1,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,nm+1,kk)=wave(j,i)
                enddo
              elseif(i .gt. One_s_state_location .and.
     &        i .le. nm+1) then
                e(i-1,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i-1,kk)=wave(j,i)
                enddo
              else
                e(i,kk)=cc(i)
                WFValue_node1=0.d0
                do j=1,2*nm
                  wave1(j,i,kk)=wave(j,i)
                enddo
              endif
            enddo
          endif
          jjj=1
          if(Manual_ncont_states)then
            do while(jjj.le. 2*nm .and. dabs(cc(jjj)).gt.up_energy)
              jjj=jjj+1
            enddo
            number_states(1,kk)=jjj
          else
            number_states(1,kk)=nm+1
            jjj=nm+1
          endif
          do while(jjj.le.2*nm .and. cc(jjj).le.up_energy)
            nstates=nstates+1
            if(((kk.gt.0).and.((kk/2)*2.eq.kk))
     &      .or.((kk.lt.0).and.((kk/2)*2.ne.kk)))then
              nste=nste+1
            else
              nsto=nsto+1
            endif
            number_states(2,kk)=jjj
            jjj=jjj+1
          enddo
        endif
      enddo

      deallocate(amat)
      deallocate(dmat)
      deallocate(cc)
      deallocate(wave)
      deallocate(b)
      deallocate(div2)
      deallocate(div1)
      deallocate(u)
      if(dkb)then
        deallocate(b_2)
        deallocate(dd)
        deallocate(d3)
        deallocate(v_2)
        deallocate(dvd)
        deallocate(dv_1)
        deallocate(dvb)
      endif
      return
      end

      subroutine b_spline_calculation_no_laser(nstates,nsto,nste,
     &nm,nu,nkap,number_states,rmin,rmax,wave1,vmat,nvmat,
     &e,up_energy,dmat1,dvdRmat_dkb1,dvdRmat_dkb2,alternate_dmat)
      include 'inc.par'
      real*8, dimension(:),allocatable:: u,cc,t
      real*8, dimension(:,:),allocatable:: amat,wave,dmat,b,b_2,dd,d3
      real*8, dimension(:,:,:),allocatable:: div1,div2,v,
     &v_2,dvd,dv_1,dvb,v2
      real*8, dimension(:,:,:,:),allocatable:: v_22b,v_22d
      real*8, dimension(:,:,:,:,:),allocatable::v_22smalla,v_22smallb
      real*8 ro(ns),ro1(ns),ro2a(ns)
      real*8 wave1(2*nm,2*nm,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat),e(2*nm,-nkap:nkap),
     &dmat1(2*nm,2*nm),dvdRmat_dkb1(nm,nm,-nkap:nkap,0:2*nkap,2),
     &dvdRmat_dkb2(nm,nm,-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &alternate_dmat(2*nm,2*nm,-nkap:nkap)
      common /weights/ w4n(4),t4n(4),w8(8),t8(8),w16(16),t16(16),
     &w32(32),t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 ttt(nuz),www(nuz),dv_dR(0:2*nkap),el(2*nm)
      common /dist/ distance,Starting_Distance
      integer number_states(2,-nkap:nkap)
      integer One_s_state_location,xi_stepslower,xi_stepsupper
      logical dkb,Manual_ncont_states
      common /neg_cont/ Manual_ncont_states
      common /common_dkb/ dkb
      common /Barycentres/ RadiusOne,RadiusTwo
      common /step_progress/ xi_stepslower,xi_stepsupper,ii_xi
      common /TargProj_prop/ Proj_mass,Targ_mass,Proj_vel,
     &b_ImpactParam

      pi4=1.d0/(2.d0)

      charge_radius=Targ_mass**(1.d0/3.d0)*1.2d0*0.0025896063d0
      Starting_R1=2.d0*Starting_Distance*Proj_mass/(Proj_mass+Targ_mass)
      Starting_R2=2.d0*Starting_Distance*Targ_mass/(Proj_mass+Targ_mass)

      nul=nu
      alternate_dmat=0.d0
      dmat1=0.d0
      dvdRmat_dkb1=0.d0
      dvdRmat_dkb2=0.d0
      e=0.d0

c b2,v2, are added for the
c calculation of the dV/dR matrix

      allocate(amat(2*nm,2*nm))
      allocate(dmat(2*nm,2*nm))
      allocate(cc(2*nm))
      allocate(wave(2*nm,2*nm))
      allocate(b(nm+2,nm+2))
c      allocate(b2(nm+2,nm+2,0:2*nkap))
      allocate(div2(nm+2,nm+2,-nkap:nkap))
      allocate(div1(nm+2,nm+2,-nkap:nkap))
      allocate(v(nm+2,nm+2,0:2*nkap))
      allocate(v2(nm+2,nm+2,0:2*nkap))
      allocate(u(0:2*nkap))
      if(dkb)then
        allocate(b_2(nm+2,nm+2))
        allocate(dd(nm+2,nm+2))
        allocate(d3(nm+2,nm+2))
        allocate(v_2(nm+2,nm+2,0:2*nkap))
        allocate(v_22b(-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(v_22d(-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(v_22smalla(-nkap:nkap,-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(v_22smallb(-nkap:nkap,-nkap:nkap,nm+2,nm+2,0:2*nkap))
        allocate(dvd(nm+2,nm+2,0:2*nkap))
        allocate(dv_1(nm+2,nm+2,0:2*nkap))
        allocate(dvb(nm+2,nm+2,0:2*nkap))
      endif

      call getIntegralWeights(ttt, www)

!      select case (nuz)
!      case(4)
!        do i=1,nuz
!          ttt(i)=t4n(i)
!          www(i)=w4n(i)
!        enddo
!      case(8)
!        do i=1,nuz
!          ttt(i)=t8(i)
!          www(i)=w8(i)
!        enddo
!      case(6)
!        do i=1,nuz
!          ttt(i)=t6(i)
!          www(i)=w6(i)
!        enddo
!      case(16)
!        do i=1,nuz
!          ttt(i)=t16(i)
!          www(i)=w16(i)
!        enddo
!      case(32)
!        do i=1,nuz
!          ttt(i)=t32(i)
!          www(i)=w32(i)
!        enddo
!      case(64)
!        do i=1,nuz
!          ttt(i)=t64(i)
!          www(i)=w64(i)
!        enddo
!      end select

      v=0.d0
      b=0.d0
      div1=0.d0
      div2=0.d0
      amat=0.d0
      dmat=0.d0
      cc=0.d0
      vmat=0.d0
      if(dkb) then
        b_2=0.d0
        dd=0.d0
        d3=0.d0
        v_2=0.d0
c      v_22=0.d0, v_22 is set to zero elsewhere.
        dvd=0.d0
        dv_1=0.d0
        dvb=0.d0
      endif
      allocate(t(nu))
      call spline_arranger(nu,rmin,rmax,t)

      if(z_nuc1.eq.z_nuc2)then
         istepforll=2
      else
         istepforll=1
      endif
c     The calculation of the matrix
      do inx1=ns,nm+2
        ac=(t(inx1)+t(inx1+1))/2.d0
        bc=(t(inx1+1)-t(inx1))/2.d0

        do j=nuz,1,-1
          xx=ac+bc*ttt(j)
          bbc=bc*www(j)

          if(dkb)then
            call ddsplines(xx,ro,ro1,ro2a,inx1,t,nu)
          else
            call dsplines(xx,ro,ro1,inx1,t,nu)
          endif
          if(z_nuc1.eq.z_nuc2)then
            call all_potential(nkap,xx,u)
          elseif(z_nuc1.ne.z_nuc2 )then
            call all_potential_baryonic(nkap,xx,u)
          endif

          do i=inx1-ns+1,inx1
            do k=inx1-ns+1,inx1

              b(i,k)=b(i,k)+ro(i-inx1+ns)*ro(k-inx1+ns)*bbc
              if(dkb)then
                b_2(i,k)=b_2(i,k)+
     &          ro(i-inx1+ns)*ro(k-inx1+ns)*bbc/xx/xx
                dd(i,k)=dd(i,k)+
     &          ro1(i-inx1+ns)*ro1(k-inx1+ns)*bbc
                d3(i,k)=d3(i,k)+
     &          ro1(i-inx1+ns)*ro2a(k-inx1+ns)*bbc
              endif
              do ll=0,2*nkap,istepforll
                v(i,k,ll)=v(i,k,ll)+
     &          ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*u(ll)

                if(dkb)then
                  v_2(i,k,ll)=v_2(i,k,ll)+
     &            ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*u(ll)/xx/xx
                  dvd(i,k,ll)=dvd(i,k,ll)+
     &            ro1(i-inx1+ns)*ro1(k-inx1+ns)*bbc*u(ll)
                  dv_1(i,k,ll)=dv_1(i,k,ll)+
     &            (ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &            ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &            bbc*u(ll)/xx
                  if((i.gt.1).and.(k.gt.1).and.(i-1.le.nm).and.
     &            (k-1.le.nm))then
                    vmat(i-1,k-1,ll,6)=vmat(i-1,k-1,ll,6)+
     &              (-ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &              bbc*u(ll)/xx/4.d0
                    vmat(i-1,k-1,ll,7)=vmat(i-1,k-1,ll,7)+
     &              (-ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &              bbc*u(ll)/2.d0
                    vmat(i-1,k-1,ll,8)=vmat(i-1,k-1,ll,8)+
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              bbc*u(ll)/xx/2.d0

                  endif
                  dvb(i,k,ll)=dvb(i,k,ll)+
     &            (ro(i-inx1+ns)*ro1(k-inx1+ns)+
     &            ro1(i-inx1+ns)*ro(k-inx1+ns))*
     &            bbc*u(ll)

                endif
              enddo
              do kk=-nkap,nkap
                if(kk.ne.0)then
                  if(dkb)then
                    div1(i,k,kk)=div1(i,k,kk)+(
     &              ro(i-inx1+ns)*ro2a(k-inx1+ns)*kk/xx-
     &              ro1(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*(kk-1)/xx/xx-
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*kk*(kk-1)/xx/xx/xx
     &              )*bbc
                    div2(i,k,kk)=div2(i,k,kk)+(
     &              ro(i-inx1+ns)*ro2a(k-inx1+ns)*kk/xx+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*(kk+1)/xx/xx-
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*
     &              kk*kk*(kk+1)/xx/xx/xx
     &              )*bbc
                  else
                    div1(i,k,kk)=div1(i,k,kk)+(
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*kk/xx-
     &              ro(i-inx1+ns)*ro1(k-inx1+ns))*bbc
                    div2(i,k,kk)=div2(i,k,kk)+(
     &              ro(i-inx1+ns)*ro(k-inx1+ns)*kk/xx+
     &              ro1(i-inx1+ns)*ro(k-inx1+ns))*bbc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      enddo


      if(dkb) then
        do i=1,nm
          do k=1,nm
            l=i+1
            m=k+1
            do ll=0,2*nkap,istepforll
              vmat(i,k,ll,1)=v(l,m,ll)
              vmat(i,k,ll,2)=dvd(l,m,ll)/4.d0
              vmat(i,k,ll,3)=v_2(l,m,ll)/4.d0
              vmat(i,k,ll,4)=dv_1(l,m,ll)/4.d0
              vmat(i,k,ll,5)=dvb(l,m,ll)/2.d0
            enddo
          enddo
        enddo
      else
        do i=1,nm
          do k=1,nm
            l=i+1
            m=k+1
            amat(i,k)=b(l,m)+v(l,m,0)*pi4

            amat(i+nm,k+nm)=-b(l,m)+v(l,m,0)*pi4

            dmat(i,k)=b(l,m)
            dmat(i+nm,k+nm)=b(l,m)
            do ll=0,2*nkap,istepforll
              vmat(i,k,ll,1)=v(l,m,ll)
            enddo
          enddo
        enddo
      endif

c****** MATRIX OF dV/dR
c This part is the calculation of my so called dV/dR matrix element,
c similar to the above, yet somewhat altered.

      v2=0.d0

      if(dkb) then
        v_22b=0.d0
        v_22d=0.d0
        v_22smalla=0.d0
        v_22smallb=0.d0
      endif

      R=2.d0*distance
      lstep=1
      if(z_nuc1.eq.z_nuc2)lstep=2

      do inx1=ns,nm+2
        ac=(t(inx1)+t(inx1+1))/2.d0
        bc=(t(inx1+1)-t(inx1))/2.d0
        do j=nuz,1,-1
          xx=ac+bc*ttt(j)
          bbc=bc*www(j)
          if(z_nuc1.eq.z_nuc2)then
            call dvdr_potential(nkap,xx,dv_dR)
          else
            call dvdr_baryonic(nkap,xx,dv_dR)
          endif
          call dsplines(xx,ro,ro1,inx1,t,nu)

          do i=inx1-ns+1,inx1
            do k=inx1-ns+1,inx1
              do ll=0,2*nkap,lstep
                v2(i,k,ll)=v2(i,k,ll)+ro(i-inx1+ns)*ro(k-inx1+ns)*bbc*
     &          dv_dR(ll)
              enddo
            enddo
          enddo
          if(dkb) then
            do i=inx1-ns+1,inx1
              do k=inx1-ns+1,inx1
                do k1=-nkap,nkap
                  if(k1.ne. 0) then
                    do ll=0,2*nkap,lstep
                      v_22b(k1,i,k,ll)=v_22b(k1,i,k,ll)+
     &                ro(k-inx1+ns)*(ro1(i-inx1+ns)-
     &                k1*ro(i-inx1+ns)/xx)*bbc*dv_dR(ll)/2.d0

                      v_22d(k1,i,k,ll)=v_22d(k1,i,k,ll)+
     &                ro(k-inx1+ns)*(ro1(i-inx1+ns)+
     &                k1*ro(i-inx1+ns)/xx)*bbc*dv_dR(ll)/2.d0
                    enddo
                  endif
                  do k22=-nkap,nkap
                    if(k1*k22.ne. 0) then
                      do ll=0,2*nkap,lstep
                        v_22smalla(k1,k22,i,k,ll)=
     &                  v_22smalla(k1,k22,i,k,ll)+
     &                  (ro1(i-inx1+ns)*ro1(k-inx1+ns)+
     &                  k1*k22*ro(i-inx1+ns)*ro(k-inx1+ns)/xx/xx-
     &                  k1*(ro(i-inx1+ns)*ro1(k-inx1+ns))/xx-
     &                  k22*(ro1(i-inx1+ns)*ro(k-inx1+ns))/xx)*
     &                  bbc*dv_dR(ll)/4.d0

                        v_22smallb(k1,k22,i,k,ll)=
     &                  v_22smallb(k1,k22,i,k,ll)+
     &                  (ro1(i-inx1+ns)*ro1(k-inx1+ns)+
     &                  k1*k22*ro(i-inx1+ns)*ro(k-inx1+ns)/xx/xx+
     &                  k1*(ro(i-inx1+ns)*ro1(k-inx1+ns))/xx+
     &                  k22*(ro1(i-inx1+ns)*ro(k-inx1+ns))/xx)*
     &                  bbc*dv_dR(ll)/4.d0
                      enddo
                    endif
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
      deallocate(t)


c******* VERY IMPORTANT  ********
c the matrix elements of dV/dR are stored in vmat at positions starting from
c element(nm+1,nm+1) this was done to because I didn't want to add any extra
c terms into the subroutine b_spline_calculation
      do i=1,nm
        do k=1,nm
          l=i+1
          m=k+1
          if(dkb) then
            do kk1=-nkap,nkap
              if(kk1.ne. 0) then
                do ll=0,2*nkap,lstep
                  dvdRmat_dkb1(i,k,kk1,ll,1)=v_22b(kk1,l,m,ll)
                  dvdRmat_dkb1(i,k,kk1,ll,2)=v_22d(kk1,l,m,ll)
                enddo
              endif
              do kk2=-nkap,nkap
                if(kk2*kk1.ne. 0) then
                  do ll=0,2*nkap,lstep
                    dvdRmat_dkb2(i,k,kk1,kk2,ll,1)=
     &              v_22smalla(kk1,kk2,l,m,ll)
                    dvdRmat_dkb2(i,k,kk1,kk2,ll,2)=
     &              v_22smallb(kk1,kk2,l,m,ll)
                  enddo
                endif
              enddo
            enddo
          endif
          do ll=0,2*nkap,lstep
            vmat(i,k,ll,nvmat)=v2(l,m,ll)
          enddo
        enddo
      enddo

c****** END OF MATRIX dV/dR

CCC      write(*,*) 'DEALLOCATING V'
      deallocate(v)
      deallocate(v2)
      nstates=0
      nsto=0
      nste=0
      number_states=0

c For the monopole case, it's a good idea to ditch the do of kk from
c -nkap to nkap, and just let it go from -nkap to -nkap

c     do kk=-nkap,-nkap
      do kk=-nkap,nkap
        if(kk.ne.0)then
          do i=1,nm
            do k=1,nm
              l=i+1
              m=k+1

              if(dkb) then
                amat(i,k)=b(l,m)+vmat(i,k,0,1)*pi4+7.5d-1*dd(l,m)+
     &          7.5d-1*kk*(kk+1)*b_2(l,m)+vmat(i,k,0,2)*pi4+
     &          kk*kk*vmat(i,k,0,3)*pi4+kk*vmat(i,k,0,4)*pi4

                amat(i+nm,k+nm)=-b(l,m)+vmat(i,k,0,1)*pi4-
     &          7.5d-1*dd(l,m)-
     &          7.5d-1*kk*(kk-1)*b_2(l,m)+vmat(i,k,0,2)*pi4+
     &          kk*kk*vmat(i,k,0,3)*pi4-kk*vmat(i,k,0,4)*pi4

                amat(i,k+nm)=vmat(i,k,0,5)*pi4+
     &          (d3(l,m)-d3(m,l))/8.d0+
     &          (div1(l,m,kk)+div2(m,l,kk))/8.d0

                amat(k+nm,i)=amat(i,k+nm)

                dmat(i,k)=b(l,m)+dd(l,m)/4.d0+kk*b_2(l,m)/
     &          4.d0*(kk+1)
                dmat(i+nm,k+nm)=b(l,m)+dd(l,m)/4.d0+
     &          kk*b_2(l,m)/4.d0*(kk-1)

c       alternate_dmat(i,k,kk)=dmat(i,k)
c       alternate_dmat(i+nm,k+nm,kk)=dmat(i+nm,k+nm)
              else
                amat(i,k+nm)=(div1(l,m,kk)+div2(l,m,kk))/2.d0

                amat(k+nm,i)=amat(i,k+nm)
              endif
            enddo
          enddo


c DMAT should now be a visible result of calling b_splines
c Use this for checking orthogonality, comment out when not in use.

          do i=1,2*nm
            do k=1,2*nm
              dmat1(i,k)=dmat(i,k)
              alternate_dmat(i,k,kk)=dmat(i,k)
            enddo
          enddo
c End of DMAT

c     Solving of the generalized problem for eigenvectors
c I am guessing that wave is vector v in eq 16 of Johnson

          call solve_equation(nm,dmat,amat,cc,wave)

          if(cc(nm+1).gt. -1.d0) then
            do i=1,2*nm
              e(i,kk)=cc(i)
              do j=1,2*nm
                wave1(j,i,kk)=wave(j,i)
              enddo
            enddo
          elseif(cc(nm+1).lt. -1.d0) then
            nrange=nm+1
            el=0.d0
            do j=1,nm+1
              if (cc(j).gt. -up_energy .and. cc(j).lt. -1.d0) then
                do l=1,nm/2+1
                  do m=1,nm/2+1
                    el(j)=el(j)+
     &              (wave(l,j)*wave(m,j)*dmat1(l,m)+wave(l+nm,j)*
     &              wave(m+nm,j)*dmat1(l+nm,m+nm))
                  enddo
                enddo
              endif
            enddo
            One_s_state_location=maxloc(el, dim=1, mask= el .lt. 1.d0)
            do i=1,2*nm
              if (i .lt. One_s_state_location) then
                e(i,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i,kk)=wave(j,i)
                enddo
              elseif(i .eq. One_s_state_location)  then
                e(nm+1,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,nm+1,kk)=wave(j,i)
                enddo
              elseif(i .gt. One_s_state_location .and.
     &        i .le. nm+1) then
                e(i-1,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i-1,kk)=wave(j,i)
                enddo
              else
                e(i,kk)=cc(i)
                do j=1,2*nm
                  wave1(j,i,kk)=wave(j,i)
                enddo
              endif
            enddo
          endif
          jjj=1
          if(Manual_ncont_states)then
            do while(jjj.le. 2*nm .and. dabs(cc(jjj)).gt.up_energy)
              jjj=jjj+1
            enddo
            number_states(1,kk)=jjj
          else
            number_states(1,kk)=nm+1
            jjj=nm+1
          endif

          do while(jjj.le.2*nm .and. cc(jjj).le.up_energy)
            if(dabs(0.5d0).le.(iabs(kk)-5.d-1))then
              nstates=nstates+1
              if(((kk.gt.0).and.((kk/2)*2.eq.kk))
     &        .or.((kk.lt.0).and.((kk/2)*2.ne.kk)))then
                nste=nste+1
              else
                nsto=nsto+1
              endif
              number_states(2,kk)=jjj
            endif
            jjj=jjj+1
          enddo
        endif
      enddo

      deallocate(amat)
      deallocate(dmat)
      deallocate(cc)
      deallocate(wave)
      deallocate(b)
      deallocate(div2)
      deallocate(div1)
      deallocate(u)
      if(dkb)then
        deallocate(b_2)
        deallocate(dd)
        deallocate(d3)
        deallocate(v_2)
        deallocate(v_22b)
        deallocate(v_22d)
        deallocate(v_22smalla)
        deallocate(v_22smallb)
        deallocate(dvd)
        deallocate(dv_1)
        deallocate(dvb)
      endif
      return
      end

      subroutine buildMultipoleBasis(nm,nkap,nstates,number_states,wave,
     &vmat,nvmat,e, ii_xi, xi_stepslower,rmin,rmax,eigval,wave_new,
     &i_even_odd_normal)
      include 'inc.par'
      integer number_states(2,-nkap:nkap), xi_stepslower
      real*8 wave(2*nm,2*nm,-nkap:nkap),vmat(nm,nm,0:2*nkap,nvmat),
     &e(2*nm,-nkap:nkap),rmin,rmax,eigval(nstates),
     &wave_new(nstates,2*nm,-nkap:nkap)
      real*8, dimension(:,:,:),allocatable:: ang
      integer, dimension(:,:), allocatable:: num_st
      real*8, dimension(:,:), allocatable:: amat, eigvec
      character*1 char_even_odd_normal

      char_even_odd_normal = 'a'
      allocate(ang(-nkap:nkap,-nkap:nkap,0:2*nkap))
      call store_angle(nkap,ang)
      allocate(num_st(-nkap:nkap,2*nm))
      allocate(amat(nstates,nstates))
      select case (i_even_odd_normal)
        case(0)
          call form_matrix(nm,nkap,nstates,number_states,num_st,
     &    amat,wave,vmat,nvmat,e,ang,1.d0)
          char_even_odd_normal = 'a'
        case(1)
          call form_matrix_even(nm,nkap,nstates,number_states,num_st,
     &    amat,wave,vmat,nvmat,e,ang,1.d0)
          char_even_odd_normal = 'e'
        case(2)
          call form_matrix_odd(nm,nkap,nstates,number_states,num_st,
     &    amat,wave,vmat,nvmat,e,ang,1.d0)
          char_even_odd_normal = 'o'
      end select
      allocate(eigvec(nstates,nstates))
      call r_diagonal(char_even_odd_normal,nstates,amat,eigval,eigvec)
      call store_new_wave(nm,nstates,nkap,num_st,wave,
     &eigvec,wave_new)
      deallocate(eigvec)
      call swapping_of_states(wave_new,nstates,nm,nkap,
     &ii_xi,xi_stepslower,rmin,rmax,eigval,char_even_odd_normal)
      call sign_of_states(wave_new,nstates,nm,nkap,
     &ii_xi,xi_stepslower,rmin,rmax,char_even_odd_normal)
      deallocate(amat)
      deallocate(num_st)
      deallocate(ang)

      return
      end
