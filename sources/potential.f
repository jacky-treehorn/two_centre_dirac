      subroutine test_potentials(nkap)
      include 'inc.par'
      real*8, dimension(:),allocatable:: u
      common /r_nuc/ r01,r02
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /nuc_mod/ nuc_model

      allocate(u(0:2*nkap))

      write(*,*) 'Enter Distance'
      read(*,*) distance
      write(*,*) 'Distance',distance
      write(*,*) 'Nuclear radii',r01,r02

      write(*,*) 'Nuclear charges',z_nuc1,z_nuc2
 1    write(*,*) 'Enter x'
      read(*,*) x
      if(x.gt.distance)then
        rg=x
        rm=distance
      else
        rm=x
        rg=distance
      endif
      rel=rm/rg
      do nuc_model=0,3
        select case(nuc_model)
        case(0)
          write(*,*) 'POINT-LIKE'
        case(1)
          write(*,*) 'Shell'
        case(2)
          write(*,*) 'Sphere'
        case(3)
          write(*,*) 'Fermi'
        end select
        if(z_nuc1.eq.z_nuc2)then
          call all_potential(nkap,x,u)
        else
          call all_potential_baryonic(nkap,x,u)
        endif
        do l=0,2*nkap
CCC            write(*,*) l,2.d0*(-az1-az2*(-1)**l)/rg*(rel)**l/(2*l+1)
          write(*,*) l,u(l)
CCC            write(*,*)
CCC     &           2.d0*(-az1-az2*(-1)**l)/rg*(rel)**l/(2*l+1)/u(l)
        enddo
        nn=5
        do j=1,nn
          xx=-1.d0+(j-1)*2.d0/(nn-1)
          vtst=0.d0
          do l=0,2*nkap
            vtst=vtst+u(l)*spharm(l,0,xx)*dsqrt(pi*(2*l+1))
          enddo
          call v_nucl(1,x,xx,vv1,vv2,v0)
          v1=vv1+v0
          call v_nucl(2,x,xx,vv1,vv2,v0)
          v2=vv2+v0
          write(*,*) j,xx
          write(*,*) (v1+v2),'analytical'
          write(*,*) vtst,'expansion'
          write(*,*) dabs((v1+v2-vtst)/vtst),'abs. accuracy'
        enddo
      enddo
      goto 1
      deallocate(u)
      return
      end

      subroutine all_potential_old(nkap,y,u)
      implicit real*8 (a-h,o-z)
      real*8 u(0:2*nkap)

      u=0.d0
      do i=0,2*nkap
        u(i)=v_nucl_l(i,y)
      enddo
      return
      end

      real*8 function v_nucl_l(l,y)
      include 'inc.par'
      common /r_nuc/ r01,r02
      external potl1,potl2
      common /dist/distance,Starting_Distance
      common /rad/ radius
      common /order/ ll
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      radius=y
      ll=l

      v_nucl_l=0.d0
      xx1=(-r01**2+y**2+distance**2)/2.d0/y/distance
      do i=1,npoints_vl
        t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
        t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
        if((xx1.le.t1).and.(xx1.ge.t0))then
          re=rint71(potl1,t0,xx1,32)+rint71(potl1,xx1,t1,32)
        else
          re=rint71(potl1,t0,t1,32)
        endif
        v_nucl_l=v_nucl_l+re
      enddo

      xx2=(r02**2-y**2-distance**2)/2.d0/y/distance
      do i=1,npoints_vl
        t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
        t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
        if((xx2.le.t1).and.(xx2.ge.t0))then
          re=rint71(potl2,t0,xx2,32)+rint71(potl2,xx2,t1,32)
        else
          re=rint71(potl2,t0,t1,32)
        endif
        v_nucl_l=v_nucl_l+re
      enddo
c      if(l.eq.0) write(*,*) 'in v_l'
c      if(l.eq.0) write(*,*) y,distance,v_nucl_l
      if(y.gt.distance)then
        v_nucl_l=v_nucl_l-2.d0*
     &  (az1+(-1)**l*az2)*distance**l/y**(l+1)/dble(2*l+1)
      else
        v_nucl_l=v_nucl_l-2.d0*
     &  (az1+(-1)**l*az2)/distance**(l+1)*y**l/dble(2*l+1)
      endif
c      if(l.eq.0) write(*,*) y,distance,v_nucl_l
c      if(l.eq.0) pause
      return
      end

      real*8 function potl1(x)
      implicit real*8(a-h,o-z)
      common /rad/ radius
      common /order/ ll
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      call v_nucl(1,radius,x,v1,v2,v0)
      potl1=v1*PLGNDR(ll,0,x)
      return
      end

      real*8 function potl2(x)
      implicit real*8(a-h,o-z)
      common /rad/ radius
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /order/ ll

      call v_nucl(2,radius,x,v1,v2,v0)
      potl2=v2*PLGNDR(ll,0,x)
      return
      end

      subroutine v_nucl(n,y,x,v1,v2,v0)
      include 'inc.par'
      common /r_nuc/ r01,r02
      common /aferm/summa1,summa2
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /dist/distance,Starting_Distance
      external r_ferma,r_ferm1a,r_fermb,r_ferm1b
      common /nuc_mod/ nuc_model

      if(dabs(x).gt.1.d0) stop '|X|>1'
      v1=0.d0
      v2=0.d0
      if((n.eq.1).and.(z_nuc1.ne.0.d0))then
        if((summa1.eq.0.d0).and.((Nuc_model.eq.3)
     &  .or.(Nuc_model.eq.4))) then
          t1a=0.d0
          t2a=r01
          re=rint(r_ferm1a,t1a,t2a,32)
          summa1=re
          do while(dabs(re/summa1).gt.1.d-16)
            t1a=t2a
            t2a=t2a*1.5d0
            re=rint(r_ferm1a,t1a,t2a,32)
            summa1=summa1+re
          enddo
        endif
        y1=dsqrt(y**2+distance**2-2.d0*y*distance*x)
        v0=-az1/y1
        select case(Nuc_model)
          case(0)
            v1=v0
          case(1)
            if(y1.gt.r01)then
              v1=-az1/y1
            else
              v1=-az1/r01
            endif
          case(2)
            if(y1.gt.r01)then
              v1=-az1/y1
            else
              v1=-az1/r01*(1.5d0-((y1/r01)**2)/2.d0)
            endif
          case(3)
            t0=0.d0
            t1=1.d-2
            nuzz=32
            v1=0.d0

            do while(t1.lt.y1)
              su=rint7(r_ferm1a,t0,t1,nuzz)
              v1=v1+su
              t0=t1
              t1=t1*1.9d0
            enddo

            su=rint7(r_ferm1a,t0,y1,nuzz)
            v1=v1+su
            t0=y1
            t1=y1*1.9d0
            su=v1

            do while(dabs(su/v1).gt.1.d-16)
              su=rint7(r_ferma,t0,t1,nuzz)
              v1=v1+su*y1
              t0=t1
              t1=t1*1.9d0
            enddo
            v1=-az1*v1/y1/summa1
        end select
        v1=v1-v0
      elseif(z_nuc2.ne.0.d0)then
        if((summa2.eq.0.d0).and.((Nuc_model.eq.3)
     &  .or.(Nuc_model.eq.4))) then
          t1a=0.d0
          t2a=r01
          re=rint(r_ferm1b,t1a,t2a,32)
          summa2=re
          do while(dabs(re/summa2).gt.1.d-16)
            t1a=t2a
            t2a=t2a*1.5d0
            re=rint(r_ferm1b,t1a,t2a,32)
            summa2=summa2+re
          enddo
        endif
        y2=dsqrt(y**2+distance**2+2.d0*y*distance*x)

        v0=-az2/y2
        select case(Nuc_model)
          case(0)
            v2=-az2/y2
          case(1)
            if(y2.gt.r02)then
              v2=-az2/y2
            else
              v2=-az2/r02
            endif
          case(2)
            if(y2.gt.r02)then
              v2=-az2/y2
            else
              v2=-az2/r02*(1.5d0-((y2/r02)**2)/2.d0)
            endif
          case(3)
            t0=0.d0
            t1=1.d-2
            nuzz=32
            v2=0.d0

            do while(t1.lt.y2)
              su=rint7(r_ferm1b,t0,t1,nuzz)
              v2=v2+su
              t0=t1
              t1=t1*1.9d0
            enddo

            su=rint7(r_ferm1b,t0,y2,nuzz)
            v2=v2+su
            t0=y2
            t1=y2*1.9d0
            su=v2

            do while(dabs(su/v2).gt.1.d-16)
              su=rint7(r_fermb,t0,t1,nuzz)
              v2=v2+su*y2
              t0=t1
              t1=t1*1.9d0
            enddo
            v2=-az2*v2/y2/summa2
        end select
        v2=v2-v0
      endif
      return
      end

      real*8 function r_ferma(x)
      implicit real*8(a-h,o-z)
      common /ferm/ c1,z1,c2,z2

      x1=x/(25.896063451984380d-4)
      if(((x1-c1)/z1).gt.500) then
        r_ferma=dexp(-(x1-c1)/z1)*x
      else
        r_ferma=1.d0/(1.d0+dexp((x1-c1)/z1))*x
      endif
      return
      end

      real*8 function r_ferm1a(x)
      implicit real*8(a-h,o-z)
      common /ferm/ c1,z1,c2,z2

      x1=x/(25.896063451984380d-4)
      if(((x1-c1)/z1).gt.500) then
        r_ferm1a=dexp(-(x1-c1)/z1)*x*x
      else
        r_ferm1a=1.d0/(1.d0+dexp((x1-c1)/z1))*x*x
      endif
      return
      end


      real*8 function r_fermb(x)
      implicit real*8(a-h,o-z)
      common /ferm/ c1,z1,c2,z2

      x1=x/(25.896063451984380d-4)
      if(((x1-c2)/z2).gt.500) then
        r_fermb=dexp(-(x1-c2)/z2)*x
      else
        r_fermb=1.d0/(1.d0+dexp((x1-c2)/z2))*x
      endif
      return
      end

      real*8 function r_ferm1b(x)
      implicit real*8(a-h,o-z)
      common /ferm/ c1,z1,c2,z2

      x1=x/(25.896063451984380d-4)
      if(((x1-c2)/z2).gt.500) then
        r_ferm1b=dexp(-(x1-c2)/z2)*x*x
      else
        r_ferm1b=1.d0/(1.d0+dexp((x1-c2)/z2))*x*x
      endif
      return
      end

      subroutine store_angle(nkap,ang)
      include 'inc.par'
      real*8 ang(-nkap:nkap,-nkap:nkap,0:2*nkap)
      common /momentum_projection/ amu,amj_max

      ang=0.d0

      do k1=-nkap,nkap
        aj1=iabs(k1)-5.d-1
        if(k1.lt.0)then
          al1=-k1-1
        else
          al1=k1
        endif
        do k2=-nkap,nkap
          aj2=iabs(k2)-5.d-1
          if(k2.lt.0)then
            al2=-k2-1
          else
            al2=k2
          endif
          if(k1*k2.ne.0)then
c$$$               write(*,*) k1,aj1,al1
c$$$               write(*,*) k2,aj2,al2
c$$$               write(*,*) idint(dabs(aj1-aj2)),idint(aj1+aj2)
c$$$               pause
            do l=idint(dabs(aj1-aj2)),idint(aj1+aj2)
              al=dble(l)
              symb1=symb6j(al1,al2,al,aj2,aj1,5.d-1)
              if(symb1.ne.0.d0)then
c$$$                     write(*,*) symb1,'SYMB1',
c$$$     &                    symb6j(al1,al2,al,aj2,aj1,5.d-1)
c$$$                     write(*,*) al1,al2,al
c$$$                     write(*,*) aj1,aj2,5.d-1
                ip=idint(aj1+aj2+al+amu+5.d-1)
                dfact=dble((-1)**ip)/2.d0
                dfact=dfact*dsqrt(u2(aj1)*u2(aj2)*u2(al1)*u2(al2))
c$$$                     write(*,*) dfact,'DFACT'
                symb2=clebsh(al1,0.d0,al2,0.d0,al,0.d0)
                symb3=clebsh(aj1,-amu,aj2,amu,al,0.d0)
c$$$                     write(*,*) symb2,'SYMB2'
c$$$                     write(*,*) symb3,'SYMB3'
                ang(k1,k2,l)=symb1*symb2*symb3*dfact
                if(dabs(ang(k1,k2,l)).lt.1.d-12) ang(k1,k2,l)=0.d0
              endif
            enddo
          endif
        enddo
      enddo
c$$$      do k1=-nkap,nkap
c$$$         do k2=-nkap,k1
c$$$            do l=0,2*nkap
c$$$               if(ang(k1,k2,l).ne.0.d0)write(*,*) k1,k2,l,ang(k1,k2,l)
c$$$            enddo
c$$$            pause
c$$$         enddo
c$$$      enddo
      return
      end

      real*8 function u2(x)
      real*8 x
      u2=2.d0*x+1.d0
      return
      end

      real*8 function s(j,l,kb,ka)
      implicit real*8(a-h,o-z)

      s=0.d0
      if((l.eq.j+1).and.(j.gt.0)) then

        s=dsqrt(dble(j+1)/dble(2*j+1))*
     &  (1.d0+dble(ka+kb)/dble(j+1))*c(j,-ka,kb)

      elseif ((l.eq.j).and.(j.gt.0)) then

        s=dble(kb-ka)/dsqrt(dble(j*(j+1)))*c(j,ka,kb)

      elseif ((l.eq.j-1).and.(j.gt.0)) then

        s=dsqrt(dble(j)/dble(2*j+1))*
     &  (-1.d0+dble(ka+kb)/dble(j))*c(j,-ka,kb)

      elseif((j.eq.0).and.(l.eq.1)) then

        s=c(j,-ka,kb)
      endif

      return
      end


      real*8 function c(j,kb,ka)
      implicit real*8(a-h,o-z)

      c=0.d0

      aja=iabs(ka)-5.d-1
      ajb=iabs(kb)-5.d-1
      if(ka.lt.0) then
        la=-ka-1
      else
        la=ka
      endif

      if(kb.lt.0) then
        lb=-kb-1
      else
        lb=kb
      endif

      i1=la+lb+j
      i2=((la+lb+j)/2)*2
      if(i1.ne.i2) return

      c=(-1)**(ajb+5.d-1)*
     &dsqrt((2.d0*aja+1.d0)*(2.d0*ajb+1.d0))*
     &symb3j(aja,dble(j),ajb,5.d-1,0.d0,-5.d-1)

      return
      end

      real*8 function symb9j(aj1,aj2,aj3,aj4,aj5,aj6,aj7,aj8,aj9)
      implicit real*8 (a-h,o-z)

      akmin=dabs(aj1-aj9)
      akmax=aj1+aj9

      akl=dabs(aj2-aj6)
      akr=aj2+aj6

      if(akmin.lt.akl) akmin=akl
      if(akmax.gt.akr) akmax=akr

      akl=dabs(aj4-aj8)
      akr=aj4+aj8

      if(akmin.lt.akl) akmin=akl
      if(akmax.gt.akr) akmax=akr

      symb9j=0.d0

      do i=int(akmin*2.d0),int(akmax*2.d0),2
      ak=dble(i)/2.d0

      symb9j=symb9j+(-1)**i*(2.d0*ak+1.d0)*
     &symb6j(aj1,aj4,aj7,aj8,aj9,ak)*symb6j(aj2,aj5,aj8,aj4,ak,aj6)*
     &symb6j(aj3,aj6,aj9,ak,aj1,aj2)

      enddo

      return
      end

      double precision function symb6j(cja,cjb,cje,cjd,cjc,cjf)
      implicit real*8(a-h,o-z)
      common /ftal/ftl(300),srftl(300)
c
c 10   ka1 = cja+cjb+cje
      ka1 = nint(cja+cjb+cje)
      ka2 = nint(cja+cjc+cjf)
      ka3 = nint(cjb+cjd+cjf)
      ka4 = nint(cjc+cjd+cje)
      kb1 = nint(cja+cjb+cjc+cjd)
      kb2 = nint(cja+cjd+cje+cjf)
      kb3 = nint(cjb+cjc+cje+cjf)
      mx = max0(ka1,ka2,ka3,ka4)
      mn = min0(kb1,kb2,kb3)
      if(mx.gt.mn) go to 3
      jmax = max0(kb1,kb2,kb3)
      jmin = min0(ka1,ka2,ka3,ka4)
      in = max0(jmax-jmin,mx+1)
c      if(in.gt.60) stop '6j-symb. error'
      a1 = (srftl(kb1-ka1+1)+srftl(kb2-ka1+1)+
     &srftl(kb3-ka1+1)-srftl(ka1+2))
      a2 = (srftl(kb1-ka2+1)+srftl(kb2-ka2+1)+
     &srftl(kb3-ka2+1)-srftl(ka2+2))
      a3 = (srftl(kb1-ka3+1)+srftl(kb2-ka3+1)+
     &srftl(kb3-ka3+1)-srftl(ka3+2))
      a4 = (srftl(kb1-ka4+1)+srftl(kb2-ka4+1)+
     &srftl(kb3-ka4+1)-srftl(ka4+2))
      s = 0.d0
      m = mx-1
2     continue
      m = m+1
      b1 = ftl(m-ka1+1)
      b2 = ftl(m-ka2+1)
      b3 = ftl(m-ka3+1)
      b4 = ftl(m-ka4+1)
      c1 = ftl(kb1-m+1)
      c2 = ftl(kb2-m+1)
      c3 = ftl(kb3-m+1)
      c = ftl(m+2)
      s = s+dexp(a1-b1-c1+a2-b2-c2+a3-b3-c3+a4+c-b4)*(-1.d0)**m
      if (m.lt.mn) go to 2
      symb6j = s
      return
3     continue
      symb6j = 0.d0
      return
      end

c******************************************************************************

      real*8 function symb3j(pj1,pj2,pj3,pm1,pm2,pm3_)
      implicit real*8 (a-h,o-z)
      parameter (Nftal = 300)
      common /ftal/ftl(Nftal),srftl(Nftal)
      data d0/0.d0/ d1/1.d0/
c
c calculation of 3j-symbol:   pj1 pj2  pj3
c                                 pm1 pm2  pm3_
c

      pm3 = -pm3_
      symb3j = d0
      pjs = pj1+ pj2+ pj3
      pk1 = dabs(pj1-pj2)
      pk2 = pj1 +pj2
      if (pj3.lt.pk1.or.pj3.gt.pk2) return
      if (dabs(pm1+pm2+pm3_).gt.1.d-10) return
      if (dabs(pm1).gt.pj1.or.dabs(pm2).gt.pj2.or.dabs(pm3).gt.pj3) then
        return
      endif
      if (pm1.eq.d0.and.pm2.eq.d0.and.pjs.eq.d0) go to 3
      k1 = 0
      k2 = idint( pj1- pj2+ pm3)
      k3 = idint( pj3- pj2+ pm1)

      l1 = idint( pj1- pj2 +pj3)
      l2 = idint( pj3+ pm3)
      l3 = idint( pj1+ pm1)

      ka = idint( pj1+ pj2+ pj3)+ 1
      ks = k1+ k2+ k3
      mx = max0(k1,k2,k3)
      mn = min0(l1,l2,l3)

      if (mx.gt.mn) return
      in = max0(max0(l1,l2,l3)- min0(k1,k2,k3),ka)

      a1 = srftl(l1-k1+1)
      a2 = srftl(l2-k1+1)
      a3 = srftl(l3-k1+1)
      a4 = srftl(l1-k2+1)
      a5 = srftl(l2-k2+1)
      a6 = srftl(l3-k2+1)
      a7 = srftl(l1-k3+1)
      a8 = srftl(l2-k3+1)
      a9 = srftl(l3-k3+1)
      a10 = srftl(ka+1)

      s = d0
      m = mx-1
 2    m = m+1
      b1 = ftl(m-k1+1)
      b2 = ftl(m-k2+1)
      b3 = ftl(m-k3+1)
      b4 = ftl(l1-m+1)
      b5 = ftl(l2-m+1)
      b6 = ftl(l3-m+1)

      ffact=dexp(a1+a2+a3+a4+a5+a6+a7+a8+a9-a10-
     &b1-b2-b3-b4-b5-b6)

      s = s+ffact*(-1.d0)**m
      if (m.lt.mn) go to 2
      s = s* (-1.d0)**ks
c     symb3j = s* dsqrt(d2*pj3+ d1)* (-1.d0)**k2
      symb3j = s

      return
3     continue
      symb3j = d1
      return
      end


      real*8 function Clebsh(pj1,pm1,pj2,pm2,pj3,pm3)
      implicit real*8 (a-h,o-z)
      parameter (Nftal = 300)
      common /ftal/ftl(Nftal),srftl(Nftal)
      data d0/0.d0/ d1/1.d0/ d2/2.d0/ !d3/3.d0/
c
c       Clebsch-Gordan coefficients are calculated
c           pj3,pm3
c         C
c           pj1,pm1;pj2,pm2
c
      if (ftl(3).eq.0.d0) stop 'Clebsh: /ftal/ is empty!'
      Clebsh = d0
      pjs = pj1+ pj2+ pj3
      pk1 = dabs(pj1-pj2)
      pk2 = pj1 +pj2
      if (pj3.lt.pk1.or.pj3.gt.pk2) return
      if (dabs(pm1+pm2-pm3).gt.1.d-10) return
      if (dabs(pm1).gt.pj1.or.dabs(pm2).gt.pj2.or.dabs(pm3).gt.pj3) then
        return
      endif
      if (pm1.eq.d0.and.pm2.eq.d0.and.pjs.eq.d0) go to 3
      k1 = 0
      k2 = idint( pj1- pj2+ pm3)
      k3 = idint( pj3- pj2+ pm1)
      l1 = idint( pj1- pj2 +pj3)
      l2 = idint( pj3+ pm3)
      l3 = idint( pj1+ pm1)
      ka = idint( pj1+ pj2+ pj3)+ 1
      ks = k1+ k2+ k3
      mx = max0(k1,k2,k3)
      mn = min0(l1,l2,l3)
      if (mx.gt.mn) return
      in = max0(max0(l1,l2,l3)- min0(k1,k2,k3),ka)
      if (in.ge.Nftal) then
        stop 'Clebsh: in > Nftal !'
      endif
      a1 = srftl(l1-k1+1)
      a2 = srftl(l2-k1+1)
      a3 = srftl(l3-k1+1)
      a4 = srftl(l1-k2+1)
      a5 = srftl(l2-k2+1)
      a6 = srftl(l3-k2+1)
      a7 = srftl(l1-k3+1)
      a8 = srftl(l2-k3+1)
      a9 = srftl(l3-k3+1)
      a10 = srftl(ka+1)
      s = d0
      m = mx-1
   2  m = m+1
      b1 = ftl(m-k1+1)
      b2 = ftl(m-k2+1)
      b3 = ftl(m-k3+1)
      b4 = ftl(l1-m+1)
      b5 = ftl(l2-m+1)
      b6 = ftl(l3-m+1)
      ffactor=dexp(a1+a2+a3+a4+a5+a6+a7+a8+a9-a10-
     &b1-b2-b3-b4-b5-b6)
      s = s+ffactor*(-1.d0)**m
      if (m.lt.mn) go to 2
      s = s* (-1.d0)**ks
      Clebsh = s* dsqrt(d2*pj3+ d1)* (-1.d0)**k2
      return

   3  continue
      Clebsh = d1
      return
      end

      subroutine factor()
c
c     the subroutine generates the table ftl of  ln of factorials
c
      implicit real*8 (a-h,o-z)
      parameter (Nftal = 300)
      common /ftal/ftl(Nftal),srftl(Nftal)
!      data d0/0.d0/ d1/1.d0/ d2/2.d0/ d3/3.d0/
      ftl(1) = 0.d0
      srftl(1) = 0.d0
      do 1 i = 1,Nftal-1
      ftl(i+1) = dlog(dble(i))+ ftl(i)
      srftl(i+1) = ftl(i+1)/2.d0
1     continue
      return
      end


      FUNCTION SPHARM(L,M,X)
C
C     SPHERICAL HARMONICS FOR PHI = 0 AND X = COS(THETA)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DATA PI /3.141592653589793D0/
      parameter (Nftal = 300)
      common /ftal/ftl(Nftal),srftl(Nftal)

      MBAR=IABS(M)
c########
c     Modified by ANA
c########
      spharm=0.d0
      if(mbar.gt.l) return
c########

      if (L+MBAR.ge.Nftal) stop 'Spharm: Nftal is too small'
      FLMM=dexp(FTL(L-MBAR+1))
      FLPM=dexp(FTL(L+MBAR+1))
      SPHARM=DSQRT(((2.D0*L+1.D0)*FLMM)/(4.D0*PI*FLPM))*
     &PLGNDR(L,MBAR,X)
      IF (M.LT.0) SPHARM=(-1)**M* SPHARM

      RETURN
      END


      subroutine all_potential(nkap,y,u)
      include 'inc.par'
      common /r_nuc/ r01,r02
      external potl1,potl2
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 u(0:2*nkap)
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_mod/ nuc_model

      u=0.d0
      radius=y
      if(z_nuc1.eq.z_nuc2)then
        istep=2
        ff=2.d0
      else
        istep=1
        ff=1.d0
      endif

      do l=0,2*nkap,istep
        ll=l
        if(y.ge.distance)then
          u(l)=
     &    -2.d0*(az1+(-1)**l*az2)*distance**l/y**(l+1)/dble(2*l+1)
        else
          u(l)=
     &    -2.d0*(az1+(-1)**l*az2)/distance**(l+1)*y**l/dble(2*l+1)
        endif
      enddo

      if(nuc_model.gt.0)then
        xx1=(-r01**2+y**2+distance**2)/2.d0/y/distance
        do i=1,npoints_vl
          t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
          t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
          if((xx1.le.t1).and.(xx1.ge.t0))then
            a=t0
            b=xx1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
            a=xx1
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          else
            a=t0
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          endif
        enddo

        if(z_nuc1.ne.z_nuc2)then
          xx1=(r02**2-y**2-distance**2)/2.d0/y/distance
          do i=1,npoints_vl
            t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
            t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
            if((xx1.le.t1).and.(xx1.ge.t0))then
              a=t0
              b=xx1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
              a=xx1
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            else
              a=t0
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            endif
          enddo
        endif
      endif
      return
      end

      subroutine all_potential_dist(nkap,y,u,distance)
      include 'inc.par'
      common /r_nuc/ r01,r02
      external potl1,potl2
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 u(0:2*nkap)
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_mod/ nuc_model

      u=0.d0
      radius=y
      if(z_nuc1.eq.z_nuc2)then
        istep=2
        ff=2.d0
      else
        istep=1
        ff=1.d0
      endif

      do l=0,2*nkap,istep
        ll=l
        if(y.ge.distance)then
          u(l)=
     &    -2.d0*(az1+(-1)**l*az2)*distance**l/y**(l+1)/dble(2*l+1)
        else
          u(l)=
     &    -2.d0*(az1+(-1)**l*az2)/distance**(l+1)*y**l/dble(2*l+1)
        endif
      enddo

      if(nuc_model.gt.0)then
        xx1=(-r01**2+y**2+distance**2)/2.d0/y/distance
        do i=1,npoints_vl
          t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
          t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
          if((xx1.le.t1).and.(xx1.ge.t0))then
            a=t0
            b=xx1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
            a=xx1
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          else
            a=t0
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          endif
        enddo

        if(z_nuc1.ne.z_nuc2)then
          xx1=(r02**2-y**2-distance**2)/2.d0/y/distance
          do i=1,npoints_vl
            t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
            t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
            if((xx1.le.t1).and.(xx1.ge.t0))then
              a=t0
              b=xx1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
              a=xx1
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            else
              a=t0
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            endif
          enddo
        endif
      endif
      return
      end

ccc***************BARYONIC CENTRE POTENTIAL**********************************

      subroutine all_potential_baryonic(nkap,y,u)
      include 'inc.par'
      common /r_nuc/ r01,r02
      external potl1,potl2
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 u(0:2*nkap)
      common /Barycentres/ RadiusOne,RadiusTwo
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_mod/ nuc_model

ccc RadiusOne= 2.d0*distance*z_nuc2/(z_nuc1+z_nuc2)
ccc RadiusTwo= 2.d0*distance*z_nuc1/(z_nuc1+z_nuc2)

      u=0.d0
      radius=y
      if(z_nuc1.eq.z_nuc2)then
        istep=2
        ff=2.d0
      else
        istep=1
        ff=1.d0
      endif

      do l=0,2*nkap,istep
        if(y.ge.RadiusTwo)then
          u(l)=-2.d0*az1*dexp(dble(l)*dlog(RadiusOne))/
     &    dexp(dble(l+1)*dlog(y))/dble(2*l+1)-
     &    2.d0*(-1)**l*az2*dexp(dble(l)*dlog(RadiusTwo))/
     &    dexp(dble(l+1)*dlog(y))/dble(2*l+1)
        elseif(y.lt.RadiusTwo .and. y.ge.RadiusOne) then
          u(l)=-2.d0*az1*dexp(dble(l)*dlog(RadiusOne))/
     &    dexp(dble(l+1)*dlog(y))/dble(2*l+1)-
     &    2.d0*(-1)**l*az2*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusTwo))/dble(2*l+1)
        else
          u(l)=-2.d0*az1*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusOne))/dble(2*l+1)-
     &    2.d0*(-1)**l*az2*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusTwo))/dble(2*l+1)
        endif
      enddo

c      do l=0,2*nkap,istep
c         if(y.ge.RadiusTwo)then
c            u(l)=-2.d0*az1*RadiusOne**l/y**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*RadiusTwo**l/y**(l+1)/(2*l+1)
c  elseif(y.lt.RadiusTwo .and. y.ge.RadiusOne) then
c     u(l)=-2.d0*az1*RadiusOne**l/y**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*y**l/RadiusTwo**(l+1)/(2*l+1)
c  else
c            u(l)=-2.d0*az1*y**l/RadiusOne**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*y**l/RadiusTwo**(l+1)/(2*l+1)
c         endif
c      enddo

      if(nuc_model.gt.0)then
        xx1=(-r01**2+y**2+distance**2)/2.d0/y/distance
        do i=1,npoints_vl
          t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
          t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
          if((xx1.le.t1).and.(xx1.ge.t0))then
            a=t0
            b=xx1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
            a=xx1
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          else
            a=t0
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          endif
        enddo

        if(z_nuc1.ne.z_nuc2)then
          xx1=(r02**2-y**2-distance**2)/2.d0/y/distance
          do i=1,npoints_vl
            t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
            t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
            if((xx1.le.t1).and.(xx1.ge.t0))then
              a=t0
              b=xx1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
              a=xx1
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            else
              a=t0
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            endif
          enddo
        endif
      endif
      return
      end

      subroutine dvdr_potential(nkap,y,u)
      include 'inc.par'
      common /r_nuc/ r01,r02
      external potl1,potl2
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 u(0:2*nkap)
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_mod/ nuc_model

      u=0.d0
      radius=y
      if(z_nuc1.eq.z_nuc2)then
        istep=2
        ff=2.d0
      else
        istep=1
        ff=1.d0
      endif

      do l=0,2*nkap,istep
        ll=l
        if(y.ge.distance)then
          u(l)=
     &    -(az1+(-1)**l*az2)*l*distance**(l-1)/y**(l+1)/2.d0
        else
          u(l)=
     &    -(az1+(-1)**l*az2)*(-l-1)*y**l/distance**(l+2)/2.d0
        endif
      enddo

      if(nuc_model.gt.0)then
        xx1=(-r01**2+y**2+distance**2)/2.d0/y/distance
        do i=1,npoints_vl
          t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
          t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
          if((xx1.le.t1).and.(xx1.ge.t0))then
            a=t0
            b=xx1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
            a=xx1
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          else
            a=t0
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          endif
        enddo

        if(z_nuc1.ne.z_nuc2)then
          xx1=(r02**2-y**2-distance**2)/2.d0/y/distance
          do i=1,npoints_vl
            t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
            t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
            if((xx1.le.t1).and.(xx1.ge.t0))then
              a=t0
              b=xx1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                 xx=ac+bc*t32(j)
                 call v_nucl(2,y,xx,v2,v1,v0)
                 v1=v1*bc
                 do l=0,2*nkap,istep
                    u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                 enddo
              enddo
              a=xx1
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                 xx=ac+bc*t32(j)
                 call v_nucl(2,y,xx,v2,v1,v0)
                 v1=v1*bc
                 do l=0,2*nkap,istep
                    u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                 enddo
              enddo
            else
              a=t0
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            endif
          enddo
        endif
      endif
      return
      end

      subroutine dvdR_baryonic(nkap,y,u)
      include 'inc.par'
      common /r_nuc/ r01,r02
      external potl1,potl2
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 u(0:2*nkap)
      common /Barycentres/ RadiusOne,RadiusTwo
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_mod/ nuc_model

ccc RadiusOne= 2.d0*distance*z_nuc2/(z_nuc1+z_nuc2)
ccc RadiusTwo= 2.d0*distance*z_nuc1/(z_nuc1+z_nuc2)
      R=RadiusOne+RadiusTwo
      u=0.d0
      radius=y
      if(z_nuc1.eq.z_nuc2)then
        istep=2
        ff=2.d0
      else
        istep=1
        ff=1.d0
      endif

      do l=0,2*nkap,istep
        if(y.ge.RadiusTwo)then
          u(l)=(-az1*dble(l)*dexp(dble(l)*dlog(RadiusOne))/
     &    dexp(dble(l+1)*dlog(y))-
     &    (-1)**l*az2*dble(l)*dexp(dble(l)*dlog(RadiusTwo))/
     &    dexp(dble(l+1)*dlog(y)))/R
        elseif(y.lt.RadiusTwo .and. y.ge.RadiusOne) then
          u(l)=-az1*dble(L)*dexp(dble(L)*dlog(RadiusOne))/
     &    dexp(dble(L+1)*dlog(y))/R-
     &    (-1)**L*az2*dble(-L-1)*dexp(dble(L)*dlog(y))/
     &    dexp(dble(L+1)*dlog(RadiusTwo))/R
        else
          u(l)=(-az1*dble(-l-1)*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusOne))-
     &    (-1)**l*az2*dble(-l-1)*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusTwo)))/R
        endif
      enddo

c      do l=0,2*nkap,istep
c         if(y.ge.RadiusTwo)then
c            u(l)=-2.d0*az1*RadiusOne**l/y**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*RadiusTwo**l/y**(l+1)/(2*l+1)
c  elseif(y.lt.RadiusTwo .and. y.ge.RadiusOne) then
c     u(l)=-2.d0*az1*RadiusOne**l/y**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*y**l/RadiusTwo**(l+1)/(2*l+1)
c  else
c            u(l)=-2.d0*az1*y**l/RadiusOne**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*y**l/RadiusTwo**(l+1)/(2*l+1)
c         endif
c      enddo

      if(nuc_model.gt.0)then
        xx1=(-r01**2+y**2+distance**2)/2.d0/y/distance
        do i=1,npoints_vl
          t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
          t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
          if((xx1.le.t1).and.(xx1.ge.t0))then
            a=t0
            b=xx1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
            a=xx1
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          else
            a=t0
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          endif
        enddo

        if(z_nuc1.ne.z_nuc2)then
          xx1=(r02**2-y**2-distance**2)/2.d0/y/distance
          do i=1,npoints_vl
            t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
            t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
            if((xx1.le.t1).and.(xx1.ge.t0))then
              a=t0
              b=xx1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
              a=xx1
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            else
              a=t0
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            endif
          enddo
        endif
      endif
      return
      end

      subroutine all_potential_baryonic_dist(nkap,y,u,dist)
      include 'inc.par'
      common /r_nuc/ r01,r02
      external potl1,potl2
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /TargProj_prop/ Proj_mass,Targ_mass,Proj_vel,
     & b_ImpactParam
      real*8 u(0:2*nkap),dist
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)
      common /nuc_mod/ nuc_model
      common /common_dkb/ dkb
      common /dist/distance,Starting_Distance
      logical dkb

      RadiusOne= 2.d0*dist*Proj_mass/(Proj_mass+Targ_mass)
      RadiusTwo= 2.d0*dist*Targ_mass/(Proj_mass+Targ_mass)

      u=0.d0
      radius=y
      if(z_nuc1.eq.z_nuc2)then
        istep=2
        ff=2.d0
      else
        istep=1
        ff=1.d0
      endif

      do l=0,2*nkap,istep
        if(y.ge.RadiusTwo)then
          u(l)=-2.d0*az1*dexp(dble(l)*dlog(RadiusOne))/
     &    dexp(dble(l+1)*dlog(y))/dble(2*l+1)-
     &    2.d0*(-1)**l*az2*dexp(dble(l)*dlog(RadiusTwo))/
     &    dexp(dble(l+1)*dlog(y))/dble(2*l+1)
        elseif(y.lt.RadiusTwo .and. y.ge.RadiusOne) then
          u(l)=-2.d0*az1*dexp(dble(l)*dlog(RadiusOne))/
     &    dexp(dble(l+1)*dlog(y))/dble(2*l+1)-
     &    2.d0*(-1)**l*az2*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusTwo))/dble(2*l+1)
        else
          u(l)=-2.d0*az1*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusOne))/dble(2*l+1)-
     &    2.d0*(-1)**l*az2*dexp(dble(l)*dlog(y))/
     &    dexp(dble(l+1)*dlog(RadiusTwo))/dble(2*l+1)
        endif
      enddo

c      do l=0,2*nkap,istep
c         if(y.ge.RadiusTwo)then
c            u(l)=-2.d0*az1*RadiusOne**l/y**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*RadiusTwo**l/y**(l+1)/(2*l+1)
c  elseif(y.lt.RadiusTwo .and. y.ge.RadiusOne) then
c     u(l)=-2.d0*az1*RadiusOne**l/y**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*y**l/RadiusTwo**(l+1)/(2*l+1)
c  else
c            u(l)=-2.d0*az1*y**l/RadiusOne**(l+1)/(2*l+1)-
c     &     2.d0*(-1)**l*az2*y**l/RadiusTwo**(l+1)/(2*l+1)
c         endif
c      enddo

      if(nuc_model.gt.0)then
        xx1=(-r01**2+y**2+distance**2)/2.d0/y/distance
        do i=1,npoints_vl
          t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
          t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
          if((xx1.le.t1).and.(xx1.ge.t0))then
            a=t0
            b=xx1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
            a=xx1
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          else
            a=t0
            b=t1
            ac=(a+b)/2.d0
            bc=(b-a)/2.d0
            do j=1,32
              xx=ac+bc*t32(j)
              call v_nucl(1,y,xx,v1,v2,v0)
              v1=v1*bc*ff
              do l=0,2*nkap,istep
                u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
              enddo
            enddo
          endif
        enddo

        if(z_nuc1.ne.z_nuc2)then
          xx1=(r02**2-y**2-distance**2)/2.d0/y/distance
          do i=1,npoints_vl
            t0=-1.d0+dble(i-1)/dble(npoints_vl)*2.d0
            t1=-1.d0+dble(i)/dble(npoints_vl)*2.d0
            if((xx1.le.t1).and.(xx1.ge.t0))then
              a=t0
              b=xx1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
              a=xx1
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            else
              a=t0
              b=t1
              ac=(a+b)/2.d0
              bc=(b-a)/2.d0
              do j=1,32
                xx=ac+bc*t32(j)
                call v_nucl(2,y,xx,v2,v1,v0)
                v1=v1*bc
                do l=0,2*nkap,istep
                  u(l)=u(l)+w32(j)*v1*PLGNDR(l,0,xx)
                enddo
              enddo
            endif
          enddo
        endif
      endif
      return
      end
