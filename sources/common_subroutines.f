
      real*8 function rint7(f,a,b,nuz)
      implicit real*8 (a-h,o-z)
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)

      rint7=0.d0
      if(a.eq.b)return
      ac=(a+b)/2.d0
      bc=(b-a)/2.d0
      select case (nuz)
        case(4)
          do i=nuz,1,-1
            rint7=rint7+w4(i)*f(ac+bc*t4(i))
          enddo

        case(8)
          do i=nuz,1,-1
            rint7=rint7+w8(i)*f(ac+bc*t8(i))
          enddo

        case(16)
          do i=nuz,1,-1
            rint7=rint7+w16(i)*f(ac+bc*t16(i))
          enddo

        case(32)
          do i=nuz,1,-1
            rint7=rint7+w32(i)*f(ac+bc*t32(i))
          enddo

        case(64)
          do i=nuz,1,-1
            rint7=rint7+w64(i)*f(ac+bc*t64(i))
          enddo
      end select

      rint7=rint7*bc

      return
      end

      real*8 function rint(f,a,b,nuz)
      implicit real*8 (a-h,o-z)
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)

      rint=0.d0
      if(a.eq.b)return
      ac=(a+b)/2.d0
      bc=(b-a)/2.d0
      select case (nuz)
        case(4)
          do i=nuz,1,-1
            rint=rint+w4(i)*f(ac+bc*t4(i))
        enddo

        case(8)
          do i=nuz,1,-1
            rint=rint+w8(i)*f(ac+bc*t8(i))
          enddo

        case(16)
          do i=nuz,1,-1
            rint=rint+w16(i)*f(ac+bc*t16(i))
          enddo

        case(32)
          do i=nuz,1,-1
            rint=rint+w32(i)*f(ac+bc*t32(i))
        enddo

        case(64)
          do i=nuz,1,-1
            rint=rint+w64(i)*f(ac+bc*t64(i))
        enddo
      end select

      rint=rint*bc

      return
      end

      real*8 function rint71(f,a,b,nuz)
      implicit real*8 (a-h,o-z)
      common /weights/ w4(4),t4(4),w8(8),t8(8),w16(16),t16(16),w32(32),
     &t32(32),w64(64),t64(64),t6(6),w6(6)

      rint71=0.d0
      if(a.eq.b)return
      ac=(a+b)/2.d0
      bc=(b-a)/2.d0

      select case (nuz)
        case(4)
          do i=nuz,1,-1
            rint71=rint71+w4(i)*f(ac+bc*t4(i))
          enddo

        case(8)
          do i=nuz,1,-1
            rint71=rint71+w8(i)*f(ac+bc*t8(i))
          enddo

        case(16)
          do i=nuz,1,-1
             rint71=rint71+w16(i)*f(ac+bc*t16(i))
          enddo

        case(32)
          do i=nuz,1,-1
            rint71=rint71+w32(i)*f(ac+bc*t32(i))
          enddo

        case(64)
          do i=nuz,1,-1
            rint71=rint71+w64(i)*f(ac+bc*t64(i))
          enddo
      end select

      rint71=rint71*bc

      return
      end

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      implicit real*8(a-h,o-z)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(*,*) 'singular matrix in function gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) then
          write(*,*) 'singular matrix in function gaussj'
        endif
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D.


      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1
      dif=abs(x-xa(1))
      do i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0) then
            do kkk=1,N
              write(*,*) kkk,xa(kkk),ya(kkk)
            enddo
            write(*,*) 'failure in function polint'
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo
      return
      END

      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      implicit real*8 (a-h,o-z)
      PARAMETER (NMAX=10,TINY=1.d-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)

      su=0.d0

      do i=1,n
        if(dabs(xa(i)-x).lt.1.d-12) then
          y = ya(i)
          return
        endif
        if(i.lt.n) su=su+dabs((ya(i+1)-ya(i))/ya(i))
      enddo
      if(su.le.1.d-12)then
         y=ya(1)
         return
      endif
      NS=1
      HH=dABS(X-XA(1))
      DO 11 I=1,N
        H=dABS(X-XA(I))
        IF (H.EQ.0.d0)THEN
          Y=YA(I)
          DY=0.d0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          if(dd.eq.0.d0) then
            do ii=1,N
              write(*,*) ii,'x,y',xa(ii),ya(ii),x
            enddo
            y = ya(N)
            if(dabs((ya(1)-ya(N))/ya(N)).gt.1.d-10) then
              write(*,*) 'x,y**',x,y
              write(*,*) '****'
              write(*,*) 'FAILURE IN FUNCTION RATINT'
            endif
            return

          endif
          IF(DD.EQ.0.d0)write(*,*) 'FAILURE IN FUNCTION RATINT'
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      FUNCTION PLGNDR(L,MBAR,X)
C
C     THIS PROGRAM CALCULATES THE LEGENDRE POLYNOMIAL.
C     WHEN USING PLGNDR() ONE HAS TO MAKE SURE THAT M>=0 AND -L<=X<=+1
c     Now the definition of Condon&Shortley is used. This differs from
c     Abramowitz's definition by factor (-1)**mbar
c
      IMPLICIT REAL*8(A-H,O-Z)
      IF(MBAR.LT.0.OR.MBAR.GT.L.OR.ABS(X).GT.1.D0) THEN
        WRITE(*,*) 'BAD ARGUMENTS TO PLGNDR'
      ENDIF
      PMM=1.D0
      IF(MBAR.GT.0) THEN
        SOMX2=DSQRT((1.D0-X)*(1.D0+X))
        FACT=1.D0
        DO 151 I=1,MBAR
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.D0
 151    CONTINUE
      ENDIF
      IF(L.EQ.MBAR) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*MBAR+1)*PMM
        IF(L.EQ.MBAR+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 152 LL=MBAR+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+MBAR-1)*PMM)/(LL-MBAR)
            PMM=PMMP1
            PMMP1=PLL
 152      CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF

      PLGNDR = (-1)**MBAR* PLGNDR
      RETURN
      END

      function coeffMakerExplicit(d_m1,d_m2,k1,k2,ll)
      implicit real*8 (a-h,o-z)
      real*8 j_1,j_2,l_1,l_2,coeffMakerExplicit
      call factor()

      call DiracAngularl(dble(k1),l_1)
      j_1=dabs(dble(k1))-0.5d0
      call DiracAngularl(dble(k2),l_2)
      j_2=dabs(dble(k2))-0.5d0

      BigL=dble(ll)
      coeffMakerExplicit=Clebsh(l_1,d_m1-0.5d0,0.5d0,0.5d0,j_1,d_m1)*
     &Clebsh(l_2,d_m2-0.5d0,0.5d0,0.5d0,j_2,d_m2)*
     &Clebsh(l_1,d_m1-0.5d0,BigL,0.d0,l_2,d_m2-0.5d0)
      coeffMakerExplicit=coeffMakerExplicit+
     &Clebsh(l_1,d_m1+0.5d0,0.5d0,-0.5d0,j_1,d_m1)*
     &Clebsh(l_2,d_m2+0.5d0,0.5d0,-0.5d0,j_2,d_m2)*
     &Clebsh(l_1,d_m1+0.5d0,BigL,0.d0,l_2,d_m2+0.5d0)
      coeffMakerExplicit=coeffMakerExplicit*sqrt(2.d0*l_1+1.d0)/
     &sqrt(2.d0*l_2+1.d0)
      coeffMakerExplicit=coeffMakerExplicit*
     &Clebsh(l_1,0.d0,BigL,0.d0,l_2,0.d0)
      if(dabs(coeffMakerExplicit).lt. 1.d-12)coeffMakerExplicit=0.d0
      return
      end function coeffMakerExplicit
C****************************************************************
C*************Added on the 110810, necessary for the angular part of the
C  matrix element shown below

C  This program calculates the matrix element <X_k,m | P_L | X_K,M >
      function Coeff_maker(mu,k1,k2,ll)
      implicit real*8 (a-h,o-z)
      integer ll
      real*8 k1,k2,SumSubject,mu,BigL
      real*8 j_1,j_2,l_1,l_2,Coeff_maker
      call factor()

      call DiracAngularl(k1,l_1)
      j_1=dabs(k1)-0.5d0
      call DiracAngularl(k2,l_2)
      j_2=dabs(k2)-0.5d0

      BigL=ll*1.d0
c     do mm = -1,1,2
c     m=mm*0.5
      Coeff_maker=SumSubject(-0.5d0,mu,l_1,BigL,l_2,j_1,j_2)+
     &SumSubject(0.5d0,mu,l_1,BigL,l_2,j_1,j_2)
c     enddo
      if(dabs(coeff_maker).lt. 1.d-12)coeff_maker=0.d0
      return
      end function Coeff_maker

      function SumSubject(m,mu,l_1,BigL,l_2,j_1,j_2)
      real*8 SumSubject,m,mu,j_1,j_2,l_1,BigL,l_2
      real*8 Clebsh, symb3j
      SumSubject= ((-1.d0)**(mu-m))*dsqrt(2.d0*l_1+1.d0)*
     &dsqrt(2.d0*l_2+1.d0)*
     &Clebsh(l_2,mu-m,0.5d0,m,j_2,mu)*
     &Clebsh(l_1,mu-m,0.5d0,m,j_1,mu)*
     &symb3j(l_1,BigL,l_2,0.d0,0.d0,0.d0)*
     &symb3j(l_1,BigL,l_2,m-mu,0.d0,mu-m)
      return
      end function SumSubject

C  This program calculates the matrix element
c     <X_k,m | sigma_x * Y_L | X_K,M >
      function Coeff_maker2(mu,k1,k2,ll)
      implicit real*8 (a-h,o-z)
      integer ll
      real*8 k1,k2,SumSubject2,mu,BigL
      real*8 j_1,j_2,l_1,l_2,Coeff_maker2
      call factor()

      call DiracAngularl(k1,l_1)
      j_1=dabs(k1)-0.5d0
      call DiracAngularl(k2,l_2)
      j_2=dabs(k2)-0.5d0

      BigL=ll*1.d0
c     do mm = -1,1,2
c     m=mm*0.5
      Coeff_maker2=SumSubject2(-0.5d0,mu,l_1,BigL,l_2,j_1,j_2)+
     &SumSubject2(0.5d0,mu,l_1,BigL,l_2,j_1,j_2)
c     enddo
      if(dabs(coeff_maker2).lt. 1.d-12)coeff_maker2=0.d0
      return
      end function Coeff_maker2

C     When mu2=mu1+1
      function Coeff_maker2plus(mu,k1,k2,ll)
      implicit real*8 (a-h,o-z)
      integer ll
      real*8 k1,k2,SumSubject2,mu,BigL
      real*8 j_1,j_2,l_1,l_2
      call factor()

      call DiracAngularl(k1,l_1)
      j_1=dabs(k1)-0.5d0
      call DiracAngularl(k2,l_2)
      j_2=dabs(k2)-0.5d0

        BigL=ll*1.d0
c     do mm = -1,1,2
c     m=mm*0.5
      Coeff_maker2plus=
     &SumSubject2(0.5d0,mu,l_1,BigL,l_2,j_1,j_2)
c     enddo

      return
      end function Coeff_maker2plus

C     When mu2=mu1-1
      function Coeff_maker2minus(mu,k1,k2,ll)
      implicit real*8 (a-h,o-z)
      integer ll
      real*8 k1,k2,SumSubject2,mu,BigL
      real*8 j_1,j_2,l_1,l_2,Coeff_maker2minus
      call factor()

      call DiracAngularl(k1,l_1)
      j_1=dabs(k1)-0.5d0
      call DiracAngularl(k2,l_2)
      j_2=dabs(k2)-0.5d0

      BigL=ll*1.d0
c     do mm = -1,1,2
c     m=mm*0.5
      Coeff_maker2minus=
     &SumSubject2(-0.5d0,mu,l_1,BigL,l_2,j_1,j_2)
c    enddo

      return
      end function Coeff_maker2minus

      function SumSubject2(m,mu,l_1,BigL,l_2,j_1,j_2)
      include 'inc.par'
      real*8 SumSubject2,m,mu,j_1,j_2,l_1,BigL,l_2
      real*8 Clebsh, symb3j
      SumSubject2=((-1.d0)**(mu+m))*dsqrt(2.d0*l_1+1.d0)*
     &dsqrt(2.d0*l_2+1.d0)*dsqrt(2.d0*BigL+1.d0)*dsqrt(1.d0/4.d0/pi)*
     &Clebsh(l_2,mu+m,0.5d0,m,j_2,mu+2.d0*m)*
     &Clebsh(l_1,mu+m,0.5d0,-m,j_1,mu)*
c     &      *clebsh(bigl,0.d0,l_2,0.d0,l_1,0.d0)
c     &      *clebsh(l_2,mu+m,bigL,0.d0,l_1,mu+m)
     &symb3j(l_1,BigL,l_2,0.d0,0.d0,0.d0)*
     &symb3j(l_1,BigL,l_2,-m-mu,0.d0,mu+m)
      return
      end function SumSubject2

c function SumSubject(m,mu,l_1,BigL,l_2,j_1,j_2)
c real*8 SumSubject,m,mu,j_1,j_2,l_1,BigL,l_2
c real*8 Clebsh, symb3j
c SumSubject= ((-1)**(mu-m))*sqrt(2.d0*l_1+1.d0)*sqrt(2.d0*l_2+1.d0)
c     & *Clebsh(l_2,mu-m,0.5d0,m,j_1,mu)
c     & *Clebsh(l_1,mu-m,0.5d0,m,j_2,mu)
c     & *symb3j(l_1,BigL,l_2,0.d0,0.d0,0.d0)
c     & *symb3j(l_1,BigL,l_2,m-mu,0.d0,mu-m)
c      return
c end function SumSubject

      subroutine ang_sigz(ang,nkap,amjf,amji)
      include 'inc.par'
      real*8 ang(-nkap:nkap,0:2*nkap,-2*nkap:2*nkap,-nkap:nkap)

      ang=0.d0
      do k1=-nkap,nkap
        call DiracAngularl(1.d0*k1,al1)
        aj1=1.d0*abs(k1)-0.5d0
        do k2=-nkap,nkap
          if(k1*k2.ne.0)then
            call DiracAngularl(1.d0*k2,al2)
            aj2=1.d0*abs(k2)-0.5d0
            do ll=0,2*nkap
              do M=-ll,ll
                a_M=dble(M)
                a_ll=dble(ll)
                ang_stor=0.d0
                do L=max(0,ll-1),ll+1
                  a_L=dble(L)
                  ang_stor=ang_stor+
     &            clebsh(a_ll,a_M,1.d0,0.d0,a_L,a_M)*
     &            dsqrt((2.d0*aj2+1.d0)*(2.d0*a_ll+1.d0)*
     &            (2.d0*aj1+1.d0)*(2.d0*a_L+1.d0)*(2.d0*al1+1.d0)*
     &            2.d0)*
     &            symb9j(al2,0.5d0,aj2,al1,0.5d0,aj1,a_ll,1.d0,a_L)*
     &            clebsh(al1,0.d0,a_ll,0.d0,al2,0.d0)*
     &            symb3j(aj2,a_L,aj1,-amjf,a_M,amji)*
     &            ((-1)**idint(aj2-amji))
                enddo
                ang(k1,ll,M,k2)=dsqrt(3.d0/4.d0/pi)*ang_stor
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine ang_sigxx(ang,nkap,amu)
      include 'inc.par'
      real*8 ang(-nkap:nkap,0:2*nkap,-nkap:nkap)

      do k1=-nkap,nkap
        call DiracAngularl(1.d0*k1,al1)
        aj1=1.d0*abs(k1)-0.5d0
        do k2=-nkap,nkap
          if(k1*k2.ne.0)then
            call DiracAngularl(1.d0*k2,al2)
            aj2=1.d0*abs(k2)-0.5d0
            do ll=0,2*nkap
              a_ll=dble(ll)
              ang_stor=0.d0
              do L=1,ll+1
                a_L=dble(L)
                ang_stor=ang_stor+
     &          clebsh(a_ll,0.d0,1.d0,-1.d0,a_L,-1.d0)*
     &          dsqrt((2.d0*aj2+1.d0)*(2.d0*a_ll+1.d0)*
     &          (2.d0*al2+1.d0)*(2.d0*a_L+1.d0)*
     &          (6.d0)/4.d0/pi)*
     &          symb9j(al1,0.5d0,aj1,al2,0.5d0,aj2,a_ll,1.d0,a_L)*
     &          clebsh(al2,0.d0,a_ll,0.d0,al1,0.d0)*
     &          (clebsh(aj2,amu+1.d0,a_L,-1.d0,aj1,amu)-
     &          ((-1)**(ll+1-L))*clebsh(aj2,amu-1.d0,a_L,1.d0,aj1,amu))
              enddo
              ang(k1,ll,k2)=ang_stor/dsqrt(2.d0)
            enddo
          endif
        enddo
      enddo

      return
      end

C     When mu2=mu1+1
      subroutine ang_sigxx_plus(ang,nkap,amu)
      include 'inc.par'
      real*8 ang(-nkap:nkap,0:2*nkap,-nkap:nkap)

      do k1=-nkap,nkap
        call DiracAngularl(1.d0*k1,al1)
        aj1=1.d0*abs(k1)-0.5d0
        do k2=-nkap,nkap
          if(k1*k2.ne.0)then
            call DiracAngularl(1.d0*k2,al2)
            aj2=1.d0*abs(k2)-0.5d0
            do ll=0,2*nkap
              a_ll=dble(ll)
              ang_stor=0.d0
              do L=1,ll+1
                a_L=dble(L)
                ang_stor=ang_stor+
     &          clebsh(a_ll,0.d0,1.d0,-1.d0,a_L,-1.d0)*
     &          dsqrt((2.d0*aj2+1.d0)*(2.d0*a_ll+1.d0)*
     &          (2.d0*al2+1.d0)*(2.d0*a_L+1.d0)*
     &          (6.d0)/4.d0/pi)*
     &          symb9j(al1,0.5d0,aj1,al2,0.5d0,aj2,a_ll,1.d0,a_L)*
     &          clebsh(al2,0.d0,a_ll,0.d0,al1,0.d0)*
     &          (clebsh(aj2,amu+1.d0,a_L,-1.d0,aj1,amu))!-
c     &         ((-1)**(ll+1-L))*clebsh(aj2,amu-1.d0,a_L,1.d0,aj1,amu))
              enddo
              ang(k1,ll,k2)=ang_stor/dsqrt(2.d0)
            enddo
          endif
        enddo
      enddo

      return
      end

C     When mu2=mu1-1
      subroutine ang_sigxx_minus(ang,nkap,amu)
      include 'inc.par'
      real*8 ang(-nkap:nkap,0:2*nkap,-nkap:nkap)

      do k1=-nkap,nkap
        call DiracAngularl(1.d0*k1,al1)
        aj1=1.d0*abs(k1)-0.5d0
        do k2=-nkap,nkap
          if(k1*k2.ne.0)then
            call DiracAngularl(1.d0*k2,al2)
            aj2=1.d0*abs(k2)-0.5d0
            do ll=0,2*nkap
              a_ll=dble(ll)
              ang_stor=0.d0
              do L=1,ll+1
                a_L=dble(L)
                ang_stor=ang_stor+
     &          clebsh(a_ll,0.d0,1.d0,-1.d0,a_L,-1.d0)*
     &          dsqrt((2.d0*aj2+1.d0)*(2.d0*a_ll+1.d0)*
     &          (2.d0*al2+1.d0)*(2.d0*a_L+1.d0)*
     &          (6.d0)/4.d0/pi)*
     &          symb9j(al1,0.5d0,aj1,al2,0.5d0,aj2,a_ll,1.d0,a_L)*
     &          clebsh(al2,0.d0,a_ll,0.d0,al1,0.d0)*
     &          (-((-1)**(ll+1-L))*
     &          clebsh(aj2,amu-1.d0,a_L,1.d0,aj1,amu))
              enddo
              ang(k1,ll,k2)=ang_stor/dsqrt(2.d0)
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine DiracAngularL_int(k,L)
      integer k,L

      if(k.gt.0)then
        L=k
      else if(k.lt.0)then
        L=-k-1
      else
        L=0
      endif
      return
      end

      subroutine DiracAngularL(kppa,L)
      real*8 kppa, L
c One enters the value of kappa firstly to determine
c the value of orbital angular momentum l
c write(*,*) 'Enter the value of kappa'
c read(*,*) kappa
c Definition of L from kappa
      if (kppa .gt. 0) then
        L=kppa
      else if (kppa .lt. 0) then
        L=-1.d0*kppa -1.d0
      endif
c write(*,*) L
      return
      end

      subroutine DiracAngularJ(kppa,J)
      real*8 kppa
      real*8 J
c One enters the value of kappa firstly to determine
c the value of angular momentum l
c write(*,*) 'Enter the value of kappa'
c read(*,*) kappa
c Definition of L from kappa
      J=dabs(kppa)*1.d0-0.5d0
c write(*,*) L
      return
      end

c********************************************************
c********************************************************
c Added on the 110830 for the rutherford trajectories and the calculation
c of the internuclear distances as a function of xi from Tupitsyn.
c The units here are Nuclear units.

      subroutine New_Elapsed_time(z1,z2,m1,m2,v,b,xi,time)
      include 'inc.par'
      real*8 z1,z2,m1,m2,v,b,xi,time,a_faktor,e_faktor

      a_faktor=z1*z2/((v**2*c)*
     &(m1*m2*pe_massratio/(m1+m2)))
      e_faktor=dsqrt(1+(b/c)**2/(a_faktor**2))
      time=(e_faktor*dsinh(xi)+ xi)*a_faktor/v
      return
      end

      subroutine NewInterNuclDist(z1,z2,m1,m2,v,b,xi,IND)
      include 'inc.par'
      real*8 z1,z2,m1,m2,v,b,xi,IND,a_faktor,e_faktor
      if(z1*z2.ne. 0.d0)then
        a_faktor=z1*z2/((v**2*c)*
     &  (m1*m2*pe_massratio/(m1+m2)))
        e_faktor=dsqrt(1+(b/c)**2/(a_faktor**2))
        IND=a_faktor*(e_faktor*dcosh(xi)+1.d0)

c IND=(sqrt(1.d0+(b/137.035999d0)**2/(z1*z2/((v**2*137.035999d0**2)*
c     & (m1*m2*1836.15267245d0/(m1+m2))))**2)*cosh(xi)+1.d0)*
c     & z1*z2/((v**2*137.035999d0**2)*(m1*m2*1836.15267245d0/(m1+m2)))
c     IND=IND*137.035999d0
      else
        IND=1.d-4
      endif
      return
      end

c subroutine NewInterNuclDist(z1,z2,m1,m2,v,b,xi,IND)
c real*8 z1,z2,m1,m2,v,b,xi,IND
c IND=(sqrt(1.d0+b**2/(z1*(z2/137.035999d0)/((v**2)*
c    & (m1*m2*1836.15267245d0/(m1+m2))))**2)*cosh(xi)+1.d0)*
c    & z1*(z2/137.035999d0)/((v**2)*(m1*m2*1836.15267245d0/(m1+m2)))
c return
c end

      subroutine DofRwrtXi(z1,z2,m1,m2,v,b,xi,dRdXi)
      include 'inc.par'
      real*8 z1,z2,m1,m2,v,b,xi,dRdXi,a_faktor,e_faktor
      if(z1*z2.ne. 0.d0)then
        a_faktor=z1*z2/((v**2*c)*
     &  (m1*m2*pe_massratio/(m1+m2)))
        e_faktor=dsqrt(1+(b/c)**2/(a_faktor**2))
        dRdXi=a_faktor*e_faktor*dsinh(xi)

c dRdXi=(sqrt(1.d0+(b/137.035999d0)**2/(z1*(z2)/((v**2*137.035999d0**2)*
c     & (m1*m2*1836.15267245d0/(m1+m2))))**2)*sinh(xi))*
c     & z1*z2/((v**2*137.035999d0)*(m1*m2*1836.15267245d0/(m1+m2)))
      else
        dRdXi=0.d0
      endif
      return
      end

c subroutine DofRwrtXi(z1,z2,m1,m2,v,b,xi,dRdXi)
c real*8 z1,z2,m1,m2,v,b,xi,dRdXi
c dRdXi=(sqrt(1.d0+b**2/(z1*(z2/137.035999d0)/((v**2)*
c    & (m1*m2*1836.15267245d0/(m1+m2))))**2)*sinh(xi))*
c    & z1*(z2/137.035999d0)/((v**2)*(m1*m2*1836.15267245d0/(m1+m2)))
c return
c end

      subroutine DofTwrtXi(z1,z2,m1,m2,v,b,xi,dTdXi)
      include 'inc.par'
      real*8 z1,z2,m1,m2,v,b,xi,dTdXi,a_faktor,e_faktor
      if(z1*z2.ne. 0.d0)then
        a_faktor=z1*z2/((v**2*c)*
     &  (m1*m2*pe_massratio/(m1+m2)))
        e_faktor=dsqrt(1+(b/c)**2/(a_faktor**2))
        dTdXi=a_faktor*(e_faktor*dcosh(xi)+1)/v

c dTdXi=(sqrt(1.d0+(b/137.035999d0)**2/(z1*z2/((v**2*137.035999d0**2)*
c     & (m1*m2*1836.15267245d0/(m1+m2))))**2)*cosh(xi)+1.d0)*
c     & z1*z2/((v**2*137.035999d0**2)*(m1*m2*1836.15267245d0/(m1+m2)))
c      dTdXi=137.035999d0*dTdXi/v
      else
        dTdXi=0.d0
      endif
      return
      end

      subroutine Axis_angle_Theta(z1,z2,m1,m2,v,b,xi,Theta)
      include 'inc.par'
      real*8 z1,z2,m1,m2,v,b,xi,Theta,a_faktor
      if(z1*z2.ne. 0.d0)then

        a_faktor=z1*z2/((v**2*c)*
     &  (m1*m2*pe_massratio/(m1+m2)))

        epslon=dsqrt(1.d0+(b/c)**2/a_faktor**2)

        Theta=2.d0*datan(dsqrt(epslon**2 - 1.d0)*(dtanh(xi/2.d0)+1.d0)/
     &  ((epslon+1.d0)-(epslon-1.d0)*dtanh(xi/2.d0)))


      else
        Theta=0.d0
      endif
      return
      end

      subroutine DofThetawrtXi(z1,z2,m1,m2,v,b,xi,dThetadXi)
      include 'inc.par'
      real*8 z1,z2,m1,m2,v,b,xi,dThetadXi
      if(z1*z2.ne. 0.d0)then

        a_faktor=z1*z2/((v**2*c)*
     &  (m1*m2*pe_massratio/(m1+m2)))
        epslon=dsqrt(1.d0+(b/c)**2/a_faktor**2)

        dThetadXi=dsqrt(-1.d0 + epslon**2)/(1.d0 + epslon*dcosh(xi))

      else
        dThetadXi=0.d0
      endif
      return
      end

c subroutine DofTwrtXi(z1,z2,m1,m2,v,b,xi,dTdXi)
c real*8 z1,z2,m1,m2,v,b,xi,dTdXi
c dTdXi=(sqrt(1.d0+b**2/(z1*(z2/137.035999d0)/((v**2)*
c    & (m1*m2*1836.15267245d0/(m1+m2))))**2)*cosh(xi)+1.d0)*
c     & z1*(z2/137.035999d0)/((v**2)*(m1*m2*1836.15267245d0/(m1+m2)))
c      dTdXi=dTdXi/v
c return
c end

      subroutine MatDimAndStateSelection_eqCharge(nkap,nm,
     &number_states,eigval1,eigval2,nsto,nste,
     &state_range_input,Start_state_e,Start_state_o,
     &Lowest_bound_e,Lowest_bound_o,End_state_e,End_state_o,
     &Mat_Dimension,Highest_bound_e,Highest_bound_o,iixi)
      include 'inc.par'
      logical state_range_input
      real*8 eigval1(nsto),eigval2(nste)
      integer nkap,number_states(2,-nkap:nkap),
     &Highest_bound_o,Highest_bound_e, xi_stepslower,xi_stepsupper,
     &Lowest_bound_e,Start_state_e,End_state_e,Mat_Dimension,iixi,
     &Lowest_bound_o,Start_state_o,End_state_o
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /momentum_projection/ amu,amj_max

      if(state_range_input)then
        open(1,file='inp.inp',status='old')
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) Manual_bound_states
        read(1,*) Manual_pcont_states
        read(1,*) Manual_ncont_states
        close(1)
      else
        open(1,file='inp.inp',status='old')
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) xi_stepslower,xi_stepsupper
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) Manual_ncont_states
        close(1)
        if(iixi.ne.xi_stepslower ) then
          open(56,file='Auto_state_allocation.dat')
          read(56,*) Manual_bound_states_e
          read(56,*) Manual_bound_states_o
          read(56,*) Manual_pcont_states_e
          read(56,*) Manual_pcont_states_o
          read(56,*) Manual_ncont_states_e
          read(56,*) Manual_ncont_states_o
          close(56)
        endif
      endif

      if((nkap/2)*2.eq.nkap)then
        nmin=-nkap+1
      else
        nmin=-nkap
      endif

      Lowest_bound_e=1
      Lowest_bound_o=1
ccc    if(z_nuc1.ne.z_nuc2)then
ccc    do kk=-nkap,nkap
ccc       if(kk.ne. 0) then
ccc          Lowest_bound=Lowest_bound+(nm+1-number_states(1,kk))
ccc       endif
ccc    enddo
ccc    else

c     if amu>=1.5, we need to reconfigure where the lowest_bound is
c     Since no bound state of 2p3/2 will ever be at risk of diving
c     into the negative continuum, we can simply search for the lowest
c     energy above -1.
      if(dabs(amu).le. 0.5d0 .and. z_nuc1+z_nuc2.gt. 174.d0)then
        do kk=nmin,-1,2
          Lowest_bound_e=Lowest_bound_e+(nm+1-number_states(1,kk))
        enddo
        do kk=2,nkap,2
          Lowest_bound_e=Lowest_bound_e+(nm+1-number_states(1,kk))
        enddo
        do kk=-nkap,-1
          Lowest_bound_o=Lowest_bound_o+(nm+1-number_states(1,kk))
        enddo
        do kk=1,nkap
          Lowest_bound_o=Lowest_bound_o+(nm+1-number_states(1,kk))
        enddo
        Lowest_bound_o=Lowest_bound_o-Lowest_bound_e+1
      else
        do while(eigval1(Lowest_bound_o).lt. -1.d0)
          Lowest_bound_o=Lowest_bound_o+1
        enddo
        do while(eigval2(Lowest_bound_e).lt. -1.d0)
          Lowest_bound_e=Lowest_bound_e+1
        enddo
      endif
ccc    endif

ccc    if(z_nuc1.ne.z_nuc2)then
ccc    i=nstates
ccc    do while(eigval_a(i).gt. 1.d0 )
ccc    i=i-1
ccc    enddo
ccc    Highest_bound=i
ccc    else
      i=nsto
      do while(eigval1(i).gt. 1.d0 )
        i=i-1
      enddo
      Highest_bound_o=i
      i=nste
      do while(eigval2(i).gt. 1.d0 )
        i=i-1
      enddo
      Highest_bound_e=i
ccc    endif

      if(state_range_input) then
        Start_state_e=Lowest_bound_e-Manual_ncont_states
        End_state_e=Highest_bound_e+
     &  Manual_pcont_states
        End_state_e=min(End_state_e,nste)

        Start_state_o=Lowest_bound_o-Manual_ncont_states
        End_state_o=Highest_bound_o+
     &  Manual_pcont_states
        End_state_o=min(End_state_o,nsto)

        Mat_Dimension=Manual_bound_states+
     &  min(Manual_pcont_states,max(0,End_state_o-
     &  Highest_bound_o,End_state_e-
     &  Highest_bound_e))+
     &  Manual_ncont_states
      elseif(.not.state_range_input .and.
     &iixi.ne.xi_stepslower)then
        Start_state_e=Lowest_bound_e-Manual_ncont_states_e
        End_state_e=Start_state_e+Manual_ncont_states_e+
     &  Manual_pcont_states_e+Manual_bound_states_e-1

        Start_state_o=Lowest_bound_o-Manual_ncont_states_o
        End_state_o=Start_state_o+Manual_ncont_states_o+
     &  Manual_pcont_states_o+Manual_bound_states_o-1

        Mat_Dimension=min(
     &  Manual_bound_states_e+Manual_pcont_states_e+
     &  Manual_ncont_states_e,
     &  Manual_bound_states_o+Manual_pcont_states_o+
     &  Manual_ncont_states_o)

      elseif(.not.state_range_input .and.
     &iixi.eq.xi_stepslower)then
        if(Manual_ncont_states.eq. 1) then
          Start_state_e=1
          Start_state_o=1
c            if(up_energy.gt. 3.d0)then
c       End_state=min(nste- 6*2*nkap,nsto- 6*2*nkap)
c       Mat_Dimension=min(nste- 6*2*nkap,nsto- 6*2*nkap)
c            else
          End_state_e=nste
          End_state_o=nsto
          Mat_Dimension=min(nste,nsto)
c            endif
        else
          Start_state_e=Lowest_bound_e
          Start_state_o=Lowest_bound_o
          End_state_e=nste
          End_state_o=nsto
          Mat_Dimension=min(nste-Start_state_e+1,nsto-Start_state_o+1)
        endif

      endif
c    pause
      write(*,*)'DIMENSION',Mat_Dimension,
     &'START STATE EVEN',Start_state_e,
     &'START STATE ODD',Start_state_o,
     &'END STATE EVEN',End_state_e,
     &'END STATE ODD',End_state_o,
     &'LOWEST BOUND EVEN', Lowest_bound_e,
     &'LOWEST BOUND ODD', Lowest_bound_o

      return
      end

      subroutine MatDimAndStateSelection_NeqCharge(nkap,nm,
     &number_states,nstates,eigval,
     &state_range_input,Start_state,Lowest_bound,End_state,
     &Mat_Dimension,Highest_bound,iixi)
      include 'inc.par'
      logical state_range_input
      real*8 eigval(nstates)
      integer nkap,number_states(2,-nkap:nkap),
     &Highest_bound,xi_stepslower,xi_stepsupper,
     &Manual_ncont_states,Manual_pcont_states,Manual_bound_states,
     &Lowest_bound,Start_state,End_state,Mat_Dimension,iixi
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      common /momentum_projection/ amu,amj_max
      common /common_dkb/ dkb
      logical dkb

      if(state_range_input)then
        open(1,file='inp.inp',status='old')
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) Manual_bound_states
        read(1,*) Manual_pcont_states
        read(1,*) Manual_ncont_states
        close(1)
      else
        open(1,file='inp.inp',status='old')
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) xi_stepslower,xi_stepsupper
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) Manual_ncont_states
        close(1)
        if(iixi.ne.xi_stepslower ) then
        open(56,file='Auto_state_allocation.dat')
        read(56,*) Manual_bound_states
        read(56,*) Manual_pcont_states
        read(56,*) Manual_ncont_states
        close(56)
        endif
      endif

c     if amu>=1.5, we need to reconfigure where the lowest_bound is
c     Since no bound state of 2p3/2 will ever be at risk of diving
c     into the negative continuum, we can simply search for the lowest
c     energy above -1.
      Lowest_bound=1
      if(dabs(amu).le. 0.5d0 .and.
     &z_nuc1+z_nuc2.gt. 174.d0)then
        do kk=-nkap,nkap
          if(kk.ne. 0) then
            Lowest_bound=Lowest_bound+(nm+1-number_states(1,kk))
          endif
        enddo
      else
        do while( eigval(Lowest_bound).lt. -1.d0)
          Lowest_bound=Lowest_bound+1
        enddo
      endif

      i=nstates
      do while(eigval(i).gt. 1.d0 )
        i=i-1
      enddo
      Highest_bound=i


      if(state_range_input) then
        Start_state=Lowest_bound-Manual_ncont_states
        End_state=min(Highest_bound+Manual_pcont_states,nstates)
        Mat_Dimension=Manual_bound_states+
     &  min(Manual_pcont_states,max(0,End_state-Highest_bound))+
     &  Manual_ncont_states
      elseif(.not.state_range_input .and. iixi.ne.
     &xi_stepslower) then
        Start_state=Lowest_bound-Manual_ncont_states
        End_state=Start_state+
     &  Manual_pcont_states+Manual_bound_states-1
        Mat_Dimension=End_state-Start_state+1
      elseif(.not.state_range_input .and.
     &iixi.eq.xi_stepslower)then
        if(Manual_ncont_states.eq. 1) then
          Start_state=1
c            if(up_energy.gt. 3.d0)then
c       End_state=nstates - 6*2*nkap
c       Mat_Dimension=nstates - 6*2*nkap
c            else
          End_state=nstates
          Mat_Dimension=nstates
c            endif
        endif
      endif
      write(*,*)'DIMENSION',Mat_Dimension,
     &'START STATE',Start_state,
     &'END STATE',End_state,
     &'LOWEST BOUND', Lowest_bound
      return
      end


c subroutine MatrixMM_Storage_DKB(wave_new,nstates,nm,nkap,
c     &  dvdRmatdkb1,dvdRmatdkb2,state_in_cont,
c     &  nvmat,wavesum_l1,Start_state,Mat_dimension,
c     &  wavesum_l2,wavesum_s1,wavesum_s2,vmat)
c include 'inc.par'
c integer Start_state,Mat_dimension,nm,nkap
c real*8 wave_new(nstates,2*nm,-nkap:nkap),
c     &  wavesum_l1(Start_state:Start_state+Mat_dimension-1,
c     &  nm,0:2*nkap,-nkap:nkap),
c     &  wavesum_l2(Start_state:Start_state+Mat_dimension-1,
c     &  nm,0:2*nkap,-nkap:nkap),
c     &  wavesum_s1(Start_state:Start_state+Mat_dimension-1,
c     &  nm,0:2*nkap,-nkap:nkap),
c     &  wavesum_s2(Start_state:Start_state+Mat_dimension-1,
c     &  nm,0:2*nkap,-nkap:nkap),
c     & dvdRmatdkb1(nm,nm,0:2*nkap,-nkap:nkap,2),
c     &  dvdRmatdkb2(nm,nm,0:2*nkap,-nkap:nkap,-nkap:nkap,2),
c     &  MatEl_Ang_comp_lge,MatEl_Ang_comp_sml,
c     &  vmat(nm,nm,0:2*nkap,nvmat)
c integer nstates,state_2,kappa_1,kappa_2,L,nvmat,i,
c     &  state_in_cont(Start_state:Start_state+Mat_dimension-1)
c logical state_range_input
c common /staterangeinput/ state_range_input
c common /momentum_projection/ amu,amj_max

c wavesum_l1=0.d0
c wavesum_l2=0.d0
c wavesum_s1=0.d0
c wavesum_s2=0.d0

c do jj=Start_state,Start_state+Mat_dimension-1!End_state
c       state_2=state_in_cont(jj)
c    do kappa_1=-nkap,nkap
c       do kappa_2=-nkap,nkap
c       if(kappa_1*kappa_2.ne. 0)then
c          do L=0,2*nkap
cc   MatEl_Ang_comp_lge=
cc     &   Coeff_maker(amu,1.d0*kappa_1,1.d0*kappa_2,L)
cc   if(MatEl_Ang_comp_lge.gt. 1.d-12) then
c      do i=1,nm
c                do j1=max(1,i-ns+1),min(i+ns-1,nm)
c
c wavesum_l1(jj,i,L,kappa_1)=
c     &  wavesum_l1(jj,i,L,kappa_1)+
c     &  vmat(i,j1,L,nvmat)*
c     &  wave_new(state_2,j1,kappa_2)+
c     &  dvdRmatdkb1(j1,i,L,kappa_2,1)*
c     &  wave_new(state_2,j1+nm,kappa_2)

c wavesum_l2(jj,i,L,kappa_1)=
c     &  wavesum_l2(jj,i,L,kappa_1)+
c     &  dvdRmatdkb1(i,j1,L,kappa_1,1)*
c     &  wave_new(state_2,j1,kappa_2)+
c     &  dvdRmatdkb2(i,j1,L,kappa_1,kappa_2,1)*
c     &  wave_new(state_2,j1+nm,kappa_2)

c         enddo
c      enddo
cc          endif
c          enddo
c       endif
c       enddo
c    enddo
c enddo

c do jj=Start_state,Start_state+Mat_dimension-1!End_state
c       state_2=state_in_cont(jj)
c    do kappa_1=-nkap,nkap
c       do kappa_2=-nkap,nkap
c       if(kappa_1*kappa_2.ne. 0)then
c          do L=0,2*nkap
cc   MatEl_Ang_comp_sml=
cc     &   Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
cc   if(MatEl_Ang_comp_sml.gt. 1.d-12) then
c      do i=1,nm
c                do j1=max(1,i-ns+1),min(i+ns-1,nm)

c wavesum_s1(jj,i,L,kappa_1)=
c     &  wavesum_s1(jj,i,L,kappa_1)+
c     &  vmat(i,j1,L,nvmat)*
c     &  wave_new(state_2,j1+nm,kappa_2)+
c     &  dvdRmatdkb1(j1,i,L,kappa_2,2)*
c     &  wave_new(state_2,j1,kappa_2)

c wavesum_s2(jj,i,L,kappa_1)=
c     & wavesum_s2(jj,i,L,kappa_1)+
c     &  dvdRmatdkb1(i,j1,L,kappa_1,2)*
c     &  wave_new(state_2,j1+nm,kappa_2)+
c     &  dvdRmatdkb2(i,j1,L,kappa_1,kappa_2,2)*
c     &  wave_new(state_2,j1,kappa_2)

c         enddo
c      enddo
cc             endif
c          enddo
c          endif
c       enddo
c    enddo
c enddo

c return
c end

      subroutine MatrixMM_Storage_DKB_1(wave_new,nm,nkap,nstates,lstep,
     &dvdRmatdkb1,nvmat,wavesum_l1,wavesum_s1,vmat,k1,i_state)
      include 'inc.par'
      real*8 wave_new(nstates,2*nm,-nkap:nkap),wavesum_l1(nm,0:2*nkap),
     &wavesum_s1(nm,0:2*nkap),dvdRmatdkb1(nm,nm,-nkap:nkap,0:2*nkap,2),
     &vmat(nm,nm,0:2*nkap,nvmat)

      wavesum_l1=0.d0
      wavesum_s1=0.d0
      do ll=0,2*nkap,lstep
        do i=1,nm
          do j1=max(1,i-ns+1),min(i+ns-1,nm)
            wavesum_l1(i,ll)=wavesum_l1(i,ll)+
     &      (vmat(i,j1,ll,nvmat)*wave_new(i_state,j1,k1)+
     &      dvdRmatdkb1(j1,i,k1,ll,1)*wave_new(i_state,j1+nm,k1))
          enddo
        enddo
      enddo

      do ll=0,2*nkap,lstep
        do i=1,nm
          do j1=max(1,i-ns+1),min(i+ns-1,nm)
            wavesum_s1(i,ll)=wavesum_s1(i,ll)+
     &      (vmat(i,j1,ll,nvmat)*wave_new(i_state,j1+nm,k1)+
     &      dvdRmatdkb1(j1,i,k1,ll,2)*wave_new(i_state,j1,k1))
          enddo
        enddo
      enddo
      return
      end

      subroutine MatrixMM_Storage_DKB_2(wave_new,nm,nkap,nstates,lstep,
     &dvdRmatdkb2,dvdRmatdkb1,wavesum_l2,wavesum_s2,k1,i_state,k2,
     &d_m1,d_m2)
      include 'inc.par'
      real*8 wave_new(nstates,2*nm,-nkap:nkap),wavesum_l2(nm,0:2*nkap),
     &wavesum_s2(nm,0:2*nkap),dvdRmatdkb1(nm,nm,-nkap:nkap,0:2*nkap,2),
     &dvdRmatdkb2(nm,nm,-nkap:nkap,-nkap:nkap,0:2*nkap,2)

      wavesum_l2=0.d0
      wavesum_s2=0.d0

      do ll=0,2*nkap,lstep
!        ang_comp_l=coeff_maker(amu,1.d0*k1,1.d0*k2,ll)
        ang_comp_l=coeffMakerExplicit(d_m1,d_m2,k1,k2,ll)
        if(ang_comp_l.ne.0.d0)then
          do i=1,nm
            do j1=max(1,i-ns+1),min(i+ns-1,nm)
              wavesum_l2(i,ll)=wavesum_l2(i,ll)+
     &        (dvdRmatdkb1(i,j1,k2,ll,1)*wave_new(i_state,j1,k1)+
     &        dvdRmatdkb2(i,j1,k2,k1,ll,1)*
     &        wave_new(i_state,j1+nm,k1))
            enddo
          enddo
        endif
      enddo

      do ll=0,2*nkap,lstep
!        ang_comp_s=coeff_maker(amu,-1.d0*k1,-1.d0*k2,ll)
        ang_comp_s=coeffMakerExplicit(d_m1,d_m2,-k1,-k2,ll)
        if(ang_comp_s.ne.0.d0)then
          do i=1,nm
            do j1=max(1,i-ns+1),min(i+ns-1,nm)
              wavesum_s2(i,ll)=wavesum_s2(i,ll)+
     &        (dvdRmatdkb1(i,j1,k2,ll,2)*wave_new(i_state,j1+nm,k1)+
     &        dvdRmatdkb2(i,j1,k2,k1,ll,2)*wave_new(i_state,j1,k1))
            enddo
          enddo
        endif
      enddo

      return
      end

      subroutine MatrixMM_Storage_DKB_even(wave_new,nstates,nm,nkap,
     &dvdRmatdkb1,dvdRmatdkb2,state_in_cont,
     &nvmat,wavesum_l1,Start_state,Mat_dimension,
     &wavesum_l2,wavesum_s1,wavesum_s2,vmat)
      include 'inc.par'
      integer Start_state,Mat_dimension,nm,nkap
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_l1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_l2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &wavesum_s1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_s2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &dvdRmatdkb1(nm,nm,0:2*nkap,-nkap:nkap,2),
     &dvdRmatdkb2(nm,nm,0:2*nkap,-nkap:nkap,-nkap:nkap,2),
     &vmat(nm,nm,0:2*nkap,nvmat),
     &MatEl_Ang_comp_lge(-nkap:nkap,0:2*nkap,-nkap:nkap),
     &MatEl_Ang_comp_sml(-nkap:nkap,0:2*nkap,-nkap:nkap)
      integer nstates,state_2,L,nvmat,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
c      logical state_range_input
c      common /staterangeinput/ state_range_input
      common /momentum_projection/ amu,amj_max


      MatEl_Ang_comp_lge=0.d0
      MatEl_Ang_comp_sml=0.d0
      do k1=-nkap,nkap
        call DiracAngularL_int(k1,L_k1)
        do k2=-nkap,nkap
          call DiracAngularL_int(k2,L_k2)
          if((1+((-1)**(L_k1))).ne.0 .and. (1+((-1)**(L_k2))).ne.0
     &    .and. k1*k2.ne.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_lge(k1,L,k2)=
     &        Coeff_maker(amu,k1*1.d0,k2*1.d0,L)
              MatEl_Ang_comp_sml(k1,L,k2)=
     &        Coeff_maker(amu,-k1*1.d0,-k2*1.d0,L)
            enddo
          endif
        enddo
      enddo


      wavesum_l1=0.d0
      wavesum_l2=0.d0
      wavesum_s1=0.d0
      wavesum_s2=0.d0
c      do kk2=-nkap,nkap
      do kk=-nkap,nkap
        call DiracAngularL_int(kk,L_kk)
c      if(kk*kk2.ne.0)then
        if(kk.ne.0 .and. (1+((-1)**(L_kk))).ne.0)then
          do jj=Start_state,Start_state+Mat_dimension-1!End_state
            state_2=state_in_cont(jj)
            do L=0,2*nkap,2
c            if(dabs(MatEl_Ang_comp_lge(kk2,L,kk)).gt.1.d-14)then
              do i=1,nm
                do j1=max(1,i-ns+1),min(i+ns-1,nm)
c             wavesum_l1(jj,i,L,kk,kk2)=wavesum_l1(jj,i,L,kk,kk2)+
                  wavesum_l1(jj,i,L,kk)=wavesum_l1(jj,i,L,kk)+
     &            (vmat(i,j1,L,nvmat)*wave_new(state_2,j1,kk)+
     &            dvdRmatdkb1(j1,i,L,kk,1)*wave_new(state_2,j1+nm,kk))
                enddo
c                  wavesum_l1(jj,i,L,kk,kk2)=wavesum_l1(jj,i,L,kk,kk2)*
c     &            MatEl_Ang_comp_lge(kk2,L,kk)
              enddo
c            endif
            enddo
          enddo
        endif
      enddo
c      enddo

      do kk2=-nkap,nkap
        call DiracAngularL_int(kk2,L_kk2)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk*kk2.ne.0 .and. (1+((-1)**(L_kk))).ne.0 .and.
     &    (1+((-1)**(L_kk2))).ne.0)then
            do jj=Start_state,Start_state+Mat_dimension-1!End_state
              state_2=state_in_cont(jj)
              do L=0,2*nkap,2
                if(dabs(MatEl_Ang_comp_lge(kk2,L,kk)).gt.1.d-14)then
                  do i=1,nm
                    do j1=max(1,i-ns+1),min(i+ns-1,nm)
                      wavesum_l2(jj,i,L,kk,kk2)=
     &                wavesum_l2(jj,i,L,kk,kk2)+
     &                (dvdRmatdkb1(i,j1,L,kk2,1)*
     &                wave_new(state_2,j1,kk)+
     &                dvdRmatdkb2(i,j1,L,kk2,kk,1)*
     &                wave_new(state_2,j1+nm,kk))
                    enddo
c                wavesum_l2(jj,i,L,kk,kk2)=wavesum_l2(jj,i,L,kk,kk2)*
c     &          MatEl_Ang_comp_lge(kk2,L,kk)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

c      do kk2=-nkap,nkap
      do kk=-nkap,nkap
        call DiracAngularL_int(kk,L_kk)
        if(kk.ne.0 .and. (1+((-1)**(L_kk))).ne.0)then
          do jj=Start_state,Start_state+Mat_dimension-1!End_state
            state_2=state_in_cont(jj)
            do L=0,2*nkap,2
c            if(dabs(MatEl_Ang_comp_sml(kk2,L,kk)).gt.1.d-14)then
              do i=1,nm
                do j1=max(1,i-ns+1),min(i+ns-1,nm)
c             wavesum_s1(jj,i,L,kk,kk2)=wavesum_s1(jj,i,L,kk,kk2)+
                  wavesum_s1(jj,i,L,kk)=wavesum_s1(jj,i,L,kk)+
     &            (vmat(i,j1,L,nvmat)*wave_new(state_2,j1+nm,kk)+
     &            dvdRmatdkb1(j1,i,L,kk,2)*wave_new(state_2,j1,kk))
                enddo
c                  wavesum_s1(jj,i,L,kk,kk2)=wavesum_s1(jj,i,L,kk,kk2)*
c     &            MatEl_Ang_comp_sml(kk2,L,kk)
              enddo
c            endif
            enddo
          enddo
        endif
      enddo
c      enddo

      do kk2=-nkap,nkap
        call DiracAngularL_int(kk2,L_kk2)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk*kk2.ne.0 .and. (1+((-1)**(L_kk))).ne.0 .and.
     &    (1+((-1)**(L_kk2))).ne.0)then
            do jj=Start_state,Start_state+Mat_dimension-1!End_state
              state_2=state_in_cont(jj)
              do L=0,2*nkap,2
                if(dabs(MatEl_Ang_comp_sml(kk2,L,kk)).gt.1.d-14)then
                  do i=1,nm
                    do j1=max(1,i-ns+1),min(i+ns-1,nm)
                      wavesum_s2(jj,i,L,kk,kk2)=
     &                wavesum_s2(jj,i,L,kk,kk2)+
     &                (dvdRmatdkb1(i,j1,L,kk2,2)*
     &                wave_new(state_2,j1+nm,kk)+
     &                dvdRmatdkb2(i,j1,L,kk2,kk,2)*
     &                wave_new(state_2,j1,kk))
                    enddo
c                  wavesum_s2(jj,i,L,kk,kk2)=wavesum_s2(jj,i,L,kk,kk2)*
c     &            MatEl_Ang_comp_sml(kk2,L,kk)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine MatrixMM_Storage_DKB_odd(wave_new,nstates,nm,nkap,
     &dvdRmatdkb1,dvdRmatdkb2,state_in_cont,
     &nvmat,wavesum_l1,Start_state,Mat_dimension,
     &wavesum_l2,wavesum_s1,wavesum_s2,vmat)
      include 'inc.par'
      integer Start_state,Mat_dimension,nm,nkap
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_l1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_l2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &wavesum_s1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_s2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &dvdRmatdkb1(nm,nm,0:2*nkap,-nkap:nkap,2),
     &dvdRmatdkb2(nm,nm,0:2*nkap,-nkap:nkap,-nkap:nkap,2),
     &vmat(nm,nm,0:2*nkap,nvmat),
     &MatEl_Ang_comp_lge(-nkap:nkap,0:2*nkap,-nkap:nkap),
     &MatEl_Ang_comp_sml(-nkap:nkap,0:2*nkap,-nkap:nkap)
      integer nstates,state_2,L,nvmat,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
c      logical state_range_input
c      common /staterangeinput/ state_range_input
      common /momentum_projection/ amu,amj_max


      MatEl_Ang_comp_lge=0.d0
      MatEl_Ang_comp_sml=0.d0
      do k1=-nkap,nkap
        call DiracAngularL_int(k1,L_k1)
        do k2=-nkap,nkap
          call DiracAngularL_int(k2,L_k2)
          if((1+((-1)**(L_k1))).eq.0 .and. (1+((-1)**(L_k2))).eq.0
     &    .and. k2*k1.ne.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_lge(k1,L,k2)=
     &        Coeff_maker(amu,k1*1.d0,k2*1.d0,L)
              MatEl_Ang_comp_sml(k1,L,k2)=
     &        Coeff_maker(amu,-k1*1.d0,-k2*1.d0,L)
            enddo
          endif
        enddo
      enddo


      wavesum_l1=0.d0
      wavesum_l2=0.d0
      wavesum_s1=0.d0
      wavesum_s2=0.d0
c      do kk2=-nkap,nkap
      do kk=-nkap,nkap
        call DiracAngularL_int(kk,L_kk)
c      if(kk*kk2.ne.0)then
        if(kk.ne.0 .and. (1+((-1)**(L_kk))).eq.0)then
          do jj=Start_state,Start_state+Mat_dimension-1!End_state
            state_2=state_in_cont(jj)
            do L=0,2*nkap,2
c            if(dabs(MatEl_Ang_comp_lge(kk2,L,kk)).gt.1.d-14)then
              do i=1,nm
                do j1=max(1,i-ns+1),min(i+ns-1,nm)
c             wavesum_l1(jj,i,L,kk,kk2)=wavesum_l1(jj,i,L,kk,kk2)+
                  wavesum_l1(jj,i,L,kk)=wavesum_l1(jj,i,L,kk)+
     &            (vmat(i,j1,L,nvmat)*wave_new(state_2,j1,kk)+
     &            dvdRmatdkb1(j1,i,L,kk,1)*wave_new(state_2,j1+nm,kk))
                enddo
c                  wavesum_l1(jj,i,L,kk,kk2)=wavesum_l1(jj,i,L,kk,kk2)*
c     &            MatEl_Ang_comp_lge(kk2,L,kk)
              enddo
c            endif
            enddo
          enddo
        endif
      enddo
c      enddo

      do kk2=-nkap,nkap
        call DiracAngularL_int(kk2,L_kk2)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk*kk2.ne.0 .and. (1+((-1)**(L_kk))).eq.0 .and.
     &    (1+((-1)**(L_kk2))).eq.0)then
            do jj=Start_state,Start_state+Mat_dimension-1!End_state
              state_2=state_in_cont(jj)
              do L=0,2*nkap,2
                if(dabs(MatEl_Ang_comp_lge(kk2,L,kk)).gt.1.d-14)then
                  do i=1,nm
                    do j1=max(1,i-ns+1),min(i+ns-1,nm)
                      wavesum_l2(jj,i,L,kk,kk2)=
     &                wavesum_l2(jj,i,L,kk,kk2)+
     &                (dvdRmatdkb1(i,j1,L,kk2,1)*
     &                wave_new(state_2,j1,kk)+
     &                dvdRmatdkb2(i,j1,L,kk2,kk,1)*
     &                wave_new(state_2,j1+nm,kk))
                    enddo
c                  wavesum_l2(jj,i,L,kk,kk2)=wavesum_l2(jj,i,L,kk,kk2)*
c     &            MatEl_Ang_comp_lge(kk2,L,kk)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

c      do kk2=-nkap,nkap
      do kk=-nkap,nkap
        call DiracAngularL_int(kk,L_kk)
        if(kk.ne.0 .and. (1+((-1)**(L_kk))).eq.0)then
          do jj=Start_state,Start_state+Mat_dimension-1!End_state
            state_2=state_in_cont(jj)
            do L=0,2*nkap,2
c            if(dabs(MatEl_Ang_comp_sml(kk2,L,kk)).gt.1.d-14)then
              do i=1,nm
                do j1=max(1,i-ns+1),min(i+ns-1,nm)
c             wavesum_s1(jj,i,L,kk,kk2)=wavesum_s1(jj,i,L,kk,kk2)+
                  wavesum_s1(jj,i,L,kk)=wavesum_s1(jj,i,L,kk)+
     &            (vmat(i,j1,L,nvmat)*wave_new(state_2,j1+nm,kk)+
     &            dvdRmatdkb1(j1,i,L,kk,2)*wave_new(state_2,j1,kk))
                enddo
c                  wavesum_s1(jj,i,L,kk,kk2)=wavesum_s1(jj,i,L,kk,kk2)*
c     &            MatEl_Ang_comp_sml(kk2,L,kk)
              enddo
c            endif
            enddo
          enddo
        endif
      enddo
c      enddo

      do kk2=-nkap,nkap
        call DiracAngularL_int(kk2,L_kk2)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk*kk2.ne.0 .and. (1+((-1)**(L_kk))).eq.0 .and.
     &    (1+((-1)**(L_kk2))).eq.0)then
            do jj=Start_state,Start_state+Mat_dimension-1!End_state
              state_2=state_in_cont(jj)
              do L=0,2*nkap,2
                if(dabs(MatEl_Ang_comp_sml(kk2,L,kk)).gt.1.d-14)then
                  do i=1,nm
                    do j1=max(1,i-ns+1),min(i+ns-1,nm)
                      wavesum_s2(jj,i,L,kk,kk2)=
     &                wavesum_s2(jj,i,L,kk,kk2)+
     &                (dvdRmatdkb1(i,j1,L,kk2,2)*
     &                wave_new(state_2,j1+nm,kk)+
     &                dvdRmatdkb2(i,j1,L,kk2,kk,2)*
     &                wave_new(state_2,j1,kk))
                    enddo
c                  wavesum_s2(jj,i,L,kk,kk2)=wavesum_s2(jj,i,L,kk,kk2)*
c     &            MatEl_Ang_comp_sml(kk2,L,kk)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end


      subroutine MatrixMM_Storage(wave_new,nstates,nm,nkap,
     &Start_state,Mat_dimension,state_in_cont,
     &nvmat,wavesum_l1,wavesum_s1,vmat)
      include 'inc.par'
      integer Start_state,Mat_dimension,nm,nkap,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_l1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_s1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat)
      integer nstates,state_2,L,nvmat
c      common /staterangeinput/ state_range_input
      common /momentum_projection/ amu,amj_max

      wavesum_l1=0.d0
      wavesum_s1=0.d0

      do jj=Start_state,Start_state+Mat_dimension-1!End_state
        state_2=state_in_cont(jj)
        do kk=-nkap,nkap
          if(kk.ne. 0)then
            do L=0,2*nkap
c          MatEl_Ang_comp_lge=
c     &    Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
c      if(MatEl_Ang_comp_lge.gt. 1.d-12) then
              do i=1,nm
                do j=max(1,i-ns+1),min(i+ns-1,nm)

                  wavesum_l1(jj,i,L,kk)=
     &            wavesum_l1(jj,i,L,kk)+
     &            (vmat(i,j,L,nvmat)*
     &            wave_new(state_2,j,kk))
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      do jj=Start_state,Start_state+Mat_dimension-1!End_state
        state_2=state_in_cont(jj)
        do kk=-nkap,nkap
          if(kk.ne. 0)then
            do L=0,2*nkap
c   MatEl_Ang_comp_sml=
c     &   Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
c   if(MatEl_Ang_comp_sml.gt. 1.d-12) then
              do i=1,nm
                do j=max(1,i-ns+1),min(i+ns-1,nm)

                  wavesum_s1(jj,i,L,kk)=
     &            wavesum_s1(jj,i,L,kk)+
     &            (vmat(i,j,L,nvmat)*
     &            wave_new(state_2,j+nm,kk))
                enddo
              enddo
c          endif
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine MatrixMM_Storage_1(wave_new,nm,nkap,nstates,lstep,
     &nvmat,wavesum_l1,wavesum_s1,vmat,k1,i_state)
      include 'inc.par'
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_l1(nm,0:2*nkap),wavesum_s1(nm,0:2*nkap),
     &vmat(nm,nm,0:2*nkap,nvmat)

      wavesum_l1=0.d0
      wavesum_s1=0.d0

      do ll=0,2*nkap,lstep
        do i=1,nm
          do j=max(1,i-ns+1),min(i+ns-1,nm)
            wavesum_l1(i,ll)=wavesum_l1(i,ll)+
     &      (vmat(i,j,ll,nvmat)*wave_new(i_state,j,k1))
          enddo
        enddo
      enddo

      do ll=0,2*nkap,lstep
        do i=1,nm
          do j=max(1,i-ns+1),min(i+ns-1,nm)
            wavesum_s1(i,ll)=wavesum_s1(i,ll)+
     &      (vmat(i,j,ll,nvmat)*wave_new(i_state,j+nm,k1))
          enddo
        enddo
      enddo

      return
      end

      subroutine MatrixMM_Storage_even(wave_new,nstates,nm,nkap,
     &Start_state,Mat_dimension,state_in_cont,
     &nvmat,wavesum_l1,wavesum_s1,vmat)
      include 'inc.par'
      integer Start_state,Mat_dimension,nm,nkap,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_l1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_s1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat)
      integer nstates,state_2,L,nvmat
c      common /staterangeinput/ state_range_input
      common /momentum_projection/ amu,amj_max


      wavesum_l1=0.d0
      wavesum_s1=0.d0

      do jj=Start_state,Start_state+Mat_dimension-1!End_state
        state_2=state_in_cont(jj)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk.ne. 0 .and. (1+((-1)**(L_kk))).ne.0)then
            do L=0,2*nkap,2
c          MatEl_Ang_comp_lge=
c     &    Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
c      if(MatEl_Ang_comp_lge.gt. 1.d-12) then
              do i=1,nm
                do j=max(1,i-ns+1),min(i+ns-1,nm)

                  wavesum_l1(jj,i,L,kk)=
     &            wavesum_l1(jj,i,L,kk)+
     &            (vmat(i,j,L,nvmat)*
     &            wave_new(state_2,j,kk))
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      do jj=Start_state,Start_state+Mat_dimension-1!End_state
        state_2=state_in_cont(jj)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk.ne. 0 .and. (1+((-1)**(L_kk))).ne.0)then
            do L=0,2*nkap,2
c   MatEl_Ang_comp_sml=
c     &   Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
c   if(MatEl_Ang_comp_sml.gt. 1.d-12) then
              do i=1,nm
                do j=max(1,i-ns+1),min(i+ns-1,nm)

                  wavesum_s1(jj,i,L,kk)=
     &            wavesum_s1(jj,i,L,kk)+
     &            (vmat(i,j,L,nvmat)*
     &            wave_new(state_2,j+nm,kk))
                enddo
              enddo
c          endif
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine MatrixMM_Storage_odd(wave_new,nstates,nm,nkap,
     &Start_state,Mat_dimension,state_in_cont,
     &nvmat,wavesum_l1,wavesum_s1,vmat)
      include 'inc.par'
      integer Start_state,Mat_dimension,nm,nkap,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_l1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_s1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &vmat(nm,nm,0:2*nkap,nvmat)
      integer nstates,state_2,L,nvmat
c      common /staterangeinput/ state_range_input
      common /momentum_projection/ amu,amj_max


      wavesum_l1=0.d0
      wavesum_s1=0.d0

      do jj=Start_state,Start_state+Mat_dimension-1!End_state
        state_2=state_in_cont(jj)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk.ne. 0 .and. (1+((-1)**(L_kk))).eq.0)then
            do L=0,2*nkap,2
c          MatEl_Ang_comp_lge=
c     &    Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
c      if(MatEl_Ang_comp_lge.gt. 1.d-12) then
              do i=1,nm
                do j=max(1,i-ns+1),min(i+ns-1,nm)

                  wavesum_l1(jj,i,L,kk)=
     &            wavesum_l1(jj,i,L,kk)+
     &            (vmat(i,j,L,nvmat)*
     &            wave_new(state_2,j,kk))
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      do jj=Start_state,Start_state+Mat_dimension-1!End_state
        state_2=state_in_cont(jj)
        do kk=-nkap,nkap
          call DiracAngularL_int(kk,L_kk)
          if(kk.ne. 0 .and. (1+((-1)**(L_kk))).eq.0)then
            do L=0,2*nkap,2
c   MatEl_Ang_comp_sml=
c     &   Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
c   if(MatEl_Ang_comp_sml.gt. 1.d-12) then
              do i=1,nm
                do j=max(1,i-ns+1),min(i+ns-1,nm)

                  wavesum_s1(jj,i,L,kk)=
     &            wavesum_s1(jj,i,L,kk)+
     &            (vmat(i,j,L,nvmat)*
     &            wave_new(state_2,j+nm,kk))
                enddo
              enddo
c          endif
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine unrolling_nondkb_1(
     &wavesum_lge1,wave_new,wavesum_sml1,
     &nkap,nm,nstates,lstep,unroll,k1,k2,j_state,d_m1,d_m2)
      include 'inc.par'
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_lge1(nm,0:2*nkap),wavesum_sml1(nm,0:2*nkap),
     &unroll(0:2*nkap,2)
      common /momentum_projection/ amu,amj_max


      unroll=0.d0
      do ll=0,2*nkap,lstep
      ! ang_comp_l=coeff_maker(amu,1.d0*k1,1.d0*k2,ll)
        ang_comp_l=coeffMakerExplicit(d_m1,d_m2,k1,k2,ll)
        if(ang_comp_l.ne.0.d0)then
          do i=1,nm
            unroll(ll,1)=unroll(ll,1)+wavesum_lge1(i,ll)*
     &      wave_new(j_state,i,k2)
          enddo
          unroll(ll,1)=unroll(ll,1)*ang_comp_l
        endif
      enddo

      do ll=0,2*nkap,lstep
!        ang_comp_s=coeff_maker(amu,-1.d0*k1,-1.d0*k2,ll)
        ang_comp_s=coeffMakerExplicit(d_m1,d_m2,-k1,-k2,ll)
        if(ang_comp_s.ne.0.d0)then
          do i=1,nm
            unroll(ll,2)=unroll(ll,2)+wavesum_sml1(i,ll)*
     &      wave_new(j_state,i+nm,k2)
          enddo
          unroll(ll,2)=unroll(ll,2)*ang_comp_s
        endif
      enddo

      return
      end

      subroutine unrolling_nondkb(state_in_cont,
     &wavesum_lge1,wave_new,wavesum_sml1,
     &Mat_dimension,Start_state,nkap,nstates,nm,unroll)
      include 'inc.par'
      integer Start_state,state_1,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_lge1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_sml1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &unroll(Start_state:Start_state+Mat_dimension-1,
     &Start_state:Start_state+Mat_dimension-1,
     &-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &MatEl_Ang_comp_sml(-nkap:nkap,0:2*nkap,-nkap:nkap),
     &MatEl_Ang_comp_lge(-nkap:nkap,0:2*nkap,-nkap:nkap)
      common /momentum_projection/ amu,amj_max

      MatEl_Ang_comp_lge=0.d0
      MatEl_Ang_comp_sml=0.d0
      do kappa_1=-nkap,nkap
        do kappa_2=-nkap,nkap
          if(kappa_1*kappa_2.ne.0)then
            do L=0,2*nkap
              MatEl_Ang_comp_lge(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
            enddo
          endif
        enddo
      enddo

      unroll=0.d0
      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              do kappa_2=-nkap,nkap
                if(kappa_1*kappa_2.ne.0)then
                  do L=0,2*nkap
                    if(dabs(MatEl_Ang_comp_lge(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,1)+
     &                  wavesum_lge1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,1)*
     &                MatEl_Ang_comp_lge(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      do kappa_1=-nkap,nkap
        do kappa_2=-nkap,nkap
          if(kappa_1*kappa_2.ne.0)then
            do L=0,2*nkap
              MatEl_Ang_comp_sml(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
            enddo
          endif
        enddo
      enddo

      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              do kappa_2=-nkap,nkap
                if(kappa_1*kappa_2.ne.0)then
                  do L=0,2*nkap
                    if(dabs(MatEl_Ang_comp_sml(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,2)+
     &                  wavesum_sml1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i+nm,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,2)*
     &                MatEl_Ang_comp_sml(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

c      if (jj.lt.Start_state+4 .and. state_1.lt.Start_state+4)then
c      write(43,*)jj,state_1,kappa_1,kappa_2,unroll
c      endif

      return
      end

      subroutine unrolling_nondkb_even(state_in_cont,
     &wavesum_lge1,wave_new,wavesum_sml1,
     &Mat_dimension,Start_state,nkap,nstates,nm,unroll)
      include 'inc.par'
      integer Start_state,state_1,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_lge1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_sml1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &unroll(Start_state:Start_state+Mat_dimension-1,
     &Start_state:Start_state+Mat_dimension-1,
     &-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &MatEl_Ang_comp_sml(-nkap:nkap,0:2*nkap,-nkap:nkap),
     &MatEl_Ang_comp_lge(-nkap:nkap,0:2*nkap,-nkap:nkap)
      common /momentum_projection/ amu,amj_max


      MatEl_Ang_comp_lge=0.d0
      MatEl_Ang_comp_sml=0.d0
      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).ne.0 .and.
     &    (1+((-1)**(L_k2))).ne.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_lge(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
            enddo
          endif
        enddo
      enddo

      unroll=0.d0
      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and.
     &          (1+((-1)**(L_k1))).ne.0 .and.
     &          (1+((-1)**(L_k2))).ne.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_lge(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,1)+
     &                  wavesum_lge1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,1)*
     &                MatEl_Ang_comp_lge(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).ne.0 .and.
     &    (1+((-1)**(L_k2))).ne.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_sml(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
            enddo
          endif
        enddo
      enddo

      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and.
     &          (1+((-1)**(L_k1))).ne.0 .and.
     &          (1+((-1)**(L_k2))).ne.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_sml(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,2)+
     &                  wavesum_sml1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i+nm,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,2)*
     &                MatEl_Ang_comp_sml(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

c      if (jj.lt.Start_state+4 .and. state_1.lt.Start_state+4)then
c      write(43,*)jj,state_1,kappa_1,kappa_2,unroll
c      endif

      return
      end

      subroutine unrolling_nondkb_odd(state_in_cont,
     &wavesum_lge1,wave_new,wavesum_sml1,
     &Mat_dimension,Start_state,nkap,nstates,nm,unroll)
      include 'inc.par'
      integer Start_state,state_1,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_lge1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_sml1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &unroll(Start_state:Start_state+Mat_dimension-1,
     &Start_state:Start_state+Mat_dimension-1,
     &-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &MatEl_Ang_comp_sml(-nkap:nkap,0:2*nkap,-nkap:nkap),
     &MatEl_Ang_comp_lge(-nkap:nkap,0:2*nkap,-nkap:nkap)
      common /momentum_projection/ amu,amj_max


      MatEl_Ang_comp_lge=0.d0
      MatEl_Ang_comp_sml=0.d0
      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).eq.0 .and.
     &    (1+((-1)**(L_k2))).eq.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_lge(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
            enddo
          endif
        enddo
      enddo

      unroll=0.d0
      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).eq.0
     &          .and.(1+((-1)**(L_k2))).eq.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_lge(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,1)+
     &                  wavesum_lge1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,1)*
     &                MatEl_Ang_comp_lge(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).eq.0 .and.
     &    (1+((-1)**(L_k2))).eq.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_sml(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
            enddo
          endif
        enddo
      enddo

      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).eq.0
     &          .and.(1+((-1)**(L_k2))).eq.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_sml(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,2)+
     &                  wavesum_sml1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i+nm,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,2)*
     &                MatEl_Ang_comp_sml(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

c      if (jj.lt.Start_state+4 .and. state_1.lt.Start_state+4)then
c      write(43,*)jj,state_1,kappa_1,kappa_2,unroll
c      endif

      return
      end


!C     *****************************************************************
!C     ORPHANED VERSION kappa_1 and kappa_2 are not initialized here!
!C     *****************************************************************
!      subroutine unrolling_nondkb_old(jj,state_1,
!     &wavesum_lge1,wave_new,wavesum_sml1,
!     &Mat_dimension,Start_state,nkap,nstates,nm,unroll)
!      include 'inc.par'
!      integer Start_state,state_1
!      real*8 wave_new(nstates,2*nm,-nkap:nkap),
!     &wavesum_lge1(Start_state:Start_state+Mat_dimension-1,
!     &nm,0:2*nkap,-nkap:nkap),
!     &wavesum_sml1(Start_state:Start_state+Mat_dimension-1,
!     &nm,0:2*nkap,-nkap:nkap),
!     &unroll(0:2*nkap,2),
!     &MatEl_Ang_comp_sml(0:2*nkap),
!     &MatEl_Ang_comp_lge(0:2*nkap)
!      common /momentum_projection/ amu,amj_max
!
!
!      MatEl_Ang_comp_lge=0.d0
!      MatEl_Ang_comp_sml=0.d0
!
!      do L=0,2*nkap
!        MatEl_Ang_comp_lge(L)=
!     &  Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
!      enddo
!
!      unroll=0.d0
!      do L=0,2*nkap
!        if(dabs(MatEl_Ang_comp_lge(L)).gt.1.d-14)then
!          do i=1,nm
!            unroll(L,1)=unroll(L,1)+
!     &      wavesum_lge1(jj,i,L,kappa_2)*
!     &      wave_new(state_1,i,kappa_1)
!          enddo
!          unroll(L,1)=unroll(L,1)*MatEl_Ang_comp_lge(L)
!        endif
!      enddo
!      do L=0,2*nkap
!        MatEl_Ang_comp_sml(L)=
!     &  Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
!      enddo
!      do L=0,2*nkap
!        if(dabs(MatEl_Ang_comp_sml(L)).gt.1.d-14)then
!          do i=1,nm
!            unroll(L,2)=unroll(L,2)+
!     &      wavesum_sml1(jj,i,L,kappa_2)*
!     &      wave_new(state_1,i+nm,kappa_1)
!          enddo
!          unroll(L,2)=unroll(L,2)*MatEl_Ang_comp_sml(L)
!        endif
!      enddo
!
!c      if (jj.lt.Start_state+4 .and. state_1.lt.Start_state+4)then
!c      write(43,*)jj,state_1,kappa_1,kappa_2,unroll
!c      endif
!
!      return
!      end

!C     *****************************************************************
!C     ORPHANED VERSION, kappa_1 and kappa_2 are not initialized here!
!C     *****************************************************************
!      subroutine unrolling_dkb_old(jj,state_1,
!     &wavesum_lge1,wavesum_lge2,wave_new,wavesum_sml1,
!     &wavesum_sml2,Mat_dimension,Start_state,nkap,nstates,nm,unroll)
!      include 'inc.par'
!      integer Start_state,state_1
!      real*8 wave_new(nstates,2*nm,-nkap:nkap),
!     &wavesum_lge1(Start_state:Start_state+Mat_dimension-1,
!     &nm,0:2*nkap,-nkap:nkap),
!     &wavesum_lge2(Start_state:Start_state+Mat_dimension-1,
!     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
!     &wavesum_sml1(Start_state:Start_state+Mat_dimension-1,
!     &nm,0:2*nkap,-nkap:nkap),
!     &wavesum_sml2(Start_state:Start_state+Mat_dimension-1,
!     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
!     &unroll(0:2*nkap,2),
!     &MatEl_Ang_comp_sml(0:2*nkap),
!     &MatEl_Ang_comp_lge(0:2*nkap)
!      common /momentum_projection/ amu,amj_max
!
!
!      MatEl_Ang_comp_lge=0.d0
!      MatEl_Ang_comp_sml=0.d0
!
!      do L=0,2*nkap
!        MatEl_Ang_comp_lge(L)=
!     &  Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
!      enddo
!
!      unroll=0.d0
!      do L=0,2*nkap
!        if(dabs(MatEl_Ang_comp_lge(L)).gt.1.d-14)then
!          do i=1,nm
!            unroll(L,1)=unroll(L,1)+
!c     & wavesum_lge1(jj,i,L,kappa_2,kappa_1)*
!     &      wavesum_lge1(jj,i,L,kappa_2)*
!     &      wave_new(state_1,i,kappa_1)+
!     &      wavesum_lge2(jj,i,L,kappa_2,kappa_1)*
!     &      wave_new(state_1,i+nm,kappa_1)
!          enddo
!          unroll(L,1)=unroll(L,1)*MatEl_Ang_comp_lge(L)
!        endif
!      enddo
!      do L=0,2*nkap
!        MatEl_Ang_comp_sml(L)=
!     &  Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
!      enddo
!
!      do L=0,2*nkap
!        if(dabs(MatEl_Ang_comp_sml(L)).gt.1.d-14)then
!          do i=1,nm
!            unroll(L,2)=unroll(L,2)+
!c     & wavesum_sml1(jj,i,L,kappa_2,kappa_1)*
!     &      wavesum_sml1(jj,i,L,kappa_2)*
!     &      wave_new(state_1,i+nm,kappa_1)+
!     &      wavesum_sml2(jj,i,L,kappa_2,kappa_1)*
!     &      wave_new(state_1,i,kappa_1)
!          enddo
!          unroll(L,2)=unroll(L,2)*MatEl_Ang_comp_sml(L)
!        endif
!      enddo
!
!c      if (jj.lt.Start_state+4 .and. state_1.lt.Start_state+4)then
!c      write(44,*)jj,state_1,kappa_1,kappa_2,unroll
!c      endif
!
!
!      return
!      end

      subroutine unrolling_dkb(wavesum_lge1,wavesum_lge2,wave_new,
     &wavesum_sml1,wavesum_sml2,nkap,nm,nstates,lstep,unroll,k1,
     &k2,j_state,d_m1,d_m2)
      include 'inc.par'
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_lge1(nm,0:2*nkap),wavesum_lge2(nm,0:2*nkap),
     &wavesum_sml1(nm,0:2*nkap),wavesum_sml2(nm,0:2*nkap),
     &unroll(0:2*nkap,2)

      unroll=0.d0
      do ll=0,2*nkap,lstep
!        ang_comp_l=coeff_maker(amu,1.d0*k1,1.d0*k2,ll)
        ang_comp_l=coeffMakerExplicit(d_m1,d_m2,k1,k2,ll)
        if(ang_comp_l.ne.0.d0)then
          do i=1,nm
            unroll(ll,1)=unroll(ll,1)+
     &      wavesum_lge1(i,ll)*wave_new(j_state,i,k2)+
     &      wavesum_lge2(i,ll)*wave_new(j_state,i+nm,k2)
          enddo
          unroll(ll,1)=unroll(ll,1)*ang_comp_l
        endif
      enddo

      do ll=0,2*nkap,lstep
!        ang_comp_s=coeff_maker(amu,-1.d0*k1,-1.d0*k2,ll)
        ang_comp_s=coeffMakerExplicit(d_m1,d_m2,-k1,-k2,ll)
        if(ang_comp_s.ne.0.d0)then
          do i=1,nm
            unroll(ll,2)=unroll(ll,2)+
     &      wavesum_sml1(i,ll)*wave_new(j_state,i+nm,k2)+
     &      wavesum_sml2(i,ll)*wave_new(j_state,i,k2)
          enddo
          unroll(ll,2)=unroll(ll,2)*ang_comp_s
        endif
      enddo
      return
      end

      subroutine unrolling_dkb_even(state_in_cont,
     &wavesum_lge1,wavesum_lge2,wave_new,wavesum_sml1,
     &wavesum_sml2,Mat_dimension,Start_state,nkap,nstates,nm,unroll)
      include 'inc.par'
      integer Start_state,state_1,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_lge1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_lge2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &wavesum_sml1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_sml2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &unroll(Start_state:Start_state+Mat_dimension-1,
     &Start_state:Start_state+Mat_dimension-1,
     &-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &MatEl_Ang_comp_sml(-nkap:nkap,0:2*nkap,-nkap:nkap),
     &MatEl_Ang_comp_lge(-nkap:nkap,0:2*nkap,-nkap:nkap)
      common /momentum_projection/ amu,amj_max


      MatEl_Ang_comp_lge=0.d0
      MatEl_Ang_comp_sml=0.d0
      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).ne.0
     &    .and. (1+((-1)**(L_k2))).ne.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_lge(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
            enddo
          endif
        enddo
      enddo

      unroll=0.d0
      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).ne.0
     &          .and.(1+((-1)**(L_k2))).ne.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_lge(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,1)+
c     & wavesum_lge1(jj,i,L,kappa_2,kappa_1)*
     &                  wavesum_lge1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i,kappa_1)+
     &                  wavesum_lge2(jj,i,L,kappa_2,kappa_1)*
     &                  wave_new(state_1,i+nm,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,1)*
     &                MatEl_Ang_comp_lge(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).ne.0
     &    .and.(1+((-1)**(L_k2))).ne.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_sml(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
            enddo
          endif
        enddo
      enddo

      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).ne.0
     &          .and.(1+((-1)**(L_k2))).ne.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_sml(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,2)+
c     & wavesum_sml1(jj,i,L,kappa_2,kappa_1)*
     &                  wavesum_sml1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i+nm,kappa_1)+
     &                  wavesum_sml2(jj,i,L,kappa_2,kappa_1)*
     &                  wave_new(state_1,i,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,2)*
     &                MatEl_Ang_comp_sml(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

c      if (jj.lt.Start_state+4 .and. state_1.lt.Start_state+4)then
c      write(44,*)jj,state_1,kappa_1,kappa_2,unroll
c      endif


      return
      end

      subroutine unrolling_dkb_odd(state_in_cont,
     &wavesum_lge1,wavesum_lge2,wave_new,wavesum_sml1,
     &wavesum_sml2,
     &Mat_dimension,Start_state,nkap,nstates,nm,unroll)
      include 'inc.par'
      integer Start_state,state_1,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap),
     &wavesum_lge1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_lge2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &wavesum_sml1(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap),
     &wavesum_sml2(Start_state:Start_state+Mat_dimension-1,
     &nm,0:2*nkap,-nkap:nkap,-nkap:nkap),
     &unroll(Start_state:Start_state+Mat_dimension-1,
     &Start_state:Start_state+Mat_dimension-1,
     &-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &MatEl_Ang_comp_sml(-nkap:nkap,0:2*nkap,-nkap:nkap),
     &MatEl_Ang_comp_lge(-nkap:nkap,0:2*nkap,-nkap:nkap)
      common /momentum_projection/ amu,amj_max


      MatEl_Ang_comp_lge=0.d0
      MatEl_Ang_comp_sml=0.d0
      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).eq.0 .and.
     &    (1+((-1)**(L_k2))).eq.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_lge(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
            enddo
          endif
        enddo
      enddo

      unroll=0.d0
      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k1))).eq.0
     &          .and.(1+((-1)**(L_k2))).eq.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_lge(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,1)+
c     & wavesum_lge1(jj,i,L,kappa_2,kappa_1)*
     &                  wavesum_lge1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i,kappa_1)+
     &                  wavesum_lge2(jj,i,L,kappa_2,kappa_1)*
     &                  wave_new(state_1,i+nm,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,1)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,1)*
     &                MatEl_Ang_comp_lge(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      do kappa_1=-nkap,nkap
        call DiracAngularL_int(kappa_1,L_k1)
        do kappa_2=-nkap,nkap
          call DiracAngularL_int(kappa_2,L_k2)
          if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k2))).eq.0 .and.
     &    (1+((-1)**(L_k1))).eq.0)then
            do L=0,2*nkap,2
              MatEl_Ang_comp_sml(kappa_1,L,kappa_2)=
     &        Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
            enddo
          endif
        enddo
      enddo

      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.lt.jj)then
            do kappa_1=-nkap,nkap
              call DiracAngularL_int(kappa_1,L_k1)
              do kappa_2=-nkap,nkap
                call DiracAngularL_int(kappa_2,L_k2)
                if(kappa_1*kappa_2.ne.0 .and. (1+((-1)**(L_k2))).eq.0
     &          .and.(1+((-1)**(L_k1))).eq.0)then
                  do L=0,2*nkap,2
                    if(dabs(MatEl_Ang_comp_sml(kappa_1,L,kappa_2))
     &              .gt.1.d-14)then
                      do i=1,nm
                        unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                  unroll(ii,jj,kappa_1,kappa_2,L,2)+
c     & wavesum_sml1(jj,i,L,kappa_2,kappa_1)*
     &                  wavesum_sml1(jj,i,L,kappa_2)*
     &                  wave_new(state_1,i+nm,kappa_1)+
     &                  wavesum_sml2(jj,i,L,kappa_2,kappa_1)*
     &                  wave_new(state_1,i,kappa_1)
                      enddo
                      unroll(ii,jj,kappa_1,kappa_2,L,2)=
     &                unroll(ii,jj,kappa_1,kappa_2,L,2)*
     &                MatEl_Ang_comp_sml(kappa_1,L,kappa_2)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

c      if (jj.lt.Start_state+4 .and. state_1.lt.Start_state+4)then
c      write(44,*)jj,state_1,kappa_1,kappa_2,unroll
c      endif


      return
      end



C     ******************************************************************
C     ORPHANED VERSION
C     ******************************************************************
      subroutine MatrixMMGenerator_neqZ_old(Mat_dimension,Start_state,
     &dTdXi,dRdXi,eigval_a,nkap,vmat,wave_new,nstates,nm,nvmat,mm,
     &dvdRmatdkb1,dvdRmatdkb2,state_in_cont)
      include 'inc.par'
      integer Mat_dimension,Start_state,state_1,
     &state_2,nstates,nkap,nm,nvmat,
     &state_in_cont(Start_state:Start_state+Mat_dimension-1)
      real*8 dTdXi,eigval_a(nstates),mm0,Cpld_chnnl_elem_lge,dRdXi,
     &Cpld_chnnl_elem_sml!MatEl_Ang_comp_lge,MatEl_Ang_comp_sml
      real*8 vmat(nm,nm,0:2*nkap,nvmat)
      real*8 wave_new(nstates,2*nm,-nkap:nkap)
      real*8 dvdRmatdkb2(nm,nm,0:2*nkap,-nkap:nkap,-nkap:nkap,2)
      real*8 dvdRmatdkb1(nm,nm,0:2*nkap,-nkap:nkap,2),
     &MatEl_Ang_comp_lge,MatEl_Ang_comp_sml
      complex*16 mm(Start_state:Start_state+Mat_dimension-1,
     &Start_state:Start_state+Mat_dimension-1)
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
c      common /staterangeinput/ state_range_input
      logical dkb!,state_range_input

      write(*,*) 'ENTER MAIN MATRIX'

      do ii=Start_state,Start_state+Mat_dimension-1!End_state
        state_1=state_in_cont(ii)
        do jj=Start_state,Start_state+Mat_dimension-1!End_state
          state_2=state_in_cont(jj)

          if (state_1.lt.state_2) then
            mm0=0.d0
            do kappa_1=-nkap,nkap
              do kappa_2=-nkap,nkap
                if (kappa_1*kappa_2.ne. 0) then
                  do L=0,2*nkap
                    Cpld_chnnl_elem_lge=0.d0
                    Cpld_chnnl_elem_sml=0.d0

                    MatEl_Ang_comp_lge=
     &              Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
                    MatEl_Ang_comp_sml=
     &              Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
                    do i=1,nm
                      do j=max(1,i-ns+1),min(i+ns-1,nm)
                        if(dkb) then
                          if(dabs(MatEl_Ang_comp_lge).gt.1.d-14)then
                            Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
     &                      (vmat(i,j,L,nvmat)*
     &                      wave_new(state_1,i,kappa_1)*
     &                      wave_new(state_2,j,kappa_2)+
     &                      dvdRmatdkb1(j,i,L,kappa_2,1)*
     &                      wave_new(state_1,i,kappa_1)*
     &                      wave_new(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb1(i,j,L,kappa_1,1)*
     &                      wave_new(state_1,i+nm,kappa_1)*
     &                      wave_new(state_2,j,kappa_2)+
     &                      dvdRmatdkb2(i,j,L,kappa_1,kappa_2,1)*
     &                      wave_new(state_2,j+nm,kappa_2)*
     &                      wave_new(state_1,i+nm,kappa_1)
     &                      )!*MatEl_Ang_comp_lge
                          endif
c Notice here that the small components are logged in wave starting from
c nm+1.
                          if(dabs(MatEl_Ang_comp_sml).gt.1.d-14)then
                            Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
     &                      (vmat(i,j,L,nvmat)*
     &                      wave_new(state_1,i+nm,kappa_1)*
     &                      wave_new(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb1(j,i,L,kappa_2,2)*
     &                      wave_new(state_2,j,kappa_2)*
     &                      wave_new(state_1,i+nm,kappa_1)+
     &                      dvdRmatdkb1(i,j,L,kappa_1,2)*
     &                      wave_new(state_1,i,kappa_1)*
     &                      wave_new(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb2(i,j,L,kappa_1,kappa_2,2)*
     &                      wave_new(state_2,j,kappa_2)*
     &                      wave_new(state_1,i,kappa_1)
     &                      )!*MatEl_Ang_comp_sml
                          endif
                        else
                          if(dabs(MatEl_Ang_comp_lge).gt.1.d-14)then
                            Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
     &                      vmat(i,j,L,nvmat)*
     &                      wave_new(state_1,i,kappa_1)*
     &                      wave_new(state_2,j,kappa_2)!*
c     &                   MatEl_Ang_comp_lge
                          endif
c Notice here that the small components are logged in wave starting from
c nm+1.
                          if(dabs(MatEl_Ang_comp_lge).gt.1.d-14)then
                            Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
     &                      vmat(i,j,L,nvmat)*
     &                      wave_new(state_1,i+nm,kappa_1)*
     &                      wave_new(state_2,j+nm,kappa_2)!*
c     &                   MatEl_Ang_comp_sml
                          endif
                        endif
                      enddo
                    enddo
c      if(dabs(MatEl_Ang_comp_lge).gt.1.d-14)write(*,*)
c     &      MatEl_Ang_comp_lge,
c     &      dabs(Cpld_chnnl_elem_lge+Cpld_chnnl_elem_sml),
c     &      kappa_1,kappa_2,L
                    mm0=mm0+(Cpld_chnnl_elem_lge*MatEl_Ang_comp_lge+
     &              Cpld_chnnl_elem_sml*MatEl_Ang_comp_sml)
                  enddo
c      write(*,*)mm0,ii-Start_state,jj-Start_state
                endif
              enddo
            enddo
            mm(ii,jj)=-(0.d0,1.d0)*dRdXi*mm0/
     &      (eigval_a(state_2)-eigval_a(state_1))
          endif
        enddo
      enddo

      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=state_in_cont(ii)
        mm(ii,ii)=dTdXi*eigval_a(state_1)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.gt. jj)then
            mm(ii,jj)=-mm(jj,ii)
          endif
          if(dkb)write(41,*)ii,jj,mm(ii,jj),eigval_a(jj)
          if(.not.dkb)write(42,*)ii,jj,mm(ii,jj),eigval_a(jj)
        enddo
      enddo
      stop
      return
      end

      subroutine projection_matrix_frozen_basis(nstates,nm,
     &nkap,dmat,wave1,wave2,projMat)
      include 'inc.par'
      real*8 projMat(nstates,nstates),wave1(nstates,2*nm,-nkap:nkap),
     &wave2(nstates,2*nm,-nkap:nkap),dmat(2*nm,2*nm,-nkap:nkap)

      projMat = 0.d0
      do i1 = 1,nstates
        do i2 = i1,nstates
          do kap = -nkap,nkap
            if (kap .eq. 0)then
              cycle
            endif
            do iw1 = 1,2*nm
              do iw2 = max(1,iw1-ns+1),min(iw1+ns-1,2*nm)
                projMat(i1,i2)=projMat(i1,i2)+wave1(i1,iw1,kap)*
     &          wave2(i2,iw2,kap)*dmat(iw1,iw2,kap)
              enddo
            enddo
          enddo
          if (i1.ne.i2)then
            projMat(i2,i1)=projMat(i1,i2)
          endif
        enddo
      enddo
      return
      end

      subroutine MatrixMMGenerator_neqZ(dTdXi,dRdXi,eigval,nkap,vmat,
     &wave_new,nstates,nm,nvmat,mm,dvdRmatdkb1,dvdRmatdkb2,
     &d_number_states_mj)
      include 'inc.par'
      real*8, dimension(:,:),allocatable:: wavesum_lge1,wavesum_sml1,
     &unroll,wavesum_lge2,wavesum_sml2
      real*8 eigval(nstates),vmat(nm,nm,0:2*nkap,nvmat),
     &wave_new(nstates,2*nm,-nkap:nkap),
     &dvdRmatdkb2(nm,nm,-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &dvdRmatdkb1(nm,nm,-nkap:nkap,0:2*nkap,2),
     &d_number_states_mj(nstates)
      complex*16 mm(nstates,nstates)
      common /common_dkb/ dkb
      logical dkb

      mm=(0.d0,0.d0)

      if(dkb) then
!$OMP PARALLEL
!$OMP&PRIVATE(id,i,d_m1,j,d_m2,k1factor,k1,ki,wavesum_lge1,
!$OMP&wavesum_sml1,k2factor,k2,kj,wavesum_lge2,wavesum_sml2,unroll,l)
!$OMP&SHARED(nstates,d_number_states_mj,nkap,nm,wave_new,dvdRmatdkb1,
!$OMP&nvmat,vmat,dvdRmatdkb2,mm,drdxi,eigval), default(none)
!$OMP DO
         do i=1,nstates
            d_m1=d_number_states_mj(i)
            do j=i+1,nstates
               d_m2=d_number_states_mj(j)
               if (d_m1 .ne. d_m2) cycle
               do ki=-nkap,nkap
                  if(ki.ne.0)then
                     allocate(wavesum_lge1(nm,0:2*nkap))
                     allocate(wavesum_sml1(nm,0:2*nkap))
                     call MatrixMM_Storage_DKB_1(wave_new,nm,nkap,
     &               nstates,1,dvdRmatdkb1,nvmat,wavesum_lge1,
     &               wavesum_sml1,vmat,ki,i)
                     do kj=-nkap,nkap
                        if(kj.ne.0)then
                           allocate(wavesum_lge2(nm,0:2*nkap))
                           allocate(wavesum_sml2(nm,0:2*nkap))
                           call MatrixMM_Storage_DKB_2(wave_new,nm,nkap,
     &                     nstates,1,dvdRmatdkb2,dvdRmatdkb1,
     &                     wavesum_lge2,wavesum_sml2,ki,i,kj,d_m1,d_m2)
                           allocate(unroll(0:2*nkap,2))
                           call unrolling_dkb(wavesum_lge1,wavesum_lge2,
     &                     wave_new,wavesum_sml1,wavesum_sml2,nkap,nm,
     &                     nstates,1,unroll,ki,kj,j,d_m1,d_m2)
                           do l=0,2*nkap
                              mm(i,j)=mm(i,j)+unroll(l,1)+
     &                        unroll(l,2)
                           enddo
                           deallocate(unroll)
                           deallocate(wavesum_lge2)
                           deallocate(wavesum_sml2)
                        endif
                     enddo
                     deallocate(wavesum_lge1)
                     deallocate(wavesum_sml1)
                  endif
               enddo
               if(dabs(eigval(j)-eigval(i)).gt.1.d-9)then
                  mm(i,j)=mm(i,j)*drdxi*(0,-1)/(eigval(j)-eigval(i))
               else
                  mm(i,j)=(0.d0,0.d0)
               endif
               mm(j,i)=-mm(i,j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      else
!$OMP PARALLEL
!$OMP&PRIVATE(i,d_m1,j,d_m2,k1factor,k1,ki,wavesum_lge1,
!$OMP&wavesum_sml1,k2factor,k2,kj,unroll,l)
!$OMP&SHARED(nstates,d_number_states_mj,nkap,nm,wave_new,dvdRmatdkb1,
!$OMP&nvmat,vmat,dvdRmatdkb2,mm,drdxi,eigval), default(none)
!$OMP DO
         do i=1,nstates
            d_m1=d_number_states_mj(i)
            do j=i+1,nstates
               d_m2=d_number_states_mj(j)
               if (d_m1 .ne. d_m2) cycle
               do ki=-nkap,nkap
                  if(ki.ne.0)then
                     allocate(wavesum_lge1(nm,0:2*nkap))
                     allocate(wavesum_sml1(nm,0:2*nkap))
                     call MatrixMM_Storage_1(wave_new,nm,nkap,nstates,
     &               1,nvmat,wavesum_lge1,wavesum_sml1,vmat,ki,i)
                     do kj=-nkap,nkap
                        if(kj.ne.0)then
                           allocate(unroll(0:2*nkap,2))
                           call unrolling_nondkb_1(wavesum_lge1,wave_new
     &                     ,wavesum_sml1,nkap,nm,nstates,1,unroll,ki,kj,
     &                     j,d_m1,d_m2)
                           do l=0,2*nkap
                              mm(i,j)=mm(i,j)+unroll(l,1)+
     &                        unroll(l,2)
                           enddo
                           deallocate(unroll)
                        endif
                     enddo
                     deallocate(wavesum_lge1)
                     deallocate(wavesum_sml1)
                  endif
               enddo
               if(dabs(eigval(j)-eigval(i)).gt.1.d-9)then
                  mm(i,j)=mm(i,j)*drdxi*(0,-1)/(eigval(j)-eigval(i))
               else
                  mm(i,j)=(0.d0,0.d0)
               endif
               mm(j,i)=-mm(i,j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      endif
      do i=1,nstates
        mm(i,i)=dtdxi*eigval(i)
      enddo
      return
      end

      subroutine MatrixMMGenerator_eqZeven(dTdXi,dRdXi,eigval,nkap,vmat,
     &wave_new,nstates,nm,nvmat,mm,dvdRmatdkb1,dvdRmatdkb2,
     &d_number_states_mj)
      include 'inc.par'
      real*8, dimension(:,:),allocatable:: wavesum_lge1,wavesum_sml1,
     &unroll,wavesum_lge2,wavesum_sml2
      real*8 eigval(nstates),vmat(nm,nm,0:2*nkap,nvmat),
     &wave_new(nstates,2*nm,-nkap:nkap),
     &dvdRmatdkb2(nm,nm,-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &dvdRmatdkb1(nm,nm,-nkap:nkap,0:2*nkap,2),
     &d_number_states_mj(nstates)
      complex*16 mm(nstates,nstates)
      common /common_dkb/ dkb
      logical dkb

      mm=(0.d0,0.d0)

      if(dkb) then
!$OMP PARALLEL
!$OMP&PRIVATE(id,i,d_m1,j,d_m2,k1factor,k1,ki,wavesum_lge1,
!$OMP&wavesum_sml1,k2factor,k2,kj,wavesum_lge2,wavesum_sml2,unroll,l)
!$OMP&SHARED(nstates,d_number_states_mj,nkap,nm,wave_new,dvdRmatdkb1,
!$OMP&nvmat,vmat,dvdRmatdkb2,mm,drdxi,eigval), default(none)
!$OMP DO
        do i=1,nstates
          d_m1=d_number_states_mj(i)
          do j=i+1,nstates
            d_m2=d_number_states_mj(j)
            if (d_m1 .ne. d_m2) cycle
            k1factor=(-1)**nkap
            do k1=nkap,1,-1
              ki=k1factor*k1
              k1factor=-k1factor
              allocate(wavesum_lge1(nm,0:2*nkap))
              allocate(wavesum_sml1(nm,0:2*nkap))
              call MatrixMM_Storage_DKB_1(wave_new,nm,nkap,nstates,
     &        2,dvdRmatdkb1,nvmat,wavesum_lge1,wavesum_sml1,vmat,ki,i)
              k2factor=(-1)**nkap
              do k2=nkap,1,-1
                kj=k2*k2factor
                k2factor=-k2factor
                allocate(wavesum_lge2(nm,0:2*nkap))
                allocate(wavesum_sml2(nm,0:2*nkap))
                call MatrixMM_Storage_DKB_2(wave_new,nm,nkap,
     &          nstates,2,dvdRmatdkb2,dvdRmatdkb1,wavesum_lge2,
     &          wavesum_sml2,ki,i,kj,d_m1,d_m2)
                allocate(unroll(0:2*nkap,2))
                call unrolling_dkb(wavesum_lge1,wavesum_lge2,
     &          wave_new,wavesum_sml1,wavesum_sml2,nkap,nm,
     &          nstates,2,unroll,ki,kj,j,d_m1,d_m2)
                do l=0,2*nkap,2
                  mm(i,j)=mm(i,j)+unroll(l,1)+unroll(l,2)
                enddo
                deallocate(unroll)
                deallocate(wavesum_lge2)
                deallocate(wavesum_sml2)
              enddo
              deallocate(wavesum_lge1)
              deallocate(wavesum_sml1)
            enddo
            if(dabs(eigval(j)-eigval(i)).gt.1.d-9)then
              mm(i,j)=mm(i,j)*drdxi*(0,-1)/(eigval(j)-eigval(i))
            else
              mm(i,j)=(0.d0,0.d0)
            endif
            mm(j,i)=-mm(i,j)
          enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
      else
!$OMP PARALLEL
!$OMP&PRIVATE(i,d_m1,j,d_m2,k1factor,k1,ki,wavesum_lge1,
!$OMP&wavesum_sml1,k2factor,k2,kj,unroll,l)
!$OMP&SHARED(nstates,d_number_states_mj,nkap,nm,wave_new,dvdRmatdkb1,
!$OMP&nvmat,vmat,dvdRmatdkb2,mm,drdxi,eigval), default(none)
!$OMP DO
        do i=1,nstates
          d_m1=d_number_states_mj(i)
          do j=i+1,nstates
            d_m2=d_number_states_mj(j)
            if (d_m1 .ne. d_m2) cycle
            k1factor=(-1)**nkap
            do k1=nkap,1,-1
              ki=k1*k1factor
              k1factor=-k1factor
              allocate(wavesum_lge1(nm,0:2*nkap))
              allocate(wavesum_sml1(nm,0:2*nkap))
              call MatrixMM_Storage_1(wave_new,nm,nkap,nstates,
     &        2,nvmat,wavesum_lge1,wavesum_sml1,vmat,ki,i)
              k2factor=(-1)**nkap
              do k2=nkap,1,-1
                kj=k2*k2factor
                k2factor=-k2factor
                allocate(unroll(0:2*nkap,2))
                call unrolling_nondkb_1(wavesum_lge1,wave_new,
     &          wavesum_sml1,nkap,nm,nstates,2,unroll,ki,kj,j,d_m1,d_m2)
                do l=0,2*nkap,2
                  mm(i,j)=mm(i,j)+unroll(l,1)+unroll(l,2)
                enddo
                deallocate(unroll)
              enddo
              deallocate(wavesum_lge1)
              deallocate(wavesum_sml1)
            enddo
            if(dabs(eigval(j)-eigval(i)).gt.1.d-9)then
              mm(i,j)=mm(i,j)*drdxi*(0,-1)/(eigval(j)-eigval(i))
            else
              mm(i,j)=(0.d0,0.d0)
            endif
            mm(j,i)=-mm(i,j)
          enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
      endif
      do i=1,nstates
        mm(i,i)=dtdxi*eigval(i)
      enddo
      return
      end

      subroutine MatrixMMGenerator_eqZodd(dTdXi,dRdXi,eigval,nkap,vmat,
     &wave_new,nstates,nm,nvmat,mm,dvdRmatdkb1,dvdRmatdkb2,
     &d_number_states_mj)
      include 'inc.par'
      real*8, dimension(:,:),allocatable:: wavesum_lge1,wavesum_sml1,
     &unroll,wavesum_lge2,wavesum_sml2
      real*8 eigval(nstates),vmat(nm,nm,0:2*nkap,nvmat),
     &wave_new(nstates,2*nm,-nkap:nkap),
     &dvdRmatdkb2(nm,nm,-nkap:nkap,-nkap:nkap,0:2*nkap,2),
     &dvdRmatdkb1(nm,nm,-nkap:nkap,0:2*nkap,2),
     &d_number_states_mj(nstates)
      complex*16 mm(nstates,nstates)
      common /common_dkb/ dkb
      logical dkb

      mm=(0.d0,0.d0)

      if(dkb) then
!$OMP PARALLEL
!$OMP&PRIVATE(id,i,d_m1,j,d_m2,k1factor,k1,ki,wavesum_lge1,
!$OMP&wavesum_sml1,k2factor,k2,kj,wavesum_lge2,wavesum_sml2,unroll,l)
!$OMP&SHARED(nstates,d_number_states_mj,nkap,nm,wave_new,dvdRmatdkb1,
!$OMP&nvmat,vmat,dvdRmatdkb2,mm,drdxi,eigval), default(none)
!$OMP DO
        do i=1,nstates
          d_m1=d_number_states_mj(i)
          do j=i+1,nstates
            d_m2=d_number_states_mj(j)
            if (d_m1 .ne. d_m2) cycle
            k1factor=(-1)**(nkap+1)
            do k1=nkap,1,-1
              ki=k1*k1factor
              k1factor=-k1factor
              allocate(wavesum_lge1(nm,0:2*nkap))
              allocate(wavesum_sml1(nm,0:2*nkap))
              call MatrixMM_Storage_DKB_1(wave_new,nm,nkap,nstates,
     &        2,dvdRmatdkb1,nvmat,wavesum_lge1,wavesum_sml1,vmat,ki,i)
              k2factor=(-1)**(nkap+1)
              do k2=nkap,1,-1
                kj=k2*k2factor
                k2factor=-k2factor
                allocate(wavesum_lge2(nm,0:2*nkap))
                allocate(wavesum_sml2(nm,0:2*nkap))
                call MatrixMM_Storage_DKB_2(wave_new,nm,nkap,
     &          nstates,2,dvdRmatdkb2,dvdRmatdkb1,wavesum_lge2,
     &          wavesum_sml2,ki,i,kj,d_m1,d_m2)
                allocate(unroll(0:2*nkap,2))
                call unrolling_dkb(wavesum_lge1,wavesum_lge2,
     &          wave_new,wavesum_sml1,wavesum_sml2,nkap,nm,
     &          nstates,2,unroll,ki,kj,j,d_m1,d_m2)
                do l=0,2*nkap,2
                  mm(i,j)=mm(i,j)+unroll(l,1)+unroll(l,2)
                enddo
                deallocate(unroll)
                deallocate(wavesum_lge2)
                deallocate(wavesum_sml2)
              enddo
              deallocate(wavesum_lge1)
              deallocate(wavesum_sml1)
            enddo
            if(dabs(eigval(j)-eigval(i)).gt.1.d-9)then
              mm(i,j)=mm(i,j)*drdxi*(0,-1)/(eigval(j)-eigval(i))
            else
              mm(i,j)=(0.d0,0.d0)
            endif
            mm(j,i)=-mm(i,j)
          enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
      else
!$OMP PARALLEL
!$OMP&PRIVATE(i,d_m1,j,d_m2,k1factor,k1,ki,wavesum_lge1,
!$OMP&wavesum_sml1,k2factor,k2,kj,unroll,l)
!$OMP&SHARED(nstates,d_number_states_mj,nkap,nm,wave_new,dvdRmatdkb1,
!$OMP&nvmat,vmat,dvdRmatdkb2,mm,drdxi,eigval), default(none)
!$OMP DO
        do i=1,nstates
          d_m1=d_number_states_mj(i)
          do j=i+1,nstates
            d_m2=d_number_states_mj(j)
            if (d_m1 .ne. d_m2) cycle
            k1factor=(-1)**(nkap+1)
            do k1=nkap,1,-1
              ki=k1*k1factor
              k1factor=-k1factor
              allocate(wavesum_lge1(nm,0:2*nkap))
              allocate(wavesum_sml1(nm,0:2*nkap))
              call MatrixMM_Storage_1(wave_new,nm,nkap,nstates,
     &        2,nvmat,wavesum_lge1,wavesum_sml1,vmat,ki,i)
              k2factor=(-1)**(nkap+1)
              do k2=nkap,1,-1
                kj=k2*k2factor
                k2factor=-k2factor
                allocate(unroll(0:2*nkap,2))
                call unrolling_nondkb_1(wavesum_lge1,wave_new,
     &          wavesum_sml1,nkap,nm,nstates,2,unroll,ki,kj,j,d_m1,d_m2)
                do l=0,2*nkap,2
                  mm(i,j)=mm(i,j)+unroll(l,1)+unroll(l,2)
                enddo
                deallocate(unroll)
              enddo
              deallocate(wavesum_lge1)
              deallocate(wavesum_sml1)
            enddo
            if(dabs(eigval(j)-eigval(i)).gt.1.d-9)then
              mm(i,j)=mm(i,j)*drdxi*(0,-1)/(eigval(j)-eigval(i))
            else
              mm(i,j)=(0.d0,0.d0)
            endif
            mm(j,i)=-mm(i,j)
          enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
      endif
      do i=1,nstates
        mm(i,i)=dtdxi*eigval(i)
      enddo
      return
      end

C     ******************************************************************
C     ORPHANED VERSION
C     ******************************************************************
      subroutine MatrixMMGenerator_eqZeven_old
     &(Mat_dimension,Start_state,dTdXi,dRdXi,
     &eigval_e,nkap,vmat,wave_new_even,nste,nm,nvmat,
     &mmeven,dvdRmatdkb1,dvdRmatdkb2,states_in_cont_e)
      include 'inc.par'
      integer Mat_dimension,Start_state,state_1,
     &state_2,nste,nkap,nm,nvmat,
     &states_in_cont_e(Start_state:Start_state+Mat_dimension-1)
      real*8 dTdXi,eigval_e(nste),mm0,Cpld_chnnl_elem_lge,dRdXi,
     &Cpld_chnnl_elem_sml,MatEl_Ang_comp_lge,MatEl_Ang_comp_sml,
     &vmat(nm,nm,0:2*nkap,nvmat),
     &wave_new_even(nste,2*nm,-nkap:nkap)
      real*8 dvdRmatdkb2(nm,nm,0:2*nkap,-nkap:nkap,-nkap:nkap,2)
      real*8 dvdRmatdkb1(nm,nm,0:2*nkap,-nkap:nkap,2)
      complex*16 mmeven(Start_state:Start_state+2*Mat_dimension-1,
     &Start_state:Start_state+2*Mat_dimension-1)
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
c      common /staterangeinput/ state_range_input
      logical dkb!,state_range_input

      mmeven=0.d0
      do ii=Start_state,Start_state+Mat_dimension-1!End_state
        state_1=states_in_cont_e(ii)
        do jj=Start_state,Start_state+Mat_dimension-1!End_state
          state_2=states_in_cont_e(jj)
          mm0=0.d0
          do kappa_1=-nkap,nkap
            do kappa_2=-nkap,nkap
              Cpld_chnnl_elem_lge=0.d0
              Cpld_chnnl_elem_sml=0.d0
              if (kappa_1*kappa_2.ne. 0) then
                do L=0,2*nkap,2
                  if (state_1.ge.state_2) then
                    Cpld_chnnl_elem_lge=0.d0
                    Cpld_chnnl_elem_sml=0.d0
                  else
                    MatEl_Ang_comp_lge=
     &              Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
                    MatEl_Ang_comp_sml=
     &              Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
                    do i=1,nm
                      do j=max(1,i-ns+1),min(i+ns-1,nm)!1,nm

                        if(dkb) then
                          if(dabs(MatEl_Ang_comp_lge).gt. 1.d-12) then
                            Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
     &                      (vmat(i,j,L,nvmat)*
     &                      wave_new_even(state_1,i,kappa_1)*
     &                      wave_new_even(state_2,j,kappa_2)+
     &                      dvdRmatdkb1(j,i,L,kappa_2,1)*
     &                      wave_new_even(state_1,i,kappa_1)*
     &                      wave_new_even(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb1(i,j,L,kappa_1,1)*
     &                      wave_new_even(state_1,i+nm,kappa_1)*
     &                      wave_new_even(state_2,j,kappa_2)+
     &                      dvdRmatdkb2(i,j,L,kappa_1,kappa_2,1)*
     &                      wave_new_even(state_2,j+nm,kappa_2)*
     &                      wave_new_even(state_1,i+nm,kappa_1)
     &                      )*
     &                      MatEl_Ang_comp_lge
                          endif
c Notice here that the small components are logged in wave starting from
c nm+1.
                          if(dabs(MatEl_Ang_comp_sml).gt. 1.d-12) then
                            Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
     &                      (vmat(i,j,L,nvmat)*
     &                      wave_new_even(state_1,i+nm,kappa_1)*
     &                      wave_new_even(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb1(j,i,L,kappa_2,2)*
     &                      wave_new_even(state_1,j,kappa_2)*
     &                      wave_new_even(state_2,i+nm,kappa_1)+
     &                      dvdRmatdkb1(i,j,L,kappa_1,2)*
     &                      wave_new_even(state_2,i,kappa_1)*
     &                      wave_new_even(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb2(i,j,L,kappa_1,kappa_2,2)*
     &                      wave_new_even(state_2,i,kappa_1)*
     &                      wave_new_even(state_1,j,kappa_2)
     &                      )*
     &                      MatEl_Ang_comp_sml
                          endif
                        else
                          if(dabs(MatEl_Ang_comp_lge).gt. 1.d-12) then
                            Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
     &                      vmat(i,j,L,nvmat)*
     &                      wave_new_even(state_1,i,kappa_1)*
     &                      wave_new_even(state_2,j,kappa_2)*
     &                      MatEl_Ang_comp_lge
                          endif
c Notice here that the small components are logged in wave starting from
c nm+1.
                          if(dabs(MatEl_Ang_comp_sml).gt. 1.d-12) then
                            Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
     &                      vmat(i,j,L,nvmat)*
     &                      wave_new_even(state_1,i+nm,kappa_1)*
     &                      wave_new_even(state_2,j+nm,kappa_2)*
     &                      MatEl_Ang_comp_sml
                          endif
                        endif
                      enddo
                    enddo
                  endif
                enddo
              endif
              mm0=mm0+(
     &        Cpld_chnnl_elem_lge+
     &        Cpld_chnnl_elem_sml)
            enddo
          enddo
          if (state_1 .ne. state_2) then
c       Off Diagonal elements of mm
            mmeven(ii,jj)=-(0.d0,1.d0)*dRdXi*
     &      mm0/(eigval_e(state_2)-eigval_e(state_1))
          endif
        enddo
      enddo
      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=states_in_cont_e(ii)
        mmeven(ii,ii)=dTdXi*eigval_e(state_1)
        do jj=Start_state,Start_state+Mat_dimension-1
          if(ii.gt. jj)then
            mmeven(ii,jj)=-mmeven(jj,ii)
          endif
c            if(dkb) write(23,*)ii,jj,mmeven(ii,jj),eigval_e(jj)
c            if(.not.dkb) write(24,*)ii,jj,mmeven(ii,jj),eigval_e(jj)
        enddo
      enddo
ccc if(ii_xi.gt. xi_stepslower+2)then
ccc deallocate(Target_evals_diff_e)
ccc endif
      return
      end

C     ******************************************************************
C     ORPHANED VERSION
C     ******************************************************************
      subroutine MatrixMMGenerator_eqZodd_old(Mat_dimension,Start_state,
     &dTdXi,dRdXi,eigval_o,nkap,vmat,wave_new_odd,nsto,nm,nvmat,
     &mmodd,dvdRmatdkb1,dvdRmatdkb2,states_in_cont_o)
      include 'inc.par'
      integer Mat_dimension,Start_state,state_1,
     &state_2,nsto,nkap,nm,nvmat,
     &states_in_cont_o(Start_state:Start_state+Mat_dimension-1)
      real*8 dTdXi,eigval_o(nsto),mm0,Cpld_chnnl_elem_lge,dRdXi,
     &Cpld_chnnl_elem_sml,MatEl_Ang_comp_lge,MatEl_Ang_comp_sml,
     &vmat(nm,nm,0:2*nkap,nvmat),
     &wave_new_odd(nsto,2*nm,-nkap:nkap)
      real*8 dvdRmatdkb2(nm,nm,0:2*nkap,-nkap:nkap,-nkap:nkap,2)
      real*8 dvdRmatdkb1(nm,nm,0:2*nkap,-nkap:nkap,2)
      complex*16 mmodd(Start_state:Start_state+2*Mat_dimension-1,
     &Start_state:Start_state+2*Mat_dimension-1)
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
c      common /staterangeinput/ state_range_input
      logical dkb!,state_range_input

      mmodd=0.d0
      do ii=Start_state,Start_state+Mat_dimension-1!End_state
        state_1=states_in_cont_o(ii)
        do jj=Start_state,Start_state+Mat_dimension-1!End_state
          state_2=states_in_cont_o(jj)
          mm0=0.d0
          do kappa_1=-nkap,nkap
            do kappa_2=-nkap,nkap
              Cpld_chnnl_elem_lge=0.d0
              Cpld_chnnl_elem_sml=0.d0
              if (kappa_1*kappa_2.ne. 0) then
                do L=0,2*nkap,2
                  if (state_1.ge.state_2) then
                    Cpld_chnnl_elem_lge=0.d0
                    Cpld_chnnl_elem_sml=0.d0
                  else
                    MatEl_Ang_comp_lge=
     &              Coeff_maker(amu,kappa_1*1.d0,kappa_2*1.d0,L)
                    MatEl_Ang_comp_sml=
     &              Coeff_maker(amu,-1.d0*kappa_1,-1.d0*kappa_2,L)
                    do i=1,nm
                      do j=max(1,i-ns+1),min(i+ns-1,nm)!1,nm

                        if(dkb) then
                          if(dabs(MatEl_Ang_comp_lge).gt. 1.d-12) then
                            Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
     &                      (vmat(i,j,L,nvmat)*
     &                      wave_new_odd(state_1,i,kappa_1)*
     &                      wave_new_odd(state_2,j,kappa_2)+
     &                      dvdRmatdkb1(j,i,L,kappa_2,1)*
     &                      wave_new_odd(state_1,i,kappa_1)*
     &                      wave_new_odd(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb1(i,j,L,kappa_1,1)*
     &                      wave_new_odd(state_1,i+nm,kappa_1)*
     &                      wave_new_odd(state_2,j,kappa_2)+
     &                      dvdRmatdkb2(i,j,L,kappa_1,kappa_2,1)*
     &                      wave_new_odd(state_2,j+nm,kappa_2)*
     &                      wave_new_odd(state_1,i+nm,kappa_1)
     &                      )*
     &                      MatEl_Ang_comp_lge
                          endif
c Notice here that the small components are logged in wave starting from
c nm+1.
                          if(dabs(MatEl_Ang_comp_sml).gt. 1.d-12) then
                            Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
     &                      (vmat(i,j,L,nvmat)*
     &                      wave_new_odd(state_1,i+nm,kappa_1)*
     &                      wave_new_odd(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb1(j,i,L,kappa_2,2)*
     &                      wave_new_odd(state_1,j,kappa_2)*
     &                      wave_new_odd(state_2,i+nm,kappa_1)+
     &                      dvdRmatdkb1(i,j,L,kappa_1,2)*
     &                      wave_new_odd(state_2,i,kappa_1)*
     &                      wave_new_odd(state_2,j+nm,kappa_2)+
     &                      dvdRmatdkb2(i,j,L,kappa_1,kappa_2,2)*
     &                      wave_new_odd(state_2,i,kappa_1)*
     &                      wave_new_odd(state_1,j,kappa_2)
     &                      )*
     &                      MatEl_Ang_comp_sml
                          endif
                        else
                          if(dabs(MatEl_Ang_comp_lge).gt. 1.d-12) then
                            Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
     &                      vmat(i,j,L,nvmat)*
     &                      wave_new_odd(state_1,i,kappa_1)*
     &                      wave_new_odd(state_2,j,kappa_2)*
     &                      MatEl_Ang_comp_lge
                          endif
c Notice here that the small components are logged in wave starting from
c nm+1.
                          if(dabs(MatEl_Ang_comp_sml).gt. 1.d-12) then
                            Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
     &                      vmat(i,j,L,nvmat)*
     &                      wave_new_odd(state_1,i+nm,kappa_1)*
     &                      wave_new_odd(state_2,j+nm,kappa_2)*
     &                      MatEl_Ang_comp_sml
                          endif
                        endif
                      enddo
                    enddo
                  endif
                enddo
              endif
              mm0=mm0+(
     &        Cpld_chnnl_elem_lge+
     &        Cpld_chnnl_elem_sml)
            enddo
          enddo
          if (state_1 .ne. state_2) then
c       Off Diagonal elements of mm
            mmodd(ii,jj)=-(0.d0,1.d0)*dRdXi*mm0/
     &      (eigval_o(state_2)-eigval_o(state_1))
          endif
        enddo
        mmodd(ii,ii)=dTdXi*eigval_o(state_1)
      enddo
      do ii=Start_state,Start_state+Mat_dimension-1
        state_1=states_in_cont_o(ii)
        mmodd(ii,ii)=dTdXi*eigval_o(state_1)
        do jj=Start_state,Start_state+Mat_dimension-1
        if(ii.gt. jj)then
          mmodd(ii,jj)=-mmodd(jj,ii)
        endif
        if(dkb)write(25,*)ii,jj,mmodd(ii,jj),eigval_o(jj)
        if(.not.dkb)write(26,*)ii,jj,mmodd(ii,jj),eigval_o(jj)
      enddo
      enddo
      stop
ccc if(ii_xi.gt. xi_stepslower+2)then
ccc deallocate(Target_evals_diff_o)
ccc endif
      return
      end

      subroutine prj_neqZ(number_states,nstates,numb_mj,
     &nstates_mj,nkap,mm,eigval,tcb_eval,
     &eigvec,i_xi, xi_stepslower)

      include 'inc.par'
      integer number_states(2,-nkap:nkap),numb_mj(4,nkap),matelm(2),
     &xi_stepslower
      real*8 mm(nstates_mj,nstates_mj),tcb_eval(nstates_mj),
     &eigval(nstates),eigvec(nstates,nstates)
      real*8, dimension(:,:), allocatable:: mm0,temp_cont_states,
     &primary_element,primary_element_prev
      real*8, dimension(:,:,:), allocatable:: contributor_states
      integer, dimension (:), allocatable:: lsize,i_holder
      integer, dimension (:,:), allocatable:: contributor_substates
      logical dkb,check!,mkdirs
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max

      !mkdirs=makedirqq('PRJ')
      CALL SYSTEM("mkdir -p PRJ")

      mj_max=nint(amj_max+0.5d0)

      if(i_xi.gt.xi_stepslower+2)then
        allocate(primary_element_prev(nstates,nstates))
        open(44,file='PRJ/proj_mat_a_monopole_prev.dat',status='old')
        read(44,*)primary_element_prev
        close(44)
      endif
      if(i_xi.gt.xi_stepslower)then
        allocate(primary_element(nstates,nstates))
        open(22,file='PRJ/proj_mat_a_monopole.dat',status='old')
        read(22,*)primary_element
        close(22)
        open(44,file='PRJ/proj_mat_a_monopole_prev.dat')
        write(44,*)primary_element
        close(44)
      endif

      allocate(mm0(nstates,nstates))
      mm=0.d0
      mm0=0.d0
      do i=1,nstates
        do j=1,nstates
          mm0(i,j)=eigvec(j,i)
        enddo
      enddo
c    do ii=1,nstates
c      j1=1
c         do kk=-nkap,nkap
c            if(kk.ne. 0)then
c               do jj=number_states(1,kk),number_states(2,kk)
c               Cpld_chnnl_elem_lge=0.d0
c               Cpld_chnnl_elem_sml=0.d0
c                  do i=1,nm
c                 do j=max(1,i-ns+1),min(i+ns-1,nm)
c                     if(dkb)then
c                 Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
c     &           wave_new(ii,i,kk)*wave(j,jj,kk)*
c     &               alternate_dmat(i,j,kk)
c                   Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
c     &            wave_new(ii,i+nm,kk)*wave(j+nm,jj,kk)*
c     &             alternate_dmat(i+nm,j+nm,kk)
c                     else
c                 Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
c     &            wave_new(ii,i,kk)*wave(j,jj,kk)*
c     &               dmat(i,j)
c                   Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
c     &           wave_new(ii,i+nm,kk)*wave(j+nm,jj,kk)*
c     &             dmat(i+nm,j+nm)
c                     endif
c                  enddo
c               enddo
c            mm0(ii,j1)=Cpld_chnnl_elem_lge+Cpld_chnnl_elem_sml
c            j1=j1+1
c            enddo
c         endif
c         enddo
c      enddo

      iz=1
      do while(eigval(iz).lt.-1.d0)
        iz=iz+1
      enddo
      if(i_xi.gt.xi_stepslower)then
        do i=1,nstates
          if(dabs(mm0(iz,i)+primary_element(iz,i)).lt.
     &      dabs(mm0(iz,i)-primary_element(iz,i)))then
            do j=1,nstates
              mm0(j,i)=-mm0(j,i)
            enddo
          endif
        enddo
        deallocate(primary_element)
      elseif(i_xi.gt.xi_stepslower+2)then
        deallocate(primary_element_prev)
      endif

      open(22,file='PRJ/proj_mat_a_monopole.dat')
      write(22,*)mm0
      close(22)

      allocate(lsize(-nkap:nkap))
      lsize=0
      do k=-nkap,nkap
        if(k.ne. 0)then
          lsize(k)=-number_states(1,k)+number_states(2,k)+1
        endif
      enddo
      lsize(0)=maxval(lsize)
      allocate(contributor_states(-nkap:nkap,nstates,lsize(0)))
      lstart=1
      do k=-nkap,nkap
        if(k.ne.0)then
          do i=1,nstates
            jl=lstart
            jk=1
            do j=number_states(1,k),number_states(2,k)
              contributor_states(k,i,jk)=dabs(mm0(i,jl))
              jl=jl+1
              jk=jk+1
            enddo
          enddo
          lstart=jl
        endif
      enddo

      allocate(contributor_substates(-nkap:nkap,lsize(0)))
      contributor_substates=0

      do k=-nkap,nkap
        allocate(temp_cont_states(nstates,lsize(0)))
        temp_cont_states=0
        if(k.ne.0)then
          do i=1,nstates
            jk=1
            do j=number_states(1,k),number_states(2,k)
              temp_cont_states(i,jk)=contributor_states(k,i,jk)
              jk=jk+1
            enddo
          enddo
          js=1
          do while(js.le.lsize(k))
            matelm=maxloc(temp_cont_states)
            js2=1
            check=.true.
            do while (js2.le. lsize(k) .and. check)
              if(contributor_substates(k,js2)-matelm(1).eq. 0)then
                check=.false.
              endif
              js2=js2+1
            enddo
            if(check)then
              contributor_substates(k,js)=matelm(1)
              js=js+1
            endif
            temp_cont_states(matelm(1),matelm(2))=0.d0
          enddo
          allocate(i_holder(lsize(0)))
          do js=1,lsize(k)
            i_holder(js)=contributor_substates(k,js)
          enddo
          do js=1,lsize(k)
            contributor_substates(k,js)=
     &      minval(i_holder,mask=i_holder.gt.0)
            i_holder(minloc(i_holder,mask=i_holder.gt.0))=-1
          enddo
          deallocate(i_holder)
        endif
        deallocate(temp_cont_states)
      enddo

      deallocate(contributor_states)

      ii=1
      lstart=1
      do mj=-mj_max,mj_max
        if(mj.ne.0)then
          i3=1
          k=-nkap
          do i2=numb_mj(1,abs(mj)),numb_mj(2,abs(mj))
            if(i3.gt. lsize(k))then
              i3=1
              k=k+1
            endif
            if(abs(mj).gt. 1)then
              i=contributor_substates(mj,i3)
            else
              i=i2
            endif
            i3=i3+1
            jj=lstart
            do j=numb_mj(1,abs(mj)),numb_mj(2,abs(mj))
              mm(ii,jj)=mm0(i,j)
              jj=jj+1
            enddo
            tcb_eval(ii)=eigval(i)
            ii=ii+1
          enddo
          lstart=jj
          i3=1
          k=1
          do i2=numb_mj(3,abs(mj)),numb_mj(4,abs(mj))
            if(i3.gt. lsize(k))then
              i3=1
              k=k+1
            endif
            if(abs(mj).gt. 1)then
              i=contributor_substates(mj,i3)
            else
              i=i2
            endif
            i3=i3+1
            jj=lstart
            do j=numb_mj(3,abs(mj)),numb_mj(4,abs(mj))
              mm(ii,jj)=mm0(i,j)
              jj=jj+1
            enddo
            tcb_eval(ii)=eigval(i)
            ii=ii+1
          enddo
          lstart=jj
        endif
      enddo

      i1=1
      do mj=-mj_max,mj_max
        if(mj.ne. 0)then
          i2=1
          do k=-nkap,nkap
            if(abs(k).ge. abs(mj))then
              jj1=1
              do jj=number_states(1,k),number_states(2,k)
                if(abs(mj).gt. 1)then
                  i3=contributor_substates(k,jj1)
                else
                  i3=i2
                endif
                jj1=jj1+1
                tcb_eval(i1)=eigval(i3)
                i1=i1+1
                i2=i2+1
              enddo
            endif
          enddo
        endif
      enddo

      deallocate(mm0)
      deallocate(lsize)
      deallocate(contributor_substates)

      return
      end

      subroutine prj_eqZ_e(number_states,nstates,numb_mj,
     &nstates_mj,nkap,mm,eigval,tcb_eval,eigvec,i_xi, xi_stepslower)

      include 'inc.par'
      integer, dimension (:), allocatable:: lsize,i_holder
      integer, dimension (:,:), allocatable:: contributor_substates
      real*8, dimension (:,:,:), allocatable:: contributor_states
      integer number_states(2,-nkap:nkap),numb_mj(2,nkap),matelm(2),
     &xi_stepslower
      real*8 mm(nstates_mj,nstates_mj),eigval(nstates),
     &tcb_eval(nstates_mj),eigvec(nstates,nstates)
      real*8, dimension(:,:), allocatable:: mm0,temp_cont_states,
     &primary_element,primary_element_prev
      logical dkb,check!,mkdirs
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max

      !mkdirs=makedirqq('PRJ')
      CALL SYSTEM("mkdir -p PRJ")

      mj_max=nint(amj_max+0.5d0)

      if(i_xi.gt.xi_stepslower+2)then
        allocate(primary_element_prev(nstates,nstates))
        open(44,file='PRJ/proj_mat_e_monopole_prev.dat',status='old')
        read(44,*)primary_element_prev
        close(44)
      endif
      if(i_xi.gt.xi_stepslower)then
        allocate(primary_element(nstates,nstates))
        open(22,file='PRJ/proj_mat_e_monopole.dat',status='old')
        read(22,*)primary_element
        close(22)
        open(44,file='PRJ/proj_mat_e_monopole_prev.dat')
        write(44,*)primary_element
        close(44)
      endif

      allocate(mm0(nstates,nstates))
      mm=0.d0
      mm0=0.d0
      do i=1,nstates
        do j=1,nstates
          mm0(i,j)=eigvec(j,i)
        enddo
      enddo

c    do ii=1,nstates
c      j1=1
c      kfact=(-1)**nkap
c         do k=nkap,1,-1
c            kk=k*kfact
c               do jj=number_states(1,kk),number_states(2,kk)
c               Cpld_chnnl_elem_lge=0.d0
c               Cpld_chnnl_elem_sml=0.d0
c                  do i=1,nm
c                 do j=max(1,i-ns+1),min(i+ns-1,nm)
c                     if(dkb)then
c                 Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
c     &            wave_new(ii,i,kk)*wave(j,jj,kk)*
c     &               alternate_dmat(i,j,kk)
c                   Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
c     &            wave_new(ii,i+nm,kk)*wave(j+nm,jj,kk)*
c     &             alternate_dmat(i+nm,j+nm,kk)
c                     else
c                 Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
c     &            wave_new(ii,i,kk)*wave(j,jj,kk)*
c     &               dmat(i,j)
c                   Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
c     &            wave_new(ii,i+nm,kk)*wave(j+nm,jj,kk)*
c     &             dmat(i+nm,j+nm)
c                     endif
c                  enddo
c               enddo
c            mm0(ii,j1)=Cpld_chnnl_elem_lge+Cpld_chnnl_elem_sml
c            j1=j1+1
c            enddo
c         kfact=-kfact
c         enddo
c      enddo

      iz=1
      do while(eigval(iz).lt.-1.d0)
        iz=iz+1
      enddo
      if(i_xi.gt.xi_stepslower)then
        do i=1,nstates
          if(dabs(mm0(iz,i)+primary_element(iz,i)).lt.
     &    dabs(mm0(iz,i)-primary_element(iz,i)))then
            do j=1,nstates
              mm0(j,i)=-mm0(j,i)
            enddo
          endif
        enddo
        deallocate(primary_element)
      elseif(i_xi.gt.xi_stepslower+2)then
        deallocate(primary_element_prev)
      endif

      open(22,file='PRJ/proj_mat_e_monopole.dat')
      write(22,*)mm0
      close(22)

      allocate(lsize(-nkap:nkap))
      lsize=0
      do k=-nkap,nkap
        if(k.ne. 0)then
          lsize(k)=-number_states(1,k)+number_states(2,k)+1
        endif
      enddo
      lsize(0)=maxval(lsize)
      allocate(contributor_states(nkap,nstates,lsize(0)))
      contributor_states=0
      kfact=(-1)**nkap
      lstart=1
      do kk=nkap,1,-1
        k=kk*kfact
        do i=1,nstates
          jl=lstart
          jk=1
            do j=number_states(1,k),number_states(2,k)
              contributor_states(kk,i,jk)=dabs(mm0(i,jl))
              jl=jl+1
              jk=jk+1
            enddo
         enddo
         lstart=jl
         kfact=-kfact
      enddo
      allocate(contributor_substates(nkap,lsize(0)))
      contributor_substates=0

      kfact=(-1)**nkap
      do kk=nkap,1,-1
        allocate(temp_cont_states(nstates,lsize(0)))
        temp_cont_states=0
        k=kk*kfact
        do i=1,nstates
          jk=1
          do j=number_states(1,k),number_states(2,k)
            temp_cont_states(i,jk)=contributor_states(kk,i,jk)
            jk=jk+1
          enddo
        enddo

        js=1
        do while(js.le.lsize(k))
          matelm=maxloc(temp_cont_states)
          js2=1
          check=.true.
          do while (js2.le. lsize(k) .and. check)
            if(contributor_substates(kk,js2)-matelm(1).eq. 0)then
              check=.false.
            endif
              js2=js2+1
            enddo
            if(check)then
              contributor_substates(kk,js)=matelm(1)
              js=js+1
            endif
            temp_cont_states(matelm(1),matelm(2))=0.d0
          enddo
          deallocate(temp_cont_states)
          allocate(i_holder(lsize(0)))
          do js=1,lsize(k)
            i_holder(js)=contributor_substates(kk,js)
          enddo
          do js=1,lsize(k)
            contributor_substates(kk,js)=minval(i_holder,
     &      mask=i_holder.gt.0)
            i_holder(minloc(i_holder,mask=i_holder.gt.0))=-1
          enddo
        deallocate(i_holder)
        kfact=-kfact
      enddo

      deallocate(contributor_states)

      ii=1
      lstart=1
      do mj=-mj_max,mj_max
        if(mj.ne.0)then
          i3=1
          kfact=(-1)**nkap
          k=nkap*kfact
          do i2=numb_mj(1,abs(mj)),numb_mj(2,abs(mj))
            if(i3.gt. lsize(k))then
              i3=1
              kfact=-kfact
              k=(abs(k)-1)*kfact
            endif
            if(abs(mj).gt. 1)then
              i=contributor_substates(abs(mj),i3)
            else
              i=i2
            endif
            i3=i3+1
            jj=lstart
            do j=numb_mj(1,abs(mj)),numb_mj(2,abs(mj))
              mm(ii,jj)=mm0(i,j)
              jj=jj+1
            enddo
            ii=ii+1
          enddo
          lstart=jj
        endif
      enddo

      i1=1
      do mj=-mj_max,mj_max
        if(mj.ne. 0)then
          i2=1
          kfact=(-1)**nkap
          do kk=nkap,1,-1
            k=kk*kfact
            if(abs(k).ge. abs(mj))then
              jj1=1
              do jj=number_states(1,k),number_states(2,k)
                if(abs(mj).gt. 1)then
                  i3=contributor_substates(abs(k),jj1)
                else
                  i3=i2
                endif
                jj1=jj1+1
                tcb_eval(i1)=eigval(i3)
                i1=i1+1
                i2=i2+1
              enddo
            endif
            kfact=-kfact
          enddo
        endif
      enddo

      deallocate(mm0)
      deallocate(lsize)
      deallocate(contributor_substates)

      return
      end

      subroutine prj_eqZ_o(number_states,nstates,numb_mj,
     &nstates_mj,nkap,mm,eigval,tcb_eval,
     &eigvec,i_xi, xi_stepslower)

      include 'inc.par'
      integer number_states(2,-nkap:nkap),numb_mj(2,nkap),matelm(2),
     &xi_stepslower
      real*8 mm(nstates_mj,nstates_mj),eigval(nstates),
     &tcb_eval(nstates_mj),eigvec(nstates,nstates)
      real*8, dimension(:,:), allocatable:: mm0,temp_cont_states,
     &primary_element,primary_element_prev
      real*8, dimension(:,:,:), allocatable:: contributor_states
      integer, dimension(:), allocatable:: lsize,i_holder
      integer, dimension(:,:), allocatable:: contributor_substates
      logical dkb,check
      common /common_dkb/ dkb!,mkdirs
      common /momentum_projection/ amu,amj_max

      !mkdirs=makedirqq('PRJ')
      CALL SYSTEM("mkdir -p PRJ")

      mj_max=nint(amj_max+0.5d0)

      if(i_xi.gt.xi_stepslower+2)then
        allocate(primary_element_prev(nstates,nstates))
        open(44,file='PRJ/proj_mat_o_monopole_prev.dat',status='old')
        read(44,*)primary_element_prev
        close(44)
      endif
      if(i_xi.gt.xi_stepslower)then
        allocate(primary_element(nstates,nstates))
        open(22,file='PRJ/proj_mat_o_monopole.dat',status='old')
        read(22,*)primary_element
        close(22)
        open(44,file='PRJ/proj_mat_o_monopole_prev.dat')
        write(44,*)primary_element
        close(44)
      endif

      allocate(mm0(nstates,nstates))
      mm=0.d0
      mm0=0.d0
      do i=1,nstates
        do j=1,nstates
          mm0(i,j)=eigvec(j,i)
        enddo
      enddo
c    do ii=1,nstates
c      j1=1
c      kfact=(-1)**(nkap+1)
c         do k=nkap,1,-1
c            kk=k*kfact
c               do jj=number_states(1,kk),number_states(2,kk)
c               Cpld_chnnl_elem_lge=0.d0
c               Cpld_chnnl_elem_sml=0.d0
c                  do i=1,nm
c                 do j=max(1,i-ns+1),min(i+ns-1,nm)
c                     if(dkb)then
c                 Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
c     &            wave_new(ii,i,kk)*wave(j,jj,kk)*
c     &               alternate_dmat(i,j,kk)
c                   Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
c     &            wave_new(ii,i+nm,kk)*wave(j+nm,jj,kk)*
c     &             alternate_dmat(i+nm,j+nm,kk)
c                     else
c                 Cpld_chnnl_elem_lge=Cpld_chnnl_elem_lge+
c     &            wave_new(ii,i,kk)*wave(j,jj,kk)*
c     &               dmat(i,j)
c                   Cpld_chnnl_elem_sml=Cpld_chnnl_elem_sml+
c     &            wave_new(ii,i+nm,kk)*wave(j+nm,jj,kk)*
c     &             dmat(i+nm,j+nm)
c                     endif
c                     enddo
c                  enddo
c               mm0(ii,j1)=Cpld_chnnl_elem_lge+Cpld_chnnl_elem_sml
c               j1=j1+1
c            enddo
c         kfact=-kfact
c         enddo
c      enddo

      iz=1
      do while(eigval(iz).lt.-1.d0)
        iz=iz+1
      enddo
      if(i_xi.gt.xi_stepslower)then
        do i=1,nstates
          if(dabs(mm0(iz,i)+primary_element(iz,i)).lt.
     &    dabs(mm0(iz,i)-primary_element(iz,i)))then
            do j=1,nstates
              mm0(j,i)=-mm0(j,i)
            enddo
          endif
        enddo
        deallocate(primary_element)
      elseif(i_xi.gt.xi_stepslower+2)then
        deallocate(primary_element_prev)
      endif

      open(22,file='PRJ/proj_mat_o_monopole.dat')
      write(22,*)mm0
      close(22)

      allocate(lsize(-nkap:nkap))
      lsize=0
      do k=-nkap,nkap
        if(k.ne. 0)then
          lsize(k)=-number_states(1,k)+number_states(2,k)+1
        endif
      enddo
      lsize(0)=maxval(lsize)
      allocate(contributor_states(nkap,nstates,lsize(0)))
      kfact=(-1)**(nkap+ 1)
      lstart=1
      do kk=nkap,1,-1
        k=kk*kfact
        do i=1,nstates
          jl=lstart
          jk=1
          do j=number_states(1,k),number_states(2,k)
            contributor_states(kk,i,jk)=dabs(mm0(i,jl))
            jl=jl+1
            jk=jk+1
          enddo
        enddo
        lstart=jl
        kfact=-kfact
      enddo

      allocate(contributor_substates(nkap,lsize(0)))
      contributor_substates=0

      kfact=(-1)**(nkap+1)
      do kk=nkap,1,-1
        allocate(temp_cont_states(nstates,lsize(0)))
        temp_cont_states=0
        k=kk*kfact
        do i=1,nstates
          jk=1
          do j=number_states(1,k),number_states(2,k)
            temp_cont_states(i,jk)=contributor_states(kk,i,jk)
            jk=jk+1
          enddo
        enddo
        js=1
        do while(js.le.lsize(k))
          matelm=maxloc(temp_cont_states)
          js2=1
          check=.true.
          do while (js2.le. lsize(k) .and. check)
            if(contributor_substates(kk,js2)-matelm(1).eq. 0)then
              check=.false.
            endif
            js2=js2+1
          enddo
          if(check)then
            contributor_substates(kk,js)=matelm(1)
            js=js+1
          endif
          temp_cont_states(matelm(1),matelm(2))=0.d0
        enddo
        deallocate(temp_cont_states)
        allocate(i_holder(lsize(0)))
        do js=1,lsize(k)
          i_holder(js)=contributor_substates(kk,js)
        enddo
        do js=1,lsize(k)
          contributor_substates(kk,js)=minval(i_holder,
     &    mask=i_holder.gt.0)
          i_holder(minloc(i_holder,mask=i_holder.gt.0))=-1
        enddo
        deallocate(i_holder)
        kfact=-kfact
      enddo

      deallocate(contributor_states)

      ii=1
      lstart=1
      do mj=-mj_max,mj_max
        if(mj.ne.0)then
          i3=1
          kfact=(-1)**(nkap+1)
          k=nkap*kfact
          do i2=numb_mj(1,abs(mj)),numb_mj(2,abs(mj))
            if(i3.gt. lsize(k))then
              i3=1
              kfact=-kfact
              k=(abs(k)-1)*kfact
            endif
            if(abs(mj).gt. 1)then
              i=contributor_substates(abs(mj),i3)
            else
              i=i2
            endif
            i3=i3+1
            jj=lstart
            do j=numb_mj(1,abs(mj)),numb_mj(2,abs(mj))
              mm(ii,jj)=mm0(i,j)
              jj=jj+1
            enddo
            ii=ii+1
          enddo
          lstart=jj
        endif
      enddo

      i1=1
      do mj=-mj_max,mj_max
        if(mj.ne. 0)then
          i2=1
          kfact=(-1)**(nkap+1)
          do kk=nkap,1,-1
            k=kk*kfact
            if(abs(k).ge. abs(mj))then
              jj1=1
              do jj=number_states(1,k),number_states(2,k)
                if(abs(mj).gt. 1)then
                  i3=contributor_substates(abs(k),jj1)
                else
                  i3=i2
                endif
                jj1=jj1+1
                tcb_eval(i1)=eigval(i3)
                i1=i1+1
                i2=i2+1
              enddo
            endif
            kfact=-kfact
          enddo
        endif
      enddo

      deallocate(mm0)
      deallocate(contributor_substates)
      deallocate(lsize)

      return
      end

      SUBROUTINE CGAMA(X,Y,KF,GR,GI)
C
C       =========================================================
C       Purpose: Compute the gamma function â(z) or ln[â(z)]
C                for a complex argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                KF --- Function code
C                       KF=0 for ln[â(z)]
C                       KF=1 for â(z)
C       Output:  GR --- Real part of ln[â(z)] or â(z)
C                GI --- Imaginary part of ln[â(z)] or â(z)
C       ========================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(10)
      PI=3.141592653589793D0
      DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
      IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
         GR=1.0D+300
         GI=0.0D0
         RETURN
      ELSE IF (X.LT.0.0D0) THEN
         X1=X
         Y1=Y
         X=-X
         Y=-Y
      ENDIF
      X0=X
      IF (X.LE.7.0) THEN
         NA=INT(7-X)
         X0=X+NA
      ENDIF
      Z1=DSQRT(X0*X0+Y*Y)
      TH=DATAN(Y/X0)
      GR=(X0-.5D0)*DLOG(Z1)-TH*Y-X0+0.5D0*DLOG(2.0D0*PI)
      GI=TH*(X0-0.5D0)+Y*DLOG(Z1)-Y
      DO K=1,10
         T=Z1**(1-2*K)
         GR=GR+A(K)*T*DCOS((2.0D0*K-1.0D0)*TH)
         GI=GI-A(K)*T*DSIN((2.0D0*K-1.0D0)*TH)
      ENDDO
      IF (X.LE.7.0) THEN
         GR1=0.0D0
         GI1=0.0D0
         DO J=0,NA-1
            GR1=GR1+.5D0*DLOG((X+J)**2+Y*Y)
            GI1=GI1+DATAN(Y/(X+J))
         ENDDO
         GR=GR-GR1
         GI=GI-GI1
      ENDIF
      IF (X1.LT.0.0D0) THEN
         Z1=DSQRT(X*X+Y*Y)
         TH1=DATAN(Y/X)
         SR=-DSIN(PI*X)*DCOSH(PI*Y)
         SI=-DCOS(PI*X)*DSINH(PI*Y)
         Z2=DSQRT(SR*SR+SI*SI)
         TH2=DATAN(SI/SR)
         IF (SR.LT.0.0D0) TH2=PI+TH2
         GR=DLOG(PI/(Z1*Z2))-GR
         GI=-TH1-TH2-GI
         X=X1
         Y=Y1
      ENDIF
      IF (KF.EQ.1) THEN
         G0=DEXP(GR)
         GR=G0*DCOS(GI)
         GI=G0*DSIN(GI)
      ENDIF
      RETURN
      END


      subroutine Ave_Potential_at_edge(v,R,z1,z2,rr)
      include 'inc.par'

      EllKFaktor1=-4.d0*rr*R*z2*(z1+z2)/((-rr*z2+R*z2-rr*z1)**2)
      EllKFaktor2=4.d0*rr*R*z1*(z1+z2)/((rr*z2+R*z1+rr*z1)**2)
      VorFaktor1=dabs((-rr*z1-rr*z2+R*z2)/(z1+z2))
      VorFaktor2=dabs((rr*z1+rr*z2+R*z1)/(z1+z2))

      call COMELP(EllKFaktor1,Ell1K,Ell1E)
      call COMELP(EllKFaktor2,Ell2K,Ell2E)

      v=(-1.d0/c/pi)*(2.d0*z1*Ell1K/VorFaktor1 +
     & 2.d0*z2*Ell2K/VorFaktor2)

      return
      end

      SUBROUTINE COMELP(HK,CK,CE)
C
C       ==================================================
C       Purpose: Compute complete elliptic integrals K(k)
C                and E(k)
C       Input  : K  --- Modulus k ( 0 ó k ó 1 )
C       Output : CK --- K(k)
C                CE --- E(k)
C       ==================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PK=1.0D0-HK*HK
      IF (HK.EQ.1.0) THEN
         CK=1.0D+300
         CE=1.0D0
      ELSE
         AK=(((.01451196212D0*PK+.03742563713D0)*PK
     &        +.03590092383D0)*PK+.09666344259D0)*PK+
     &        1.38629436112D0
         BK=(((.00441787012D0*PK+.03328355346D0)*PK+
     &        .06880248576D0)*PK+.12498593597D0)*PK+.5D0
         CK=AK-BK*DLOG(PK)
         AE=(((.01736506451D0*PK+.04757383546D0)*PK+
     &        .0626060122D0)*PK+.44325141463D0)*PK+1.0D0
         BE=(((.00526449639D0*PK+.04069697526D0)*PK+
     &        .09200180037D0)*PK+.2499836831D0)*PK
         CE=AE-BE*DLOG(PK)
      ENDIF
      RETURN
      END

      subroutine xc_am05_labertw(z,val)
! !INPUT/OUTPUT PARAMETERS:
!   z   : function argument (in,real)
!   val : value of lambert W function of z (out,real)
! !DESCRIPTION:
!   Lambert $W$-function using the method of Corless, Gonnet, Hare, Jeffrey and
!   Knuth, {\it Adv. Comp. Math.} {\bf 5}, 329 (1996). The approach is based
!   loosely on that in GNU Octave by N. N. Schraudolph, but this implementation
!   is for real values and the principal branch only.
!
! !REVISION HISTORY:
!   Created April 2005 (RAR)
!EOP
!BOC
      implicit none
! arguments
      real(8), intent(in) :: z
      real(8), intent(out) :: val
!  local variables
      real(8) e,t,p
      integer i
! if z too low, go with the first term of the power expansion, z
      if (z.lt.1.d-20) then
        val=z
        return
      end if
      e=exp(1.d0)
! inital guess
      if (abs(z+1.d0/e).gt.1.45d0) then
! asymptotic expansion at 0 and Inf
        val=log(z)
        val=val-log(val)
      else
! series expansion about -1/e to first order
        val=1.d0*sqrt(2.d0*e*z+2.d0)-1.d0
      end if
! find val through iteration
      do i=1,10
        p=exp(val)
        t=val*p-z
        if (val.ne.-1.d0) then
          t=t/(p*(val+1.d0)-0.5d0*(val+2.d0)*t/(val+1.d0))
        else
          t=0.d0
        end if
        val=val-t
        if (abs(t).lt.(2.48d0*1.d-14)*(1.d0+abs(val))) return
      end do
! this should never happen!
      write(*,*)
      write(*,*)'("Error(xc_am05_labertw): iteration limit reached")'
      write(*,*)'(" Likely cause: improper numbers (INFs, NaNs)'
      write(*,*)'in density")'
      write(*,*)
      stop
      end subroutine

      SUBROUTINE ROOT(A,FA,B,FB,U,FU,KOUNT,IFLAG,IERROR,EPMACH)
C
C***********************************************************************
C
C  ROOT SEEKS A ROOT OF THE EQUATION F(X)=0.0,
C  GIVEN A STARTING INTERVAL (A,B) ON WHICH F CHANGES SIGN.
C  ON FIRST CALL TO ROOT, THE INTERVAL AND FUNCTION VALUES FA AND FB
C  ARE INPUT AND AN APPROXIMATION U FOR THE ROOT IS RETURNED.
C  BEFORE EACH SUBSEQUENT CALL, THE USER EVALUATES FU=F(U), AND THE
C  PROGRAM TRIES TO RETURN A BETTER APPROXIMATION U.
C
C  SEE THE BOOK BY BRENT LISTED in the documentation.
C
C  VARIABLES
C
C  A     = IS ONE ENDPOINT OF AN INTERVAL IN WHICH F CHANGES SIGN.
C  FA    = THE VALUE OF F(A).  THE USER MUST EVALUATE F(A) BEFORE
C          FIRST CALL ONLY.  THEREAFTER THE PROGRAM SETS FA.
C  B     = IS THE OTHER ENDPOINT OF THE INTERVAL IN WHICH
C          F CHANGES SIGN.  NOTE THAT THE PROGRAM WILL RETURN
C          IMMEDIATELY WITH AN ERROR FLAG IF FB*FA.GT.0.0.
C  FB    = THE VALUE OF F(B).  THE USER MUST EVALUATE F(B) BEFORE
C          FIRST CALL ONLY.  THERAFTER THE PROGRAM SETS FB.
C  U     = ON FIRST CALL, U SHOULD NOT BE SET BY THE USER.
C          ON SUBSEQUENT CALLS, U SHOULD NOT BE CHANGED
C          FROM ITS OUTPUT VALUE, THE CURRENT APPROXIMANT
C          TO THE ROOT.
C  FU    = ON FIRST CALL, FU SHOULD NOT BE SET BY THE USER.
C          THEREAFTER, THE USER SHOULD EVALUATE THE FUNCTION
C          AT THE OUTPUT VALUE U, AND RETURN FU=F(U).
C  KOUNT = A COUNTER FOR THE NUMBER OF CALLS TO ROOT.  KOUNT
C          SHOULD BE SET TO ZERO ON THE FIRST CALL FOR A GIVEN
C          ROOT PROBLEM.
C  IFLAG = PROGRAM RETURN FLAG
C          IFLAG=-1  THE CURRENT BRACKETING INTERVAL WITH
C                    ENDPOINTS  STORED IN A AND B IS LESS THAN
C                    4*EPMACH*ABS(U)+EPMACH
C                    HENCE U SHOULD BE ACCEPTED AS THE ROOT.
C                    THE FUNCTION VALUE F(U) IS STORED IN FU.
C          IFLAG= 0  THE INPUT VALUE FU IS EXACTLY ZERO, AND
C                    U SHOULD BE ACCEPTED AS THE ROOT.
C          IFLAG>0   THE CURRENT APPROXIMATION TO THE ROOT IS
C                    CONTAINED IN U.  IF A BETTER APPROXIMATION IS
C                    DESIRED, SET FU=F(U) AND CALL ROOT AGAIN.
C                    THE VALUE OF IFLAG INDICATES
C                    THE METHOD THAT WAS USED TO PRODUCE U.
C
C          IFLAG= 1  BISECTION WAS USED.
C          IFLAG= 2  LINEAR INTERPOLATION (SECANT METHOD).
C          IFLAG= 3  INVERSE QUADRATIC INTERPOLATION.
C
C  IERROR= GLOBAL ERROR FLAG.
C          IERROR=0 NO ERROR FOUND
C          IERROR=7  ON FIRST CALL, FA*FB.GT.0.0. AND HENCE
C                    THE GIVEN INTERVAL IS UNACCEPTABLE.
C  EPMACH= THE RELATIVE MACHINE PRECISION
C
      DOUBLE PRECISION EIGHT
      DOUBLE PRECISION HALF
      DOUBLE PRECISION ONE
      DOUBLE PRECISION ONEP5
      DOUBLE PRECISION TWO
      DOUBLE PRECISION ZERO
C
      PARAMETER (EIGHT=8.0)
      PARAMETER (HALF=0.5)
      PARAMETER (ONE=1.0)
      PARAMETER (ONEP5=1.5)
      PARAMETER (TWO=2.0)
      PARAMETER (ZERO=0.0)
C
      INTRINSIC ABS
      INTRINSIC SIGN
C
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION EPMACH
      DOUBLE PRECISION FA
      DOUBLE PRECISION FB
      DOUBLE PRECISION FU
      DOUBLE PRECISION HALFUB
      INTEGER   IERROR
      INTEGER   IFLAG
      INTEGER   KOUNT
      DOUBLE PRECISION P
      DOUBLE PRECISION Q
      DOUBLE PRECISION R
      DOUBLE PRECISION S
      DOUBLE PRECISION SDEL1
      DOUBLE PRECISION SDEL2
      DOUBLE PRECISION SDEL3
      DOUBLE PRECISION SDEL4
      DOUBLE PRECISION STEP
      DOUBLE PRECISION TOLER
      DOUBLE PRECISION U
C
C  SEGMENT 1   FIRST CALL HANDLED SPECIALLY.  DO BOOKKEEPING.
C
      IF(KOUNT.LE.0)THEN
        IF((FA.GT.ZERO.AND.FB.GT.ZERO)
     1  .OR.(FA.LT.ZERO.AND.FB.LT.ZERO))THEN
          IERROR=7
          KOUNT=0
          RETURN
          ENDIF
        KOUNT=1
        SDEL1=TWO*ABS(B-A)
        SDEL2=TWO*SDEL1
        SDEL3=TWO*SDEL2
        U=B
        B=A
        FU=FB
        FB=FA
      ELSE
C
C  INCREMENT COUNTER AND CHECK WHETHER F(U)=0.
C
        KOUNT=KOUNT+1
        IF(FU.EQ.ZERO)THEN
          IFLAG=0
          RETURN
          ENDIF
C
C  IF FU AND FB HAVE THE SAME SIGN OVERWRITE B WITH A
C
        IF(SIGN(ONE,FU).EQ.SIGN(ONE,FB))THEN
          B=A
          FB=FA
          ENDIF
        ENDIF
C
C  SEGMENT 2   REARRANGE POINTS AND FUNCTION VALUES IF
C  NECESSARY SO THAT  ABS(FU).LT.ABS(FB)
C
      IF(ABS(FB).LT.ABS(FU)) THEN
        A=U
        U=B
        B=A
        FA=FU
        FU=FB
        FB=FA
        ENDIF
C
C  SEGMENT 3   CHECK FOR ACCEPTANCE BECAUSE OF SMALL INTERVAL
C  CURRENT CHANGE IN SIGN INTERVAL IS (B,U) OR (U,B).
C
      TOLER=TWO*EPMACH*ABS(U)+EPMACH
      HALFUB=HALF*(B-U)
      SDEL4=SDEL3
      SDEL3=SDEL2
      SDEL2=SDEL1
      SDEL1=ABS(B-U)
      IF(ABS(HALFUB).LE.TOLER) THEN
        IFLAG=-1
        A=U
        FA=FU
        RETURN
        ENDIF
C
C  SEGMENT 4   COMPUTE NEW APPROXIMANT TO ROOT OF THE FORM
C  U(NEW)=U(OLD)+STEP.
C  METHODS AVAILABLE ARE LINEAR INTERPOLATION
C  INVERSE QUADRATIC INTERPOLATION
C  AND BISECTION.
C
      IF(ABS(FU).GE.ABS(FA)) THEN
        IFLAG=1
        STEP=HALFUB
        GO TO 80
        ENDIF
C
C  ATTEMPT LINEAR INTERPOLATION IF ONLY TWO POINTS AVAILABLE
C  COMPUTE P AND Q FOR APPROXIMATION U(NEW)=U(OLD)+P/Q
C
      IF(A.EQ.B) THEN
        IFLAG=2
        S=FU/FA
        P=TWO*HALFUB*S
        Q=ONE-S
C
C  ATTEMPT INVERSE QUADRATIC INTERPOLATION IF THREE POINTS AVAILABLE
C  COMPUTE P AND Q FOR APPROXIMATION U(NEW)=U(OLD)+P/Q
C
      ELSE
        IFLAG=3
        S=FU/FA
        Q=FA/FB
        R=FU/FB
        P=S*(TWO*HALFUB*Q*(Q-R)-(U-A)*(R-ONE))
        Q=(Q-ONE)*(R-ONE)*(S-ONE)
        ENDIF
C
C  CORRECT THE SIGNS OF P AND Q
C
      IF(P.GT.ZERO)Q=-Q
      P=ABS(P)
C
C  IF P/Q IS TOO LARGE, GO BACK TO BISECTION
C
      IF((EIGHT*SDEL1.GT.SDEL4)
     1 .OR.(P.GE.ONEP5*ABS(HALFUB*Q)-ABS(TOLER*Q))) THEN
        IFLAG=1
        STEP=HALFUB
        GO TO 80
        ENDIF
      STEP=P/Q
C
C  SEGMENT 5   VALUE OF STEP HAS BEEN COMPUTED.
C  UPDATE INFORMATION  A =U, FA =FU, U =U+STEP.
C  CHANGE IN SIGN INTERVAL IS NOW (A,B) OR (B,A).
C
80    CONTINUE
      A=U
      FA=FU
      IF(ABS(STEP).LE.TOLER) STEP=SIGN(TOLER,HALFUB)
      U=U+STEP
      RETURN
      END

      SUBROUTINE RPOLY(OP, DEGREE, ZEROR, ZEROI,
     *  FAIL)
C FINDS THE ZEROS OF A REAL POLYNOMIAL
C OP  - DOUBLE PRECISION VECTOR OF COEFFICIENTS IN
C       ORDER OF DECREASING POWERS.
C DEGREE   - INTEGER DEGREE OF POLYNOMIAL.
C ZEROR, ZEROI - OUTPUT DOUBLE PRECISION VECTORS OF
C                REAL AND IMAGINARY PARTS OF THE
C                ZEROS.
C FAIL  - OUTPUT LOGICAL PARAMETER, TRUE ONLY IF
C         LEADING COEFFICIENT IS ZERO OR IF RPOLY
C         HAS FOUND FEWER THAN DEGREE ZEROS.
C         IN THE LATTER CASE DEGREE IS RESET TO
C         THE NUMBER OF ZEROS FOUND.
C TO CHANGE THE SIZE OF POLYNOMIALS WHICH CAN BE
C SOLVED, RESET THE DIMENSIONS OF THE ARRAYS IN THE
C COMMON AREA AND IN THE FOLLOWING DECLARATIONS.
C THE SUBROUTINE USES SINGLE PRECISION CALCULATIONS
C FOR SCALING, BOUNDS AND ERROR CALCULATIONS. ALL
C CALCULATIONS FOR THE ITERATIONS ARE DONE IN DOUBLE
C PRECISION.
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      REAL ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION OP(101), TEMP(101),
     * ZEROR(100), ZEROI(100), T, AA, BB, CC, DABS,
     * FACTOR, L
      REAL PT(101), LO, MAX, MIN, XX, YY, COSR,
     * SINR, XXX, X, SC, BND, XM, FF, DF, DX, INFIN,
     * SMALNO, BASE
      INTEGER DEGREE, CNT, NZ, I, J, JJ, NM1
      LOGICAL FAIL, ZEROK
C THE FOLLOWING STATEMENTS SET MACHINE CONSTANTS USED
C IN VARIOUS PARTS OF THE PROGRAM. THE MEANING OF THE
C FOUR CONSTANTS ARE...
C ETA     THE MAXIMUM RELATIVE REPRESENTATION ERROR
C         WHICH CAN BE DESCRIBED AS THE SMALLEST
C         POSITIVE FLOATING POINT NUMBER SUCH THAT
C         1.D0+ETA IS GREATER THAN 1.
C INFINY  THE LARGEST FLOATING-POINT NUMBER.
C SMALNO  THE SMALLEST POSITIVE FLOATING-POINT NUMBER
C         IF THE EXPONENT RANGE DIFFERS IN SINGLE AND
C         DOUBLE PRECISION THEN SMALNO AND INFIN
C         SHOULD INDICATE THE SMALLER RANGE.
C BASE    THE BASE OF THE FLOATING-POINT NUMBER
C         SYSTEM USED.
C THE VALUES BELOW CORRESPOND TO THE BURROUGHS B6700
      BASE = 8.
      ETA = .5*BASE**(1-26)
      INFIN = HUGE(0.0)!4.3E68
      SMALNO = 1.0E-37
C ARE AND MRE REFER TO THE UNIT ERROR IN + AND *
C RESPECTIVELY. THEY ARE ASSUMED TO BE THE SAME AS
C ETA.
      ARE = ETA
      MRE = ETA
      LO = SMALNO/ETA
C INITIALIZATION OF CONSTANTS FOR SHIFT ROTATION
      XX = .70710678
      YY = -XX
      COSR = -.069756474
      SINR = .99756405
      FAIL = .FALSE.
      N = DEGREE
      NN = N + 1
C ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO.
      IF (OP(1).NE.0.D0) GO TO 10
      FAIL = .TRUE.
      DEGREE = 0
      RETURN
C REMOVE THE ZEROS AT THE ORIGIN IF ANY
   10 IF (OP(NN).NE.0.0D0) GO TO 20
      J = DEGREE - N + 1
      ZEROR(J) = 0.D0
      ZEROI(J) = 0.D0
      NN = NN - 1
      N = N - 1
      GO TO 10
C MAKE A COPY OF THE COEFFICIENTS
   20 DO 30 I=1,NN
        P(I) = OP(I)
   30 CONTINUE
C START THE ALGORITHM FOR ONE ZERO
   40 IF (N.GT.2) GO TO 60
      IF (N.LT.1) RETURN
C CALCULATE THE FINAL ZERO OR PAIR OF ZEROS
      IF (N.EQ.2) GO TO 50
      ZEROR(DEGREE) = -P(2)/P(1)
      ZEROI(DEGREE) = 0.0D0
      RETURN
   50 CALL QUAD(P(1), P(2), P(3), ZEROR(DEGREE-1),
     * ZEROI(DEGREE-1), ZEROR(DEGREE), ZEROI(DEGREE))
      RETURN
C FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
   60 MAX = 0.
      MIN = INFIN
      DO 70 I=1,NN
        X = ABS(SNGL(P(I)))
        IF (X.GT.MAX) MAX = X
        IF (X.NE.0. .AND. X.LT.MIN) MIN = X
   70 CONTINUE
C SCALE IF THERE ARE LARGE OR VERY SMALL COEFFICIENTS
C COMPUTES A SCALE FACTOR TO MULTIPLY THE
C COEFFICIENTS OF THE POLYNOMIAL. THE SCALING IS DONE
C TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
C INTERFERING WITH THE CONVERGENCE CRITERION.
C THE FACTOR IS A POWER OF THE BASE
      SC = LO/MIN
      IF (SC.GT.1.0) GO TO 80
      IF (MAX.LT.10.) GO TO 110
      IF (SC.EQ.0.) SC = SMALNO
      GO TO 90
   80 IF (INFIN/SC.LT.MAX) GO TO 110
   90 L = ALOG(SC)/ALOG(BASE) + .5
      FACTOR = (BASE*1.0D0)**L
      IF (FACTOR.EQ.1.D0) GO TO 110
      DO 100 I=1,NN
        P(I) = FACTOR*P(I)
  100 CONTINUE
C COMPUTE LOWER BOUND ON MODULI OF ZEROS.
  110 DO 120 I=1,NN
        PT(I) = ABS(SNGL(P(I)))
  120 CONTINUE
      PT(NN) = -PT(NN)
C COMPUTE UPPER ESTIMATE OF BOUND
      X = EXP((ALOG(-PT(NN))-ALOG(PT(1)))/FLOAT(N))
      IF (PT(N).EQ.0.) GO TO 130
C IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
      XM = -PT(NN)/PT(N)
      IF (XM.LT.X) X = XM
C CHOP THE INTERVAL (0,X) UNTIL FF .LE. 0
  130 XM = X*.1
      FF = PT(1)
      DO 140 I=2,NN
        FF = FF*XM + PT(I)
  140 CONTINUE
      IF (FF.LE.0.) GO TO 150
      X = XM
      GO TO 130
  150 DX = X
C DO NEWTON ITERATION UNTIL X CONVERGES TO TWO
C DECIMAL PLACES
  160 IF (ABS(DX/X).LE..005) GO TO 180
      FF = PT(1)
      DF = FF
      DO 170 I=2,N
        FF = FF*X + PT(I)
        DF = DF*X + FF
  170 CONTINUE
      FF = FF*X + PT(NN)
      DX = FF/DF
      X = X - DX
      GO TO 160
  180 BND = X
C COMPUTE THE DERIVATIVE AS THE INTIAL K POLYNOMIAL
C AND DO 5 STEPS WITH NO SHIFT
      NM1 = N - 1
      DO 190 I=2,N
        K(I) = FLOAT(NN-I)*P(I)/FLOAT(N)
  190 CONTINUE
      K(1) = P(1)
      AA = P(NN)
      BB = P(N)
      ZEROK = K(N).EQ.0.D0
      DO 230 JJ=1,5
        CC = K(N)
        IF (ZEROK) GO TO 210
C USE SCALED FORM OF RECURRENCE IF VALUE OF K AT 0 IS
C NONZERO
        T = -AA/CC
        DO 200 I=1,NM1
          J = NN - I
          K(J) = T*K(J-1) + P(J)
  200   CONTINUE
        K(1) = P(1)
        ZEROK = DABS(K(N)).LE.DABS(BB)*ETA*10.
        GO TO 230
C USE UNSCALED FORM OF RECURRENCE
  210   DO 220 I=1,NM1
          J = NN - I
          K(J) = K(J-1)
  220   CONTINUE
        K(1) = 0.D0
        ZEROK = K(N).EQ.0.D0
  230 CONTINUE
C SAVE K FOR RESTARTS WITH NEW SHIFTS
      DO 240 I=1,N
        TEMP(I) = K(I)
  240 CONTINUE
C LOOP TO SELECT THE QUADRATIC  CORRESPONDING TO EACH
C NEW SHIFT
      DO 280 CNT=1,20
C QUADRATIC CORRESPONDS TO A DOUBLE SHIFT TO A
C NON-REAL POINT AND ITS COMPLEX CONJUGATE. THE POINT
C HAS MODULUS BND AND AMPLITUDE ROTATED BY 94 DEGREES
C FROM THE PREVIOUS SHIFT
        XXX = COSR*XX - SINR*YY
        YY = SINR*XX + COSR*YY
        XX = XXX
        SR = BND*XX
        SI = BND*YY
        U = -2.0D0*SR
        V = BND
C SECOND STAGE CALCULATION, FIXED QUADRATIC
        CALL FXSHFR(20*CNT, NZ)
        IF (NZ.EQ.0) GO TO 260
C THE SECOND STAGE JUMPS DIRECTLY TO ONE OF THE THIRD
C STAGE ITERATIONS AND RETURNS HERE IF SUCCESSFUL.
C DEFLATE THE POLYNOMIAL, STORE THE ZERO OR ZEROS AND
C RETURN TO THE MAIN ALGORITHM.
        J = DEGREE - N + 1
        ZEROR(J) = SZR
        ZEROI(J) = SZI
        NN = NN - NZ
        N = NN - 1
        DO 250 I=1,NN
          P(I) = QP(I)
  250   CONTINUE
        IF (NZ.EQ.1) GO TO 40
        ZEROR(J+1) = LZR
        ZEROI(J+1) = LZI
        GO TO 40
C IF THE ITERATION IS UNSUCCESSFUL ANOTHER QUADRATIC
C IS CHOSEN AFTER RESTORING K
  260   DO 270 I=1,N
          K(I) = TEMP(I)
  270   CONTINUE
  280 CONTINUE
C RETURN WITH FAILURE IF NO CONVERGENCE WITH 20
C SHIFTS
      FAIL = .TRUE.
      DEGREE = DEGREE - N
      RETURN
      END

      SUBROUTINE FXSHFR(L2, NZ)
C COMPUTES UP TO  L2  FIXED SHIFT K-POLYNOMIALS,
C TESTING FOR CONVERGENCE IN THE LINEAR OR QUADRATIC
C CASE. INITIATES ONE OF THE VARIABLE SHIFT
C ITERATIONS AND RETURNS WITH THE NUMBER OF ZEROS
C FOUND.
C L2 - LIMIT OF FIXED SHIFT STEPS
C NZ - NUMBER OF ZEROS FOUND
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      REAL ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION SVU, SVV, UI, VI, S
      REAL BETAS, BETAV, OSS, OVV, SS, VV, TS, TV,
     * OTS, OTV, TVV, TSS
      INTEGER L2, NZ, TYPE, I, J, IFLAG
      LOGICAL VPASS, SPASS, VTRY, STRY
      NZ = 0
      BETAV = .25
      BETAS = .25
      OSS = REAL(SR)
      OVV = REAL(V)
C EVALUATE POLYNOMIAL BY SYNTHETIC DIVISION
      CALL QUADSD(NN, U, V, P, QP, A, B)
      CALL CALCSC(TYPE)
      DO 80 J=1,L2
C CALCULATE NEXT K POLYNOMIAL AND ESTIMATE V
        CALL NEXTK(TYPE)
        CALL CALCSC(TYPE)
        CALL NEWEST(TYPE, UI, VI)
        VV = REAL(VI)
C ESTIMATE S
        SS = 0.
        IF (K(N).NE.0.D0) SS = REAL(-P(NN)/K(N))
        TV = 1.
        TS = 1.
        IF (J.EQ.1 .OR. TYPE.EQ.3) GO TO 70
C COMPUTE RELATIVE MEASURES OF CONVERGENCE OF S AND V
C SEQUENCES
        IF (VV.NE.0.) TV = ABS((VV-OVV)/VV)
        IF (SS.NE.0.) TS = ABS((SS-OSS)/SS)
C IF DECREASING, MULTIPLY TWO MOST RECENT
C CONVERGENCE MEASURES
        TVV = 1.
        IF (TV.LT.OTV) TVV = TV*OTV
        TSS = 1.
        IF (TS.LT.OTS) TSS = TS*OTS
C COMPARE WITH CONVERGENCE CRITERIA
        VPASS = TVV.LT.BETAV
        SPASS = TSS.LT.BETAS
        IF (.NOT.(SPASS .OR. VPASS)) GO TO 70
C AT LEAST ONE SEQUENCE HAS PASSED THE CONVERGENCE
C TEST. STORE VARIABLES BEFORE ITERATING
        SVU = U
        SVV = V
        DO 10 I=1,N
          SVK(I) = K(I)
   10   CONTINUE
        S = SS
C CHOOSE ITERATION ACCORDING TO THE FASTEST
C CONVERGING SEQUENCE
        VTRY = .FALSE.
        STRY = .FALSE.
        IF (SPASS .AND. ((.NOT.VPASS) .OR.
     *   TSS.LT.TVV)) GO TO 40
   20   CALL QUADIT(UI, VI, NZ)
        IF (NZ.GT.0) RETURN
C QUADRATIC ITERATION HAS FAILED. FLAG THAT IT HAS
C BEEN TRIED AND DECREASE THE CONVERGENCE CRITERION.
        VTRY = .TRUE.
        BETAV = BETAV*.25
C TRY LINEAR ITERATION IF IT HAS NOT BEEN TRIED AND
C THE S SEQUENCE IS CONVERGING
        IF (STRY .OR. (.NOT.SPASS)) GO TO 50
        DO 30 I=1,N
          K(I) = SVK(I)
   30   CONTINUE
   40   CALL REALIT(S, NZ, IFLAG)
        IF (NZ.GT.0) RETURN
C LINEAR ITERATION HAS FAILED. FLAG THAT IT HAS BEEN
C TRIED AND DECREASE THE CONVERGENCE CRITERION
        STRY = .TRUE.
        BETAS = BETAS*.25
        IF (IFLAG.EQ.0) GO TO 50
C IF LINEAR ITERATION SIGNALS AN ALMOST DOUBLE REAL
C ZERO ATTEMPT QUADRATIC INTERATION
        UI = -(S+S)
        VI = S*S
        GO TO 20
C RESTORE VARIABLES
   50   U = SVU
        V = SVV
        DO 60 I=1,N
          K(I) = SVK(I)
   60   CONTINUE
C TRY QUADRATIC ITERATION IF IT HAS NOT BEEN TRIED
C AND THE V SEQUENCE IS CONVERGING
        IF (VPASS .AND. (.NOT.VTRY)) GO TO 20
C RECOMPUTE QP AND SCALAR VALUES TO CONTINUE THE
C SECOND STAGE
        CALL QUADSD(NN, U, V, P, QP, A, B)
        CALL CALCSC(TYPE)
   70   OVV = VV
        OSS = SS
        OTV = TV
        OTS = TS
   80 CONTINUE
      RETURN
      END
      SUBROUTINE QUADIT(UU, VV, NZ)
C VARIABLE-SHIFT K-POLYNOMIAL ITERATION FOR A
C QUADRATIC FACTOR CONVERGES ONLY IF THE ZEROS ARE
C EQUIMODULAR OR NEARLY SO.
C UU,VV - COEFFICIENTS OF STARTING QUADRATIC
C NZ - NUMBER OF ZERO FOUND
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      REAL ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION UI, VI, UU, VV, DABS
      REAL MP, OMP, EE, RELSTP, T, ZM
      INTEGER NZ, TYPE, I, J
      LOGICAL TRIED
      NZ = 0
      TRIED = .FALSE.
      U = UU
      V = VV
      J = 0
C MAIN LOOP
   10 CALL QUAD(1.D0, U, V, SZR, SZI, LZR, LZI)
C RETURN IF ROOTS OF THE QUADRATIC ARE REAL AND NOT
C CLOSE TO MULTIPLE OR NEARLY EQUAL AND  OF OPPOSITE
C SIGN
      IF (DABS(DABS(SZR)-DABS(LZR)).GT..01D0*
     * DABS(LZR)) RETURN
C EVALUATE POLYNOMIAL BY QUADRATIC SYNTHETIC DIVISION
      CALL QUADSD(NN, U, V, P, QP, A, B)
      MP = REAL(DABS(A-SZR*B) + DABS(SZI*B))
C COMPUTE A RIGOROUS  BOUND ON THE ROUNDING ERROR IN
C EVALUTING P
      ZM = SQRT(ABS(SNGL(V)))
      EE = 2.*ABS(SNGL(QP(1)))
      T = REAL(-SZR*B)
      DO 20 I=2,N
        EE = EE*ZM + ABS(SNGL(QP(I)))
   20 CONTINUE
      EE = EE*ZM + ABS(SNGL(A)+T)
      EE = (5.*MRE+4.*ARE)*EE - (5.*MRE+2.*ARE)*
     * (ABS(SNGL(A)+T)+ABS(SNGL(B))*ZM) +
     * 2.*ARE*ABS(T)
C ITERATION HAS CONVERGED SUFFICIENTLY IF THE
C POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND
      IF (MP.GT.20.*EE) GO TO 30
      NZ = 2
      RETURN
   30 J = J + 1
C STOP ITERATION AFTER 20 STEPS
      IF (J.GT.20) RETURN
      IF (J.LT.2) GO TO 50
      IF (RELSTP.GT..01 .OR. MP.LT.OMP .OR. TRIED)
     * GO TO 50
C A CLUSTER APPEARS TO BE STALLING THE CONVERGENCE.
C FIVE FIXED SHIFT STEPS ARE TAKEN WITH A U,V CLOSE
C TO THE CLUSTER
      IF (RELSTP.LT.ETA) RELSTP = ETA
      RELSTP = SQRT(RELSTP)
      U = U - U*RELSTP
      V = V + V*RELSTP
      CALL QUADSD(NN, U, V, P, QP, A, B)
      DO 40 I=1,5
        CALL CALCSC(TYPE)
        CALL NEXTK(TYPE)
   40 CONTINUE
      TRIED = .TRUE.
      J = 0
   50 OMP = MP
C CALCULATE NEXT K POLYNOMIAL AND NEW U AND V
      CALL CALCSC(TYPE)
      CALL NEXTK(TYPE)
      CALL CALCSC(TYPE)
      CALL NEWEST(TYPE, UI, VI)
C IF VI IS ZERO THE ITERATION IS NOT CONVERGING
      IF (VI.EQ.0.D0) RETURN
      RELSTP = REAL(DABS((VI-V)/VI))
      U = UI
      V = VI
      GO TO 10
      END
      SUBROUTINE REALIT(SSS, NZ, IFLAG)
C VARIABLE-SHIFT H POLYNOMIAL ITERATION FOR A REAL
C ZERO.
C SSS   - STARTING ITERATE
C NZ    - NUMBER OF ZERO FOUND
C IFLAG - FLAG TO INDICATE A PAIR OF ZEROS NEAR REAL
C         AXIS.
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      REAL ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION PV, KV, T, S, SSS, DABS
      REAL MS, MP, OMP, EE
      INTEGER NZ, IFLAG, I, J, NM1
      NM1 = N - 1
      NZ = 0
      S = SSS
      IFLAG = 0
      J = 0
C MAIN LOOP
   10 PV = P(1)
C EVALUATE P AT S
      QP(1) = PV
      DO 20 I=2,NN
        PV = PV*S + P(I)
        QP(I) = PV
   20 CONTINUE
      MP = REAL(DABS(PV))
C COMPUTE A RIGOROUS BOUND ON THE ERROR IN EVALUATING
C P
      MS = REAL(DABS(S))
      EE = (MRE/(ARE+MRE))*ABS(SNGL(QP(1)))
      DO 30 I=2,NN
        EE = EE*MS + ABS(SNGL(QP(I)))
   30 CONTINUE
C ITERATION HAS CONVERGED SUFFICIENTLY IF THE
C POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND
      IF (MP.GT.20.*((ARE+MRE)*EE-MRE*MP)) GO TO 40
      NZ = 1
      SZR = S
      SZI = 0.D0
      RETURN
   40 J = J + 1
C STOP ITERATION AFTER 10 STEPS
      IF (J.GT.10) RETURN
      IF (J.LT.2) GO TO 50
      IF (DABS(T).GT..001*DABS(S-T) .OR. MP.LE.OMP)
     * GO TO 50
C A CLUSTER OF ZEROS NEAR THE REAL AXIS HAS BEEN
C ENCOUNTERED RETURN WITH IFLAG SET TO INITIATE A
C QUADRATIC ITERATION
      IFLAG = 1
      SSS = S
      RETURN
C RETURN IF THE POLYNOMIAL VALUE HAS INCREASED
C SIGNIFICANTLY
   50 OMP = MP
C COMPUTE T, THE NEXT POLYNOMIAL, AND THE NEW ITERATE
      KV = K(1)
      QK(1) = KV
      DO 60 I=2,N
        KV = KV*S + K(I)
        QK(I) = KV
   60 CONTINUE
      IF (DABS(KV).LE.DABS(K(N))*10.*ETA) GO TO 80
C USE THE SCALED FORM OF THE RECURRENCE IF THE VALUE
C OF K AT S IS NONZERO
      T = -PV/KV
      K(1) = QP(1)
      DO 70 I=2,N
        K(I) = T*QK(I-1) + QP(I)
   70 CONTINUE
      GO TO 100
C USE UNSCALED FORM
   80 K(1) = 0.0D0
      DO 90 I=2,N
        K(I) = QK(I-1)
   90 CONTINUE
  100 KV = K(1)
      DO 110 I=2,N
        KV = KV*S + K(I)
  110 CONTINUE
      T = 0.D0
      IF (DABS(KV).GT.DABS(K(N))*10.*ETA) T = -PV/KV
      S = S + T
      GO TO 10
      END
      SUBROUTINE CALCSC(TYPE)
C THIS ROUTINE CALCULATES SCALAR QUANTITIES USED TO
C COMPUTE THE NEXT K POLYNOMIAL AND NEW ESTIMATES OF
C THE QUADRATIC COEFFICIENTS.
C TYPE - INTEGER VARIABLE SET HERE INDICATING HOW THE
C CALCULATIONS ARE NORMALIZED TO AVOID OVERFLOW
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      REAL ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION DABS
      INTEGER TYPE
C SYNTHETIC DIVISION OF K BY THE QUADRATIC 1,U,V
      CALL QUADSD(N, U, V, K, QK, C, D)
      IF (DABS(C).GT.DABS(K(N))*100.*ETA) GO TO 10
      IF (DABS(D).GT.DABS(K(N-1))*100.*ETA) GO TO 10
      TYPE = 3
C TYPE=3 INDICATES THE QUADRATIC IS ALMOST A FACTOR
C OF K
      RETURN
   10 IF (DABS(D).LT.DABS(C)) GO TO 20
      TYPE = 2
C TYPE=2 INDICATES THAT ALL FORMULAS ARE DIVIDED BY D
      E = A/D
      F = C/D
      G = U*B
      H = V*B
      A3 = (A+G)*E + H*(B/D)
      A1 = B*F - A
      A7 = (F+U)*A + H
      RETURN
   20 TYPE = 1
C TYPE=1 INDICATES THAT ALL FORMULAS ARE DIVIDED BY C
      E = A/C
      F = D/C
      G = U*E
      H = V*B
      A3 = A*E + (H/C+G)*B
      A1 = B - A*(D/C)
      A7 = A + G*D + H*F
      RETURN
      END

      SUBROUTINE NEXTK(TYPE)
C COMPUTES THE NEXT K POLYNOMIALS USING SCALARS
C COMPUTED IN CALCSC
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      REAL ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION TEMP, DABS
      INTEGER TYPE
      IF (TYPE.EQ.3) GO TO 40
      TEMP = A
      IF (TYPE.EQ.1) TEMP = B
      IF (DABS(A1).GT.DABS(TEMP)*ETA*10.) GO TO 20
C IF A1 IS NEARLY ZERO THEN USE A SPECIAL FORM OF THE
C RECURRENCE
      K(1) = 0.D0
      K(2) = -A7*QP(1)
      DO 10 I=3,N
        K(I) = A3*QK(I-2) - A7*QP(I-1)
   10 CONTINUE
      RETURN
C USE SCALED FORM OF THE RECURRENCE
   20 A7 = A7/A1
      A3 = A3/A1
      K(1) = QP(1)
      K(2) = QP(2) - A7*QP(1)
      DO 30 I=3,N
        K(I) = A3*QK(I-2) - A7*QP(I-1) + QP(I)
   30 CONTINUE
      RETURN
C USE UNSCALED FORM OF THE RECURRENCE IF TYPE IS 3
   40 K(1) = 0.D0
      K(2) = 0.D0
      DO 50 I=3,N
        K(I) = QK(I-2)
   50 CONTINUE
      RETURN
      END

      SUBROUTINE NEWEST(TYPE, UU, VV)
C COMPUTE NEW ESTIMATES OF THE QUADRATIC COEFFICIENTS
C USING THE SCALARS COMPUTED IN CALCSC.
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      REAL ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION A4, A5, B1, B2, C1, C2, C3,
     * C4, TEMP, UU, VV
      INTEGER TYPE
C USE FORMULAS APPROPRIATE TO SETTING OF TYPE.
      IF (TYPE.EQ.3) GO TO 30
      IF (TYPE.EQ.2) GO TO 10
      A4 = A + U*B + H*F
      A5 = C + (U+V*F)*D
      GO TO 20
   10 A4 = (A+G)*F + H
      A5 = (F+U)*C + V*D
C EVALUATE NEW QUADRATIC COEFFICIENTS.
   20 B1 = -K(N)/P(NN)
      B2 = -(K(N-1)+B1*P(N))/P(NN)
      C1 = V*B2*A1
      C2 = B1*A7
      C3 = B1*B1*A3
      C4 = C1 - C2 - C3
      TEMP = A5 + B1*A4 - C4
      IF (TEMP.EQ.0.D0) GO TO 30
      UU = U - (U*(C3+C2)+V*(B1*A1+B2*A7))/TEMP
      VV = V*(1.+C4/TEMP)
      RETURN
C IF TYPE=3 THE QUADRATIC IS ZEROED
   30 UU = 0.D0
      VV = 0.D0
      RETURN
      END

      SUBROUTINE QUADSD(NN, U, V, P, Q, A, B)
C DIVIDES P BY THE QUADRATIC  1,U,V  PLACING THE
C QUOTIENT IN Q AND THE REMAINDER IN A,B
      DOUBLE PRECISION P(NN), Q(NN), U, V, A, B, C
      INTEGER I
      B = P(1)
      Q(1) = B
      A = P(2) - U*B
      Q(2) = A
      DO 10 I=3,NN
        C = P(I) - U*A - V*B
        Q(I) = C
        B = A
        A = C
   10 CONTINUE
      RETURN
      END

      SUBROUTINE QUAD(A, B1, C, SR, SI, LR, LI)
C CALCULATE THE ZEROS OF THE QUADRATIC A*Z**2+B1*Z+C.
C THE QUADRATIC FORMULA, MODIFIED TO AVOID
C OVERFLOW, IS USED TO FIND THE LARGER ZERO IF THE
C ZEROS ARE REAL AND BOTH ZEROS ARE COMPLEX.
C THE SMALLER REAL ZERO IS FOUND DIRECTLY FROM THE
C PRODUCT OF THE ZEROS C/A.
      DOUBLE PRECISION A, B1, C, SR, SI, LR, LI, B,
     * D, E, DABS, DSQRT
      IF (A.NE.0.D0) GO TO 20
      SR = 0.D0
      IF (B1.NE.0.D0) SR = -C/B1
      LR = 0.D0
   10 SI = 0.D0
      LI = 0.D0
      RETURN
   20 IF (C.NE.0.D0) GO TO 30
      SR = 0.D0
      LR = -B1/A
      GO TO 10
C COMPUTE DISCRIMINANT AVOIDING OVERFLOW
   30 B = B1/2.D0
      IF (DABS(B).LT.DABS(C)) GO TO 40
      E = 1.D0 - (A/B)*(C/B)
      D = DSQRT(DABS(E))*DABS(B)
      GO TO 50
   40 E = A
      IF (C.LT.0.D0) E = -A
      E = B*(B/DABS(C)) - E
      D = DSQRT(DABS(E))*DSQRT(DABS(C))
   50 IF (E.LT.0.D0) GO TO 60
C REAL ZEROS
      IF (B.GE.0.D0) D = -D
      LR = (-B+D)/A
      SR = 0.D0
      IF (LR.NE.0.D0) SR = (C/LR)/A
      GO TO 10
C COMPLEX CONJUGATE ZEROS
   60 SR = -B/A
      LR = SR
      SI = DABS(D/A)
      LI = -SI
      RETURN
      END

      subroutine gamma_inc_values ( n_data, a, x, fx )

c*********************************************************************72
c
cc GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
c
c  Discussion:
c
c    The (normalized) incomplete Gamma function P(A,X) is defined as:
c
c      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
c
c    With this definition, for all A and X,
c
c      0 <= PN(A,X) <= 1
c
c    and
c
c      PN(A,INFINITY) = 1.0
c
c    In Mathematica, the function can be evaluated by:
c
c      1 - GammaRegularized[A,X]
c
c  Modified:
c
c    19 January 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision A, the parameter of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision a
      double precision a_vec ( n_max )
      double precision fx
      double precision fx_vec ( n_max )
      integer n_data
      double precision x
      double precision x_vec ( n_max )

      data a_vec /
     &   0.10D+00,
     &   0.10D+00,
     &   0.10D+00,
     &   0.50D+00,
     &   0.50D+00,
     &   0.50D+00,
     &   0.10D+01,
     &   0.10D+01,
     &   0.10D+01,
     &   0.11D+01,
     &   0.11D+01,
     &   0.11D+01,
     &   0.20D+01,
     &   0.20D+01,
     &   0.20D+01,
     &   0.60D+01,
     &   0.60D+01,
     &   0.11D+02,
     &   0.26D+02,
     &   0.41D+02 /
      data fx_vec /
     &   0.7382350532339351D+00,
     &   0.9083579897300343D+00,
     &   0.9886559833621947D+00,
     &   0.3014646416966613D+00,
     &   0.7793286380801532D+00,
     &   0.9918490284064973D+00,
     &   0.9516258196404043E-01,
     &   0.6321205588285577D+00,
     &   0.9932620530009145D+00,
     &   0.7205974576054322E-01,
     &   0.5891809618706485D+00,
     &   0.9915368159845525D+00,
     &   0.01018582711118352D+00,
     &   0.4421745996289254D+00,
     &   0.9927049442755639D+00,
     &   0.4202103819530612E-01,
     &   0.9796589705830716D+00,
     &   0.9226039842296429D+00,
     &   0.4470785799755852D+00,
     &   0.7444549220718699D+00 /
      data x_vec /
     &   0.30E-01,
     &   0.30D+00,
     &   0.15D+01,
     &   0.75D-01,
     &   0.75D+00,
     &   0.35D+01,
     &   0.10D+00,
     &   0.10D+01,
     &   0.50D+01,
     &   0.10D+00,
     &   0.10D+01,
     &   0.50D+01,
     &   0.15D+00,
     &   0.15D+01,
     &   0.70D+01,
     &   0.25D+01,
     &   0.12D+02,
     &   0.16D+02,
     &   0.25D+02,
     &   0.45D+02 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end

      SUBROUTINE amat_inverse(a,a_dim,a_inv)
      include 'inc.par'
      integer a_dim
      complex*16 a(a_dim,a_dim),a_inv(a_dim,a_dim)
      complex*16,dimension(:),allocatable:: WORK_proj
      integer,dimension(:),allocatable:: IPIV_proj


      allocate(WORK_proj(26*a_dim))
      allocate(IPIV_proj(a_dim))

      a_inv=a

      CALL ZSYTRF ('U',a_dim,a_inv,a_dim,
     &IPIV_proj, WORK_proj,-1, INFO_proj)
      LWORK_prj=nint(dble(WORK_proj(1)))
      deallocate(WORK_proj)
      allocate(WORK_proj(LWORK_prj))
      CALL ZSYTRF ('U',a_dim,a_inv,a_dim,
     &IPIV_proj, WORK_proj,LWORK_prj,INFO_proj)
      CALL ZSYTRI ('U',a_dim,a_inv,a_dim,
     &IPIV_proj, WORK_proj, INFO_proj )

      deallocate(IPIV_proj)
      deallocate(WORK_proj)

      return
      end

      SUBROUTINE amat_inverse_general(a,a_dim,a_inv)
      include 'inc.par'
      integer a_dim
      complex*16 a(a_dim,a_dim),a_inv(a_dim,a_dim)
      complex*16,dimension(:),allocatable:: WORK_proj
      integer,dimension(:),allocatable:: IPIV_proj

      allocate(WORK_proj(26*a_dim))
      allocate(IPIV_proj(a_dim))

      a_inv=a

      CALL ZGETRF (a_dim,a_dim,a_inv,a_dim,
     &IPIV_proj, INFO_proj)
      if(info_proj.ne.0)then
        write(*,*)'amat_inverse_general 1: INFO NON-ZERO',INFO_proj
        stop
      endif
      CALL ZGETRI (a_dim,a_inv,a_dim,
     &IPIV_proj, WORK_proj, -1,INFO_proj)
      LWORK=nint(dble(Work_proj(1)))
      deallocate(WORK_proj)
      allocate(WORK_proj(LWORK))
      if(info_proj.ne.0)then
        write(*,*)'amat_inverse_general 2: INFO NON-ZERO',INFO_proj
        stop
      endif
      CALL ZGETRI (a_dim,a_inv,a_dim,
     &IPIV_proj, WORK_proj,LWORK,INFO_proj)
      if(info_proj.ne.0)then
        write(*,*)'amat_inverse_general 3: INFO NON-ZERO',INFO_proj
        stop
      endif

      deallocate(IPIV_proj)
      deallocate(WORK_proj)

      return
      end

      SUBROUTINE amat_inverse_general_alternate(a,a_dim,a_inv)
      include 'inc.par'
      integer a_dim
      complex*16 a(a_dim,a_dim),a_inv(a_dim,a_dim)
      complex*16,dimension(:),allocatable:: WORK,det
      integer,dimension(:),allocatable:: IPIV

      allocate(IPIV(a_dim))
      allocate(work(a_dim))
      allocate(det(2))
      a_inv=a

      CALL zgefa(a_inv,a_dim,a_dim,ipiv,info)
      if(info.ne.0)then
        write(*,*)'amat_inverse_general_alternate: INFO NON-ZERO',INFO
        stop
      endif
      CALL zgedi(a_inv,a_dim,a_dim,ipiv,det,work,1)

      deallocate(IPIV)
      deallocate(WORK)
      deallocate(det)

      return
      end

      subroutine skew_diagonal(a_dim,aa,eignum,
     &eigvec)
      include 'inc.par'
      integer a_dim
      real*8 eignum(a_dim),aa(a_dim,a_dim),e(a_dim),trace
      complex*16 eigvec(a_dim,a_dim)
      real*8, dimension(:,:),allocatable:: a,z
      logical skew

      eignum=0.d0
      eigvec=(0.d0,0.d0)

      allocate(a(a_dim,a_dim))
      a=aa
      skew=.false.

      allocate(z(a_dim,a_dim))

      call trizd(a_dim, a_dim, a, e)

      trace=0.d0
      do i=1,a_dim
        trace=trace+a(i,i)
      enddo

      if(trace.eq.0.d0) skew=.true.

      do i=1,a_dim
        do j=1,a_dim
          if(i.ne.j)then
            if(a(i,j).ne. -a(j,i)) skew=.false.
          endif
        enddo
      enddo

      call imzd(a_dim, e, .true., skew, a_dim, z, ierr)
c e=e/dsqrt(2.d0)
      call tbakzd( a_dim, a_dim, a, a_dim, a_dim, z )

      eignum=e

      z_norm=0.d0

      do i=1,a_dim
        if(i+1.le.a_dim)then
          if(e(i).eq.-e(i+1))then
            do k=1,a_dim
              eigvec(k,i)=z(k,i)*(1.d0,0.d0)+z(k,i+1)*(0.d0,1.d0)
            enddo
            do k=1,a_dim
              eigvec(k,i+1)=z(k,i)*(1.d0,0.d0)+z(k,i+1)*(0.d0,-1.d0)
            enddo
          endif
        endif
      enddo
      do i=1,a_dim
        if(i.ge. max(2,a_dim) .and. i.le. max(1,a_dim-1))then
          if(e(i-1).ne.-e(i) .and. e(i).ne.-e(i+1))then
            do k=1,a_dim
              eigvec(k,i)=z(k,i)*(1,0)
            enddo
          endif
        elseif(i.eq.1)then
          if(e(i).ne.-e(i+1))then
            do k=1,a_dim
              eigvec(k,i)=z(k,i)*(1,0)
            enddo
          endif
        elseif(i.eq.a_dim)then
          if(e(i).ne.-e(i-1))then
            do k=1,a_dim
              eigvec(k,i)=z(k,i)*(1,0)
            enddo
          endif
        endif
      enddo

      do i=1,a_dim
        z_norm=z_norm+cdabs(eigvec(1,i)*conjg(eigvec(1,i)))
      enddo

      eigvec=eigvec/dsqrt(z_norm)

      deallocate(a)
      deallocate(z)
      return
      end

      subroutine coeff_reader(coefficient_strt,coefficient_end,
     &end_of_filename,coeff)
      include 'inc.par'
      integer coefficient_strt,coefficient_end
      complex*16 coeff(coefficient_strt:coefficient_end)
      real*8, dimension(:), allocatable:: Line_value
      character*50 end_of_filename
      character*3 en_eval1_fileoutput
      character*80 Line_readin
      integer, parameter:: Linemax=4

      do i=coefficient_strt,coefficient_end
        if(i.lt. 0) then
          write(en_eval1_fileoutput, 73)i
        else
          write(en_eval1_fileoutput, 78)i
        endif
        open(513+i, file=
     &  'Coeff_Raw_'//en_eval1_fileoutput//end_of_filename,
     &  status='old')
      enddo

      allocate(Line_value(Linemax))
      do j=coefficient_strt,coefficient_end
        do
          read (513+j, '(A)', end=999) line_readin
          read (line_readin, *, end=888) (Line_value(ii),ii=1,linemax)
 888      continue
          coeff(j)=line_value(2)*(1,0)+line_value(3)*(0,1)
        enddo
 999    continue
      enddo
      deallocate(Line_value)

  73  format(I3.2)
  78  format(I3.3)

      do i=coefficient_strt,coefficient_end
        close(513+i)
      enddo

      return
      end

      subroutine diving_sorter(eigval_e,nste,
     &wave_new_even,Lowest_bound_e,nu,alternate_dmat,dmat,
     &nkap,nm,dkb,tknot)
      include 'inc.par'
      real*8 eigval_e(nste),wave_new_even(nste,2*nm,-nkap:nkap),
     & alternate_dmat(2*nm,2*nm,-nkap:nkap),dmat(2*nm,2*nm),tknot(nu)
      integer One_s_state_location
      real*8, dimension(:),allocatable :: el,tempeigval_e
      real*8, dimension(:,:,:),allocatable:: wave_new_even_shift
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      logical dkb

      allocate(el(nste))
      allocate(wave_new_even_shift(nste,2*nm,-nkap:nkap))
      allocate(tempeigval_e(nste))
      tempeigval_e=eigval_e
      wave_new_even_shift=wave_new_even
      el=0.d0
      if(dkb)then
        do j=1,nste
          if(eigval_e(j).gt. -2.5d0 .and. eigval_e(j).lt.-1.d0)then
            l=1
              do while(tknot(l).lt. min(1.d0/(az1+az2),tknot(nu)))
                l=l+1
                m=1
                do while(tknot(m).lt. min(1.d0/(az1+az2),tknot(nu)))
                  m=m+1
c       do l=1,nm-10!nu/2!nm/2+1
c          do m=1,nm-10!nu/2!nm/2+1
                  el(j)=el(j)+(wave_new_even(j,l,-1)*
     &            wave_new_even(j,m,-1)*
     &            alternate_dmat(l,m,-1)+
     &            wave_new_even(j,l+nm,-1)*
     &            wave_new_even(j,m+nm,-1)*
     &            alternate_dmat(l+nm,m+nm,-1))
               enddo
             enddo
c         write(*,*)j,el(j)
          endif
        enddo
      else
        do j=1,nste
          if(eigval_e(j).gt. -2.5d0 .and. eigval_e(j).lt.-1.d0)then
            l=1
            do while(tknot(l).lt. min(1.d0/(az1+az2),tknot(nu)))
              l=l+1
              m=1
              do while(tknot(m).lt. min(1.d0/(az1+az2),tknot(nu)))
                m=m+1
c       do l=1,nm-10!nu/2!nm/2+1
c          do m=1,nm-10!nu/2!nm/2+1
                el(j)=el(j)+(wave_new_even(j,l,-1)*
     &          wave_new_even(j,m,-1)*
     &          dmat(l,m)+
     &          wave_new_even(j,l+nm,-1)*
     &          wave_new_even(j,m+nm,-1)*
     &          dmat(l+nm,m+nm))
              enddo
            enddo
c         write(*,*)j,el(j)
          endif
        enddo
      endif

      One_s_state_location=maxloc(el, dim=1, mask= el .lt. 1.d0)
      actual_1sEn=eigval_e(One_s_state_location)
      deallocate(el)
      do j=One_s_state_location,Lowest_bound_e-1
        eigval_e(j)=tempeigval_e(j+1)
        do l=1,2*nm
          do kk=-nkap,nkap
            if(kk.ne. 0) wave_new_even(j,l,kk)=
     &      wave_new_even_shift(j+1,l,kk)
          enddo
        enddo
      enddo
      eigval_e(Lowest_bound_e)=actual_1sEn
      do l=1,2*nm
        do kk=-nkap,nkap
          if(kk.ne. 0) wave_new_even(Lowest_bound_e,l,kk)=
     &      wave_new_even_shift(One_s_state_location,l,kk)
        enddo
      enddo
      deallocate(wave_new_even_shift)
      deallocate(tempeigval_e)
      write(*,*)'EXIT DIVING SORTER'
      return
      end

      subroutine spline_arranger_old(nu,rmin,rmax,t)
      include 'inc.par'
      common /Barycentres/ RadiusOne,RadiusTwo
      common /TargProj_prop/ Proj_mass,Targ_mass,Proj_vel,
     &b_ImpactParam
      common /dist/distance,Starting_Distance
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2
      real*8 t(nu)

      hh=dexp((dlog(rmax)-dlog(rmin))/dble(nu-2*ns))

      do i=1,ns
        t(i)=0.d0
      enddo
      do i=nu-ns+1,nu
        t(i)=rmax
      enddo

      if(z_nuc1.eq. z_nuc2) then
        ii=ns+1
        dist_min=dabs(t(ii)-distance)
        do i=ns+2,nu-ns
          if(dabs(t(i)-distance).lt.dist_min) then
            ii=i
            dist_min=dabs(t(i)-distance)
          endif
        enddo
        t(ii)=distance
      else
        t(ns+1)=rmin
        do i=ns+2,nu-ns
          t(i)=hh*t(i-1)
        enddo

        ii=ns+1
        dist_min=dabs(t(ii)-RadiusOne)
        do i=ns+2,nu-ns
          if(dabs(t(i)-RadiusOne).lt.dist_min) then
            ii=i
            dist_min=dabs(t(i)-RadiusOne)
          endif
        enddo
        t(ii)=RadiusOne

        ii=ns+1
        dist_min=dabs(t(ii)-RadiusTwo)
        do i=ns+2,nu-ns
          if(dabs(t(i)-RadiusTwo).lt.dist_min) then
            ii=i
            dist_min=dabs(t(i)-RadiusTwo)
          endif
        enddo
        t(ii)=RadiusTwo
      endif

      return
      end

      subroutine spline_arranger(nu,rmin,rmax,t)
      include 'inc.par'
      real*8 t(nu)

      hh=dexp((dlog(rmax)-dlog(rmin))/dble(nu-2*ns))

      do i=1,ns
        t(i)=0.d0
      enddo
      do i=nu-ns+1,nu
        t(i)=rmax
      enddo

      t(ns+1)=rmin

      do i=ns+2,nu-ns
        t(i)=hh*t(i-1)
      enddo

      return
      end

      subroutine sign_of_states_old(start_state,mat_dimension,
     &wave_new,nstates,nm,nkap,ii_xi,
     &xi_stepslower,rmin,rmax,alternate_dmat)

      include 'inc.par'
      integer start_state,mat_dimension,xi_stepslower,
     &kappa(nstates)
      real*8 signn(Start_state:Start_state+Mat_dimension-1),
     &gradient(Start_state:Start_state+Mat_dimension-1)
      real*8 wave_new(nstates,2*nm,-nkap:nkap)
      real*8 ro(ns),ro1(ns),ro2(ns),el(nstates,-nkap:nkap),
     &el_temp(-nkap:nkap),alternate_dmat(2*nm,2*nm,-nkap:nkap)
      real*8,dimension(:),allocatable::t
      character*80 Line_readin
      real*8, dimension(:), allocatable:: Line_value
      logical dkb!,mkdirs
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      common /Barycentres/ RadiusOne,RadiusTwo

      !mkdirs=makedirqq('WF')
      CALL SYSTEM("mkdir -p WF")

      do i_nst=1,nstates
        do kk=-nkap,nkap
          if(kk.ne.0)then
            b_norm=0.d0
            do i=1,nm
              do j=max(1,i-ns+1),min(i+ns-1,nm)
                b_norm=b_norm+
     &          wave_new(i_nst,i,kk)*wave_new(i_nst,j,kk)*
     &          alternate_dmat(i,j,kk)+
     &          wave_new(i_nst,i+nm,kk)*wave_new(i_nst,j+nm,kk)*
     &          alternate_dmat(i+nm,j+nm,kk)
              enddo
            enddo
            el(i_nst,kk)=dabs(b_norm)
          endif
        enddo
      enddo

      do i_nst=1,nstates
        do kk=-nkap,nkap
          if(kk.ne.0)then
            el_temp(kk)=el(i_nst,kk)
          endif
        enddo
c      write(*,*)maxloc(el_temp,dim=1)-nkap-1,eigval(i_nst)
        kappa(i_nst)=maxloc(el_temp,dim=1)-nkap-1
c      write(*,*)kappa(i_nst),
      enddo

      nu=nm+ns+2
      allocate(t(nu))

      if(ii_xi.gt.xi_stepslower)then
        open(22,file='WF/D_of_wf_at_small_r.dat',status='old')
        allocate(line_value(2))
        do i=start_state,start_state+Mat_dimension-1
          read(22,'(A)',end=300)Line_readin
          read(line_readin,*,end=400)(Line_value(ii),ii=1,2)
 400      continue
          signn(i)=Line_value(1)
          gradient(i)=Line_value(2)
c      write(*,*)i-start_state+1,signn(i)
        enddo
 300    continue
        deallocate(line_value)
        close(22)
      endif

      call spline_arranger(nu,rmin,rmax,t)
      sampl_point=1.01d0*rmin

      i2=1
      do while(t(i2).lt. sampl_point)
        i2=i2+1
      enddo
      if (t(i2).gt. sampl_point) then
        i2=i2-1
      endif

      call ddsplines(sampl_point,ro,ro1,ro2,i2,t,nu)
      Kapp_range=nint(dabs(amu)+0.5d0)
      open(22,file='WF/D_of_wf_at_small_r.dat',status='old')
      do i=start_state,start_state+Mat_dimension-1
        WFValue_node1=0.d0      !Value of G
        WFValue_node2=0.d0      !Value of DG
        WFValue_node3=0.d0      !Value of F
        WFValue_node4=0.d0      !Value of DF
c     MAYBE IT'S BETTER TO SAMPLE THE DERIVATIVE OF THE FUNCTION
c     NEAR 0.
      do kapp=-nkap,nkap
        if(abs(kapp).ge.Kapp_range)then
          do i_spl=max(2,i2-ns+1),min(i2,nm+1)
            if(dkb)then
              WFValue_node2=WFValue_node2+
     &        wave_new(i,i_spl-1,Kapp)*ro1(i_spl-i2+ns)+
     &        wave_new(i,i_spl+nm-1,Kapp)*0.5d0*(
     &        ro2(i_spl-i2+ns)-Kapp*ro1(i_spl-i2+ns)/sampl_point+
     &        Kapp*ro(i_spl-i2+ns)/sampl_point/sampl_point)
              WFValue_node1=WFValue_node1+
     &        wave_new(i,i_spl-1,Kapp)*ro(i_spl-i2+ns)+
     &        wave_new(i,i_spl+nm-1,Kapp)*0.5d0*(
     &        ro1(i_spl-i2+ns)-Kapp*ro(i_spl-i2+ns)/sampl_point)
              WFValue_node3=WFValue_node3+
     &        wave_new(i,i_spl+nm-1,Kapp)*ro(i_spl-i2+ns)+
     &        wave_new(i,i_spl-1,Kapp)*0.5d0*(
     &        ro1(i_spl-i2+ns)+Kapp*ro(i_spl-i2+ns)/sampl_point)
              WFValue_node4=WFValue_node4+
     &        wave_new(i,i_spl-1+nm,Kapp)*ro1(i_spl-i2+ns)+
     &        wave_new(i,i_spl-1,Kapp)*0.5d0*(
     &        ro2(i_spl-i2+ns)+Kapp*ro1(i_spl-i2+ns)/sampl_point-
     &        Kapp*ro(i_spl-i2+ns)/sampl_point/sampl_point)
            else
              WFValue_node2=WFValue_node2+
     &        wave_new(i,i_spl-1,Kapp)*ro1(i_spl-i2+ns)
              WFValue_node1=WFValue_node1+
     &        wave_new(i,i_spl-1,Kapp)*ro(i_spl-i2+ns)
              WFValue_node3=WFValue_node3+
     &        wave_new(i,i_spl-1+nm,Kapp)*ro(i_spl-i2+ns)
              WFValue_node4=WFValue_node4+
     &        wave_new(i,i_spl-1+nm,Kapp)*ro1(i_spl-i2+ns)
            endif
          enddo
        endif
      enddo
c         DofFonG=WFValue_node4/WFValue_node1-
c     &      WFValue_node3*WFValue_node2/WFValue_node1**2
c      write(*,*)i,DofFonG,WFValue_node1,WFValue_node2

      if(ii_xi.eq.xi_stepslower+2)then
        if(WFValue_node2*signn(i).lt. 0.d0)then
          do kap=-nkap,nkap
            if(kap.ne. 0) then
              do j=1,2*nm
                wave_new(i,j,kap)=-wave_new(i,j,kap)
              enddo
            endif
          enddo
          WFValue_node2=-WFValue_node2
        endif
      endif

      if(ii_xi.gt.xi_stepslower+3)then
        prediction=gradient(i)+signn(i)
        write(*,*)'PREDICTION',prediction,'ACTUAL',WFValue_node2,
     &  'PREVIOUS',signn(i),'VARIANCE',variance,'GRADIENT',gradient
        variance=dabs(WFValue_node2-prediction)
c            if(WFValue_node2*signn(i).lt. 0.d0)then
        if(variance.gt. dabs(10.d0*gradient(i)))then
          do kap=-nkap,nkap
            if(kap.ne. 0) then
              do j=1,2*nm
                wave_new(i,j,kap)=-wave_new(i,j,kap)
              enddo
            endif
          enddo
          WFValue_node2=-WFValue_node2
        endif
      endif

      if(ii_xi.eq.xi_stepslower)then
        delta_value=0.d0
      else
        delta_value=(WFValue_node2-signn(i))
      endif

      write(22,*)WFValue_node2,delta_value

      enddo
      close(22)
      deallocate(t)
c      pause

      return
      end

C     FULL WAVEFUNCTION SAMPLING
      subroutine sign_of_states(wave_new,nstates,nm,nkap,
     & ii_xi,xi_stepslower,rmin,rmax,typ)

      include 'inc.par'
      integer xi_stepslower
      real*8 plus_minus(nstates),plus_plus(nstates),
     &wfsummation_g(nstates),wfsummation_dg(nstates),
     &wfsummation_f(nstates),wfsummation_df(nstates),
     &sampl_point(nm,nstates),
     &wave_new(nstates,2*nm,-nkap:nkap),ro(ns),ro1(ns),ro2(ns),
     &abstand(nstates)
      real*8,dimension(:),allocatable::t,Line_value
      character*80 Line_readin
      character typ
      logical dkb,swapover
      common /common_dkb/ dkb!,mkdirs
      common /momentum_projection/ amu,amj_max
      common /Barycentres/ RadiusOne,RadiusTwo
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      !mkdirs=makedirqq('WF')
      CALL SYSTEM("mkdir -p WF")

      abstand=rmax*0.5d0
      lstep=1
      if(z_nuc1.eq.z_nuc2)lstep=2
      sampl_point=0.d0
      swapover=.false.

      if(ii_xi.gt.xi_stepslower)then
        open(22,file='WF/wf_at_small_r'//typ//'.dat',status='old')
        allocate(line_value(4))
        do i=1,nstates
          read(22,'(A)',end=300)Line_readin
          read(line_readin,'(4E17.10)',end=400)(Line_value(ii),ii=1,4)
 400      continue
          wfsummation_g(i)=Line_value(1)
          wfsummation_dg(i)=Line_value(2)
          wfsummation_f(i)=Line_value(3)
          wfsummation_df(i)=Line_value(4)
c      swapped_q(i)=Line_value(5)
c      plus_minus(i)=Line_value(2)
c      plus_plus(i)=Line_value(3)
c      write(*,*)i-start_state+1,signn(i)
        enddo
 300    continue
        deallocate(line_value)
        close(22)
c      write(*,*)wfsummation
c      pause
      endif

      nu=nm+ns+2
      allocate(t(nu))

      call spline_arranger(nu,rmin,rmax,t)

      do jstate=1,nstates
        ii=1
        do i=1,nm
          sampl_point(ii,jstate)=(abstand(jstate)/
     &    dabs(abstand(jstate)-1.d0))*
     &    dabs(abstand(jstate)**(dble((ii))/dble(nm))-1.d0)
          ii=ii+1
        enddo
      enddo

      Kapp_range=nint(dabs(amu)+0.5d0)
      open(22, file='WF/wf_at_small_r'//typ//'.dat')
      do i=1,nstates

        sampl_sum_g=0.d0
        sampl_sum_dg=0.d0
        sampl_sum_f=0.d0
        sampl_sum_df=0.d0

        do ji=1,ii-1
          if(ji.ne.0)then
            Sample_point=sampl_point(abs(ji),i)
            i2=1
            do while(t(i2).lt. Sample_point)
              i2=i2+1
            enddo
            if (t(i2).gt. Sample_point) then
              i2=i2-1
            endif
c         write(*,*)Sample_point,sampl_point(abs(ji),i),i2
c         pause
            call ddsplines(Sample_point,ro,ro1,ro2,i2,t,nu)

            wfgtotal=0.d0
            wfdgtotal=0.d0
            wfftotal=0.d0
            wfdftotal=0.d0

            do kk=-nkap,nkap
              if(abs(kk).ge. abs(Kapp_range))then

                WFVal_G1=0.d0
                WFVal_G=0.d0
                WFVal_DG=0.d0
                WFVal_F=0.d0
                WFVal_DF=0.d0

                call DiracAngularL(-1.d0*kk,al1l)
                aj1j=dabs(1.d0*kk)-0.5d0
                ffakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
                ffakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
                call DiracAngularL(1.d0*kk,al1l)
                gfakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
                gfakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
                do i_spl=max(2,i2-ns+1),min(i2,nm+1)
                  if(dkb)then

                    WFVal_G=WFVal_G+
     &              wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1+nm,kk)*0.5d0*
     &              (ro1(i_spl-i2+ns)-kk*ro(i_spl-i2+ns)/Sample_point)
                    WFVal_DG=WFVal_DG+
     &              wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1+nm,kk)*(ro2(i_spl-i2+ns)-
     &              kk*ro1(i_spl-i2+ns)/Sample_point+
     &              kk*ro(i_spl-i2+ns)/Sample_point**2)
     &              /2.d0
                    WFVal_F=WFVal_F+
     &              wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1,kk)*0.5d0*
     &              (ro1(i_spl-i2+ns)+kk*ro(i_spl-i2+ns)/Sample_point)
                    WFVal_DF=WFVal_DF+
     &              wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1,kk)*(ro2(i_spl-i2+ns)+
     &              kk*ro1(i_spl-i2+ns)/Sample_point-
     &              kk*ro(i_spl-i2+ns)/Sample_point**2)
     &              /2.d0

                  else

                    WFVal_G=WFVal_G+
     &              wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)
                    WFVal_DG=WFVal_DG+
     &              wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)
                    WFVal_F=WFVal_F+
     &              wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)
                    WFVal_DF=WFVal_DF+
     &              wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)

                  endif
                enddo
                sampl_sum_g=sampl_sum_g+WFVal_G*gfakt1!*
c     &      (Sample_point**2+1.d0)
                sampl_sum_dg=sampl_sum_dg+WFVal_DG*gfakt1!*
c     &      (Sample_point**2+1.d0)
                wfgtotal=wfgtotal+WFVal_G*gfakt1
                wfdgtotal=wfdgtotal+WFVal_DG*gfakt1
                sampl_sum_f=sampl_sum_f+WFVal_F*ffakt1!*
c     &      (Sample_point**2+1.d0)
                sampl_sum_df=sampl_sum_df+WFVal_DF*ffakt1!*
c     &      (Sample_point**2+1.d0)
                wfftotal=wfftotal+WFVal_F*ffakt1
                wfdftotal=wfdftotal+WFVal_DF*ffakt1
              endif
            enddo
          endif
        enddo

        if(ii_xi.gt.xi_stepslower)then
          if(((dabs(sampl_sum_g-wfsummation_g(i)) .gt.
     &    dabs(sampl_sum_g+wfsummation_g(i))).and.
     &    (dabs(sampl_sum_dg-wfsummation_dg(i)) .gt.
     &    dabs(sampl_sum_dg+wfsummation_dg(i)))
     &    ).or.
     &    ((dabs(sampl_sum_f-wfsummation_f(i)) .gt.
     &    dabs(sampl_sum_f+wfsummation_f(i))).and.
     &    (dabs(sampl_sum_df-wfsummation_df(i)) .gt.
     &    dabs(sampl_sum_df+wfsummation_df(i))))
     &    )then
            swapover=.true.
            do kp=-nkap,nkap
              if(kp.ne.0)then
                do ji1=1,2*nm
                  wave_new(i,ji1,kp)=-wave_new(i,ji1,kp)
                enddo
              endif
            enddo
          else
            swapover=.false.
          endif
        endif

        if(swapover)then
          sampl_sum_g=-sampl_sum_g
          sampl_sum_dg=-sampl_sum_dg
          sampl_sum_f=-sampl_sum_f
          sampl_sum_df=-sampl_sum_df
        endif

        if(ii_xi.eq.xi_stepslower)then
          plus_minus=0.d0
          plus_plus=0.d0
          swapover=.false.
        else
          plus_minus(i)=dabs(sampl_sum_g-wfsummation_g(i))
          plus_plus(i)=dabs(sampl_sum_g+wfsummation_g(i))
        endif
        if(swapover)then
          swpvr=1.d0
        else
          swpvr=0.d0
        endif
        write(22,27)sampl_sum_g,sampl_sum_dg,sampl_sum_f,sampl_sum_df
 27     format(4E17.10)
      enddo
      close(22)
      deallocate(t)
      return
      end

      subroutine sign_of_states_monopole(number_states,
     &wave,nstates,nm,nkap,ii_xi,xi_stepslower,rmin,rmax,typ)

      include 'inc.par'
      integer xi_stepslower,number_states(2,-nkap:nkap)
      real*8 wfsummation_g(nstates),wfsummation_dg(nstates),
     &wfsummation_f(nstates),wfsummation_df(nstates),
     &sampl_point(nm,nstates),wave(2*nm,2*nm,-nkap:nkap),ro(ns),
     &ro1(ns),ro2(ns),abstand(nstates)
      real*8,dimension(:),allocatable::t,Line_value
      character*80 Line_readin
      character*3 env
      character typ
      logical dkb,swapover!,mkdirs
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      common /Barycentres/ RadiusOne,RadiusTwo
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      !mkdirs=makedirqq('WF')
      CALL SYSTEM("mkdir -p WF")

      sampl_point=0.d0
      swapover=.false.
      abstand=1.d0/(az1+az2)

      if(ii_xi.gt.xi_stepslower)then
        open(22,file='WF/wf_at_small_r_'//typ//'_monopole.dat',
     &  status='old')
        allocate(line_value(4))
        do i=1,nstates
          read(22,'(A)',end=300)Line_readin
          read(line_readin,'(4E17.10)',end=400)(Line_value(ii),ii=1,4)
 400      continue
          wfsummation_g(i)=Line_value(1)
          wfsummation_dg(i)=Line_value(2)
          wfsummation_f(i)=Line_value(3)
          wfsummation_df(i)=Line_value(4)
        enddo
 300    continue
        deallocate(line_value)
        close(22)
      endif

      nu=nm+ns+2
      allocate(t(nu))

      call spline_arranger(nu,rmin,rmax,t)

      do jstate=1,nstates
        ii=1
        do i=1,nm
          sampl_point(ii,jstate)=(abstand(jstate)/
     &    dabs(abstand(jstate)-1.d0))*
     &    dabs(abstand(jstate)**(dble((ii))/dble(nm))-1.d0)
          ii=ii+1
        enddo
      enddo
      Kapp_range=nint(dabs(amu)+0.5d0)
      open(22, file='WF/wf_at_small_r_'//typ//'_monopole.dat')
      i9=1
      do k=1,2*nkap
        kk=(k+(1-(-1)**k)/2)*(-1)**k/2
        do i=number_states(1,kk),number_states(2,kk)
          write(env,24)i9
 24       format(I3.3)

          sampl_sum_g=0.d0
          sampl_sum_dg=0.d0
          sampl_sum_f=0.d0
          sampl_sum_df=0.d0

          do ji=1,ii-1
            Sample_point=sampl_point(abs(ji),i9)
            i2=1
            do while(t(i2).lt. Sample_point)
              i2=i2+1
            enddo
            if (t(i2).gt. Sample_point) then
              i2=i2-1
              call ddsplines(Sample_point,ro,ro1,ro2,i2,t,nu)

              WFVal_G=0.d0
              WFVal_DG=0.d0
              WFVal_F=0.d0
              WFVal_DF=0.d0

              do i_spl=max(2,i2-ns+1),min(i2,nm+1)
                if(dkb)then

                  WFVal_G=WFVal_G+
     &            wave(i_spl-1,i,kk)*ro(i_spl-i2+ns)+
     &            wave(i_spl-1+nm,i,kk)*0.5d0*
     &            (ro1(i_spl-i2+ns)-kk*ro(i_spl-i2+ns)/Sample_point)
                  WFVal_DG=WFVal_DG
     &            +wave(i_spl-1,i,kk)*ro1(i_spl-i2+ns)+
     &            wave(i_spl-1+nm,i,kk)*(ro2(i_spl-i2+ns)
     &            -kk*ro1(i_spl-i2+ns)/Sample_point+
     &            kk*ro(i_spl-i2+ns)/Sample_point**2)
     &            /2.d0
                  WFVal_F=WFVal_F
     &            +wave(i_spl-1+nm,i,kk)*ro(i_spl-i2+ns)+
     &            wave(i_spl-1,i,kk)*0.5d0*
     &            (ro1(i_spl-i2+ns)+kk*ro(i_spl-i2+ns)/Sample_point)
                  WFVal_DF=WFVal_DF
     &            +wave(i_spl-1+nm,i,kk)*ro1(i_spl-i2+ns)+
     &            wave(i_spl-1,i,kk)*(ro2(i_spl-i2+ns)
     &            +kk*ro1(i_spl-i2+ns)/Sample_point
     &            -kk*ro(i_spl-i2+ns)/Sample_point**2)
     &            /2.d0

                else

                  WFVal_G=WFVal_G+
     &            wave(i_spl-1,i,kk)*ro(i_spl-i2+ns)
                  WFVal_DG=WFVal_DG+
     &            wave(i_spl-1,i,kk)*ro1(i_spl-i2+ns)
                  WFVal_F=WFVal_F+
     &            wave(i_spl-1+nm,i,kk)*ro(i_spl-i2+ns)
                  WFVal_DF=WFVal_DF+
     &            wave(i_spl-1+nm,i,kk)*ro1(i_spl-i2+ns)

                endif
              enddo
              sampl_sum_g=sampl_sum_g+WFVal_G
              sampl_sum_dg=sampl_sum_dg+WFVal_DG
              sampl_sum_f=sampl_sum_f+WFVal_F
              sampl_sum_df=sampl_sum_df+WFVal_DF
            endif
          enddo

          if(ii_xi.gt.xi_stepslower)then
            if(((dabs(sampl_sum_g-wfsummation_g(i9)) .gt.
     &      dabs(sampl_sum_g+wfsummation_g(i9))).and.
     &      (dabs(sampl_sum_dg-wfsummation_dg(i9)) .gt.
     &      dabs(sampl_sum_dg+wfsummation_dg(i9)))
     &      ).or.
     &      ((dabs(sampl_sum_f-wfsummation_f(i9)) .gt.
     &      dabs(sampl_sum_f+wfsummation_f(i9))).and.
     &      (dabs(sampl_sum_df-wfsummation_df(i9)) .gt.
     &      dabs(sampl_sum_df+wfsummation_df(i9)))
     &      )
     &      )then
              swapover=.true.
              do ji1=1,2*nm
                wave(ji1,i,kk)=-wave(ji1,i,kk)
              enddo
            else
              swapover=.false.
            endif
          endif

          if(swapover)then
            sampl_sum_g=-sampl_sum_g
            sampl_sum_dg=-sampl_sum_dg
            sampl_sum_f=-sampl_sum_f
            sampl_sum_df=-sampl_sum_df
          endif
          write(22,27)sampl_sum_g,sampl_sum_dg,sampl_sum_f,sampl_sum_df
 27       format(4E17.10)
          i9=i9+1
        enddo
      enddo
      close(22)
      deallocate(t)
      return
      end

c     SINGLE POINT SAMPLING
      subroutine sign_of_states_sps(wave_new,nstates,nm,nkap,rmin,rmax,
     &alternate_dmat)
      include 'inc.par'
      integer kappa(nstates)
      real*8 alternate_dmat(2*nm,2*nm,-nkap:nkap),
     &el(nstates,-nkap:nkap),el_temp(-nkap:nkap)
      real*8 wave_new(nstates,2*nm,-nkap:nkap)
      real*8 ro(ns),ro1(ns),ro2(ns)
      real*8,dimension(:),allocatable::t
      logical dkb
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      common /Barycentres/ RadiusOne,RadiusTwo
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      do i_nst=1,nstates
        do kk=-nkap,nkap
          if(kk.ne.0)then
            b_norm=0.d0
            do i=1,nm
              do j=max(1,i-ns+1),min(i+ns-1,nm)
                b_norm=b_norm+
     &          wave_new(i_nst,i,kk)*wave_new(i_nst,j,kk)*
     &          alternate_dmat(i,j,kk)+
     &          wave_new(i_nst,i+nm,kk)*wave_new(i_nst,j+nm,kk)*
     &          alternate_dmat(i+nm,j+nm,kk)
              enddo
            enddo
            el(i_nst,kk)=dabs(b_norm)
          endif
        enddo
      enddo

      do i_nst=1,nstates
        do kk=-nkap,nkap
          if(kk.ne.0)then
            el_temp(kk)=el(i_nst,kk)
          endif
        enddo
c      write(*,*)maxloc(el_temp,dim=1)-nkap-1,eigval(i_nst)
        kappa(i_nst)=maxloc(el_temp,dim=1)-nkap-1
c      write(*,*)kappa(i_nst),
      enddo

      nu=nm+ns+2
      allocate(t(nu))

      call spline_arranger(nu,rmin,rmax,t)
      sampl_point=RadiusOne+1.d0/az1

      i2=1
      do while(t(i2).lt. sampl_point)
        i2=i2+1
      enddo
      if (t(i2).gt. sampl_point) then
        i2=i2-1
      endif

      call ddsplines(sampl_point,ro,ro1,ro2,i2,t,nu)
      Kapp_range=nint(dabs(amu)+0.5d0)
c      kk=-kapp_range
      do i=1,nstates
        kk1=kappa(i)

        do kk=-nkap,nkap
          if(abs(kk).ge. abs(Kapp_range))then

            WFValue_G=0.d0
            WFValue_DG=0.d0
            WFValue_F=0.d0
            WFValue_DF=0.d0

            do i_spl=max(2,i2-ns+1),min(i2,nm+1)
              if(dkb)then
                WFValue_G=WFValue_G+
     &          wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)+
     &          wave_new(i,i_spl-1+nm,kk)*0.5d0*
     &          (ro1(i_spl-i2+ns)-kk*ro(i_spl-i2+ns)/sampl_point)
                WFValue_DG=WFValue_DG+
     &          wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)+
     &          wave_new(i,i_spl-1+nm,kk)*(ro2(i_spl-i2+ns)-
     &          kk*ro1(i_spl-i2+ns)/sampl_point+
     &          kk*ro(i_spl-i2+ns)/sampl_point**2)
     &          /2.d0
                WFValue_F=WFValue_F+
     &          wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)+
     &          wave_new(i,i_spl-1,kk)*0.5d0*
     &          (ro1(i_spl-i2+ns)+kk*ro(i_spl-i2+ns)/sampl_point)
                WFValue_DF=WFValue_DF+
     &          wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)+
     &          wave_new(i,i_spl-1,kk)*(ro2(i_spl-i2+ns)+
     &          kk*ro1(i_spl-i2+ns)/sampl_point-
     &          kk*ro(i_spl-i2+ns)/sampl_point**2)
     &          /2.d0
              else
                WFValue_G=WFValue_G+
     &          wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)
                WFValue_DG=WFValue_DG+
     &          wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)
                WFValue_F=WFValue_F+
     &          wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)
                WFValue_DF=WFValue_DF+
     &          wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)
              endif
            enddo

c       if(WFValue_G*kk.gt. 0.d0 )then!.or.
cc     &      WFValue_node1*kk.lt. 0.d0)then
c          do j=1,nm!2*nm
c               wave_new(i,j,kk)=-wave_new(i,j,kk)
c               enddo
c       endif
            if(WFValue_F.gt. 0.d0 )then
              do j=1,2*nm!nm+1,2*nm
                wave_new(i,j,kk)=-wave_new(i,j,kk)
              enddo
            endif

          endif
        enddo
      enddo

      deallocate(t)

      return
      end

      subroutine Rotating_INaxis(nstates,nm,nkap,dthetadt,
     &nsts,bb_mjj,dmat,alt_dmat,wave_new,d_mjMax, d_number_states_mj,
     &n_jstates)
      include 'inc.par'
      complex*16 bb_mjj(nsts,nsts)
      real*8 dmat(2*nm,2*nm),alt_dmat(2*nm,2*nm,-nkap:nkap),
     &wave_new(nstates,2*nm,-nkap:nkap)
      real*8, dimension(:,:,:),allocatable:: norm_mat
      real*8, dimension(:,:,:,:),allocatable:: shifted_norm_plus,
     & shifted_norm_minus
      common /common_dkb/ dkb
      logical dkb
      real*8 d_number_states_mj(2*n_jstates*nstates)

      mj_max=nint(d_mjMax+0.5d0)

      allocate(norm_mat(nstates,nstates,-nkap:nkap))
      norm_mat=0.d0
!$OMP PARALLEL
!$OMP&PRIVATE(i,j,ki,stor_norm,ii,jj)
!$OMP&SHARED(nstates,nkap,norm_mat,dkb,nm,wave_new,alt_dmat,dmat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          do ki=-nkap,nkap
            if(ki.ne.0)then
              stor_norm=0.d0
              if(dkb)then
                do ii=1,nm
                  do jj=max(1,i-ns+1),min(i+ns-1,nm)
                    stor_norm=stor_norm+wave_new(i,ii,ki)*
     &              wave_new(j,jj,ki)*alt_dmat(ii,jj,ki)+
     &              wave_new(i,ii+nm,ki)*wave_new(j,jj+nm,ki)*
     &              alt_dmat(ii+nm,jj+nm,ki)
                  enddo
                enddo
              else
                do ii=1,nm
                  do jj=max(1,i-ns+1),min(i+ns-1,nm)
                    stor_norm=stor_norm+wave_new(i,ii,ki)*
     &              wave_new(j,jj,ki)*dmat(ii,jj)+
     &              wave_new(i,ii+nm,ki)*wave_new(j,jj+nm,ki)*
     &              dmat(ii+nm,jj+nm)
                  enddo
                enddo
              endif
            endif
            norm_mat(i,j,ki)=stor_norm
            if(i.ne.j)norm_mat(j,i,ki)=norm_mat(i,j,ki)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      allocate(shifted_norm_plus(nstates,nstates,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))
      allocate(shifted_norm_minus(nstates,nstates,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))

      shifted_norm_plus=0.d0
      shifted_norm_minus=0.d0
      bb_mjj=0.d0
!$OMP PARALLEL
!$OMP&PRIVATE(i,j,mj,a_mj,mj2,a_mj2,k,a_jk)
!$OMP&SHARED(nstates,mj_max,nkap,shifted_norm_minus,norm_mat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          do mj=-mj_max+1,mj_max-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            do k=-nkap,nkap
              a_jk=1.d0*abs(k)-0.5d0
              if(dabs(a_mj2).le.a_jk .and. dabs(a_mj).le.a_jk
     &        .and. k.ne.0)then
                shifted_norm_minus(i,j,mj,mj2)=
     &          shifted_norm_minus(i,j,mj,mj2)+
     &          dsqrt((a_jk+a_mj+1.d0)*(a_jk-a_mj))*norm_mat(i,j,k)
              endif
            enddo
            shifted_norm_minus(j,i,mj,mj2)=
     &      shifted_norm_minus(i,j,mj,mj2)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP&PRIVATE(i,j,mj,a_mj,mj2,a_mj2,k,a_jk)
!$OMP&SHARED(nstates,mj_max,nkap,shifted_norm_plus,norm_mat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          do mj=-mj_max,mj_max-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            do k=-nkap,nkap
              a_jk=1.d0*abs(k)-0.5d0
              if(dabs(a_mj2).le.a_jk .and. dabs(a_mj).le.a_jk
     &        .and. k.ne.0)then
                shifted_norm_plus(i,j,mj,mj2)=
     &          shifted_norm_plus(i,j,mj,mj2)+
     &          dsqrt((a_jk-a_mj+1.d0)*(a_jk+a_mj))*norm_mat(i,j,k)
              endif
            enddo
            shifted_norm_plus(j,i,mj,mj2)=
     &      shifted_norm_plus(i,j,mj,mj2)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      deallocate(norm_mat)

      j_start=1
      i_count=0
      do mj1=-mj_max,mj_max-1
        d_number_states_mj((i_count*nstates+1):(i_count+1)*nstates) =
     &  -d_mjMax+dble(i_count)
        i_count = i_count + 1
        i_col=1
        do mj2=-mj_max,mj_max-1
          do i=1,nstates
            j_row=j_start
            do j=1,nstates
              bb_mjj(i_col,j_row)=(0.5d0,0.d0)*
     &        (shifted_norm_minus(i,j,mj1,mj2)-
     &        shifted_norm_plus(i,j,mj1,mj2))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
        enddo
        j_start=j_row
      enddo

      bb_mjj=bb_mjj*(0,-1)*dthetadt
      deallocate(shifted_norm_minus)
      deallocate(shifted_norm_plus)

      return
      end

      subroutine Rotating_INaxis_even(nstates,nm,nkap,dthetadt,
     &nsts,bb_mjj,dmat,alt_dmat,wave_new,d_mjMax,d_number_states_mj,
     &n_jstates)
      include 'inc.par'
      complex*16 bb_mjj(nsts,nsts)
      real*8 dmat(2*nm,2*nm),alt_dmat(2*nm,2*nm,-nkap:nkap),
     &wave_new(nstates,2*nm,-nkap:nkap)
      real*8, dimension(:,:,:),allocatable:: norm_mat
      real*8, dimension(:,:,:,:),allocatable:: shifted_norm_plus,
     &shifted_norm_minus
      common /common_dkb/ dkb
      logical dkb
      real*8 d_number_states_mj(2*n_jstates*nstates)

      mj_max=nint(d_mjMax+0.5d0)

      allocate(norm_mat(nstates,nstates,-nkap:nkap))
      norm_mat=0.d0
!$OMP PARALLEL
!$OMP&PRIVATE(i,j,ki,stor_norm,ii,jj,k1fact)
!$OMP&SHARED(nstates,nkap,norm_mat,dkb,nm,wave_new,alt_dmat,dmat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          k1fact=(-1)**nkap
          do k1=nkap,1,-1
            ki=k1*k1fact
            k1fact=-k1fact
            stor_norm=0.d0
            if(dkb)then
              do ii=1,nm
                do jj=max(1,i-ns+1),min(i+ns-1,nm)
                  stor_norm=stor_norm+wave_new(i,ii,ki)*
     &            wave_new(j,jj,ki)*alt_dmat(ii,jj,ki)+
     &            wave_new(i,ii+nm,ki)*wave_new(j,jj+nm,ki)*
     &            alt_dmat(ii+nm,jj+nm,ki)
                enddo
              enddo
            else
              do ii=1,nm
                do jj=max(1,i-ns+1),min(i+ns-1,nm)
                  stor_norm=stor_norm+wave_new(i,ii,ki)*
     &            wave_new(j,jj,ki)*dmat(ii,jj)+
     &            wave_new(i,ii+nm,ki)*wave_new(j,jj+nm,ki)*
     &            dmat(ii+nm,jj+nm)
                enddo
              enddo
            endif
            norm_mat(i,j,ki)=stor_norm
            if(i.ne.j)norm_mat(j,i,ki)=norm_mat(i,j,ki)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      allocate(shifted_norm_plus(nstates,nstates,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))
      allocate(shifted_norm_minus(nstates,nstates,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))

      shifted_norm_plus=0.d0
      shifted_norm_minus=0.d0
      bb_mjj=0.d0
!$OMP PARALLEL
!$OMP&PRIVATE(i,j,mj,a_mj,mj2,a_mj2,k,a_jk,kfact)
!$OMP&SHARED(nstates,mj_max,nkap,shifted_norm_minus,norm_mat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          do mj=-mj_max+1,mj_max-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            kfact=(-1)**nkap
            do kk=nkap,1,-1
              a_jk=1.d0*kk-0.5d0
              k=kk*kfact
              kfact=-kfact
              if(dabs(a_mj2).le.a_jk .and. dabs(a_mj).le.a_jk)then
                shifted_norm_minus(i,j,mj,mj2)=
     &          shifted_norm_minus(i,j,mj,mj2)+
     &          dsqrt((a_jk+a_mj+1.d0)*(a_jk-a_mj))*norm_mat(i,j,k)
              endif
            enddo
            shifted_norm_minus(j,i,mj,mj2)=
     &      shifted_norm_minus(i,j,mj,mj2)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP&PRIVATE(i,j,mj,a_mj,mj2,a_mj2,k,a_jk,kfact)
!$OMP&SHARED(nstates,mj_max,nkap,shifted_norm_plus,norm_mat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          do mj=-mj_max,mj_max-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            kfact=(-1)**nkap
            do kk=nkap,1,-1
              a_jk=1.d0*kk-0.5d0
              k=kk*kfact
              kfact=-kfact
              if(dabs(a_mj2).le.a_jk .and. dabs(a_mj).le.a_jk)then
                shifted_norm_plus(i,j,mj,mj2)=
     &          shifted_norm_plus(i,j,mj,mj2)+
     &          dsqrt((a_jk-a_mj+1.d0)*(a_jk+a_mj))*norm_mat(i,j,k)
              endif
            enddo
            shifted_norm_plus(j,i,mj,mj2)=
     &      shifted_norm_plus(i,j,mj,mj2)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      deallocate(norm_mat)


      i_col=1
      i_count=0
      do mj1=-mj_max,mj_max-1
        d_number_states_mj((i_count*nstates+1):(i_count+1)*nstates) =
     &  -d_mjMax+dble(i_count)
        i_count = i_count + 1
        do i=1,nstates
          j_row=1
          do mj2=-mj_max,mj_max-1
            do j=1,nstates
              bb_mjj(i_col,j_row)=(0.5d0,0.d0)*
     &        (shifted_norm_minus(i,j,mj1,mj2)-
     &        shifted_norm_plus(i,j,mj1,mj2))
              j_row=j_row+1
            enddo
          enddo
          i_col=i_col+1
        enddo
      enddo
!      j_start=1
!      do mj1=-mj_max,mj_max-1
!        i_col=1
!        do mj2=-mj_max,mj_max-1
!          do i=1,nstates
!            j_row=j_start
!            do j=1,nstates
!              bb_mjj(i_col,j_row)=(0.5d0,0.d0)*
!     &        (shifted_norm_minus(i,j,mj1,mj2)-
!     &        shifted_norm_plus(i,j,mj1,mj2))
!              j_row=j_row+1
!            enddo
!            i_col=i_col+1
!          enddo
!        enddo
!        j_start=j_row
!      enddo

      bb_mjj=bb_mjj*(0,-1)*dthetadt

      deallocate(shifted_norm_minus)
      deallocate(shifted_norm_plus)

      return
      end

      subroutine Rotating_INaxis_odd(nstates,nm,nkap,dthetadt,
     &nsts,bb_mjj,dmat,alt_dmat,wave_new,d_mjMax,d_number_states_mj,
     &n_jstates)
      include 'inc.par'
      complex*16 bb_mjj(nsts,nsts)
      real*8 dmat(2*nm,2*nm),alt_dmat(2*nm,2*nm,-nkap:nkap),
     &wave_new(nstates,2*nm,-nkap:nkap)
      real*8, dimension(:,:,:),allocatable:: norm_mat
      real*8, dimension(:,:,:,:),allocatable:: shifted_norm_plus,
     &shifted_norm_minus
      common /common_dkb/ dkb
      logical dkb
      real*8 d_number_states_mj(2*n_jstates*nstates)

      mj_max=nint(d_mjMax+0.5d0)

      allocate(norm_mat(nstates,nstates,-nkap:nkap))
      norm_mat=0.d0
!$OMP PARALLEL
!$OMP&PRIVATE(i,j,ki,stor_norm,ii,jj,k1fact)
!$OMP&SHARED(nstates,nkap,norm_mat,dkb,nm,wave_new,alt_dmat,dmat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          k1fact=(-1)**(nkap+1)
          do k1=nkap,1,-1
            ki=k1*k1fact
            k1fact=-k1fact
            stor_norm=0.d0
            if(dkb)then
              do ii=1,nm
                do jj=max(1,i-ns+1),min(i+ns-1,nm)
                  stor_norm=stor_norm+wave_new(i,ii,ki)*
     &            wave_new(j,jj,ki)*alt_dmat(ii,jj,ki)+
     &            wave_new(i,ii+nm,ki)*wave_new(j,jj+nm,ki)*
     &            alt_dmat(ii+nm,jj+nm,ki)
                enddo
              enddo
            else
              do ii=1,nm
                do jj=max(1,i-ns+1),min(i+ns-1,nm)
                  stor_norm=stor_norm+wave_new(i,ii,ki)*
     &            wave_new(j,jj,ki)*dmat(ii,jj)+
     &            wave_new(i,ii+nm,ki)*wave_new(j,jj+nm,ki)*
     &            dmat(ii+nm,jj+nm)
                enddo
              enddo
            endif
            norm_mat(i,j,ki)=stor_norm
            if(i.ne.j)norm_mat(j,i,ki)=norm_mat(i,j,ki)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      allocate(shifted_norm_plus(nstates,nstates,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))
      allocate(shifted_norm_minus(nstates,nstates,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))

      shifted_norm_plus=0.d0
      shifted_norm_minus=0.d0
      bb_mjj=0.d0
!$OMP PARALLEL
!$OMP&PRIVATE(i,j,mj,a_mj,mj2,a_mj2,k,a_jk,kfact)
!$OMP&SHARED(nstates,mj_max,nkap,shifted_norm_minus,norm_mat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          do mj=-mj_max+1,mj_max-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            kfact=(-1)**(nkap+1)
            do kk=nkap,1,-1
              a_jk=1.d0*kk-0.5d0
              k=kk*kfact
              kfact=-kfact
              if(dabs(a_mj2).le.a_jk .and. dabs(a_mj).le.a_jk)then
                shifted_norm_minus(i,j,mj,mj2)=
     &          shifted_norm_minus(i,j,mj,mj2)+
     &          dsqrt((a_jk+a_mj+1.d0)*(a_jk-a_mj))*norm_mat(i,j,k)
              endif
            enddo
            shifted_norm_minus(j,i,mj,mj2)=
     &      shifted_norm_minus(i,j,mj,mj2)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP&PRIVATE(i,j,mj,a_mj,mj2,a_mj2,k,a_jk,kfact)
!$OMP&SHARED(nstates,mj_max,nkap,shifted_norm_plus,norm_mat)
!$OMP&default(none)
!$OMP DO
      do i=1,nstates
        do j=i,nstates
          do mj=-mj_max,mj_max-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            kfact=(-1)**(nkap+1)
            do kk=nkap,1,-1
              a_jk=1.d0*kk-0.5d0
              k=kk*kfact
              kfact=-kfact
              if(dabs(a_mj2).le.a_jk .and. dabs(a_mj).le.a_jk)then
                shifted_norm_plus(i,j,mj,mj2)=
     &          shifted_norm_plus(i,j,mj,mj2)+
     &          dsqrt((a_jk-a_mj+1.d0)*(a_jk+a_mj))*norm_mat(i,j,k)
              endif
            enddo
            shifted_norm_plus(j,i,mj,mj2)=
     &      shifted_norm_plus(i,j,mj,mj2)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL


      deallocate(norm_mat)

      j_start=1
      i_count=0
      do mj1=-mj_max,mj_max-1
        d_number_states_mj((i_count*nstates+1):(i_count+1)*nstates) =
     &  -d_mjMax+dble(i_count)
        i_count = i_count + 1
        i_col=1
        do mj2=-mj_max,mj_max-1
          do i=1,nstates
            j_row=j_start
            do j=1,nstates
              bb_mjj(i_col,j_row)=(0.5d0,0.d0)*
     &        (shifted_norm_minus(i,j,mj1,mj2)-
     &        shifted_norm_plus(i,j,mj1,mj2))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
        enddo
        j_start=j_row
      enddo

      bb_mjj=bb_mjj*(0,-1)*dthetadt

      deallocate(shifted_norm_minus)
      deallocate(shifted_norm_plus)

      return
      end

      subroutine EM_unrolling(wave_n,v2sbj1,unrolling1,nkap,nm,nstates,
     &k1,i_state,ll_start,ll_step)
      include 'inc.par'
      real*8 wave_n(nstates,2*nm,-nkap:nkap),v2sbj1(nm,nm,0:2*nkap),
     &unrolling1(nm,0:2*nkap)

      unrolling1=0.d0
      do L=ll_start,2*nkap,ll_step
        do i=1,nm
          do j=max(1,i-ns+1),min(i+ns-1,nm)
            unrolling1(i,L)=unrolling1(i,L)+wave_n(i_state,j+nm,k1)*
     &      v2sbj1(i,j,L)
          enddo
        enddo
      enddo
      return
      end

      subroutine EM_unrolling_dkb_1(wave_n,v2sbj1,v2sbj2,
     &unrolling1,nkap,nm,nstates,k1,i_state,ll_start,ll_step)
      include 'inc.par'
      real*8 wave_n(nstates,2*nm,-nkap:nkap),
     &v2sbj1(nm,nm,0:2*nkap),v2sbj2(-nkap:nkap,nm,nm,0:2*nkap),
     &unrolling1(nm,0:2*nkap)

      unrolling1=0.d0
      do L=ll_start,2*nkap,ll_step
        do i=1,nm
          do j=max(1,i-ns+1),min(i+ns-1,nm)
            unrolling1(i,L)=unrolling1(i,L)+
     &      wave_n(i_state,j+nm,k1)*v2sbj1(i,j,L)+
     &      wave_n(i_state,j,k1)*v2sbj2(k1,j,i,L)
          enddo
        enddo
      enddo

      return
      end

      subroutine EM_unrolling_dkb_2(wave_n,v2sbj3,v2sbj4,unrolling2,
     &nkap,nm,nstates,k1,i_state,k2,ll_start,ll_step)
      include 'inc.par'
      real*8 wave_n(nstates,2*nm,-nkap:nkap),
     &v2sbj3(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj4(-nkap:nkap,-nkap:nkap,nm,nm,0:2*nkap),
     &unrolling2(nm,0:2*nkap)

      unrolling2=0.d0
      do L=ll_start,2*nkap,ll_step
        do i=1,nm
          do j=max(1,i-ns+1),min(i+ns-1,nm)
            unrolling2(i,L)=unrolling2(i,L)+
     &      wave_n(i_state,j+nm,k1)*v2sbj3(k2,i,j,L)+
     &      wave_n(i_state,j,k1)*v2sbj4(k1,k2,j,i,L)
          enddo
        enddo
      enddo

      return
      end

      subroutine EM_storage_mat(wave_n,stor_em_int_mat1,nkap,nm,nstates,
     &k2,j_state,unrolling1,ll_start,ll_step)
      real*8 wave_n(nstates,2*nm,-nkap:nkap),
     &stor_em_int_mat1(0:2*nkap),unrolling1(nm,0:2*nkap)

      stor_em_int_mat1=0.d0
      do L=ll_start,2*nkap,ll_step
        do i=1,nm
          stor_em_int_mat1(L)=stor_em_int_mat1(L)+wave_n(j_state,i,k2)*
     &    unrolling1(i,L)
        enddo
      enddo

      return
      end

      subroutine EM_storage_mat_dkb(wave_n,stor_em_int_mat1,nkap,nm,
     &nstates,k2,j_state,unrolling1,unrolling2,ll_start,ll_step)
      real*8 wave_n(nstates,2*nm,-nkap:nkap),
     &stor_em_int_mat1(0:2*nkap),unrolling1(nm,0:2*nkap),
     &unrolling2(nm,0:2*nkap)

      stor_em_int_mat1=0.d0
      do L=ll_start,2*nkap,ll_step
        do i=1,nm
          stor_em_int_mat1(L)=stor_em_int_mat1(L)+
     &    wave_n(j_state,i,k2)*unrolling1(i,L)+
     &    wave_n(j_state,i+nm,k2)*unrolling2(i,L)
        enddo
      enddo

      return
      end

      subroutine Electro_Mag_Field(wave_n,
     &nm,nkap,dtdxi,nsts,nsts_mj,bb_mjj,v2sbj1,
     &v2sbj2,v2sbj3,v2sbj4,Amplitude,wavenumber,elapsed_time)
      include 'inc.par'
      real*8 wave_n(nsts,2*nm,-nkap:nkap),
     &v2sbj1(nm,nm,0:2*nkap),v2sbj2(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj3(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj4(-nkap:nkap,-nkap:nkap,nm,nm,0:2*nkap)
      complex*16 bb_mjj(nsts_mj,nsts_mj)
      complex*16, dimension(:,:,:,:),allocatable:: em_int_mat1,
     &em_int_mat2
      real*8, dimension(:,:,:,:), allocatable:: ang_bit_plus,
     &ang_bit_minus
      real*8, dimension(:,:), allocatable:: unrolling1,unrolling2
      real*8, dimension(:), allocatable:: stor_em_int_mat1
      real*8, dimension(:,:,:), allocatable:: ang_sigx_temp_p,
     &ang_sigx_temp_m
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      logical dkb

      mj_max=nint(amj_max+0.5d0)

c      write(*,*)stat_in_cont

      allocate(ang_bit_plus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))
      allocate(ang_bit_minus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))

      allocate(ang_sigx_temp_p(-nkap:nkap,0:2*nkap,-nkap:nkap))
      allocate(ang_sigx_temp_m(-nkap:nkap,0:2*nkap,-nkap:nkap))

      ang_bit_minus=0.d0
      ang_bit_plus=0.d0

      do mu=-mj_max,mj_max-1
        aamu=dble(mu)+0.5d0
        call ang_sigxx_plus(ang_sigx_temp_p,nkap,aamu)
        call ang_sigxx_minus(ang_sigx_temp_m,nkap,aamu)
        do k1=-nkap,nkap
          do k2=-nkap,nkap
            if(k1*k2.ne.0)then
              do L=0,2*nkap
                ang_bit_plus(k1,L,k2,mu)=ang_sigx_temp_p(k1,L,k2)
                ang_bit_minus(k1,L,k2,mu)=ang_sigx_temp_m(k1,L,k2)
              enddo
            endif
          enddo
        enddo
      enddo

      deallocate(ang_sigx_temp_p)
      deallocate(ang_sigx_temp_m)

      allocate(em_int_mat1(nsts,nsts,-mj_max:mj_max-1,-mj_max:mj_max-1))
      em_int_mat1=0.d0

      do i=1,nsts
        do j=1,nsts
          do mj=-mj_max+1,mj_max-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            do ki=-nkap,nkap
              a_jk1=1.d0*abs(ki)-0.5d0
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb .and. ki.ne.0)then
                call EM_unrolling_dkb_1(wave_n,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nsts,ki,i,0,1)
              elseif(.not.dkb .and. ki.ne.0)then
                call EM_unrolling(wave_n,v2sbj1,
     &          unrolling1,nkap,nm,nsts,ki,i,0,1)
              endif
              do kj=-nkap,nkap
                a_jk2=1.d0*abs(kj)-0.5d0
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1 .and. ki*kj.ne.0)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_n,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nsts,ki,i,kj,0,1)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_n,stor_em_int_mat1,
     &              nkap,nm,nsts,kj,j,unrolling1,unrolling2,0,1)
                  else
                    call EM_storage_mat(wave_n,stor_em_int_mat1,nkap,
     &              nm,nsts,kj,j,unrolling1,0,1)
                  endif
                  do L=0,2*nkap
                    if(dabs(ang_bit_minus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_minus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat1=em_int_mat1*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_minus)
      allocate(em_int_mat2(nsts,nsts,-mj_max:mj_max-1,-mj_max:mj_max-1))
      em_int_mat2=0.d0

      do i=1,nsts
        do j=1,nsts
          do mj=-mj_max,mj_max-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            do ki=-nkap,nkap
              a_jk1=1.d0*abs(ki)-0.5d0
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb .and. ki.ne.0)then
                call EM_unrolling_dkb_1(wave_n,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nsts,ki,i,0,1)
              elseif(.not.dkb .and. ki.ne.0)then
                call EM_unrolling(wave_n,v2sbj1,
     &          unrolling1,nkap,nm,nsts,ki,i,0,1)
              endif
              do kj=-nkap,nkap
                a_jk2=1.d0*abs(kj)-0.5d0
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1 .and. ki*kj.ne.0)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_n,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nsts,ki,i,kj,0,1)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_n,stor_em_int_mat1,
     &              nkap,nm,nsts,kj,j,unrolling1,unrolling2,0,1)
                  else
                    call EM_storage_mat(wave_n,stor_em_int_mat1,nkap,
     &              nm,nsts,kj,j,unrolling1,0,1)
                  endif
                  do L=0,2*nkap
                    if(dabs(ang_bit_plus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_plus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat2=em_int_mat2*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_plus)

      j_start=1
      do mj1=-mj_max,mj_max-1
        i_col=1
        do mj2=-mj_max,mj_max-1
          do i=1,nsts
            j_row=j_start
            do j=1,nsts
              bb_mjj(i_col,j_row)=
     &        (em_int_mat1(i,j,mj1,mj2)+
     &        em_int_mat2(i,j,mj1,mj2))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
        enddo
        j_start=j_row
      enddo

      bb_mjj=bb_mjj*dtdxi

      deallocate(em_int_mat1)
      deallocate(em_int_mat2)

      return
      end

      subroutine Electro_Mag_Field_even(wave_n,
     &nm,nkap,dtdxi,nsts,nsts_mj,bb_mjj,v2sbj1,
     &v2sbj2,v2sbj3,v2sbj4,Amplitude,wavenumber,elapsed_time)
      include 'inc.par'
      real*8 wave_n(nsts,2*nm,-nkap:nkap),
     &v2sbj1(nm,nm,0:2*nkap),v2sbj2(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj3(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj4(-nkap:nkap,-nkap:nkap,nm,nm,0:2*nkap)
      complex*16 bb_mjj(nsts_mj,nsts_mj)
      complex*16, dimension(:,:,:,:),allocatable:: em_int_mat1,
     &em_int_mat2
      real*8, dimension(:,:,:,:), allocatable:: ang_bit_plus,
     &ang_bit_minus
      real*8, dimension(:,:), allocatable:: unrolling1,unrolling2
      real*8, dimension(:), allocatable:: stor_em_int_mat1
      real*8, dimension(:,:,:), allocatable:: ang_sigx_temp_p,
     &ang_sigx_temp_m
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      logical dkb

      mj_max=nint(amj_max+0.5d0)

c      write(*,*)stat_in_cont

      allocate(ang_bit_plus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))
      allocate(ang_bit_minus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))

      allocate(ang_sigx_temp_p(-nkap:nkap,0:2*nkap,-nkap:nkap))
      allocate(ang_sigx_temp_m(-nkap:nkap,0:2*nkap,-nkap:nkap))

      ang_bit_minus=0.d0
      ang_bit_plus=0.d0

      do mu=-mj_max,mj_max-1
        aamu=dble(mu)+0.5d0
        call ang_sigxx_plus(ang_sigx_temp_p,nkap,aamu)
        call ang_sigxx_minus(ang_sigx_temp_m,nkap,aamu)
        do k1=-nkap,nkap
          do k2=-nkap,nkap
            if(k1*k2.ne.0)then
              do L=1,2*nkap,2
                ang_bit_plus(k1,L,k2,mu)=ang_sigx_temp_p(k1,L,k2)
                ang_bit_minus(k1,L,k2,mu)=ang_sigx_temp_m(k1,L,k2)
              enddo
            endif
          enddo
        enddo
      enddo
      deallocate(ang_sigx_temp_p)
      deallocate(ang_sigx_temp_m)

      allocate(em_int_mat1(nsts,nsts,-mj_max:mj_max-1,-mj_max:mj_max-1))
      em_int_mat1=0.d0

      do i=1,nsts
        do j=1,nsts
          do mj=-mj_max+1,mj_max-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            k1fact=(-1)**nkap
            do k1=nkap,1,-1
              a_jk1=1.d0*k1-0.5d0
              ki=k1*k1fact
              k1fact=-k1fact
              k2fact=(-1)**nkap
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb)then
                call EM_unrolling_dkb_1(wave_n,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              else
                call EM_unrolling(wave_n,v2sbj1,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              endif
              do k2=nkap,1,-1
                a_jk2=1.d0*k2-0.5d0
                kj=k2*k2fact
                k2fact=-k2fact
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_n,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nsts,ki,i,kj,1,2)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_n,stor_em_int_mat1,
     &              nkap,nm,nsts,kj,j,unrolling1,unrolling2,1,2)
                  else
                    call EM_storage_mat(wave_n,stor_em_int_mat1,nkap,
     &              nm,nsts,kj,j,unrolling1,1,2)
                  endif
                  do L=1,2*nkap,2
                    if(dabs(ang_bit_minus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_minus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat1=em_int_mat1*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_minus)
      allocate(em_int_mat2(nsts,nsts,-mj_max:mj_max-1,-mj_max:mj_max-1))
      em_int_mat2=0.d0

      do i=1,nsts
        do j=1,nsts
          do mj=-mj_max,mj_max-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            k1fact=(-1)**nkap
            do k1=nkap,1,-1
              a_jk1=1.d0*k1-0.5d0
              ki=k1*k1fact
              k1fact=-k1fact
              k2fact=(-1)**nkap
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb)then
                call EM_unrolling_dkb_1(wave_n,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              else
                call EM_unrolling(wave_n,v2sbj1,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              endif
              do k2=nkap,1,-1
                a_jk2=1.d0*k2-0.5d0
                kj=k2*k2fact
                k2fact=-k2fact
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_n,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nsts,ki,i,kj,1,2)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_n,stor_em_int_mat1,
     &              nkap,nm,nsts,kj,j,unrolling1,unrolling2,1,2)
                  else
                    call EM_storage_mat(wave_n,stor_em_int_mat1,nkap,
     &              nm,nsts,kj,j,unrolling1,1,2)
                  endif
                  do L=1,2*nkap,2
                    if(dabs(ang_bit_plus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_plus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat2=em_int_mat2*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_plus)

      j_start=1
      do mj1=-mj_max,mj_max-1
        i_col=1
        do mj2=-mj_max,mj_max-1
          do i=1,nsts
            j_row=j_start
            do j=1,nsts
              bb_mjj(i_col,j_row)=
     &        (em_int_mat1(i,j,mj1,mj2)+
     &        em_int_mat2(i,j,mj1,mj2))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
        enddo
        j_start=j_row
      enddo

      bb_mjj=bb_mjj*dtdxi

      deallocate(em_int_mat1)
      deallocate(em_int_mat2)

      return
      end

      subroutine Electro_Mag_Field_odd(wave_n,
     &nm,nkap,dtdxi,nsts,nsts_mj,bb_mjj,v2sbj1,
     &v2sbj2,v2sbj3,v2sbj4,Amplitude,wavenumber,elapsed_time)
      include 'inc.par'
      real*8 wave_n(nsts,2*nm,-nkap:nkap),
     &v2sbj1(nm,nm,0:2*nkap),v2sbj2(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj3(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj4(-nkap:nkap,-nkap:nkap,nm,nm,0:2*nkap)
      complex*16 bb_mjj(nsts_mj,nsts_mj)
      complex*16, dimension(:,:,:,:),allocatable:: em_int_mat1,
     &em_int_mat2
      real*8, dimension(:,:,:,:), allocatable:: ang_bit_plus,
     &ang_bit_minus
      real*8, dimension(:,:), allocatable:: unrolling1,unrolling2
      real*8, dimension(:), allocatable:: stor_em_int_mat1
      real*8, dimension(:,:,:), allocatable:: ang_sigx_temp_p,
     &ang_sigx_temp_m
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      logical dkb

      mj_max=nint(amj_max+0.5d0)

c      write(*,*)stat_in_cont

      allocate(ang_bit_plus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))
      allocate(ang_bit_minus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))

      allocate(ang_sigx_temp_p(-nkap:nkap,0:2*nkap,-nkap:nkap))
      allocate(ang_sigx_temp_m(-nkap:nkap,0:2*nkap,-nkap:nkap))

      ang_bit_minus=0.d0
      ang_bit_plus=0.d0

      do mu=-mj_max,mj_max-1
        aamu=dble(mu)+0.5d0
        call ang_sigxx_plus(ang_sigx_temp_p,nkap,aamu)
        call ang_sigxx_minus(ang_sigx_temp_m,nkap,aamu)
        do k1=-nkap,nkap
          do k2=-nkap,nkap
            if(k1*k2.ne.0)then
              do L=1,2*nkap,2
                ang_bit_plus(k1,L,k2,mu)=ang_sigx_temp_p(k1,L,k2)
                ang_bit_minus(k1,L,k2,mu)=ang_sigx_temp_m(k1,L,k2)
              enddo
            endif
          enddo
        enddo
      enddo

      deallocate(ang_sigx_temp_p)
      deallocate(ang_sigx_temp_m)

      allocate(em_int_mat1(nsts,nsts,-mj_max:mj_max-1,-mj_max:mj_max-1))
      em_int_mat1=0.d0

      do i=1,nsts
        do j=1,nsts
          do mj=-mj_max+1,mj_max-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            k1fact=(-1)**(nkap+1)
            do k1=nkap,1,-1
              a_jk1=1.d0*k1-0.5d0
              ki=k1*k1fact
              k1fact=-k1fact
              k2fact=(-1)**(nkap+1)
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb)then
                call EM_unrolling_dkb_1(wave_n,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              else
                call EM_unrolling(wave_n,v2sbj1,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              endif
              do k2=nkap,1,-1
                a_jk2=1.d0*k2-0.5d0
                kj=k2*k2fact
                k2fact=-k2fact
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_n,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nsts,ki,i,kj,1,2)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_n,stor_em_int_mat1,
     &              nkap,nm,nsts,kj,j,unrolling1,unrolling2,1,2)
                  else
                    call EM_storage_mat(wave_n,stor_em_int_mat1,nkap,
     &              nm,nsts,kj,j,unrolling1,1,2)
                  endif
                  do L=1,2*nkap,2
                    if(dabs(ang_bit_minus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_minus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat1=em_int_mat1*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_minus)
      allocate(em_int_mat2(nsts,nsts,-mj_max:mj_max-1,-mj_max:mj_max-1))
      em_int_mat2=0.d0

      do i=1,nsts
        do j=1,nsts
          do mj=-mj_max,mj_max-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            k1fact=(-1)**(nkap+1)
            do k1=nkap,1,-1
              a_jk1=1.d0*k1-0.5d0
              ki=k1*k1fact
              k1fact=-k1fact
              k2fact=(-1)**(nkap+1)
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb)then
                call EM_unrolling_dkb_1(wave_n,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              else
                call EM_unrolling(wave_n,v2sbj1,
     &          unrolling1,nkap,nm,nsts,ki,i,1,2)
              endif
              do k2=nkap,1,-1
                a_jk2=1.d0*k2-0.5d0
                kj=k2*k2fact
                k2fact=-k2fact
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_n,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nsts,ki,i,kj,1,2)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_n,stor_em_int_mat1,
     &              nkap,nm,nsts,kj,j,unrolling1,unrolling2,1,2)
                  else
                    call EM_storage_mat(wave_n,stor_em_int_mat1,nkap,
     &              nm,nsts,kj,j,unrolling1,1,2)
                  endif
                  do L=1,2*nkap,2
                    if(dabs(ang_bit_plus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_plus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat2=em_int_mat2*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_plus)

      j_start=1
      do mj1=-mj_max,mj_max-1
        i_col=1
        do mj2=-mj_max,mj_max-1
          do i=1,nsts
            j_row=j_start
            do j=1,nsts
              bb_mjj(i_col,j_row)=
     &        (em_int_mat1(i,j,mj1,mj2)+
     &        em_int_mat2(i,j,mj1,mj2))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
        enddo
        j_start=j_row
      enddo

      bb_mjj=bb_mjj*dtdxi

      deallocate(em_int_mat1)
      deallocate(em_int_mat2)

      return
      end

      subroutine Electro_Mag_Field_EO(wave_ne,wave_no,
     &nm,nkap,dtdxi,nstse,nstse_mj,nstso,nstso_mj,bb_mjj,v2sbj1,
     &v2sbj2,v2sbj3,v2sbj4,Amplitude,wavenumber,elapsed_time)
      include 'inc.par'
      real*8 wave_ne(nstse,2*nm,-nkap:nkap),
     &wave_no(nstso,2*nm,-nkap:nkap),
     &v2sbj1(nm,nm,0:2*nkap),v2sbj2(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj3(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj4(-nkap:nkap,-nkap:nkap,nm,nm,0:2*nkap)
      complex*16 bb_mjj(nstse_mj,nstso_mj)
      complex*16, dimension(:,:,:,:),allocatable:: em_int_mat1,
     &em_int_mat2
      real*8, dimension(:,:,:,:), allocatable:: ang_bit_plus,
     &ang_bit_minus
      real*8, dimension(:,:), allocatable:: unrolling1,unrolling2
      real*8, dimension(:), allocatable:: stor_em_int_mat1
      real*8, dimension(:,:,:), allocatable:: ang_sigx_temp_p,
     &ang_sigx_temp_m
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      logical dkb

      mj_max=nint(amj_max+0.5d0)

c      write(*,*)stat_in_cont

      allocate(ang_bit_plus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))
      allocate(ang_bit_minus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-mj_max:mj_max-1))

      allocate(ang_sigx_temp_p(-nkap:nkap,0:2*nkap,-nkap:nkap))
      allocate(ang_sigx_temp_m(-nkap:nkap,0:2*nkap,-nkap:nkap))

      ang_bit_minus=0.d0
      ang_bit_plus=0.d0

      do mu=-mj_max,mj_max-1
        aamu=dble(mu)+0.5d0
        call ang_sigxx_plus(ang_sigx_temp_p,nkap,aamu)
        call ang_sigxx_minus(ang_sigx_temp_m,nkap,aamu)
        do k1=-nkap,nkap
          do k2=-nkap,nkap
            if(k1*k2.ne.0)then
              do L=0,2*nkap,2
                ang_bit_plus(k1,L,k2,mu)=ang_sigx_temp_p(k1,L,k2)
                ang_bit_minus(k1,L,k2,mu)=ang_sigx_temp_m(k1,L,k2)
              enddo
            endif
          enddo
        enddo
      enddo

      deallocate(ang_sigx_temp_p)
      deallocate(ang_sigx_temp_m)

      allocate(em_int_mat1(nstse,nstso,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))
      em_int_mat1=0.d0

      do i=1,nstse
        do j=1,nstso
          do mj=-mj_max+1,mj_max-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            k1fact=(-1)**(nkap)
            do k1=nkap,1,-1
              a_jk1=1.d0*k1-0.5d0
              ki=k1*k1fact
              k1fact=-k1fact
              k2fact=(-1)**(nkap+1)
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb)then
                call EM_unrolling_dkb_1(wave_ne,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nstse,ki,i,0,2)
              else
                call EM_unrolling(wave_ne,v2sbj1,
     &          unrolling1,nkap,nm,nstse,ki,i,0,2)
              endif
              do k2=nkap,1,-1
                a_jk2=1.d0*k2-0.5d0
                kj=k2*k2fact
                k2fact=-k2fact
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_ne,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nstse,ki,i,kj,0,2)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_no,stor_em_int_mat1,
     &              nkap,nm,nstso,kj,j,unrolling1,unrolling2,0,2)
                  else
                    call EM_storage_mat(wave_no,stor_em_int_mat1,
     &              nkap,nm,nstso,kj,j,unrolling1,0,2)
                  endif
                  do L=0,2*nkap,2
                    if(dabs(ang_bit_minus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_minus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat1(i,j,mj,mj2)=
     &                em_int_mat1(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_minus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat1=em_int_mat1*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_minus)
      allocate(em_int_mat2(nstse,nstso,-mj_max:mj_max-1,
     &-mj_max:mj_max-1))
      em_int_mat2=0.d0

      do i=1,nstse
        do j=1,nstso
          do mj=-mj_max,mj_max-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            k1fact=(-1)**(nkap)
            do k1=nkap,1,-1
              a_jk1=1.d0*k1-0.5d0
              ki=k1*k1fact
              k1fact=-k1fact
              k2fact=(-1)**(nkap+1)
              allocate(unrolling1(nm,0:2*nkap))
              if(dkb)then
                call EM_unrolling_dkb_1(wave_ne,v2sbj1,v2sbj2,
     &          unrolling1,nkap,nm,nstse,ki,i,0,2)
              else
                call EM_unrolling(wave_ne,v2sbj1,
     &          unrolling1,nkap,nm,nstse,ki,i,0,2)
              endif
              do k2=nkap,1,-1
                a_jk2=1.d0*k2-0.5d0
                kj=k2*k2fact
                k2fact=-k2fact
                if(dabs(a_mj2).le.a_jk2 .and.
     &          dabs(a_mj).le.a_jk1)then
                  if(dkb)then
                    allocate(unrolling2(nm,0:2*nkap))
                    call EM_unrolling_dkb_2(wave_ne,v2sbj3,v2sbj4,
     &              unrolling2,nkap,nm,nstse,ki,i,kj,0,2)
                  endif
                  allocate(stor_em_int_mat1(0:2*nkap))
                  if(dkb)then
                    call EM_storage_mat_dkb(wave_no,stor_em_int_mat1,
     &              nkap,nm,nstso,kj,j,unrolling1,unrolling2,0,2)
                  else
                    call EM_storage_mat(wave_no,stor_em_int_mat1,
     &              nkap,nm,nstso,kj,j,unrolling1,0,2)
                  endif
                  do L=0,2*nkap,2
                    if(dabs(ang_bit_plus(ki,L,-kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(ki,L,-kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,1))*(0,1)**L
                    endif
                    if(dabs(ang_bit_plus(-ki,L,kj,mj)).ne.
     &              0.d0)then
                      em_int_mat2(i,j,mj,mj2)=
     &                em_int_mat2(i,j,mj,mj2)+
     &                (stor_em_int_mat1(L)*
     &                ang_bit_plus(-ki,L,kj,mj)*
     &                dsqrt(2.d0*L+1.d0)*(0,-1))*(0,1)**L
                    endif
                  enddo
                  deallocate(stor_em_int_mat1)
                  if(dkb)deallocate(unrolling2)
                endif
              enddo
              deallocate(unrolling1)
            enddo
          enddo
        enddo
      enddo

      em_int_mat2=em_int_mat2*dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*
     &elapsed_time)*Amplitude/dsqrt(c)
      deallocate(ang_bit_plus)

      j_start=1
      do mj1=-mj_max,mj_max-1
        i_col=1
        do mj2=-mj_max,mj_max-1
          do i=1,nstse
            j_row=j_start
            do j=1,nstso
              bb_mjj(i_col,j_row)=
     &        (em_int_mat1(i,j,mj1,mj2)+
     &        em_int_mat2(i,j,mj1,mj2))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
        enddo
        j_start=j_row
      enddo

      bb_mjj=bb_mjj*dtdxi

      deallocate(em_int_mat1)
      deallocate(em_int_mat2)

      return
      end

      subroutine Electro_Mag_Field_OE(wave_n,wave_o,
     &nm,nkap,Start_stat,Start_stat_o,Mat_Dimension,
     &stat_in_cont,stat_in_cont_o,dtdxi,nsts,nste,bb_mjj,v2sbj1,v2sbj2,
     &v2sbj3,v2sbj4,Amplitude,wavenumber,elapsed_time)
      include 'inc.par'
      integer Start_stat,Start_stat_o
      real*8 wave_n(nsts,2*nm,-nkap:nkap),
     &wave_o(nste,2*nm,-nkap:nkap),
     &v2sbj1(nm,nm,0:2*nkap),v2sbj2(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj3(-nkap:nkap,nm,nm,0:2*nkap),
     &v2sbj4(-nkap:nkap,-nkap:nkap,nm,nm,0:2*nkap)
      integer stat_in_cont(Start_stat:
     &Start_stat+Mat_dimension-1),stat_in_cont_o(Start_stat_o:
     &Start_stat_o+Mat_dimension-1)
      complex*16 bb_mjj(Mat_dimension*2*nkap,
     &Mat_dimension*2*nkap)
      complex*16, dimension(:,:,:,:),allocatable:: em_int_mat1,
     &em_int_mat2
      real*8, dimension(:,:,:,:,:),allocatable:: stor_em_int_mat1,
     &unrolling2,unrolling2_back,stor_em_int_mat2
      real*8, dimension(:,:,:,:), allocatable:: ang_bit_plus,
     &ang_bit_minus,unrolling1,unrolling1_back
      real*8, dimension(:,:,:), allocatable:: ang_sigx_temp_p,
     &ang_sigx_temp_m
      common /common_dkb/ dkb
      logical dkb

c      write(*,*)stat_in_cont

      allocate(ang_bit_plus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-nkap:nkap-1))
      allocate(ang_bit_minus(-nkap:nkap,0:2*nkap,-nkap:nkap,
     &-nkap:nkap-1))

      allocate(ang_sigx_temp_p(-nkap:nkap,0:2*nkap,-nkap:nkap))
      allocate(ang_sigx_temp_m(-nkap:nkap,0:2*nkap,-nkap:nkap))

      ang_bit_minus=0.d0
      ang_bit_plus=0.d0

      do mu=-nkap,nkap-1
        amu=dble(mu)+0.5d0
        call ang_sigxx_plus(ang_sigx_temp_p,nkap,amu)
        call ang_sigxx_minus(ang_sigx_temp_m,nkap,amu)
        do k1=-nkap,nkap
          do k2=-nkap,nkap
            if(k1*k2.ne.0)then
              do L=0,2*nkap!,2
                ang_bit_plus(k1,L,k2,mu)=ang_sigx_temp_p(k1,L,k2)
                ang_bit_minus(k1,L,k2,mu)=ang_sigx_temp_m(k1,L,k2)
              enddo
            endif
          enddo
        enddo
      enddo

      deallocate(ang_sigx_temp_p)
      deallocate(ang_sigx_temp_m)

      allocate(unrolling1(Mat_dimension,-nkap:nkap,nm,0:2*nkap))
      allocate(unrolling1_back(Mat_dimension,-nkap:nkap,nm,0:2*nkap))
      unrolling1=0.d0
      unrolling1_back=0.d0
      if(dkb)then
        allocate(unrolling2(Mat_dimension,-nkap:nkap,-nkap:nkap,
     &  nm,0:2*nkap))
        allocate(unrolling2_back(Mat_dimension,-nkap:nkap,-nkap:nkap,
     &  nm,0:2*nkap))
        unrolling2=0.d0
        unrolling2_back=0.d0
      endif


      do jj_ev=Start_stat,Start_stat+Mat_dimension-1
        j_ev_begin=jj_ev-Start_stat+1
        j_ev=stat_in_cont(jj_ev)
        do k=-nkap,nkap
          call DiracAngularL_int(k,L_k)
          if(k.ne.0 )then!.and. (1+((-1)**(L_k))).ne.0)then
            do L=0,2*nkap!,2
              do i=1,nm
                if(dkb)then
                  do j=max(1,i-ns+1),min(i+ns-1,nm)
                    unrolling1(j_ev_begin,k,i,L)=
     &              unrolling1(j_ev_begin,k,i,L)+
     &              wave_n(j_ev,j+nm,k)*
     &              v2sbj1(i,j,L)+
     &              wave_n(j_ev,j,k)*
     &              v2sbj2(k,j,i,L)
                  enddo
                else
                  do j=max(1,i-ns+1),min(i+ns-1,nm)
                    unrolling1(j_ev_begin,k,i,L)=
     &              unrolling1(j_ev_begin,k,i,L)+
     &              wave_n(j_ev,j+nm,k)*
     &              v2sbj1(i,j,L)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo
      if(dkb)then
        do jj_ev=Start_stat,Start_stat+Mat_dimension-1
          j_ev_begin=jj_ev-Start_stat+1
          j_ev=stat_in_cont(jj_ev)
          do k1=-nkap,nkap
            call DiracAngularL_int(k1,L_k1)
            do k2=-nkap,nkap
              call DiracAngularL_int(k2,L_k2)
              if(k1*k2.ne.0 )then!.and. (1+((-1)**(L_k2))).ne.0 .and.
c     &      (1+((-1)**(L_k1))).eq.0)then
                do L=0,2*nkap!,2
                  do i=1,nm
                    do j=max(1,i-ns+1),min(i+ns-1,nm)
                      unrolling2(j_ev_begin,k1,k2,i,L)=
     &                unrolling2(j_ev_begin,k1,k2,i,L)+
     &                wave_n(j_ev,j+nm,k2)*
     &                v2sbj3(k1,i,j,L)+
     &                wave_n(j_ev,j,k2)*
     &                v2sbj4(k2,k1,j,i,L)
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      endif
      do jj_ev=Start_stat,Start_stat+Mat_dimension-1
        j_ev_begin=jj_ev-Start_stat+1
        j_ev=stat_in_cont(jj_ev)
        do k=-nkap,nkap
          call DiracAngularL_int(k,L_k)
          if(k.ne.0 )then!.and. (1+((-1)**(L_k))).eq.0)then
            do L=0,2*nkap!,2
              do i=1,nm
                if(dkb)then
                  do j=max(1,i-ns+1),min(i+ns-1,nm)
                    unrolling1_back(j_ev_begin,k,i,L)=
     &              unrolling1_back(j_ev_begin,k,i,L)+
     &              wave_n(j_ev,j,k)*
     &              v2sbj1(i,j,L)+
     &              wave_n(j_ev,j+nm,k)*
     &              v2sbj3(k,j,i,L)
                  enddo
                else
                  do j=max(1,i-ns+1),min(i+ns-1,nm)
                    unrolling1_back(j_ev_begin,k,i,L)=
     &              unrolling1_back(j_ev_begin,k,i,L)+
     &              wave_n(j_ev,j,k)*
     &              v2sbj1(i,j,L)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo
      if(dkb)then
        do jj_ev=Start_stat,Start_stat+Mat_dimension-1
          j_ev_begin=jj_ev-Start_stat+1
          j_ev=stat_in_cont(jj_ev)
          do k1=-nkap,nkap
            call DiracAngularL_int(k1,L_k1)
            do k2=-nkap,nkap
              call DiracAngularL_int(k2,L_k2)
              if(k1*k2.ne.0 )then!.and. (1+((-1)**(L_k2))).eq.0 .and.
c     &      (1+((-1)**(L_k1))).ne.0)then
                do L=0,2*nkap!,2
                  do i=1,nm
                    do j=max(1,i-ns+1),min(i+ns-1,nm)
                      unrolling2_back(j_ev_begin,k1,k2,i,L)=
     &                unrolling2_back(j_ev_begin,k1,k2,i,L)+
     &                wave_n(j_ev,j,k2)*
     &                v2sbj2(k1,i,j,L)+
     &                wave_n(j_ev,j+nm,k2)*
     &                v2sbj4(k1,k2,i,j,L)
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      endif

      allocate(stor_em_int_mat1(Mat_dimension,Mat_dimension,
     &-nkap:nkap,-nkap:nkap,0:2*nkap))

      stor_em_int_mat1=0.d0
      do ii_ev=Start_stat_o,Start_stat_o+Mat_dimension-1
        i_ev_begin=ii_ev-Start_stat_o+1
        i_ev=stat_in_cont_o(ii_ev)
        do jj_ev=Start_stat,Start_stat+Mat_dimension-1
          j_ev_begin=jj_ev-Start_stat+1
          j_ev=stat_in_cont(jj_ev)
          do kk1=-nkap,nkap
            call DiracAngularL_int(kk1,L_kk1)
            do kk2=-nkap,nkap
              call DiracAngularL_int(kk2,L_kk2)
              if(kk1*kk2.ne.0 )then!.and. (1+((-1)**(L_kk1))).eq.0 .and.
c     &         (1+((-1)**(L_kk2))).ne.0)then
                do L=0,2*nkap!,2
                  if(dkb)then
                    do i=1,nm
                      stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)=
     &                stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)+
     &                wave_o(i_ev,i,kk1)*
     &                (unrolling1(j_ev_begin,kk2,i,L))+
     &                wave_o(i_ev,i+nm,kk1)*(
     &                unrolling2(j_ev_begin,kk1,kk2,i,L))
                    enddo
                  else
                    do i=1,nm
                      stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)=
     &                stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)+
     &                wave_o(i_ev,i,kk1)*unrolling1(j_ev_begin,kk2,i,L)
                    enddo
                  endif
c                  if(.not. dkb)then
c                  write(988,*)i_ev_begin,j_ev_begin,kk1,kk2,L
c                  write(988,*)
c     &            stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)
c                  else
c                  write(989,*)i_ev_begin,j_ev_begin,kk1,kk2,L
c                  write(989,*)
c     &            stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)
c                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo

      allocate(stor_em_int_mat2(Mat_dimension,Mat_dimension,
     & -nkap:nkap,-nkap:nkap,0:2*nkap))

      stor_em_int_mat2=0.d0
      do ii_ev=Start_stat_o,Start_stat_o+Mat_dimension-1
        i_ev_begin=ii_ev-Start_stat_o+1
        i_ev=stat_in_cont_o(ii_ev)
        do jj_ev=Start_stat,Start_stat+Mat_dimension-1
          j_ev_begin=jj_ev-Start_stat+1
          j_ev=stat_in_cont(jj_ev)
          do kk1=-nkap,nkap
            call DiracAngularL_int(kk1,L_kk1)
            do kk2=-nkap,nkap
              call DiracAngularL_int(kk2,L_kk2)
              if(kk1*kk2.ne.0 )then!.and. (1+((-1)**(L_kk1))).ne.0 .and.
c     &         (1+((-1)**(L_kk2))).eq.0)then
                do L=0,2*nkap!,2
                  if(dkb)then
                    do i=1,nm
                      stor_em_int_mat2(i_ev_begin,j_ev_begin,kk1,kk2,L)=
     &                stor_em_int_mat2(i_ev_begin,j_ev_begin,kk1,kk2,L)+
     &                wave_o(i_ev,i+nm,kk1)*(unrolling1_back
     &                (j_ev_begin,kk2,i,L))
     &                +
     &                wave_o(i_ev,i,kk1)*(
     &                unrolling2_back(j_ev_begin,kk1,kk2,i,L))
                    enddo
                  else
                    do i=1,nm
                      stor_em_int_mat2(i_ev_begin,j_ev_begin,kk1,kk2,L)=
     &                stor_em_int_mat2(i_ev_begin,j_ev_begin,kk1,kk2,L)+
     &                wave_o(i_ev,i+nm,kk1)*unrolling1_back
     &                (j_ev_begin,kk2,i,L)
                    enddo
                  endif
c                  if(.not. dkb)then
c                  write(988,*)i_ev_begin,j_ev_begin,kk1,kk2,L
c                  write(988,*)
c     &            stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)
c                  else
c                  write(989,*)i_ev_begin,j_ev_begin,kk1,kk2,L
c                  write(989,*)
c     &            stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)
c                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo

      deallocate(unrolling1)
      deallocate(unrolling1_back)
      if(dkb)then
        deallocate(unrolling2)
        deallocate(unrolling2_back)
      endif

      allocate(em_int_mat1(Mat_dimension,Mat_dimension,
     &-nkap:nkap-1,-nkap:nkap-1))
      allocate(em_int_mat2(Mat_dimension,Mat_dimension,
     &-nkap:nkap-1,-nkap:nkap-1))

      em_int_mat1=0.d0
      em_int_mat2=0.d0
      bb_mjj=0.d0

      do ii_ev=Start_stat_o,Start_stat_o+Mat_dimension-1
        i_ev_begin=ii_ev-Start_stat_o+1
        i_ev=stat_in_cont_o(ii_ev)
        do jj_ev=Start_stat,Start_stat+Mat_dimension-1
          j_ev_begin=jj_ev-Start_stat+1
          j_ev=stat_in_cont(jj_ev)
          do mj=-nkap+1,nkap-1
            a_mj=1.d0*mj+0.5d0
            mj2=mj-1
            a_mj2=1.d0*mj2+0.5d0
            do kk1=-nkap,nkap
              call DiracAngularL_int(kk1,L_kk1)
              a_jk1=1.d0*abs(kk1)-0.5d0
              do kk2=-nkap,nkap
                call DiracAngularL_int(kk2,L_kk2)
                a_jk2=1.d0*abs(kk2)-0.5d0
                if(kk1*kk2.ne.0 .and. dabs(a_mj2).le.a_jk2
     &          .and. dabs(a_mj).le.a_jk1)then !.and.
c     &            (1+((-1)**(L_kk1))).eq.0 .and.
c     &            (1+((-1)**(L_kk2))).ne.0)then
                  do L=0,2*nkap!,2
c                    write(*,*)kk1,L,kk2,mj,ang_bit_minus(kk1,L,-kk2,mj)
                    em_int_mat1(i_ev_begin,j_ev_begin,mj,mj2)=
     &              em_int_mat1(i_ev_begin,j_ev_begin,mj,mj2)+
     &              (stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)*
     &              ang_bit_minus(kk1,L,-kk2,mj)*
     &              dsqrt(2.d0*L+1.d0)*(0,1)+
     &              stor_em_int_mat2(i_ev_begin,j_ev_begin,kk1,kk2,L)*
     &              ang_bit_minus(-kk1,L,kk2,mj)*
     &              dsqrt(2.d0*L+1.d0)*(0,-1))
     &              *(0,1)**L
                   enddo
                 endif
               enddo
             enddo
c               if(jj_ev.gt.ii_ev)then
c               em_int_mat1(j_ev_begin,i_ev_begin,mj,mj2)=
c     &         em_int_mat1(i_ev_begin,j_ev_begin,mj,mj2)
c               endif
c            write(*,*)i_ev_begin,j_ev_begin,mj,mj2
c            write(*,*)em_int_mat1(i_ev_begin,j_ev_begin,mj,mj2)
c            pause
          enddo
        enddo
      enddo


      em_int_mat1=em_int_mat1*
     &dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*elapsed_time)*
     &Amplitude

      do ii_ev=Start_stat_o,Start_stat_o+Mat_dimension-1
        i_ev_begin=ii_ev-Start_stat_o+1
        i_ev=stat_in_cont_o(ii_ev)
        do jj_ev=Start_stat,Start_stat+Mat_dimension-1
          j_ev_begin=jj_ev-Start_stat+1
          j_ev=stat_in_cont(jj_ev)
          do mj=-nkap,nkap-2
            a_mj=1.d0*mj+0.5d0
            mj2=mj+1
            a_mj2=1.d0*mj2+0.5d0
            do kk1=-nkap,nkap
              call DiracAngularL_int(kk1,L_kk1)
              a_jk1=1.d0*abs(kk1)-0.5d0
              do kk2=-nkap,nkap
                call DiracAngularL_int(kk2,L_kk2)
                a_jk2=1.d0*abs(kk2)-0.5d0
                if(kk1*kk2.ne.0 .and. dabs(a_mj2).le.a_jk2
     &          .and. dabs(a_mj).le.a_jk1 )then!.and.
c     &            (1+((-1)**(L_kk1))).eq.0 .and.
c     &            (1+((-1)**(L_kk2))).ne.0)then
                  do L=0,2*nkap!,2
                    em_int_mat2(i_ev_begin,j_ev_begin,mj,mj2)=
     &              em_int_mat2(i_ev_begin,j_ev_begin,mj,mj2)+
     &              (stor_em_int_mat1(i_ev_begin,j_ev_begin,kk1,kk2,L)*
     &              ang_bit_plus(kk1,L,-kk2,mj)*(0,1)*
     &              dsqrt(2.d0*L+1.d0)+
     &              stor_em_int_mat2(i_ev_begin,j_ev_begin,kk1,kk2,L)*
     &              ang_bit_plus(-kk1,L,kk2,mj)*(0,-1)*
     &              dsqrt(2.d0*L+1.d0))*
     &              (0,1)**L
                  enddo
                endif
              enddo
            enddo
c               if(jj_ev.gt.ii_ev)then
c               em_int_mat2(j_ev_begin,i_ev_begin,mj,mj2)=
c     &         em_int_mat2(i_ev_begin,j_ev_begin,mj,mj2)
c               endif
          enddo
        enddo
      enddo

      em_int_mat2=em_int_mat2*
     &dsqrt(4.d0*pi)*cdexp((0,-1)*wavenumber*elapsed_time)*
     &Amplitude

      deallocate(ang_bit_plus)
      deallocate(ang_bit_minus)

      lstart=1
      lstep=1

c      if(z_nuc1.eq.z_nuc2)then
c      Mat_dimension22=2*Mat_dimension
c      endif

c      i=1
      do i=1,Mat_dimension
c      j=1
        do j=1,Mat_dimension
          i_col=lstart+(i-1)*2*nkap*lstep
          do mj1=-nkap,nkap-1
            j_row=lstart+(j-1)*2*nkap*lstep
            do mj2=-nkap,nkap-1
              bb_mjj(i_col,j_row)=
     &        (em_int_mat1(i,j,mj1,mj2)+
     &        em_int_mat2(i,j,mj1,mj2))
c               if(.not. dkb)then
c               write(788,*)i_col,j_row,bb_mjj(i_col,j_row)
c               else
c               write(790,*)i_col,j_row,bb_mjj(i_col,j_row)
c               endif
c               pause
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
c         j=j+1
        enddo
c      i=i+1
      enddo
cc      write(*,*)i,j
c      write(*,*)i_col,j_row!,Mat_Dimension*2*nkap*lstep

      bb_mjj=bb_mjj*dtdxi

      deallocate(em_int_mat1)
      deallocate(em_int_mat2)

      return
      end

      subroutine dvdr_matrix_inkl_mj(mm,nste,nstates,mm_mj,
     &nmj)
      include 'inc.par'
      complex*16 mm(nste,nste,nmj),mm_mj(nstates,nstates)
      common /momentum_projection/ amu,amj_max

      mj_max=nint(amj_max+0.5d0)

      mm_mj=(0,0)
      i_step=1
      i_col=1
      j_row=1
      do mj1=-mj_max,mj_max
        mjj=abs(mj1)
        if(mj1.ne.0)then
          do i=1,nste
            j_row=i_step
            do j=1,nste
              mm_mj(i_col,j_row)=mm(i,j,abs(mj1))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
          i_step=i_col
        endif
      enddo

      return
      end

      subroutine dvdr_matrix_inkl_mj_real(mm,nste,nstates,numb_mj,
     &nkap,mm_mj,nmj)
      include 'inc.par'
      integer numb_mj(4,nkap)
      real*8 mm(nste,nste,nmj),mm_mj(nstates,nstates)
      common /momentum_projection/ amu,amj_max

      mj_max=nint(amj_max+0.5d0)

      mm_mj=0.d0
      i_step=1
      i_col=1
      do mj1=-mj_max,mj_max
        mjj=abs(mj1)
        if(mj1.ne. 0)then
          do i=numb_mj(1,mjj),numb_mj(2,mjj)
            j_row=i_step
            do j=numb_mj(1,mjj),numb_mj(2,mjj)
              mm_mj(i_col,j_row)=mm(i,j,abs(mj1))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
          i_step=i_col
          do i=numb_mj(3,mjj),numb_mj(4,mjj)
            j_row=i_step
            do j=numb_mj(3,mjj),numb_mj(4,mjj)
              mm_mj(i_col,j_row)=mm(i,j,abs(mj1))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
          i_step=j_row
        endif
      enddo

      return
      end

      subroutine dvdr_matrix_inkl_mj_eo(mm,nste,nstates,
     &mm_mj,nmj)
      include 'inc.par'
      complex*16 mm(nste,nste,nmj),mm_mj(nstates,nstates)
      common /momentum_projection/ amu,amj_max

      mj_max=nint(amj_max+0.5d0)

      mm_mj=(0,0)
      i_step=1
      i_col=1
      j_row=1
      do mj1=-mj_max,mj_max
        mjj=abs(mj1)
        if(mj1.ne.0)then
          do i=1,nste
            j_row=i_step
            do j=1,nste
              mm_mj(i_col,j_row)=mm(i,j,abs(mj1))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
        endif
        i_step=i_col
      enddo

      return
      end

      subroutine dvdr_matrix_inkl_mj_eo_real(mm,nste,nstates,numb_mj,
     &nkap,mm_mj,nmj)
      include 'inc.par'
      integer numb_mj(2,nkap)
      real*8 mm(nste,nste,nmj),mm_mj(nstates,nstates)
      common /momentum_projection/ amu,amj_max

      mj_max=nint(amj_max+0.5d0)

      mm_mj=0.d0
      i_step=1
      i_col=1
      do mj1=-mj_max,mj_max
        mjj=abs(mj1)
        if(mj1.ne. 0)then
          do i=numb_mj(1,mjj),numb_mj(2,mjj)
            j_row=i_step
            do j=numb_mj(1,mjj),numb_mj(2,mjj)
              mm_mj(i_col,j_row)=mm(i,j,abs(mj1))
              j_row=j_row+1
            enddo
            i_col=i_col+1
          enddo
          i_step=j_row
        endif
      enddo

      return
      end

      subroutine spherical_bessel_j(x,sbj,nkap)
      include 'inc.par'
      real*8 sbj(0:2*nkap)

      Lmax=2*nkap
      if(x.ge.1.d-3)then
        sbj(0)=dsin(x)/x
        sbj(1)=dsin(x)/x/x-dcos(x)/x
        do ll=1,Lmax-1
          sbj(ll+1)=((2.d0*ll+1.d0)/x)*sbj(ll)-sbj(ll-1)
        enddo
      else
        do ll=0,Lmax
          ifk=int_factorial(2*ll+1)
          if(ifk.le.6)then
            sbj(ll)=x**ll/int_factorial(ifk)
          else
            sbj(ll)=0.d0
          endif
        enddo
      endif

      return
      end

      subroutine c_diagonal_general(nm,a,eignum,eigvectr)
      implicit complex*16(a-h,o-z)
      integer nm,ilo,ihi,info,ldvl,lwork
      complex*16 a(nm,nm),eigvectr(nm,nm),eignum(nm)
      real*8 RWORK(2*nm),SCALES(nm,nm),RCONDE(nm),RCONDV(nm),ABRNM
      character*1 JOBVL,JOBVR,SENSE
      complex*16, dimension(:),allocatable:: WORK
      complex*16, dimension(:,:), allocatable:: vl

      LDVL=1
      allocate(vl(ldvl,ldvl))
      allocate(WORK(1))
      SENSE='V'
      JOBVL='N'
      JOBVR='V'
      LWORK=-1

c      CALL ZGEEV(JOBVL, JOBVR, nm, A, nm, eignum, VL, LDVL, eigvectr,
c     $             nm, WORK, LWORK, RWORK, INFO)
c      if(INFO.ne.0)then
c      write(*,*)'INFO NON-ZERO',INFO
c      pause
c      endif

      CALL ZGEEVX('S', JOBVL, JOBVR, SENSE, nm, A, nm, eignum, vl,
     &LDVL,eigvectr,nm, ILO, IHI, SCALES, ABRNM, RCONDE, RCONDV, WORK,
     &LWORK,RWORK, INFO )
      if(INFO.ne.0)then
        write(*,*)'c_diagonal_general 1: INFO NON-ZERO',INFO
        stop
      endif

      LWORK=nint(dble(WORK(1)))
      deallocate(WORK)
      allocate(WORK(LWORK))

      CALL ZGEEVX('S', JOBVL, JOBVR, SENSE, nm, A, nm, eignum, vl,
     &LDVL,eigvectr,nm, ILO, IHI, SCALES, ABRNM, RCONDE, RCONDV, WORK,
     &LWORK, RWORK, INFO )
      if(INFO.ne.0)then
        write(*,*)'c_diagonal_general 2: INFO NON-ZERO',INFO
        stop
      endif
      deallocate(WORK)
      deallocate(vl)

      return
      end

      subroutine c_diagonal_general_option(nm,a,eignum,eigvectr)
      implicit complex*16(a-h,o-z)
      integer nm
      complex*16 a(nm,nm),eigvectr(nm,nm),eignum(nm)
      real*8 RWORK(2*nm)
      character*1 JOBVL,JOBVR
      complex*16, dimension(:),allocatable:: WORK
      complex*16, dimension(:,:), allocatable:: vl

      LDVL=1
      allocate(vl(ldvl,ldvl))
      allocate(WORK(1))
      JOBVL='N'
      JOBVR='V'
      LWORK=-1

      CALL ZGEEV(JOBVL, JOBVR, nm, A, nm, eignum, VL, LDVL, eigvectr,
     &nm, WORK, LWORK, RWORK, INFO)
      if(INFO.ne.0)then
        write(*,*)'c_diagonal_general_option 1: INFO NON-ZERO',INFO
        stop
      endif

      LWORK=nint(dble(WORK(1)))
      deallocate(WORK)
      allocate(WORK(LWORK))

      CALL ZGEEV(JOBVL, JOBVR, nm, A, nm, eignum, VL, LDVL, eigvectr,
     &nm, WORK, LWORK, RWORK, INFO)
      if(INFO.ne.0)then
        write(*,*)'c_diagonal_general_option 2: INFO NON-ZERO',INFO
        stop
      endif
      deallocate(WORK)
      deallocate(vl)

      return
      end

      function int_factorial(n)
      integer int_factorial, n, p
      p = 1
      do i = 1, n
        p = p * i
      end do
      int_factorial = p
      end

      subroutine swapping_of_states(wave_new,nstates,nm,
     &nkap,ii_xi,ixi_stepslower,rmin,rmax,eigval,typ)

      include 'inc.par'
      real*8 wave_new(nstates,2*nm,-nkap:nkap),eigval(nstates)
      integer, dimension(:),allocatable::nearest_neighbours_stor,
     &nn,nn2
      integer, dimension(:,:), allocatable::nearest_neighbours
      real*8, dimension(:),allocatable::t,ro,ro1,ro2,abstand,
     &eigval_sort,eigval_old,eigval_old1
      real*8, dimension(:,:), allocatable::sampl_point,
     &sample_differences_stor,Psi_sum,Psi_sum1,
     &Psi_sum_old,g_sum,g_sum1,g_sum_old,f_sum,f_sum1,f_sum_old,dg_sum,
     &dg_sum1,dg_sum_old,df_sum,df_sum1,df_sum_old,g_sum_old1,
     &f_sum_old1,dg_sum_old1,df_sum_old1,Psi_sum_old1,g_on_f,dg_on_df,
     &g_on_f1,dg_on_df1,g_on_f_old,dg_on_df_old,g_on_f_old1,
     &dg_on_df_old1
      real*8, dimension(:,:,:), allocatable::sample_differences,
     &wave_new_sort
      character*1 typ
      logical dkb!,mkdirs
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      common /Barycentres/ RadiusOne,RadiusTwo
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      !mkdirs=makedirqq('WF')
      CALL SYSTEM("mkdir -p WF")

      a_IND=RadiusOne+RadiusTwo
      allocate(abstand(nstates))
      abstand=rmax*0.5d0/(z_nuc1+z_nuc2)**0.25d0
      if(ii_xi.gt.ixi_stepslower+2)then
        allocate(g_sum_old1(nstates,2*nm))
        allocate(f_sum_old1(nstates,2*nm))
        allocate(dg_sum_old1(nstates,2*nm))
        allocate(df_sum_old1(nstates,2*nm))
        allocate(Psi_sum_old1(nstates,2*nm))
        allocate(eigval_old1(nstates))
        allocate(g_on_f_old1(nstates,2*nm))
        allocate(dg_on_df_old1(nstates,2*nm))
        open(33,file='WF/Psi_sum_prev_'//typ//'_1.dat',status='old')
        open(34,file='WF/Psi_sum_prev_'//typ//'_2.dat',status='old')
        open(35,file='WF/Psi_sum_prev_'//typ//'_3.dat',status='old')
        do i=1,nstates
          read(34,*)eigval_old1(i)
        enddo
        do i=1,nstates
          do n=1,2*nm
            read(33,*)g_sum_old1(i,n),f_sum_old1(i,n),dg_sum_old1(i,n)
            read(34,*)df_sum_old1(i,n),Psi_sum_old1(i,n)
            read(35,*)g_on_f_old1(i,n),dg_on_df_old1(i,n)
          enddo
        enddo
        close(33)
        close(34)
        close(35)
      endif
      if(ii_xi.gt.ixi_stepslower)then
        allocate(g_sum_old(nstates,2*nm))
        allocate(f_sum_old(nstates,2*nm))
        allocate(dg_sum_old(nstates,2*nm))
        allocate(df_sum_old(nstates,2*nm))
        allocate(Psi_sum_old(nstates,2*nm))
        allocate(eigval_old(nstates))
        allocate(g_on_f_old(nstates,2*nm))
        allocate(dg_on_df_old(nstates,2*nm))
        open(33,file='WF/Psi_sum_'//typ//'_1.dat',status='old')
        open(34,file='WF/Psi_sum_'//typ//'_2.dat',status='old')
        open(35,file='WF/Psi_sum_'//typ//'_3.dat',status='old')
        do i=1,nstates
        read(34,*)eigval_old(i)
        enddo
        do i=1,nstates
          do n=1,2*nm
            read(33,*)g_sum_old(i,n),f_sum_old(i,n),dg_sum_old(i,n)
            read(34,*)df_sum_old(i,n),Psi_sum_old(i,n)
            read(35,*)g_on_f_old(i,n),dg_on_df_old(i,n)
          enddo
        enddo
        close(33)
        close(34)
        close(35)
        open(33,file='WF/Psi_sum_prev_'//typ//'_1.dat')
        open(34,file='WF/Psi_sum_prev_'//typ//'_2.dat')
        open(35,file='WF/Psi_sum_prev_'//typ//'_3.dat')
        do i=1,nstates
        write(34,*)eigval_old(i)
        enddo
        do i=1,nstates
          do n=1,2*nm
            write(33,*)g_sum_old(i,n),f_sum_old(i,n),dg_sum_old(i,n)
            write(34,*)df_sum_old(i,n),Psi_sum_old(i,n)
            write(35,*)g_on_f_old(i,n),dg_on_df_old(i,n)
          enddo
        enddo
        close(33)
        close(34)
        close(35)
      endif
      allocate(eigval_sort(nstates))
      allocate(wave_new_sort(nstates,2*nm,-nkap:nkap))
      eigval_sort=eigval
      wave_new_sort=wave_new

      nu=nm+ns+2
      allocate(t(nu))
      call spline_arranger(nu,rmin,rmax,t)
      allocate(sampl_point(2*nm,nstates))
      do jstate=1,nstates
        ii=1
        do i=1,2*nm
          sampl_point(ii,jstate)=(abstand(jstate)/
     &    dabs(abstand(jstate)-1.d0))*
     &    dabs(abstand(jstate)**(dble(ii)/dble(2*nm))-1.d0)+RadiusOne
          ii=ii+1
        enddo
      enddo
      deallocate(abstand)

      Kapp_range=nint(dabs(amu)+0.5d0)

      allocate(Psi_sum(nstates,2*nm))
      allocate(g_sum(nstates,2*nm))
      allocate(f_sum(nstates,2*nm))
      allocate(dg_sum(nstates,2*nm))
      allocate(df_sum(nstates,2*nm))
      allocate(g_on_f(nstates,2*nm))
      allocate(dg_on_df(nstates,2*nm))

      Psi_sum=0.d0
      g_sum=0.d0
      f_sum=0.d0
      dg_sum=0.d0
      df_sum=0.d0
      g_on_f=0.d0
      dg_on_df=0.d0
      do i=1,nstates
        sampl_sum_g=0.d0
        sampl_sum_f=0.d0
c      pause
        do ji=1,2*nm
          Sample_point=sampl_point(ji,i)
          i2=1
          do while(t(i2).lt. Sample_point)
            i2=i2+1
          enddo
          if (t(i2).gt. Sample_point)i2=i2-1
c         write(*,*)Sample_point,sampl_point(abs(ji),i),i2,i
c         pause
          allocate(ro(ns))
          allocate(ro1(ns))
          allocate(ro2(ns))
          call ddsplines(Sample_point,ro,ro1,ro2,i2,t,nu)

          wfgtotal=0.d0
          wfdgtotal=0.d0
          wfgtotal1=0.d0
          wfftotal=0.d0
          wfdftotal=0.d0
          wfftotal1=0.d0
          wfgtotal_tst=0.d0
          wfftotal_tst=0.d0
          wfgtotal1_tst=0.d0
          wfftotal1_tst=0.d0

          do kk=-nkap,nkap
            if(abs(kk).ge. abs(Kapp_range))then

              WFVal_G=0.d0
              WFVal_F=0.d0
              WFVal_DG=0.d0
              WFVal_DF=0.d0

              call DiracAngularL(-1.d0*kk,al1l)
              aj1j=dabs(1.d0*kk)-0.5d0
              ffakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &        spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
              ffakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &        spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
              call DiracAngularL(1.d0*kk,al1l)
              gfakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &        spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
              gfakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &        spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
              call DiracAngularL_int(kk,llk)
              ang_test=coeff_maker(amu,-1.d0*kk,1.d0*kk,llk)+
     &        coeff_maker(amu,1.d0*kk,1.d0*kk,llk)
              do i_spl=max(2,i2-ns+1),min(i2,nm+1)
                if(dkb)then

                  WFVal_G=WFVal_G+
     &            wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)+
     &            wave_new(i,i_spl-1+nm,kk)*0.5d0*
     &            (ro1(i_spl-i2+ns)-kk*ro(i_spl-i2+ns)/Sample_point)
                  WFVal_DG=WFVal_DG+
     &            wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)+
     &            wave_new(i,i_spl-1+nm,kk)*(ro2(i_spl-i2+ns)-
     &            kk*ro1(i_spl-i2+ns)/Sample_point+
     &            kk*ro(i_spl-i2+ns)/Sample_point**2)
     &            /2.d0
                  WFVal_F=WFVal_F+
     &            wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)+
     &            wave_new(i,i_spl-1,kk)*0.5d0*
     &            (ro1(i_spl-i2+ns)+kk*ro(i_spl-i2+ns)/Sample_point)
                  WFVal_DF=WFVal_DF+
     &            wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)+
     &            wave_new(i,i_spl-1,kk)*(ro2(i_spl-i2+ns)+
     &            kk*ro1(i_spl-i2+ns)/Sample_point-
     &            kk*ro(i_spl-i2+ns)/Sample_point**2)
     &            /2.d0

                else

                  WFVal_G=WFVal_G+
     &            wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)
                  WFVal_DG=WFVal_DG+
     &            wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)
                  WFVal_F=WFVal_F+
     &            wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)
                  WFVal_DF=WFVal_DF+
     &            wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)

                endif
              enddo
              wfgtotal=wfgtotal+WFVal_G*gfakt1
              wfdgtotal=wfdgtotal+WFVal_DG*gfakt1
              wfgtotal1=wfgtotal1+WFVal_G*gfakt2
              wfftotal=wfftotal+WFVal_F*ffakt1
              wfdftotal=wfdftotal+WFVal_DF*ffakt1
              wfftotal1=wfftotal1+WFVal_F*ffakt2

              wfgtotal_tst=wfgtotal_tst+WFVal_G*ang_test*gfakt1
              wfftotal_tst=wfftotal_tst+WFVal_F*ang_test*gfakt1
              wfgtotal1_tst=wfgtotal1_tst+WFVal_G*ang_test*gfakt2
              wfftotal1_tst=wfftotal1_tst+WFVal_F*ang_test*gfakt2
            endif
          enddo

          Psi_sum(i,ji)=wfgtotal**2+wfftotal**2+wfgtotal1**2+
     &    wfftotal1**2
c     &      (wfgtotal_tst**2+wfftotal_tst**2+wfgtotal1_tst**2+
c     &      wfftotal1_tst**2)*dsin(Sample_point**2)
          g_sum(i,ji)=dabs(wfgtotal)
          f_sum(i,ji)=dabs(wfftotal)
          dg_sum(i,ji)=dabs(wfdgtotal)
          df_sum(i,ji)=dabs(wfdftotal)
          g_on_f(i,ji)=wfgtotal/wfftotal
          dg_on_df(i,ji)=wfdgtotal/wfftotal-
     &    wfdftotal*wfgtotal/wfftotal**2

          deallocate(ro)
          deallocate(ro1)
          deallocate(ro2)
        enddo
c         g_sum(i)=dabs(g_sum(i))
c         f_sum(i)=dabs(f_sum(i))
c         dg_sum(i)=dabs(dg_sum(i))
c         df_sum(i)=dabs(df_sum(i))
      enddo
      deallocate(t)
      deallocate(sampl_point)

      if(ii_xi.gt.ixi_stepslower)then
        allocate(sample_differences(nstates,nstates,16))
        sample_differences=0.d0
        do i=1,nstates
          do j=i,nstates
            sample_differences(j,i,1)=dabs(eigval(i)-eigval_old(j))**2
            if(ii_xi.gt.ixi_stepslower+2)then
              grad_o=eigval_old(j)-eigval_old1(j)
              sample_differences(j,i,13)=
     &        dabs(eigval(i)-(eigval_old(j)+grad_o))**2
            endif
            do n=1,2*nm
              sample_differences(j,i,2)=sample_differences(j,i,2)+
     &        dabs(Psi_sum(i,n)-Psi_sum_old(j,n))**2
              sample_differences(j,i,3)=sample_differences(j,i,3)+
     &        dabs(g_sum(i,n)-g_sum_old(j,n))**2
              sample_differences(j,i,5)=sample_differences(j,i,5)+
     &        dabs(g_on_f(i,n)-g_on_f_old(j,n))**2
              sample_differences(j,i,4)=sample_differences(j,i,4)+
     &        dabs(f_sum(i,n)-f_sum_old(j,n))**2
              sample_differences(j,i,10)=sample_differences(j,i,10)+
     &        dabs(dg_sum(i,n)-dg_sum_old(j,n))**2
              sample_differences(j,i,11)=sample_differences(j,i,11)+
     &        dabs(df_sum(i,n)-df_sum_old(j,n))**2
              sample_differences(j,i,12)=sample_differences(j,i,12)+
     &        dabs(dg_on_df(i,n)-dg_on_df_old(j,n))**2
              if(ii_xi.gt.ixi_stepslower+2)then
                sample_differences(j,i,6)=sample_differences(j,i,6)+
     &          dabs(Psi_sum(i,n)-Psi_sum_old1(j,n))**2
                sample_differences(j,i,7)=sample_differences(j,i,7)+
     &          dabs(g_sum(i,n)-g_sum_old1(j,n))**2
                sample_differences(j,i,8)=sample_differences(j,i,8)+
     &          dabs(g_on_f(i,n)-g_on_f_old1(j,n))**2
                sample_differences(j,i,9)=sample_differences(j,i,9)+
     &          dabs(f_sum(i,n)-f_sum_old1(j,n))**2
                sample_differences(j,i,14)=sample_differences(j,i,14)+
     &          dabs(dg_sum(i,n)-dg_sum_old1(j,n))**2
                sample_differences(j,i,15)=sample_differences(j,i,15)+
     &          dabs(df_sum(i,n)-df_sum_old1(j,n))**2
                sample_differences(j,i,16)=sample_differences(j,i,16)+
     &          dabs(dg_on_df(i,n)-dg_on_df_old1(j,n))**2
              endif
            enddo
            sample_differences(i,j,:)=sample_differences(j,i,:)
          enddo
        enddo
        allocate(nearest_neighbours(nstates,16))

        do k=1,16
          allocate(sample_differences_stor(nstates,nstates))
          sample_differences_stor=sample_differences(:,:,k)
          allocate(nearest_neighbours_stor(nstates))
          nearest_neighbours_stor=minloc(sample_differences_stor,
     &    dim=1,mask=sample_differences_stor.ge. 0.d0)
          deallocate(sample_differences_stor)
          nearest_neighbours(:,k)=nearest_neighbours_stor
          deallocate(nearest_neighbours_stor)
        enddo

        allocate(g_sum1(nstates,2*nm))
        allocate(f_sum1(nstates,2*nm))
        allocate(dg_sum1(nstates,2*nm))
        allocate(df_sum1(nstates,2*nm))
        allocate(Psi_sum1(nstates,2*nm))
        allocate(g_on_f1(nstates,2*nm))
        allocate(dg_on_df1(nstates,2*nm))
        g_sum1=g_sum
        f_sum1=f_sum
        dg_sum1=dg_sum
        df_sum1=df_sum
        g_on_f1=g_on_f
        dg_on_df1=dg_on_df
        Psi_sum1=Psi_sum

        i7=4
        do i=1,nstates
          allocate(nn(i7))
          allocate(nn2(nstates))
          nn2=0
          nn=nearest_neighbours(i,:)
          do k=1,i7
            nn2(nn(k)) = nn2(nn(k)) + 1
          enddo

          i1=maxloc(nn2,dim=1)
          deallocate(nn)
          deallocate(nn2)

          g_sum1(i1,:)=g_sum(i,:)
          f_sum1(i1,:)=f_sum(i,:)
          dg_sum1(i1,:)=dg_sum(i,:)
          df_sum1(i1,:)=df_sum(i,:)
          g_on_f1(i1,:)=g_on_f(i,:)
          dg_on_df1(i1,:)=dg_on_df(i,:)
          Psi_sum1(i1,:)=Psi_sum(i,:)
          eigval_sort(i1)=eigval(i)
          wave_new_sort(i1,:,:)=wave_new(i,:,:)
        enddo

        g_sum=g_sum1
        f_sum=f_sum1
        dg_sum=dg_sum1
        df_sum=df_sum1
        g_on_f=g_on_f1
        dg_on_df=dg_on_df1
        Psi_sum=Psi_sum1
        eigval=eigval_sort
        wave_new=wave_new_sort

        deallocate(nearest_neighbours)
        deallocate(sample_differences)
        deallocate(wave_new_sort)
        deallocate(g_sum1)
        deallocate(f_sum1)
        deallocate(g_on_f1)
        deallocate(dg_on_df1)
        deallocate(dg_sum1)
        deallocate(df_sum1)
        deallocate(Psi_sum1)
        deallocate(g_sum_old)
        deallocate(f_sum_old)
        deallocate(g_on_f_old)
        deallocate(dg_on_df_old)
        deallocate(dg_sum_old)
        deallocate(df_sum_old)
        deallocate(Psi_sum_old)
        deallocate(eigval_old)
        if(ii_xi.gt.ixi_stepslower+2)then
          deallocate(g_sum_old1)
          deallocate(f_sum_old1)
          deallocate(g_on_f_old1)
          deallocate(dg_on_df_old1)
          deallocate(dg_sum_old1)
          deallocate(df_sum_old1)
          deallocate(Psi_sum_old1)
          deallocate(eigval_old1)
        endif
      endif

      open(33,file='WF/Psi_sum_'//typ//'_1.dat')
      open(34,file='WF/Psi_sum_'//typ//'_2.dat')
      open(35,file='WF/Psi_sum_'//typ//'_3.dat')
      do i=1,nstates
        write(34,*)eigval(i)
      enddo
      do i=1,nstates
        do n=1,2*nm
          write(33,*)g_sum(i,n),f_sum(i,n),dg_sum(i,n)
          write(34,*)df_sum(i,n),Psi_sum(i,n)
          write(35,*)g_on_f(i,n),dg_on_df(i,n)
        enddo
      enddo
      close(33)
      close(34)
      close(35)
      deallocate(g_sum)
      deallocate(f_sum)
      deallocate(g_on_f)
      deallocate(dg_on_df)
      deallocate(dg_sum)
      deallocate(df_sum)
      deallocate(Psi_sum)
      deallocate(eigval_sort)
      return
      end

      subroutine swapping_of_states_abandoned_2(up_energy,
     &wave_new,nstates,nm,nkap,ii_xi,
     &ixi_stepslower,rmin,rmax,eigval,typ)

      include 'inc.par'
      real*8 Psi_sum(nstates),Psi_sum1(nstates),Psi_sum_old(nstates),
     &wave_new_sort(nstates,2*nm,-nkap:nkap),
     &wave_new(nstates,2*nm,-nkap:nkap),g_sum(nstates),g_sum1(nstates),
     &g_sum_old(nstates),f_sum(nstates),f_sum1(nstates),
     &f_sum_old(nstates),dg_sum(nstates),dg_sum1(nstates),
     &dg_sum_old(nstates),df_sum(nstates),df_sum1(nstates),
     &df_sum_old(nstates),ro(ns),ro1(ns),ro2(ns),
     &abstand(nstates),
     &sampl_point(nm,nstates),eigval(nstates),
     &eigval_old(nstates),eigval_sort(nstates),
     &energy_differences(nstates,nstates)
      integer nearest_neighbours(nstates,2),
     &nearest_neighbours_stor(nstates),flipped_old(nstates),
     &flipped(nstates)
      real*8,dimension(:),allocatable::t
      character*1 typ
      logical dkb,en_min(4),skipq(nstates)!,mkdirs
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      common /Barycentres/ RadiusOne,RadiusTwo
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      !mkdirs=makedirqq('WF')
      CALL SYSTEM("mkdir -p WF")

      a_IND=RadiusOne+RadiusTwo
      abstand=rmax*0.5d0/(z_nuc1+z_nuc2)**0.25d0
      sampl_point=0.d0
      if(ii_xi.gt.ixi_stepslower)then
        open(33,file='WF/Psi_sum_'//typ//'.dat',status='old')
        open(66,file='WF/Backflipped_'//typ//'.dat',status='old')
        do i=1,nstates
          read(33,'(6E13.6)')g_sum_old(i),f_sum_old(i),dg_sum_old(i),
     &    df_sum_old(i),Psi_sum_old(i),eigval_old(i)
          read(66,*)flipped_old(i)
        enddo
        close(33)
        close(66)
      else
        nearest_neighbours=0
      endif
      eigval_sort=eigval
      wave_new_sort=wave_new

      if(ii_xi.gt.ixi_stepslower)then
        do i=1,nstates
          do j=1,nstates
            energy_differences(j,i)=dabs(eigval(i)-eigval_old(j))
          enddo
        enddo
        nearest_neighbours_stor=minloc(energy_differences,dim=1)
        do i=1,nstates
          nearest_neighbours(i,1)=nearest_neighbours_stor(i)
        enddo
        do i=1,nstates
          j=nearest_neighbours_stor(i)
          energy_differences(j,i)=up_energy
        enddo
        nearest_neighbours_stor=minloc(energy_differences,dim=1,
     &  mask=energy_differences.ne.up_energy)
        do i=1,nstates
          nearest_neighbours(i,2)=nearest_neighbours_stor(i)
        enddo
      endif

      nu=nm+ns+2
      allocate(t(nu))

      call spline_arranger(nu,rmin,rmax,t)

      do jstate=1,nstates
        ii=1
        do i=1,nm
          sampl_point(ii,jstate)=(abstand(jstate)/
     &    dabs(abstand(jstate)-1.d0))*
     &    dabs(abstand(jstate)**(dble(ii)/dble(nm))-1.d0)
          ii=ii+1
        enddo
      enddo

      Kapp_range=nint(dabs(amu)+0.5d0)

      Psi_sum=0.d0
      g_sum=0.d0
      f_sum=0.d0
      dg_sum=0.d0
      df_sum=0.d0
      do i=1,nstates
        sampl_sum_g=0.d0
        sampl_sum_f=0.d0
c      pause
        do ji=1,ii-1
          if(ji.ne.0)then
            Sample_point=sampl_point(abs(ji),i)
            i2=1
            do while(t(i2).lt. Sample_point)
              i2=i2+1
            enddo
            if (t(i2).gt. Sample_point) then
              i2=i2-1
            endif
c         write(*,*)Sample_point,sampl_point(abs(ji),i),i2,i
c         pause
            call ddsplines(Sample_point,ro,ro1,ro2,i2,t,nu)

            wfgtotal=0.d0
            wfdgtotal=0.d0
            wfgtotal1=0.d0
            wfftotal=0.d0
            wfdftotal=0.d0
            wfftotal1=0.d0
            wfgtotal_tst=0.d0
            wfftotal_tst=0.d0
            wfgtotal1_tst=0.d0
            wfftotal1_tst=0.d0

            do kk=-nkap,nkap
              if(abs(kk).ge. abs(Kapp_range))then

                WFVal_G=0.d0
                WFVal_F=0.d0
                WFVal_DG=0.d0
                WFVal_DF=0.d0

                call DiracAngularL(-1.d0*kk,al1l)
                aj1j=dabs(1.d0*kk)-0.5d0
                ffakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
                ffakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
                call DiracAngularL(1.d0*kk,al1l)
                gfakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
                gfakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
                call DiracAngularL_int(kk,llk)
                ang_test=coeff_maker(amu,-1.d0*kk,1.d0*kk,llk)
                write(*,*)'spline_arranger: ', kk,llk,ang_test
                do i_spl=max(2,i2-ns+1),min(i2,nm+1)
                  if(dkb)then

                    WFVal_G=WFVal_G+
     &              wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1+nm,kk)*0.5d0*
     &              (ro1(i_spl-i2+ns)-kk*ro(i_spl-i2+ns)/Sample_point)
                    WFVal_DG=WFVal_DG+
     &              wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1+nm,kk)*(ro2(i_spl-i2+ns)-
     &              kk*ro1(i_spl-i2+ns)/Sample_point+
     &              kk*ro(i_spl-i2+ns)/Sample_point**2)
     &              /2.d0
                    WFVal_F=WFVal_F+
     &              wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1,kk)*0.5d0*
     &              (ro1(i_spl-i2+ns)+kk*ro(i_spl-i2+ns)/Sample_point)
                    WFVal_DF=WFVal_DF+
     &              wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1,kk)*(ro2(i_spl-i2+ns)+
     &              kk*ro1(i_spl-i2+ns)/Sample_point-
     &              kk*ro(i_spl-i2+ns)/Sample_point**2)
     &              /2.d0

                  else

                    WFVal_G=WFVal_G+
     &              wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)
                    WFVal_DG=WFVal_DG+
     &              wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)
                    WFVal_F=WFVal_F+
     &              wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)
                    WFVal_DF=WFVal_DF+
     &              wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)

                  endif
                enddo
                wfgtotal=wfgtotal+WFVal_G*gfakt1
                wfdgtotal=wfdgtotal+WFVal_DG*gfakt1
                wfgtotal1=wfgtotal1+WFVal_G*gfakt2
                wfftotal=wfftotal+WFVal_F*ffakt1
                wfdftotal=wfdftotal+WFVal_DF*ffakt1
                wfftotal1=wfftotal1+WFVal_F*ffakt2

                wfgtotal_tst=wfgtotal_tst+WFVal_G*ang_test*gfakt1
                wfftotal_tst=wfftotal_tst+WFVal_F*ang_test*gfakt1
                wfgtotal1_tst=wfgtotal1_tst+WFVal_G*ang_test*gfakt2
                wfftotal1_tst=wfftotal1_tst+WFVal_F*ang_test*gfakt2
              endif
            enddo

            Psi_sum(i)=Psi_sum(i)+
c     &      wfgtotal**2+wfftotal**2+wfgtotal1**2+wfftotal1**2
     &      wfgtotal_tst**2+wfftotal_tst**2+wfgtotal1_tst**2+
     &      wfftotal1_tst**2
            g_sum(i)=g_sum(i)+wfgtotal
            f_sum(i)=f_sum(i)+wfftotal
            dg_sum(i)=dg_sum(i)+wfdgtotal
            df_sum(i)=df_sum(i)+wfdftotal

          endif
        enddo
        g_sum(i)=dabs(g_sum(i))
        f_sum(i)=dabs(f_sum(i))
        dg_sum(i)=dabs(dg_sum(i))
        df_sum(i)=dabs(df_sum(i))
      enddo

      ii=0
      flipped=0
      skipq=.true.
      if(ii_xi.gt.ixi_stepslower)then
        do while(ii.le.nstates-1)
          ii=ii+1
          iii1=nearest_neighbours(ii,1)
          iii2=nearest_neighbours(ii,2)

          en_min=.false.
          if(ii+1.le.nstates)then
            if(dabs(eigval(ii)-eigval(iii1)).le.
     &      dabs(eigval(ii)-eigval(ii+1)))then
              en_min(1)=.true.
            endif
          else
            en_min(1)=.true.
          endif

          if(ii-1.ge.1)then
            if(dabs(eigval(ii)-eigval(iii1)).le.
     &      dabs(eigval(ii)-eigval(ii-1)))then
              en_min(2)=.true.
            endif
          else
            en_min(2)=.true.
          endif

          if(ii+1.le.nstates)then
            if(dabs(eigval(ii)-eigval(iii2)).le.
     &      dabs(eigval(ii)-eigval(ii+1)))then
              en_min(3)=.true.
            endif
          else
            en_min(3)=.true.
          endif

          if(ii-1.ge.1)then
            if(dabs(eigval(ii)-eigval(iii2)).le.
     &      dabs(eigval(ii)-eigval(ii-1)))then
              en_min(4)=.true.
            endif
          else
            en_min(4)=.true.
          endif

          if(flipped_old(ii).eq.iii1)then
            en_min(1)=.true.
            en_min(2)=.true.
          elseif(flipped_old(ii).eq.iii2)then
            en_min(3)=.true.
            en_min(4)=.true.
          endif

          if(en_min(1) .and. en_min(2).and. skipq(ii).and.
     &    dabs(Psi_sum(ii)-Psi_sum_old(iii1)).lt.
     &    dabs(Psi_sum(ii)-Psi_sum_old(iii2)).and.
     &    dabs(g_sum(ii)-g_sum_old(iii1)).lt.
     &    dabs(g_sum(ii)-g_sum_old(iii2)) .and.
     &    dabs(f_sum(ii)-f_sum_old(iii1)).lt.
     &    dabs(f_sum(ii)-f_sum_old(iii2)).and.
     &    dabs(dg_sum(ii)-dg_sum_old(iii1)).lt.
     &    dabs(dg_sum(ii)-dg_sum_old(iii2)) .and.
     &    dabs(df_sum(ii)-df_sum_old(iii1)).lt.
     &    dabs(df_sum(ii)-df_sum_old(iii2))
     &    )then

            if(flipped_old(iii1).ne.0)then
              if(
     &        dabs(Psi_sum(iii1)-Psi_sum_old(flipped_old(iii1))).lt.
     &        dabs(Psi_sum(iii1)-Psi_sum_old(iii1)).and.
     &        dabs(g_sum(iii1)-g_sum_old(flipped_old(iii1))).lt.
     &        dabs(g_sum(iii1)-g_sum_old(iii1)) .and.
     &        dabs(f_sum(iii1)-f_sum_old(flipped_old(iii1))).lt.
     &        dabs(f_sum(iii1)-f_sum_old(iii1)).and.
     &        dabs(dg_sum(iii1)-dg_sum_old(flipped_old(iii1))).lt.
     &        dabs(dg_sum(iii1)-dg_sum_old(iii1)) .and.
     &        dabs(df_sum(iii1)-df_sum_old(flipped_old(iii1))).lt.
     &        dabs(df_sum(iii1)-df_sum_old(iii1))
     &        )then
                flipped(iii1)=flipped_old(iii1)
                skipq(iii1)=.false.
                eigval_sort(iii1)=eigval(flipped_old(iii1))
                do k=-nkap,nkap
                  if(k.ne.0)then
                    do j=1,2*nm
                      wave_new_sort(iii1,j,k)=
     &                wave_new(flipped_old(iii1),j,k)
                    enddo
                  endif
                enddo
                g_sum1(iii1)=g_sum(flipped_old(iii1))
                f_sum1(iii1)=f_sum(flipped_old(iii1))
                dg_sum1(iii1)=dg_sum(flipped_old(iii1))
                df_sum1(iii1)=df_sum(flipped_old(iii1))
                Psi_sum1(iii1)=Psi_sum(flipped_old(iii1))
              endif
            endif

            flipped(ii)=iii1
            eigval_sort(ii)=eigval(iii1)
            do k=-nkap,nkap
              if(k.ne.0)then
                do j=1,2*nm
                  wave_new_sort(ii,j,k)=wave_new(iii1,j,k)
                enddo
              endif
            enddo
            g_sum1(ii)=g_sum(iii1)
            f_sum1(ii)=f_sum(iii1)
            dg_sum1(ii)=dg_sum(iii1)
            df_sum1(ii)=df_sum(iii1)
            Psi_sum1(ii)=Psi_sum(iii1)

          elseif(en_min(3) .and. en_min(4).and. skipq(ii).and.
     &    dabs(Psi_sum(ii)-Psi_sum_old(iii2)).lt.
     &    dabs(Psi_sum(ii)-Psi_sum_old(iii1)).and.
     &    dabs(g_sum(ii)-g_sum_old(iii2)).lt.
     &    dabs(g_sum(ii)-g_sum_old(iii1)) .and.
     &    dabs(f_sum(ii)-f_sum_old(iii2)).lt.
     &    dabs(f_sum(ii)-f_sum_old(iii1)).and.
     &    dabs(dg_sum(ii)-dg_sum_old(iii2)).lt.
     &    dabs(dg_sum(ii)-dg_sum_old(iii1)) .and.
     &    dabs(df_sum(ii)-df_sum_old(iii2)).lt.
     &    dabs(df_sum(ii)-df_sum_old(iii1))
     &    )then

            if(flipped_old(iii2).ne.0)then
              if(
     &        dabs(Psi_sum(iii2)-Psi_sum_old(flipped_old(iii2))).lt.
     &        dabs(Psi_sum(iii2)-Psi_sum_old(iii2)).and.
     &        dabs(g_sum(iii2)-g_sum_old(flipped_old(iii2))).lt.
     &        dabs(g_sum(iii2)-g_sum_old(iii2)) .and.
     &        dabs(f_sum(iii2)-f_sum_old(flipped_old(iii2))).lt.
     &        dabs(f_sum(iii2)-f_sum_old(iii2)).and.
     &        dabs(dg_sum(iii2)-dg_sum_old(flipped_old(iii2))).lt.
     &        dabs(dg_sum(iii2)-dg_sum_old(iii2)) .and.
     &        dabs(df_sum(iii2)-df_sum_old(flipped_old(iii2))).lt.
     &        dabs(df_sum(iii2)-df_sum_old(iii2))
     &        )then
                flipped(iii2)=flipped_old(iii2)
                skipq(iii2)=.false.
                eigval_sort(iii2)=eigval(flipped_old(iii2))
                do k=-nkap,nkap
                  if(k.ne.0)then
                    do j=1,2*nm
                      wave_new_sort(iii2,j,k)=
     &                wave_new(flipped_old(iii2),j,k)
                    enddo
                  endif
                enddo
                g_sum1(iii2)=g_sum(flipped_old(iii2))
                f_sum1(iii2)=f_sum(flipped_old(iii2))
                dg_sum1(iii2)=dg_sum(flipped_old(iii2))
                df_sum1(iii2)=df_sum(flipped_old(iii2))
                Psi_sum1(iii2)=Psi_sum(flipped_old(iii2))
              endif
            endif

            flipped(ii)=iii2
            eigval_sort(ii)=eigval(iii2)
            do k=-nkap,nkap
              if(k.ne.0)then
                do j=1,2*nm
                  wave_new_sort(ii,j,k)=wave_new(iii2,j,k)
                enddo
              endif
            enddo
            g_sum1(ii)=g_sum(iii2)
            f_sum1(ii)=f_sum(iii2)
            dg_sum1(ii)=dg_sum(iii2)
            df_sum1(ii)=df_sum(iii2)
            Psi_sum1(ii)=Psi_sum(iii2)

          endif
        enddo
        eigval=eigval_sort
        wave_new=wave_new_sort
        g_sum=g_sum1
        f_sum=f_sum1
        dg_sum=dg_sum1
        df_sum=df_sum1
        Psi_sum=Psi_sum1
      endif
      open(33,file='WF/Psi_sum_'//typ//'.dat')
      open(66,file='WF/Backflipped_'//typ//'.dat')
      do i=1,nstates
        write(33,'(6E13.6)')g_sum(i),f_sum(i),dg_sum(i),df_sum(i),
     &  Psi_sum(i),eigval(i)
        write(66,*)flipped(i)
      enddo
      deallocate(t)
      return
      end

      subroutine swapping_of_states_abandoned(up_energy,
     &wave_new,nstates,nm,nkap,ii_xi,
     &ixi_stepslower,rmin,rmax,eigval,typ)

      include 'inc.par'
      real*8 Psi_sum(nstates),Psi_sum_old(nstates),
     &wave_new(nstates,2*nm,-nkap:nkap),g_sum(nstates),
     &g_sum_old(nstates),f_sum(nstates),f_sum_old(nstates),
     &dg_sum(nstates),dg_sum_old(nstates),
     &df_sum(nstates),df_sum_old(nstates),ro(ns),ro1(ns),ro2(ns),
     &abstand(nstates),coeff_storage(-nkap:nkap,2*nm),
     &sampl_point(nm,nstates),eigval(nstates),eigval_1(nstates),
     &eigval_sort(nstates)
      integer backflip(nstates),nearest_neighbours(nstates)
      real*8,dimension(:),allocatable::t
      character*1 typ
      logical dkb!,mkdirs
      common /common_dkb/ dkb
      common /momentum_projection/ amu,amj_max
      common /Barycentres/ RadiusOne,RadiusTwo
      common /nuc_charge/ z_nuc1,az1,z_nuc2,az2

      !mkdirs=makedirqq('WF')
      CALL SYSTEM("mkdir -p WF")

      a_IND=RadiusOne+RadiusTwo
      abstand=rmax*0.5d0/(z_nuc1+z_nuc2)**0.25d0
      sampl_point=0.d0
      if(ii_xi.gt.ixi_stepslower)then
        open(33,file='WF/Psi_sum_'//typ//'.dat',status='old')
        open(74,file='WF/Backflip'//typ//'.dat',status='old')
        do i=1,nstates
          read(33,'(5E16.9)')g_sum_old(i),f_sum_old(i),dg_sum_old(i),
     &    df_sum_old(i),Psi_sum_old(i)
          read(74,*)backflip(i)
        enddo
        close(33)
        close(74)
      else
        do i=1,nstates
          backflip(i)=i
        enddo
      endif
      eigval_sort=eigval

      do i=1,nstates
        if(backflip(i).ne. i)then
          eval_storage=eigval(i)
          eigval(i)=eigval(backflip(i))
          eigval(backflip(i))=eval_storage
          do k=-nkap,nkap
            if(k.ne.0)then
              do j=1,2*nm
                coeff_storage(k,j)=wave_new(i,j,k)
              enddo
            endif
          enddo
          do k=-nkap,nkap
            if(k.ne.0)then
              do j=1,2*nm
                wave_new(i,j,k)=wave_new(backflip(i),j,k)
              enddo
            endif
          enddo
          do k=-nkap,nkap
            if(k.ne.0)then
              do j=1,2*nm
                wave_new(backflip(i),j,k)=coeff_storage(k,j)
              enddo
            endif
          enddo
        endif
      enddo
      eigval_1=eigval
      do i=1,nstates
        nearest_neighbours(i)=minloc(eigval_1,dim=1,mask=eigval_1.ne.
     &  up_energy)
        eigval_1(nearest_neighbours(i))=up_energy
      enddo
      eigval_1=eigval

      nu=nm+ns+2
      allocate(t(nu))

      call spline_arranger(nu,rmin,rmax,t)

      do jstate=1,nstates
        ii=1
        do i=1,nm
          sampl_point(ii,jstate)=(abstand(jstate)/
     &    dabs(abstand(jstate)-1.d0))*
     &    dabs(abstand(jstate)**(dble(ii)/dble(nm))-1.d0)
          ii=ii+1
        enddo
      enddo

      Kapp_range=nint(dabs(amu)+0.5d0)

      Psi_sum=0.d0
      g_sum=0.d0
      f_sum=0.d0
      dg_sum=0.d0
      df_sum=0.d0
      do i=1,nstates
        sampl_sum_g=0.d0
        sampl_sum_f=0.d0
c      pause
        do ji=1,ii-1
          if(ji.ne.0)then
            Sample_point=sampl_point(abs(ji),i)
            i2=1
            do while(t(i2).lt. Sample_point)
              i2=i2+1
            enddo
            if (t(i2).gt. Sample_point) then
              i2=i2-1
            endif
c         write(*,*)Sample_point,sampl_point(abs(ji),i),i2,i
c         pause
            call ddsplines(Sample_point,ro,ro1,ro2,i2,t,nu)

            wfgtotal=0.d0
            wfdgtotal=0.d0
            wfgtotal1=0.d0
            wfftotal=0.d0
            wfdftotal=0.d0
            wfftotal1=0.d0
            wfgtotal_tst=0.d0
            wfftotal_tst=0.d0
            wfgtotal1_tst=0.d0
            wfftotal1_tst=0.d0

            do kk=-nkap,nkap
              if(abs(kk).ge. abs(Kapp_range))then

                WFVal_G=0.d0
                WFVal_F=0.d0
                WFVal_DG=0.d0
                WFVal_DF=0.d0

                call DiracAngularL(-1.d0*kk,al1l)
                aj1j=dabs(1.d0*kk)-0.5d0
                ffakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
                ffakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
                call DiracAngularL(1.d0*kk,al1l)
                gfakt1=clebsh(al1l,amu-5.d-1,5.d-1,5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu-5.d-1),-1.d0)
                gfakt2=clebsh(al1l,amu+5.d-1,5.d-1,-5.d-1,aj1j,amu)*
     &          spharm(idint(al1l),idint(amu+5.d-1),-1.d0)
                ang_test=coeff_maker(amu,-1.d0*kk,1.d0*kk,1)
                do i_spl=max(2,i2-ns+1),min(i2,nm+1)
                  if(dkb)then

                    WFVal_G=WFVal_G+
     &              wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1+nm,kk)*0.5d0*
     &              (ro1(i_spl-i2+ns)-kk*ro(i_spl-i2+ns)/Sample_point)
                    WFVal_DG=WFVal_DG+
     &              wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1+nm,kk)*(ro2(i_spl-i2+ns)-
     &              kk*ro1(i_spl-i2+ns)/Sample_point+
     &              kk*ro(i_spl-i2+ns)/Sample_point**2)
     &              /2.d0
                    WFVal_F=WFVal_F+
     &              wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1,kk)*0.5d0*
     &              (ro1(i_spl-i2+ns)+kk*ro(i_spl-i2+ns)/Sample_point)
                    WFVal_DF=WFVal_DF+
     &              wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)+
     &              wave_new(i,i_spl-1,kk)*(ro2(i_spl-i2+ns)+
     &              kk*ro1(i_spl-i2+ns)/Sample_point-
     &              kk*ro(i_spl-i2+ns)/Sample_point**2)
     &              /2.d0

                  else

                    WFVal_G=WFVal_G+
     &              wave_new(i,i_spl-1,kk)*ro(i_spl-i2+ns)
                    WFVal_DG=WFVal_DG+
     &              wave_new(i,i_spl-1,kk)*ro1(i_spl-i2+ns)
                    WFVal_F=WFVal_F+
     &              wave_new(i,i_spl-1+nm,kk)*ro(i_spl-i2+ns)
                    WFVal_DF=WFVal_DF+
     &              wave_new(i,i_spl-1+nm,kk)*ro1(i_spl-i2+ns)

                  endif
                enddo
                wfgtotal=wfgtotal+WFVal_G*gfakt1
                wfdgtotal=wfdgtotal+WFVal_DG*gfakt1
                wfgtotal1=wfgtotal1+WFVal_G*gfakt2
                wfftotal=wfftotal+WFVal_F*ffakt1
                wfdftotal=wfdftotal+WFVal_DF*ffakt1
                wfftotal1=wfftotal1+WFVal_F*ffakt2

                wfgtotal_tst=wfgtotal_tst+WFVal_G*ang_test*gfakt1
                wfftotal_tst=wfftotal_tst+WFVal_F*ang_test*gfakt1
                wfgtotal1_tst=wfgtotal1_tst+WFVal_G*ang_test*gfakt2
                wfftotal1_tst=wfftotal1_tst+WFVal_F*ang_test*gfakt2
              endif
            enddo

            Psi_sum(i)=Psi_sum(i)+
c     &      wfgtotal**2+wfftotal**2+wfgtotal1**2+wfftotal1**2
     &      wfgtotal_tst**2+wfftotal_tst**2+wfgtotal1_tst**2+
     &      wfftotal1_tst**2
            g_sum(i)=g_sum(i)+wfgtotal
            f_sum(i)=f_sum(i)+wfftotal
            dg_sum(i)=dg_sum(i)+wfdgtotal
            df_sum(i)=df_sum(i)+wfdftotal

          endif
        enddo
        g_sum(i)=dabs(g_sum(i))
        f_sum(i)=dabs(f_sum(i))
        dg_sum(i)=dabs(dg_sum(i))
        df_sum(i)=dabs(df_sum(i))
      enddo

      do i=1,nstates
        jj=1
        jjj=nstates
        do while(jjj.gt.1 .and. jj.le.nstates)
          if(eigval_1(jj)-up_energy.eq. 0.d0)jjj=jjj-1
          jj=jj+1
        enddo
        ii=minloc(eigval_1,dim=1, mask=eigval_1.ne. up_energy)
        nn1=1
        do while(eigval_sort(nn1)-eigval_1(ii).ne. 0.d0)
          nn1=nn1+1
        enddo
        eigval_1(ii)=up_energy
        iii=minloc(eigval_1,dim=1, mask=eigval_1.ne. up_energy)
        nn2=1
        do while(eigval_sort(nn2)-eigval_1(iii).ne. 0.d0)
          nn2=nn2+1
        enddo
        if(jjj.gt.1 .and.
     &  dabs(Psi_sum(ii)-Psi_sum_old(iii)).lt.
     &  dabs(Psi_sum(ii)-Psi_sum_old(ii)).and.
     &  dabs(g_sum(ii)-g_sum_old(iii)).lt.
     &  dabs(g_sum(ii)-g_sum_old(ii)) .and.
     &  dabs(f_sum(ii)-f_sum_old(iii)).lt.
     &  dabs(f_sum(ii)-f_sum_old(ii)).and.
     &  dabs(dg_sum(ii)-dg_sum_old(iii)).lt.
     &  dabs(dg_sum(ii)-dg_sum_old(ii)) .and.
     &  dabs(df_sum(ii)-df_sum_old(iii)).lt.
     &  dabs(df_sum(ii)-df_sum_old(ii)).and.
     &  ii_xi.gt.ixi_stepslower) then

          if(backflip(iii).eq.iii)then
            backflip(iii)=ii
          else
            backflip(ii)=nearest_neighbours(ii+1)
          endif
          eigval_1(iii)=up_energy
          eval_storage=eigval(ii)
          eigval(ii)=eigval(iii)
          eigval(iii)=eval_storage
          do k=-nkap,nkap
            if(abs(k).ge. abs(Kapp_range))then
              do j=1,2*nm
                coeff_storage(k,j)=wave_new(ii,j,k)
              enddo
            endif
          enddo
          do k=-nkap,nkap
            if(abs(k).ge. abs(Kapp_range))then
              do j=1,2*nm
                wave_new(ii,j,k)=wave_new(iii,j,k)
              enddo
            endif
          enddo
          do k=-nkap,nkap
            if(abs(k).ge. abs(Kapp_range))then
              do j=1,2*nm
                wave_new(iii,j,k)=coeff_storage(k,j)
              enddo
            endif
          enddo
          Psi_storage=Psi_sum(ii)
          Psi_sum(ii)=Psi_sum(iii)
          Psi_sum(iii)=Psi_storage
          g_storage=g_sum(ii)
          g_sum(ii)=g_sum(iii)
          g_sum(iii)=g_storage
          f_storage=f_sum(ii)
          f_sum(ii)=f_sum(iii)
          f_sum(iii)=f_storage
          dg_storage=dg_sum(ii)
          dg_sum(ii)=dg_sum(iii)
          dg_sum(iii)=dg_storage
          df_storage=df_sum(ii)
          df_sum(ii)=df_sum(iii)
          df_sum(iii)=df_storage
c              Psi_storage=Psi_sum_old(ii)
c              Psi_sum_old(ii)=Psi_sum_old(iii)
c              Psi_sum_old(iii)=Psi_storage
c              g_storage=g_sum_old(ii)
c              g_sum_old(ii)=g_sum_old(iii)
c              g_sum_old(iii)=g_storage
c              f_storage=f_sum_old(ii)
c              f_sum_old(ii)=f_sum_old(iii)
c              f_sum_old(iii)=f_storage
c              dg_storage=dg_sum_old(ii)
c              dg_sum_old(ii)=dg_sum_old(iii)
c              dg_sum_old(iii)=dg_storage
c              df_storage=df_sum_old(ii)
c              df_sum_old(ii)=df_sum_old(iii)
c              df_sum_old(iii)=df_storage
        endif
      enddo
      open(33,file='WF/Psi_sum_'//typ//'.dat')
      open(74,file='WF/Backflip'//typ//'.dat')
      do i=1,nstates
        write(33,'(5E16.9)')g_sum(i),f_sum(i),dg_sum(i),df_sum(i),
     &  Psi_sum(i)
        write(74,*)backflip(i)
      enddo
      deallocate(t)
      return
      end

      subroutine mj_proj_accounting_neqz(nkap,nstates_mj,number_states,
     &numb_mj)
      include 'inc.par'
      integer number_states(2,-nkap:nkap),numb_mj(4,nkap)
      common /momentum_projection/ amu,amj_max

      mj_max=nint(amj_max+0.5d0)

      nstates_mj=0
      do mj=mj_max,1,-1
        do k=-nkap,-mj
          do i=number_states(1,k),number_states(2,k)
            nstates_mj=nstates_mj+1
          enddo
        enddo
        do k=mj,nkap
          do i=number_states(1,k),number_states(2,k)
            nstates_mj=nstates_mj+1
          enddo
        enddo
      enddo
      nstates_mj=2*nstates_mj

      do mj=mj_max,1,-1
        iz=0
        numb_mj(1,mj)=iz+1
        do kk=-nkap,-mj
          do j=number_states(1,kk),number_states(2,kk)
            iz=iz+1
          enddo
        enddo
        numb_mj(2,mj)=iz
        do kk=-mj+1,mj-1
          if(kk.ne.0)then
            do j=number_states(1,kk),number_states(2,kk)
              iz=iz+1
            enddo
          endif
        enddo
        numb_mj(3,mj)=iz+1
        do kk=mj,nkap
          do j=number_states(1,kk),number_states(2,kk)
            iz=iz+1
          enddo
        enddo
        numb_mj(4,mj)=iz
      enddo

      return
      end

      subroutine mj_proj_accounting_eqz(nkap,nste_mj,nsto_mj,
     &number_states,numbe_mj,numbo_mj)
      include 'inc.par'
      integer number_states(2,-nkap:nkap),numbe_mj(2,nkap),
     &numbo_mj(2,nkap)
      common /momentum_projection/ amu,amj_max

      mj_max=nint(amj_max+0.5d0)

      nste_mj=0
      do mj=1,mj_max
        kfactor=(-1)**(mj+1)
        do kk=-mj,-nkap,-1
          k=kk*kfactor
          do i=number_states(1,k),number_states(2,k)
            nste_mj=nste_mj+1
          enddo
          kfactor=-kfactor
        enddo
      enddo

      nsto_mj=0
      do mj=1,mj_max
        kfactor=(-1)**(mj+1)
        do kk=mj,nkap
          k=kk*kfactor
          do i=number_states(1,k),number_states(2,k)
            nsto_mj=nsto_mj+1
          enddo
          kfactor=-kfactor
        enddo
      enddo
      nste_mj=nste_mj*2
      nsto_mj=nsto_mj*2

      do mj=mj_max,1,-1
        ize=0
        izo=0
        numbe_mj(1,mj)=ize+1
        numbo_mj(1,mj)=izo+1
        if((-1)**nkap.gt. 0)then
          kfactor=1
          do k=-nkap,-mj
            kk=k*kfactor
            do j=number_states(1,kk),number_states(2,kk)
              izo=izo+1
            enddo
            do j=number_states(1,-kk),number_states(2,-kk)
              ize=ize+1
            enddo
            kfactor=-kfactor
          enddo
          numbo_mj(2,mj)=izo
          numbe_mj(2,mj)=ize
        else
          kfactor=1
          do k=-nkap,-mj
            kk=k*kfactor
            do j=number_states(1,kk),number_states(2,kk)
              ize=ize+1
            enddo
            do j=number_states(1,-kk),number_states(2,-kk)
              izo=izo+1
            enddo
            kfactor=-kfactor
          enddo
          numbe_mj(2,mj)=ize
          numbo_mj(2,mj)=izo
        endif
      enddo

      return
      end

      subroutine redefineEigvalWaveNew(n_jstates,eigval,eigval_mj,
     &wave_new,wave_new_mj,nstates,nm,nkap)
      include 'inc.par'
      real*8 eigval(2*n_jstates*nstates),eigval_mj(nstates,-n_jstates:
     &n_jstates),
     &wave_new(2*n_jstates*nstates,2*nm,-nkap:nkap),
     &wave_new_mj(nstates,2*nm,-nkap:nkap,-n_jstates:n_jstates)

      i_count = 0
      do n_jstate=-n_jstates,n_jstates
        if (n_jstate .ne. 0)then
          eigval((i_count*nstates+1):(i_count+1)*nstates) =
     &    eigval_mj(:,n_jstate)
          wave_new((i_count*nstates+1):(i_count+1)*nstates,:,
     &    -nkap:nkap) = wave_new_mj(:,:,:,n_jstate)
          i_count = i_count+1
        endif
      enddo
      nstates = 2*n_jstates*nstates
      return
      end

      subroutine getIntegralWeights(ttt, www)
      include 'inc.par'
      real*8 ttt(nuz), www(nuz)
      common /weights/ w4n(4),t4n(4),w8(8),t8(8),w16(16),t16(16),
     &w32(32),t32(32),w64(64),t64(64),t6(6),w6(6)

      select case (nuz)
        case(4)
          do i=1,4
            ttt(i)=t4n(i)
            www(i)=w4n(i)
          enddo
        case(8)
          do i=1,8
            ttt(i)=t8(i)
            www(i)=w8(i)
          enddo
        case(6)
          do i=1,6
            ttt(i)=t6(i)
            www(i)=w6(i)
          enddo
        case(16)
          do i=1,16
            ttt(i)=t16(i)
            www(i)=w16(i)
          enddo
        case(32)
          do i=1,32
            ttt(i)=t32(i)
            www(i)=w32(i)
          enddo
        case(64)
          do i=1,64
            ttt(i)=t64(i)
            www(i)=w64(i)
          enddo
      end select
      return
      end
