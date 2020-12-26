      subroutine io(n,de)
      implicit real*8(a-h,o-z)
      character*60 podp(50)
      real*8 anum(50)
      
      podp(1)='zero order energy'
      podp(2)='1st order additional potential'
      podp(3)='2nd order additional potential non_qed'
      podp(4)='SES irred. 1s shell'
      podp(5)='SES irred. 2s shell'
      podp(6)='SES red. 1s shell'
      podp(7)='SES red. 2s shell'
      podp(8)='SES VR_0 1s shell'
      podp(9)='SES VR_0 2s shell'
      podp(10)='SES VR_many 1s shell'
      podp(11)='SES VR_many 2s shell'
      podp(12)='SE'
      podp(13)='SE + add. pot.'
      podp(14)='SE + add. pot. VR_0'
      podp(15)='SE + add. pot. VR_many'
      podp(16)='1ph exch. 1s shell'
      podp(17)='1ph exch. 2s shell'
      podp(18)='2ph exch. ladd. dir 1s shell'
      podp(19)='2ph exch. ladd. exch 1s shell'
      podp(20)='2ph exch. ladd. dir 2s shell'
      podp(21)='2ph exch. ladd. exch 2s shell'
      podp(22)='2ph exch. cros. dir 1s shell'
      podp(23)='2ph exch. cros. exch 1s shell'
      podp(24)='2ph exch. cros. dir 2s shell'
      podp(25)='2ph exch. cros. exch 2s shell'
      podp(26)='2ph exch. red. 1s shell'
      podp(27)='2ph exch. red. 2s shell'
      podp(28)='2ph exch. 3 el.'
      podp(29)='1ph exhc. + add. potl. 1s shell'
      podp(30)='1ph exhc. + add. potl. 2s shell'
      podp(31)='VPS 1s shell'
      podp(32)='VPS 2s shell'
      podp(33)='VPS loop 1s shell'
      podp(34)='VPS loop 2s shell'
      podp(35)='VP'
      podp(36)='VP + add. potl. line'
      podp(37)='VP + add. potl. loop'
      podp(44)='2nd order additional potential qed'
      
      
      open(1,file='total.dat')
      read(1,*) anum
      close(1)

      anum(n)=de

      open(1,file='total.dat')
      write(1,*) anum
      close(1)

      su=0.d0
      open(1,file='total.read')
      do i=1,50
         if(anum(i).ne.0.d0)then
            write(1,*) podp(i)
            write(1,*) anum(i)
            su=su+anum(i)
         endif
      enddo

      write(1,*)'**********************'
      write(1,*) su,'Total'
      close(1)
      return
      end

