      program iterator
! Iterate brtest on BR read from file
! Extrapolate from 118 and 121 in logspace
! + impose random condition on b_CH3  at 140 nm
! + separate CH from others
      implicit real*8 (a-h,o-z)
      character*14 chain
      character*4  count
      dimension b_1(4), b_2(4), b_3(4), b(4)
      dimension ab(4), bb(4), cb(4)
      dimension b_1_CH(2), b_2_CH(2), db_CH(2), b_CH(2)
      real*8 l1, l2, l3
      common /bran/ iseed

      iseed=-12547

      nrun=60
      nskip1=036 ! Lines to skip in sample 1
      nskip2=128 ! id
      eps=1e-6   ! Thresholding factor

      print *,'Monte Carlo Uncertainty Propagation'
      print *,'-----------------------------------'
      print *
      print *,'Nombre de runs = ',nrun
      print *,"Facteurs d'incertitude lus dans INPUT/uncert.dat"
      print *

      l1=118.2
      l2=121.6
      l3=140.0

      dl_CH=123.6-105.8
      
      open(10,file='rb-118.dat')
      open(11,file='rb-121.dat')
      read(10,*) ! Order in files : CH3,CH2a,CH2X,CH 
      do iskip=1,nskip1
         read(10,*,iostat=ios) 
      enddo
      read(11,*) ! Skip header
      do iskip=1,nskip2
         read(11,*,iostat=ios) 
      enddo
      if=20
      do i= 100, 140 , 1
         write(count,'(i4.4)') i
         open(if+i-100,file='wlBR/br_CH4'//count//'.dat')
      enddo
      do irun= 1, nrun         
         read(10,*,iostat=ios) b_1
         if(ios.ne.0) exit
         read(11,*,iostat=ios) b_2
         if(ios.ne.0) exit

         b_3(1)=unidev_01()
         b_3(2)=1-b_3(1)
         b_3(3)=0

         ! Thresholding to avoid numerical problems
         b_1=min(max(b_1,eps),1-eps)
         b_2=min(max(b_2,eps),1-eps)
         b_3=min(max(b_3,eps),1-eps)

         ! Treat CH at the first level
         b_1_CH(1)=0.23+gasdev(3d-2) !b_1(4)
         b_2_CH(1)=0.06+gasdev(5d-3) !b_2(4)
         b_1_CH(2)=1-b_1_CH(1)
         b_2_CH(2)=1-b_2_CH(1)
         ! Centered logratio transform
         b_1_CH=log(b_1_CH/product(b_1_CH)**0.5d0)
         b_2_CH=log(b_2_CH/product(b_2_CH)**0.5d0)

         ! Renormalize remaining components
         b_1=b_1/sum(b_1(1:3))
         b_2=b_2/sum(b_2(1:3))
         b_3=b_3/sum(b_3(1:3))

         ! Centered logratio transform
         b_1=log(b_1/product(b_1(1:3))**(1/3.d0))
         b_2=log(b_2/product(b_2(1:3))**(1/3.d0))
         b_3=log(b_3/product(b_3(1:3))**(1/3.d0))
        
         ! Slope for Linear interpolation of CH
         db_CH =(b_2_CH-b_1_CH)/dl_CH

         ! Solve quadratic equation for other components
         del=l1**2*(l2-l3)
     >      -l2**2*(l1-l3)
     >      +l3**2*(l1-l2)
         ab=(b_1*(l2-l3)
     >      -b_2*(l1-l3)
     >      +b_3*(l1-l2))/del
         bb=(l1**2*(b_2-b_3)
     >      -l2**2*(b_1-b_3)
     >      +l3**2*(b_1-b_2))/del
         cb=(l1**2*(l2*b_3-l3*b_2)
     >      -l2**2*(l1*b_3-l3*b_1)
     >      +l3**2*(l1*b_2-l2*b_1))/del

         if=20
         do i= 100, 140 , 1

            ! Linear inter/extra-polation in Cl space
            b_CH = b_1_CH  + (i-105.8)*db_CH
            ! Inverse CLT
            b_CH= exp(b_CH)/sum(exp(b_CH))

            ! Quadratic ...
            b = max(min(ab*i**2  + bb*i + cb,1.d2),-1e2)
            ! Inverse CLT and normalize to 1-CH
            b= b_CH(2)*exp(b)/sum(exp(b(1:3)))

            ! Output and congregate CH2 states
            write(if+i-100,'(4e15.3e3)') dble(i),b(2)+b(3),b_CH(1),b(1)
         enddo
 
      enddo
      close(10)
      close(11)
        if=20
         do i= 100, 140 , 1
            close(if+i-100)
         enddo

      end
