      programme q
      implicit double precision (a-h,o-z)
c-----
      parameter (nd=4000)
      character*30 filename
c      character*3 gate,mass
c      character*7 dope
      common /var/ f(nd),v(nd),s(nd),cn(nd),s0(nd),
     *             ep(nd),ca(nd),v0(nd), cband(nd), af(nd)
      dimension a(nd),b(nd),c(nd),d(nd),r(nd)
      parameter (qel=1.60219e-19,   !electronic charge
     *           bkb=1.38062e-23,   !Boltamann constant
     *           am0=0.910956e-30,  !electron mass
     *           hbr=1.05459e-32,   !h bar
     *           e00=8.85419e-14,   !permitivity of free space
     *           esr=11.9,    	    !rel. permitivity of Si
     *           eor=3.9,    	    !rel. permitivity of oxide
     *           afs=4.05,    	    !electron affinity of Si
     *           afo=0.95,    	    !electron affinity of oxide
     *           tab=300,     	    !temperature
     *           ft=bkb*tab/qel,    !kT/q
     *           cni=1.45e10) 	    !intrinsic concentration
         
      parameter (vzero=1e-6)

c-----Read input file
      open (unit=20,file='1D.inp',status='old')
      read (20,*) icl_or_dg
      read (20,*) Vgate
      read (20,*) cna
      read (20,*) tox
      read (20,*) amr
      close (20)
c-----
      amrox=0.22d0	  !relative effective mass
      aml=am0*amr	 !absolute effective mass
      amox=am0*amrox	 !absolute effective mass
      
      bn=((hbr/aml)*hbr)/(12d0*qel)
      box=((hbr/amox)*hbr)/(12d0*qel)
      bn1=1.0/bn
      
      xp=hbr/sqrt(2d0*amox*(afs-afo)*qel)     !char. penetration depth
      print*,'xp = ',xp
c-----
      niter=5000   	      !max number of iterations
      niterp=10   	      !max number of iterations
      niterdg=5   	      !max number of iterations
      
      epsr=1d-4
      epsr_big=1d-8
c-----
      open (unit=10,file='d.d')
!      vg0=2.0	  !1.4282
      do vg=Vgate,Vgate,0.1   !set Vg
c      do vg=1.5d0,1.5d0	!,0.5d0   !set Vg
      dx=1e-8	  	      !grid spacing
c      tox=1.0001e-7	      !oxide thickness
      tox=tox*1d-7	      !oxide thickness
      nox=tox/dx+1	      !grid line for interface
      print*,nox
      n=800+nox 	  	      !number of grid nodes in device
!      cna=5e17	  	      !substrate doping conc
      fn=ft*log(cna/cni)      !fermi level
      
      grad_at_interface = -(box/bn)*dx/xp - 1d0
      
      print*,'Interface gradient = ',grad_at_interface
      
      filename='cn_dg.dat'
 
      do i=1,n
       if (i.lt.nox) then
        ep(i)=e00*eor	      !permitivity in oxide
        af(i)=afo	      !electron affinity in oxide
       else
        ep(i)=e00*esr	      !permitivity in Si
        af(i)=afs 	      !electron affinity in Si
        ca(i)=cna 	      !doping in Si
       end if
       v(i)=0.0   	      !initial potential = 0
      end do
      v(1)=vg	  	      !potential on gate = Vg
c-----First Poisson solution
      do iter=1,niter
       do i=2,n-1
        if (i.lt.nox) then
         bb=0.0
         bj=0.0
        else
c        calculate RHS charge
         bb=qel*(-ca(i)+cni*(-exp((v(i)-fn)/ft)+exp((fn-v(i))/ft)))*dx
c        calculate Jacobian
         bj=-(qel*cni/ft)*(exp((fn-v(i))/ft)+exp((v(i)-fn)/ft))*dx
        end if
	
	if (i.eq.nox) then   ! Half the volume at the interface
	  bb = bb/2d0
	  bj = bj/2d0
        end if
	
c       calculate discretisation	
        a(i-1)=ep(i-1)/dx
        b(i-1)=-((ep(i-1)+ep(i))/dx)+bj
        c(i-1)=ep(i)/dx
        d(i-1)=(-bb+bj*v(i))
	
       end do
c-----
c      fix boundary conditions
       d(1)=d(1)-a(1)*v(1)
       a(1)=0.0
       d(n-2)=d(n-2)-c(n-2)*v(n)
       c(n-2)=0.0
c-----
c      solve with Gaussian elimination
       call gaus3d(a,b,c,d,r,n-2)
c-----
c      check error between this solution and previous value of v
       error=0.0
       ecount=0d0
       do i=2,n-1
c        if (abs(r(i-1)).ge.vzero) then
         error=error+abs((v(i)-r(i-1))/(vzero+r(i-1)))
	 ecount=ecount+1d0
c        end if
        v(i)=r(i-1)
       end do
       print*,error,ecount
       error=error/ecount
       print*,'P1==>',iter,error
       if (error.le.epsr) goto 100  !if met tolerance then break out
      end do
      print*,'No convergence'
  100 continue
      print*,'===>',v(16)
      do i=1,n
       v0(i)=v(i)
       if (i.ge.nox) then
         s0(i)=cni*exp((v(i)-fn)/ft)
         cn(i)=s0(i)
       end if
      end do
      
      if (icl_or_dg.eq.0) goto 400
      
c---- Start Density Gradient solution

      do ibig=1,niter
      do i=nox,n
       cn(i)=cni*exp((v(i)-fn)/ft)
       s(i)=sqrt(cn(i)/cni)
      end do
c      s(nox)=0.0
      
c-----
      do iter=1,niterdg
      
c          bb=-bn1*(s(nox)*(fn-v(nox))/2.d0+s(nox)*ft*log(s(nox)))
c c        calculate Jacobian
c          bj=-bn1*((fn-v(nox))/2.0+ft*log(s(nox))+ft)
c          
c c       disctretisation at nox
c         a(1)=0d0
c         b(1)=-1.0/(dx*dx)+bj
c         c(1)=1.0/(dx*dx)
c         d(1)=(-bb+bj*s(nox)) - (box/bn*dx/xp)*s(nox)

c       disctretisation at nox
        a(1)=0d0
        b(1)=grad_at_interface
        c(1)=1d0
        d(1)=0d0
      
        do i=nox+1,n-1
c        calculate RHS 
         bb=-bn1*(s(i)*(fn-v(i))/2.d0+s(i)*ft*log(s(i)))*dx
c        calculate Jacobian
         bj=-bn1*((fn-v(i))/2.0+ft*log(s(i))+ft)*dx
         
c        calculate discretisation
         a(i-nox+1)=1.0/dx
         b(i-nox+1)=-2.0/dx+bj
         c(i-nox+1)=1.0/dx
         d(i-nox+1)=(-bb+bj*s(i))
        end do
c-----
c       fix boundary conditions
        d(n-nox)=d(n-nox)-c(n-nox)*s(n)
        c(n-nox)=0.0
c-----
c       solve with Gaussian elimination
        call gaus3d(a,b,c,d,r,n-nox)
c-----
        error=0.0
        ecount=0d0
c-cor   nox > nox+1
        do i=nox,n-1
c         if (abs(s(i)).ge.vzero) then
          error=error+abs((s(i)-r(i-nox+1))/(vzero+s(i)))
          ecount=ecount+1d0
c         end if
         s(i)=r(i-nox+1)
        end do
        error=error/ecount
        print*,'S===>',iter,error
        if (error.le.epsr) goto 200
      end do
      print*,'No convergence'
 200  continue
 
c     reconstitute cn based on s(i) for use in Poisson equation
      do i=nox,n
       cn(i)=s(i)*s(i)*cni
      end do
      
c-----Next Poisson solutions
      do iter=1,niterp
       do i=2,n-1
        if (i.lt.nox) then
         bb=0.0
         bj=0.0
        else
c        calculate RHS charge density
         bb=qel*(-ca(i)-cn(i)+cni*exp((fn-v(i))/ft))*dx
c        calculate Jacobian
         bj=-(qel*cni/ft)*exp((fn-v(i))/ft)*dx
        end if
	
	if (i.eq.nox) then   ! Half the volume at the interface
	  bb = bb/2d0
	  bj = bj/2d0
        end if
	
c       calculate discretisation	
        a(i-1)=ep(i-1)/dx
        b(i-1)=-((ep(i-1)+ep(i))/dx)+bj
        c(i-1)=ep(i)/dx
        d(i-1)=(-bb+bj*v(i))
	
       end do
c-----
       d(1)=d(1)-a(1)*v(1)
       a(1)=0.0
       d(n-2)=d(n-2)-c(n-2)*v(n)
       c(n-2)=0.0
c-----
       call gaus3d(a,b,c,d,r,n-2)
c-----
       error=0.0
       do i=2,n-1
        if (abs(r(i)).ge.vzero) then
         error=error+abs((v(i)-r(i-1))/r(i))
        end if
        v(i)=r(i-1)
       end do
       error=error/(n-1)
       print*,'P2==>',iter,error
       if (error.le.epsr) goto 300
      end do
      print*,'No convergence'
  300 continue
c-----
       do i=2,n-1
        if (abs(v(i)).ge.vzero) then
         error=error+abs((v(i)-v0(i))/v(i))
        end if
        v(i)=v0(i)+0.02*(v(i)-v0(i))	  ! damp the change in potential
        v0(i)=v(i)
       end do
       error=error/(n-1)
       print*,iter,error
       if (error.le.epsr_big) goto 400
       print*,'Big===>',ibig,error
      end do  ! ibig
c-----
      print*,'No big convergence'
  400 continue
  
      do i=1,n
         cband(i) = - v(i) - af(i)
      end do  
      
      open (unit=16, file='cband.dat', status='unknown')
      do i=1,n
         write (16,*) (i-nox)*dx/1e-7,cband(i)
      end do
      close (unit=16)
      
      open (unit=17, file='pot.dat', status='unknown')
      do i=1,n
         write (17,*) (i-nox)*dx/1e-7,v(i)
      end do
      close (unit=17)
      
      open (unit=32,file=filename)
      do i=nox,n
c       print*,i-nox,(i-nox)*dx/1e-7,cn(i),s0(i)
       write (32,1234)(i-nox)*dx/1e-7,cn(i),s0(i),v(i)
      end do
      
      q1=0.0
      q2=0.0
      do i=nox,n-1
       q1=q1+(s0(i)+s0(i+1))*dx/2.0
       q2=q2+(cn(i)+cn(i+1))*dx/2.0
      end do
      
      eox = (v(nox)-v(nox-1))/dx
      esi = (v(nox+1)-v(nox))/dx
      print*,ep(nox-1)/e00,eox,ep(nox-1)/e00*eox
      print*,ep(nox)/e00,esi,ep(nox)/e00*esi
      
      write (6,1234) vg,q1,q2
 1234 format (f10.5,1p2e12.5,0pf10.5)
c      write (10,*) vg,q1,q2
      write (10,*) vg,q2
      print*,filename
      
      end do
c-----
      end 
C======================================================================C
      SUBROUTINE GAUS3D(A,B,C,D,R,K)
C----------------------------------------------------------------------C
C     A. ASENOV    17.04.1984
C----------------------------------------------------------------------C
      implicit double precision (a-h,o-z)
      DIMENSION A(K),B(K),C(K),D(K),R(K)
C----------------------------------------------------------------------C
      D(1)=D(1)/B(1)
      DO 10001 I=2,K 
         B(I)=B(I)-A(I)*C(I-1)/B(I-1)
         D(I)=(D(I)-A(I)*D(I-1))/B(I)
10001 CONTINUE
      R(K)=D(K)
      DO 10002 I=K-1,1,-1
         R(I)=D(I)-C(I)*R(I+1)/B(I)
10002 CONTINUE
C----------------------------------------------------------------------C
  999 END
