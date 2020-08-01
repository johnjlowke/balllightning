c trich6 is trich5c
c trich5c has correction to make positive ions absorbed at glass
c trich5b has zero volts at i=1 for j gt. 80
c trich5 has radial conductiin of negative ins at i=32
c trich4 has cathode conducting of incident positive ions.
c diffmeta changed from 0.2 to 0.4 following value from Findlay and Snedding
c trichel3 has glass nonabsorbing or emitting of positive and negative particles through field at glass surface set to zero
c trichel2 has conducting surface at i=32 for particles from plasma flowingupstream to z=0 ie positive
c    ions from plasma flowing to glass are absorbed by glass.
c  trichel2 has metastable diffusion added
c  trichel is ballbirthmeta8.f
c  Ions and electrons impinging on cathode are absorbed;  line 239
c  ballbirthmeta8 has linear volts to zero from i > 20 to the outer edge of top plate boundary
c  ballbirthmeta8 has zero voltage at sides
c  if statement rne .ge. 0.0 then rion = 0.0   removed
c  next label 788
c  Ionization zero if ne > 1.0E14
c  ballbirth7 has far side of insulating sheet at i = 32 conducting
c  ballbirth5 has field next to glass notaveraged with high field
c  if ez.ge 0.0 alpha = 0.0
c  ballbirth5 has weav limited to 1.5E6
c  ballbirth5 has rn1 = rn - ne -ni -np -nm
c  ballbirth5 corrects read of rnm, also halves max alpha
c  ballbirth4.f has 900 MV
c  ballbirthmeta4.f has radial velocities at glass i=32 set to zero
c  ballbirthmeta3 has no input sphere of metastables, ions or electrons
c  ballbirthmeta2 has initial metastable density of 5*10^17 /cc and recombination coeff of 10^-9 cm3/s and other densities icrresed by factor 5
c  ballbirthmeta1 allows radial convection on surface of glass at i=32
c  ballbirth7 has insert of sphere of positive and negative ions
c  ballbirth.f is calculation in 2 dimensions for a stream of ions
c     on a glass sheet to produce a ball on the far side of the glass.
c  ballbirth6c has only odd points in MATLAB files for radius so number of radial points are less than 100
c  ballbirth6b has glass at i = 32
c  ballbirth6a has imposed voltage and ion density at z = 0, no current density input 
c  ballbirth6,f has drift velocity limit of 1.0E6 for electrons only
c  ballbirth5a.f has an(1) in volume divided by 2 not 8
c  ballbirth5.f has minimum ne = 1000
c  ballbirth4.f has radial velocities at glass sheet zeo
c  ballbirth3.f has boundary for V fixed at nr proportional to z.
c  ballbirth2.f has glass sheet at i = 52; insulating earth of BallEarth5c removed
c  It is a combination of BallEarth5c.f and step12b.f 
c  BallEarth5c has no ionization and treats a stream of ions impinging on the 
c      earth in two dimensions. 
c      Has plot from -R to R Space charges included in field calculation
c      Calculates electricfield and charge density of a cylindrical      
c       column of negative ions from a stepped leader being attracted   
c       to the earth by the electric field from a cloud of potential  V 
c     33 volume    used in meshls
c     34 cellpos   used in meshls 
c     36 readarcin reads input parameters from arcinls
c     37 readfn reads intial set and writes final solution
c     38 cont      writes to plot files if plot is set true
c step12b treats a stream of ions hitting glass with 
c     ionization on the far side of the glass but in one dimension.
c     STEP12b has minimum rne = 0.0000001
c     iion = 0 no ionization
c     glass at I = 781
c     electrons and positive ions included  	
c----------------------------------------------------------------------
      program ballbirth
      implicit double precision (a-h,o-z)
      logical test
c----------------------------------------------------------------------
      parameter(mr=200, mz=200, nsh=21)
      dimension rr(0:200)
      real ion
      common/var/
     r  volts,ez(0:mr,0:mz),er(0:mr,0:mz),phi(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rreabs(0:mr,0:mz),eabs(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rni(0:mr,0:mz),rjniz(0:mr,0:mz),rjnir(0:mr,0:mz),dni(0:mr,0:mz),
     r  rne(0:mr,0:mz),rjnez(0:mr,0:mz),rjner(0:mr,0:mz),dne(0:mr,0:mz),
     r  rnp(0:mr,0:mz),rjnpz(0:mr,0:mz),rjnpr(0:mr,0:mz),dnp(0:mr,0:mz),
     r  rrnp(0:mr,0:mz),rrne(0:mr,0:mz),rrni(0:mr,0:mz),rphi(0:mr,0:mz),
     r  dnm(0:mr,0:mz),rrphi(0:mr,0:mz),rrnet(0:mr,0:mz),
     r  rrnetp(0:mr,0:mz),rrnetn(0:mr,0:mz),rrnm(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlognm(0:mr,0:mz),rlognetp(0:mr,0:mz),rlognetn(0:mr,0:mz)  
c  rni,rnp,rne and rnm are negative and positive ions, electrons and metastables
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),
     r  an(0:mr), ae(0:mr,0:mz), vol(0:mr,0:mz),
     i  nr, nz
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)     
      common/plot/ iplot(50),jplot(50),iplotn,jplotn,ib,iend,jend
      common/dummy/resphi,relphi,rlxphi,mm,jscreen,iphi,jphi,iiphi,
     r	  jjphi
c---------------------------------------------------------------------
      data pid,pi/6.283185307d+00,3.141592654d+00/
c---------------------------------------------------------------------
      e=1.6E-19
      epsi=8.85E-14
      detach = 2.0E-10
      qench = 2.22E-18
      gamma=1.0E-9
      rn=2.5E19
c      diffmeta = 0.2
      diffmeta = 0.4
c Findlay and Snedding value 
      time = 0.0
      resphi = 1.0
c  ion mobility rmob is 2 where wi = rmob*e cm/s; e in V/cm
      rmob = 2.0 
      open(3,file='arcin')
      rewind 3
      read(3,*)
c volts is volts at top of ion column, ifread=1 reads arcfn,delt=timestep,
c     nr and nz are number of r and z points, mm is number of time steps
c     if iion=0 no ionization
c     nz3 is the i value beyond which currents are set to zero in earth surface
      read(3,*) volts,ifread,rlxphi,delt,nr,nz,iscreen,jscreen,mm,iion
      nz3 = nz-3
      write(*,29) mm,volts
      read(3,*) (dr(j), j=1,nr)
      read(3,*) (dz(i), i=1,nz)
      z(1) = dz(1)/2
      do 10 i=1,nz-1
          z(i+1) = z(i) + (dz(i)+dz(i+1))/2
 10   continue
      r(1) = dr(1)/2      
      do 11 j=1,nr-1
          r(j+1) = r(j) + (dr(j)+dr(j+1))/2
 11   continue
c save mesh data in file mesh
      open(20,file='mesh')
      nzp = nz+1
      nrp = nr+1
      rewind 20
        write(20,*) '     i,      z,      dz '
        do 780 i=0,nzp
            write(20,*) i,     z(i),      dz(i)
  780   continue
        write(20,*) '     j,      r,      dr'
        do 788 j=0,nrp
            write(20,*) j,     r(j),       dr(j)
  788	continue
      close(20)
      read(3,*) test
c if test true reads elements for vector plot; 
      if (test) then
c ib, iend, jend first and last ii and endj for contour plot; iplot,jplotn arrow numbers
          read(3,*) ib,iend,jend
c          read(3,*) (iplot(i), i=1,iplotn)
c          read(3,*) (jplot(j), j=1,jplotn)
      endif
      close(3)
c-----------------------------------------------------------------
c     initial conditions:
c if ifread = 1 start from file arcfn
c if ifread = 0 start from scratch.
c----------------------------------------------------------------
      if (ifread .eq. 1) then
        open(2,file='arcfn')
         read(2,*) time
         read(2,*) ((rni(j,i),j=1,nr),i=1,nz)
         read(2,*) ((rne(j,i),j=1,nr),i=1,nz)
         read(2,*) ((rnp(j,i),j=1,nr),i=1,nz)
         read(2,*) ((phi(j,i),j=1,nr),i=1,nz)
         read(2,*) ((rnm(j,i),j=1,nr),i=1,nz)
        close(2)
        endif
c initial values in ball of positive and negative ions electrons and metastables, Fig7 J Phys D 1992
c      do 205 i = 1, nz
c      do 206 j = 1,nr
c        Radd = 10.0
c	xx=r(j)
c	xx2 = xx**2
c	zz = z(i)-310.0
c	zz2 = zz**2
c	dist=sqrt(xx2+zz2)
c	if(dist.le.Radd) then
c	   rnp(j,i)=5.0e10
c	   rni(j,i)=2.5e10
c	   rne(j,i)=2.5e10
c	   rnm(j,i)=5.0e17
c	endif
c 206  continue
c 205  continue
c input ion density at z = 0
        open(4,file='rni0')
        rewind 4
        read(4,*) (rni(j,1),j=0,nr)
        close(4)
      write(*,181) time
      call volume
      open(7,file='recaxes')
        write(7,30) 
        write(7,511) (z(i),i=0,nz)
        write(7,40) 
        write(7,511) (r(j),j=0,nr)
      close(7)
c----------------------------------------------------------------
c     time iterations start here
c      time = 0.0
      do 188 j=0,nr
          rjniz(j,0) = 0.0
          rjnez(j,0) = 0.0
          rjnpz(j,0) = 0.0
 188  continue
      do 187 i=0,nz
c	  rni(0,i) = rni(1,i)
c	  phi(0,i) = phi(1,i)
          rjnir(0,i) = 0.0
          rjner(0,i) = 0.0
          rjnpr(0,i) = 0.0
 187  continue
      write(*,990)
      do 200 m=1,mm
 990    format(1x,' <solving for phi>')
      call poiss
c      write(*,850) m,resphi,jphi,iphi,relphi,jjphi,iiphi
c Ion current densities
      rmaxne = 0.0
      do 193 i=1,nz
      do 194 j=1,nr
          ion = 0.0
          att = 0.0
          rece = 0.0
          reci = 0.0
c Radial current densities
          wir = -rmob*er(j,i)
c Radial Ion and electron currents now not permitted at glass at i = 32
c          if (i.eq.32) wir = 0.0
          wire = 100*wir
          if (wire.ge.1000000.) wire = 1000000.
          if (wire.le.-1000000.) wire = -1000000.
          if(er(j,i).le.0.0) then
              rjnir(j,i) = wir*rni(j,i)*(1+(dr(j)/(2*r(j))))
              rjner(j,i) = wire*rne(j,i)*(1+(dr(j)/(2*r(j))))
              rjnpr(j,i) = -wir*rnp(j+1,i)*(1+(dr(j)/(2*r(j))))
          else
              rjnir(j,i) = wir*rni(j+1,i)*(1+(dr(j)/(2*r(j))))
              rjner(j,i) = wire*rne(j+1,i)*(1+(dr(j)/(2*r(j))))
              rjnpr(j,i) = -wir*rnp(j,i)*(1+(dr(j)/(2*r(j))))
          endif
          wiz = -rmob*ez(j,i)
c Axial Ion and electron currents zero at glass at i = 32
          if (i.eq.32) wiz = 0.0
          wize = 100*wiz
          if (wize.ge.1000000.) wize = 1000000.
          if (wize.le.-1000000.) wize = -1000000.
          if(ez(j,i).le.0.0) then
              rjniz(j,i) = wiz*rni(j,i)  
              rjnez(j,i) = wize*rne(j,i)  
              rjnpz(j,i) = -wiz*rnp(j,i+1)  
          else
              rjniz(j,i) = wiz*rni(j,i+1) 
              rjnez(j,i) = wize*rne(j,i+1) 
              rjnpz(j,i) = -wiz*rnp(j,i) 
          endif
c Increment particle densities        
          reci = gamma*rni(j,i)*rnp(j,i)
          rece = gamma*rne(j,i)*rnp(j,i)
          erav = (er(j,i)+er(j-1,i))/2
          ezav = (ez(j,i)+ez(j,i-1))/2
          eabs(j,i) = sqrt(erav**2 + ezav**2)
          if (i.eq.33) eabs(j,i) = sqrt(ez(j,i)**2)
          eav = eabs(j,i)
          eon = eav/rn
          weav = 100*rmob*eav
          if (weav.ge.1500000.0) weav = 1500000.0
c Enhanced ionization from Ball Birth paper for metastables	     
c          if (eon.ge.10.0E-17) then
c              alphaon= (eon-10.0E-17)/200.0
c          else
c              alphaon = 0.0
c          endif
c          if (eon.le.24.0E-17) then
c	      etaon = 4.0E-19 - 0.001*eon 
c          else
c              alphaon = 7.0E-19
c              etaon = 1.6E-19
c          endif
c Tam representtion of alpha/n (E,ne) from Petrova; eont = E/N(Tam)
          rn1 = rn - rne(j,i) - rnp(j,i) - rni(j,i) -rnm(j,i)
          if (rn1.le.0.0) rn1 = 0.0
          rn2 = rne(j,i)
          eont = eon/1.0E4
c omit tam = 7.66893E-19/(eont+20*eont*log((rn1+rn2)/rn1)**0.25)
          if (iion.eq.0) then
             tam = 0.0
             alphaon = 0.0
             etaon = 0.d0
          else
             tam = 7.66893E-19/eont
             alphaon = 2.42013*1.0E-16*exp(-tam)
             etaon = 0.078E-17 - 0.0029*eon
             if (etaon.le.1.6E-19) etaon = 1.6E-19
          endif
          if (eon.ge.1.0E-15) rmetaon = 1.0E-17
          if (eon.ge.5.0E-16.and.eon.le.1.0E-15) 
     r       rmetaon = 0.018*eon-0.8E-17 
          if (eon.ge.1.0E-16.and.eon.le.5.0E-16) rmetaon = 1.0E-18
          if (eon.le.1.0E-16) rmetaon = 1.0E-2*eon
          eta = etaon*rn1
          alpha = alphaon*rn1
          if (ez(j,i).ge.0.0) alpha = 0.0
          rmeta = rmetaon*rn1
c	  if (alpha.ge.1.0E4) alpha = 1.0E4
          if (alpha.ge.0.5E4) alpha = 0.5E4
          att = rne(j,i)*eta*weav
c correction to make attachment rate constant near zero electric field
          if (eon.le.3.0E-17) att = rne(j,i)*3.75E6
          rion = rne(j,i)*alpha*weav
c	  if (rne(j,i).ge.0.0) rion = 0.0
          rmet = rne(j,i)*rmeta*weav
          det = detach*rnm(j,i)*rni(j,i)
          qen = qench*rnm(j,i)*rn/5
          dni(j,i) = -(rjnir(j,i)-rjnir(j-1,i))/dr(j) - 
     r        (rjniz(j,i)-rjniz(j,i-1))/dz(i)+att-reci - det
          if (i.ge.2) rni(j,i) = rni(j,i) + dni(j,i)*delt
          dne(j,i) = -(rjner(j,i)-rjner(j-1,i))/dr(j) - 
     r        (rjnez(j,i)-rjnez(j,i-1))/dz(i)+rion-att-rece + det
          if (i.ge.2) rne(j,i) = rne(j,i) + dne(j,i)*delt
          if (rne(j,i).le.1.0E-7.and.i.gt.32) rne(j,i) =1.0E-7
          if (i.le.32) rne(j,i) =1.0E-7
          if (iion.eq.0) rne(j,i) =1.d0
          if (rne(j,i).ge.rmaxne) rmaxne = rne(j,i)
c Make under surface of cathode absorb positive ions - corrected
c          if (i.eq.33.and.rjnpz(j,i).le.0.0) rjnpz(j,i-1)=rjnpz(j,i) 
          if (i.eq.33.and.rjnpz(j,i).le.0.0) then
             wiz = -rmob*ez(j,32)     
             rjnpz(j,i-1)= -wiz*rnp(j,i) 
          endif   
          dnp(j,i) = -(rjnpr(j,i)-rjnpr(j-1,i))/dr(j) - 
     r        (rjnpz(j,i)-rjnpz(j,i-1))/dz(i)+rion-rece-reci
          if (i.ge.2) rnp(j,i) = rnp(j,i) + dnp(j,i)*delt
c Diffusion of metastables; diff coefficient = diffmeta
c	  dnm(j,i) = rmet -det -qen
          diffmet = diffmeta*(((rnm(j,i+1)-rnm(j,i))/(z(i+1)-z(i))
     r	  -(rnm(j,i)-rnm(j,i-1))/(z(i)-z(i-1)))/dz(i)
     r    +((rnm(j+1,i)-rnm(j,i))/(r(j+1)-r(j))
     r    -(rnm(j,i)-rnm(j-1,i))/(r(j)-r(j-1)))/dr(j))
	  dnm(j,i) = rmet -det -qen +diffmet
	  rnm(j,i) = rnm(j,i) + dnm(j,i)*delt
          if (iion.eq.0) rnm(j,i) =1.d0
	  rnet(j,i) = rnp(j,i) - rne(j,i) - rni(j,i)
 194  continue
 193  continue
      do 196 j=1,nr
          rni(j,nr+1) = rni(j,nr) 
 196  continue
      time = time + delt
      write(*,851) m,time,rmaxne
 200  continue
c     end internal iteration
      open(304,file='arcfnout')
        rewind 304
         write(304,520) time
         write(304,520) ((rni(j,i),j=1,nr),i=1,nz)
         write(304,520) ((rne(j,i),j=1,nr),i=1,nz)
         write(304,520) ((rnp(j,i),j=1,nr),i=1,nz)
         write(304,520) ((phi(j,i),j=1,nr),i=1,nz)
         write(304,520) ((rnm(j,i),j=1,nr),i=1,nz)
      close(304)
  25  format(1x,'rlxn. factors =',6(f6.2,','),' no. of iterations  
     r         =',i5)
  29  format(' mm =',i6,' volts =',e10.3)
  30  format(' axial values of z')
  40  format(' radial values of r')
  50  format(1p6e13.5)
  53  format(5(e11.4))
 181  format(' time= ',e10.4,' dt = ',e10.4,'  mm = ',i6)
 463  format(2x,5(3x,f6.3),3x,i5)
 495  format(2e15.5)
 496  format(2e15.5)
 500  format(8f9.0)
 501  format(/)
 502  format(8f9.0)
 503  format(/'  rjniz')
 504  format(/'  rjnez')
 505  format(/'  rjnpz')
 507  format(/'  potential')
 508  format(/'  z  electric field')
 509  format(/'  r  electric field')
 510  format('  i =',i4,'   z =',1pe12.3)
 511  format(1p10e10.2)
 512  format(14f6.0)
 513  format(11f8.0)
 514  format(/'  Negative Ion Density')
 522  format(/'  Electron Density')
 523  format(/'  Positive Ion Density')
 516  format(9f7.0)
 518  format (3e24.8)
 520  format (5e16.8)
 521  format (1pe10.1,1pe10.1,1pe10.1,1pe10.1,1pe10.1,1pe10.1,
     1       1pe10.1,1pe10.1,1pe10.1,1pe10.1)
 524  format (/' Net Charge Density')
 527  format (/' Metastable Density')
 572  format(1p6e12.4)
 620  format(1i3,6e12.3)
 796  format(200e11.3)
 797  format(200f9.3)
 850  format(' iteration no. ',i6,'   resphi =',1pe9.2,2i3,
     r    ' relphi =',e9.2,2i3/)
 851  format(1x,'iteration no.',i6,'  time =',e12.4,
     r '  max electron density =',1pe9.2)
 859  format('     res    je     ie      phi        rel    jje  ',
     r 'iie    phi') 
 861  format(1pe12.2,2i4,1pe12.2,1pe12.2,2i4,e12.2)
 899     format(1x,'iteration no. ',i4)
 909  format(1x,i3,6(2pe12.4))
 918  format(2(1p,i4),2(2p,f12.4))
 920  format((1p,i4),2(2p,f12.4))
 978  format(i3.3)
 995  format(1p6e13.5)
 996  format(1pe10.2)
       open(6,file='record')
       rewind 6
        write(6,181) time,delt,mm
        write(6,859)
        write(6,861) resphi,jphi,iphi,phi(jphi,iphi),relphi,jjphi,
     1        iiphi,phi(jjphi,iiphi)
        write(6,40) 
        write(6,511) (r(j),j=0,nr)
         write(6,507)
      do 778 i=0,nz+1
         write(6,510) i,z(i)
         write(6,511) (phi(j,i),j=0,nr+1)
 778  continue
         write(6,501)
         write(6,514)
      do 785 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rni(j,i),j=0,nr+1)
 785  continue
         write(6,501)
         write(6,522)
      do 786 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rne(j,i),j=0,nr+1)
 786  continue
         write(6,501)
         write(6,523)
      do 787 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rnp(j,i),j=0,nr+1)
 787  continue
         write(6,501)
         write(6,524)
      do 772 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rnet(j,i),j=0,nr+1)
 772  continue
         write(6,501)
         write(6,527)
      do 798 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rnm(j,i),j=0,nr+1)
 798  continue
         write(6,501)
         write(6,508)
      do 779 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (ez(j,i),j=0,nr+1)
 779  continue
         write(6,501)
         write(6,509)
      do 781 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (er(j,i),j=0,nr+1)
 781  continue
         write(6,501)
         write(6,503)
      do 782 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rjniz(j,i),j=0,nr+1)
 782  continue
         write(6,501)
         write(6,504)
      do 783 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rjnez(j,i),j=0,nr+1)
 783  continue
         write(6,501)
         write(6,505)
      do 784 i=0,nz
         write(6,510) i,z(i)
         write(6,511) (rjnpz(j,i),j=0,nr+1)
 784  continue
       close(6)
c files for MATLAB
      ib = 1
      iend = nz
      jend = nr
      open(21,file='z.dat')
      write(21,797) (z(i),i=ib,iend)
      close(21)
c err and rr files for full plot from r = -R to r = +R
      i = ib
      jjend = 2*jend
      n = 0
      do 789 j=1,jjend
	  n = n + 1
          jj = jend-j+1
          if (jj.le.0) then
	      jj = j - jend
              rr(n) = r(jj) 
          else
              rr(n) = -r(jj)
          endif
 789   continue
      rr(0) = 0.d0   
       open(23,file='rr.dat')
       rewind 23
      write(23,996) (rr(j),j=1,jjend)
      close(23)
      open(28,file='rrni.dat')
      do 799 i=ib,iend
      do 800 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j-jend
          rrni(j,i) = rni(jj,i)
          if (rrni(j,i).le.1.0E-7) rrni(j,i) = 1.0E-7
          rlogni(j,i) = log10(rrni(j,i)) 
 800  continue
 799  continue
      do 793 i=ib,iend
          write(28,796) (rlogni(j,i),j=1,jjend)
 793  continue
      open(30,file='rrnp.dat')
      do 806 i=ib,iend
      do 807 j=1,jjend
          jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rrnp(j,i) = rnp(jj,i)
          if (rrnp(j,i).le.1.0D-7) rrnp(j,i)= 1.0D-7
          rlognp(j,i) = log10(rrnp(j,i))
  807 continue
  806 continue
      do 808 i=ib,iend
          write(30,796) (rlognp(j,i),j=1,jjend)
  808 continue      
      open(26,file='rrne.dat')
      do 794 i=ib,iend
      do 801 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rrne(j,i) = rne(jj,i) 
	  if (rrne(n,i).le.1.0D-7) rrne(n,i)=1.0D-7
	  rlogne(j,i) = log10(rrne(j,i))
 801  continue
 794  continue
      do 802 i=ib,iend
          write(26,796) (rlogne(j,i),j=1,jjend)
 802  continue
       open(25,file='rrmeta.dat')
      do 810 i=ib,iend
      do 811 j=1,jjend
          jj = jend-j+1
          if (jj .le.0) jj = j-jend
          rrnm(j,i) = rnm(jj,i)
          if (rrnm(j,i).le.1.0E-7) rrnm(j,i)=1.0E-7
          rlognm(j,i) = log10(rrnm(j,i))
 811  continue
 810  continue
      do 812 i = ib,iend
          write(25,796) (rlognm(j,i),j=1,jjend)
 812  continue 
      close (25)
      close (26)
      close (27)
      close (28)
      close (30)
      do 795 i=ib,iend
      do 809 j=1,jjend 
          jj = jend-j+1
          if (jj .le.0) jj = j-jend
          rrnet(j,i) = rnet(jj,i)
          if (rrnet(j,i).le.1.0D-7.and.rrnet(j,i).ge.-1.0D-7)
     r        rrnet(j,i) = 1.0D-7
          if (rrnet(j,i).gt.0.0) then
              rrnetp(j,i) = rrnet(j,i)
              rlognetp(j,i) = log10(rrnetp(j,i))
          else
              rrnetn(j,i) = -rrnet(j,i)
              rlognetn(j,i) = log10(rrnetn(j,i))
          endif
 809  continue
 795  continue
      open(31,file='netp.dat')
      open(32,file='netn.dat') 
      do 813 i = ib,iend
          write(31,796) (rlognetp(j,i),j=1,jjend)
 813  continue
      do 814 i = ib,iend
          write(32,796) (rlognetn(j,i),j=1,jjend)
 814  continue
      close (31)
      close (32)       
      open(33,file='rreabs.dat')
      do 815 i=ib,iend
      do 816 j=1,jjend
          jj = jend-j+1
          if (jj .le.0) jj = j-jend
          rreabs(j,i) = eabs(jj,i)
 816  continue
 815  continue
      do 817 i = ib,iend
          write(33,796) (rreabs(j,i),j=1,jjend)
 817  continue
      close (33) 
      open(34,file='rrphi.dat')
      do 818 i=ib,iend
      do 819 j=1,jjend
          jj = jend-j+1
          if (jj .le.0) jj = j-jend
          rrphi(j,i) = phi(jj,i)
 819  continue
 818  continue
      do 820 i = ib,iend
          write(34,796) (rrphi(j,i),j=1,jjend)
 820  continue
      close (34)   
c      pause
      end 
c-------------------------------------------------
c  end main
c-------------------------------------------------------------------------------------------------------
c     subroutine poiss
c
c     this subroutine solves Poisson's equation to calculate the
c     electric potential (phi) for a column of ions attracted to      
c     the earth by a potential V on the thunder cloud. The axial and
c     radial electric fields are also calculated(cdr)
c     poiss uses the subroutines fill for the coeff for the phi 
c     equation and solve to solve the phi equation line by line.
c-----------------------------------------------------------------------
      subroutine poiss
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=200, nsh=21)
      common/var/
     r  volts,ez(0:mr,0:mz),er(0:mr,0:mz),phi(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rreabs(0:mr,0:mz),eabs(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rni(0:mr,0:mz),rjniz(0:mr,0:mz),rjnir(0:mr,0:mz),dni(0:mr,0:mz),
     r  rne(0:mr,0:mz),rjnez(0:mr,0:mz),rjner(0:mr,0:mz),dne(0:mr,0:mz),
     r  rnp(0:mr,0:mz),rjnpz(0:mr,0:mz),rjnpr(0:mr,0:mz),dnp(0:mr,0:mz),
     r  rrnp(0:mr,0:mz),rrne(0:mr,0:mz),rrni(0:mr,0:mz),rphi(0:mr,0:mz),
     r  dnm(0:mr,0:mz),rrphi(0:mr,0:mz),rrnet(0:mr,0:mz),
     r  rrnetp(0:mr,0:mz),rrnetn(0:mr,0:mz),rrnm(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlognm(0:mr,0:mz),rlognetp(0:mr,0:mz),rlognetn(0:mr,0:mz)  
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),
     r  an(0:mr), ae(0:mr,0:mz), vol(0:mr,0:mz),
     i  nr, nz
      common/dummy/resphi,relphi,rlxphi,mm,jscreen,iphi,jphi,iiphi,
     r	  jjphi
c---calculate coefficients  
      call fill
c---relaxation
      do 15 i=1,nz
      do 15 j=1,nr
         dc(j,i) = dc(j,i)/rlxphi
         df(j,i) = df(j,i) + (1.d0-rlxphi)*dc(j,i)*phi(j,i)
  15  continue
c---solve matrix
  17      call solve(1,nr,1,nz,de,dw,dn,ds,dc,df,phi,0,jscreen)
c--reset boundries
      do 30 i=0,nz+1
         phi(nr+1,i)=phi(nr,i)
         phi(0,i)=phi(1,i)
  30  continue

c--calculate residual
      resphi = 0.d0
      relphi = 0.d0
      do 40 i=2,nz-1
      do 40 j=2,nr
          east= de(j,i)*phi(j+1,i)/dc(j,i)
          west= dw(j,i)*phi(j-1,i)/dc(j,i)
          rnorth= dn(j,i)*phi(j,i+1)/dc(j,i)
          south= ds(j,i)*phi(j,i-1)/dc(j,i)
          centre= df(j,i)/dc(j,i)
	  ph = dabs(phi(j,i))
          res = east+west+south+rnorth+centre- phi(j,i)
          res = dabs(res)
        if( abs(phi(j,i)) .ge. 1e-6 ) then 
		rel=res/ph
        else
		rel = res/1e-6
        endif
        if (res .gt. resphi) then
           resphi = res
           iphi = i
           jphi = j
        end if
        if (rel .gt.relphi) then
            relphi = rel
            iiphi = i
            jjphi = j
        endif
c----end of residu calculation
 40   continue

c--- calculate the electric field;
      do 70 i=1,nz-1
      do 70 j=1,nr
         ez(j,i) = 2*(phi(j,i)-phi(j,i+1))/(dz(i)+dz(i+1))
  70  continue
      do 71 j=1,nr
         ez(j,nz) = 2*phi(j,nz)/dz(nz)
  71  continue
c----calculate the electric field
      do 80 i=1,nz
      do 80 j=1,nr
         er(j,i) = 2*(phi(j,i)-phi(j+1,i))/(dr(j)+dr(j+1))
  80  continue
      do 75 i=0,nz
         er(0,i) = 0.0
  75  continue
      do 76 j=0,nr
         ez(j,0) = ez(j,1)
  76  continue
 100  continue
      return
      end
c------------------------------------------------------------------
c end of subroutine poiss
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c     subroutine fill
c     fill calculates coefficients needed in the calculation 
c     of the electric potential in poiss. 
c----------------------------------------------------------------------
      subroutine fill
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=200, nsh=21)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),
     r  an(0:mr), ae(0:mr,0:mz), vol(0:mr,0:mz),
     i  nr,nz
      common/var/
     r  volts,ez(0:mr,0:mz),er(0:mr,0:mz),phi(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rreabs(0:mr,0:mz),eabs(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rni(0:mr,0:mz),rjniz(0:mr,0:mz),rjnir(0:mr,0:mz),dni(0:mr,0:mz),
     r  rne(0:mr,0:mz),rjnez(0:mr,0:mz),rjner(0:mr,0:mz),dne(0:mr,0:mz),
     r  rnp(0:mr,0:mz),rjnpz(0:mr,0:mz),rjnpr(0:mr,0:mz),dnp(0:mr,0:mz),
     r  rrnp(0:mr,0:mz),rrne(0:mr,0:mz),rrni(0:mr,0:mz),rphi(0:mr,0:mz),
     r  dnm(0:mr,0:mz),rrphi(0:mr,0:mz),rrnet(0:mr,0:mz),
     r  rrnetp(0:mr,0:mz),rrnetn(0:mr,0:mz),rrnm(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlognm(0:mr,0:mz),rlognetp(0:mr,0:mz),rlognetn(0:mr,0:mz)  
c----------------------------------------------------------------------
c     setting up interpolated values
c     electron charge e(C) and physical constant epsilon
      e  = 1.60206d-19
      epsilon = 8.85d-14
      do 10 i=1,nz
      do 10 j=1,nr
         dife = 2.*ae(j,i)/(dr(j)+dr(j+1))
         de(j,i) = dife
         dw(j+1,i) = dife
         difn = 2.*an(j)/(dz(i)+dz(i+1))
         dn(j,i) = difn
         ds(j,i+1) = difn
c inclusion of charge density in Poisson Equation
         df(j,i) = (e/epsilon)*(rnp(j,i)-rne(j,i)-rni(j,i))*vol(j,i)
  10  continue
c     setting up boundary values
      do 20 i=1,nz
         dw(1,i) = 0.d0
c         de(nr,i) = 0.d0
 20   continue
      do 40 i=2,nz
      do 39 j=1,nr
         dc(j,i) = dn(j,i) + ds(j,i) + de(j,i) + dw(j,i)
 39   continue
 40   continue
      do 30 j=1,nr
c boundary voltage on top plate beyong j = 20 to scale linearly to zero         
          voltr = volts
          if (j.gt.20) voltr = volts*(r(nr)-r(j))/(r(nr)-r(20)) 
          if (j.gt.80) voltr = 0.0 
          call zerocoef(dc,df,de,dw,dn,ds,voltr,j,1)
          call zerocoef(dc,df,de,dw,dn,ds,0.d0,j,nz)
 30   continue
       do 21 i=1,nz
         call zerocoef(dc,df,de,dw,dn,ds,0.d0,nr,i)
 21   continue
c       do 21 i=1,nz
c         vv = volts*(1-(z(i)-z(1))/(z(nz)-z(1)))
c         call zerocoef(dc,df,de,dw,dn,ds,vv,nr,i)
c 21   continue
 520     format (5e16.8)
 501  format(/)
c	 write(305,501)
c         write(305,520) ((df(j,i),j=1,nr),i=1,nz)
c      close(305)
      return
      end
c end of subroutine fill
c-----------------------------------------------------------------
c     subroutine solve
c     v.1 ls 18.9.98
c
c     solve is the solver for the poisson equation
c     solve uses the line by line technique in conjonction with 
c     the block correction technique.
c-----------------------------------------------------------------

      subroutine solve(m,nr,n,nz,de,dw,dn,ds,dc,df,phi,flag1,jscreen)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=200, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1   df(mr,mz),phi(0:mr,0:mz),v(mz),x(0:mr,0:mz) 
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz) 
      integer flag1
c------------------------------------------------------------------  
c      write(*,29) jscreen
      nr1 = nr-1
      nz1 = nz-1
      do 5 i=0,mz
      do 5 j=0,mr
   5     x(j,i) = 0.d0
      do 10 i=0,nz+1
      do 10 j=0,nr+1
         x(j,i) = phi(j,i)
  10  continue
      rel = 1.85
  29  format(' jscreen =',i4)
      do 200 iter=1,40
c
c      first step: block correction along lines of constant j
c
      do 100 j=m,nr1
         a(j) = 0.d0
         b(j) = 0.d0
         c(j) = 0.d0
         d(j) = 0.d0
 100   continue
      do 110 i=n,nz
      do 110 j=m,nr1
         a(j) = a(j) + dc(j,i) - dn(j,i) - ds(j,i)
         b(j) = b(j) + de(j,i)
         c(j) = c(j) + dw(j,i)
         d(j) = d(j) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 110  continue
      call tdma(m,nr1,v)
      do 120 i=n,nz
      do 120 j=m,nr1
         x(j,i) = x(j,i) + v(j)
 120  continue
c
c     second step: block corrections along lines of constant i
c
      do 130 i=n,nz1
         a(i) = 0.d0
         b(i) = 0.d0
         c(i) = 0.d0
         d(i) = 0.d0
 130  continue
      do 140 i=n,nz1
      do 140 j=m,nr
         a(i) = a(i) + dc(j,i) - de(j,i) - dw(j,i)
         b(i) = b(i) + dn(j,i)
         c(i) = c(i) + ds(j,i)
         d(i) = d(i) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 140  continue
      call tdma(n,nz1,v)
      do 150 i=n,nz1
      do 150 j=m,nr
         x(j,i) = x(j,i) + v(i)
 150  continue
c
c     third step: sweep from i=1 to nz
c
      do 30 i=n,nz
         do 35 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*dn(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*(x(j,i+1) - (rel-1.)*x(j,i))
     1              + ds(j,i)*x(j,i-1)
  35     continue
         call tdma(m,nr,v)
         do 40 j=m,nr
            x(j,i) = v(j)
  40     continue
  30  continue
c
c     fourth step: sweep from j=nr to 1
c
      do 70 j=nr,m,-1
         do 75 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*dw(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + de(j,i)*x(j+1,i)
     1             + dw(j,i)*(x(j-1,i) - (rel-1.)*x(j,i))
  75     continue
         call tdma(n,nz,v)
         do 80 i=n,nz
            x(j,i) = v(i)
  80     continue
  70  continue
c
c     fifth step: sweep from i=nz to 1
c
      do 50 i=nz,n,-1
         do 55 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*ds(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*x(j,i+1) + ds(j,i)*(x(j,i-1) -
     1             (rel-1.)*x(j,i))
  55     continue
         call tdma(m,nr,v)
         do 60 j=m,nr
            x(j,i) = v(j)
  60     continue
  50  continue
c
c     sixth step: sweep from j=1 to nr
c
      do 90 j=m,nr
         do 95 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*de(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + dw(j,i)*x(j-1,i) + de(j,i)*(x(j+1,i) -
     1             (rel-1.)*x(j,i))
  95     continue
         call tdma(n,nz,v)
         do 96 i=n,nz
            x(j,i) = v(i)
  96     continue
  90  continue
c
c     check convergence
c
       res = 0.d0
       do 160 i=2,nz-1
       do 160 j=1,nr
	   rres =(de(j,i)*x(j+1,i)+dw(j,i)*x(j-1,i)+dn(j,i)*x(j,i+1)
     r    +ds(j,i)*x(j,i-1)+df(j,i)-dc(j,i)*x(j,i))/dc(j,i)
	  rres = dabs(rres)
	  if ( rres.gt.res) res = rres
  160  continue     
       reldiff = 0.0
         do 260 i=2,nz
         do 260 j=2,nr
            if (x(j,i) .ne. 0.) then
               diff = (x(j,i)-phi(j,i))/x(j,i)
            else
               diff = -phi(j,i)
            end if
	    diff = dabs(diff)
            if (diff.gt.reldiff) then
		reldiff = diff
		ii = i
		jj = j
	    endif
  260    continue
  998    format(1x,'iter =',i4,'    res =',1pe12.3,
     1            '   reldiff  =',1pe12.3,2i4)
         do 270 i=n,nz
         do 270 j=m,nr
            phi(j,i) = x(j,i)
  270    continue
      if (jscreen.eq.1) write(*,998) iter,res,reldiff,jj,ii
      if (reldiff .lt. 1.e-3) return
 200  continue
      write(*,300)
c      write(5,300)
 300  format(1x,'convergence not attained in subroutine solve')
      return
      end
c------------------------------------------------------------------
c end of subroutine solve
c-----------------------------------------------------------------------

c     subroutine tdma
c
c     tdma solves equations of the form
c     a(i)*v(i) = b(i)*v(i+1) + c(i)*v(i-1) + d(i)  ,i=1,....,n
c     with c(1) = 0 and b(n) = 0
c----------------------------------------------------------------------
      subroutine tdma(if,l,v)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=200, nsh=21)
      dimension v(mz),p(mz),q(mz)
      common/tri/ a(mz),b(mz),c(mz),d(mz)  
c----------------------------------------------------------------------
      l1 = l-1
      p(if) = b(if)/a(if)
      q(if) = d(if)/a(if)
      ifp1 = if + 1
      do 10 i=ifp1,l
         den = 1./(a(i) - c(i)*p(i-1))
         p(i) = b(i)*den
         q(i) = (d(i)+c(i)*q(i-1))*den
  10  continue
      v(l) = q(l)
      do 20 i=if,l1
         k = l1 -i + if
         v(k) = p(k)*v(k+1) + q(k)
  20  continue
      return
      end
c-------------------------------------------------------
c end subroutine tdma
c-----------------------------------------------------------
      subroutine zerocoef(dc,df,de,dw,dn,ds,dfval,j,i) 
c sets value at j,i to be dfval
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=200, nsh=21)
      dimension  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1    df(mr,mz)
         dc(j,i)=1.d0
         df(j,i)=dfval
         de(j,i)=0.d0
         dw(j,i)=0.d0
         dn(j,i)=0.d0
         ds(j,i)=0.d0
        end
c end of subroutine zerocoef
c-----------------------------------------------------
c     subroutine volume
c 
c calculates the volume and the surfaces of each elementary
c cell. calculates also the ponderation factor used for the 
c geometric average calculation.
c an top surface of the cell
c ae external radial surface of the cell
c vol volume of the cell
c---------------------------------------------------- 

      subroutine volume
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=200, nsh=21)
        common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),
     r  an(0:mr), ae(0:mr,0:mz), vol(0:mr,0:mz),
     i  nr,nz
      an(0) = 0.d0
      an(1) = 0.5*dr(1)*dr(1)
      do 30 j=2,nr
         xa = r(j)*dr(j)
         an(j) = xa
  30  continue
      do 50 i=0,nz
      do 50 j=0,nr
         xa = (r(j)+0.5*dr(j))*dz(i)
         ae(j,i) = xa
  50  continue
      do 90 i=1,nz
      do 90 j=1,nr
        vol(j,i) = an(j)*dz(i)
  90  continue

      return
      end
c------------------------------------------------------
c end subroutine volume
c------------------------------------------------------
     
