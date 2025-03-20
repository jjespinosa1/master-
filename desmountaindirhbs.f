c======================================================================c
c 
c     PROGRAM DIRHB-Spherical  
c
c======================================================================c
c     Relativistic Hartree-Bogoliubov theory in a spherical basis
c     Main part of the code
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
c
c---- sets data
      call default
c
c---- reads in data     
      call reader
c
c---- force-parameters
      call forces(.true.)
c
c---- Gauss-Hermite mesh points
      call gaush(.true.)
c
c---- preparations
      call prep(.false.)
c
c---- initialization of the potentials
      call inout(1,.true.)
      call start(.true.)
c
c---- oscillator basis for single particle states
      call base(.true.)
c
c---- initialization of the pairing field
      call dinout(1,.false.)
c
c---- wavefunctions at Gauss-Meshpoints (RN d^r Rn , Rb)
      call gaupol(.false.)
c
c---- single-particle matix elements
      call singf(.false.)
c
c---- pairing matix elements
      call singd(.false.)
c      
c---- coulomb and meson propagators
      call greecou(.false.)
      call greemes(.false.)
c
c---- iteration
C      call iter(.true.)
c
c---- transformation to the canonical basis
C      call canon(.true.)
c
c---- center of the mass correction
C      call centmas(.false.)
c
c---- results
C      call resu(.true.)
c
c---- densities in the coordinate space
C      call plot(.false.)
c---- punching of the potentials to the tape  dirhb.wel
C      call inout(2,.false.)
c---- punching of the pairing field to the tape  dirhb.wel      
C      call dinout(2,.false.)
      stop ' FINAL STOP OF DIRHBS'
c-end-DIRHBS
      end



c======================================================================c

      subroutine default()

c======================================================================c
c
c     Default for Relativistic Mean Field spherical
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tp*1,tis*1,tit*8,tl*1                   ! textex

      common /baspar/ hom,hb0,b0
      common /couplf/ ff(0:ngh,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /fermi / ala(2),tz(2)
      common /initia/ inin,inink
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /broyde2/ ibroyd
c

c---- fixed texts
      data tp/'+','-'/,tis/'n','p'/,tit/'Neutron:','Proton: '/
      data tl/'s','p','d','f','g','h','i','j','k','l','m',
     &            'n','o','P','q','r','S','t','u','v','w',
     &            'x','y','z','0','0','0','0','0','0','0'/
c
c---- physical constants
      data hbc/197.328284d0/,r0/1.2d0/,alphi/137.03602/
c======================================================================c
c---- signs and factorials
c-----------------------------------------------------------------------
      call gfv()
c======================================================================c
c---- Neutrons and Protons
c-----------------------------------------------------------------------
c     number of isospins 
c     itx = 1   only neutrons
c           2   neutrons and protons
c
      itx = 2
c
c======================================================================c
 


c======================================================================c
c     center off mass correction: 
c-----------------------------------------------------------------------
c
c     a) the kinetic energy in the variation can be
c           <T>  = hbc**2/(2*amu)*Delta   (no correction)
c        or <T>' = hbc**2/(2*amu)*(1-1/A)*Delta
c                  the diagonal part of ecm=<P**2>/2Am is included
c     b) the total energy is 
c        <T> + <V> - <P**2>/2Am
c
c     icm: 0   variation of <T>
c              hb0 = hbc**2/(2*amu)       ecm = 3/4*hom
c          1   variation of <T>'
c              hb0 = hb0*(1-1/A)          ecm = 3/4*hom
c          2   variation of <T>
c              hb0 = hbc**2/(2*amu)       ecm = <P**2>/2M 
c
      icm = 0       
  
c
c======================================================================c



c======================================================================c
c     Coulomb-Field
c-----------------------------------------------------------------------
c     icou: Coulomb-field:  0   not at all 
c                           1   only direct term  
c                           2   plus exchange 
      icou = 1
c======================================================================c

C======================================================================c
c     pairing
c----------------------------------------------------------------------c
      do it = 1,2                   
	     spk0(it) = zero
	     del(it)  = zero
	     ala(it)  = -7.0
      enddo   ! it
c======================================================================c
c     iteration
c----------------------------------------------------------------------c
c
      maxi = 500               ! maximal number of iteration
      si   = one             ! actual error in the main iteration
      epsi = 1.d-6       ! accuracy for the main iteration
      iaut  = 1              ! automatic change of xmix: 0 (no) 1 (yes)
      inxt  = 3              ! number of iterations until next question
      xmix  = 0.5d0          ! starting xmix
      xmax  = 0.7d0          ! maximal value for xmix
      xmix0 = xmix           ! basic xmix, where it returns to
      ibroyd= 1
c======================================================================c
c
c
c---- parameters of the initial potentials
c     inin = 0: fields read, 1: saxon-wood,
      inin  = 1         
c     inink = 0: pairing potential read, 1: pairing potential monopol
      inink = 1
c
c     oscillator length b0 (is calcuated for b0 <= 0)
      b0 = -2.320
c
c


c======================================================================c
c---- tapes
c----------------------------------------------------------------------c
      l6   = 10
      lin  = 3
      lou  = 6
      lwin = 1
      lwou = 2
      lplo = 11
      laka = 12
      lvpp = 13
c======================================================================c
c---- preparation of density dependence
c----------------------------------------------------------------------c
      do m = 1,4
        f=0
         lmes(m) = 0
         ggmes(m) = zero
         do i = 0,ngh
            ff(i,m,1) = one
            !f = f + 1
           ! print*, 'Counter ',f
            !print*, 'ff(i,m,1) = ',ff(i,m,1)    
	    ff(i,m,2) = zero
        !print*, 'len',f
         enddo   ! ngh
      enddo   ! m

      return
c-end-DEFAULT
      end

C=======================================================================

      subroutine gfv

C=======================================================================
C
C     Calculates sign, sqrt, factorials, etc. of integers and half int.
C
c     iv(n)  =  (-1)**n
c     sq(n)  =  sqrt(n)
c     sqi(n) =  1/sqrt(n)
c     sqh(n) =  sqrt(n+1/2)
c     shi(n) =  1/sqrt(n+1/2)
c     fak(n) =  n!
c     fad(n) =  (2*n+1)!!
c     fdi(n) =  1/(2*n+1)!!
c     fi(n)  =  1/n!
c     wf(n)  =  sqrt(n!)
c     wfi(n) =  1/sqrt(n!)
c     wfd(n) =  sqrt((2*n+1)!!)
c     gm2(n) =  gamma(n+1/2)
c     gmi(n) =  1/gamma(n+1/2)
c     wg(n)  =  sqrt(gamma(n+1/2))
c     wgi(n) =  1/sqrt(gamma(n+1/2))
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
c     include 'diz.par'
      parameter (igfv = 100)
c
      common /gfviv / iv(0:igfv)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvfak/ fak(0:igfv)
      common /gfvfad/ fad(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvfdi/ fdi(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfvwfd/ wfd(0:igfv)
      common /gfvgm2/ gm2(0:igfv)
      common /gfvgmi/ gmi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
c
c---- mathemathical constants
c     data zero/0.0d0/,one/1.d0/,two/2.d0/
c     data half/0.5d0/,third/0.333333333333333333d0/
c     data pi/3.141592653589793d0/
c
      zero  = 0.d0
      one   = 1.d0
      two   = 2.d0
      half  = one/two
      third = one/3.d0
      pi    = 4*atan(one)
c
      iv(0)  = +1
      sq(0)  =  zero
      sqi(0) =  1.d30
      sqh(0) =  sqrt(half)
      shi(0) =  1/sqh(0)
      fak(0) =  one
      fad(0) =  one
      fi(0)  =  one
      fdi(0) =  one
      wf(0)  =  one
      wfi(0) =  one
      wfd(0)=  one
c     gm2(0) = Gamma(1/2) = sqrt(pi)
      gm2(0) =  sqrt(pi)
      gmi(0) =  1/gm2(0)
      wg(0)  =  sqrt(gm2(0))
      wgi(0) =  1/wg(0)
     
      do i = 1,igfv
        iv(i)  = -iv(i-1)
        sq(i)  = dsqrt(dfloat(i))
        sqi(i) = one/sq(i)
        sqh(i) = sqrt(i+half)
        shi(i) = one/sqh(i)
        fak(i) = i*fak(i-1)
        fad(i) = (2*i+1)*fad(i-1)
        fi(i)  = one/fak(i)
        fdi(i) = one/fad(i)
        wf(i)  = sq(i)*wf(i-1)
        wfi(i) = one/wf(i)
        wfd(i) = sqrt(fad(i))
        gm2(i) = (i-half)*gm2(i-1)
        gmi(i) = one/gm2(i)
        wg(i)  = sqh(i-1)*wg(i-1)
        wgi(i) = one/wg(i)
        enddo
c
c     write(6,*) ' ****** END GFV *************************************' 
      return
c-end-GFV
      end

c=======================================================================
c======================================================================c

      subroutine reader

c======================================================================c
          implicit real*8 (a-h,o-z)
c
          include 'dirhb.par'
c 
          character parname*10                                      ! partyp
          character tp*1,tis*1,tit*8,tl*1                           ! textex
          character nucnam*2                                        ! nucnuc
    
          CHARACTER*8 date
          CHARACTER*10 time
          CHARACTER*5 zone
          INTEGER*4 VALUES(8)
c
          common /basnnn/ n0f,n0b
          common /baspar/ hom,hb0,b0
          common /fermi / ala(2),tz(2)
          common /initia/ inin,inink
          common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
          common /mathco/ zero,one,two,half,third,pi
          common /nucnuc/ amas,nama,nneu,npro,nucnam
          common /optopt/ itx,icm,icou,ipc,inl,idd
          common /partyp/ parname
          common /pair  / del(2),spk(2),spk0(2)
          common /physco/ hbc,alphi,r0
          common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
          common /textex/ tp(2),tis(2),tit(2),tl(0:30)
          common /wplot/ itpl,jpl,ippl,kpl      
c............................
caf ... reading from dirhb.dat
c............................
c
          if (lin.eq.0) return
          open(lin,file='dirhb.dat',status='old')
c
c
c---- Basisparameters:            
          read(lin,'(10x,2i5)') n0f,n0b
c
c---- Initialization of wavefunctions:
          read(lin,'(10x,2i5)') inin,inink
c
c---- Nucleus under consideration
          read(lin,'(a2,i4)') nucnam,nama
    
c---- Initial Gap Parameters
          read(lin,*)
          read(lin,102) del
c
c---- Parameterset of the Lagrangian
          read(lin,*)
          read(lin,'(12x,a10)') parname	
c      read(lin,'(/12x,i4)') itpl
c      read(lin,'(12x,i4)') jpl
c      read(lin,'(12x,i4)') ippl
c      read(lin,'(12x,i4)') kpl
c
c..................................
caf ... end of reading from DIS.dat
c..................................
c
    
          close(lin)
          call nucleus(2,npro,nucnam)
c
          nneu = nama - npro
          amas = nama 
          
c..............................................................
caf ... begin of output into .OUT file
c..............................................................
c
      if (l6.ne.6) open(l6,file='dirhb.out',status='unknown')
c
      call date_and_time( date, time, zone, values )
      write(l6,'(a)')
      write(l6,'(a)') '  ******************************************  '
      write(l6,'(a)') '  *           Fortran 77 Code              *  '
      write(l6,'(a)') '  *         Spherical H.O. basis            *  '
      write(l6,'(a)') '  * Dirac-Hartree-Bogoliubov calculation   *  '
      write(l6,'(a)') '  *        with density-dependent          *  '
      write(l6,'(a)') '  * meson-exchange or point coupling force *  '
      write(l6,'(a)') '  *          and separable pairing         *  '
      write(l6,'(a)') '  * -------------------------------------- *  '
      write(l6,'(a)') '  *      Niksic, Paar, Vretenar, Ring      *  '
      write(l6,'(a,i2,a,i2,a,i4,a,a2,a,a2,a,a2,a)')
     &               '  *          ',values(3),'/',values(2),'/',
     &          values(1),'/',time(1:2),':',time(3:4),':',time(5:6),
     &          '           *'
      write(l6,'(a,a2,i4,a,i3,a,i3)')
     &              '  *       ',
     &    nucnam,nama,'  N = ',nama-npro,'  Z = ',npro
      write(6,'(a,16x,a10)') '  *',parname
      write(l6,'(a)') '  ******************************************  '

      write(6,'(a)')
      write(6,'(a)') '  ******************************************  '
      write(6,'(a)') '  *           Fortran 77 Code              *  '
      write(6,'(a)') '  *         Spherical H.O. basis            *  '
      write(6,'(a)') '  * Dirac-Hartree-Bogoliubov calculation   *  '
      write(6,'(a)') '  *        with density-dependent          *  '
      write(6,'(a)') '  * meson-exchange or point coupling force *  '
      write(6,'(a)') '  *          and separable pairing         *  '
      write(6,'(a)') '  * -------------------------------------- *  '
      write(6,'(a)') '  *      Niksic, Paar, Vretenar, Ring      *  '
      write(6,'(a,i2,a,i2,a,i4,a,a2,a,a2,a,a2,a)')
     &               '  *          ',values(3),'/',values(2),'/',
     &          values(1),'/',time(1:2),':',time(3:4),':',time(5:6),
     &          '           *'
      write(6,'(a,a2,i4,a,i3,a,i3)')
     &              '  *       ',
     &    nucnam,nama,'  N = ',nama-npro,'  Z = ',npro
      write(6,'(a,16x,a10)') '  *',parname
      write(6,'(a)') '  ******************************************  '

      write(l6,*) ' ****** BEGIN READER *******************************'
c
c---- Basisparameters:            
      write(l6,101) ' Number of oscillator shells : ',n0f,n0b
c
c---- Initialization of wavefunctions:
      write(l6,'(a,2i5)') ' Initialization inin,inink   : ',inin,inink  
c
c---- Nucleus under consideration
      write(l6,'(a,a,i4,i6,i4)') ' Nucleus: ',nucnam,nama,nneu,npro
c
c---- Initial Gap Parameters
      write(l6,103) ' Initial Gap Parameters      : ',del
c---- Parameterset of the Lagrangian
      write(l6,106) ' Parameter set               : ',parname
c.............................................................
caf ... end of output into .OUT file
c.............................................................
      tz(1) = nneu
      tz(2) = npro
c
c
  100 format(10x,2i5)
  101 format(a,2i5)
  102 format(10x,2f10.4) 
  103 format(a,2f10.4) 
  106 format(a,'   ',a10) 
  110 format(10x,4f10.4)
  111 format(9x,i1,4f10.4)

      write(l6,*) ' ****** END READER *********************************'
c-end-READER 
      end

C=======================================================================

      subroutine nucleus(is,npro,te)

C=======================================================================
C
C     is = 1 determines the symbol for a given proton number npro
c          2 determines the proton number for a given symbol te
c
C-----------------------------------------------------------------------
C
      PARAMETER (MAXZ=140)
C
      CHARACTER TE*2,T*(2*MAXZ+2)
C
      T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
      T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
      T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
      T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
      T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
      T(201:240) = 'FmMdNoLrRfHaSgNsHsMr10111213141516171819'
      T(241:280) = '2021222324252627282930313233343536373839'
      T(281:282) = '40'
c
c ... Rf is called also as Ku (kurchatovium)
c ... Ha: IUPAC calls it as dubnium (Db). J.Chem.Educ. 1997, 74, 1258
c ... Ha is called also as Db (Dubnium)
c
      if (is.eq.1) then
        if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong' 
            te = t(2*npro+1:2*npro+2)
             return
            else
c
            do np = 0,maxz
               if (te.eq.t(2*np+1:2*np+2)) then
                  npro = np
               if (npro.gt.maxz) write(6,100) TE
                   return
                endif
            enddo
c
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
         endif
c
          stop
C-END-NUCLEUS
      END

c=====================================================================c

      subroutine forces(lpr)

c=====================================================================c
c
c---- options
c     center of mass:
c        icm = 0   hb0 = hbc**2/(2*amu) ecm = 3/4*hom
c              1   hb0 = hb0*(1-1/A)    ecm = 3/4*hom
c              2   hb0 = hb0            ecm = <P**2>/2M
c
c     model-type   DD     density-dependent meson-coupling
c                  PC     point-coupling
c
c---------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      character parname*10                     ! common partyp

      common /dforce/ a_s,a_v,a_ts,a_tv,b_s,b_v,b_ts,b_tv,
     &                c_s,c_v,c_ts,c_tv,d_s,d_v,d_ts,d_tv,dsat
      common /masses/ amu,amsig,amome,amdel,amrho
      common /coupld/ ddsig,ddome,dddel,ddrho
      common /couplg/ ggsig,ggome,ggdel,ggrho,lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /partyp/ parname
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/  gl(2),gal
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN FORCES *******************************'
c---- forces
c===============================================================
       if (parname.eq.'DD-ME2') then
         amu    =  939.d0                 ! MeV
         amsig  =  550.12380d0            ! MeV
           amome  =  783.d0                 ! MeV
           amrho  =  763.d0                 ! MeV
         gsig   =   10.5396d0
         gome   =   13.0189d0
         gdel   =   zero
           grho   =    3.6836d0
         b_s    =    1.0943d0
           c_s    =    1.7057d0
         c_v    =    1.4620d0
         a_tv   =    0.5647d0
         dsat   =    0.152d0
         d_s  = one/sqrt(3.d0*c_s)
         a_s  = (one+c_s*(one+d_s)**2)/(one+b_s*(one+d_s)**2)
         d_v  = one/sqrt(3.d0*c_v)
         facs =two*a_s*(b_s-c_s)*(one-3.d0*c_s*(one+d_s)**2)/
     &                          (one+c_s*(1+d_s)**2)**3
         faco = (one-3.d0*c_v*(one+d_v)**2)/(one+c_v*(1+d_v)**2)**3
         x    = facs/(two*faco)
         fac1 = x+c_v*(one+c_v*(one+d_v)**2)
         fac2 = one+c_v*(one+d_v)**2-x*(one+d_v)**2
         b_v = fac1/fac2
         a_v =(one+c_v*(one+d_v)**2)/(one+b_v*(one+d_v)**2)

         a_ts   = zero
         b_ts   = zero
         c_ts   = zero
         d_ts   = zero

         ipc    =  0
       icm    =  2
       idd    =  2
c===============================================================
       elseif (parname.eq.'DD-PC1') then
c---------------------------------------------------------------
c        G(x) = a + (b + c*x) * exp(-d*x)
c---------------------------------------------------------------

         dsat   =  0.152d0              ! fm^-3
         amu    =  939.d0               ! MeV
c
c        scalar-isoscalar
         a_s    = -10.0462d0           ! fm^-2
         b_s    =  -9.1504d0           ! fm^-2
         c_s    =  -6.4273d0           ! fm^-2
         d_s    =  +1.3724d0
c
c        vector-isoscalar
         a_v    =  +5.9195d0           ! fm^-2
         b_v    =  +8.8637d0           ! fm^-2
         c_v    =   0.00000d0           ! fm^-2
         d_v    =  +0.6584d0
c
c        scalar-isovector
         a_ts   =   zero
         b_ts   =   zero
         c_ts   =   zero
         d_ts   =   zero
c
c        vector-isovector
         a_tv   =   0.0000d0            ! fm^-2
         b_tv   =   1.8360d0           ! fm^-2
         c_tv   =   0.0000d0            ! fm^-2
         d_tv   =   0.6403d0
c
c----- derivative terms                 ! MeV^-4
         ddsig  =  -0.8149d0
         ddome  =   zero
         dddel  =   zero
         ddrho  =   zero
c
c----------------------------------------------------
c
c        gg = 1
         ggsig = one
         ggome = one
         ggdel = one
         ggrho = one
         amsig = zero
         amome = zero
         amdel = zero
         amrho = zero

         icm    = 2
         idd    = 2
         ipc    = 1
      else
          stop 'This type of force is not defined'
      endif
c===============================================================
      amu      = amu/hbc
      amsig    = amsig/hbc
      amome    = amome/hbc
      amdel    = zmdel/hbc
      amrho    = amrho/hbc

      if (ipc.eq.0) then
         ggsig = -(gsig/(amsig+1.d-10))**2
         ggome = +(gome/(amome+1.d-10))**2
         ggdel = -(gdel/(amdel+1.d-10))**2
         ggrho = +(grho/(amrho+1.d-10))**2
         if (abs(gsig).gt.1.d-5) lmes(1) = 1
         if (abs(gome).gt.1.d-5) lmes(2) = 1
         if (abs(gdel).gt.1.d-5) lmes(3) = 1
         if (abs(grho).gt.1.d-5) lmes(4) = 1
      endif

c----------------------------------------------------------------------
c---- Separable pairing force

      gl(1) = -728.d0
      gl(2) = -728.d0
      gal   = 0.415d0
c
c---- printout of force:
      if (lpr) call pripar
c
   10 if (lpr)
     &write(l6,*) ' ****** END FORCES *********************************'
c
      return
c-end-FORCES
      end




c=====================================================================c

      subroutine pripar

c=====================================================================c
c
c     prints parameters of the Lagrangian
c
c---------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
c
      character parname*10                     ! common partyp
c
      common /partyp/ parname
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,amsig,amome,amdel,amrho
      common /pair  / del(2),spk(2),spk0(2)
      common /coupld/ ddmes(4)
      common /couplg/ ggmes(4),lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/  gl(2),gal
c
      write(l6,'(/,a,a10)') ' NEDF Parameters: ',parname
c
      if (ipc.eq.0) then
	 write(l6,'(a,f8.3)') ' AMU   = ',amu*hbc
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' msig  = ',amsig*hbc,
     &                 '  gsig = ',gsig,'  Gsig = ',ggmes(1)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mome  = ',amome*hbc,
     &                 '  gome = ',gome,'  Gome = ',ggmes(2)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mdel  = ',amdel*hbc,
     &                 '  gdel = ',gdel,'  Gdel = ',ggmes(3)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mrho  = ',amrho*hbc,
     &                 '  grho = ',grho,'  Grho = ',ggmes(4)
c
      elseif (ipc.eq.1) then
         write(l6,'(11x,a,4x,a,4x,a,4x,a)') 'scasca','scavec',
     &                                         'vecsca','vecvec'
         write(l6,'(a,4f10.6)') ' GG   = ',ggmes
         write(l6,'(a,4f10.6)') ' DD   = ',ddmes
      else
         stop 'in PRIPAR: ipc not properly defined'
      endif   ! ipc=1
      write(l6,*) ' '
      write(l6,'(a)') ' Density dependence parameters:'
      write(l6,'(11x,a,4x,a,4x,a,4x,a)') 'scasca','scavec',
     &                                   'vecsca','vecvec'
      write(l6,'(a,4f10.6)') ' a    = ',(a_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' b    = ',(b_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' c    = ',(c_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' d    = ',(d_m(m),m=1,4)
      write(l6,'(a,f10.6)') ' dsat = ',dsat

      write(l6,*) ' '
      write(l6,'(a)') ' TMR pairing: Tian,Ma,Ring, PRB 676, 44 (2009)'
      write(l6,'(a,f11.5,a)') ' TMR pairing strength    gl =: ',gl(1),
     &              ' [MeV*fm^3]'
      write(l6,'(a,f11.5,a)') ' TMR pairing width        a =: ',
     &                           sqrt(gal),' [fm]'

      return
C-end-PRIPAR
      end


c======================================================================c

       subroutine gaush(lpr)

c======================================================================c
c
c     Gauss-Hermite integration data
c     ------------------------------
c     for integration from minus to plus infinity
c                     or from 0 to infinity
c
c     ph  =  wh * exp(-xh**2)
c
c     whh=ph
c
c     \int_-\infty^+\infty  f(z) exp(-z**2) dz  =   \sum_i f(xh(i)) ph(i)
c    possible alternative
c     \int_-\infty^+\infty  f(z) dz             =   \sum_i f(xh(i)) wh(i)
c
c     fak = 1  --- for integration from 0 to +infinity
c     fak = 2  --- for integration from -infinity to +infinity
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      dimension xt(2*ngh),pt(2*ngh)

      if (lpr)
     &write(l6,*) ' ****** BEGIN GAUSH *******************************'
      call gauher(xt,pt,2*ngh)
      do i=1,ngh
         xh(i)=xt(ngh+1-i)
         ph(i)=pt(ngh+1-i)
         wh(i) = ph(i)*dexp(xh(i)*xh(i))
      enddo
      xh(0)=1.d-10
      ph(0)=1.d-10
      wh(0)=1.d-10
c
c      write(l6,100) ngh
C
      if (.not.lpr) return
      write(l6,101)
c
      write(l6,*) ' xh'
      write(l6,105) (xh(i),i=1,ngh)
      if (ngh.eq.8)  write(l6,102) (xh(i),i=1,ngh)
      if (ngh.eq.12) write(l6,103) (xh(i),i=1,ngh)
      if (ngh.eq.16) write(l6,104) (xh(i),i=1,ngh)
c
      write(l6,*) ' wh'
      write(l6,105) (wh(i),i=1,ngh)
      if (ngh.eq.8)  write(l6,102) (wh(i),i=1,ngh)
      if (ngh.eq.12) write(l6,103) (wh(i),i=1,ngh)
      if (ngh.eq.16) write(l6,104) (wh(i),i=1,ngh)
c
      write(l6,*) ' ph'
      write(l6,105) (ph(i),i=1,ngh)
      if (ngh.eq.8)  write(l6,102) (ph(i),i=1,ngh)
      if (ngh.eq.12) write(l6,103) (ph(i),i=1,ngh)
      if (ngh.eq.16) write(l6,104) (ph(i),i=1,ngh)
c
  100 format('  GAUSH:  G-H-Integration  ngh =',i3)
  101 format(1x,36(1h-))
  102 format(2(3e19.11,/),2(e19.11))
  103 format(4(3e19.11,/))
  104 format(5(3e19.11,/),e19.11,/)
  105 format(3e19.11)
c
      if (lpr)
     &write(l6,*) ' ****** END GAUSH *********************************'
      return
c-end-GAUSH
      end


c=======================================================================
c
      SUBROUTINE gauher(x,w,n)
c     Calculates gaussian mesh points
c
c=======================================================================
c
      INTEGER n,MAXIT
      DOUBLE PRECISION w(n),x(n)
      DOUBLE PRECISION EPS,PIM4
      PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=10)
      INTEGER i,its,j,m
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      m=(n+1)/2
      do 13 i=1,m
        if(i.eq.1)then
          z=sqrt(dble(2*n+1))-1.85575d0*(2*n+1)**(-.16667d0)
        else if(i.eq.2)then
          z=z-1.14d0*n**.426d0/z
        else if (i.eq.3)then
          z=1.86d0*z-.86d0*x(1)
        else if (i.eq.4)then
          z=1.91d0*z-.91d0*x(2)
        else
          z=2.d0*z-x(i-2)
        endif
        do 12 its=1,MAXIT
          p1=PIM4
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
11        continue
          pp=sqrt(2.d0*n)*p2
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        stop 'too many iterations in gauher'
1       x(i)=z
        x(n+1-i)=-z
        w(i)=2.d0/(pp*pp)
        w(n+1-i)=w(i)
13    continue
      return
      END


c======================================================================c

      subroutine prep(lpr)

c======================================================================c
c
c     preparations
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

c
      character tp*1,tis*1,tit*8,tl*1                    ! common textex
      character nucnam*2                                 ! common nucnuc
      character tb*5
      logical lpr

c
      common /baspar/ hom,hb0,b0
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c
      write(l6,*) ' ****** BEGIN PREP *********************************'
c
c---- basis parameters
      hb0 = hbc/(two*amu)
      hom = 41.0*amas**(-third)
      if (icm.eq.1) hb0 = hb0*(one - one/amas)
      b0 = sqrt(two*hb0/hom)
c
      write(l6,*) ' '
      write(l6,100) ' hom =                         ',hom
      write(l6,100) ' hb0 =                         ',hb0
      write(l6,100) ' b0  =                         ',b0
c
      do ih = 0,ngh
         rb(ih) = xh(ih)*b0
c        metric element for three-dimensional integration
         wdcor(ih) = b0**3 * 4*pi * xh(ih)**2 * wh(ih)
      enddo
c
c---- printout pairing:
c      write(l6,*)
  100 format(a,4f11.6)
c      write(l6,100) ' Initial Gap   = ',del
c
   10 if (itx.eq.1) icou = 0
      if (icou.eq.0) write(l6,100) ' without Coulomb force'
      if (icou.eq.1) write(l6,100) ' with Coulomb force'
      if (icou.eq.2) write(l6,100) ' with Coulomb force with exchange'
      if(icm.eq.2) write(l6,*) 'With microscopic c.m. correction'
  101 format(a,i4)


      write(l6,*) ' ****** END PREP ***********************************'
      return
c-end PREP
      end

c======================================================================c

      subroutine inout(is,lpr)

c======================================================================c
c
c     is = 1: reads meson-fields from tape
c          2: writes meson-fields  to tape
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
c
      common /fermi / ala(2),tz(2)
      common /initia/ inin,inink
      common /potpot/ vps(0:ngh,2),vms(0:ngh,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (is.eq.1.and.inin.ne.0) return
c
      if(lpr)
     &write(l6,*) ' ****** BEGIN INOUT ********************************'
c
c---- reading of meson fields from tape:
c-------------------------------------
      if (is.eq.1) then
         open(lwin,file='dirhb.wel',status='old',form='unformatted')
         read(lwin) ng0
         if (ng0.ne.ngh) stop 'in INOUT: ngh wrong'
         read(lwin) ala

c------- reading of the potentials
         read(lwin) vms

         read(lwin) vps
         close(lwin)

         write(l6,*) ' potentials read from tape ','dirhb.wel'
         write(l6,*) ' vms(0,1) = ',vms
         write(l6,*) ' vps(0,1) = ',vps 
c---- writing to tape:
      else
         open(lwou,file='dirhb.wel',status='unknown',form='unformatted')
         write(lwou) ngh
         write(lwou) ala
         write(lwou) vms
         write(lwou) vps
         close(lwou)

         if(lpr) write(l6,*) ' potentials written to tape dirhb.wel'
      endif
c
      if(lpr)
     &write(l6,*) ' ****** END INOUT **********************************'
      return
c-end-INOUT
      end

c======================================================================c

      subroutine start(lpr)

c======================================================================c
c
c     initializes potentials
c     inin = 0:   reads fields from tape lwin
c            1:   saxon-woods
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      character*2 nucnam
c
      dimension vso(2),r0v(2),av(2),rso(2),aso(2)
c
      common /baspar/ hom,hb0,b0
      common /coulmb/ cou(0:ngh),drvp(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /initia/ inin,inink
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(0:ngh,2),vms(0:ngh,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c=======================================================================
c     Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
      data v0/-71.28/,akv/0.4616/
      data r0v/1.2334,1.2496/,av/0.615,0.6124/
      data vso/11.1175,8.9698/
      data rso/1.1443,1.1401/,aso/0.6476,0.6469/
c---- potentials are read in INOUT
      if (inin.eq.0) return
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN START ********************************'
c
c
      if (lpr) then
         write(l6,'(a,f10.4)') ' v0     = ',v0
         write(l6,'(a,f10.4)') ' kappa  = ',akv
         write(l6,'(a,2f10.4)') ' lambda = ',vso
         write(l6,'(a,2f10.4)') ' r0     = ',r0v
         write(l6,'(a,2f10.4)') ' a      = ',av
         write(l6,'(a,2f10.4)') ' r0-so  = ',rso
         write(l6,'(a,2f10.4)') ' a-so   = ',aso
      endif
      do ih = 0,ngh
         r = rb(ih)
c
c------- Woods-Saxon potential
         do it = 1,itx
            ita = 3-it
            rav = r0v(it)*amas**third
            rao = rso(it)*amas**third
            vp  = v0*(one - akv*(npr(it)-npr(ita))/amas)
            vls = vp * vso(it)
c
            argv = (r - rav)/av(it)
            if (argv.le.65.d0) then
               u = vp/(one + exp(argv))
            else
               u = zero
            endif
            argo = (r - rao)/aso(it)
            if (argo.le.65.d0) then
               w = -vls/(one + exp(argo))
            else
               w = zero
            endif
            vps(ih,it) = u
            vms(ih,it) = w
         enddo   ! it
           if (itx.eq.1) then
              vms(ih,2) = vms(ih,1)
              vps(ih,2) = vps(ih,1)
         endif   ! itx=1
c
c------- Coulomb potential
         cou(ih) = zero
         if (icou.ne.0) then
            rc = r0v(2)*amas**third
            if (r.lt.rc) then
               c = half*(3/rc - r*r/rc**3)
            else
               c = one/r
            endif
            cou(ih)   = c*npr(2)/alphi
              vps(ih,2) = vps(ih,2) + cou(ih)*hbc
              vms(ih,2) = vms(ih,2) + cou(ih)*hbc
           endif   ! icou.ne.0
      enddo   ! ih
      if (lpr)
     &write(l6,'(/,a)') ' Initial potentials of Saxon-Woods shape '

      if (lpr)
     &write(l6,*) ' ****** END START **********************************'
      return
c-end START
      end


c======================================================================c

      subroutine base(lpr)

c======================================================================c
c
c     determines the basis in spherical oscillators for Dirac solution
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      common /basnnn/ n0f,n0b
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2),ic(nbx),iz(nbx)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)
      common /bosqua/ no
      common /gfviv / iv(0:igfv)
      common /sdimos/ nrm,nlm,nrbm
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /vvvikf/ mv,ipos(nbx),nib(mvx),nni(2,mvx)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN BASE *********************************'

c
      !write(l6,*) ' n0fx0', n0fx0 
      n0fx0 = n0fx
      if (n0f.gt.n0fx0) stop ' in BASE: n0f too large'
      if (is.eq.2.and.n0fx0.lt.20) stop 'in BASE: n0fx < 20'
c
      nrm = 1 + (n0f+1)/2
      nlm = n0f+1
      !write(l6,*)'nrm', nrm
      !write(l6,*)'nlm',nlm
      if (nrm.gt.ndx) stop ' in BASE: ndx too small'
      if (nrm.gt.nrx) stop ' in BASE: nrx too small '
      if (nlm.gt.nlx) stop ' in BASE: nlx too small '
c
c----------------------------------------------------------------------c
c     construction of the different kappa-blocks
c----------------------------------------------------------------------c
      kapmax = n0f + 1
      ib  = 0
      il  = 0
      ik  = 0
      ik1 = 0
      ntz = 0
      nfm = 0
      ngm = 0
c
c
      do jj = 1,n0f+1
            do ipp = 1,2
            ibb = 2*(jj-1) + ipp
            if (mod(jj,2).eq.ipp-1) then
                  ik1 = +jj
                  l = jj
                  la = jj-1
            else
                  ik1 = -jj
                  l = jj-1
                  la = jj
            endif
            ijb(ibb) = jj
            !write(l6,*)ijb(ibb)
            ilb(ibb,1) = l
            !write(l6,*)ilb(ibb,1)
            ilb(ibb,2) = la
            !write(l6,*)ilb(ibb,2)
            ipb(ibb) = ipp
            !write(l6,*)ipb(ibb)
            ikb(ibb) = ik1
            !write(l6,*)ikb(ibb)
            enddo	! ip
      enddo	! j
      !write(l6,*)ijb(ibb),ilb(ibb,1),ilb(ibb,2),ipb(ibb),ikb(ibb)
      
      
      do j = 1,kapmax
         do ipar = 1,-1,-2
              lf = j - (1-ipar*iv(j))/2
            kappa = j * (2*(lf-j)+1)
              lg    = 2*j - lf -1
            !write(l6,*)'lf',ipar, lf ,kappa
            !write(l6,*)'lg',lg
            if (lf.le.n0f) then
               ip = 1 + mod(lf,2)
               nf = (n0f-lf)/2 + 1
               ng = (n0f+1-lg)/2 + 1
               ib = ib + 1
               !write(l6,*)'nf :',nf
               if (ib.gt.nbx) stop ' in BASE: nbx too small'
               write(tb(ib),'(a1,i2,2h/2)') tl(lf),j+j-1
c
                 kb(ib)    = kappa
                 mb(ib)    = 2*iabs(kappa)
               id(ib,1)  = nf
               id(ib,2)  = ng
               if (ib.eq.0) then
                    ic(ib) = 0
               else
                    ic(ib) = ic(ib - 1) + id(ib-1,1)
               endif
               if (ib.eq.0) then
                  iz(ib) = 0
             else
                  iz(ib) = iz(ib - 1) + id(ib-1,1)
             endif
               !ic(0) = 0 

                !ic(ib) = ic(ib - 1) +  id(ib,2)
                 ia(ib,1)  = il

                 ia(ib,2)  = il + id(ib,1)
c
               ntz = ntz + 2*j
c
               do ifg = 1,2
                  ne = id(ib,ifg)
                  l  = lfgkap(kappa,ifg)
                  do n = 1,ne
                     il = il + 1
                     if (il .gt.ntx) stop ' in BASE: ntx too small'
                     nr(il) = n
                     nl(il) = l
                     nj(il) = j
                     kk(il) = kappa
                         np(il) = ifg
                     write(tt(il),'(i2,a1,i3,2h/2)') n,tl(l),j+j-1
                     nn     = 2*(n-1)+l
                  enddo   ! n
               enddo   ! ifg
                 ik = ik + max(nf,ng)
                 nfm = max(nfm,nf)
                 ngm = max(ngm,ng)
            endif
c
   10    enddo   ! ivkap
      enddo   ! j
      nb  = ib
      nt  = il
      nk  = ik
      no  = n0b/2 + 1
      if (nk.gt.nkx)   stop 'in BASE: nkx too small'
      if (nfm.gt.nfx)  stop 'in BASE: nfx too small '
      if (ngm.gt.ngx)  stop 'in BASE: ngx too small '
c
c
c
c
c----------------------------------------------------------------------c
c     Printout
c----------------------------------------------------------------------c
      if (lpr) then
         do i = 1,nt
            nn     =  2*(nr(i)-1) + nl(i)
            write(l6,'(i3,i3,1x,a,i4)') i,nn,tt(i),np(i)
         enddo
         do ib = 1,nb
            nf = id(ib,1)
            ng = id(ib,2)
            nh = nf + ng
            j  = iabs(kb(ib))
            write(l6,'(/,a,5i4)') tb(ib),id(ib,1),id(ib,2),
     &                                   ic(ib),ia(ib,1),ia(ib,2)
            do i = ia(ib,1)+1,ia(ib,1)+nh

               nn     =  2*(nr(i)-1) + nl(i)
               write(l6,102) i,'   NN = ',nn,
     &         '   nr = ',nr(i),'   l =',nl(i),'   j =',j+j-1,tt(i)
               if (i.eq.ia(ib,1)+nf) write(l6,'(3x,61(1h-))')

            enddo   ! i
         enddo   ! ib
c
         write(l6,'(/,a,2i4)') ' Number of blocks: nb  = ',nb,nbx
         write(l6,100) ' Number of levels  nt  = ',nt,ntx
         write(l6,100) ' Number of levels  nk  = ',nk,nkx
         write(l6,100) ' Maximal n:        nrm = ',nrm,nrx
         write(l6,100) ' Maximal l:        nlm = ',nlm,nlx
         write(l6,100) ' Maximal nf            = ',nfm,nfx
         write(l6,100) ' Maximal ng            = ',ngm,ngx
      endif
c
c
c----------------------------------------------------------------------c
c     Construction of two-body pairs (i1,i2)
c----------------------------------------------------------------------c
c---- only f-pairs
      if (lpr)
     &write(l6,*) ' ****** PAIRING ***********************************'
      il = 0
      do ib = 1,nb
         ipos(ib) = il
         nf = id(ib,1)
         do n2 =  1,nf
         do n1 = n2,nf
            il = il + 1
            nib(il)   = ib
            nni(1,il) = n1
            nni(2,il) = n2
         enddo  ! n1
         enddo  ! n2
      enddo  ! ib
      mv = il
      if (lpr) write(l6,100) ' number of pairs f  mv = ',mv,mvx
      if (mv.gt.mvx) then
       write(6,*) 'mv =',mv,' mvx = ',mvx
       stop ' in BASE: mvx too small'
      endif
c
  100 format(a,2i6)
  102 format(i4,a,i2,a,i2,a,i2,a,i2,3h/2 ,a)
c
      if (lpr)
     &write(l6,*) ' ****** END BASE ***********************************'
c
c-end-BASE
      end

c======================================================================c
      integer function lfgkap(kappa,is)
c======================================================================c
      if (is.eq.1) then
         if (kappa.gt.0) then
            lfgkap = kappa
         else
            lfgkap = - kappa - 1
         endif
      else
         if (kappa.gt.0) then
            lfgkap = kappa - 1
         else
            lfgkap = - kappa
         endif
      endif
      return
c-end-LFGKAP
      end

c======================================================================c

      subroutine dinout(is,lpr)

c======================================================================c
c
c     IS = 1 : for ININK = 0  reads pairing field Delta from tape  
c                  ININK = 1  calculates pairing field Delta 
c     IS = 2 : writes pairing field Delta to tape                  
c
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tt*8                                            ! bloqua
      character tb*5                                            ! blokap
c
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
      common /deldel/ de(nhhx,nb2x)
      common /initia/ inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /vvvikf/ mv,ipos(nbx),nib(mvx),nni(2,mvx)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN DINOUT ****************************'
c
c
c==== initialize or read pairing potential
      if (is.eq.1) then
         do ib = 1,nb2x
            do i = 1,nhhx
	       de(i,ib) = zero
	    enddo   ! i
         enddo   ! ib
c        
c---- reads pairing potential
      if (inink.eq.0) then
c
         open(laka,file='dirhb.del',form='unformatted',status='unknown')
         read(laka) mv0
         if (mv0.ne.mv) stop 'in DINOUT: mv wrong'
         do it = 1,itx
            do ib = 1,nb
	           nf = id(ib,1)
	           ng = id(ib,2)
	           nh = nf + ng
	           m  = ib + (it-1)*nbx
	           do n2 = 1,nf
 	              read(laka) (de(n1+(n2-1)*nh,m),n1=n2,nf)
	              do n1 = n2,nf
                     de(n2+(n1-1)*nh,m) = de(n1+(n2-1)*nh,m)
                  enddo   ! n1
               enddo   ! n2
            enddo   ! ib
         enddo   ! it
         close(laka)
c     
c     
c---- initial constant pairing field
      elseif (inink.eq.1) then
         do it = 1,itx
            do ib = 1,nb
               nf = id(ib,1)
               ng = id(ib,2)
               nh = nf + ng
               m  = ib + (it-1)*nbx
               do n = 1,nf
                  de(n+(n-1)*nh,m) = del(it)
               enddo   ! n
            enddo   ! ib
         enddo   ! it
      else
         stop 'in DINOUT: inink wrong'
      endif   ! inink
      if (lpr) write(l6,*) ' Pairing field has been calculated'
c
c==== writing of the pairing potential
      elseif (is.eq.2) then
c
c
         open(laka,file='dirhb.del',form='unformatted',status='unknown')
         write(laka) mv
         do it = 1,itx
            do ib = 1,nb
	           nf = id(ib,1)
	           ng = id(ib,2)
	           nh = nf + ng
	           m  = ib + (it-1)*nbx
	           do n2 =  1,nf
 	              write(laka) (de(n1+(n2-1)*nh,m),n1=n2,nf)
               enddo   ! n2
            enddo   ! ib
         enddo   ! it
c
      else
         stop 'in DINOUT: is wrong'
      endif   ! is
c
      if (lpr) then
         do it = 1,itx
            do ib = 1,nb
	           nf = id(ib,1) 
	           ng = id(ib,2)
	           nh = nf + ng
	           m  = ib + (it-1)*nbx
	           k0 = ia(ib,1)+1
	        enddo   ! ib
	     enddo   ! it
         write(l6,*) ' ****** END DINOUT ****************************'
      endif
      return
c-end-DINOUT
      end

c======================================================================c

      subroutine gaupol(lpr)

c======================================================================c
c
c     calculates the radial functions for the spherical oscillator
c
c     the wave function phi(nlj) of the spherical oscillator are: 
c
c     phi(r,Omega) = b^(-3/2) * R_nl(r) * Y_ljm(Omega) 
c     
c     R_nl(r) = N_nl * r^l  * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
c
c     N_nl    = sqrt(2 * (n-1)!/(n+l-1/2)!)     and    x=r/b
c
c     the contribution to the density from the shell j is
c
c     rho_j(r)= 1/(4*pi*b0**3) * (2*j+1) * R_nl(r)^2
c
c     the radial function at meshpoint xh(ih) is stored in RNL(n,l,ih)
c     in the following way: RNL is normalized in such way that the
c     norm integral reads
c
c     \int d^3r |phi(r)|^2 = 1 = \sum_i RNL(n,l,i)**2
c
c     this means, that RNL contains the following factors:
c
c     a)  the radial part of the wavefunction r * R_nl(r)
c     b)  the length units factor  b ** (3/2)
c     c)  the Gaussian weight sqrt( WH(i) ): 
c         \inf_0^inf f(x) dx = \sum_i f(x_i) * WH(i)
c
c     having RNL(n,l,i) we get the radial wavefunction:
c
c     R_nl(r) =  RNL(n,l,i) / ( x_i * sqrt(WH(i)) )  
c
c     and the density contribution from the shell j
c
c     rho_j(r) = (2*j+1) * RNL(n,l,i)**2 / ( 4 * pi x_i**2 * WH(i) * b**3)   
c
c----------------------------------------------------------------------c
c

c     RNL1 contains the radial derivatives in the following form:
c
c     d/dr R_nl(r) = 1/b * RNL1(n,l,i) / (x_i * sqrt(WH(i) )
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tt*8                                            ! bloqua
c
      common /bosqua/ no
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /radosc/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /bospol/ rnb(1:nox,0:ngh)
      common /sdimos/ nrm,nlm,nrbm
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN GAUPOL ************************'
c
c     f = 2/pi**0.25
      f = 2*wgi(0)
c
      do ih = 0,ngh
         r  = xh(ih)
         rr = r*r
         ri = one/r 
         fe = f*exp(-half*rr)
c
c------------------------------------
c        the functions rnl,rnl1 contain already the measure
         u1 = fe*sqrt(wh(ih)*rr)
c------------------------------------
c
c------- basis for fermions
         do l = 0,nlm
            rnl(1,l,ih)  = u1
            rnl(2,l,ih)  = u1*(l+1.5d0-rr)*shi(l+1)
            u1           = u1*r*shi(l+1)
            rnl1(1,l,ih) =    (l-rr)*rnl(1,l,ih)*ri
            rnl1(2,l,ih) = ((2+l-rr)*rnl(2,l,ih) - 
     &                       2*sqh(l+1)*rnl(1,l,ih))*ri
c
            do n = 3,nrm
               rnl(n,l,ih)  = ((2*n+l-2.5d0-rr)*rnl(n-1,l,ih) -
     &           sq(n-2)*sqh(n-2+l)*rnl(n-2,l,ih))*sqi(n-1)*shi(n-1+l)
               rnl1(n,l,ih) = ((2*n+l-2-rr)*rnl(n,l,ih) -
     &           2*sq(n-1)*sqh(n-1+l)*rnl(n-1,l,ih))*ri
            enddo
         enddo
         !write(l6,*)'base for fermions',rnl, rnl1
c
c------- basis for bosons
         rnb(1,ih)  = fe
         rnb(2,ih)  = fe*(1.5d0-rr)*shi(1)
         do n = 3,no
            rnb(n,ih) = ((2*n-2.5d0-rr)*rnb(n-1,ih) -
     &           sq(n-2)*sqh(n-2)*rnb(n-2,ih))*sqi(n-1)*shi(n-1)
         enddo
c
      enddo 
      !write(l6,*)'base for fermions',rnl, rnl1
c
c
c---- Test of orthogonality
      if (lpr) then
         do 40 l = 0,nlm
            write(l6,'(/,80(1h*))')
            do 41 n = 1,nrm
               write(l6,'(a,2i3)') 
     &         ' Radial function and derivative for n,l =',n,l
               ix = 5
               write(l6,'(5f15.8)') (rnl(n,l,ih),ih=1,ix)
               write(l6,'(5f15.8)') (rnl1(n,l,ih),ih=1,ix)
   41       continue
            do 50 n2 = 1,nrm
            do 50 n1 = n2,nrm 
               s1 = zero
               s2 = zero
               s3 = zero
               sb = zero
               do ih = 0,ngh
                  rr = xh(ih)**2
                  s0 = rnl(n1,l,ih)*rnl(n2,l,ih)
                  s1 = s1 + s0
                  s2 = s2 + rr*s0
                  s3 = s3 + (rnl1(n1,l,ih)*rnl1(n2,l,ih)
     &                       + rnl1(n1,l,ih)*rnl(n2,l,ih)/xh(ih)
     &                       + rnl(n1,l,ih)*rnl1(n2,l,ih)/xh(ih)
     &                       + s0*(1+l*(l+1))/rr)
		  if (l.eq.0) then
		     sb = sb + rnb(n1,ih)*rnb(n2,ih)*rr*wh(ih)
                  endif
               enddo
               write(l6,'(a,2i3,4f12.8)') 
     &                  ' RNL(n,l) test ',n1,n2,s1,s2,s3,sb
   50       continue
   40    continue
         write(l6,'(/,80(1h*))')
      endif
c
      if (lpr)
     &write(l6,*) ' ****** END GAUPOL *************************'
      return
c-end-GAUPOL
      end

c======================================================================c

      subroutine singf(lpr)

c======================================================================c
c
c     calculates single particle matrix elements for Fermions       
c     in the spherical oscillator basis
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tb*5                                            ! blokap
c
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /single/ sp(nfgx,nbx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN SINGF ********************************'
c
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
c
c        SIGMA*P
c----------------
         call sigp(nf,ng,ib,sp(1,ib),lpr)
c
      enddo
C
      if (lpr)
     &write(l6,*) ' ****** END SINGF **********************************'
      return
c-end-SINGF
      end  
c=====================================================================c

      subroutine sigp(nf,ng,ib,aa,lpr)

c=====================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      dimension aa(ng,nf)
c
      common /baspar/ hom,hb0,b0
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /radosc/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      kappa = kb(ib)
      lf    = lfkap(kappa)
      lg    = lgkap(kappa)
c
      kkk = - kappa - 1
      do n2 = 1,nf
      do n1 = 1,ng
         s = zero
         do ih = 1,ngh
            s = s + rnl(n1,lg,ih) * 
     &              ( - rnl1(n2,lf,ih) + kkk*rnl(n2,lf,ih)/xh(ih))    
         enddo
         aa(n1,n2) = s
      enddo
      enddo
      !write(l6,*)'s ', kkk , 'ih', ih, 'ss', s
      return
c-end-SIGP
      end

c======================================================================c
      integer function lfkap(kappa)
c======================================================================c
      if (kappa.gt.0) then
         lfkap = kappa
      else
         lfkap = - kappa - 1
      endif
      return
c-end-LFKAP
      end
c======================================================================c
      integer function lgkap(kappa)
c======================================================================c
      if (kappa.gt.0) then
         lgkap = kappa - 1
      else
         lgkap = - kappa
      endif
      return
c-end-LGKAP
      end

c======================================================================c

      subroutine singd(lpr)

c======================================================================c
c
c     calculates single particle matrix V_nn' for separable-pairing
c     in the spherical oscillator basis
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      dimension oscn(0:n0fx) 
c
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      common /baspar/ hom,hb0,b0
      common /physco/ hbc,alphi,r0
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),
     &                tt(ntx)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /vvvsep/ vn(mvx,0:n0fx)
      common /tmrpar/  gl(2),gal
      common /gfviv / iv(0:igfv)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN SINGD ********************************'
c	
      call setwig
      
      call vn0(oscn)
      do i = 1,mvx
	     do nn = 0,n0fx
		    vn(i,nn) = zero
	     enddo    ! nn
	  enddo    ! i
      il = 0
	  do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
         m  = ib + (it-1)*nbx
         kappa = kb(ib)
         j     = iabs(kappa)
         l     = lfkap(kappa)
         do n2 = 1,nf
            do n1 = n2,nf
               il = il + 1
               do nn = 0,n1-1+n2-1+l
                  nx = n1-1+n2-1+l-nn
                  vn(il,nn) = iv(l)*talmos(n1-1,l,n2-1,l,nn,
     &                        0,nx,0,0)*oscn(nx)/sqrt(two*l+one)
               enddo                  ! nn
               
            enddo                     ! n1
         enddo                     ! n2 
      enddo                     ! ib
      if (lpr)
     &write(l6,*) ' ****** END SINGD **********************************'
      return
c-end-SINGD
      end 

c======================================================================c

      subroutine setwig

c======================================================================c
c
c     computes and stores in an efficient way couplings for Moshinsky
c
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
c
      include 'dirhb.par'
c
      integer*2 locs
c
      common /gfvsqh/ sqh(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /wigwig/ wig(jdim),locs(0:lx,0:lx,0:lx)
c
      r2 = sqh(0)
      do 10 l1 = 0,lx
      do 10 l2 = 0,lx
      do 10 l3 = 0,lx
   10    locs(l1,l2,l3) = jdim

      ico = 0
      do 20 l1 = 0,lx
      do 20 l2 = l1,lx
         kmin = l2+mod(l1,2)
         ktop = min0(lx,l1+l2)
         if (kmin.gt.ktop) goto 20
         do l3 = kmin,ktop,2
            ico = ico+1
            if (ico.gt.jdim)  stop ' in SETWIG: jdim too small'
            ip = (l1+l2+l3)/2
            wig(ico) = wg(ip-l1)*wfi(ip-l1)*wg(ip-l2)*wfi(ip-l2)*
     &                 wg(ip-l3)*wfi(ip-l3)*wf(ip)*wgi(ip+1)*r2
            locs(l1,l2,l3) = ico
            locs(l1,l3,l2) = ico
            locs(l2,l3,l1) = ico
            locs(l2,l1,l3) = ico
            locs(l3,l1,l2) = ico
            locs(l3,l2,l1) = ico
         enddo
   20 continue
      wig(jdim) = 0.0d0
c
      return
c-end-SETWIG
      end
c======================================================================c

      subroutine vn0(oscn)

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
c
      logical lpr
c
      dimension gll(4),oscn(0:n0fx)
      character tk*8                                            ! blolev
      character tp*1,tis*1,tit*8,tl*1                           ! textex
c
      common /baspar/ hom,hb0,b0
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /wavefg/ fg(nhx,nkx,4)
      common /gfvfak/ fak(0:igfv) 
      common /gfviv / iv(0:igfv) 
      common /gfvwfi/ wfi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /tmrpar/  gl(2),gal
c
	  a=sqrt(gal)
	  f1=one/sqrt(2*pi*a)
	  f1=f1**3
	  f2=a**3/b0**1.5
	  fac=f1*f2
c
      do i = 0,n0fx
         f1 = sqrt(2.d0)/pi/2.d0**0.75d0
         f2 = (one/b0)**1.5d0
         f3 = wg(i+1)
         f4 = wfi(i)
         f5 = (b0**2/(b0**2+a**2))**1.5d0
         f6 = ((b0**2-a**2)/(b0**2+a**2))**i
         s2 = f1*f2*f3*f4*f5*f6
	     oscn(i) = s2 
      enddo    ! i
      
      return      
c
      end
c======================================================================c
c
      function talmos(m1,k1,m2,k2,m3,k3,m4,k4,ilam)
c
c======================================================================c
C
c
C     Talmi - Moshinsky Bracket:( lambda n1 l1 n2 l2 n3 l3 n4 l4 )
C
c     < n1 l1, n2 l2, lambda | n3 l3, n4 l4, lambda >
C
c     see T.A. Brody & M. Moshinsky, "tablas de parentesis de
c     transformacion" (monografias del instituto de fisica, Mexico
c     1960)
c
c     zero is returned for cases that violate the triangle rules
c     energy consevation.
c
c     radial quantum number start from zero:   n1 = 0,1,2,....
c
c     the coefficients are calculated here in the phase convention
c     of Baranger. 
c     <n1 l1 n2 l2 n3 l3 n4 l4, lambda>(Brody-Moszkowsky)
c     = (-)^(b3+n3-lambda) * <n1 l1 n2 l2 n4 l4 n3 l3, lambda>(Baranger)
c       
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical testj
      integer*2 locs
c
      common /gfviv / iv(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvgmi/ gmi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /wigwig/ wig(jdim),locs(0:lx,0:lx,0:lx)
c
c     itestj is true if any triangle condition is violated
      testj(l1,l2,l3) = l1+l2.lt.l3 .or. iabs(l1-l2).gt.l3
c
      talmos = zero
c
c     check triangular rules and energy
      if (testj(k1,k2,ilam).or.testj(k3,k4,ilam)) return
      nn1 = 2*m1+k1
      nn2 = 2*m2+k2
      nn3 = 2*m3+k3
      nn4 = 2*m4+k4
      ichi = nn1+nn2
      if (ichi.ne.nn3+nn4) return
c
c---- reordering of indices
      if (min0(k1,k2).lt.min0(k3,k4)) then
         if (k1.gt.k2) then
            n3 = m2
            l3 = k2
            n4 = m1
            l4 = k1
            if (nn3.gt.nn4) then
               n1 = m4
               l1 = k4
               n2 = m3
               l2 = k3
               iph = iv(l1+l4)
            else
               n1 = m3
               l1 = k3
               n2 = m4
               l2 = k4
               iph = iv(l1+ilam)
            endif
         else
            n3 = m1
            l3 = k1
            n4 = m2
            l4 = k2
            if (nn3.gt.nn4) then
               n1 = m4
               l1 = k4
               n2 = m3
               l2 = k3
               iph = iv(l3+ilam)
            else
               n1 = m3
               l1 = k3
               n2 = m4
               l2 = k4
               iph = +1
            endif
         endif
      else
         if (k3.gt.k4) then
            n3 = m4
            l3 = k4
            n4 = m3
            l4 = k3
            if (nn1.gt.nn2) then
               n1 = m2
               l1 = k2
               n2 = m1
               l2 = k1
               iph = iv(l2+l3)
            else
               n1 = m1
               l1 = k1
               n2 = m2
               l2 = k2
               iph = iv(l1+ilam)
            endif
         else
            n3 = m3
            l3 = k3
            n4 = m4
            l4 = k4
            if (nn1.gt.nn2) then
               n1 = m2
               l1 = k2
               n2 = m1
               l2 = k1
               iph = iv(l3+ilam)
            else
               n1 = m1
               l1 = k1
               n2 = m2
               l2 = k2
               iph = +1
            endif
         endif
      endif
c
c---------------------------------------------------------------
      nn1 = 2*n1+l1
      nn2 = 2*n2+l2
      nn3 = 2*n3+l3
      nn4 = 2*n4+l4
      prout = sqh(l1)*wf(n1)*wg(n1+l1+1)*sqh(l2)*wf(n2)*wg(n2+l2+1)*
     &        sqh(l3)*wf(n3)*wg(n3+l3+1)*sqh(l4)*wf(n4)*wg(n4+l4+1)
      xf = sqh(0)**ichi
      iat = min0(nn1,nn3)
      if (l1.gt.lx.or.l2.gt.lx) goto 910
      if (iat.gt.lx)  goto 910
c
c---- main loop
      t = zero
      do 10 la = 0,iat
         lbl = iabs(la-l3)
         ibl = lbl+mod(la+lbl+l3,2)
         lbt = min0(nn2,nn3-la,la+l3)
         ibt = lbt-mod(la+lbt+l3,2)
         if (ibl.gt.ibt) goto 10
         if (ibt.gt.lx)  goto 910
c-------------
         do 20 lb = ibl,ibt,2
            lcl = iabs(la-l1)
            icl = lcl+mod(la+lcl+l1,2)
            lct = min0(nn1-la,nn4,la+l1)
            ict = lct-mod(la+lct+l1,2)
            if (icl.gt.ict) goto 20
            lo  = locs(la,lb,l3)
            c3  = wig(lo)
            if (ict.gt.lx) goto 910
c----------------
            do 30 lc = icl,ict,2
               ldl = max0(iabs(lb-l2),iabs(lc-l4))
               idl = ldl+mod(lb+ldl+l2,2)
               ldt = min0(nn2-lb,nn4-lc,lb+l2,lc+l4)
               idt = ldt-mod(lb+ldt+l2,2)
               if (idl.gt.idt) goto 30
               lo = locs(la,lc,l1)
               c1 = wig(lo)
               if (idt.gt.lx)  goto 910
c-------------------
               do 40 ld = idl,idt,2
                  mm1 = nn4-lc-ld
                  mm2 = nn2-lb-ld
                  ndt = min0(mm1/2,mm2/2)
                  if (ndt.lt.0) goto 40
                  iphd = iv(ld)
                  lo = locs(lb,ld,l2)
                  c2 = wig(lo)
                  lo = locs(lc,ld,l4)
                  c4 = wig(lo)
                  if (min0(ilam,la).gt.l3) then
                     qj = sninj0(l1,lc,la,l2,ld,lb,ilam,l4,l3)
                     qj = qj*iv(l1+l2+ilam)
                  else
                     qj = sninj0(la,lc,l1,lb,ld,l2,l3,l4,ilam)
                  endif
c----------------------
                  do 50 nd = 0,ndt
                     nnd = 2*nd+ld
                     na = (nn3-nn2+nnd-la)/2
                     if (na.lt.0) goto 50 
                     nb = (nn2-nnd-lb)/2
                     if (nb.lt.0) goto 50 
                     nc = (nn4-nnd-lc)/2
                     if (nc.lt.0) goto 50 
                     ll = (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1)
                     gg = fi(na)*gmi(na+la+1)*fi(nb)*gmi(nb+lb+1)*
     &                    fi(nc)*gmi(nc+lc+1)*fi(nd)*gmi(nd+ld+1)
                     t  = t + iphd*ll*c1*c2*c3*c4*qj*gg*prout
   50             continue
   40          continue
   30       continue
   20    continue
   10 continue
      talmos = iph*t*xf/pi
      if (dabs(talmos).lt.1.0d-10) talmos = zero
      return
c
 900  stop  'in TALMOS: n_j too large '
 910  stop ' in TALMOS: lx  too small '
c-end-TALMOS
      end

C=======================================================================

      function sninj0(j11,j12,j13,j21,j22,j23,j31,j32,j33)

C=======================================================================
C
C     calculates 9j-Symbol for integer j-values
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c 
      sninj0 = 0.0
c
c---- check for triangel rules (can be skipped if know before) 
      if (j12.lt.abs(j11-j13).or.j12.gt.j11+j13) return
      if (j22.lt.abs(j21-j23).or.j22.gt.j21+j23) return
      if (j32.lt.abs(j31-j33).or.j32.gt.j31+j33) return
      if (j21.lt.abs(j11-j31).or.j21.gt.j11+j31) return
      if (j22.lt.abs(j12-j32).or.j22.gt.j12+j32) return
      if (j23.lt.abs(j13-j33).or.j23.gt.j13+j33) return
c--------------------------------------------------------------
c
      k1 = max0(iabs(j21-j32),iabs(j11-j33),iabs(j12-j23))
      k2 = min0(j21+j32,j12+j23,j11+j33)
      do k = k1,k2
         sninj0 = sninj0 + (k+k+1)*
     &            racah0(j32,j21,k,j11,j33,j31)*
     &            racah0(j12,j22,j32,j21,k,j23)*
     &            racah0(k,j23,j12,j13,j11,j33)
      enddo
c
      return
c-end-SNINJ0
      end

C=======================================================================

      double precision function racah0(j1,j2,j3,l1,l2,l3)

C=======================================================================
C
C     Calculates 6j-symbol (notation of Edmonds) for integer ang.momenta
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
      include 'dirhb.par'
c
      common /gfviv / iv(0:igfv)
      common /gfvfak/ fak(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
c
      dsq(i,k,l) = wf(i+k-l)*wfi(i+k+l+1)*wf(i-k+l)*wf(k+l-i)
c
      racah0 = 0.
      i1 = j1+j2+j3
      i2 = j1+l2+l3
      i3 = l1+j2+l3
      i4 = l1+l2+j3
      i5 = j1+j2+l1+l2
      i6 = j2+j3+l2+l3
      i7 = j3+j1+l3+l1
      n1 = max0(i1,i2,i3,i4)
      n2 = min0(i5,i6,i7)
      if (n1.gt.n2) return
      do 10 n = n1,n2
   10 racah0 = racah0+iv(n)*fak(n+1)*
     &         fi(n-i1)*fi(n-i2)*fi(n-i3)*fi(n-i4)*
     &         fi(i5-n)*fi(i6-n)*fi(i7-n)
   20 racah0 = dsq(j1,j2,j3)*dsq(j1,l2,l3)*
     &         dsq(l1,j2,l3)*dsq(l1,l2,j3)*racah0
c
      return
c-end-RACAH0
      end

c======================================================================c

      subroutine greecou(lpr)

c======================================================================c
C
C     calculation of the meson and Coulomb-propagator
c
c     imes = 1:   sigma
c            2:   omega
c            3:   delta
c            4:   rho
c            0:   coulomb
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tt*8                                            ! bloqua
c
      dimension gi(nox,nox),rr(nox,ngh)
c
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /procou/ ggc(0:ngh,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (icou.eq.0) return
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN GREECOU ******************************'
c
c
c---- Coulomb-propagator
      f0 = one/(6*alphi)
      do kh = 0,ngh
         f  = f0*wdcor(kh)
         do ih = 0,ngh
            r1 = xh(ih)
            r2 = xh(kh)
            rg = dmax1(r1,r2)
            rk = dmin1(r1,r2)
            ggc(ih,kh) =  f * ( 3*rg + rk**2/rg)
         enddo   ! ih
      enddo    ! ik
c      if (lpr) call aprint(1,1,6,ngh,ngh,ngh,ggc,' ',' ','VC')
c
      if (lpr)
     &write(l6,*) ' ****** END GREECOU ********************************'
c
      return
c-end-GREECOU
      end

c======================================================================c

      subroutine greemes(lpr)

c======================================================================c
C
C     calculation of the meson and Coulomb-propagator
c
c     imes = 1:   sigma
c            2:   omega
c            3:   delta
c            4:   rho
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
c
      dimension gi(nox,nox),rr(nox,ngh)
c
      common /baspar/ hom,hb0,b0
      common /bosqua/ no
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,ames(4)
      common /physco/ hbc,alphi,r0
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /propag/ gg(ngh,ngh,4)
      common /bospol/ rnb(1:nox,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN GREEMES ******************************'
c
      if (ipc.eq.1) return
c
c     meson-propagators
      do imes = 1,4
         f = (one/((ames(imes)+1.d-10)*b0))**2
         do i = 1,no
         do k = 1,no
            gi(i,k) = zero
         enddo   ! k
         enddo   ! i
         do n = 1,no
            gi(n,n) = one + f*(2*n-half) 
            if (n.lt.no) then
               gi(n,n+1) = f*sq(n)*sqh(n)
               gi(n+1,n) = gi(n,n+1)
            endif
            do ih = 1,ngh
               rr(n,ih) = rnb(n,ih)
            enddo   ! ih
         enddo   ! n
         call lingd(nox,nox,no,ngh,gi,rr,d,ifl)
         f0 = one/(4*pi*b0**3)
         do kh = 1,ngh
	        f = f0*wdcor(kh)
            do ih = 1,ngh
               s = zero
               do n = 1,no
                  s = s + rnb(n,ih)*rr(n,kh)
               enddo
               gg(ih,kh,imes) = f * s
            enddo
         enddo      
      enddo   ! imes
c
      if (lpr)
     &write(l6,*) ' ****** END GREEMES ********************************'
c
      return
c-end-GREEMES
      end

C=======================================================================
  
      subroutine lingd(ma,mx,n,m,a,x,d,ifl)

C=======================================================================
C
C     solves the system of linear equations A*X = B 
C     at the beginning the matrix B is stored in X
C     during the calculation it will be overwritten
C     D is the determinant of A
C
C-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
C
c     solves the system A*X=B, where B is at the beginning on X
c     it will be overwritten lateron, d is the determinant
C
      dimension  a(ma,m),x(mx,m)
C
      data tollim/1.d-10/,one/1.d0/,zero/0.d0/
C
      ifl=1 
      p=zero 
      do 10 i=1,n   
         q=zero         
         do 20 j=1,n
 20      q=q+ abs(a(i,j)) 
         if (q.gt.p)   p=q 
 10   continue         
      tol=tollim*p
      d=one
      do 30 k=1,n     
         p=zero           
         do 40 j=k,n   
            q = abs(a(j,k))
            if (q.lt.p) goto 40
            p=q 
            i=j 
 40      continue          
         if (p.gt.tol) goto 70
         write (6,200) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
  200    format (/1x,80a1/' *****  ERROR IN LINGD , TOLERANZ =',e10.4,
     1 ' VALUE OF A(',i3,',',i3,') IS ',e10.4/1x,80a1)
         ifl=-1                                         
         return
   70    cp=one/a(i,k)
         if (i.eq.k) goto 90
         d=-d
         do 81 l=1,m
            cq=x(i,l)
            x(i,l)=x(k,l) 
   81       x(k,l)=cq
         do 80 l=k,n
            cq=a(i,l)
            a(i,l)=a(k,l) 
   80       a(k,l)=cq
   90       d=d*a(k,k)
            if (k.eq.n) goto 1
            k1=k+1
            do 120 i=k1,n 
               cq=a(i,k)*cp
               do 106 l=1,m
  106             x(i,l)=x(i,l)-cq*x(k,l) 
               do 120 l=k1,n 
  120             a(i,l)=a(i,l)-cq*a(k,l) 
   30 continue
    1 do 126 l=1,m
  126    x(n,l)=x(n,l)*cp
         if (n.eq.1) return
         n1=n-1
         do 140 k=1,n1 
            cp=one/a(n-k,n-k)
            do 140 l=1,m
               cq=x(n-k,l)
               do 141 i=1,k
  141             cq=cq-a(n-k,n+1-i)*x(n+1-i,l)
  140          x(n-k,l)=cq*cp
c
c
      return
c-end-LINGD
      end 


