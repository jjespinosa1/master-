c     Esse programa faz :
c       - que o codigo delta2 imprima na tela as
c         normas de cada estado, e os coefs. da exp. do oscilador no 
c         arquivo dis.out. Conti 08/10/99
c       - inclui termos nao lineares para os mesons sigma e omega.
c         Conti 05/00
c     
c======================================================================c

c     PROGRAM SNU

c======================================================================c
c     Relativistic mean field theory in a spherical basis
CL--------------------------------------------------------------------CL
CL    This program introduces the mixing coupling for rho meson and
CL    extracts the contact term to the tensor coupling - 04/99
c----------------------------------------------------------------------c

c---- reads in data     
      call reader
c
c---- preparations
      call prep
c
c---- initialization of the potentials
C      call inout(1,.false .)
      call start(.true.)
c
c---- oscillator basis for single particle states
      call base(.true.)
c
c---- wavefunctions at Gauss-Meshpoints
      call radgh(.true.)

      !call gfv
      !a=4

      !write ( *, '(A, F10.2)' ) 'Oi', fak(8)
c
c---- single-particle matix elements
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!! Now this subroutine is iterative                             !!       
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CL      call singf(.false.)
c
c---- iteration
C      call iter(.true.)
c
c---- results
C      call resu
c
c---- punching of potentials
C      call inout(2,.true.)
c
c---- plotting of densities in coordinate space
c      call plotd(.true.)
c
c---- plotting of wavefunctions in coordinate space
C      call plw
      it = 1
      ib = 1
      k  = 1
c      call plotw(it,ib,k,.true.)
c
c---- equivalent Schroedinger equation
c     call schroed(.true.) 
c
c
c
      stop ' FINAL STOP'
c-end-DIS
      end


c======================================================================c

      blockdata default

c======================================================================c
c     
c     Default for Relativistic Mean Field spherical
c
c----------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
      include 'paramet'
c
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor
c
      common /baspar/ hom,hb0,b0f,b0b,bc,cf,cb
      common /dimens/ n0f,n0b,nrm,nlm
      common /fermi / ala(2),tz(2)
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /fixocc/ ioc(nbx,2)
      common /initia/ vin,rin,ain,inin,iplot
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /extrct/ gprime      
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,nneu,npro,jmax
      common /optopt/ icm,icou,it1,it2,ncut
      common /pair  / ga(2),gg(2),del(2),spk(2),dec(2),pwi
      common /physco/ amu,hqc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /temper/ temp
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
      common /ugugug/ itbl(2),jbl(2),ipbl(2),nbl(2),nrbl(2)
      common /woodsa/ v0,akv,vso(2),r0v(2),av(2),rso(2),aso(2)
c
c
c---- nucleus
      data nucnam/' O'/,nama/16/
      data jmax/8/
c
c---- blocking structure (for neutrons and protons)
c     itbl = 0: no blockung    = 1: one orbit is blocked
c     jbl, ipbl, nbl: quantum numbers j,ip,n  of the blocked level
      data itbl/0,0/,jbl/3,1/,ipbl/1,1/,nbl/1,1/,nrbl/0,0/
c
c---- options
c     center of mass:
      data icm/1/
c     protons or neutrons of both
      data it1/1/,it2/2/
c
c---- iteration
      data maxi/500/,si/1.0d0/,epsi/0.000001/
      data iaut/0/,inxt/3/
      data xmix0/0.1d0/,xmix/0.1d0/,xmax/0.8d0/
c
c
c---- force-parameters
c     Horotitz and Serot:
c     data txtfor/'HS'/
c     data amu/939.0/
c     data amsig/520.0/,amome/783.0/,amrho/770.0/
c     data gsig/10.47/,gome/13.80/,grho/8.07/,fr/0.0/
c     data g2/0.0/,g3/0.0/
c---------------------------------------------------------------
c     L1:
c     data txtfor/'L1'/
c     data amu/938.0/
c     data amsig/550.0/,amome/783.0/,amrho/763.0/
c     data gsig/10.30/,gome/12.60/,grho/0.0/,fr/0.0/
c     data g2/0.0/,g3/0.0/,w3/0.0/
c---------------------------------------------------------------
c     NL1:
c       data txtfor/'NL1'/
c       data amu/938.0/
c       data amsig/492.25/,amome/795.359/,amrho/763.0/
c       data gsig/10.138/,gome/13.285/,grho/4.9755/,fr/0.00/
c       data g2/-12.172/,g3/-36.265/,w3/0.0/
c---------------------------------------------------------------
c     NL1+frho:
c       data txtfor/'NL1+frho'/
c       data amu/938.0/
c       data amsig/492.25/,amome/795.359/,amrho/763.0/
c       data gsig/10.138/,gome/13.285/,grho/4.9755/,fr/9.95/
c       data g2/-12.172/,g3/-36.265/,w3/0.0/
c---------------------------------------------------------------
c     NL2:
c      data txtfor/'NL2'/
c      data amu/938.0/
c      data amsig/504.89/,amome/780.0/,amrho/763.0/
c      data gsig/9.111/,gome/11.493/,grho/5.507/,fr/18.41/
c      data g2/-2.304/,g3/13.783/,w3/0.0/
c---------------------------------------------------------------
c     NL3: 
c      data txtfor/'NL3'/
c      data amu/939.0/
c      data amsig/508.194/,amome/782.501/,amrho/763.0/
c      data gsig/10.217/,gome/12.868/,grho/4.474/,fr/32.84/
c      data g2/-10.431/,g3/-28.885/,w3/0.0/
c---------------------------------------------------------------
c     NLSH: (good version)
c     data txtfor/'NLSH'/
c     data amu/939.0/
c     data amsig/526.059/,amome/783.0/,amrho/763.0/
c     data gsig/10.444/,gome/12.945/,grho/4.383/,fr/0.0/
c     data g2/-6.9099/,g3/-15.8337/,w3/0.0/
c---------------------------------------------------------------
c     Bouyssy1:
c      data txtfor/'Bouyssy1'/
c      data amu/938.9/
c      data amsig/440.0/,amome/783.0/,amrho/770.0/
c      data gsig/5.341/,gome/11.21/,grho/2.63/,fr/17.358/
c      data g2/-12.172/,g3/-36.265/,w3/0.0/
c---------------------------------------------------------------
c     Bouyssy2:
c      data txtfor/'Bouyssy2'/
c      data amu/938.9/
c      data amsig/440.0/,amome/783.0/,amrho/770.0/
c      data gsig/7.230/,gome/11.853/,grho/2.629/,fr/9.727/
c      data g2/-12.172/,g3/-36.265/,w3/0.0/
c---------------------------------------------------------------
c    NL1+Bouyssy2:
c      data txtfor/'NL1+Bouyssy2'/
c      data amu/938.0/
c      data amsig/492.25/,amome/795.359/,amrho/770.0/
c      data gsig/10.138/,gome/13.285/,grho/2.629/,fr/9.727/
c      data g2/-12.172/,g3/-36.265/,w3/0.0/
c---------------------------------------------------------------
c     NL-VT1:
c       data txtfor/'NL-VT1'/
c       data amu/938.9/
c       data amsig/484.307/,amome/780.000/,amrho/763.000/
c       data gsig/9.81307/,gome/12.6504/,grho/4.63432/,fr/21.8343/
c       data g2/-13.2808/,g3/-38.0773/,w3/0.0/
c---------------------------------------------------------------
c     TM1:
c      data txtfor/'TM1'/
c      data amu/938./
c      data amsig/511.198/,amome/783.000/,amrho/770.000/
c      data gsig/10.0289/,gome/12.6139/,grho/4.6322/,fr/0.00000/
c      data g2/-7.2325/,g3/0.6183/,w3/71.3075/
c---------------------------------------------------------------
c     TM1 + frho:
      data txtfor/'TM1+frho'/
      data amu/938./
      data amsig/511.198/,amome/783.000/,amrho/770.000/
      data gsig/10.0289/,gome/12.6139/,grho/4.6322/,fr/17.14/
      data g2/-7.2325/,g3/0.6183/,w3/71.3075/
c---------------------------------------------------------------
c
c---- gprime = 1/3 extracts the contact part from derivative coupling  
c                  to the rho meson 
      data gprime/0.d0/
c      data gprime/0.333333333333333333d0/  
c---------------------------------------------------------------
c      
c---- Coulomb-field: not at all (0), direct term (1), plus exchange (2)
      data icou/1/
c---------------------------------------------------------------
c
c
c---- pairing
      data dec/0.d0,0.d0/
      data del/6.d0,6.d0/
      data ala/-7.d0,-7.d0/
c     data ga/22.,27./
      data ga/0.0,0.0/
c
c---- temperture
      data temp/0.0/
c
c---- parameters of the initial potentials
c     inin = 0: fields read, 1: default, 2: saxon-wood, 3:oscillator
      data inin/2/
c
c     Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
c      data v0/-71.28/,akv/0.4616/
c      data r0v/1.2334,1.2496/,av/0.615,0.6124/
c      data vso/11.1175,8.9698/
c      data rso/1.1443,1.1401/,aso/0.6476,0.6469/
c---- Saxon-Woods parameter for test without Coulomb
      data v0/-71.28/,akv/0.4616/
      data r0v/1.2334,1.2334/,av/0.6150,0.6150/
      data vso/11.1175,11.1175/
      data rso/1.1443,1.1443/,aso/0.6476,0.6476/
c
c     data vin/-55.0/,rin/1.0/,qin/1.3/,ain/0.6/
c     Woods-Saxon Potential Dudek
c     data wsv0/-49.6/,wsv1/0.86/,wsr/1.347,1.275/,wsa/0.7/
c
c---- basis parameters:
c     number of major oscillator shells 
      data n0f/2/,n0b/20/
c     oscillator length b0f (is calcuated for b0f <= 0)
      data b0f/-2.320/
c
c
c---- tapes
      data l6/10/,lin/3/,lou/6/,lwin/1/,lwou/2/,lplo/11/
c
c---- fixed texts
      data tp/'+','-'/,tis/'n','p'/
      data tl/'s','p','d','f','g','h','i','j','k','l','m',
     &            'n','o','P','q','r','S','t','u','v','w'/
c
c
c---- physical constants
      data hqc/197.328284d0/,r0/1.2/,alphi/137.03602/
c
c---- mathemathical constants
c     are determined in GFV
c      data zero/0.0d0/,one/1.d0/,two/2.d0/
c      data half/0.5d0/,third/0.333333333333333333d0/
c      data pi/3.141592653589793d0/
c
c
c---- fixed occupation patterns
c     O-16
c     data ioc/1,1,0,1,0,0,0,0,0,0, 32*0,
c    &         1,1,0,1,0,0,0,0,0,0, 32*0/
c
c     Pb-208
c     data ioc/3,3,2,3,2,2,1,2,1,1, 0,1,1,29*0,
c    &         3,2,2,2,2,1,1,1,1,0, 0,1,0,29*0/
c
c     GG-298
c     data ioc/4,3,3,3,3,2,2,2,2,1, 1,1,1,0,0,1,26*0,
c    &         3,2,2,3,2,2,1,2,1,0, 0,1,1,0,0,0,26*0/
c
c---- fields O-16 for NL1, ngh = 12
c     data sig/0.d0,
c    1 -.197896732795d+00, -.203336593050d+00, -.175178199706d+00,
c    2 -.110773903644d+00, -.492821497342d-01, -.155056088258d-01,
c    3 -.370467731601d-02, -.752571143561d-03, -.128392213020d-03,
c    4 -.178840524301d-04, -.123043868927d-06, -.138401023191d-05/
c     data ome/0.d0,
c    1  .122884775525d+00,  .128313360162d+00,  .111273919274d+00,
c    2  .678286964552d-01,  .271388311726d-01,  .717766421530d-02,
c    3  .140426026681d-02,  .243999373057d-03,  .378195438601d-04,
c    4  .312885710667d-05,  .507556868833d-06, -.331104601547d-06/
c     data rho/0.d0,
c    1 -.402717077863d-03, -.401598408644d-03, -.327222122144d-03,
c    2 -.134890782160d-03,  .226171401831d-04,  .537452482103d-04,
c    3  .302903348538d-04,  .103608532225d-04,  .228698677638d-05,
c    4  .197570927351d-06,  .514723017016d-07, -.380613801692d-07/
c     data cou/0.d0,
c    1  .280242732120d-01,  .268490299578d-01,  .243582917371d-01,
c    2  .209044213183d-01,  .173601945696d-01,  .143915626933d-01,
c    3  .121166187436d-01,  .103808124163d-01,  .901598395252d-02,
c    4  .790128300969d-02,  .694924512594d-02,  .607533630278d-02/
c
c-end-DEFAULT
      end
c======================================================================c

      blockdata gauss

c======================================================================c
c
c     12 Meshpoints 
c     ph  =  wh * exp(-xh)
c
c     \int_0^\infty  f(z) exp(-z^2) dz   =   \sum_i f(xh(i)) ph(i)
c     \int_0^\infty  f(z) dz             =   \sum_i f(xh(i)) wh(i) 
c
c----------------------------------------------------------------------c
      include 'paramet'
      implicit real*8 (a-h,o-z)
c
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
c
      DATA XH  /0.000001E+00,
     1   0.22441453833E+00,   0.67417108003E+00,   0.11267607739E+01,
     2   0.15842499524E+01,   0.20490035025E+01,   0.25238809358E+01,
     3   0.30125460488E+01,   0.35200067194E+01,   0.40536643061E+01,
     4   0.46256626595E+01,   0.52593828320E+01,   0.60159254688E+01/
      DATA WH  /0.00000000001E+00,
     1   0.44898282712E+00,   0.45084608343E+00,   0.45467546127E+00,
     2   0.46069389583E+00,   0.46928487425E+00,   0.48107219405E+00,
     3   0.49707668668E+00,   0.51904386381E+00,   0.55020754585E+00,
     4   0.59739170078E+00,   0.67857466439E+00,   0.86887556926E+00/
      DATA PH  /0.00000000001E+00,
     1   0.42693114831E+00,   0.28617953473E+00,   0.12773962990E+00,
     2   0.37445476328E-01,   0.70483576945E-02,   0.82369280545E-03,
     3   0.56886946065E-04,   0.21582471124E-05,   0.40189743017E-07,
     4   0.30462570024E-09,   0.65846268900E-12,   0.16643703593E-15/
c
c-end BLOCKDATA GAUSS
      end

c======================================================================c

      subroutine reader

c======================================================================c
c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor
c
      common /baspar/ hom,hb0,b0f,b0b,bc,cf,cb
      common /dimens/ n0f,n0b,nrm,nlm
      common /initia/ vin,rin,ain,inin,iplot
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,nneu,npro,jmax
      common /pair  / ga(2),gg(2),del(2),spk(2),dec(2),pwi
      common /physco/ amu,hqc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
      common /temper/ temp
c
c
c
      if (l6.ne.6) open(l6,file='dis.out',status='unknown')
      if (lin.eq.0) return
      open(lin,file='snu.dat',status='old')
c
c
c---- Output-File:            
      read(lin,'(10x,i5)') l6
      write(l6,*) ' ****** BEGIN READER *******************************'
      write(l6,'(a,i5)') ' Output file                 : ',l6 
c
c---- Basisparameters:            
      read(lin,'(10x,2i5)') n0f,n0b
      write(l6,'(a,2i5)') ' Number of oscillator shells : ',n0f,n0b
      read(lin,'(10x,f10.3)') b0f
      write(l6,'(a,f9.3)') ' Oscillator length b0f (fm)  : ',b0f
c
c---- Parameter for the iteration:
      read(lin,'(10x,i5)') maxi
      write(l6,'(a,i5)') ' Maximal number of iterations: ',maxi
      read(lin,'(10x,f10.3)') xmix
      write(l6,'(a,f9.3)') ' Mixing parameter            : ',xmix
      xmix0 = xmix
c
c---- Initialization of wavefunctions:
      read(lin,'(10x,i5)') inin
      write(l6,'(a,i5)') ' Initial wavefunctions       :  ',inin
c
c---- Nucleus under consideration
      read(lin,'(a2,i4)') nucnam,nama
      write(l6,'(a,20x,a2,i4)') ' Nucleus:      ',nucnam,nama
c
c---- Gap-Parameters
      read(lin,'(10x,2f10.3)') dec
      write(l6,'(a,2f10.3)') ' Gap Parameters             : ',dec
c
c---- Pairing-Constants
      read(lin,'(10x,2f10.3)') ga
      write(l6,'(a,2f10.3)') ' Pairing-Constants          : ',ga
c
c---- Temperature    
      read(lin,'(10x,f10.3)') temp
      write(l6,'(a,2f10.3)') ' Temperature                : ',temp
c---- Plot of Wavefunctions
      read(lin,'(10x,i4)') iplot
      write(l6,'(a,i5)') ' Plot of Wavefunctions      : ',iplot
c
c      close(lin)
c
      write(l6,*) ' ****** END READER *********************************'
c-end-reader 
      end

c======================================================================c

      subroutine start(lpr)

c======================================================================c
c
c     initializes potentials
c     inin = 0:   reads fields from tape lwin
c            1:   saxon-woods
c            2:   uses default fields
c            3:   oscillator
c
c----------------------------------------------------------------------c
c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /initia/ vin,rin,ain,inin,iplot
      common /mathco/ zero,one,two,half,third,pi
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /nucnuc/ amas,nama,npr(2),jmax
      common /optopt/ icm,icou,it1,it2,ncut
      common /physco/ amu,hqc,alphi,r0
      common /potpot/ ss(0:ngh,1:2),vv(0:ngh,1:2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /woodsa/ v0,akv,vso(2),r0v(2),av(2),rso(2),aso(2)
c
c---- potentials are read in INOUT
      if (inin.eq.0) return
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN START ********************************'
c
c---- default values
      if (inin.eq.2) then
         write(l6,'(/,a)') 
     &      ' uses default values as initial fields '
      endif

c---- saxon-woods potential
      if (inin.eq.1) then
c
         write(l6,'(a,f10.4)') ' v0     = ',v0
         write(l6,'(a,f10.4)') ' kappa  = ',akv
         write(l6,'(a,2f10.4)') ' lambda = ',vso
         write(l6,'(a,2f10.4)') ' r0     = ',r0v
         write(l6,'(a,2f10.4)') ' a      = ',av
         write(l6,'(a,2f10.4)') ' r0-so  = ',rso
         write(l6,'(a,2f10.4)') ' a-so   = ',aso
         do 10 it = 1,2
            ita = 3-it
            rav = r0v(it)*amas**third
            rao = rso(it)*amas**third
            vp  = v0*(1 - akv*(npr(it)-npr(ita))/amas)
c           vls = half*(hqc/amu)**2 * vp * vso(it)
            vls = vp * vso(it)
            do ih = 0,ngh
               r = xh(ih)*b0c
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
c              hb(ih,it) = hb0
c              v(ih,it)  = u
c              vs(ih,it) = w
               ss(ih,it) = half * (u - w)
               vv(ih,it) = half * (u + w)
c              ss(ih,it) = half*(u-w/(1-w/(2*amu)))
c              vv(ih,it) = half*(u+w/(1-w/(2*amu)))
            enddo
   10    continue
c
c------- Coulomb potenital
         if (icou.eq.0) then
            do ih =0,ngh
               cou(ih) = zero
            enddo
         else
            rc = r0v(2)*amas**third
            do ih = 0,ngh
               r = xh(ih)*b0c
               if (r.lt.rc) then
                  c = half*(3/rc - r*r/rc**3)
               else
                  c = one/r
               endif
               cou(ih) = c*npr(2)/alphi
            enddo
         endif
         write(l6,'(/,a)') 
     &      ' Initial potentials of Saxon-Woods shape '
      endif
c
c
c---- calculation of the meson-fields
      if (inin.ge.2) then
         fsig = half/(hqc*gsig)
         fome = half/(hqc*gome)
         frho = half/(hqc*grho)
         do ih = 0,ngh
            sig(ih) = fsig * ( ss(ih,1) + ss(ih,2))
            ome(ih) = fome * ( vv(ih,1) + vv(ih,2))
            rho(ih) = frho * (-vv(ih,1) + vv(ih,2))
         enddo
      endif
c
      if (lpr) then
            call prigh(0,sig,b0c,'Mesh ')
            call prigh(1,sig,one,' SIG')
            call prigh(1,ome,one,' OME')
            call prigh(1,rho,one,' RHO')
            call prigh(1,cou,one,' COU')
      endif    
c
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
c
c----------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor
      character*8 tb
      character*25 txb
c
      common /baspar/ hom,hb0,b0f,b0b,bc,cf,cb
      common /dimens/ n0f,n0b,nrm,nlm
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)
      common /nucnuc/ amas,nama,nneu,npro,jmax
      common /quaqua/ nt,nr(ntx),nl(ntx),nj(ntx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /texblo/ tb(ntx),txb(nbx)
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
      common /ugugug/ itbl(2),jbl(2),ipbl(2),nbl(2),nrbl(2)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN BASE *********************************'
c
      write(6,*) 'n0f,n0fx',n0f,n0fx
      if (n0f.gt.n0fx) stop ' in BASE: n0f too large'
      je = min0(n0f+1,jmax)
c
      ib     = 0
      ilauf  = 0
      nrm    = 0
      nlm    = 0
c
c---- loop over j-quantum number
      do 10 j = 1,je
c
c---- loop over parity
      do 10 ip = 1,2
         l = j - mod(j+ip+1,2)   
         nlm = max0(nlm,l)
c
c        write(l6,'(/,a,i3,a,a1)') ' j = ',j+j-1,'/2',tp(ip)
         ne = (n0f+1-l)/2 + 1
         jlauf = ilauf
         do 20 ir = 1,ne
            ilauf  = ilauf + 1
            if (ilauf.gt.ntx) stop ' in BASE: ntx too small'
            nr(ilauf) = ir
            nl(ilauf) = l
            nj(ilauf) = j
            mu(ilauf) = j
            write(tb(ilauf),'(i2,a1,i3,2h/2)') ir,tl(l),j+j-1
            nn        = 2*(ir-1)+l
c           write(l6,'(i4,a,i2,a,i2,a,i2)') 
c    &           ilauf,'  N = ',nn,'   n = ',ir,'   l = ',l
            nrm = max0(nrm,ir)
c
c           determination of blocked level
            do it = 1,2
               if ( itbl(it).eq.1  .and. jbl(it).eq.j .and. 
     &              ipbl(it).eq.ip .and. nbl(it).eq.ir) nrbl(it)=ilauf
            enddo
   20    continue
c          
         if (ilauf.gt.jlauf) then
            ib = ib + 1
            if (ib.gt.nbx)  stop ' in BASE: nbx too small'
            ia(ib)  = jlauf 
            id(ib)  = ilauf-jlauf
            if (id(ib).gt.ndx) stop ' in BASE: ndx too small'
            ijb(ib) = j
            ilb(ib) = l
            write(txb(ib),'(i3,a,i2,a,a1)') 
     &            ib,'. block:  j = ',j+j-1,'/2',tp(ip)
         endif
c
   10 continue
      nb  = ib 
      nt  = ilauf
      if (nrm.gt.nrx) stop 'in BASE: nrx too small '
      if (nlm.gt.nlx) stop 'in BASE: nlx too small '
c
c---- determination of the corresponding small components quantum numbers
      do ib = 1,nb
         ibq = ib - 1 + 2*mod(ib,2)
         idq(ib) = id(ibq)
         iaq(ib) = ia(ibq)
      enddo
c
c---- printout
      if (lpr) then
         do 40 ib = 1,nb
            j  = ijb(ib)
            ip = mod(ilb(ib),2)+1
            i1 = ia(ib)+1
            i2 = ia(ib)+id(ib)
            write(l6,'(/,i3,a,i3,a,a1)') 
     &                ib,'. block:  j = ',j+j-1,'/2',tp(ip)
            do i = i1,i2
               n = nr(i)
               l = nl(i)
               nn = 2*(n-1)+l
               write(l6,'(i4,a,i2,a,i2,a,i2)') 
     &               i,'  N = ',nn,'   n = ',n,'   l = ',l
            enddo
            i1 = iaq(ib)+1
            i2 = iaq(ib)+idq(ib)
            do i = i1,i2
               n = nr(i)
               l = nl(i)
               nn = 2*(n-1)+l
               write(l6,'(i4,a,i2,a,i2,a,i2)') 
     &               i,'  N = ',nn,'   n = ',ir,'   l = ',l
            enddo
   40    continue
         write(l6,'(/,a,2i4)') ' Number of blocks: nb  = ',nb,nbx
         write(l6,100) ' Number of levels  nt  = ',nt,ntx 
         write(l6,100) ' Maximal n:        nrm = ',nrm,nrx
         write(l6,100) ' Maximal l:        nlm = ',nlm,nlx
         write(l6,100) ' Blocked levels:   nrbl= ',nrbl
  100    format(a,2i4)
      endif
c    
      if (lpr)
     &write(l6,*) ' ****** END BASE ***********************************'
      return
c-end-BASE
      end

c=====================================================================c

      subroutine prigh(is,ff,f,text)

c=====================================================================c
C
c     is=0  prints gauss-meshpoints * b0  (b0 = f)
c        1  prints f*ff(x) at gauss-meshpoints 
c        2  prints f*ff(x)/(wh(ih)*xh(i)**2)
c
c---------------------------------------------------------------------c
      include 'paramet'
c     
      implicit real*8 (a-h,o-z)
      character text*(*)
c
      dimension ff(0:ngh)
c
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      data ix/6/
c
      if (is.eq.0) then
         write(l6,100) text,(f*xh(ih),ih=0,ix)
  100    format(/,1x,a6,12f10.4) 
      elseif (is.eq.1) then  
         write(l6,101) text,(f*ff(ih),ih=0,ix)
  101    format(1x,a6,12f10.4)
      elseif (is.eq.2) then  
         write(l6,101) text,(f*ff(ih)/(wh(ih)*xh(ih)**2),ih=0,ix)
      endif
c
      return
c-end-PRIGH
      end

c======================================================================c

      subroutine prep

c======================================================================c
c
c     preparations
c
c----------------------------------------------------------------------c
c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /gaucor/ rb(0:ngh),wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /nucnuc/ amas,nmas,nneu,npro,jmax
      common /optopt/ icm,icou,it1,it2,ncut
      common /pair  / ga(2),gg(2),del(2),spk(2),dec(2),pwi
      common /physco/ amu,hqc,alphi,r0
      common /potpot/ ss(0:ngh,1:2),vv(0:ngh,1:2)
      common /rhorho/ ros(0:ngh,2),rov(0:ngh,2),dro(0:ngh)     
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
c
c
cdiz  rtest = sqrt(0.22441453833d+00**2 + 0.11572211736d+00)
      rtest = sqrt(3*xh(1)**2)
c
      write(l6,*) ' ****** BEGIN PREP *********************************'
c
c---- signs and factorials
      call gfv
c
c---- nuclear parameters
      call nucleus(2,npro,nucnam)
      nneu = nmas - npro
      amas = nmas
      if (nneu.eq.0) it1 = 2
      if (npro.eq.0) it2 = 1
      write(l6,'(a,a,i4,i6,i4)') ' Nucleus: ',nucnam,nmas,nneu,npro
c
c---- basis parameters
      hb0 = hqc**2/(two*amu)
      hom = 41.0*amas**(-third)    
      if (icm.gt.0) hb0 = hb0*(one - one/amas)
      if (b0f.le.0.0) then
          b0f = sqrt(two*hb0/hom)
          write(l6,*) ' b0f is calculated: ',b0f
      endif
c      hom = one
c      hb0 = half
c      b0  = one
c      xh(0) = rtest
      wh(0) = one
      b0b   = b0f/sqrt(two)
      b0c   = one/sqrt(one/b0f**2+half/b0b**2)
      cf    = b0c/b0f
      cb    = b0c/b0b
      do ih = 0,ngh
         rb(ih) = xh(ih)*b0f
c
c        metric element for three-dimensional integration 
         wdcor(ih) = b0f**3 * 4*pi * xh(ih)**2 * wh(ih)
      enddo
c
c---- pairing force
c     dec(1) = 12.d0/sqrt(nneu*one)
c     dec(2) = 12.d0/sqrt(npro*one)
c      dec(1) = 12.d0/sqrt(amas)
c      dec(2) = 12.d0/sqrt(amas)
      do it = 1,2
         gg(it) = ga(it)/amas
         if (gg(it).eq.zero)  del(it) = dec(it)
      enddo
      pwi = 2*hom
c     pwi = zero
c
      write(l6,*) ' hom = ',hom
      write(l6,*) ' hb0 = ',hb0
      write(l6,'(a,3f15.8)') ' b0  = ',b0f,b0b,b0c 
c
c
c---- printout of force:
      csig = gsig*amu/amsig
      come = gome*amu/amome
      crho = grho*amu/amrho
      write(l6,'(/,a,a)') ' Meson-Parameters: ',txtfor
      write(l6,'(a,f10.4,a,f10.4,a,f10.4)') 
     &          ' msig  = ',amsig,'  gsig = ',gsig,'  Csig',csig
      write(l6,'(a,f10.4,a,f10.4,a,f10.4)') 
     &          ' mome  = ',amome,'  gome = ',gome,'  Come',come
      write(l6,'(a,f10.4,a,f10.4,a,f10.4,a,f8.4)') 
     &          ' mrho  = ',amrho,'  grho = ',grho,'  Crho',crho,
     &          '  frho =',fr
      write(l6,'(a,f10.4,a,f10.4)') ' g2    = ',g2  
      write(l6,'(a,f10.4,a,f10.4)') ' g3    = ',g3
      write(l6,'(a,f10.4,a,f10.4)') ' c3    = ',w3  
c
c---- printout pairing:
      write(l6,'(a,4f10.6)') '  Gap parameter = ',dec,del
      write(l6,'(a,2(f8.4,2h/A),2f10.6)') '  Pairing const.= ',ga,gg
      write(l6,'(a,4f10.6)') '  Pairing window= ',pwi/hom

c
c
      if (icou.eq.0) write(l6,*) ' without Coulomb force'
      if (icou.eq.1) write(l6,*) ' with Coulomb force'
      if (icou.eq.2) write(l6,*) ' with Coulomb force with exchange'
      if (icm.eq.1)  write(l6,*) ' with center of mass correction '
      write(l6,*) ' Mixing-Parameter xmix: ',xmix
c
c
      do it = 1,2
      do ih = 0,ngh
            ss(ih,it)  = zero
            vv(ih,it)  = zero
            ros(ih,it) = zero
            rov(ih,it) = zero
      enddo
      enddo
c
c
c
      write(l6,*) ' ****** END PREP ***********************************'
      return
c-end PREP
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
c     fi(n)  =  1/n!
c     wf(n)  =  sqrt(n!)
c     wfi(n) =  1/sqrt(n!)
c     gm2(n) =  gamma(n+1/2)
c     gmi(n) =  1/gamma(n+1/2)
c     wg(n)  =  sqrt(gamma(n+1/2))
c     wgi(n) =  1/sqrt(gamma(n+1/2))
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      include 'paramet'
c
      common /gfviv / iv(0:igfv)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvfak/ fak(0:igfv)
      common /gfvfad/ fad(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfvgm2/ gm2(0:igfv)
      common /gfvgmi/ gmi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
c
c---- mathemathical constants
      data zero/0.0d0/,one/1.d0/,two/2.d0/
      data half/0.5d0/,third/0.333333333333333333d0/
      data pi/3.141592653589793d0/
C
      third = one/3.d0
      pi = 4*atan(one)
c
      iv(0)  = +1
      sq(0)  =  zero
      sqi(0) =  1.d30
      sqh(0) =  sqrt(half)
      shi(0) =  1/sqh(0)
      fak(0) =  1
      fad(0) =  1
      fi(0)  =  1
      wf(0)  =  1
      wfi(0) =  1
c     gm2(0) = Gamma(1/2) = sqrt(pi)
      gm2(0) =  sqrt(pi)
      gmi(0) =  1/gm2(0)
      wg(0)  =  sqrt(gm2(0))
      wgi(0) =  1/wg(0)
      do i = 1,igfv
         iv(i)  = -iv(i-1)
         sq(i)  = dsqrt(dfloat(i))
         sqi(i) = 1/sq(i)
         sqh(i) = sqrt(i+half)
         shi(i) = 1/sqh(i)
         fak(i) = i*fak(i-1)
	 fad(i) = (2*i+1)*fad(i-1)
         fi(i)  = 1/fak(i)
         wf(i)  = sq(i)*wf(i-1)
         wfi(i) = 1/wf(i)
         gm2(i) = (i-half)*gm2(i-1)
         gmi(i) = 1/gm2(i)
         wg(i)  = sqh(i-1)*wg(i-1)
         wgi(i) = 1/wg(i)
      enddo
c
c     write(6,*) ' ****** END GFV *************************************' 
      return
c-end-GFV
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
      PARAMETER (MAXZ=110)
C
      CHARACTER TE*2,T*(2*MAXZ+2)
C
      T(  1: 40) = '   HHELIBE B C N O FNENAMGALSI P SCLAR K'
      T( 41: 80) = 'CASCTI VCRMNFECONICUZNGAGEASSEBRCRRBSR Y'
      T( 81:120) = 'ZRNBMOTCRORHPDAGCDINSNSBTE OXECSBA;ACEPR'
      T(121:160) = 'NDPMSMEUGDTBDYHOERTMYBLUHFTA WREODIRPTAU'
      T(161:200) = 'HGTLPBBIPOATRNFRRAACTHPA UNPPUAMCMBKCFES'
      T(201:222) = 'FMMDNOLR040506070809GG'
C
      if (is.eq.1) then
         if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong' 
         te = t(2*npro+1:2*npro+2)
         return
      else
         do np = 0,maxz
            if (te.eq.t(2*np+1:2*np+2)) then
               npro = np
	       if (npro.eq.110) npro = 114
               return
            endif
         enddo
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
      endif
c
      stop
C-END-NUCLEUS
      END

c======================================================================c

      subroutine radgh(lpr)

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
c
c     RNB(n,i) is the radial function for the expansion of the mesonfields
c     differnently defined form RNL: 
c
c     RNB(n,i) = R_n0(r) 
c
c----------------------------------------------------------------------c
c   
c
c
c----------------------------------------------------------------------c
c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /dimens/ n0f,n0b,nrm,nlm
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /rados1/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /rados2/ rnl2(1:nrx,0:nlx,0:ngh),rnl12(1:nrx,0:nlx,0:ngh)
      common /radbos/ rnb(1:nox,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN RADGH ********************************'
c
      f = 2*wgi(0)
      nbo = n0b/2+1
c
      do 10 ih = 0,ngh
         r  = xh(ih)
         rr = r*r
         ri = one/r 
         fe = f*exp(-half*rr)
c
         r2  = xh(ih)*cf
         rr2 = r2*r2
         ri2 = one/r2
         fe2 = f*exp(-half*rr2)        
c
c------- basis for fermions
c------------------------------------
c        renormalization for fermions
         u1 = fe*sqrt(wh(ih)*rr)
         u2 = fe2*sqrt(wh(ih)*rr2*cf)
c------------------------------------
         do l = 0,nlm
            rnl(1,l,ih)  = u1
            rnl(2,l,ih)  = u1*(l+1.5d0-rr)*shi(l+1)
            u1           = u1*r*shi(l+1)
            rnl1(1,l,ih) =    (l-rr)*rnl(1,l,ih)*ri
            rnl1(2,l,ih) = ((2+l-rr)*rnl(2,l,ih) - 
     &                       2*sqh(l+1)*rnl(1,l,ih))*ri
c
            rnl2(1,l,ih) = u2
            rnl2(2,l,ih) = u2*(l+1.5d0-rr2)*shi(l+1)
            u2           = u2*r2*shi(l+1)            
            rnl12(1,l,ih)=    (l-rr2)*rnl2(1,l,ih)*ri2
            rnl12(2,l,ih)= ((2+l-rr2)*rnl2(2,l,ih) - 
     &                       2*sqh(l+1)*rnl(1,l,ih))*ri2
c
            do n = 3,nrm
               rnl(n,l,ih)  = ((2*n+l-2.5d0-rr)*rnl(n-1,l,ih) -
     &           sq(n-2)*sqh(n-2+l)*rnl(n-2,l,ih))*sqi(n-1)*shi(n-1+l)
               rnl1(n,l,ih) = ((2*n+l-2-rr)*rnl(n,l,ih) -
     &           2*sq(n-1)*sqh(n-1+l)*rnl(n-1,l,ih))*ri
               rnl2(n,l,ih) = ((2*n+l-2.5d0-rr2)*rnl2(n-1,l,ih) -
     &           sq(n-2)*sqh(n-2+l)*rnl2(n-2,l,ih))*sqi(n-1)*shi(n-1+l)               
               rnl12(n,l,ih)= ((2*n+l-2-rr2)*rnl2(n,l,ih) -
     &           2*sq(n-1)*sqh(n-1+l)*rnl2(n-1,l,ih))*ri2
            enddo
         enddo
c
c
c---- basis for bosons
         r  = xh(ih)*cb
         rr = r*r
         fe = f*exp(-half*rr)
c         
         rnb(1,ih) = fe
         rnb(2,ih) = rnb(1,ih)*(1.5d0-rr)*shi(1)
         do n = 3,nbo
            rnb(n,ih)  = ((2*n-2.5d0-rr)*rnb(n-1,ih) -
     &          sq(n-2)*sqh(n-2)*rnb(n-2,ih))*sqi(n-1)*shi(n-1)
         enddo
c
   10 continue
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
               do ih = 1,ngh
                  rr = xh(ih)**2
                  s0 = rnl(n1,l,ih)*rnl(n2,l,ih)
                  s1 = s1 + s0
                  s2 = s2 + rr*s0
                  s3 = s3 + (rnl1(n1,l,ih)*rnl1(n2,l,ih)
     &                       + rnl1(n1,l,ih)*rnl(n2,l,ih)/xh(ih)
     &                       + rnl(n1,l,ih)*rnl1(n2,l,ih)/xh(ih)
     &                       + s0*(1+l*(l+1))/rr)
 
                  if (l.eq.0) then
                     sb = sb + rnb(n1,ih)*rnb(n2,ih)*rr*wh(ih)*cb**3
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
     &write(l6,*) ' ****** END RADGH **********************************'
      return
c-end-RADGH
      end
