c     Esse programa faz :
c       - que o codigo delta2 imprima na tela as
c         normas de cada estado, e os coefs. da exp. do oscilador no 
c         arquivo dis.out. Conti 08/10/99
c       - inclui termos nao lineares para os mesons sigma e omega.
c         Conti 05/00
c     
c======================================================================c

      PROGRAM SNU

c======================================================================c
c     Relativistic mean field theory in a spherical basis
CL--------------------------------------------------------------------CL
CL    This program introduces the mixing coupling for rho meson and
CL    extracts the contact term to the tensor coupling - 04/99
c----------------------------------------------------------------------c
      
      

c
      !call gfv
      

c---- reads in data     
      call reader
c
c---- preparations
      call prep
c
c---- initialization of the potentials
      call inout(1,.false .)
      call start(.false.)
c
c---- oscillator basis for single particle states
      call base(.false.)
c
c---- wavefunctions at Gauss-Meshpoints
      call radgh(.false.)
c
c---- single-particle matix elements
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!! Now this subroutine is iterative                             !!       
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CL      call singf(.false.)
c
c---- iteration
      call iter(.true.)
c
c---- results
      call resu
c
c---- punching of potentials
      call inout(2,.true.)
c
c---- plotting of densities in coordinate space
      call plotd(.false.)
c
c---- plotting of wavefunctions in coordinate space
      call plw
      it = 1
      ib = 1
      k  = 1
      call plotw(it,ib,k,.true.)
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

      subroutine aprint(is,it,ns,ma,n1,n2,a,t1,t2,text)

c======================================================================c
C
C     IS = 1    Full matrix  
C          2    Lower diagonal matrix    
c          3    specially stored symmetric matrix
C 
C     IT = 1    numbers for rows and columns
C          2    text for rows and numbers for columns
C          3    text for rows and columns
C
C     NS = 1     FORMAT  8F8.4      80 Coulums
C     NS = 2     FORMAT  8f8.2      80 Coulums
C     NS = 3     FORMAT 17F4.1      80 Coulums
C     NS = 4     FORMAT 30F4.1     120 Coulums
C     NS = 5     FORMAT  5F12.8     80 Coulums
C     NS = 6     FORMAT 10F12.8    130 Coulums
C     NS = 7     FORMAT  4E13.6     80 Coulums
C     NS = 8     FORMAT  8E15.8    130 Coulums
c
c----------------------------------------------------------------------c
      implicit double precision (a-h,o-z)
C
      character*8 t1(1),t2(1)
      character text*(*)
C
      dimension a(1)
C
      character*30 fmt1,fmt2
      character*20 fti,ftt,fmt(8),fmti(8),fmtt(8)
      dimension nsp(8)
c
      common /tapes / l6,lin,lou,lwin,lwou,lplo
C
      data nsp/8,8,17,30,5,10,4,8/
      data fmt /'8f8.4)',            '8F8.2)',
     &          '17f4.1)',           '30f4.1)',
     &          '5f12.8)',           '10f12.8)',
     &          '4e13.6)',           '8e15.8)'/
      data fmti/'(11x,8(i4,4x))',    '(11x,8(i4,4x))',
     &          '(11x,17(1x,i2,1x))','(11x,30(1x,i2,1x))',
     &          '(11x,6(i4,8x))',    '(11x,10(i4,8x))',
     &          '(11x,5(i4,9x))',    '(11x,8(i4,11x))'/
      data fmtt/'(11x,8a8)',         '(11x,8a8)',
     &          '(11x,17a4)',        '(11x,30a4)',
     &          '(11x,6(a8,2x))',    '(11x,10(a8,4x)',
     &          '(11x,5(a8,5x))',    '(11x,8(a8,7x))'/
C
      fmt1   = '(4x,i3,4x,' // fmt(ns)
      fmt2   = '(1x,a8,2x' // fmt(ns)
      fti    = fmti(ns)
      ftt    = fmtt(ns)
      nspalt = nsp(ns)

C
      write(l6,'(//,3x,a)') text
C
      ka = 1
      ke = nspalt
      nteil = n2/nspalt
      if (nteil*nspalt.ne.n2) nteil = nteil + 1
C
      do  10  nt = 1,nteil
      if (n2.gt.nspalt)  write(L6,100)  nt
  100 format(//, 10x,'Part',i5,' of the Matrix',/)
      if(nt.eq.nteil) ke = n2
      if (it.lt.3) then
        write(L6,fti) (k,k=ka,ke)
      else
        write(L6,ftt) (t2(k),k=ka,ke)
      endif
C
      do 20  i=1,n1
         kee=ke
         if (is.eq.2.and.ke.gt.i) kee=i
         if (ka.gt.kee) goto 20
         if (is.eq.3) then
            if (it.eq.1) then
               write(L6,fmt1) i,(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
            else
               write(L6,fmt2) t1(i),(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
            endif
         else
            if (it.eq.1) then
               write(L6,fmt1) i,(a(i+(k-1)*ma),k=ka,kee)
            else
               write(L6,fmt2) t1(i),(a(i+(k-1)*ma),k=ka,kee)
            endif
         endif
   20 continue
c
      ka=ka+nspalt
      ke=ke+nspalt
   10 continue
C
      return
C-end-APRINT
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
c======================================================================c
  
      function itestc()
  
c======================================================================c
c
C    yields 1, if interrupted
c           2, if convergence
c
c    the iteration is determined by the parameter XMIX
c    it can be fixed, or automatically adjusted according to
c
c     IAUT = 0 fixed value for XMIX
c          = 1 automatic adjustment of XMIX
c           
c     INXT = 0 with out inputs from the console
c          > 0 after INXT iterations INX and XMIX are read
c
c     INX  = 0 immediate interruption of the iteration
c          > 0 further INX steps with fixed XMIX, which is read  
c          < 0 further ABS(INX) steps with automatic change of XMIX  
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      integer itestc
c
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
c
      if (ii.ge.2.and.si.lt.epsi) then
         itestc = 2
         return
      endif
c
c     change of XMIX by reading from the console:
      if (inxt.eq.ii) then
         write(6,*) 
     &   'next stop? (0 right now, >0 fixed xmix, <0 autom. xmix)'
         read(*,*)  inx
         if (inx.eq.0) then
            itestc = 1
            return
         endif
         if (inx.lt.0) then
            iaut = 1
         endif
         if (inx.gt.0) then
            iaut = 0
         endif
         inxt = ii+iabs(inx)
         write(6,*) 'new value for xmix?'
         read(*,*) xmix
         write(6,*) inxt,xmix
         xmix0 = xmix
      endif
c
c     automatic change of XMIX:
      if ((si.lt.siold).and.iaut.eq.1) THEN
         xmix = xmix * 1.04
         if (xmix.gt.xmax) xmix = xmax
      else
         xmix = xmix0
      endif
      siold  = si
      itestc = 0
c
      return
c-end-ITESTC
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
      dimension  a(ma,1),x(mx,m)
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
c=======================================================================
  
      subroutine ordi(n,e,mu)
  
c=======================================================================
c     
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n),mu(n)
c  
      do 10 i = 1,n
         k  = i 
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n 
               if (e(j).lt.p) then 
                  k = j 
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
               mk    = mu(k)
               mu(k) = mu(i)
               mu(i) = mk
            endif
         endif
   10 continue
c
      return
c-end-ORDI
      end 
C=======================================================================

      subroutine sdiag(nmax,n,a,d,x,e,is)

C=======================================================================
C
C     A   matrix to be diagonalized
C     D   eigenvalues    
C     X   eigenvectors
C     E   auxiliary field
C     IS = 1  eigenvalues are ordered and major component of X is positiv
C          0  eigenvalues are not ordered            
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension a(nmax,nmax),x(nmax,nmax),e(n),d(n)
C
      data tol,eps/1.e-32,1.e-10/                           
C
      if (n.eq.1) then
         d(1)=a(1,1)  
         x(1,1)=1.
         return
      endif
c
      do 10 i=1,n 
      do 10 j=1,i 
   10    x(i,j)=a(i,j)
c
ccc   householder-reduktion
      i=n
   15 if (i-2) 200,20,20
   20 l=i-2
      f=x(i,i-1)
      g=f            
      h=0  
      if (l) 31,31,32
   32 do 30 k=1,l
   30 h=h+x(i,k)*x(i,k)
   31 s=h+f*f         
      if (s-tol) 33,34,34              
   33 h=0                             
      goto 100                       
   34 if (h) 100,100,40             
   40 l=l+1                        
      g= dsqrt(s)
      if (f.ge.0.) g=-g        
      h=s-f*g                 
      hi=1.d0/h                
      x(i,i-1)=f-g          
      f=0.0                 
      if (l) 51,51,52     
   52 do 50 j=1,l        
      x(j,i)=x(i,j)*hi  
      s=0.0             
      do 55 k=1,j     
   55 s=s+x(j,k)*x(i,k)                      
      j1=j+1                                
      if (l-j1) 57,58,58                   
   58 do 59 k=j1,l                        
   59 s=s+x(k,j)*x(i,k)                  
   57 e(j)=s*hi                         
   50 f=f+s*x(j,i)                     
   51 f=f*hi*.5d0                      
c                                    
      if (l) 100,100,62             
   62 do 60 j=1,l                  
      s=x(i,j)                    
      e(j)=e(j)-f*s              
      p=e(j)                    
      do 65 k=1,j              
   65 x(j,k)=x(j,k)-s*e(k)-x(i,k)*p        
   60 continue                            
  100 continue                           
      d(i)=h                            
      e(i-1)=g                         
      i=i-1                           
      goto 15            
c            
ccc   Bereitstellen der Transformationmatrix 
  200 d(1)=0.0                               
      e(n)=0.0                              
      b=0.0                                
      f=0.0                               
      do 210 i=1,n                      
      l=i-1                            
      if (d(i).eq.0.) goto 221        
      if (l) 221,221,222             
  222 do 220 j=1,l                  
      s=0.0                         
      do 225 k=1,l                
  225 s=s+x(i,k)*x(k,j)          
      do 226 k=1,l              
  226 x(k,j)=x(k,j)-s*x(k,i)   
  220 continue                
  221 d(i)=x(i,i)            
      x(i,i)=1              
      if (l) 210,210,232   
  232 do 230 j=1,l        
      x(i,j)=0.0          
  230 x(j,i)=0.0         
  210 continue         
c
ccc   Diagonalisieren der Tri-Diagonal-Matrix
      DO 300 L=1,N                     
      h=eps*( abs(d(l))+ abs(e(l)))
      if (h.gt.b) b=h             
c
ccc   Test fuer Splitting        
      do 310 j=l,n              
      if ( abs(e(j)).le.b) goto 320
  310 continue                 
c
ccc   test fuer konvergenz    
  320 if (j.eq.l) goto 300   
  340 p=(d(l+1)-d(l))/(2*e(l))          
      r= dsqrt(p*p+1.d0)
      pr=p+r                           
      if (p.lt.0.) pr=p-r             
      h=d(l)-e(l)/pr                 
      do 350 i=l,n                  
  350 d(i)=d(i)-h                  
      f=f+h                       
c
ccc   QR-transformation          
      p=d(j)                    
      c=1.d0                     
      s=0.0                    
      i=j                    
  360 i=i-1                 
      if (i.lt.l) goto 362 
      g=c*e(i)            
      h=c*p              
      if ( abs(p)- abs(e(i))) 363,364,364
  364 c=e(i)/p                          
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*p*r                     
      s=c/r                           
      c=1.d0/r                         
      goto 365                      
  363 c=p/e(i)                     
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*e(i)*r             
      s=1.d0/r                      
      c=c/r                     
  365 p=c*d(i)-s*g             
      d(i+1)=h+s*(c*g+s*d(i)) 
      do 368 k=1,n           
         h=x(k,i+1)            
         x(k,i+1)=x(k,i)*s+h*c
  368    x(k,i)=x(k,i)*c-h*s 
      goto 360           
  362 e(l)=s*p          
      d(l)=c*p         
      if ( abs(e(l)).gt.b) goto 340
c
ccc   konvergenz      
  300 d(l)=d(l)+f    
c
      if (is.eq.0) return
ccc   ordnen der eigenwerte    
      do 400 i=1,n            
      k=i                    
      p=d(i)                
      j1=i+1               
      if (j1-n) 401,401,400   
  401 do 410 j=j1,n          
      if (d(j).ge.p) goto 410 
      k=j                    
      p=d(j)                
  410 continue             
  420 if (k.eq.i) goto 400
      d(k)=d(i)          
      d(i)=p            
      do 425 j=1,n     
      p=x(j,i)        
      x(j,i)=x(j,k)  
  425 x(j,k)=p      
  400 continue     
c                 
c     signum
      do  71 k=1,n
      s=0.0
      do 72 i=1,n
      h= abs(x(i,k))
      if (h.gt.s) then
         s=h
         im=i
      endif
   72 continue
      if (x(im,k).lt.0.0) then
         do 73 i=1,n
   73    x(i,k)=-x(i,k)
      endiF
   71 continue
c 
      return
c-end-SDIAG
      end 
c======================================================================c

      subroutine spline(x,y,y2,n,yp0,ypn)

c======================================================================c
c
c     SPLINE-Routine of "Numerical Recipies" p.88
c
c input:
c     X,Y       tabulated function (0..N)
c     YP0,YPN   first derivatives at point 0 and N 
c               of the interpolating spline-function
c               (if larger then 1.e30, natural spline: y''= 0)   
c output:
c     Y2        second derivatives of the interpolating function
c               is used as input for function SPLINT
c
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      parameter (nmax=500)
      dimension x(0:n),y(0:n),y2(0:n),u(0:nmax)
c
      if (nmax.lt.n) stop ' in SPLINE: nmax too small'
      if (yp0.gt.999d30) then
         y2(0) = 0.0
         u(0)  = 0.0
      else
         y2(0) = -0.5d0
         u(0)  = (3.d0/(x(1)-x(0))) * ((y(1)-y(0))/(x(1)-x(0))-yp0)
      endif
      do 11 i = 1,n-1
         sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p     = sig*y2(i-1) + 2.d0
         y2(i) = (sig - 1.d0)/p
         u(i)  = (6.d0*( (y(i+1)-y(i))/(x(i+1)-x(i)) -
     &                   (y(i)-y(i-1))/(x(i)-x(i-1)) )/
     &                   (x(i+1)-x(i-1)) - sig*u(i-1))/p
   11    continue
      if (ypn.gt..999d30) then
         qn = 0.0
         un = 0.0
      else
         qn = 0.5d0
         un = (3.d0/(x(n)-x(n-1))) * (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k = n-1,0,-1
         y2(k) = y2(k)*y2(k+1)+u(k)
   12 continue
c
      return
c-end-SPLINE
      end
c======================================================================c

      subroutine splint(is,xa,ya,y2a,n,x,y,y1,y2)

c======================================================================c
c
c     SPLINT-Routine of "Numerical Recipies" p.89
c
c input:
c     XA,YA     tabulated function (0:N)
c     Y2A       first derivatives (output von SPLINE) 
c     X         given value on the abscissa
c
c output:
c   is = 0:  Y  interpolated value  Y(x)
c   is = 1;  y1 in addition interpolated value of the derivativ dY/dx  
c   is = 2;  y2 in addition interpolated value of 2nd derivativ d2Y/dx2
c  
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      dimension xa(0:n),ya(0:n),y2a(0:n)
      data sixth/0.1666666666666666667d0/
c
      klo = 0
      khi = n
    1 if (khi-klo.gt.1) then
         k = (khi+klo)/2
         if (xa(k).gt.x) then
            khi = k
         else 
            klo = k
         endif
         goto 1
      endif
      h = xa(khi)-xa(klo)
      if (h.eq.0.0) pause ' in SPINT: bad xa input '
      hi = 1.d0/h
      a  = (xa(khi)-x)*hi
      b  = (x-xa(klo))*hi
c
c     value of the function
      y = a*ya(klo)+b*ya(khi)+
     &    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)*sixth
c 
c     first derivative 
      if (is.lt.1) return
      y1 = hi*(-ya(klo)+ya(khi)) + 
     &    (-(3*a**2-1)*y2a(klo)+(3*b**2-1)*y2a(khi))*h*sixth
c
c     second derivative
      if (is.lt.2) return
      y2 = a*y2a(klo) + b*y2a(khi)
c
      return
c-end-SPLINT
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
      
C      common /basnnn/ n0f,n0b
C      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
C      common /bloosc/ ia(nbx,2),id(nbx,2)
C      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
C      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)
C      common /bosqua/ no
C      common /gfviv / iv(0:igfv)
C      common /sdimos/ nrm,nlm,nrbm
C      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
C      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
C      common /vvvikf/ mv,ipos(nbx),nib(mvx),nni(2,mvx)





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

      subroutine densit(lpr)

c=====================================================================c
C
C     density at the radius r = xh(ih)*b0 is given by
C     b0**(-3) * rv(ih) / ( 4*pi * r**2 ) in units of fm**(-3)
C
c---------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor
c
      dimension frs(2),frv(2),frr(2)
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)
      common /gaucor/ rb(0:ngh),wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),jmax
      common /optopt/ icm,icou,it1,it2,ncut
      common /rados1/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /rados2/ rnl2(1:nrx,0:nlx,0:ngh),rnl12(1:nrx,0:nlx,0:ngh)
      common /rhoshe/ rrs(nqx,nb2x),rrv(nqx,nb2x),rrr(nqx,nb2x)
      common /rhorho/ rs(0:ngh,1:2),rv(0:ngh,1:2),dro(0:ngh)
      common /rhoro2/ rs2(0:ngh,1:2),rv2(0:ngh,1:2),rr2(0:ngh,1:2),
     &                dr(0:ngh,1:2),drold(0:ngh,1:2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
c
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN DENSIT *******************************'
c
c
CL----09/04/99
c      drold = dr
CL----
      do it = it1,it2
         do ih = 0,ngh
            dro(ih)    = zero
            rs(ih,it)  = zero
            rv(ih,it)  = zero
            rs2(ih,it) = zero
            rv2(ih,it) = zero
            rr2(ih,it) = zero
CL----09/04/99
            drold(ih,it) = dr(ih,it)
CL----
CL  25/03/99 
            dr(ih,it)  = zero
         enddo
      enddo
c
c     loop over j-blocks
      do 10 ib = 1,nb
         ibg = ib - 1 + 2*mod(ib,2)
         j   = ijb(ib)
         l   = ilb(ib)
         lg  = ilb(ibg)
         nd  = id(ib)
         ll  = l*(l+1)
c
         do 20 n2 =  1,nd
         do 20 n1 = n2,nd
            i12 = 2 - n2/n1
c
CL--------------------------------------------------------------
CL             frs, frv e frr in units of fm**(-3)
CL--------------------------------------------------------------
            do it = it1,it2
               frs(it) = rrs(n1+(n2-1)*nd,ib+(it-1)*nbx)*i12
               frv(it) = rrv(n1+(n2-1)*nd,ib+(it-1)*nbx)*i12
               frr(it) = rrr(n1+(n2-1)*nd,ib+(it-1)*nbx)*i12
            enddo   
c
            nn = 2*(n1+n2+l)-1
            do ih = 0,ngh
               r  = xh(ih)*b0c
               s  = rnl(n1,l,ih)*rnl(n2,l,ih)
               s2 = rnl2(n1,l,ih)*rnl2(n2,l,ih)
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!! The sign minus bellow take into account the difference       !! 
C !!!!!!! between convection Ring's and Serot&Walecka's Dirac spinors  !!       
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               s3 = -( 2*rnl2(n1,l,ih)*rnl2(n2,lg,ih)/r
     &            +  ( rnl12(n1,l,ih)*rnl2(n2,lg,ih)
     &            +    rnl2(n1,l,ih)*rnl12(n2,lg,ih) )/b0c )
CL  25/03/99
               s4 = -rnl(n1,l,ih)*rnl(n2,lg,ih)
CL
               do it = it1,it2
                   rs(ih,it) = rs(ih,it) + frs(it)*s
                   rv(ih,it) = rv(ih,it) + frv(it)*s
                  rs2(ih,it) = rs2(ih,it) + frs(it)*s2
                  rv2(ih,it) = rv2(ih,it) + frv(it)*s2
                  rr2(ih,it) = rr2(ih,it) + frr(it)*s3
CL------------ 25/03/99
CL------------ Delta-rho for extracting contact term
CL                dr in units of fm**(-3)
CL
                  dr(ih,it) = dr(ih,it) + frr(it)*s4
CL------------
               enddo
c
c------------- Delta-rho for calculation of Coulomb field
               if (icou.gt.0) then
                  xx = xh(ih)*xh(ih)
                  dro(ih) = dro(ih) + 2*frv(2) * ( s*(xx+ll/xx-nn) 
     &                              + rnl1(n1,l,ih)*rnl1(n2,l,ih)) 
               endif
            enddo
   20    continue
   10 continue
c
c
c---- check, whether integral over dro vanishes
      s = zero
      do ih = 1,ngh
         s = s + dro(ih)
      enddo
      if (lpr) write(l6,*) 'integral over dro',s
c
c
c---- normalization and renormalization to particle number
      do it = it1,it2
         s  = zero
         do ih = 1,ngh
            s  =  s + rv(ih,it)
         enddo
CL------------------------------------------------------
CL        s in units of fm**(-3)
CL------------------------------------------------------
         if (lpr) write(l6,'(a,i3,2f15.8)') 
     &                  ' norm of the vector density = ',it,s
         s = npr(it)/s                                  ! [s]=fm**3
         bi2 = one/(b0f*b0f)
         do ih = 0,ngh
            f  = s/wdcor(ih)                            ! [f]=1
            rs(ih,it)   = f*rs(ih,it)
            rv(ih,it)   = f*rv(ih,it)
            rs2(ih,it)  = s*rs2(ih,it)
            rv2(ih,it)  = s*rv2(ih,it)
            rr2(ih,it)  = s*rr2(ih,it)
CL---------- 25/03/99
            dr(ih,it)   = f*dr(ih,it)                   ! [dr]=fm**(-3)
CL----------
            if (it.eq.2) dro(ih) = f*dro(ih)*bi2
         enddo
         if (lpr) then
            call prigh(0,ro,b0f,'x(fm) ')
            call prigh(1,rs(0,it),one,'ROS '//tis(it))
            call prigh(1,rv(0,it),one,'ROV '//tis(it))
            if (it.eq.2) call prigh(1,dro,one,'DRO  ')
            call prigh(0,ro,b0c,'x(fm) ')
            f = one/(4.d0*pi*b0c**3)
            call prigh(2,rs2(0,it),f,'ROS2'//tis(it))
            call prigh(2,rv2(0,it),f,'ROV2'//tis(it))
         endif
      enddo
c
c
      if (lpr)
     &write(l6,*) ' ****** END DENSIT *********************************'
      return
C-end-DENSIT
      end
c======================================================================c

      subroutine denssh(lpr)

c======================================================================c
c
c     calculates densities in oscillator basis 
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
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ icm,icou,it1,it2,ncut
      common /rhoshe/ rrs(nqx,nb2x),rrv(nqx,nb2x),rrr(nqx,nb2x)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /texblo/ tb(ntx),txb(nbx)
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
      common /wavefg/ fg(nq2x,nb2x)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN DENSSH *******************************'
c
c---- loop over the j-blocks
      do 10 ib = 1,nb
         ibg = ib - 1 + 2*mod(ib,2)
         nf  = id(ib)
         ng  = id(ibg)
         nd  = nf + ng
         imf = ia(ib)
         img = ia(ibg)
c         nfg = nf
c         if (ng.lt.nf) nfg = ng
c
c------- loop over neutron and proton
         do 20 it = it1,it2
            mf  = ib  + (it-1)*nbx
            mg  = ibg + (it-1)*nbx
c
c---------- loop over the n-quantum numbers
            do 30 n2 =  1,nf
            do 30 n1 = n2,nf
               i12 = 2 - n2/n1
c
               sf = zero
               do k = 1,nf
                  sf = sf + fg(n1+(k-1)*nd,mf)*fg(n2+(k-1)*nd,mf) 
     &                    * vv(imf+k,it)
               enddo
               sg = zero
               do k = 1,ng
                  sg = sg + fg(ng+n1+(k-1)*nd,mg)*fg(ng+n2+(k-1)*nd,mg) 
     &                    * vv(img+k,it)
               enddo
C                             
C
CL------------------------------------------
CL            sfg in units of fm**(-3)
CL------------------------------------------
               sfg = zero
               do k = 1,nf
                  sfg = sfg + fg(n1+(k-1)*nd,mf)*fg(nf+n2+(k-1)*nd,mf) 
     &                      * vv(imf+k,it)
               enddo
CL-----------------------------------------------
CL            rrs, rrv e rrr in units of fm**(-3)
CL-----------------------------------------------            
               rrs(n1+(n2-1)*nf,mf) = sf-sg
               rrv(n1+(n2-1)*nf,mf) = sf+sg
               rrr(n1+(n2-1)*nf,mf) = 2*sfg
   30       continue
c
            if (lpr) then
               write(l6,'(/,a,1x,a)') txb(ib),tis(it)
               call aprint(2,2,1,nf,nf,nf,rrs(1,mf),
     &                     tb(imf+1),' ','RRS ')
               call aprint(2,2,1,nf,nf,nf,rrv(1,mf),
     &                     tb(imf+1),' ','RRV ')
            endif
   20    continue
   10 continue
c
      if (lpr)
     &write(l6,*) ' ****** END DENSSH *********************************'
      return
c-end-DENSSH
      end      
c======================================================================c

      subroutine dirac(lpr)

c======================================================================c
c
c     solves the Dirac-Equation in spherical oscillator basis c
c     units:    fields and Hamiltonian in fm^(-1)
c               eigenvalues in MeV
c 
c----------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical  lpr
      character*8 tb,tbb(nd2x)
      character*25 txb
c
      dimension e(nd2x),ez(nd2x)
      dimension hh(nq4x)
c
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /dimens/ n0f,n0b,nrm,nlm
      common /eeeeee/ ee(ntx,2),v2(ntx,2),mu(ntx)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ icm,icou,it1,it2,ncut
      common /physco/ amu,hqc,alphi,r0
      common /potpot/ vps(0:ngh,1:2),vms(0:ngh,1:2)
      common /single/ sp(nq2x,nb2x)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /texblo/ tb(ntx),txb(nbx)
      common /wavefg/ fg(nq2x,nb2x)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN DIRAC ********************************'
c
      emcc2 = 2*amu/hqc
c
c
c     loop over protons and neutrons 
C-----------------------------------
      do 10 it = it1,it2
c
c
c------- loop over the different j-blocks
         do 30 ib = 1,nb
            ibg = ib - 1 + 2*mod(ib,2)
            nf  = id(ib)
            ng  = id(ibg)
            nd  = nf + ng
            imf = ia(ib)
            img = ia(ibg)
            lf  = ilb(ib)
            lg  = ilb(ibg)
            m   = ib + (it-1)*nbx
c
c           calculation of the Dirac-Matrix:
C-------------------------------------------
            do i2 = 1,nf
            do i1 = 1,ng
c               hh(nf+i1+(i2-1)*nd) = sp(i1+(i2-1)*ng,ib)
               hh(nf+i1+(i2-1)*nd) = sp(i1+(i2-1)*ng,m)
            enddo
            enddo
            call pot(nf,nd,lf,vps(1,it),hh)
            call pot(ng,nd,lg,vms(1,it),hh(nf+1+nf*nd))
            do i = nf+1,nd
               hh(i+(i-1)*nd) = hh(i+(i-1)*nd) - emcc2
            enddo
c
c           cut off large components with highest N
            if (2*(nf-1)+lf.gt.n0f) then
               do i = 1,nf
                  hh(nf+(i-1)*nd) = zero
               enddo
               do i = nf+1,nd
                  hh(i+(nf-1)*nd) = zero
               enddo
               hh(nf+(nf-1)*nd) = 1000.0
            endif 
            if (lpr) then
               do i = 1,nf
                  tbb(i) = tb(imf+i)
               enddo
               do i = 1,ng
                  tbb(nf+i) = tb(img+i)
               enddo
               write(l6,'(/,a)') txb(ib)
               call aprint(2,2,1,nd,nd,nd,hh,tbb,' ','HH') 
            endif
c
c---------- Diagonalization:
            call sdiag(nd,nd,hh,e,hh,ez,+1)
c           call aprint(1,1,1,1,1,nd,e,' ',' ','E')
c           call aprint(1,2,1,nd,nd,nd,hh,tbb,' ','XX')
            do k = 1,nf
               ee(imf+k,it) = e(ng+k)*hqc 
               do i = 1,nd
                  fg(i+(k-1)*nd,m) = hh(i+(ng+k-1)*nd)
               enddo
            enddo
            if (lpr) then
               call aprint(1,1,1,1,1,nf,ee(imf+1,it),' ',' ','E')
               call aprint(1,2,1,nd,nd,nf,fg(1,m),tbb,' ','FG')
            endif
   30    continue
   10 continue
c
      if (lpr)
     &write(l6,*) ' ****** END DIRAC **********************************'
      return
C-end-DIRAC
      end
c======================================================================c

      subroutine pot(n,nd,l,v,tt)

c======================================================================c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
c
      dimension tt(nd,nd),v(ngh)
c
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /rados2/ rnl2(1:nrx,0:nlx,0:ngh),rnl12(1:nrx,0:nlx,0:ngh)
c
      do 10 n2 = 1,n
      do 10 n1 = n2,n
         s = 0.0
         do 20 ih = 1,ngh
            s = s + v(ih)*rnl2(n1,l,ih)*rnl2(n2,l,ih)
   20    continue
         tt(n1,n2) = s
   10 continue
c
      return
c-end-POT
      end
c======================================================================c

      subroutine expect(lpr)

c======================================================================c
c
c     calculates expectation values
c
c----------------------------------------------------------------------c
      include 'paramet'
      parameter (ndwork = nwork-4*ngh-4)
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      dimension xn(3),r2(3)
      dimension ept(3),ekt(3),epart(3),entro(3),vir(3)
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)
      common /erwar / ea,rms,qp
      common /fermi / ala(2),tz(2)
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /gaucor/ rb(0:ngh),wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /nucnuc/ amas,nama,npr(2),jmax 
      common /optopt/ icm,icou,it1,it2,ncut
      common /pair  / ga(2),gg(2),del(2),spk(2),dec(2),pwi
      common /physco/ amu,hqc,alphi,r0
      common /potpot/ vps(0:ngh,1:2),vms(0:ngh,1:2)
      common /quaqua/ nt,nr(ntx),nl(ntx),nj(ntx)
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2),dro(0:ngh)
      common /rhoro2/ rs2(0:ngh,1:2),rv2(0:ngh,1:2),rr2(0:ngh,1:2),
     &                dr(0:ngh,1:2),drold(0:ngh,1:2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /temper/ temp
      common /work  / dsig(0:ngh),dome(0:ngh),drho(0:ngh),dcou(0:ngh),
     &                dwork(ndwork)
c     
      !print *, 'Factorial of 5 is', fak(5)
      if (lpr)
     &write(l6,*) ' ****** BEGIN EXPECT *******************************'
c     
c
      do it = 1,3
         xn(it) = zero
         r2(it) = zero
      enddo
c
c---- particle number and radii
      do ih = 1,ngh
         r  = rb(ih)
         rr = r*r
         wx = wdcor(ih)
         do it = it1,it2
            x      = wx*rv(ih,it)
            xn(it) = xn(it) + x
            r2(it) = r2(it) + x*rr
         enddo
      enddo
c
      do it = it1,it2
         r2(it) = sqrt(r2(it)/xn(it))
      enddo
c
      xn(3) = xn(1) + xn(2)
      r2(3) = sqrt((npr(1)*r2(1)**2+npr(2)*r2(2)**2)/amas)
      rc    = sqrt(r2(2)**2 + 0.64)
      rms   = r2(3)
c
c
c---- pairing-energy
      do it = it1,it2
         del(it) = gg(it)*spk(it) + dec(it)
         ept(it) = - del(it)*spk(it)
      enddo
      ept(3) = ept(1) + ept(2)
c
c
      do it = it1,it2
c
c------- single particle energy
         epart(it) = zero
         do i = 1,nt
            epart(it) = epart(it) + ee(i,it)*vv(i,it)
         enddo
c
c------- entropy
	   entro(it) = zero
         if (temp.gt.zero) then
            do i = 1,nt
	           fn = vv(i,it)/mu(i)
	           if (fn.gt.zero) entro(it) = entro(it)
     %	           - mu(i)*fn*dlog(fn)
	           fn = one - fn
	           if (fn.gt.zero) entro(it) = entro(it)
     %	            - mu(i)*fn*dlog(fn)
            enddo
         endif
      enddo
      epart(3) = epart(1) + epart(2)
      entro(3) = entro(1) + entro(2)
c
c---- field energies
      b2    = g2/3
      b3    = g3/2
      bo3   = w3/2
      esig  = zero
      eome  = zero
      erho  = zero
      ecou  = zero
      ecoex = zero
      esnl   = zero
      eonl   = zero
      do ih = 1,ngh
         wx   = wdcor(ih)
         esig = esig - sig(ih)*( rs2(ih,1)+rs2(ih,2))
         eome = eome - ome(ih)*( rv2(ih,1)+rv2(ih,2))
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!! rr2 is in fm^-1                                        !!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         erho = erho - rho(ih)*(-rv2(ih,1)+rv2(ih,2))
     +               - rho(ih)*(-rr2(ih,1)+rr2(ih,2))
     +               *  fr*hqc/2/amu/grho
         if (icou.ge.1) ecou = ecou - cou(ih)*rv2(ih,2)
         if (icou.eq.2) ecoex = ecoex + rv(ih,2)**(4*third)*wx
         esnl  = esnl  - (b2*sig(ih)**3+b3*sig(ih)**4)*wx*cf**3
         eonl  = eonl  + (bo3*ome(ih)**4)*wx*cf**3
      enddo
      esig  = hqc*gsig*esig/2
      eome  = hqc*gome*eome/2
      erho  = hqc*grho*erho/2
      ecou  = hqc*ecou/2
      ecoex =  -(0.75d0*hqc/alphi)*(3/pi)**third*ecoex
      esnl   = hqc * esnl/2
      eonl   = hqc * eonl/2
c
      ecm  = -0.75d0*hom
c
      ekt(1) = zero
      ekt(2) = zero
      ekt(3) = epart(3) + two*(esig + eome + erho + ecou)
c
      etot = epart(3) + esig + eome + erho + ecou + ecoex +
     &                  esnl + eonl + ept(3) + ecm
      ea = etot/amas
c
c---- virial theorem: J.Rafelski; Phys.Rev.D16(1977)1890
c
c---- interpolate the fields for calculation of derivative
c     call spline(rb,sig,yp,ngh,zero,zero)
c     dsig(0) = zero
c     do ih = 1,ngh
c        call splint(1,xx,yy,yp,ngh,rb(ih),y,dsig(ih),y2)
c     enddo
c     call spline(rb,ome,yp,ngh,zero,zero)
c     dome(0) = zero
c     do ih = 1,ngh
c        call splint(1,xx,yy,yp,ngh,rb(ih),y,dome(ih),y2)
c     enddo
c     call spline(rb,rho,yp,ngh,zero,zero)
c     drho(0) = zero
c     do ih = 1,ngh
c        call splint(1,xx,yy,yp,ngh,rb(ih),y,drho(ih),y2)
c     enddo
c     yp2 = -npr(2)/(alphi*rb(ngh)**2)
c     call spline(rb,cou,yp,ngh,zero,yp2)
c     dcou(0) = zero
c     do ih = 1,ngh
c        call splint(1,xx,yy,yp,ngh,rb(ih),y,dcou(ih),y2)
c     enddo
c---- integration
c     vsig = zero
c     vome = zero
c     vrho = zero
c     vcou = zero
c     do ih = 1,ngh
c        wx = rb(ih)*wdcor(ih)
c        vsig = vsig + dsig(ih)*(rs(ih,1)+rs(ih,2))*wx
c        vome = vome + dome(ih)*(rv(ih,1)+rv(ih,2))*wx
c        vrho = vrho + drho(ih)*(-rv(ih,1)+rv(ih,2))*wx
c        vcou = vcou + dcou(ih)*rv(ih,2)*wx
c     enddo
c     vsig = hqc*gsig*vsig
c     vome = hqc*gome*vome
c     vrho = hqc*grho*vrho
c     vcou = hqc*vcou
c     vir(3) = - ekt(3) + vsig + vome + vrho + vcou
c
c      call csplin(ngh,rb(1),sig(1),bb,cc,dd)
c      do ip = 1,ngh
c	 uu = rb(ip)
c	 call cseval(ngh,uu,1,rb(1),sig(1),bb,cc,dd,dsig(ip))
c      enddo 
c      call csplin(ngh,rb(1),ome(1),bb,cc,dd) 
c      do ip = 1,ngh
c         uu = rb(ip)
c	 call cseval(ngh,uu,1,rb(1),ome(1),bb,cc,dd,dome(ip))
c      enddo 
c      call csplin(ngh,rb(1),rho(1),bb,cc,dd)
c      do ip = 1,ngh
c	 uu = rb(ip)
c	 call cseval(ngh,uu,1,rb(1),rho(1),bb,cc,dd,drho(ip))
c      enddo
c      call csplin(ngh,rb(1),cou(1),bb,cc,dd)
c      do ip = 1,ngh
c	 uu = rb(ip)
c	 call cseval(ngh,uu,1,rb(1),cou(1),bb,cc,dd,dcou(ip))
c      enddo 
c---- calculate the integral
c      vm   = zero
c      vsig = zero
c      vome = zero
c      vrho = zero
c      vcou = zero
c      do ip = 1,ngh
c	 wx = wdcor(ip)
c	 vm = vm + amu*(rs(ip,1)+rs(ip,2))*wx
c	 vsig = vsig + rb(ip) *gsig*dsig(ip)*(rs(ip,1)+rs(ip,2))*wx
c	 vome = vome + rb(ip) *gome*dome(ip)*(rv(ip,1)+rv(ip,2))*wx
c	 vrho = vrho + rb(ip) *grho*drho(ip)*(-rv(ip,1)+rv(ip,2))*wx
c	 vcou = vcou + rb(ip)*      dcou(ip)*rv(ip,2)*wx
c      enddo
c      vsig=hqc*vsig
c      vome=hqc*vome
c      vrho=hqc*vrho
c      vcou=hqc*vcou
c
c---- the left term tl
c      tl=epartv(3)-vm+2.0*esig+2.0*eome+2.0*erho+2.0*(ecou+ecoex)
c      write(l6,'(a,5x,F10.3)') 'epartv(3) ',epartv(3)
c      write(l6,'(a,5x,F10.3)') 'vm       ',vm
c      write(l6,'(a,5x,F10.3)') '2.0*esig ',2.0*esig
c      write(l6,'(a,5x,F10.3)') '2.0*eome ',2.0*eome
c      write(l6,'(a,5x,F10.3)') '2.0*erho ',2.0*erho
c      write(l6,'(a,5x,F10.3)') '2.0*(ecou+ecoex) ',2.0*(ecou+ecoex)
c---- the right term tr 
c      tr= vsig + vome + vrho + vcou
c      write(l6,'(a,5x,F10.3)') 'vsig   ',vsig
c      write(l6,'(a,5x,F10.3)') 'vome   ',vome
c      write(l6,'(a,5x,F10.3)') 'vrho   ',vrho
c      write(l6,'(a,5x,F10.3)') 'vcou   ',vcou 
c---- the virial theorem
c      tvir = tl - tr
    
c
c---- printout
      if (.not.lpr) return
c
      write(l6,'(/,28x,a,8x,a,9x,a)') 'neutron','proton','total'
c
c     particle number
      write(l6,'(a,6x,3f15.6)') ' particle number',xn
c
c     Lambda
      write(l6,'(a,7x,2f15.6)') ' lambda        ',ala
c
c     Delta
      write(l6,'(a,7x,2f15.6)') ' Delta         ',del
c
c     trace of kappa
      write(l6,'(a,7x,2f15.6)') ' spk           ',spk
c
c     rms-Radius    
      write(l6,'(a,7x,3f15.6)') ' rms-Radius    ',r2
c
c     charge-Radius    
      write(l6,'(a,22x,f15.6)') ' charge-Radius ',rc
c
c
c     kinetic energy
      write(l6,'(/,a,7x,3f15.6)') ' Kinetic Energy',ekt
c
c     virial theorem
      write(l6,'(a,22x,f15.6)') ' Virial Theorem',vir(3)
c
c     single-particle energy
      write(l6,'(/,a,7x,3f15.6)') ' Particle Energ',epart
c 
c     sigma energy 
      write(l6,'(a,37x,f15.6)') ' E-sigma       ',esig 
c 
c     nonlinear part sigma energy 
      write(l6,'(a,37x,f15.6)') ' E-sigma n.lin ',esnl 
c 
c     omega energy  
      write(l6,'(a,37x,f15.6)') ' E-omega       ',eome
c 
c     nonlinear part omega energy 
      write(l6,'(a,37x,f15.6)') ' E-omega n.lin ',eonl 
c 
c     rho-energy       
      write(l6,'(a,37x,f15.6)') ' E-rho         ',erho 
c
c     Coulomb energy (direct part)
      write(l6,'(a,37x,f15.6)') ' Coulomb direct',ecou 
c
c     Coulomb energy (exchange part in Slater approx.)
      write(l6,'(a,37x,f15.6)') ' Coulomb exch. ',ecoex 
c
c     pairing energy
      write(l6,'(a,7x,3f15.6)') ' Pairing Energy',ept
c
c     center of mass correction
      write(l6,'(a,37x,f15.6)') ' E-cm          ',ecm
c
c     total energy
      write(l6,'(a,37x,f15.6)') ' Total Energy  ',etot
c
c     energy per particle
      write(l6,'(a,37x,f15.6)') ' E/A           ',ea 
c
c     entropy 
      write(l6,'(/,a,7x,3f15.6)') ' Entropy       ',entro
c
      write(l6,*) ' ****** END EXPECT *********************************'
      return
c-end-EXPECT
      end
c======================================================================c

      subroutine field(lpr)

c======================================================================c
c
c     calculation of the meson-fields in the oscillator basis
c     the fields are given in (fm^-1)
c
c----------------------------------------------------------------------c
      include 'paramet'
      parameter (ndwork = nwork-3*ngh-3)
c
      implicit real*8 (a-h,o-z)
      logical lpr
      dimension pn(nox)
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /gaucor/ rb(0:ngh),wdcor(0:ngh)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /optopt/ icm,icou,it1,it2,ncut
      common /physco/ amu,hqc,alphi,r0
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2),dro(0:ngh)
      common /rhoro2/ rs2(0:ngh,1:2),rv2(0:ngh,1:2),rr2(0:ngh,1:2),
     &                dr(0:ngh,1:2),drold(0:ngh,1:2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /work  / so(0:ngh),phi(0:ngh),sig1(0:ngh),dwork(ndwork)
c
      data maxs/50/
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN FIELD ********************************'
c
c
c---- Mixing of new and old fields: 
      xmi = xmix
      si  = zero
c
c
c---- sigma-meson
      ami2 = (hqc/amsig)**2
      f    = -gsig*ami2
      sb   = -g2*ami2
      sc   = -g3*ami2
      epss = epsi*0.01d0
      do ih = 0,ngh
         sig1(ih) = sig(ih)
      enddo
c
c     loop of sigma-iteration
      do is = 0,maxs
         do ih = 0,ngh
            si1 = sig1(ih)
            si2 = si1*si1
            si3 = si2*si1
            w   = wdcor(ih)*cf**3
            so(ih) = f*(rs2(ih,1)+rs2(ih,2)) + w*(sb*si2 + sc*si3)
         enddo
         call gordon(1,so,phi)
         ss = zero
         do ih = 0,ngh
            ss = dmax1(ss,abs(phi(ih)-sig1(ih)))
            sig1(ih) = phi(ih)
         enddo
c        write(6,100) is,ss,(phi(i),i=1,4)
c        write(l6,100) is,ss,(phi(i),i=1,4)
c 100    format(i4,'.te s-iteration:',f10.8,4f10.5)
         if (ss.lt.epss) goto 10
      enddo
      stop ' in FIELD: sigma-iteration has not converged'
   10 do ih = 0,ngh
         sv       = phi(ih) - sig(ih)
         sig(ih)  = sig(ih) + xmi*sv
         si = max(si,abs(sv))
      enddo
c
c
c---- omega-meson
      ami2 = (hqc/amome)**2
      f    = +gome*ami2
      sc   = -w3*ami2
      epso = epsi*1.d-10
      do ih = 0,ngh
         sig1(ih) = ome(ih)
      enddo
c
c     loop of omega-iteration
      do is = 0,maxs
        do ih = 0,ngh
          si1    = sig1(ih)
          si2    = si1*si1
          si3    = si2*si1
         w      = wdcor(ih)*cf**3
         so(ih) = f*(rv2(ih,1)+rv2(ih,2)) + w*sc*si3
        enddo
        call gordon(2,so,phi)
        ss = zero
         do ih = 0,ngh
           ss = dmax1(ss,abs(phi(ih) - sig1(ih)))
           sig1(ih) = phi(ih)
         enddo
         if (ss.lt.epso) goto 20
      enddo
      stop ' in FIELD: omega-iteration has not converged'
   20 do ih = 0,ngh
         sv       = phi(ih) - ome(ih)
         ome(ih)  = ome(ih) + xmi*sv
         si = max(si,abs(sv))
      enddo
c
c
c---- rho-meson
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!! f1 is in fm^3 and rr2 fm^-1                            !!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ami2 = (hqc/amrho)**2
      f    = +grho*ami2
      f1   = +fr*ami2*hqc/2/amu
CL----09/04/99
      df   =  fr*hqc*third/2/amu
      delr = 0
      do ih = 0,ngh
           delr = delr + xh(ih)
      enddo
           delr = delr/ngh
           df   = df*delr/2
CL----
      do ih = 0,ngh
         so(ih)  = f*(-rv2(ih,1)+rv2(ih,2) ) +
     &            f1*(-rr2(ih,1)+rr2(ih,2) )
CL------  31/03/99
         dr(ih,1) = rdens(1,dr(0,1),pn,xh(ih))
         dr(ih,2) = rdens(1,dr(0,2),pn,xh(ih))
CL------
      enddo
      call gordon(3,so,phi)
      do ih = 0,ngh
         sv       = phi(ih) - rho(ih)
         rho(ih)  = rho(ih) + xmi*sv
         si = max(si,abs(sv))
CL------  31/03/99
         sv1      = dr(ih,1) - drold(ih,1)
         dr(ih,1) = dr(ih,1) + xmi*sv1
         sv1      = sv1*df
         si = max(si,abs(sv1))
         sv2      = dr(ih,2) - drold(ih,2)
         dr(ih,2) = dr(ih,2) + xmi*sv2
         sv2      = sv2*df  
         si = max(si,abs(sv2))
CL------
      enddo
C
C
C---- photon
      if (icou.ge.1) then
         call coulom(phi)
         do ih = 0,ngh
            sv       = phi(ih) - cou(ih)
            cou(ih)  = cou(ih) + xmi*sv
            si = max(si,abs(sv))
         enddo
      endif
c
      si = si*hqc
c
c
      if (lpr) then
         call prigh(0,sig,b0c,'x(fm) ')
         call prigh(1,sig,one,'Sigma ')
         call prigh(1,ome,one,'Omega ')
         call prigh(1,rho,one,'Rho   ')
         call prigh(1,cou,one,'Coulom')
      endif
C   
      if (lpr)
     &write(l6,*) ' ****** END FIELD **********************************'
      return
C-end-FIELD
      end
c======================================================================c

      subroutine gordon(imes,so,phi)

c======================================================================c
C
C     SOLUTION OF THE KLEIN-GORDON-EQU. BY EXPANSION IN OSCILLATOR
c     imes: number of the meson
c     so:   source
c     phi:  meson field
c
c----------------------------------------------------------------------c
      include 'paramet'
      parameter (nox2 = nox*nox)
c
      implicit real*8 (a-h,o-z)
c
      dimension so(0:ngh),phi(0:ngh)
      dimension gi(nox,nox),rr(nox,ngh)
C
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /physco/ amu,hqc,alphi,r0
      common /dimens/ n0f,n0b,nrm,nlm
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /greens/ gg(0:ngh,1:ngh,1:3),igreen(3)
      common /gfvsq / sq(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /mespar/ ames(3),gsig,gome,grho,fr,g2,g3,w3
      common /radbos/ rnb(1:nox,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      data igreen/0,0,0/
c
c     write(l6,*) ' ****** BEGIN GORDON *******************************'
c
c
c---- calculation of the Greens function
      if (igreen(imes).eq.0) then
         f  = (hqc/(ames(imes)*b0b))**2
         no = n0b/2 + 1
         do i = 1,no
            do k = 1,no
               gi(i,k) = zero
            enddo
         enddo
         do n = 1,no
            gi(n,n) = one + f*(2*n-half) 
            if (n.lt.no) then
               gi(n,n+1) = f*sq(n)*sqh(n)
               gi(n+1,n) = gi(n,n+1)
            endif
            do ih = 1,ngh
               rr(n,ih)=rnb(n,ih)
            enddo
         enddo
         call lingd(nox,nox,no,ngh,gi,rr,d,ifl)
         f = one/(4*pi*b0b**3)
         do ih = 0,ngh
            do kh = 1,ngh
               s = zero
               do n = 1,no
                  s = s + rnb(n,ih)*rr(n,kh)
               enddo
               gg(ih,kh,imes) = f*s
            enddo
         enddo      
         igreen(imes) = 1
      endif
c
c
c---- multiplication of source with Greens function
      do ih = 0,ngh
         s = zero
         do kh = 1,ngh
            s = s + gg(ih,kh,imes)*so(kh)
         enddo
         phi(ih) = s
      enddo
c
c     write(l6,*) ' ****** END GORDON *********************************'
      return
c-end-GORDON
      end
c======================================================================c
 
      subroutine coulom(cou)

c======================================================================c
c
c     Coulom-Field (direct part)
c
c----------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
c
      dimension cou(0:ngh)
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /coucal/ vc(0:ngh,1:ngh),icacou
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /gaucor/ rb(0:ngh),wdcor(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ amu,hqc,alphi,r0
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2),dro(0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      data icacou/0/
c
c     write(l6,*) ' ****** BEGIN COULOM *******************************'
c
c
c---- calculation of the Coulomb interaction
      if (icacou.eq.0) then
         f = one/(6*alphi)
         do i1 = 0,ngh
            r1 = xh(i1)*b0c
            do i2 = 1,ngh
               r2 = xh(i2)*b0f
               rg = dmax1(r1,r2)
               rk = dmin1(r1,r2)
               vc(i1,i2) = f * ( 3*rg + rk**2/rg) * wdcor(i2) 
            enddo
         enddo
         icacou = 1
      endif
c
c---- calculation of the Coulomb field
      do ih = 0,ngh
         s = zero
         do kh = 1,ngh
             s = s + vc(ih,kh) * dro(kh) 
         enddo
         cou(ih) = s
      enddo
c
c     write(l6,*) ' ****** END COULOM *********************************'
      return
c-end-COULOM
      end      
c======================================================================c

      subroutine inout(is,lpr)

c======================================================================c
c
c     is = 1: reads meson-fields from tape
c          2: writes meson-fields  to tape
c
c----------------------------------------------------------------------c
      include 'paramet'
      parameter (ndwork = nwork-3*ngh-3)
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      character*1 tp,tl,tis
      character*2 nucnam,nucnam1
      character*10 txtfor,txtfor1
c
      dimension ga1(2),gg1(2),del1(2),spk1(2),dec1(2),tz1(2)
      dimension npr1(2)
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /dimens/ n0f,n0b,nrm,nlm
      common /erwar / ea,rms,qp
      common /fermi / ala(2),tz(2)
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /gaucor/ rb(0:ngh),wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /initia/ vin,rin,ain,inin,iplot
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /nucnuc/ amas,nama,npr(2),jmax
      common /optopt/ icm,icou,it1,it2,ncut
      common /pair  / ga(2),gg(2),del(2),spk(2),dec(2),pwi
      common /physco/ amu,hqc,alphi,r0
      common /quaqua/ nt,nr(ntx),nl(ntx),nj(ntx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
      common /work  / xx(0:ngh),yy(0:ngh),yp(0:ngh),dwork(ndwork)
      common /rhoro2/ rs2(0:ngh,1:2),rv2(0:ngh,1:2),rr2(0:ngh,1:2),
     &                dr(0:ngh,1:2),drold(0:ngh,1:2)
c
      if (is.eq.1.and.inin.ne.0) return
c
      write(l6,*) ' ****** BEGIN INOUT ********************************'
c
c
c
c---- reading of meson fields from tape:
c-------------------------------------
      if (is.eq.1) then
         open(lwin,file='del.wel',status='old')
         read(lwin,100) 
     &        nucnam1,nama1,npr1,ngh1,n0f1,n0b1,nb1,nt1
  100    format(1x,a2,8i4)
         read(lwin,'(5x,2f12.6,6x,f12.9,2f12.6)') 
     &        b01,b02,si,ea,rms
         read(lwin,'(a10,5f10.4)') txtfor1,amsig1,amome1,amrho1
         read(lwin,'(10x,7f10.4)') gsig1,gome1,grho1,fr1,g21,g31,w31
         read(lwin,103) ga1,gg1,pwi1
  103    format(10x,5f12.6)
         read(lwin,103) del1,dec1
         read(lwin,103) spk1
         read(lwin,103) ala,tz1
c
         read(lwin,101) sig
         read(lwin,101) ome
         read(lwin,101) rho
         read(lwin,101) cou
CL-------06/04/99
         read(lwin,101) (dr(i,1),i=0,ngh)
         read(lwin,101) (dr(i,2),i=0,ngh)
CL-------
  101    format(4e20.12)
C
c
         close(lwin)
         write(l6,*) ' potentials read from tape del.wel:'
         write(l6,100) nucnam1,nama1,npr1,ngh1,n0f1,nb1,nt1 
         write(l6,102) b01,b02,si
  102    format(5h b0 =,2f12.6,6h  si =,f12.6) 
c
c
c     writing of potentials to tape:
c-----------------------------------
      else
         open(lwou,file='del.wel',status='unknown')
         write(lwou,100) nucnam,nama,npr,ngh,n0f,n0b,nb,nt
         write(lwou,'(5h b0 =,2f12.6,6h  si =,f12.9,2f12.6)') 
     &              b0f,b0b,si,ea,rms
         write(lwou,'(a10,5f10.4)') txtfor,amsig,amome,amrho
         write(lwou,'(10x,7f10.4)') gsig,gome,grho,fr,g2,g3,w3
         write(lwou,104) 'Pairing:  ',ga,gg,pwi
  104    format(a,5f12.6)
         write(lwou,104) 'Delta:    ',del,dec
         write(lwou,104) 'Spk:      ',spk
         write(lwou,104) 'Lambda:   ',ala,tz
         write(lwou,101) sig
         write(lwou,101) ome
         write(lwou,101) rho
         write(lwou,101) cou
CL-------06/04/99
         write(lwou,101) (dr(i,1),i=0,ngh)
         write(lwou,101) (dr(i,2),i=0,ngh)
CL-------

         close(lwou)
         write(l6,*) ' potentials written to tape del.wel'
      endif
c
      if (lpr) then
         call prigh(0,sig,b0c,'Mesh ')
         call prigh(1,sig,one,'SIG ')
         call prigh(1,ome,one,'OME ')
         call prigh(1,rho,one,'RHO ')
         call prigh(1,cou,one,'COU ')
      endif
c
      close(lwou)
c
      write(l6,*) ' ****** END INOUT **********************************'
      return
c-end-INOUT
      end      
c======================================================================c

      subroutine iter(lpr)

c======================================================================c
c
c     main iteration for the spherical Dirac program
c
c----------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
      logical lpr,lprx
c
      common /erwar / ea,rms,qp
      common /optopt/ icm,icou,it1,it2,ncut
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      write(l6,*) ' ****** BEGIN ITER *********************************'
c
      do 10 ite = 1,maxi
         ii = ite
c
         if (lpr) then
            write(l6,'(i3,a,f12.8,2(a,f7.3),a,f5.2)') 
     &      ii,'.Iteration:  si = ',
     &      si,'  E/A = ',ea,'  R = ',rms,'  mix =',xmix  
            if (l6.ne.6) write(6,'(i3,a,f12.8,2(a,f7.3),a,f5.2)') 
     &      ite,'. Iteration:  si = ',si,'  E/A = ',ea,'  R = ',rms, 
     &      '  mix =',xmix
         endif
c------- potentials in oscillator space
         call potgh(.false.)
c
c------- single-particle matrix elements
         call singf(.false.)
c
c------- diagonalization of Dirac equation in the oscillatorbasis 
         call dirac(.false.) 
c
c------- occupation of single-particle levels
         call occup(2,.false.)
c
c------- calculation of new densities in oscillator basis
         call denssh(.false.)
c
c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- calculation of expectation values
         lprx = mod(ii,10).eq.1
         call expect(lprx)
c
c------- calculation of new fields
         call field(.false.)
c
c------- check for convergence
         ic = itestc()
         if (ic.eq.1) goto 20
         if (ic.eq.2) goto 30
c
   10 continue
   20 write(l6,101) ii,si
      if (l6.ne.6) write(6,101) ii,si
  101 format(//,1x,68(1h*),/,' *   Iteration interrupted after',i4,
     &             ' steps   si =',f17.10,'  *',/,1x,68(1h*))
      goto 40
c
   30 continue
      write(l6,100) ii,si
      if (l6.ne.6) write(6,100) ii,si
  100 format(//,1x,68(1h*),/,' *   Iteration converged after',i4,
     &             ' steps   si =',f17.10,'    *',/,1x,68(1h*))
c
   40 write(l6,*) ' ****** END ITER ***********************************'
      return
c-end-ITER
      end
c======================================================================c
c
      subroutine occup(is,lpr)
c
c======================================================================c
C
c     IS=1  occupation in fixed j-blocks
c        2  occupation by lambda-iteration (from the bottom)
c        3  occupation by particle number projected BCS     
c        4  occupation by temperature     
c
c     nt:    number of levels
c     ibl:   number of blocked level
c     tz:    particle number
c     al:    chemical potential 
c     gg:    pairing strength
c     dec:   fixed external pairing field
c     del:   gap parameter 
c     pwi:   pairing window
c
c     ee:    single particle energies
c     vv:    occupation probabilities (0 < vv < 1)
c     mu:    multiplicities of j-shells (j+1/2)
c----------------------------------------------------------------------c
      include 'paramet'
      parameter (ndwork = nwork-2*ntx)
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      common /bcsbcs/ zz,al,g,dc,sp,ex(ntx),vx(ntx),uv(ntx),mx(ntx),kt
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)
      common /fermi / ala(2),tz(2)
      common /nucnuc/ amas,nama,npr(2),jmax
      common /fixocc/ ioc(nbx,2)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ icm,icou,it1,it2,ncut
      common /pair  / ga(2),gg(2),del(2),spk(2),dec(2),pwi
      common /quaqua/ nt,nr(ntx),nl(ntx),nj(ntx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /ugugug/ itbl(2),jpnbl(3,2),nrbl(2)
      common /work  / exx(ntx),ix(ntx),mxx(ntx),dwork(ndwork)
c
c     external enbcs
c
      data maxd/100/,epsd/1.d-6/
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN OCCUP ********************************'
c
c
c---- loop over neutron-proton
      do 10 it = it1,it2
c
c
c---- fixed occupation within the differten j-blocks
      if (is.eq.1) then
         do ib = 1,nb
            im  = ia(ib)
            nd  = id(ib)
            ifx = ioc(ib,it)
            do i = 1,ifx
               vv(im+i,it) = 2*ijb(ib)      
            enddo
            do i = ifx+1,nd
               vv(im+i,it) = zero
            enddo
         enddo
      endif
c
c
c---- BCS
      if (is.gt.1) then
         al  = ala(it)
         de  = del(it)
         g   = gg(it)
         ibl = nrbl(it)
c
c------- pairing-window and blocking
         if (ibl.eq.0) then
            nz = npr(it)
         else  
            nz = npr(it) - one
         endif
         zz = nz
         kt = 0
	 nz0 = 0
         do i = 1,nt
            m = mu(i)
            if (i.eq.ibl) m = m-1
            if (ee(i,it).le.pwi.and.m.gt.0) then
               kt = kt+1
               ix(kt) = i
               ex(kt) = ee(i,it)
               mx(kt) = m
	       nz0    = nz0 + m*2 
            else
               vv(i,it) = zero
            endif
         enddo
	 if (nz0.lt.nz) stop ' in OCCUP: nz0 too small'
c
c------- simple BCS with lambda iteration
         if (is.eq.2) then
c
c           Delta-iteration
            do itd = 0,maxd
               call bcslam(kt,mx,ex,vx,uv,de,zz,al,sp,mxx,.false.)
               d0 = de
               de = g*sp + dec(it)
c               write(6,100) itd,de,d0,de-d0
c               write(l6,100) itd,de,d0,de-d0
c  100          format(i3,13h. D-Iteration,3f12.8)
               if (abs(de-d0).le.epsd) goto 20
               if (itd.eq.0) then
                  x0 = d0
                  y0 = de
               else
                  x1 = d0
                  y1 = de
                  x  = x0-x1
                  y  = y0-y1
                  de = (y1*x-x1*y)/(x-y)
                  x0 = x1
                  y0 = y1
               endif
            enddo
            stop 'in OCCUP: Delta-Iteration not converged'
   20       continue
         endif
c
c
c------- Numberprojected BCS
c        if (is.eq.3) then
c           stop ' BCSPNP not mounted '
c
c           bracketing the minimum
c           ax = de
c           bx = ax + 0.2d0
c           call mnbrak(ax,bx,cx,fa,fb,fc,enbcs)
c
c           Brent's search for minimum
c           emin = brent(ax,bx,cx,enbcs,epsd,xmin)
c           de = xmin
c           write(6,*) ' Numb.Proj. D =',xmin,emin 
c        endif
c
         ala(it) = al
         de = de + dec(it)
         del(it) = de
         spk(it) = sp
         tz(it)  = zz
         do k = 1,kt
            i = ix(k)
            vv(i,it) = 2*mu(i)*vx(k)
         enddo
         if (ibl.ne.0) vv(ibl,it) = vv(ibl,it) + one
c
         if (lpr) then
            write(l6,'(/,a)') '   k         e(k)         vv(k)'
            sn = zero
            es = zero
            do k = 1,nt
               write(l6,'(i4,i5,2f15.8)') k,mu(k),ee(k,it), 
     &                                            vv(k,it)/(2*mu(k))
               sn = sn + vv(k,it)
               es = es + vv(k,it)*ee(k,it)
            enddo
            epair = -g*sp*sp
            write(l6,'(/,a,f15.8)') ' OCCUP: delta             =',de
            write(l6,'(a,f15.8)')   '        trace of kappa    =',sp
            write(l6,'(a,f15.8)')   '        chemical potential=',al
            write(l6,'(a,f15.8)')   '        particle number   =',sn
            write(l6,'(a,f15.8)')   '        particle energy   =',es
            write(l6,'(a,f15.8)')   '        pairing energy    =',epair 
         endif
c
c
c
c     temperature   
c----------------
c     if (is.eq.4) then
c        zz = npr(it)
c        al  = ala(it)
c        stop ' in OCCUP: lamtem not initialized '
c        call lamtem(nt,zz,al,ee(1,it),mu,vv(1,it),ex,mx,lpr)
c        ala(it) = al
c        tz(it)  = zz
c     endif
      endif
   10 continue
c
      if (lpr)
     &write(l6,*) ' ****** END OCCUP **********************************'
      return
c-end-OCCUP
      end
c======================================================================c

      subroutine bcslam(n,mu,ee,vv,uv,delta,tz,al,spk,mx,lpr)

c======================================================================c
c
c    determines chemical potential for fixed Gap DELTA
c
c    N     dimension of fields ee, vv
c    TZ    wanted particle number 
c    ALAM  chemical potential
c    DELTA gap parameter 
c    SPK   trace of kappa
c    ee    single particle energies
c    vv    BCS occupation numbers (0 < vv < 1)
c    uv    BCS amplitudes u*v
c    mu    multiplicities of each shell (j+1/2)
c    mx    auxilary field
C
c----------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      dimension mu(n),mx(n),ee(n),vv(n),uv(n)
c
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      data maxl/50/,epsl/1.d-10/
c
c
      if (lpr) write(l6,*) 'BCSLAM: Delta = ',delta
c
      if (abs(delta).lt.0.0001) then
c
c------- Delta = 0: counting levels from the bottom of the well
         ntz = tz + 0.1
         do i = 1,n
            vv(i) = ee(i)
            mx(i) = mu(i)
         enddo
         call ordi(n,vv,mx)
         nt = 0
         do i = 1,n
            nt = nt + 2*mx(i)
            if (nt.ge.ntz) goto 10
         enddo             
   10    al  = vv(i)
         vtz = 2*mx(i) + ntz - nt 
         do i = 1,n
            if (ee(i).lt.al) vv(i) = one
            if (ee(i).gt.al) vv(i) = zero
            if (ee(i).eq.al) vv(i) = vtz/(2*mu(i))
            uv(i) = sqrt(vv(i)*(one-vv(i)))
         enddo
         spk = zero
      else
c
c------- Delta > 0: lambda iteration
         del2  = delta**2
         dx    = 100.d0
         dxold = dx
         do lit = 1,maxl
c
c           calculation of particle number and derivative with lambda
            dn = -tz
            dd = zero
            do i = 1,n 
               el = ee(i) - al
               e2 = one/(el**2+del2)
               e1 = sqrt(e2)
               dn = dn + mu(i)*(one-el*e1)
               dd = dd + mu(i)*e1*e2
            enddo
            dd = dd*del2
c           improvement: for lit=1 here determination of dx by bracketing
            if (dn.lt.zero) then
               xl = al
               if (lit.eq.1) xh = al + dx
            else
               xh = al
               if (lit.eq.1) xl = al - dx
            endif
            if (((al-xh)*dd-dn)*((al-xl)*dd-dn).ge.zero 
     &         .or. abs(2*dn).gt.abs(dxold*dd) ) then
c  
c              bisection because newton out of range or too slow
           write(6,*) ' Bisection'
           write(l6,*) ' Bisection'
               dxold = dx
               dx    = half*(xh-xl)
               al    = xh - dx
            else
c
c              newton 
               dxold = dx
               dx    = dn/dd
               al    = al - dx
            endif
            if (abs(dx).lt.epsl) goto 20
c   
            if (lpr.or.lit.gt.10) then
               write(l6,100) lit,'. Lambda-Iteration:',al+dx,dn,al    
               write(6,100)  lit,'. Lambda-Iteration:',al+dx,dn,al    
  100          format(i4,a,3f15.8) 
            endif
c
         enddo
         write(l6,'(a,i4,a)') 
     &        ' Lambda-Iteration interupted after',lit,' steps'
         stop
   20    if (lpr) write(l6,'(i4,a,2f15.8)')
     &      lit,'. Lambda-Iteration successful:    ',al,dn
c
c        calculation of occupation factors for fixed lambda             
         spk = zero
         do i = 1,n
            el    = ee(i) - al
            eki   = one/sqrt(el**2+del2)
            vv(i) = half*(one-el*eki)
            uv(i) = half*delta*eki
            spk   = spk + mu(i)*uv(i)
         enddo
      endif
c
      return
c-end-BCSLAM
      end
c=====================================================================c

      subroutine plotd(lpr)

c=====================================================================c
C
C     prepares plot of densities in coordinate space
C
c---------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      dimension pn(nox)
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ icm,icou,it1,it2,ncut
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2),dro(0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN PLOTD ********************************'
c
c
c-------------------------------------------------
c     test:
c     do ih=1,ngh
c        ro = rdens(1,rs(1,1),pn,xh(ih))
c        write(l6,*) 'rs neut',ih,xh(ih)*b0f,ro
c     enddo
c     do ih=1,ngh
c        ro = rdens(1,rv(1,1),pn,xh(ih))
c        write(l6,*) 'rv neut',ih,xh(ih)*b0f,ro
c     enddo
c-------------------------------------------------
c
c     number of points for the plot
      mxpl = 80
c
c     plot step in (fm)
      stpl = 0.1

c
c     plot for densities:
c------------------------
      open(lplo,file='dis.plo',status='unknown')
      do 10 it = it1,it2
c
         write(lplo,'(/,a,i3)') ' scalar density it =',it
         if (lpr) then
            write(l6,'(/,a,i3)') ' scalar density it =',it
            ro = rdens(1,rs(1,it),pn,xh(0))
            write(l6,100) b0f*xh(0),ro
         endif
         r = zero
         do ist = 0,mxpl
            x  = r/b0f
            ro = rdens(1,rs(1,it),pn,x)
            write(lplo,100) r,ro
  100       format(f10.3,f15.6) 
            if (lpr) write(l6,100) r,ro 
            r = r + stpl
         enddo
c
         write(lplo,'(/,a,i3)') ' vector density it =',it
         if (lpr) write(l6,'(/,a,i3)') ' vector density it =',it
         r = zero
         do ist = 0,mxpl
            x  = r/b0f
            ro = rdens(1,rv(1,it),pn,x)
            write(lplo,100) r,ro 
            if (lpr) write(l6,100) r,ro 
            r = r + stpl
         enddo
   10 continue
c
c     plot of specific wavefunction f(r) und g(r)

      close(lplo)
c
      if (lpr)
     &write(l6,*) ' ****** END PLOTD **********************************'
      return
C-end-PLOT
      end
c======================================================================c

      double precision function rdens(is,ro,pn,x)

c======================================================================c
c
c     calculation of the density ro at arbitrary meshpoint x
c     x is given in units of the oscillator lenght: x = r/b0f 
c
c     is = 1: density given at Gauss-Meshpoints ro(ih)
c          2: density given through oscillator expansion pn(n)
c
c----------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
c
      dimension ro(ngh),pn(nox),rnr(nox)
C
      common /dimens/ n0f,n0b,nrm,nlm
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /radbos/ rnb(1:nox,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
c
      no = n0b/2 + 1
c
c---- calculation of the oscillator expansion for ro
      if (is.eq.1) then
         do n = 1,no
            s = zero
            do ih = 1,ngh
               s = s + ro(ih)*rnb(n,ih)*wh(ih)*xh(ih)**2
            enddo
            pn(n) = s
         enddo
      endif
c
      xx = x*x
      rnr(1) = 2*wgi(0)*exp(-half*xx)
      rnr(2) = rnr(1)*(1.5d0-xx)*shi(1)
      do n = 3,no
         rnr(n)  = ((2*n-2.5d0-xx)*rnr(n-1) -
     &              sq(n-2)*sqh(n-2)*rnr(n-2))*sqi(n-1)*shi(n-1)
      enddo
c
      s = zero
      do n = 1,no
         s = s + pn(n)*rnr(n)
      enddo
      rdens = s
c
      return
c-end-RDENS
      end
c=====================================================================c

      subroutine plotw(it,ib,k,lpr)

c=====================================================================c
C
C     prepares plot of specific wafefunctions f(r) and g(r)
c     it = 1 for neutrons,  it = 2 for protons
c     ib number of block
c     k  number of specific wavefunction in this block
C
c---------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor           
c
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)      
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /initia/ vin,rin,ain,inin,iplot
      common /optopt/ icm,icou,it1,it2,ncut
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
cllllll 06/99
      common /match/ mxpl,stpl
cllllll 06/99
claudio 04/2000
      ibg = ib - 1 + 2*mod(ib,2)
      nf  = id(ib)
      ng  = id(ibg)
claudio 04/2000
c
c
c      if (lpr)
c     &write(l6,*) ' ****** BEGIN PLOTW ********************************'
       k1=k-ia(ib)
c
c     number of points for the plot
clll      mxpl = 180
c      mxpl=1
c
c     plot step in (fm)
clll      stpl = 0.1
c
c     plot for wavefunctions:
c----------------------------
      open(lplo,file='dis.wplo',status='unknown')
      ip = 2-mod(ib,2)
c      write(lplo,)'El negocio es aqui'
      write(lplo,111)  it,k1,2*ijb(ib)-1,tp(ip),ee(k,it)
111   format(/,3i3,'/2',a1,f10.3)     
c         write(lplo,'(/,a,3i3)') ' wavefunction f(r)',it,ib,k
         if (lpr)  then
          if (iplot.eq.1 .or. iplot.eq.2)
     +     write(l6,111) it,k1,2*ijb(ib)-1,tp(ip),ee(k,it)
claudio 04/2000
          if (iplot.eq.1) write(l6,'(1x,2i3)') nf,ng
claudio 04/2000
         endif
c         r = zero
c         s = zero
c         do ist = 0,mxpl
c            call rwave(it,ib,k1,r,f,g)
c           write(lplo,100) r,f,g
c            write(lplo,100) r,f
  100       format(f10.3,2f15.6) 
c           if (lpr) write(l6,100) r,f,g 
c            if (lpr) write(l6,100) r,f 
c            s = s + (f*f+g*g)*r*r
c            r = r + stpl
c         enddo
c         write(6,*) ' check norm of f and g',s*stpl
c         write(l6,*) ' check norm of f and g',s*stpl
c
c         write(lplo,'(/,a,3i3)') ' wavefunction g(r)',it,ib,k
c         if (lpr) then
c            write(l6,'(/,a,3i3)') ' wavefunction f(r) and g(r)',it,ib,k
c         endif
         r = zero
         s = zero
         s1 = zero
         s2 = zero
         do ist = 0,mxpl
            !write(l6,*)'ib',ib
           call rwave(it,ib,k1,r,f,g,ist)
           write(lplo,100) r,f,g
c            write(lplo,100) r,g
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Coloquei o comentario abaixo para evitar a impressao das funcoes 
C     de onda no arquivo dis.out 06/98
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (lpr .and. iplot.eq.2 ) write(l6,100) r,f,g 
c            if (lpr) write(l6,100) r,g 
            s = s + (f*f+g*g)*r*r
            s1 = s1 + f*f*r*r
            s2 = s2 + g*g*r*r
            r = r + stpl
         enddo
c         write(6,*) ' check norm of f and g',s*stpl
c         write(l6,*) ' check norm of f and g',s*stpl
        ! if (lpr .and. iplot.eq.1 ) 
    ! +    write(l6,*) ' integral de (f*f+g*g)*r*r --> ',s*stpl
    !     if (lpr .and. iplot.eq.1 ) 
    ! +    write(l6,*) ' integral de  f*f*r*r      --> ',s1*stpl
    !     if (lpr .and. iplot.eq.1 ) 
    ! +    write(l6,*) ' integral de  g*g*r*r      --> ',s2*stpl                           
      close(lplo)
c
c      if (lpr)
c     &write(l6,*) ' ****** END PLOTW **********************************'
      return
C-end-PLOT
      end
c======================================================================c

      subroutine rwave(it,ib,k,r,f,g,imc)

c======================================================================c
c
c     calculation of the wavefunctions f(r) and g(r) at point x
c     x is given in units of the oscillator lenght: x = r/b0f 
c     it = 1:  neutron,  it = 2: proton
c     ib is the block charakterized by j,l
c     k  is the number of the wavefunction within this block
c
c----------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
c
      dimension rnl(nrx,2)
C
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /dimens/ n0f,n0b,nrm,nlm
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /initia/ vin,rin,ain,inin,iplot
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /wavefg/ fg(nq2x,nb2x)
c
      if (r.eq.zero) r = 0.0000001
      ibg = ib - 1 + 2*mod(ib,2)
      nf  = id(ib)
      ng  = id(ibg)
      nd  = nf + ng
      mf  = ib + (it-1)*nbx
      lf  = ilb(ib)
      lg  = ilb(ibg)
      
C---------------------      
      k1  = k
C---------------------      
c     
      write(l6,*) 'NG',ng,'NF' ,nf,'k',k,'lf',lf,'lg',lg
      x  = r/b0f
      xx = x*x
      f  = b0f**(-1.5d0)
      !write(l6,*) x, xx, f
      rnl(1,1) = sq(2)*f*wgi(lf+1)*x**lf*exp(-half*xx)
      rnl(2,1) = rnl(1,1)*(lf+1.5d0-xx)*shi(lf+1)
      rnl(1,2) = sq(2)*f*wgi(lg+1)*x**lg*exp(-half*xx)
      rnl(2,2) = rnl(1,2)*(lg+1.5d0-xx)*shi(lg+1)
      do n = 3,nrm
         rnl(n,1)  = ((2*n+lf-2.5d0-xx)*rnl(n-1,1) -
     &           sq(n-2)*sqh(n-2+lf)*rnl(n-2,1))*sqi(n-1)*shi(n-1+lf)
         rnl(n,2)  = ((2*n+lg-2.5d0-xx)*rnl(n-1,2) -
     &           sq(n-2)*sqh(n-2+lg)*rnl(n-2,2))*sqi(n-1)*shi(n-1+lg)
      enddo
c
      sf = zero
      sg = zero
C------------------      
      snf = zero
      sng = zero
C------------------
      if (imc.eq.0 .and. iplot.eq.1) then 
      write(l6,*) ' Coeficientes de f(r)'           
      end if
      do n = 1,nf
           ! write(l6,*) 'MF' ,mf
         sf = sf + fg(n+(k-1)*nd,mf)*rnl(n,1)
         write(l6,*)'fgcan1',fg(n+(k-1)*nd,mf)
C------------------         
         snf = snf + fg(n+(k-1)*nd,mf)*fg(n+(k-1)*nd,mf)
         if (imc.eq.0 .and. iplot.eq.1) then 
         write(l6,'(f15.6)') fg(n+(k-1)*nd,mf)
         end if
C------------------                
      enddo
      if (imc.eq.0 .and. iplot.eq.1) then       
      write(l6,*) ' Coeficientes de g(r)'      
      end if
      do n = 1,ng
         sg = sg + fg(nf+n+(k-1)*nd,mf)*rnl(n,2)
         
C------------------         
cll         sng = sng + fg(nf+n+(k-1)*nd,mf)*rnl(n,2)
         sng = sng + fg(nf+n+(k-1)*nd,mf)*fg(nf+n+(k-1)*nd,mf)
         if (imc.eq.0 .and. iplot.eq.1) then 
         write(l6,'(f15.6)') fg(nf+n+(k-1)*nd,mf)
         end if
C------------------                        
      enddo
      !write(l6,*) 'sf,sg',sf,sg
      f = sf
      g = sg
C------------------
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Coloquei o comentario abaixo para evitar a impressao das normas 
C     \int f^2 + g^2 na tela 06/98
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      write(*,*) 'isospin,n_r,j,norm',it,k1,2*ijb(ib)-1,snf+sng
         if (imc.eq.0 .and. iplot.eq.1) then 
      write(l6,*) 'Sum f,Sum g',snf,sng
         end if
C------------------     
c
      return
c-end-RWAVE
      end
c======================================================================c

      subroutine potgh(lpr)

c======================================================================c
C
C     CALCULATION OF THE POTENTIALS AT GAUSS-MESHPOINTS
C
c----------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /physco/ amu,hqc,alphi,r0
      common /potpot/ vps(0:ngh,1:2),vms(0:ngh,1:2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN POTGH ********************************'
c
      do ih = 0,ngh
         s  = gsig*sig(ih)
         go = gome*ome(ih)
         gr = grho*rho(ih) 
         v1 = go - gr
         v2 = go + gr + cou(ih) 
         vps(ih,1) = v1 + s
         vms(ih,1) = v1 - s
         vps(ih,2) = v2 + s
         vms(ih,2) = v2 - s
      enddo
C
      if (lpr) then
         call prigh(0,sig,b0c,'X(FM) ')
         call prigh(1,vps(0,1),hqc,'V+S  n')
         call prigh(1,vms(0,1),hqc,'V-S  n')
         call prigh(1,vps(0,2),hqc,'V+S  p')
         call prigh(1,vms(0,2),hqc,'V-S  p')
      endif

c
      if (lpr)
     &write(l6,*) ' ****** END POTGH **********************************'
      return
c-end-POTGH
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
cb0 = sqrt(two*hb0/hom)
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

      subroutine resu

c======================================================================c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
c
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2),dro(0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
c
      write(l6,*) ' ****** BEGIN RESU *********************************'
c
c
c---- printing of densities
      call prigh(0,rs,b0f,'x(fm) ')
      do it = 1,2
         call prigh(1,rs(0,it),one,'ROS '//tis(it))
         call prigh(1,rv(0,it),one,'ROV '//tis(it))
      enddo
      call prigh(1,dro,one,'DRO  ')
c
c---- printing of fields
      call prigh(0,sig,b0c,'x(fm) ')
      call prigh(1,sig,one,'Sigma ')
      call prigh(1,ome,one,'Omega ')
      call prigh(1,rho,one,'Rho   ')
      call prigh(1,cou,one,'Coulom')
c
c
c
c---- single particle energies
      write(l6,100) 
  100 format(//,' Single-particle Energies',/,1x,24(1h-),/,
     $20x,'neutrons',35x,'protons')
      en0 = ee(1,1)
      ep0 = ee(1,2)
      do ib = 1,nb
         nf  = id(ib)
         im  = ia(ib)+1
         j2  = 2*ijb(ib)-1
         ip  = 2-mod(ib,2)
c
         ie = im+nf-1
         do 20 i = im,ie 
            v1 = vv(i,1)/(2*mu(i))
            v2 = vv(i,2)/(2*mu(i))
            if (v1.lt.0.0001.and.v2.lt.0.0001) goto 20
            if (i.eq.im) write(l6,*)
            write(l6,101) j2,tp(ip),ee(i,1),ee(i,1)-en0,v1,
     &                               ee(i,2),ee(i,2)-ep0,v2
 101        format(' j =',i2,'/2',a1,2f10.3,f8.3,8x,2f10.3,f8.3) 
  20     continue
      enddo
c
c
      call expect(.true.)
c
      write(l6,*) ' ****** END RESU ***********************************'
      return
c-end-RESU
      end 
c======================================================================c

      subroutine singf(lpr)

c======================================================================c
c
c     calculates single particle matrix elements for Fermions       
c     in the spherical oscillator basis
c
c----------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /single/ sp(nq2x,nb2x)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN SINGF ********************************'
c
      do  5 it = 1,2
      do 10 ib = 1,nb,2
         np = id(ib)
         nm = id(ib+1)
         mp = ib + (it-1)*nbx         
         mm = ib + 1 + (it-1)*nbx
c
c        SIGMA*P
c----------------
         call sigp(nm,np,it,ib,sp(1,mp),lpr)
         if (ib+1.gt.nbx) stop ' in SINGF: nbx too small'
         do k = 1,nm
         do i = 1,np
cl            sp(i+(k-1)*np,ib+1) = sp(k+(i-1)*nm,ib)
            sp(i+(k-1)*np,mm) = sp(k+(i-1)*nm,mp)
         enddo
         enddo
c
   10 continue
    5 continue
C
      if (lpr)
     &write(l6,*) ' ****** END SINGF **********************************'
      return
c-end-SINGF
      end  
c=====================================================================c

      subroutine sigp(nm,np,it,ib,tt,lpr)

c=====================================================================c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
      logical lpr
      character*8 tb
      character*25 txb
c
      dimension tt(nm,np),yp(0:ngh),rb(0:ngh)
CL      dimension drb(0:ngh),pn(nox)
c
      common /baspar/ hom,hb0,b0f,b0b,b0c,cf,cb
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ amu,hqc,alphi,r0
      common /rados1/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /texblo/ tb(ntx),txb(nbx)      
C      
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3      
      common /fields/ sig(0:ngh),ome(0:ngh),rho(0:ngh),cou(0:ngh)
      common /extrct/ gprime
      common /rhoro2/ rs2(0:ngh,1:2),rv2(0:ngh,1:2),rr2(0:ngh,1:2),
     &                dr(0:ngh,1:2),drold(0:ngh,1:2)
      common /gfviv / iv(0:igfv)
c
      j  = ijb(ib)
      lp = j - mod(j,2)
      lm = j - mod(j+1,2)
      if (lp.eq.j) then
         kk =-j
      else
         kk = j
      endif
c
      kk = kk - 1    
      f  = one/b0f
C
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
C !!!!!!! The field rho has dimension fm^-1 and it is in the common    !! 
C !!!!!!! called field at the radius rb(ih) = xh(ih)*b0c               !!
C !!!!!!! xh(ih) is dimensionless                                      !!
C !!!!!!! gt is in fm^2  and vt is dimensionless                       !!
C-----------------------------------------------------------------------C
C !!!!!!! The densit drb has dimension fm**(-3)                        !!
C !!!!!!! dgt is in fm^3  and dvt is dimensionless                     !!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
CL      gt = fr*hqc**2/2/amu/amrho
      gt  = fr*hqc/2/amu/f
CL-------- 25/03/99
      dgt = (fr*hqc/2/amu)**2/f
ccccc      dgt = zero
CL--------
      do ih = 0,ngh
         rb(ih)  = xh(ih)*b0c
CL         drb(ih) = rdens(1,dr(0,it),pn,xh(ih))
      enddo   
      do n2 = 1,np
      do n1 = 1,nm
         s   = zero
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!! The subroutines spline and splint calculate rho and theirs   !! 
C !!!!!!! first derivative at rb(ih)                                   !!       
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call spline(rb,rho,yp,ngh,zero,zero)
         do ih =1,ngh        
         call splint(1,rb,rho,yp,ngh,rb(ih),rh,drh,y2)
CL EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
CL EEEEEEE Foi encontrado um erro nessa modificacao 15/11/98                         
CL           vt = iv(it)*gt*( drh - (kk+1)*rh/rb(ih) )
CL EEEEEEE A modificacao correta esta na linha abaixo
CL EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
           vt  = iv(it)*gt*drh
CL--------- 25/03/99
CL           dvt = third*dgt*drb(ih)
           dvt = gprime*dgt*dr(ih,it)
           s = s + rnl(n1,lm,ih) * 
     &              ( - rnl1(n2,lp,ih) + kk*rnl(n2,lp,ih)/xh(ih)
     &                + (-vt + dvt)*rnl(n2,lp,ih)  )
CL---------
         enddo
         tt(n1,n2) = f*s
      enddo
      enddo
c
      if (lpr) then
         write(l6,'(a)') txb(2*j-1)
         iap = ia(2*j-1) + 1
         iam = ia(2*j) + 1
         call aprint(1,3,1,nm,nm,np,tt,tb(iam),tb(iap),'Sigma * P')
      endif
c
      return
c-end-SIGP
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
c=====================================================================c

      subroutine plw

c=====================================================================c
C
C     plotes the  wafefunctions f(r) and g(r)
C
c---------------------------------------------------------------------c
      include 'paramet'
c
      implicit real*8 (a-h,o-z)
cllllll 06/99
      character*1 tp,tl,tis
      character*2 nucnam
      character*10 txtfor
cllllll 06/99
c
c
      common /bloblo/ nb,ijb(nbx),ilb(nbx),
     &                id(nbx),idq(nbx),ia(nbx),iaq(nbx)
      common /eeeeee/ ee(ntx,2),vv(ntx,2),mu(ntx)     
cllllll 06/99      
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /physco/ amu,hqc,alphi,r0
      common /nucnuc/ amas,nama,nneu,npro,jmax
      common /textex/ nucnam,tp(2),tis(2),tl(0:20),txtfor
      common /mespar/ amsig,amome,amrho,gsig,gome,grho,fr,g2,g3,w3
      common /extrct/ gprime            
      common /mathco/ zero,one,two,half,third,pi
      common /initia/ vin,rin,ain,inin,iplot
      common /match/ mxpl,stpl
cllllll 06/99
cllllll 06/99
c     number of points for the plot
cll      mxpl = 180
      mxpl = 180
c
c     plot step in (fm)
      stpl = 0.1
        ii = 0
cllllll 06/99
      do it=1,2
         do ib=1,nb
            nf  = id(ib)
            im  = ia(ib)+1
            ie = im+nf-1
            do 10 k=im,ie
c            v1 = vv(i,1)/(2*mu(i))
c            v2 = vv(i,2)/(2*mu(i))
c            if (v1.lt.0.0001.and.v2.lt.0.0001) goto 10
             if (ee(k,it).gt.13.)               goto 10
             ii = ii + 1
10          continue
         enddo
      enddo
      if (iplot.eq.1 .or. iplot.eq.2) then
      write(l6,'(1X,i4,2F6.2,i4,F6.2)') ii,zero,stpl,mxpl+1,mxpl/10.
      write(l6,*) ' ****** END PLW ************************************'
      write(l6,'(1X,a2,i3)') nucnam,nama
      write(l6,'(1X,a,a10)') 'force ',txtfor      
      write(l6,'(1X,4F8.3)') amsig,amome,amrho,amu      
      write(l6,'(1X,3F10.3)') gsig,gome,grho
      write(l6,'(1X,3F10.3)') g2,g3,w3
      write(l6,'(1X,F7.3)') fr
      write(l6,'(1X,F7.3)') gprime      
c      
      end if

      write(l6,*) ' ****** END PLW ************************************'
      do it=1,2
         do ib=1,nb
            nf  = id(ib)
            im  = ia(ib)+1
            ie = im+nf-1
            do 20 k=im,ie

                  !write(l6,*)'nf', nf, 'ib', ib, 'ie', ie, 'k', k
c            v1 = vv(i,1)/(2*mu(i))
c            v2 = vv(i,2)/(2*mu(i))
c            if (v1.lt.0.0001.and.v2.lt.0.0001) goto 20
                  if (ee(k,it).gt.13.)               goto 20
                  call plotw(it,ib,k,.true.)
20          continue
         enddo
      enddo
      write(l6,*) ' ****** END PLW ************************************'
      return
C-end-PLW
      end
c======================================================================c


