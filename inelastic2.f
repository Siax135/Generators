      program inelastic2


C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...PYTHIA Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYDATR/MRPY(6),RRPY(100)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      COMMON/PYMSRV/RVLAM(3,3,3), RVLAMP(3,3,3), RVLAMB(3,3,3)
      COMMON/PYTCSM/ITCM(0:99),RTCM(0:99)
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      COMMON/PYLH3P/MODSEL(200),PARMIN(100),PAREXT(200),RMSOFT(0:100),
     &     AU(3,3),AD(3,3),AE(3,3)
      COMMON/PYLH3C/CPRO(2),CVER(2)
      CHARACTER CPRO*12,CVER*12

c...Required varibales for running
      integer I,MJ,NUMEV,MELECLOC,NGENEV
      integer NPRINT,NSEED
      integer MCHARGE(-2212:2212)
      double precision Q2,W
      logical NEU,PIONPL,PIONMI,PRO,ELEC,VERBOSE,HREAC,DREAC
      character ARG*32,OUTPUT*32,PAIRSTR*10,TARGET*2,PAIROP*1
      real ELECANG
      double precision P(4000,5),V(4000,5)
      integer N,NPAD,K(4000,5)
      COMMON/PYJETS/N,NPAD,K,P,V
      COMMON/PIDCHRG/MCHARGE

      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYDATR/,/PYSUBS/,
     &/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,/PYINT4/,/PYINT5/,
     &/PYINT6/,/PYINT7/,/PYMSSM/,/PYSSMT/,/PYMSRV/,/PYTCSM/,
     &/PYBINS/,/PYLH3P/,/PYLH3C/,/PIDCHRG/

c...Set start values for variables
      NUMEV = 0
      NEU = .FALSE.
      PIONPL = .FALSE.
      PIONMI = .FALSE.
      PRO = .FALSE.
      ELEC = .FALSE.

c...Set defaults for options
      OUTPUT = 'out.dat'
      TARGET = 'p+'
      NGENEV = 20
      NPRINT = 5
      PAIROP = '0'
      NSEED = 19780503  ! default set by PYTHIA
      VERBOSE = .FALSE.
      HREAC = .FALSE.
      DREAC  = .FALSE.
      PAIRSTR = 'MSTJ(12)=0'

c...Set up format for the LUND format
100   format(11X,I2,     ! Number of particles
     1 2X,I1,            ! Num of target nucleons
     2 2X,I1,            ! Num of target protons
     3 2X,F5.3,          ! Target polarization
     4 2X,F5.3,          ! Beam polarization
     5 2X,F5.3,          ! x
     6 1X,F5.3,          ! y
     7 2X,F8.3,          ! W
     8 1X,F8.3,          ! Q^2
     9 2X,F5.3)          ! nu

200   format(4X,I2,      ! Index
     1 3X,I2,            ! Charge
     2 8X,I1,            ! Type
     3 8X,I7,            ! PID
     4 7X,I2,            ! Parent Index
     5 7X,I2,            ! Daughter Index
     6 8X,F7.4,          ! Px
     7 3X,F7.4,          ! Py
     8 3X,F7.4,          ! Pz
     9 3X,F7.4,          ! E
     1 3X,F7.4,          ! mass
     2 3X,F8.4,          ! Vertex x
     3 3X,F8.4,          ! Vertex y
     4 3X,F8.4)          ! Vertex z

300   format('Run Parameters',/,
     1 'Output file:',13X,A32,/,
     2 'Num events:',6X,I10,/,
     3 'Num events/print:',I10,/,
     4 'Target:',18X,A2,/,
     5 'Pair prodution:',9X,A1,/,
     6 'RNG seed:',15X,I9,/)

c...Parse given arguments
      J = 1
      do
        call get_command_argument(J,ARG)
        if(LEN_TRIM(ARG) .EQ. 0) exit

        if(TRIM(ARG) .EQ. '-o') then  ! parse output file
          J = J+1
          call get_command_argument(J,OUTPUT)
        else if(TRIM(ARG) .EQ. '-n') then  ! parse num of events to generate
          J = J+1
          call get_command_argument(J,ARG)
          read(ARG,'(I10)') NGENEV
        else if(TRIM(ARG) .EQ. '-n_print') then  ! parse num of events between prints
          J = J+1
          call get_command_argument(J,ARG)
          read(ARG,'(I10)') NPRINT
        else if(TRIM(ARG) .EQ. '-seed') then  ! parse RNG seed
          J = J+1
          call get_command_argument(J,ARG)
          read(ARG,'(I9)') NSEED
        else if(TRIM(ARG) .EQ. '-pair_pro') then  ! parse pair production settings
          J = J+1
          call get_command_argument(J,ARG)
          PAIROP = ARG
          PAIRSTR = 'MSTJ(12)='//ARG
        else if(TRIM(ARG) .EQ. '-hreac') then  ! parse e,pi-,neutron option
          HREAC = .TRUE.
          print*,"Looking for events with an electron, pi-, and neutron 
     +in the final state."
        else if(TRIM(ARG) .EQ. '-dreac') then  ! parse e,pi-,pi+,proton,neutron option
          DREAC = .TRUE.
          print*,"Looking for events with an electron, proton, pi+, pi-,
     + and neutron in the final state."
        else if(TRIM(ARG) .EQ. '-v') then  ! parse verbose option
          VERBOSE = .TRUE.
        else if(TRIM(ARG) .EQ. '-neutron') then ! parse option to change target to neutron
          TARGET = 'n0'
        else if(TRIM(ARG) .EQ. '-h') then  ! parse help message
          write(*,*) 'Options:'
          write(*,*) '-o           Output file name (default: out.dat)'
          write(*,*) '-n           Number of events to generate (default
     +: 20)'
          write(*,*) '-n_print     Number of events between print statem
     +ents (default: 5)'
          write(*,*) '-neutron     Change target from proton to neutron,
     + this flag doesn''t take a following argument'
          write(*,*) '-hreac       Look for events containing an electro
     +n, pi-, and neutron in the final state, this flag doesn''t take a 
     +following argument'
          write(*,*) '-dreac       Look for events containing an electro
     +n, proton, pi+, pi-, and neutron in the final state, this flag doe
     +sn''t take a following argument'
          write(*,*) '-seed        Seed for RNG, allowed values 0 <= see
     +d <= 900000000 (default: 19780503)'
          write(*,*) '-v           Add this flag to set the output to be
     + more verbose, this flag doesn''t take a following argument'
          write(*,*) '-pair_pro    Set how diquark-antidiquark pair prod
     +uction is handled. Options are as follows:'
          write(*,*) '              0 - no diquark-antidiquark pair prod
     +ution at all (default)'
          write(*,*) '              1 - diquark-antidiquark pair product
     +ion allowed, diquark treated as unit'
          write(*,*) '              2 - diquark-antidiquark pair product
     +ion allowed, with possibility for diquark to be split'
          goto 40
          else
            write(*,*) 'Flag ',ARG,' not recognized'
        endif
        J = J+1
      enddo

c...Pass parameters to pythia
      call PYGIVE('PARP(2)=2D0')
      call PYGIVE(PAIRSTR)

c...Print run parameters
      write(*,300) OUTPUT,NGENEV,NPRINT,TARGET,PAIROP,NSEED

c...Open output file
      open(2,file=OUTPUT)

c...Set seed
      MRPY(1) = NSEED

      write(*,*) 'Starting Initialization'

c...Initialize pythia
      call pyinit('FIXT','gamma/e-',TARGET,11D0)

      write(*,*) 'Initialized'

c...Start event loop
      do I=0,NGENEV

        ELEC = .FALSE.
        NEU = .FALSE.
        PRO = .FALSE.
        PIONPL = .FALSE.
        PIONMI = .FALSE.

c...Print event number
        if( MOD(I,NPRINT) .EQ. 0) then
          write(*,*) 'Event: ',I
        endif

c...Generate event
        call pyevnt

        if(VERBOSE) then
          write(*,*)
          write(*,*) 'Event: ',I
          call pylist(2)
        endif

c...Look for electron
        do MJ=3,N
          if (K(MJ,2) .EQ. 11) then
            ELEC = .TRUE.
            MELECLOC = MJ
            ELECANG = ACOS((P(MELECLOC,3)/P(MELECLOC,4)))*(180/3.14159)
          endif
        enddo

c...Calculate Q2 and W
        Q2 = 4*11.00051*P(MELECLOC,4)*(SIN((ELECANG*(3.14159/180))/2)
     + **2)
        W = SQRT((P(2,5)**2)+(2*(11.00051-P(MELECLOC,4))*P(2,5))-Q2)

c...Trim event list to only include particle in the final state
        call pyedit(1)

c...Look for events that have at least an electron, pi+, and neutron in the final state
        if (HREAC) then
          do MJ=1,N
            if (K(MJ,2) .EQ. 2112) then ! look for neutron
              NEU = .TRUE.
            else if (K(MJ,2) .EQ. 211) then ! look for pi+
              PIONPL = .TRUE.
            endif
          enddo

          if (ELEC .AND. NEU .AND. PIONPL) then
            goto 50
          else
            goto 10
          endif
        endif

c...Look for events that have at least an electron, pi+, pi-, proton, and neutron in the final state
        if(DREAC) then
          do MJ=1,N
            if (K(MJ,2) .EQ. 2112) then ! look for neutron
              NEU = .TRUE.
            else if (K(MJ,2) .EQ. 211) then ! look for pi+
              PIONPL = .TRUE.
            else if (K(MJ,2) .EQ. -211) then ! look for pi-
              PIONMI = .TRUE.
            else if (K(MJ,2) .EQ. 2212) then ! look for proton
              PRO = .TRUE.
            endif
          enddo

          if(ELEC .AND. PRO .AND. PIONPL .AND. PIONMI .AND. NEU) then
            goto 50
          else
            goto 10
          endif
        endif

c...Write event out in LUND format
50      write(2,100) N,1,1,0.000,0.000,0.000,0.000,W,Q2,0.000
          do MJ=1,N
            write(2,200) MJ,MCHARGE(K(MJ,2)),1,K(MJ,2),K(MJ,3),0,
     +          P(MJ,1),P(MJ,2),P(MJ,3),P(MJ,4),
     +          P(MJ,5),V(MJ,1)/10,V(MJ,2)/10,
     +          V(MJ,3)/10
          enddo

          NUMEV = NUMEV+1
10    enddo
      
      if(VERBOSE) then
        call pystat(1)
      endif

      print*,"Total events in ",TRIM(OUTPUT),": ",NUMEV

40    stop
      end