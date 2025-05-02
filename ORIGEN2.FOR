      PROGRAM ORIGEN2
C
C-----------------------------------------------------------------------
C
C  OOOOOOO    RRRRRR     IIIIII     GGGGGG     EEEEEE    NN    NN    2222
C  OO   OO    RR   RR      II      GG          EE        NNN   NN   22  22
C  OO   OO    RRRRRR       II      GG    GGG   EEEEE     NN NN NN      22
C  OO   OO    RR   RR      II      GG     GG   EE        NN   NNN    22
C  OOOOOOO    RR   RR    IIIIII     GGGGGGG    EEEEEE    NN    NN   222222
C
C-------------------------------------------------------------------------------
C     ORIGEN2 V2.1    8/1/91    by:    Scott B. Ludwig
C                                      Chemical Technology Division
C                                      Oak Ridge National Laboratory
C                                      (615) 574-7916   FTS 624-7916
C
C     This is a preliminary update of the ORIGEN2 code package.  The primary 
C     purpose of this update is to include the revised set of PWR and BWR
C     cross section libraries for standard- and extended-burnups, as documented
C     in ORNL/TM-11018 (December, 1989).  These libraries were used in updating
C     the Characteristics Data Base (see DOE/RW-0184, Rev. 1).
C
C     This update to ORIGEN2 also provides a single version for both mainframe
C     and PC-based hardware.  The PC version is limited to use on hardware 
C     with 32-bit processors (80386 or 80486 with coprocessors).  The source
C     code used to create the PC version executable may be uploaded to the
C     mainframe computer and (should) compile and execute correctly.
C
C     This code compiles and executes correctly on VAX hardware (VMS operating
C     system).  Any problems with using this code on other machines should be
C     reported to S. B. Ludwig.  IBM (MVS) users may with to replace the OPEN
C     statements below with appropriate JCL.
C
C-------------------------------------------------------------------------------
C
C     ORIGEN2 - OAK RIDGE ISOTOPE GENERATION AND DEPLETION CODE
C
C     First converted to the IBM PC in 1986 by J. Marable
C
C     First full conversion (mainframe -> PC)
C                                     by Rob Tayloe, Battelle 21-Apr-1991
C
C-----------------------------------------------------------------------
C
C     ORIGEN2 DOCUMENTATION MAY BE OBTAINED FROM
C
C     RADIATION SHIELDING INFORMATION CENTER
C     OAK RIDGE NATIONAL LABORATORY
C     P.O. BOX 2008
C     BLDG. 6025, MS-6362
C     OAK RIDGE, TENNESSEE  37381
C     (615) 574-6176     FTS 624-6176
C
C     DOCUMENTATION
C     *************
C
C     ORNL-4628     "ORIGEN:  THE ORNL ISOTOPE GENERATION AND DEPLETION
C                    CODE" BY M. J. BELL, OAK RIDGE NATIONAL LABORATORY
C
C     ORNL-5621     "ORIGEN2 - A REVISED AND UPDATED VERSION OF THE
C                    OAK RIDGE ISOTOPE GENERATION AND DEPLETION CODE"
C                    BY A. G. CROFF, OAK RIDGE NATIONAL LABORATORY
C
C     ORNL/TM-7175  "A USER'S MANUAL FOR THE ORIGEN2 COMPUTER CODE"
C                    BY A. G. CROFF, OAK RIDGE NATIONAL LABORATORY
C
C-----------------------------------------------------------------------
C
C     This is the new BIG version, using bigger arrays in PARAMS.O2
C
      INCLUDE 'PARAMS.O2'
      LOGICAL   LONG
      INTEGER*2 LOCA,NON0,KD,LOC,NGF,NGN,NGR,NYIELD,NONP,NQ,MMAX,KAP,
     $LOCP,NFUDFP
      DOUBLE PRECISION CIMN,CSUM
      DIMENSION XNEW(AAAA,BBBB),COEFF(CCCC,BBBB),NPROD(CCCC,BBBB),
     $MMAX(BBBB),KAP(BBBB)
C  THE FIRST DIMENSION OF STTFPB,ISTOTI,IS, AND RSTOTI IS THE # OF
C  STORAGE VECTORS, WHICH IS LX            （这段注释说明了STTFPB、ISTOTI、IS 和 RSTOTI 这些数组的第一个维度是储存向量的数量，即 LX。）
      DIMENSION STTFPB(JJJJ,10),ISTOTI(JJJJ,03),IS(JJJJ),RSTOTI(JJJJ)
C  ARRAYS DR, ER, AND FR ARE USED TO INITIALIZE RMULV, THE VARIABLE
C  MULTIPLIER ARRAY.  THEIR DIMENSION SHOULD BE THE SAME AS THE FIRST
C  DIMENSION OF RMULV, WHICH IS LAM.（使用 DR、ER 和 FR 数组来初始化 RMULV 数组。具体来说，这些数组的维度应与 RMULV 数组的第一个维度相同，即 LAM。）
      DIMENSION DR(NNNN),ER(NNNN),FR(NNNN)
      DIMENSION YIELD(EEEE),NYIELD(FFFF),RMULV(NNNN,3)
      DIMENSION ALPHN(GGGG),NUCAN(GGGG),NUCSFU(GGGG),NY(GGGG),YY(GGGG),
     $FFSF(GGGG),FFA(GGGG)
      COMMON /JUNK/ERR,IDM(1),ILITE,IACT,IFP,ITOT,ILMAX,IAMAX,IFMAX,
     $ITMAX,IZMAX,AXN,QXN,FLUX,POWER,INDEX,TFPBAV(4),IPHMAX
      COMMON /MAIN03/NSTP,ANMUL,ANEXP,NABMAX,      ICNMAX,IAPMAX,IFYMAX
C 1766 WORDS ARE NECESSARY IN /NUDSCR/ BEGINNING WITH S
C /NUDSCR/ IS USED FOR MULTIPLE PURPOSES.（在名为/NUDSCR/的公共块中，需要1766个单词，并且这些单词的名称以S开头。）
      COMMON /NUDSCR/DUM1(CCCC,BBBB),DUM2(HHHH,BBBB),S(2),CIMN(BBBB),
     $ CSUM(BBBB),NONP(BBBB),NQ(BBBB),XP(BBBB),XPAR(BBBB),XTEMP(BBBB),
     $ D(BBBB),AP(IIII),LOCP(IIII),LONG(BBBB)
C  /BIG/ CONTAINS ARRAYS THAT ARE TO BE KEPT INTACT FOR THE ENTIRE RUN  （/BIG/ 包含整个运行过程中需要保持完整的数组）
      COMMON /BIG/NUCL(BBBB),Q(BBBB),FG(0004),TOCAP(BBBB),GENNEU(GGGG),
     $ALPHAN(GGGG),SPONF(GGGG),SFNU(GGGG),FISS(GGGG),NUCAB(BBBB),
     $AMPC(BBBB),WMPC(BBBB),XSTORE(JJJJ,BBBB),DIS(BBBB),B(BBBB),
     $ABUND(KKKK),NON0(BBBB),KD(BBBB),LOC(DDDD),NGF(BBBB),NGN(BBBB),
     $NGR(LLLL),GGR(LLLL),A(DDDD),LOCA(DDDD),NFUDFP(FFFF,OOOO)
C DR,ER, AND FR PROVIDE A CONVENIENT MECHANISM FOR INITIALIZING VARIABLE
C MULTIPLIER ARRAY RMULV.（DR, ER 和 FR 提供了一种便捷的机制来初始化变量倍乘数组 RMULV。）
      EQUIVALENCE (DR(1),RMULV(1,1)),(ER(1),RMULV(1,2)),
     $(FR(1),RMULV(1,3))
      EQUIVALENCE (DUM1(1,1),COEFF(1,1)),(DUM2(1,1),NPROD(1,1)),
     $(NONP(1),MMAX(1)),(KAP(1),NQ(1)),(XNEW(1,1),DUM1(1,1))
      EQUIVALENCE (XP(1),ALPHN(1)),(ALPHN(GGGG),NUCAN(1)),(NUCAN(GGGG),
     $NUCSFU(1)),(NUCSFU(GGGG),NY(1)),(NY(GGGG),YY(1)),(YY(GGGG),
     $FFSF(1)),(FFSF(GGGG),YIELD(1)),(YIELD(EEEE),NYIELD(1))
      CHARACTER*8 CODE
C
C     THE FOLLOWING FILES ARE USED FOR INPUT DATA（下面的文件用于输入数据）
C
      OPEN(UNIT=5  ,FILE='TAPE5.INP'  ,STATUS='unknown')
      OPEN(UNIT=3  ,FILE='TAPE3.INP'  ,STATUS='unknown')
      OPEN(UNIT=4  ,FILE='TAPE4.INP'  ,STATUS='unknown')
      OPEN(UNIT=9  ,FILE='TAPE9.INP'  ,STATUS='unknown')
      OPEN(UNIT=10 ,FILE='TAPE10.INP' ,STATUS='unknown')
C
C     THE FOLLOWING FILES ARE USED FOR OUTPUT DATA（下面的文件用于输出数据）
C
      OPEN(UNIT=6  ,FILE='TAPE6.OUT'  ,STATUS='unknown')
      OPEN(UNIT=7  ,FILE='TAPE7.OUT'  ,STATUS='unknown')
      OPEN(UNIT=11 ,FILE='TAPE11.OUT' ,STATUS='unknown')
      OPEN(UNIT=12 ,FILE='TAPE12.OUT' ,STATUS='unknown')
      OPEN(UNIT=13 ,FILE='TAPE13.OUT' ,STATUS='unknown')
      OPEN(UNIT=15 ,FILE='TAPE15.OUT' ,STATUS='unknown')
      OPEN(UNIT=16 ,FILE='TAPE16.OUT' ,STATUS='unknown')
      OPEN(UNIT=50 ,FILE='TAPE50.OUT' ,STATUS='unknown')
C
C     PRINT HEADER QAPAGE DATA AS FIRST PART OF OUTPUT（将QAPAGE数据作为输出的第一部分）
C
      CALL QAPAGE (6)
      CALL QAPAGE (11)
C
C  Q105F(8) IS USED TO PRODUCE MICROFICHE AT ORNL; A DUMMY ROUTINE IS
C  PRESENT IN THE RSIC VERSION TO ACCOMODATE THIS  （Q105F(8)用于在ORNL生成缩微胶片；在RSIC版本中存在一个虚拟例程以适应这一点）
      CALL Q105F(8)
C INITIALIZE PAGE COUNTER （ 初始化页计数器）
      NPAGE=IPAGE(0)
C  ERRSET SUPRESSES ERROR MESSAGES FROM UNDERFLOWS, WHICH ARE NORMAL
C  OCCURRENCES IN ORIGEN2 CALCULATIONS（ERRSET 抑制来自UNDERFLOWS的错误消息，这在 ORIGEN2 计算中是正常现象）
      CALL ERRSET(208,300,-1,1)
C  NUDOC PRINTS ORIGEN2 DOCUMENTATION DATA（NUDOC 打印 ORIGEN2 文档数据）
      CALL NUDOC(6)
C  NUMBER OF STORAGE VECTORS（存储向量的数量）
      LX= JJJJ
C  NUMBER OF OUTPUT VECTORS+1（输出向量的数量 + 1）
      MX= AAAA
C  MX- THE NUMBER OF CROSS-SECTION REACTIONS PER NUCLIDE (USUALLY 6)（ MX - 每个核素的截面反应数量（通常为 6））
      LC= CCCC
C  MAX # OF ACTIVATION PRODUCTS（最大活化产物数量）
      ILMAX= MMMM
C  MAX # OF ACTINIDES（最大锕系元素数量）
      IAMAX= GGGG
C  MAX # OF FISSION PRODUCTS（裂变产物的最大数量）
      IFMAX= FFFF
C  MAX TOTAL # OF NUCLIDES（核素的最大总数）
      ITMAX= BBBB
C  MAX # OF REACTIONS IN MATRIX A（矩阵 A 中反应的最大数量）
      IZMAX= DDDD
C  MAX # OF PHOTONS（光子的最大数量）
      IPHMAX=LLLL
C  MAX # OF REACTIONS IN LONG-LIVED MATRIX AP（长寿命矩阵中反应的最大数量）
      IAPMAX=IIII
C  MAX # OF FISSION PRODUCT YIELDS（裂变产物产额的最大数量。）
      IFYMAX=EEEE
C  MAX # OF NATURAL ABUNDANCES（自然丰度的最大数量。）
      NABMAX=KKKK
C  MAX # OF ACTINIDES WITH VARIABLE CROSS SECTIONS AND FISSION
C  PRODUCT YIELDS(变截面锕系核素和裂变产物产额的最大值)
      ICNMAX=OOOO
C  MAX # OF FISSION PRODUCTS ATTACHED TO ACTINIDES WITH VARIABLE
C  FISSION CROSS SECTIONS(具有可变裂变截面的锕系同位素附加的裂变产物的最大数量。)
      IFD= FFFF
C  MAX # OF VARIABLE MULTIPLIERS IN A SET IN RMULV(Rmulv中集合中变量乘子的最大值)
      LAM=NNNN
C  NEUTRONS PER NEUTRON-INDUCED FISSION:
C  0= THERMAL SPECTRUM, 1= FAST SPECTRUM(每个中子诱导裂变所产生的中子数。1是快光谱0是热光谱)
      NYTF=1
      NYTF=0
C  LISTIT READS CARD INPUT FROM UNIT 5, LISTS IT ON UNIT 6, AND
C  WRITES IT ON UNIT 50.  UNIT 50 IS THEN REWOUND AND ORIGEN2
C  READS THE DATA FROM UNIT 50.（这段注释描述了LISTIT功能的作用。它从输入单元5读取卡片输入，将其列印在输出单元6上，并将其写入单元50。然后，单元50被倒带，ORIGEN2从单元50读取数据。）
      CALL LISTIT(5,6,50)
      REWIND 50
C MAIN1 HANDLES THE MISCELLANEOUS INITIALIZATION DATA（主程序MAIN1负责处理各种杂项初始化数据）
    1 CALL MAIN1(NYTF,SFNU,ALPHN,NUCAN,NUCSFU,NY,YY,ANMUL,ANEXP)
C MAIN2 READS THE ORIGEN2 COMMANDS（程序MAIN2负责读取ORIGEN2命令）
    2 CALL MAIN2(NSTP)
C MAIN3 EXECUTES THE ORIGEN2 COMMANDS(主程序MAIN3执行ORIGEN2命令。)
    3 CALL MAIN3(
     $ LONG,STTFPB,ISTOTI,IS,RSTOTI,  LX, MX, LC,IFD,
     $NUCAB,NON0,KD,LOC,NGF,NGN,NGR,NYIELD,NONP,NQ,LOCP,MMAX,KAP,
     $LOCA,NFUDFP,    CIMN,CSUM,                   S,
     $NUCL,Q,FG,TOCAP,GENNEU,ALPHAN,SPONF,SFNU,FISS,AMPC,WMPC,XSTORE,
     $DIS,B,GGR,YIELD,A     ,XP,XPAR,XTEMP,D,AP,COEFF,NPROD,   XNEW,
     $ALPHN,NUCAN,NUCSFU,NY,YY,FFSF,FFA,ABUND,RMULV,LAM)
C  THIS "GO TO" PROVIDES THE MECHANISM FOR EXECUTING MULTIPLE
C  PROBLEMS WITHIN A SINGLE JOB OR EXTENDING THE CALCULATIONS AD
C  INFINITUM（这个 "GO TO" 语句提供了在单个作业中执行多个问题或无限扩展计算的机制）
      GO TO (1,2,3,4),NSTP
    4 CONTINUE
      CALL Q105F(8)
      STOP 'ORIGEN2 V2.1 (8-1-91) Execution Completed'
      END
      SUBROUTINE MAIN1(NYTF, SFNU, ALPHN, NUCAN, NUCSFU, NY, YY, ANMUL, MAI   10
     *     ANEXP)                                                       MAI   20
      DOUBLE PRECISION UNIT, TENS, DNUM                                 MAI   30
      DIMENSION SFNU(1), ALPHN(1), NUCAN(1), NUCSFU(1), NY(1), YY(1)    MAI   40
      COMMON /BDCOM/ NUCLA(100), TNPNIF(100), NUCLB(020), ALN(020),     MAI   50
     *     NUCLC(100), SFNDAT(100), CNMUL, CNEXP, FNPNIF(100)           MAI   60
      COMMON /FRECOV/ FRECG(10,20), NGROUP(100), FREC(10,100)           MAI   70
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  MAI   80
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           MAI   90
     *     TFPBAV(4), IPHMAX                                            MAI  100
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    MAI  110
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  MAI  120
     *     ITITLE(18)                                                   MAI  130
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      MAI  140
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           MAI  150
C     LEVEL 2,SFNU,ALPHN,NUCAN,NUCSFU,NY,YY,FRECG（这些参数级别为2级）                       MAI  160
      DO 10 I=1,IAMAX                                                   MAI  170
         NY(I) = 0                                                      MAI  180
         NUCAN(I) = 0                                                   MAI  190
         NUCSFU(I) = 0                                                  MAI  200
   10 CONTINUE                                                          MAI  210
      INPT = IUNIT                                                      MAI  220
      KO = KOUT                                                         MAI  230
      IF (IAMAX.LT.10) GO TO 110                                        MAI  240
C    TRANSFER NEUTRON YIELDS PER NEUTRON-INDUCED FISSION(每次中子诱发裂变的转移中子产额)                MAI  250
      DO 20 I=1,100                                                     MAI  260
         IF (NUCLA(I).LT.0) GO TO 30                                    MAI  270
         NY(I) = NUCLA(I)                                               MAI  280
         YY(I) = FNPNIF(I)                                              MAI  290
         IF (NYTF.NE.1) YY(I) = TNPNIF(I)                               MAI  300
   20 CONTINUE                                                          MAI  310
   30 CONTINUE                                                          MAI  320
      II = I - 1                                                        MAI  330
      IYY = (II/6) + 1                                                  MAI  340
      IF (MOD(II,6).EQ.0) IYY = IYY - 1                                 MAI  350
      IF (I2.GT.II) I2 = II                                             MAI  360
      WRITE (IOUT,99995) IOUT, NPAGE                                    MAI  370
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 01, 01, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI  380
      NPAGE = IPAGE(1)                                                  MAI  390
      WRITE (KOUT,99994)                                                MAI  400
      DO 40 I=1,IYY                                                     MAI  410
         I1 = 6*(I-1) + 1                                               MAI  420
         I2 = I1 + 5                                                    MAI  430
         WRITE (KOUT,99993) (NY(IY),YY(IY),IY=I1,I2)                    MAI  440
   40 CONTINUE                                                          MAI  450
C    TRANSFER EMPIRICAL (ALPHA,N) DATA   (转移经验数据)                               MAI  470
      DO 50 I=1,50                                                      MAI  480
         IF (NUCLB(I).LT.0) GO TO 60                                    MAI  490
         NUCAN(I) = NUCLB(I)                                            MAI  500
         ALPHN(I) = ALN(I)                                              MAI  510
   50 CONTINUE                                                          MAI  520
   60 IJ = I - 1                                                        MAI  530
      ANMUL = CNMUL                                                     MAI  540
      ANEXP = CNEXP                                                     MAI  550
      IF(IJ.LT.1)GO TO 71                                               MAI  551
      IAN = IJ/6 + 1                                                    MAI  560
      IF (MOD(IJ,6).EQ.0) IAN = IAN - 1                                 MAI  570
      WRITE (IOUT,99995) IOUT, NPAGE                                    MAI  590
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 01, 02, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI  600
      NPAGE = IPAGE(1)                                                  MAI  610
      WRITE (KO,99984)                                                  MAI  620
      DO 70 I=1,IAN                                                     MAI  630
         I1 = 6*(I-1) + 1                                               MAI  640
         I2 = I1 + 5                                                    MAI  650
         IF (I2.GT.IJ) I2 = IJ                                          MAI  660
         WRITE (KO,99992) (NUCAN(IA),ALPHN(IA),IA=I1,I2)                MAI  670
   70 CONTINUE                                                          MAI  680
   71 CONTINUE                                                          MAI  700
C    TRANSFER EMPIRICAL SF NEUTRON YIELDS     （转移经验数据中的裂变中子产额）                          MAI  710
      DO 80 I=1,100                                                     MAI  720
         IF (NUCLC(I).LT.0) GO TO 90                                    MAI  730
         NUCSFU(I) = NUCLC(I)                                           MAI  740
         SFNU(I) = SFNDAT(I)                                            MAI  750
   80 CONTINUE                                                          MAI  760
   90 IK = I - 1                                                        MAI  770
      ISF = IK/6 + 1                                                    MAI  780
      IF (MOD(IK,6).EQ.0) ISF = ISF - 1                                 MAI  790
      WRITE (IOUT,99995) IOUT, NPAGE                                    MAI  810
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 01, 03, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI  820
      NPAGE = IPAGE(1)                                                  MAI  830
      WRITE (KO,99991)                                                  MAI  840
      DO 100 I=1,ISF                                                    MAI  850
         I1 = 6*(I-1) + 1                                               MAI  860
         I2 = I1 + 5                                                    MAI  870
         IF (I2.GT.IK) I2 = IK                                          MAI  880
         WRITE (KO,99990) (NUCSFU(IV),SFNU(IV),IV=I1,I2)                MAI  890
  100 CONTINUE                                                                  
      WRITE (KO,99998)                                                  MAI  920
  110 CONTINUE                                                          MAI  930
C    READ INDIVIDUAL ELEMENT FRACTIONAL RECOVERIES TO OVERRIDE DEFAULT  MAI  940
C    VALUES（读取单个元素的分数回收率以覆盖默认值）                                                             MAI  950
      DO 120 I=1,1000                                                   MAI  960
         NELEM = IREAD(1)                                               MAI  970
         IF (NELEM.LT.0) GO TO 130                                      MAI  980
         NSET = IREAD(0)                                                MAI  990
         RECOV = READ(0)                                                MAI 1000
         FREC(NSET,NELEM) = RECOV                                       MAI 1010
  120 CONTINUE                                                          MAI 1020
  130 CONTINUE                                                          MAI 1030
      WRITE (IOUT,99995) IOUT, NPAGE                                    MAI 1050
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 01, 04, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI 1060
      NPAGE = IPAGE(1)                                                  MAI 1070
      WRITE (KO,99989)                                                  MAI 1080
      DO 150 NE=1,99                                                    MAI 1090
         IF (MOD(NE,50).NE.0) GO TO 140                                 MAI 1100
         WRITE (IOUT,99995) IOUT, NPAGE                                 MAI 1110
      CALL IDENTIFY (IOUT, 0)                                           PC
         NPAGE = IPAGE(1)                                               MAI 1120
         WRITE (KOUT,99989)                                             MAI 1130
  140    WRITE (KOUT,99988) NE, (FREC(ISET,NE),ISET=1,10)               MAI 1140
  150 CONTINUE                                                          MAI 1150
C    READ GROUP FRAC RECOV TO OVERRIDE DEFAULT VALUES                   MAI 1160
      DO 160 I=1,200                                                    MAI 1170
         NGRP = IREAD(1)                                                MAI 1180
         IF (NGRP.LT.0) GO TO 170                                       MAI 1190
         NSET = IREAD(0)                                                MAI 1200
         RECOV = READ(0)                                                MAI 1210
         FRECG(NSET,NGRP) = RECOV                                       MAI 1220
  160 CONTINUE                                                          MAI 1230
  170 CONTINUE                                                          MAI 1240
      WRITE (IOUT,99995) IOUT, NPAGE                                    MAI 1260
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 01, 05, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI 1270
      NPAGE = IPAGE(1)                                                  MAI 1280
      WRITE (KO,99987)                                                  MAI 1290
      DO 180 I=1,20                                                     MAI 1300
         WRITE (KO,99986) I, (FRECG(IGFR,I),IGFR=1,10)                  MAI 1310
  180 CONTINUE                                                          MAI 1320
C    READ GROUP ASSIGNMENTS TO OVERRIDE DEFAULT VALUES                  MAI 1330
      DO 190 I=1,100                                                    MAI 1340
         NELEM = IREAD(1)                                               MAI 1350
         IF (NELEM.LT.0) GO TO 200                                      MAI 1360
         NGRP = IREAD(0)                                                MAI 1370
         NGROUP(NELEM) = NGRP                                           MAI 1380
  190 CONTINUE                                                          MAI 1390
  200 CONTINUE                                                          MAI 1400
      WRITE (IOUT,99995) IOUT, NPAGE                                    MAI 1420
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 01, 06, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI 1430
      NPAGE = IPAGE(1)                                                  MAI 1440
      WRITE (KO,99997)                                                  MAI 1450
      DO 210 NROW=1,11                                                  MAI 1460
         I1 = 9*(NROW-1) + 1                                            MAI 1470
         I2 = I1 + 8                                                    MAI 1480
         WRITE (KOUT,99996) (I,NGROUP(I),I=I1,I2)                       MAI 1490
  210 CONTINUE                                                          MAI 1500
C    WRITE CHEMICAL TOXICITY DATA                                       MAI 1510
      WRITE (IOUT,99995) IOUT, NPAGE                                    MAI 1530
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 01, 07, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI 1540
      NPAGE = IPAGE(1)                                                  MAI 1550
      WRITE (KOUT,99985)                                                MAI 1560
      DO 220 I=1,17                                                     MAI 1570
         I1 = 6*(I-1) + 1                                               MAI 1580
         I2 = I1 + 5                                                    MAI 1590
         IF (I2.GT.99) I2 = 99                                          MAI 1600
         WRITE (KOUT,99992) (IA,CEMTOX(IA),IA=I1,I2)                    MAI 1610
  220 CONTINUE                                                          MAI 1620
      ISKIP = 0                                                         MAI 1640
      RETURN                                                            MAI 1650
99998 FORMAT (54H0THE REMAINING NEUTRON YIELDS ARE CALCULATED FROM THE ,MAI 1700
     *     44HEQUATION: NEUT/FISS=(2.84+0.1225*(AT WT-244))             MAI 1710
99997 FORMAT (50H ELEMENT ASSIGNMENTS TO FRACTIONAL RECOVERY GROUPS, /, MAI 1720
     *     1H0, 9(14H ELEM GROUP   ))                                   MAI 1730
99996 FORMAT (1H , 9(2X, I2, 3X, I2, 5X))                               MAI 1740
99995 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        MAI 1750
99994 FORMAT (42H NEUTRON YIELD PER NEUTRON-INDUCED FISSION, /, /,      MAI 1760
     *     6(4X, 8H NUCLIDE, 5X, 5HYIELD))                              MAI 1770
99993 FORMAT (1H0, 6(6X, I6, 4X, F6.4))                                 MAI 1780
99992 FORMAT (1H , 6(I6, 3X, 1PE8.2, 3X))                               MAI 1790
99991 FORMAT (48H SPONTANEOUS FISSION NEUTRON YIELD, NEUT/FISSION, /,   MAI 1800
     *     /, 1H , 6(7HNUCLIDE, 13X))                                   MAI 1810
99990 FORMAT (1H , 6(I6, 4X, F6.4, 4X))                                 MAI 1820
99989 FORMAT (41H INDIVIDUAL ELEMENT FRACTIONAL RECOVERIES, /, /,       MAI 1830
     *     12H   ELEMENT  , 10(12H SET NUMBER ), /, 1H , T18, 1H1, T30, MAI 1840
     *     1H2, T42, 1H3, T54, 1H4, T66, 1H5, T78, 1H6, T90, 1H7, T102, MAI 1850
     *     1H8, T114, 1H9, T126, 2H10)                                  MAI 1860
99988 FORMAT (1H , I6, 5X, 10(1X, F9.7, 2X))                            MAI 1870
99987 FORMAT (28H GROUP FRACTIONAL RECOVERIES, /, /, 13H    GROUP    ,  MAI 1880
     *     10(12H SET NUMBER ), /, 1H , T18, 1H1, T30, 1H2, T42, 1H3,   MAI 1890
     *     T54, 1H4, T66, 1H5, T78, 1H6, T90, 1H7, T102, 1H8, T114, 1H9,MAI 1900
     *     T126, 2H10)                                                  MAI 1910
99986 FORMAT (1H , I6, 5X, 10(1X, F9.7, 2X))                            MAI 1920
99985 FORMAT (1H , T40, 41HCHEMICAL TOXICITIES, GRAMS PER M**3 WATER,   MAI 1930
     *     /, /, 1H0, 6(7HNUCLIDE, 3X, 8HTOXICITY, 2X))                 MAI 1940
99984 FORMAT (39H (ALPHA,N) NEUTRON YIELD IN OXIDE FUELS, /, 1H , T10,  MAI 1950
     *     25HUNITS ARE NEUT/SEC-G ATOM, /, /, 1H , 6(7HNUCLIDE, 13X))  MAI 1960
      END                                                               MAI 1970
      SUBROUTINE MAIN2(NSTP)                                            MAI   10
      LOGICAL WRDFLG, LSCAN, LALLC, LNEG                                MAI   20
      INTEGER AR, AREAD                                                 MAI   30
      DOUBLE PRECISION DNUM, UNIT, TENS                                 MAI   40
      DIMENSION MAXCOM(32),CL(1),BP(1),BS(1)                            MAI   50
      COMMON /BDAT/ ZAP(130), KZAP(12), JZAP(6,26), NTB(3,4),           MAI   60
     *     INSTR0(32), RRCUT(30), NDI(11,7), NDL(6,5), NDO(9,5)         MAI   70
      COMMON /FUDG/ NFUD(7,20), FUDSIG(32,20), CBUP(32), NXSEC, NSM,    MAI   80
     *     NPT, NCO(2,20), FPBST(050,3), VECT(050), FOLD(20), FNEW(20)  MAI   90
      COMMON /INST/ INSTR(300), NINSTR(300), IFIND(300), NNTYP(32),     MAI  100
     *     LAB(300)                                                     MAI  110
      COMMON /MANE4I/ NINP(15,06), NIRR(150,04), NPRO(20,05),           MAI  120
     *     NCON(20,01), NDOL(02,04), NMOV(99,03), NADD(30,03),          MAI  130
     *     NHED(50,01), NWAC(09,02), NLPU(09,100), NPCH(54,03),         MAI  140
     *     NBUP(20), INUCL(5), NKEQ(03,05), NFAC(20,4), NGTO(10,2),     MAI  150
     *     NOFF(5)                                                              
      COMMON /MANE4R/ RADD(30), RIRR(150,02), RCUT(03,30), RMOV(99),    MAI  160
     *     XCOMP(5), RKEQ(3), RFAC(20), RGTO(10), ROFF(5)               MAI  170
      COMMON /MANE2I/ NOPT(20,80), NLIB(5,11), NLIP(5,3), NPHO(5,4)     MAI  180
      COMMON /MANE1L/ NTIT(20,18), NBAS(10,18), KHED(50,3)              MAI  190
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    MAI  200
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  MAI  210
     *     ITITLE(18)                                                   MAI  220
      COMMON /OUTOPT/ NOUT(20,4), NTO(80),MCOUT                         MAI  230
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      MAI  240
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           MAI  250
      EQUIVALENCE (BL,CL(1)),(AP,BP(1)),(AS,BS(1))                      MAI  251
C     LEVEL 2,NOPT,NINP                                                 MAI  260
C                                                                       MAI  270
C       MAXIMUM NUMBER OF EACH INSTRUCTION TYPE FOR OVERFLOW TEST       MAI  280
      DATA MAXCOM /1000,20,10,20,20,15,50,1,3,3,2,99,30,20,54,5,9,5,5,9,MAI  290
     *     150,0,0,20,20,20,20,20,1,1,10,5/                             MAI  300
      DATA CL /4H    /                                                  MAI  310
      DATA BP /4H'   /                                                  MAI  320
      DATA BS /4H*   /                                                  MAI  330
C                                                                       MAI  340
      MSUM = 1                                                          MAI  350
C    MAXIMUM ALLOWABLE NUMBER OF COMMANDS                               MAI  360
      INSMAX = 300                                                      MAI  370
      NINTYP = 32                                                       MAI  380
C    READ STRING OF INSTRUCTIONS DEFINING FLOWSHEET                     MAI  390
C                                                                       MAI  400
      NC = NTAB + 4                                                     MAI  410
      IF (NREC.LE.0) GO TO 20                                           MAI  420
      NTAB2 = NC - 1                                                    MAI  430
      DO 10 J=1,NTAB2                                                   MAI  440
         RCUT(1,J) = RRCUT(J)                                           MAI  450
   10 CONTINUE                                                          MAI  460
      ERR = RRCUT(NC)                                                   MAI  470
   20 CONTINUE                                                          MAI  480
      DO 30 N=1,INSMAX                                                  MAI  490
         NINSTR(N) = 0                                                  MAI  500
   30 CONTINUE                                                          MAI  510
      L = 1                                                             MAI  520
      WRITE (IOUT,99996) IOUT, NPAGE                                    MAI  530
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 02, 00, 00, 00, ITITLE, IOUT, NTOC, NREC)         MAI  540
      NPAGE = IPAGE(1)                                                  MAI  550
      DO 40 I=1,NINTYP                                                  MAI  560
         NNTYP(I) = 0                                                   MAI  570
   40 CONTINUE                                                          MAI  580
   50 AR = AREAD(1)                                                     MAI  590
      IF (L.GT.INSMAX) WRITE (15,99999) L, INSMAX                       MAI  600
      IF (MOD(L,50).NE.0) GO TO 60                                      MAI  610
      WRITE (IOUT,99996) IOUT, NPAGE                                    MAI  620
      CALL IDENTIFY (IOUT, 0)                                           PC
      NPAGE = IPAGE(1)                                                  MAI  630
   60 INSTR(L) = AR                                                     MAI  640
      LAB(L) = 0                                                        MAI  650
      NINSTR(L) = -1                                                    MAI  660
      WRITE (IOUT,99997) L, LAB(L), INSTR(L), CARD                      MAI  670
      DO 70 IJ=1,NINTYP                                                 MAI  680
         IF (INSTR(L).EQ.INSTR0(IJ)) NINSTR(L) = IJ                     MAI  690
   70 CONTINUE                                                          MAI  700
      IF (NINSTR(L).LT.0) WRITE (IOUT,99995) L, INSTR(L)                MAI  710
      MDUM = NINSTR(L)                                                  MAI  720
      M = MDUM                                                          MAI  730
      IF (M.EQ.22 .OR. M.EQ.23) M = 21                                  MAI  740
      IF (NNTYP(M).GE.MAXCOM(M)) WRITE (15,99998) L, INSTR(L),          MAI  750
     *     MAXCOM(MDUM), NNTYP(M), CARD                                 MAI  760
      GO TO (80, 90, 110, 130, 150, 170, 190, 240, 250, 290, 310, 320,  MAI  770
     *     340, 360, 370, 390, 410, 430, 870, 890, 920, 950, 980, 1000, MAI  780
     *     1020, 1030, 1040, 1090, 1100, 1120, 1110, 1111), MDUM        MAI  790
C    RDA--READ DESCRIPTIVE TEXT                                         MAI  800
   80 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI  810
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI  820
      L = L + 1                                                         MAI  830
      GO TO 50                                                          MAI  840
C    TIT--READ CASE TITLE                                               MAI  850
   90 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI  860
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI  870
      ICHR = 8                                                          MAI  880
      DO 100 J=1,18                                                     MAI  890
         NTIT(IFIND(L),J) = AREAD(0)                                    MAI  900
  100 CONTINUE                                                          MAI  910
      L = L + 1                                                         MAI  920
      GO TO 50                                                          MAI  930
C    BAS--READ BASIS FOR THIS CALCULATION                               MAI  940
  110 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI  950
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI  960
      ICHR = 8                                                          MAI  970
      DO 120 J=1,18                                                     MAI  980
         NBAS(IFIND(L),J) = AREAD(0)                                    MAI  990
  120 CONTINUE                                                          MAI 1000
      L = L + 1                                                         MAI 1010
      GO TO 50                                                          MAI 1020
C    FAC--CALCULATE MULTIPLICATION FACTORS FOR MOV AND ADD              MAI 1030
  130 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1040
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1050
      DO 140 J=1,4                                                      MAI 1060
         NFAC(IFIND(L),J) = IREAD(0)                                    MAI 1070
  140 CONTINUE                                                          MAI 1080
      RFAC(IFIND(L)) = READ(0)                                          MAI 1090
      L = L + 1                                                         MAI 1100
      GO TO 50                                                          MAI 1110
C    OUT--OUTPUT RESULTS                                                MAI 1120
  150 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1130
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1140
      DO 160 J=1,4                                                      MAI 1150
         NOUT(IFIND(L),J) = IREAD(0)                                    MAI 1160
  160 CONTINUE                                                          MAI 1170
      L = L + 1                                                         MAI 1180
      GO TO 50                                                          MAI 1190
C    INP--READ INITIAL COMPOSITIONS,CONT. REMOVAL,CONT. FEED            MAI 1200
  170 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1210
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1220
      DO 180 J=1,6                                                      MAI 1230
         NINP(IFIND(L),J) = IREAD(0)                                    MAI 1240
  180 CONTINUE                                                          MAI 1250
      L = L + 1                                                         MAI 1260
      GO TO 50                                                          MAI 1270
C    HED--ALLOWS VECTOR HEADINGS TO BE READ IN SEPARATELY               MAI 1280
  190 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1290
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1300
      NHED(IFIND(L),1) = IREAD(0)                                       MAI 1310
      J = ICHR                                                          MAI 1320
  200 IF (CARD(J).EQ.BL) GO TO 220                                      MAI 1330
      IF (CARD(J).EQ.AP .OR. CARD(J).EQ.AS) J = J + 1                   MAI 1340
      ICHR = J - 1                                                      MAI 1350
      DO 210 JK=1,3                                                     MAI 1360
         KHED(IFIND(L),JK) = AREAD(0)                                   MAI 1370
  210 CONTINUE                                                          MAI 1380
      GO TO 230                                                         MAI 1390
  220 J = J + 1                                                         MAI 1400
      IF (J.GT.68) GO TO 230                                            MAI 1410
      GO TO 200                                                         MAI 1420
  230 L = L + 1                                                         MAI 1430
      GO TO 50                                                          MAI 1440
C    REC--RECYCLE COUNTER                                               MAI 1450
  240 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1460
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1470
      NREC = 0                                                          MAI 1480
      L = L + 1                                                         MAI 1490
      GO TO 50                                                          MAI 1500
C    CUT--READ CUTOFF VALUES FOR SUMMARY TABLES & TRUNCATION ERROR      MAI 1510
  250 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1520
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1530
C    CUT IN MAIN2                                                       MAI 1540
      DO 260 J=1,NC                                                     MAI 1550
         RCUT(IFIND(L),J) = RRCUT(J)                                    MAI 1560
  260 CONTINUE                                                          MAI 1570
      DO 270 J=1,NC                                                     MAI 1580
         IICHR = ICHR                                                   MAI 1590
         NCT = IREAD(0)                                                 MAI 1600
         IF (ICHR.LE.IICHR) WRITE (IOUT,99994) CARD                     MAI 1610
         IF (NCT.LE.0) GO TO 280                                        MAI 1620
         RCUT(IFIND(L),NCT) = READ(0)                                   MAI 1630
  270 CONTINUE                                                          MAI 1640
  280 L = L + 1                                                         MAI 1650
      GO TO 50                                                          MAI 1660
C    KEQ-- MATCH K-INFINITIES                                           MAI 1670
  290 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1680
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1690
      DO 300 J=1,5                                                      MAI 1700
         NKEQ(IFIND(L),J) = IREAD(0)                                    MAI 1710
  300 CONTINUE                                                          MAI 1720
      RKEQ(IFIND(L)) = READ(0)                                          MAI 1730
      L = L + 1                                                         MAI 1740
      GO TO 50                                                          MAI 1750
C    DOL--DO LOOP                                                       MAI 1760
  310 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1770
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1780
      NDOL(IFIND(L),1) = IREAD(0)                                       MAI 1790
      NDOL(IFIND(L),2) = IREAD(0)                                       MAI 1800
      NDOL(IFIND(L),3) = 0                                              MAI 1810
      L = L + 1                                                         MAI 1820
      GO TO 50                                                          MAI 1830
C    MOV--COPY A VECTOR OVER ANOTHER VECTOR                             MAI 1840
  320 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1850
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1860
      DO 330 J=1,3                                                      MAI 1870
         NMOV(IFIND(L),J) = IREAD(0)                                    MAI 1880
  330 CONTINUE                                                          MAI 1890
      RMOV(IFIND(L)) = READ(0)                                          MAI 1900
      L = L + 1                                                         MAI 1910
      GO TO 50                                                          MAI 1920
C    ADD--ADD ONE VECTOR TO ANOTHER VECTOR                              MAI 1930
  340 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 1940
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 1950
      DO 350 J=1,3                                                      MAI 1960
         NADD(IFIND(L),J) = IREAD(0)                                    MAI 1970
  350 CONTINUE                                                          MAI 1980
      RADD(IFIND(L)) = READ(0)                                          MAI 1990
      L = L + 1                                                         MAI 2000
      GO TO 50                                                          MAI 2010
C    BUP--DETERMINES BASIS FOR BURNUP CALCULATION                       MAI 2020
  360 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 2030
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 2040
      L = L + 1                                                         MAI 2050
      GO TO 50                                                          MAI 2060
C    PCH-- PUNCH OUTPUT VECTOR                                          MAI 2070
  370 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 2080
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 2090
      DO 380 J=1,3                                                      MAI 2100
         NPCH(IFIND(L),J) = IREAD(0)                                    MAI 2110
  380 CONTINUE                                                          MAI 2120
      L = L + 1                                                         MAI 2130
      GO TO 50                                                          MAI 2140
C    LIP--SUPRESS PRINTING OF DATA LIBRARIES                            MAI 2150
  390 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 2160
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 2170
      DO 400 J=1,3                                                      MAI 2180
         NLIP(IFIND(L),J) = IREAD(0)                                    MAI 2190
  400 CONTINUE                                                          MAI 2200
      L = L + 1                                                         MAI 2210
      GO TO 50                                                          MAI 2220
C    WAC--PLACE CONC VECTOR X CONT REMOVAL RATE IN B(I) T0 ACCUM WASTE  MAI 2230
  410 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 2240
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 2250
      DO 420 J=1,2                                                      MAI 2260
         NWAC(IFIND(L),J) = IREAD(0)                                    MAI 2270
  420 CONTINUE                                                          MAI 2280
      L = L + 1                                                         MAI 2290
      GO TO 50                                                          MAI 2300
C    LIB-- READ DECAY AND X-SEC DATA                                    MAI 2310
  430 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 2320
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 2330
      DO 440 J=1,11                                                     MAI 2340
         NLIB(IFIND(L),J) = IREAD(0)                                    MAI 2350
  440 CONTINUE                                                          MAI 2360
C    STORE VARIABLE XSEC INFO FOR REACTOR OF INTEREST                   MAI 2370
      LFUD = IREAD(0)                                                   MAI 2380
      NXSEC = 0                                                         MAI 2390
      IF (LFUD.GT.0) GO TO 450                                          MAI 2400
      L = L + 1                                                         MAI 2410
      GO TO 50                                                          MAI 2420
  450 CONTINUE                                                          MAI 2430
      GO TO (460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560,     MAI 2440
     *     570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680,  MAI 2450
     *     690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800,  MAI 2460
     *     810, 820, 830, 840, 850, 851, 852), LFUD                     ORIGEN21
  460 CALL XSEC01(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2480
      GO TO 860                                                         MAI 2490
  470 CALL XSEC02(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2500
      GO TO 860                                                         MAI 2510
  480 CALL XSEC03(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2520
      GO TO 860                                                         MAI 2530
  490 CALL XSEC04(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2540
      GO TO 860                                                         MAI 2550
  500 CALL XSEC05(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2560
      GO TO 860                                                         MAI 2570
  510 CALL XSEC06(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2580
      GO TO 860                                                         MAI 2590
  520 CALL XSEC07(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2600
      GO TO 860                                                         MAI 2610
  530 CALL XSEC08(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2620
      GO TO 860                                                         MAI 2630
  540 CALL XSEC09(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2640
      GO TO 860                                                         MAI 2650
  550 CALL XSEC10(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2660
      GO TO 860                                                         MAI 2670
  560 CALL XSEC11(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2680
      GO TO 860                                                         MAI 2690
  570 CALL XSEC12(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2700
      GO TO 860                                                         MAI 2710
  580 CALL XSEC13(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2720
      GO TO 860                                                         MAI 2730
  590 CALL XSEC14(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2740
      GO TO 860                                                         MAI 2750
  600 CALL XSEC15(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2760
      GO TO 860                                                         MAI 2770
  610 CALL XSEC16(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2780
      GO TO 860                                                         MAI 2790
  620 CALL XSEC17(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2800
      GO TO 860                                                         MAI 2810
  630 CALL XSEC18(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2820
      GO TO 860                                                         MAI 2830
  640 CALL XSEC19(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2840
      GO TO 860                                                         MAI 2850
  650 CALL XSEC20(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2860
      GO TO 860                                                         MAI 2870
  660 CALL XSEC21(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2880
      GO TO 860                                                         MAI 2890
  670 CALL XSEC22(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2900
      GO TO 860                                                         MAI 2910
  680 CALL XSEC23(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2920
      GO TO 860                                                         MAI 2930
  690 CALL XSEC24(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2940
      GO TO 860                                                         MAI 2950
  700 CALL XSEC25(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2960
      GO TO 860                                                         MAI 2970
  710 CALL XSEC26(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 2980
      GO TO 860                                                         MAI 2990
  720 CALL XSEC27(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3000
      GO TO 860                                                         MAI 3010
  730 CALL XSEC28(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3020
      GO TO 860                                                         MAI 3030
  740 CALL XSEC29(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3040
      GO TO 860                                                         MAI 3050
  750 CALL XSEC30(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3060
      GO TO 860                                                         MAI 3070
  760 CALL XSEC31(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3080
      GO TO 860                                                         MAI 3090
  770 CALL XSEC32(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3100
      GO TO 860                                                         MAI 3110
  780 CALL XSEC33(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3120
      GO TO 860                                                         MAI 3130
  790 CALL XSEC34(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3140
      GO TO 860                                                         MAI 3150
  800 CALL XSEC35(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3160
      GO TO 860                                                         MAI 3170
  810 CALL XSEC36(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3180
      GO TO 860                                                         MAI 3190
  820 CALL XSEC37(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       MAI 3200
      GO TO 860                                                         MAI 3210
  830 CALL XSEC38(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       ORIGEN21
      GO TO 860                                                         ORIGEN21
  840 CALL XSEC39(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       ORIGEN21
      GO TO 860                                                         ORIGEN21
  850 CALL XSEC40(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       ORIGEN21
      GO TO 860                                                         ORIGEN21
  851 CALL XSEC41(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       ORIGEN21
      GO TO 860                                                         ORIGEN21
  852 CALL XSEC42(NFUD, FUDSIG, CBUP, NXSEC, NPT)                       ORIGEN21
  860 CONTINUE                                                          MAI 3270
      L = L + 1                                                         MAI 3280
      GO TO 50                                                          MAI 3290
C    PHO--READ PHOTON LIBRARY                                           MAI 3300
  870 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 3310
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 3320
      DO 880 J=1,4                                                      MAI 3330
         NPHO(IFIND(L),J) = IREAD(0)                                    MAI 3340
  880 CONTINUE                                                          MAI 3350
      L = L + 1                                                         MAI 3360
      GO TO 50                                                          MAI 3370
C    LPU--READ SUBSTITUTE DATA CARDS                                    MAI 3380
  890 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 3390
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 3400
      DO 900 I=1,100                                                    MAI 3410
         IICHR = ICHR                                                   MAI 3420
         NLPU(IFIND(L),I) = IREAD(0)                                    MAI 3430
         IF (ICHR.LE.IICHR) WRITE (IOUT,99994) CARD                     MAI 3440
         IF (NLPU(IFIND(L),I).LT.0) GO TO 910                           MAI 3450
  900 CONTINUE                                                          MAI 3460
  910 L = L + 1                                                         MAI 3470
      GO TO 50                                                          MAI 3480
C    IRF IRRADIATION; FLUX GIVEN                                        MAI 3490
  920 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 3500
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 3510
      DO 930 J=1,2                                                      MAI 3520
         RIRR(IFIND(L),J) = READ(0)                                     MAI 3530
  930 CONTINUE                                                          MAI 3540
      DO 940 J=1,4                                                      MAI 3550
         NIRR(IFIND(L),J) = IREAD(0)                                    MAI 3560
  940 CONTINUE                                                          MAI 3570
      L = L + 1                                                         MAI 3580
      GO TO 50                                                          MAI 3590
C    IRP--IRRADIATION; SPECIFIC POWER GIVEN                             MAI 3600
  950 NNTYP(21) = NNTYP(21) + 1                                         MAI 3610
      IFIND(L) = NNTYP(21)                                              MAI 3620
      DO 960 J=1,2                                                      MAI 3630
         RIRR(IFIND(L),J) = READ(0)                                     MAI 3640
  960 CONTINUE                                                          MAI 3650
      DO 970 J=1,4                                                      MAI 3660
         NIRR(IFIND(L),J) = IREAD(0)                                    MAI 3670
  970 CONTINUE                                                          MAI 3680
      L = L + 1                                                         MAI 3690
      GO TO 50                                                          MAI 3700
C    DEC--DECAY                                                         MAI 3710
  980 NNTYP(21) = NNTYP(21) + 1                                         MAI 3720
      IFIND(L) = NNTYP(21)                                              MAI 3730
      RIRR(IFIND(L),1) = READ(0)                                        MAI 3740
      DO 990 J=1,4                                                      MAI 3750
         NIRR(IFIND(L),J) = IREAD(0)                                    MAI 3760
  990 CONTINUE                                                          MAI 3770
      L = L + 1                                                         MAI 3780
      GO TO 50                                                          MAI 3790
C    PRO--REPROCESSING                                                  MAI 3800
 1000 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 3810
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 3820
      DO 1010 J=1,4                                                     MAI 3830
         NPRO(IFIND(L),J) = IREAD(0)                                    MAI 3840
 1010 CONTINUE                                                          MAI 3850
      L = L + 1                                                         MAI 3860
      GO TO 50                                                          MAI 3870
C    OPTL-- PRINTOUT OPTIONS FOR LITE ELEMENTS                          MAI 3880
 1020 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 3890
      N1 = 1                                                            MAI 3900
      N2 = NTAB                                                         MAI 3910
      GO TO 1050                                                        MAI 3920
C    OPTA-- PRINTOUT OPTIONS FOR ACTINIDE ELEMENTS                      MAI 3930
 1030 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 3940
      N1 = NTAB + 1                                                     MAI 3950
      N2 = 2*NTAB                                                       MAI 3960
      GO TO 1050                                                        MAI 3970
C    OPTF-- PRINTOUT OPTIONS FOR FISSION PRODUCT ELEMENTS               MAI 3980
 1040 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 3990
      N1 = 2*NTAB + 1                                                   MAI 4000
      N2 = 3*NTAB                                                       MAI 4010
 1050 IFIND(L) = NNTYP(NINSTR(L))                                       MAI 4020
      DO 1060 J=N1,N2                                                   MAI 4030
      NOP=IREAD(0)                                                      MAI 4031
      IF(NOP.GT.-2)GO TO 1061                                           MAI 4032
      MCOUT=NOP                                                         MAI 4033
      NOP=IREAD(0)                                                      MAI 4034
 1061 NOPT(IFIND(L),J)=NOP                                              MAI 4035
 1060 CONTINUE                                                          MAI 4050
      IF (NOPT(IFIND(L),N1).NE.-1.AND. MSUM.GT.0) GO TO 1080            MAI 4060
      MSUM = -1                                                         MAI 4070
      DO 1070 J=N1,N2                                                   MAI 4080
         NOPT(IFIND(L),J) = 8                                           MAI 4090
 1070 CONTINUE                                                          MAI 4100
      NOPT(IFIND(L),N1+4) = 7                                           MAI 4110
 1080 CONTINUE                                                          MAI 4120
      L = L + 1                                                         MAI 4130
      GO TO 50                                                          MAI 4140
C    CON--CONTINUE CARD                                                 MAI 4150
 1090 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 4160
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 4170
      LAB(L) = IREAD(0)                                                 MAI 4180
      NCON(IFIND(L),1) = LAB(L)                                         MAI 4190
      L = L + 1                                                         MAI 4200
      GO TO 50                                                          MAI 4210
C    STP--TEMPORARY END OF INSTRUCTION SET                              MAI 4220
 1100 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 4230
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 4240
      NSTP = IREAD(0)                                                   MAI 4250
      L = L + 1                                                         MAI 4260
      GO TO 1120                                                        MAI 4270
C    GTO--GO TO A COMMAND & EXECUTE THOSE THAT FOLLOW DOWN TO A         MAI 4280
C    SPECIFIED CON COMMAND                                              MAI 4290
 1110 NNTYP(NINSTR(L)) = NNTYP(NINSTR(L)) + 1                           MAI 4300
      IFIND(L) = NNTYP(NINSTR(L))                                       MAI 4310
      NGTO(IFIND(L),1) = IREAD(0)                                       MAI 4320
      NGTO(IFIND(L),2) = IREAD(0)                                       MAI 4330
      RGTO(IFIND(L)) = READ(0)                                          MAI 4340
      L = L + 1                                                         MAI 4350
      GO TO 50                                                          MAI 4360
C    TOF--OFFSET TIME BETWEEN DISCHARGE AND REPROCESSING (MAKES                 
C    WASTE DECAY TIMES RELATIVE TO TIME OF DISCHARGE)                           
 1111 NNTYP(NINSTR(L))=NNTYP(NINSTR(L))+1                                       
      IFIND(L)=NNTYP(NINSTR(L))                                                 
      ROFF(IFIND(L))=READ(0)                                                    
      NOFF(IFIND(L))=IREAD(0)                                                   
      L=L+1                                                                     
      GO TO 50                                                                  
C    END--END OF PROGRAM                                                MAI 4370
 1120 CONTINUE                                                          MAI 4380
      RETURN                                                            MAI 4390
99999 FORMAT (38H **WARNING 9612** NUMBER OF COMMANDS =, I5, 8H IS GREA,MAI 4400
     *     20HTER THAN ALLOWABLE =, I5)                                 MAI 4410
99998 FORMAT (1H , 45H**WARNING 9613: TOO MANY OF A COMMAND TYPE** ,    MAI 4420
     *     10HCOMMAND #=, I4, 14H COMMAND TYPE=, A4, 11H # ALLOWED=,    MAI 4430
     *     I4, 13H # SPECIFIED=, I4, /1H , 20X, 80A1)                   MAI 4440
99997 FORMAT (1H , I4, 1X, I4, 1X, A4, 1X, 1H*, 80A1, 1H*)              MAI 4450
99996 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        MAI 4460
99995 FORMAT (1H0, 38HINCORRECT INSTRUCTION; INSTRUCTION # =, I5,       MAI 4470
     *     18HINSTRUCTION TYPE =, A4)                                   MAI 4480
99994 FORMAT (1H , 15X, 1H*, 80A1, 1H*)                                 MAI 4490
      END                                                               MAI 4500
      SUBROUTINE MAIN3(LONG, STTFPB, ISTOTI, IS, RSTOTI, LX, MX, LC,    MAI   10
     *     IFD, NUCAB, NON0, KD, LOC, NGF, NGN, NGR, NYIELD, NONP, NQ,  MAI   20
     *     LOCP, MMAX, KAP, LOCA, NFUDFP, CIMN, CSUM, S, NUCL, Q, FG,   MAI   30
     *     TOCAP, GENNEU, ALPHAN, SPONF, SFNU, FISS, AMPC, WMPC,        MAI   40
     *     XSTORE, DIS, B, GGR, YIELD, A, XP, XPAR, XTEMP, D, AP,       MAI   50
     *     COEFF, NPROD, XNEW, ALPHN, NUCAN, NUCSFU, NY, YY, FFSF, FFA, MAI   60
     *     ABUND, RMULV, LAM)                                           MAI   70
      LOGICAL LONG                                                      MAI   80
      INTEGER*2       NON0,KD,LOC,NGF,NGN,NGR,NYIELD,NONP,NQ,LOCP,MMAX, MAI   90
     *     KAP,LOCA,NFUDFP                                              MAI  100
      INTEGER BASIS                                                     MAI  110
      DIMENSION LOCA(1), NFUDFP(IFD,1)                                  MAI  120
      DOUBLE PRECISION CIMN, CSUM                                       MAI  130
      DOUBLE PRECISION DNUM, UNIT, TENS                                 MAI  140
      DIMENSION S(1),INFE(3)                                            MAI  150
      DIMENSION STTFPB(LX,1), ISTOTI(LX,1), IS(1), RSTOTI(1)            MAI  160
      DIMENSION FLSTOR(150), AA(4), TOCS(150), NVECST(150), POSTOR(150) MAI  170
      DIMENSION YY(1), NY(1), FFA(1), FFSF(1), ABUND(1), NUCSFU(1),     MAI  180
     *     NUCAN(1), ALPHN(1), RMULV(LAM,3)                             MAI  190
      DIMENSION LONG(1), TITLE(1), NUCAB(1), NON0(1), KD(1), LOC(1),    MAI  200
     *     NGF(1), NGN(1), NGR(1), NYIELD(1), NONP(1), NQ(1), LOCP(1),  MAI  210
     *     MMAX(1), KAP(1), CIMN(1), CSUM(1), NUCL(1), Q(1), FG(1),     MAI  220
     *     TOCAP(1), GENNEU(1), ALPHAN(1), SPONF(1), SFNU(1), FISS(1),  MAI  230
     *     AMPC(1), WMPC(1), XSTORE(LX,1), XNEW(MX,1), DIS(1), B(1),    MAI  240
     *     GGR(1), YIELD(1), A(1), XP(1), XPAR(1), XTEMP(1), D(1),      MAI  250
     *     AP(1), COEFF(LC,1), NPROD(LC,1)                              MAI  260
      DIMENSION FAC(20)                                                 MAI  270
      COMMON /FRECOV/ FRECG(10,20), NGROUP(100), FREC(10,100)           MAI  280
      COMMON /FUDG/ NFUD(7,20), FUDSIG(32,20), CBUP(32), NXSEC, NSM,    MAI  290
     *     NPT, NCO(2,20), FPBST(050,3), VECT(050), FOLD(20), FNEW(20)  MAI  300
      COMMON /INST/ INSTR(300), NINSTR(300), IFIND(300), NNTYP(32),     MAI  310
     *     LAB(300)                                                     MAI  320
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  MAI  330
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           MAI  340
     *     TFPBAV(4), IPHMAX                                            MAI  350
      COMMON /MAIN03/ NSTP, ANMUL, ANEXP, NABMAX, ICNMAX, IAPMAX, IFYMAXMAI  360
      COMMON /MANE4I/ NINP(15,06), NIRR(150,04), NPRO(20,05),           MAI  370
     *     NCON(20,01), NDOL(02,04), NMOV(99,03), NADD(30,03),          MAI  380
     *     NHED(50,01), NWAC(09,02), NLPU(09,100), NPCH(54,03),         MAI  390
     *     NBUP(20), INUCL(5), NKEQ(03,05), NFAC(20,4), NGTO(10,2),     MAI  400
     *     NOFF(5)                                                              
      COMMON /MANE4R/ RADD(30), RIRR(150,02), RCUT(03,30), RMOV(99),    MAI  410
     *     XCOMP(5), RKEQ(3), RFAC(20), RGTO(10), ROFF(5)               MAI  420
      COMMON /MANE2I/ NOPT(20,80), NLIB(5,11), NLIP(5,3), NPHO(5,4)     MAI  430
      COMMON /MANE1L/ NTIT(20,18), NBAS(10,18), KHED(50,3)              MAI  440
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            MAI  450
     *     INCYLD(8,4), MYLD(4)                                         MAI  460
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    MAI  470
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  MAI  480
     *     ITITLE(18)                                                   MAI  490
      COMMON /OUTOPT/ NOUT(20,4), NTO(80),MCOUT                         MAI  500
      COMMON /PROCSS/ MPROS, PRATE(8), NOPROS(8), NZPROS(8,20)          MAI  510
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      MAI  520
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           MAI  530
      COMMON /SZ/ NSIZE(25)                                             MAI  540
      COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)   MAI  550
      EXTERNAL GAMMA                                                    MAI  551
      DATA INFE/4HINFE,4HASIB,4HLE  /                                           
C      LEVEL 2,NUCAB,NON0,KD,LOC,NGF,NGN,NGR,NYIELD,NONP,NQ,LOCP,MMAX,  MAI  560
C    1  KAP,CIMN,CSUM,S,NUCL,Q,FG,TOCAP,GENNEU,ALPHAN,SPONF,SFNU,FISS,  MAI  570
C    2  AMPC,WMPC,XSTORE,DIS,B,GGR,YIELD,XP,XPAR,XTEMP,D,AP,COEFF,NPROD,MAI  580
C    4  XNEW,ALPHN,NUCAN,NUCSFU,NY,YY,FFSF,ABUND,FRECG,NOPT,NINP        MAI  590
C    5  ,A,LOCA,NFUDFP                                                  MAI  600
      IF (NREC.GE.0) GO TO 30                                           MAI  610
      DO 20 I=1,LX                                                      MAI  620
         IS(I) = 1                                                      MAI  630
         DO 10 J=1,10                                                   MAI  640
            STTFPB(I,J) = 0.0                                           MAI  650
   10    CONTINUE                                                       MAI  660
   20 CONTINUE                                                          MAI  670
   30 CONTINUE                                                          MAI  680
      NB = NTAB*3                                                       MAI  690
      IDM(1) = 0                                                        MAI  700
      DO 40 JJ=1,ITMAX                                                  MAI  710
         B(JJ) = 0.                                                     MAI  720
   40 CONTINUE                                                          MAI  730
      MPROS = 0                                                         MAI  740
      WRITE (16,99986)                                                  MAI  750
      MIRR = 0                                                          MAI  760
      RPF = 1.0                                                         MAI  770
      ANMO=1.0/ANMUL                                                    MAI  771
      ANEO=1.0/ANEXP                                                    MAI  772
      NJUMP = 0                                                         MAI  780
      BU = 1.                                                           MAI  790
      NREC = 0                                                          MAI  800
      TOFF=0.0                                                                  
      XOFF=0.0                                                                  
      L = 1                                                             MAI  810
   50 MDUM = NINSTR(L)                                                  MAI  820
      IF = IFIND(L)                                                     MAI  830
      WRITE (15,99984) L, IFIND(L), NNTYP(MDUM), INSTR(L)               MAI  840
      GO TO (60, 70, 90, 110, 270, 330, 720, 770, 780, 800, 860, 900,   MAI  850
     *     910, 920, 1020, 1100, 1120, 1140, 1180, 1200, 1210, 1230,    MAI  860
     *     1250, 1350, 1500, 1510, 1520, 1550, 1590, 1650, 1600,        MAI  870
     *     1641), MDUM                                                          
C    RDA--EXECUTION                                                     MAI  880
   60 L = L + 1                                                         MAI  890
      GO TO 50                                                          MAI  900
C    TIT--EXECUTION                                                     MAI  910
   70 DO 80 J=1,18                                                      MAI  920
         ITITLE(J) = NTIT(IFIND(L),J)                                   MAI  930
   80 CONTINUE                                                          MAI  940
      L = L + 1                                                         MAI  950
      GO TO 50                                                          MAI  960
C    BAS--EXECUTION                                                     MAI  970
   90 DO 100 J=1,18                                                     MAI  980
         BASIS(J) = NBAS(IFIND(L),J)                                    MAI  990
  100 CONTINUE                                                          MAI 1000
      L = L + 1                                                         MAI 1010
      GO TO 50                                                          MAI 1020
C    FAC--EXECUTION                                                     MAI 1030
  110 FA = 0.0                                                          MAI 1040
      FB = 0.0                                                          MAI 1050
      ILO = ILITE + 1                                                   MAI 1060
      IHI = ITOT                                                        MAI 1070
      LA = NFAC(IFIND(L),2)                                             MAI 1080
      LB = NFAC(IFIND(L),3)                                             MAI 1090
      LD = NFAC(IFIND(L),1)                                             MAI 1100
      LE = NFAC(IFIND(L),4)                                             MAI 1110
C    FOR VECTOR LA                                                      MAI 1120
      IF (RFAC(IFIND(L)).LT.ERR) GO TO 120                              MAI 1130
      FA = RFAC(IFIND(L))                                               MAI 1140
      GO TO 140                                                         MAI 1150
  120 DO 130 I=ILO,IHI                                                  MAI 1160
         KM = MOD(NUCL(I)/10,1000)                                      MAI 1170
         IF (LA.GT.0) FA = FA + XNEW(LA,I)*KM*0.001                     MAI 1180
         IF (LA.LT.0) FA = FA + XSTORE(-LA,I)*KM*0.001                  MAI 1190
  130 CONTINUE                                                          MAI 1200
  140 IF (RFAC(IFIND(L)).GT.(-ERR)) GO TO 150                           MAI 1210
      FB = (-RFAC(IFIND(L)))                                            MAI 1220
      GO TO 170                                                         MAI 1230
  150 DO 160 I=ILO,IHI                                                  MAI 1240
         KM = MOD(NUCL(I)/10,1000)                                      MAI 1250
         IF (LB.GT.0) FB = FB + XNEW(LB,I)*KM*0.001                     MAI 1260
         IF (LB.LT.0) FB = FB + XSTORE(-LB,I)*KM*0.001                  MAI 1270
  160 CONTINUE                                                          MAI 1280
  170 CONTINUE                                                          MAI 1290
      GO TO (180, 190, 200, 210, 220, 230, 240, 250, 251), LE           MAI 1300
  180 FAC(LD) = FA + FB                                                 MAI 1310
      GO TO 260                                                         MAI 1320
  190 FAC(LD) = FA - FB                                                 MAI 1330
      GO TO 260                                                         MAI 1340
  200 FAC(LD) = FA*FB                                                   MAI 1350
      GO TO 260                                                         MAI 1360
  210 FAC(LD) = FA/FB                                                   MAI 1370
      GO TO 260                                                         MAI 1380
  220 FAC(LD) = FA                                                      MAI 1390
      GO TO 260                                                         MAI 1400
  230 FAC(LD) = FB                                                      MAI 1410
      GO TO 260                                                         MAI 1420
  240 FAC(LD) = 1.0/FA                                                  MAI 1430
      GO TO 260                                                         MAI 1440
  250 FAC(LD) = 1.0/FB                                                  MAI 1450
      GO TO 260                                                                 
  251 DO 252 I=1,ITOT                                                           
      XP(I)=0.0                                                                 
      IF(LA.GT.0.AND.XNEW(LA,I).NE.0.0)XP(I)=1.0/XNEW(LA,I)                     
      IF(LA.LT.0.AND.XSTORE(-LA,I).NE.0.0)XP(I)=1.0/XSTORE(-LA,I)               
  252 CONTINUE                                                                  
      XP(ITOT+1)=XP(1)                                                          
      FAC(LD)=XP(1)                                                             
  260 L = L + 1                                                         MAI 1460
      WRITE (15,99999) LD, FAC(LD)                                      MAI 1470
      GO TO 50                                                          MAI 1480
C    OUT--EXECUTION                                                     MAI 1490
  270 NTEM = NOUT(IFIND(L),2)                                           MAI 1500
      IF (NREC.LT.2) GO TO 280                                          MAI 1510
      IF ((NREC/NTEM)*NTEM.NE.NREC) GO TO 320                           MAI 1520
  280 MOUT = NOUT(IFIND(L),1)                                           MAI 1530
      ISW = 0                                                           MAI 1540
      IF(MOUT.GT.0.AND.MCOUT.EQ.-2)GO TO 320                            MAI 1541
      IF(MOUT.LT.0.AND.MCOUT.EQ.-3)GO TO 320                            MAI 1542
      IF (MOUT.GT.0) GO TO 290                                          MAI 1550
      NTOC = NTOCB                                                      MAI 1560
      ISW = IOUT                                                        MAI 1570
      NPAGE = IPAGE(2)                                                  MAI 1580
      IOUT = 11                                                         MAI 1590
      JOUT = 11                                                         MAI 1600
      KOUT = 11                                                         MAI 1610
      MOUT = -MOUT                                                      MAI 1620
  290 CONTINUE                                                          MAI 1630
      CALL OUTPUT(L, NUCL, Q, FG, AMPC, WMPC, TOCAP, FISS, DIS, XNEW,   MAI 1640
     *     MX, S(1), S(1302), S(1318), S(1334), S(1350), S(1366), FFA,  MAI 1650
     *     S(1382))                                                     MAI 1660
      CALL NUTRON(XNEW, MX, NUCL, ALPHAN,        SFNU, S(1), S(16),     MAI 1670
     *     S(32), S(48), S(64), S(80), S(96))                           MAI 1680
      IF (NSIZE(20).LT.1) GO TO 300                                     MAI 1690
      MSTAR = NOUT(IFIND(L),4)                                          MAI 1700
      CALL GAMMA(NUCL, DIS, NGF, NGN, NGR, GGR, XNEW, MX, S(1), S(16),  MAI 1710
     *     S(252), S(268), S(284), S(300), MSTAR)                       MAI 1720
  300 CONTINUE                                                          MAI 1730
      IF (ISW.LT.1) GO TO 310                                           MAI 1740
      IOUT = ISW                                                        MAI 1750
      JOUT = ISW                                                        MAI 1760
      KOUT = ISW                                                        MAI 1770
      NTOC = NTOCA                                                      MAI 1780
      NPAGE = IPAGE(2)                                                  MAI 1790
  310 CONTINUE                                                          MAI 1800
  320 L = L + 1                                                         MAI 1810
      GO TO 50                                                          MAI 1820
C    INP--EXECUTION                                                     MAI 1830
  330 IF (NINP(IFIND(L),2).EQ.0) GO TO 520                              MAI 1840
      IUSW = 1                                                          MAI 1850
      IF (NINP(IFIND(L),2).GT.0) GO TO 340                              MAI 1860
      IUSW = 0                                                          MAI 1870
      NINP(IFIND(L),2) = -NINP(IFIND(L),2)                              MAI 1880
      INPT = 4                                                          MAI 1890
  340 CONTINUE                                                          MAI 1900
C    READ INITIAL CONCENTRATION VECTOR                                  MAI 1910
C    IF NEXT=0  END OF READ                                             MAI 1920
C    IF NEXT=1  LIGHT ELEMENT(ACTIVATION PRODUCTS)                      MAI 1930
C    IF NEXT=2  HEAVY ELEMENT(ACTINIDES+DAUGHTERS)                      MAI 1940
C    IF NEXT=3  FISSION PRODUCT                                         MAI 1950
C    IF NEXT=4,5,6  MULTIPLY BY ISOTOPIC ABUNDANCE                      MAI 1960
C     IF NEXT=7  USE COMPOSITION FOR (ALPHA,N) ADJUSTMENT               MAI 1961
      DO 350 JJ=1,ITOT                                                  MAI 1970
         XNEW(13,JJ) = 0.                                               MAI 1980
  350 CONTINUE                                                          MAI 1990
      DO 360 J=1,10                                                     MAI 2000
         OUTFPB(13,J) = 0.0                                             MAI 2010
  360 CONTINUE                                                          MAI 2020
  370 NEXT = IREAD(1)                                                   MAI 2030
      IF (NEXT.NE.0) GO TO 380                                          MAI 2040
      IF (IUSW.NE.0) GO TO 510                                          MAI 2050
C    READ BURNUP                                                        MAI 2060
      OUTFPB(13,08) = READ(0)                                           MAI 2070
C    READ FLUX                                                          MAI 2080
      OUTFPB(13,09) = READ(0)                                           MAI 2090
C    READ SPECIFIC POWER                                                MAI 2100
      OUTFPB(13,10) = READ(0)                                           MAI 2110
      GO TO 510                                                         MAI 2120
  380 CONTINUE                                                          MAI 2130
      MEXT=NEXT                                                         MAI 2131
      KINP = NINP(IFIND(L),2)                                           MAI 2140
      DO 390 I=1,4                                                      MAI 2150
         INUCL(I) = IREAD(0)                                            MAI 2160
         XCOMP(I) = READ(0)                                             MAI 2170
         IF (INUCL(I).EQ.0) GO TO 400                                   MAI 2180
  390 CONTINUE                                                          MAI 2190
  400 CONTINUE                                                          MAI 2200
      GO TO (410, 420, 430, 410, 420, 430, 410), NEXT                   MAI 2210
  410 ILO = 1                                                           MAI 2220
      IHI = ILITE                                                       MAI 2230
      GO TO 440                                                         MAI 2240
  420 ILO = ILITE + 1                                                   MAI 2250
      IHI = ILITE + IACT                                                MAI 2260
      GO TO 440                                                         MAI 2270
  430 ILO = ILITE + IACT + 1                                            MAI 2280
      IHI = ITOT                                                        MAI 2290
  440 IF (NEXT.GT.3) GO TO 470                                          MAI 2300
      DO 460 I=1,4                                                      MAI 2310
         NUCLI = INUCL(I)                                               MAI 2320
         IF (NUCLI.EQ.0) GO TO 370                                      MAI 2330
         IF (XCOMP(I).LT.ERR) GO TO 460                                 MAI 2340
         DO 450 I1=ILO,IHI                                              MAI 2350
            IF (NUCL(I1).NE.NUCLI) GO TO 450                            MAI 2360
            XNEW(13,I1) = XNEW(13,I1) + XCOMP(I)/RMASS(NUCLI,KINP)      MAI 2370
  450    CONTINUE                                                       MAI 2380
  460 CONTINUE                                                          MAI 2390
      GO TO 370                                                         MAI 2400
  470 DO 500 I=1,4                                                      MAI 2410
         NUCLI = INUCL(I)/10000                                         MAI 2420
         IF (NUCLI.EQ.0) GO TO 370                                      MAI 2430
         RMW = 0.0                                                      MAI 2440
         DO 480 I1=ILO,IHI                                              MAI 2450
            NZ = NUCL(I1)/10000                                         MAI 2460
            IF (NUCLI.NE.NZ) GO TO 480                                  MAI 2470
            RMW = RMW + ABUND(NUCAB(I1))*0.01*RMASS(NUCL(I1),KINP)      MAI 2480
  480    CONTINUE                                                       MAI 2490
         DO 490 I1=ILO,IHI                                              MAI 2500
            NZ = NUCL(I1)/10000                                         MAI 2510
            IF (NUCLI.NE.NZ) GO TO 490                                  MAI 2520
            NA = NUCAB(I1)                                              MAI 2530
            XNEWA = 0.                                                  MAI 2540
            IF (RMW.NE.0.) XNEWA = XCOMP(I)*ABUND(NA)*0.01/RMW          MAI 2550
            XNEW(13,I1) = XNEW(13,I1) + XNEWA                           MAI 2560
  490    CONTINUE                                                       MAI 2570
  500 CONTINUE                                                          MAI 2580
      GO TO 370                                                         MAI 2590
C   SET UP (ALPHA,N) AND SF MULTIPLIERS                                         
  510 CONTINUE                                                                  
      IF(MEXT.NE.7)GO TO 515                                                    
      ANMUL=1.0                                                                 
      ANEXP=1.0                                                                 
      DO 513 I=1,ILITE                                                          
      NZ=NUCL(I)/10000                                                          
      YNEWI=XNEW(13,I)/1000000.0                                                
      IF(NZ.EQ.03)ANMUL=ANMUL+40.0*YNEWI                                        
      IF(NZ.EQ.04)ANMUL=ANMUL+1298.0*YNEWI                                      
      IF(NZ.EQ.05)ANMUL=ANMUL+302.0*YNEWI                                       
      IF(NZ.EQ.06)ANMUL=ANMUL+1.74*YNEWI                                        
      IF(NZ.EQ.09)ANMUL=ANMUL+178.0*YNEWI                                       
      IF(NZ.EQ.11)ANMUL=ANMUL+23.1*YNEWI                                        
      IF(NZ.EQ.12)ANMUL=ANMUL+20.5*YNEWI                                        
      IF(NZ.EQ.01.AND.XNEW(13,I).NE.0.0)                                        
     $ ANEXP=XNEW(13,I)/(ABUND(NUCAB(I))*0.01)                                  
  513 CONTINUE                                                                  
      DO 514 I=1,IACT                                                           
      ALPHAN(I)=ALPHAN(I)*ANMUL*ANMO*ANEXP*ANEO                                 
      SFNU(I)=SFNU(I)*ANEXP*ANEO                                                
  514 CONTINUE                                                                  
      ANMO=ANMO/ANMUL                                                           
      ANEO=ANEO/ANEXP                                                           
C     CALL ANSF(NUCL, Q, ALPHAN,        SFNU, ALPHN, NUCAN, NUCSFU,             
C    *     FFSF, FFA, DIS, ANMUL, ANEXP)                                        
  515 CONTINUE                                                                  
C    SBL - ADDED XP TO ADDMOV LIST                                             
      CALL ADDMOV(2, 13, NINP(IFIND(L),1), ITOT, 1.0, XSTORE, LX,       MAI 2600
     *     XTEMP, XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI, XP)             MAI 2610
      IF (IUSW.LT.1) INPT = IUNIT                                       MAI 2620
  520 IF (NINP(IFIND(L),3).LE.0) GO TO 670                              MAI 2630
      DO 530 JJ=1,ITOT                                                  MAI 2640
         B(JJ) = 0.                                                     MAI 2650
  530 CONTINUE                                                          MAI 2660
C    READ CONTINUOUS FEED PARAMETERS                                    MAI 2670
      ITM = NINP(IFIND(L),5)                                            MAI 2680
      IMS = NINP(IFIND(L),3)                                            MAI 2690
  540 NEXT = IREAD(1)                                                   MAI 2700
      IF (NEXT.EQ.0) GO TO 670                                          MAI 2710
      DO 550 I=1,4                                                      MAI 2720
         INUCL(I) = IREAD(0)                                            MAI 2730
         XCOMP(I) = READ(0)                                             MAI 2740
         IF (INUCL(I).EQ.0) GO TO 560                                   MAI 2750
  550 CONTINUE                                                          MAI 2760
  560 CONTINUE                                                          MAI 2770
      GO TO (570, 580, 590, 570, 580, 590), NEXT                        MAI 2780
  570 ILO = 1                                                           MAI 2790
      IHI = ILITE                                                       MAI 2800
      GO TO 600                                                         MAI 2810
  580 ILO = ILITE + 1                                                   MAI 2820
      IHI = ILITE + IACT                                                MAI 2830
      GO TO 600                                                         MAI 2840
  590 ILO = ILITE + IACT + 1                                            MAI 2850
      IHI = ITOT                                                        MAI 2860
  600 IF (NEXT.GT.3) GO TO 630                                          MAI 2870
      DO 620 I=1,4                                                      MAI 2880
         NUCLI = INUCL(I)                                               MAI 2890
         IF (NUCLI.EQ.0) GO TO 540                                      MAI 2900
         IF (XCOMP(I).LT.ERR) GO TO 620                                 MAI 2910
         DO 610 I1=ILO,IHI                                              MAI 2920
            IF (NUCL(I1).NE.NUCLI) GO TO 610                            MAI 2930
            B(I1) = B(I1) + XCOMP(I)/(RTIME(ITM)*RMASS(NUCLI,IMS))      MAI 2940
  610    CONTINUE                                                       MAI 2950
  620 CONTINUE                                                          MAI 2960
      GO TO 540                                                         MAI 2970
  630 DO 660 I=1,4                                                      MAI 2980
         NUCLI = INUCL(I)/10000                                         MAI 2990
         IF (NUCLI.EQ.0) GO TO 540                                      MAI 3000
         RMW = 0.0                                                      MAI 3010
         DO 640 I1=ILO,IHI                                              MAI 3020
            NZ = NUCL(I1)/10000                                         MAI 3030
            IF (NUCLI.NE.NZ) GO TO 640                                  MAI 3040
            RMW = RMW + ABUND(NUCAB(I1))*0.01*RMASS(NUCL(I1),IMS)       MAI 3050
  640    CONTINUE                                                       MAI 3060
         DO 650 I1=ILO,IHI                                              MAI 3070
            NZ = NUCL(I1)/10000                                         MAI 3080
            IF (NUCLI.NE.NZ) GO TO 650                                  MAI 3090
            NA = NUCAB(I1)                                              MAI 3100
            BI1 = 0.                                                    MAI 3110
            IF (RMW.NE.0.) BI1 = XCOMP(I)*ABUND(NA)*0.01/(RMW*RTIME(ITM)MAI 3120
     *           )                                                      MAI 3130
            B(I1) = B(I1) + BI1                                         MAI 3140
  650    CONTINUE                                                       MAI 3150
  660 CONTINUE                                                          MAI 3160
      GO TO 540                                                         MAI 3170
C    READ CONTINUOUS REPROCESSING PARAMETERS                            MAI 3180
  670 MPROS = NINP(IFIND(L),4)                                          MAI 3190
      IF (MPROS.LT.1) GO TO 710                                         MAI 3200
C    MPROS.GT.0  READ CONT REP                                          MAI 3210
      DO 680 M=1,MPROS                                                  MAI 3220
C    .EQ.0  PROPOGATE PREVIOUS CONT REP, NO READ                        MAI 3230
C    .LT.0  ZERO PREVIOUS CONT REP,NO READ                              MAI 3240
         PRATE(M) = READ(0)                                             MAI 3250
         NOPROS(M) = IREAD(0)                                           MAI 3260
  680 CONTINUE                                                          MAI 3270
      DO 700 M=1,MPROS                                                  MAI 3280
         PRATE(M) = PRATE(M)/RTIME(NINP(IFIND(L),6))                    MAI 3290
         NOPR = NOPROS(M)                                               MAI 3300
         DO 690 N=1,NOPR                                                MAI 3310
            NZPROS(M,N) = IREAD(0)                                      MAI 3320
  690    CONTINUE                                                       MAI 3330
  700 CONTINUE                                                          MAI 3340
  710 L = L + 1                                                         MAI 3350
      GO TO 50                                                          MAI 3360
C    HED--EXECUTION                                                     MAI 3370
  720 IF (NHED(IFIND(L),1).LT.0) GO TO 740                              MAI 3380
      DO 730 J=1,3                                                      MAI 3390
         IOUTTI(NHED(IFIND(L),1),J) = KHED(IFIND(L),J)                  MAI 3400
  730 CONTINUE                                                          MAI 3410
      IO(NHED(IFIND(L),1)) = -1                                         MAI 3420
      GO TO 760                                                         MAI 3430
  740 DO 750 J=1,3                                                      MAI 3440
         ISTOTI(-NHED(IFIND(L),1),J) = KHED(IFIND(L),J)                 MAI 3450
  750 CONTINUE                                                          MAI 3460
      IS(-NHED(IFIND(L),1)) = -1                                        MAI 3470
  760 L = L + 1                                                         MAI 3480
      GO TO 50                                                          MAI 3490
C    REC--EXECUTION                                                     MAI 3500
  770 NREC = NREC + 1                                                   MAI 3510
      L = L + 1                                                         MAI 3520
      GO TO 50                                                          MAI 3530
C    CUT-- EXECUTION                                                    MAI 3540
  780 ERR = RCUT(IFIND(L),NC)                                           MAI 3550
      NTAB2 = NC - 1                                                    MAI 3560
      DO 790 J=1,NTAB2                                                  MAI 3570
         CUTOFF(J) = RCUT(IFIND(L),J)                                   MAI 3580
  790 CONTINUE                                                          MAI 3590
      L = L + 1                                                         MAI 3600
      GO TO 50                                                          MAI 3610
C    KEQ--EXECUTION                                                     MAI 3620
  800 QAPRO = 0.0                                                       MAI 3630
      RIMF = RKEQ(IFIND(L))                                             MAI 3640
      QBPRO = 0.0                                                       MAI 3650
      QCPRO = 0.0                                                       MAI 3660
      QADES = 0.0                                                       MAI 3670
      QBDES = 0.0                                                       MAI 3680
      QCDES = 0.0                                                       MAI 3690
      JECA = NKEQ(IFIND(L),1)                                           MAI 3700
      JECB = NKEQ(IFIND(L),2)                                           MAI 3710
      JECC = NKEQ(IFIND(L),3)                                           MAI 3720
      JECD = NKEQ(IFIND(L),4)                                           MAI 3730
      JECE = NKEQ(IFIND(L),5)                                           MAI 3740
C    PUT ORIGINAL ACT XSEC IN TOCAP AND SAVE PRESENT ACT XSEC           MAI 3750
      DO 810 I=1,150                                                    MAI 3760
         T = TOCAP(I+ILITE)                                             MAI 3770
         TOCAP(I+ILITE) = TOCS(I)                                       MAI 3780
         TOCS(I) = T                                                    MAI 3790
  810 CONTINUE                                                          MAI 3800
      DO 820 I=1,IACT                                                   MAI 3810
         IF (JECA.LT.0) QAPRO = QAPRO + XSTORE(-JECA,I+ILITE)*GENNEU(I)*MAI 3820
     *        0.6023                                                    MAI 3830
         IF (JECB.LT.0) QBPRO = QBPRO + XSTORE(-JECB,I+ILITE)*GENNEU(I)*MAI 3840
     *        0.6023                                                    MAI 3850
         IF (JECC.LT.0) QCPRO = QCPRO + XSTORE(-JECC,I+ILITE)*GENNEU(I)*MAI 3860
     *        0.6023                                                    MAI 3870
         IF (JECA.GT.0) QAPRO = QAPRO + XNEW(JECA,I+ILITE)*GENNEU(I)*   MAI 3880
     *        0.6023                                                    MAI 3890
         IF (JECB.GT.0) QBPRO = QBPRO + XNEW(JECB,I+ILITE)*GENNEU(I)*   MAI 3900
     *        0.6023                                                    MAI 3910
         IF (JECC.GT.0) QCPRO = QCPRO + XNEW(JECC,I+ILITE)*GENNEU(I)*   MAI 3920
     *        0.6023                                                    MAI 3930
  820 CONTINUE                                                          MAI 3940
      DO 830 I=1,ITOT                                                   MAI 3950
         IF (JECA.LT.0) QADES = QADES + XSTORE(-JECA,I)*TOCAP(I)*0.6023 MAI 3960
         IF (JECB.LT.0) QBDES = QBDES + XSTORE(-JECB,I)*TOCAP(I)*0.6023 MAI 3970
         IF (JECC.LT.0) QCDES = QCDES + XSTORE(-JECC,I)*TOCAP(I)*0.6023 MAI 3980
         IF (JECA.GT.0) QADES = QADES + XNEW(JECA,I)*TOCAP(I)*0.6023    MAI 3990
         IF (JECB.GT.0) QBDES = QBDES + XNEW(JECB,I)*TOCAP(I)*0.6023    MAI 4000
         IF (JECC.GT.0) QCDES = QCDES + XNEW(JECC,I)*TOCAP(I)*0.6023    MAI 4010
  830 CONTINUE                                                          MAI 4020
C    PUT PRESENT ACT XSEC IN TOCAP AND SAVE ORIGINAL ACT XSEC IN TOCS   MAI 4030
      DO 840 I=1,150                                                    MAI 4040
         T = TOCAP(I+ILITE)                                             MAI 4050
         TOCAP(I+ILITE) = TOCS(I)                                       MAI 4060
         TOCS(I) = T                                                    MAI 4070
  840 CONTINUE                                                          MAI 4080
      AIMF = RKEQ(IFIND(L))                                             MAI 4090
      IF (AIMF.LE.0.0 .AND. QADES.GT.ERR) AIMF = QAPRO/QADES            MAI 4100
      IF (RKEQ(IFIND(L)).EQ.0.0) AIMF = RMULV(NREC,1)                   MAI 4110
      IF (QBDES.GT.ERR) BIMF = QBPRO/QBDES                              MAI 4120
      IF (QCDES.GT.ERR) CIMF = QCPRO/QCDES                              MAI 4130
      IF (ABS(AIMF-CIMF).GT.ERR) FRC = (BIMF-AIMF)/(AIMF-CIMF)          MAI 4140
      IF (FRC.GE.0.0) GO TO 850                                         MAI 4150
      WRITE (06,99998) FRC, AIMF, JECA, BIMF, JECB, CIMF, JECC          MAI 4160
      FRC = 0.0                                                         MAI 4170
  850 CONTINUE                                                          MAI 4180
      IF (ABS(QCDES).GT.0.0) FRD = FRC*QBDES/QCDES                      MAI 4190
      GRD = 1.0 - FRD                                                   MAI 4200
      IF (GRD.LT.0.0) GRD = 0.0                                         MAI 4210
      WRITE (15,99997) QAPRO, QADES, AIMF, QBPRO, QBDES, BIMF, QCPRO,   MAI 4220
     *     QCDES, CIMF, FRC, FRD                                        MAI 4230
C    SBL - ADDED XP TO ADDMOV LIST                                             
      CALL ADDMOV(2, JECB, JECD, ITOT, 1.0, XSTORE, LX, XTEMP, XNEW,    MAI 4240
     *     MX, STTFPB, ISTOTI, IS, RSTOTI,XP)                           MAI 4250
      CALL ADDMOV(1, JECC, JECD, ITOT, FRD, XSTORE, LX, XTEMP, XNEW,    MAI 4260
     *     MX, STTFPB, ISTOTI, IS, RSTOTI,XP)                           MAI 4270
      CALL ADDMOV(2, JECC, JECE, ITOT, GRD, XSTORE, LX, XTEMP, XNEW,    MAI 4280
     *     MX, STTFPB, ISTOTI, IS, RSTOTI,XP)                           MAI 4290
      L = L + 1                                                         MAI 4300
      GO TO 50                                                          MAI 4310
C    D0L--EXECUTION                                                     MAI 4320
  860 IF (NDOL(IFIND(L),3).GT.0) GO TO 890                              MAI 4330
      DO 870 LL=1,300                                                   MAI 4340
         IF (NDOL(IFIND(L),1).EQ.LAB(LL)) GO TO 880                     MAI 4350
  870 CONTINUE                                                          MAI 4360
      WRITE (IOUT,99985) IFIND(L)                                       MAI 4370
  880 NDOL(IFIND(L),4) = LL                                             MAI 4380
      NCON(IFIND(LL),1) = L                                             MAI 4390
  890 L = L + 1                                                         MAI 4400
      GO TO 50                                                          MAI 4410
C    MOV--EXECUTION                                                     MAI 4420
  900 IN = IFIND(L)                                                     MAI 4430
      RMULT = RMOV(IN)                                                  MAI 4440
      IF (NMOV(IN,3).LT.0) RMULT = RMULT*FAC(-NMOV(IN,3))               MAI 4450
      IF (NMOV(IN,3).GT.0) RMULT = RMULV(NREC,NMOV(IN,3))*RMOV(IN)      MAI 4460
      CALL ADDMOV(2, NMOV(IN,1), NMOV(IN,2), ITOT, RMULT, XSTORE, LX,   MAI 4470
     *     XTEMP, XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI, XP)             MAI 4480
      L = L + 1                                                         MAI 4490
      GO TO 50                                                          MAI 4500
C    ADD--EXECUTION                                                     MAI 4510
  910 IN = IFIND(L)                                                     MAI 4520
      RMULT = RADD(IN)                                                  MAI 4530
      IF (NADD(IN,3).LT.0) RMULT = RMULT*FAC(-NADD(IN,3))               MAI 4540
      IF (NADD(IN,3).GT.0) RMULT = RMULV(NREC,NADD(IN,3))*RADD(IN)      MAI 4550
      CALL ADDMOV(1, NADD(IN,1), NADD(IN,2), ITOT, RMULT, XSTORE, LX,   MAI 4570
     *     XTEMP, XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI, XP)             MAI 4580
      L = L + 1                                                         MAI 4590
      GO TO 50                                                          MAI 4600
C    BUP--EXECUTION                                                     MAI 4610
  920 BU = -BU                                                          MAI 4620
      IF (BU.LT.0.0) GO TO 1000                                         MAI 4630
      OUTFPB(13,4) = POWER*DELT/86400.0                                 MAI 4640
      TT = 0.0                                                          MAI 4650
      IF (MIRR.EQ.0) GO TO 960                                          MAI 4660
      DO 930 MV=1,MIRR                                                  MAI 4670
         IF (FPBST(MV,2).GT.ERR) TT = TT + FPBST(MV,1)                  MAI 4680
  930 CONTINUE                                                          MAI 4690
      DO 950 MV=1,MIRR                                                  MAI 4700
         DO 940 J=2,3                                                   MAI 4710
            IF (TT.GT.ERR) AA(J) = AA(J) + FPBST(MV,J)*FPBST(MV,1)/TT   MAI 4720
  940    CONTINUE                                                       MAI 4730
  950 CONTINUE                                                          MAI 4740
  960 CONTINUE                                                          MAI 4750
      AA(1) = AA(3)*TT/86400.0                                          MAI 4760
      IF (MIRR.EQ.0) GO TO 990                                          MAI 4770
      DO 980 MV=1,MIRR                                                  MAI 4780
         DO 970 J=1,3                                                   MAI 4790
            NVE = NVECST(MV)                                            MAI 4800
            IF (NVE.LT.0) STTFPB(-NVE,J+7) = AA(J)                      MAI 4810
            IF (NVE.GT.0) OUTFPB(NVE,J+7) = AA(J)                       MAI 4820
  970    CONTINUE                                                       MAI 4830
  980 CONTINUE                                                          MAI 4840
  990 CONTINUE                                                          MAI 4850
 1000 L = L + 1                                                         MAI 4860
      MIRR = 0                                                          MAI 4870
      IF (BU.GT.0.0) GO TO 50                                           MAI 4880
      DO 1010 J=8,10                                                    MAI 4890
         AA(J-7) = 0.0                                                  MAI 4900
 1010 CONTINUE                                                          MAI 4910
      GO TO 50                                                          MAI 4920
C    PCH--EXECUTION                                                     MAI 4930
 1020 NEXT = 0                                                          MAI 4940
      DO 1080 J=1,3                                                     MAI 4950
         MPCH = NPCH(IFIND(L),J)                                        MAI 4960
         IF (IDM(J+1).LT.1) GO TO 1080                                  MAI 4970
         IF (MPCH.EQ.0) GO TO 1080                                      MAI 4980
         I1 = 0                                                         MAI 4990
         DO 1030 JG=1,J                                                 MAI 5000
            I1 = I1 + IDM(JG)                                           MAI 5010
 1030    CONTINUE                                                       MAI 5020
         I2 = I1 + IDM(J+1)                                             MAI 5030
         I1 = I1 + 1                                                    MAI 5040
         KPCH = IABS(MPCH)                                              MAI 5050
         K = 1                                                          MAI 5060
         DO 1060 I=I1,I2                                                MAI 5070
            INUCL(K) = NUCL(I)                                          MAI 5080
            IF (MPCH.GT.0) XCOMP(K) = XNEW(KPCH,I)                      MAI 5090
            IF (MPCH.LT.0) XCOMP(K) = XSTORE(KPCH,I)                    MAI 5100
            IF (XCOMP(K).LT.ERR) GO TO 1060                             MAI 5110
            IF (K.NE.4) GO TO 1050                                      MAI 5120
            WRITE (7,99983) J, (INUCL(KJ),XCOMP(KJ),KJ=1,4)             MAI 5130
            DO 1040 KK=1,4                                              MAI 5140
               INUCL(KK) = 0                                            MAI 5150
               XCOMP(KK) = 0.0                                          MAI 5160
 1040       CONTINUE                                                    MAI 5170
            K = 0                                                       MAI 5180
 1050       K = K + 1                                                   MAI 5190
 1060    CONTINUE                                                       MAI 5200
         IF (K.EQ.1) GO TO 1080                                         MAI 5210
         DO 1070 KK=K,4                                                 MAI 5220
            INUCL(KK) = 0                                               MAI 5230
            XCOMP(KK) = 0.0                                             MAI 5240
 1070    CONTINUE                                                       MAI 5250
         WRITE (7,99983) J, (INUCL(KJ),XCOMP(KJ),KJ=1,4)                MAI 5260
 1080 CONTINUE                                                          MAI 5270
      DO 1090 J=1,3                                                     MAI 5280
         IF (MPCH.GT.0) AA(J) = OUTFPB(KPCH,J+7)                        MAI 5290
         IF (MPCH.LT.0) AA(J) = STTFPB(KPCH,J+7)                        MAI 5300
 1090 CONTINUE                                                          MAI 5310
      WRITE (7,99996) NEXT, (AA(J),J=1,3)                               MAI 5320
C    NEXT POWER BURNUP FLUX                                             MAI 5330
      L = L + 1                                                         MAI 5340
      GO TO 50                                                          MAI 5350
C    LIP--EXECUTION                                                     MAI 5360
 1100 DO 1110 J=1,3                                                     MAI 5370
         MLIP(J) = NLIP(IFIND(L),J)                                     MAI 5380
 1110 CONTINUE                                                          MAI 5390
      L = L + 1                                                         MAI 5400
      GO TO 50                                                          MAI 5410
C    WAC--EXECUTION                                                     MAI 5420
C    SBL ADDED XP TO ADDMOV LIST                                               
 1120 CALL ADDMOV(2, NWAC(IFIND(L),2), 0, ITOT, 1.0, XSTORE, LX, XTEMP, MAI 5430
     *     XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI,XP)                     MAI 5440
      NW = NWAC(IFIND(L),1)                                             MAI 5450
      DO 1130 I=1,ITOT                                                  MAI 5460
         NELEM = (NUCL(I)-MOD(NUCL(I),10000))/10000                     MAI 5470
         IF (NW.GT.0) B(I) = XTEMP(I)*FREC(NW,NELEM)                    MAI 5480
         IF (NW.LT.0) B(I) = XTEMP(I)*FRECG(-NW,NGROUP(NELEM))          MAI 5490
 1130 CONTINUE                                                          MAI 5500
      L = L + 1                                                         MAI 5510
      GO TO 50                                                          MAI 5520
C    LIB--EXECUTION                                                     MAI 5530
 1140 DO 1150 J=1,11                                                    MAI 5540
         MLIB(J) = NLIB(IFIND(L),J)                                     MAI 5550
 1150 CONTINUE                                                          MAI 5560
      NNSORS = MYLD(MLIB(11))                                           MAI 5570
      IFIRST = 0                                                        MAI 5580
      IJKO = 0                                                          MAI 5590
      FPRATN = 1.0                                                      MAI 5600
      IIUT = IOUT                                                       MAI 5610
      OLDFLX = 0.0                                                      MAI 5620
      CALL NUDAT1(NUCL, Q, FG, GENNEU, AMPC, WMPC, TOCAP, FISS, DIS,    MAI 5630
     *     YIELD, COEFF, LC, NPROD, FFSF, FFA, NY, YY, NUCAB, NYIELD,   MAI 5640
     *     LOC, NON0, KD, MMAX, KAP, AP, ABUND, NNSORS, LPUCNT)         MAI 5650
      CALL NUDAT2(NUCL, Q, FG, GENNEU, AMPC, WMPC, TOCAP, FISS, DIS,    MAI 5660
     *     YIELD, COEFF, LC, NPROD, FFSF, FFA, NY, YY, NUCAB, NYIELD,   MAI 5670
     *     LOC, NON0, KD, MMAX, KAP, AP, ABUND, NNSORS, LPUCNT, ICNMAX) MAI 5680
      CALL NUDAT3(NUCL, Q, FG, GENNEU, AMPC, WMPC, TOCAP, A, FISS, DIS, MAI 5690
     *     YIELD, COEFF, LC, NPROD, FFSF, FFA, NY, YY, NUCAB, NYIELD,   MAI 5700
     *     LOC, NON0, KD, MMAX, KAP, AP, ABUND, NNSORS, NFUDFP, IFD)    MAI 5710
      IF (NSIZE(13).GT.IZMAX) WRITE (15,99995) NSIZE(13), IZMAX         MAI 5720
      IF (ILITE.GT.ILMAX) WRITE (15,99994) ILITE, ILMAX                 MAI 5730
      IF (IACT.GT.IAMAX) WRITE (15,99993) IACT, IAMAX                   MAI 5740
      IF (IFP.GT.IFMAX) WRITE (15,99992) IFP, IFMAX                     MAI 5750
      IF (ITOT.GT.ITMAX) WRITE (15,99991) ITOT, ITMAX                   MAI 5760
      IF (NSIZE(19).GT.NABMAX) WRITE (15,99990) NSIZE(19), NABMAX       MAI 5770
      IF (NSIZE(14).GT.IFYMAX) WRITE (15,99989) NSIZE(14), IFYMAX       MAI 5780
      NUDUN = MLIB(8)                                                   MAI 5790
      IF (NUDUN.NE.IUNIT) REWIND NUDUN                                  MAI 5800
      CALL ANSF(NUCL, Q, ALPHAN,        SFNU, ALPHN, NUCAN, NUCSFU,     MAI 5810
     *     FFSF, FFA, DIS, ANMUL, ANEXP)                                MAI 5820
      DO 1160 I=2,ITOT                                                  MAI 5830
         NON0(I) = NON0(I) + NON0(I-1)                                  MAI 5840
         KD(I) = KD(I) + NON0(I-1)                                      MAI 5850
 1160 CONTINUE                                                          MAI 5860
      DO 1170 I=1,150                                                   MAI 5870
         TOCS(I) = TOCAP(I+ILITE)                                       MAI 5880
 1170 CONTINUE                                                          MAI 5890
      QXN = 0.001                                                       MAI 5900
      AXN = -ALOG(QXN)                                                  MAI 5910
      L = L + 1                                                         MAI 5920
      GO TO 50                                                          MAI 5930
C    PHO--EXECUTION                                                     MAI 5940
 1180 DO 1190 J=1,4                                                     MAI 5950
         MPHO(J) = NPHO(IFIND(L),J)                                     MAI 5960
 1190 CONTINUE                                                          MAI 5970
      CALL PHOLIB(NUCL, NGF, NGN, NGR, GGR, AP(1), AP(20))              MAI 5980
      IF (NSIZE(20).GT.IPHMAX) WRITE (15,99988) NSIZE(20), IPHMAX       MAI 5990
      NPHUN = MPHO(4)                                                   MAI 6000
      IF (NPHUN.NE.IUNIT) REWIND NPHUN                                  MAI 6010
      INPT = IUNIT                                                      MAI 6020
      L = L + 1                                                         MAI 6030
      GO TO 50                                                          MAI 6040
C    LPU--EXECUTION                                                     MAI 6050
 1200 L = L + 1                                                         MAI 6060
      GO TO 50                                                          MAI 6070
C    IRF,IRP,AND DEC--EXECUTION                                         MAI 6080
 1210 INDEX = 1                                                         MAI 6090
      NEGF = 1                                                          MAI 6100
      FLUX = RIRR(IFIND(L),2)                                           MAI 6110
      IF (FLUX.GT.0.0) FLUX = FLUX*RPF                                  MAI 6120
      N4 = NIRR(IFIND(L),4)                                             MAI 6130
      IF (N4.EQ.1 .OR. N4.EQ.2 .OR. N4.EQ.4) MIRR = 0                   MAI 6140
      IF (FLUX.GT.0.0) GO TO 1220                                       MAI 6150
      FLUX = FLSTOR(MIRR+1)*(-FLUX)*RPF                                 MAI 6160
      NEGF = -1                                                         MAI 6170
 1220 CONTINUE                                                          MAI 6180
      POWER = 0.                                                        MAI 6190
      GO TO 1260                                                        MAI 6200
 1230 INDEX = 0                                                         MAI 6210
      POWER = RIRR(IFIND(L),2)                                          MAI 6220
      IF (POWER.GT.0.0) POWER = POWER*RPF                               MAI 6230
      IF (POWER.GT.0.0) GO TO 1240                                      MAI 6240
      POWER = POSTOR(MIRR+1)*(-POWER)*RPF                               MAI 6250
      NEGF = -1                                                         MAI 6260
 1240 CONTINUE                                                          MAI 6270
      FLUX = 0.                                                         MAI 6280
      GO TO 1260                                                        MAI 6290
 1250 INDEX = 2                                                         MAI 6300
      POWER = 0.                                                        MAI 6310
      FLUX = 0.                                                         MAI 6320
 1260 CONTINUE                                                          MAI 6330
      IF (FLUX.EQ.0.0 .AND. POWER.EQ.0.0) INDEX = 2                     MAI 6340
      N4 = NIRR(IFIND(L),4)                                             MAI 6350
      IF (N4.EQ.0) GO TO 1290                                           MAI 6360
      IF (N4.NE.1) NSIZE(21) = 0                                        MAI 6370
      IF (N4.EQ.1 .OR. N4.EQ.2 .OR. N4.EQ.4) TMO = 0.0                  MAI 6380
      IF (N4.EQ.1 .OR. N4.EQ.2 .OR. N4.EQ.4) MIRR = 0                   MAI 6390
      IF (N4.NE.2 .AND. N4.NE.3 .AND. N4.NE.4) GO TO 1290               MAI 6400
      IF(MIRR.EQ.0) XOFF=TOFF                                                   
      KKK = 7                                                           MAI 6410
      IF (N4.EQ.2) KKK = 10                                             MAI 6420
      DO 1280 K=1,13                                                    MAI 6430
         DO 1270 KK=1,KKK                                               MAI 6440
            OUTFPB(K,KK) = 0.0                                          MAI 6450
 1270    CONTINUE                                                       MAI 6460
 1280 CONTINUE                                                          MAI 6470
C    SBL   ADDED XP TO ADDMOV LIST                                             
 1290 IF (NIRR(IFIND(L),1).NE.0) CALL ADDMOV(2, NIRR(IFIND(L),1), 0,    MAI 6480
     *     ITOT, 1.0, XSTORE, LX, XTEMP, XNEW, MX, STTFPB, ISTOTI, IS,  MAI 6490
     *     RSTOTI,XP)                                                   MAI 6500
      OUTFPB(13,6) = 0.0                                                MAI 6510
      OUTFPB(13,7) = 0.0                                                MAI 6520
      TSEC = RIRR(IFIND(L),1)*RTIME(NIRR(IFIND(L),3))                   MAI 6530
      OUTFPB(13,1) = TSEC                                               MAI 6540
      DELT=TSEC-TMO-XOFF                                                        
      IF(DELT.GT.0.0)XOFF=0.0                                                   
      IF(DELT.GT.0.0)GO TO 1292                                                 
C     SBL ADDED XP TO ADDMOV LIST                                              
      CALL ADDMOV(2,NIRR(IFIND(L),1),NIRR(IFIND(L),2),ITOT,1.0,XSTORE,          
     * LX,XTEMP,XNEW,MX,STTFPB,ISTOTI,IS,RSTOTI,XP)                            
      NB=NIRR(IFIND(L),2)                                                       
      IF(NB.EQ.0)GO TO 1291                                                     
      IF(NB.LT.0)GO TO 1293                                                     
      IO(NB)=-1                                                                 
      IOUTTI(NB,1)=INFE(1)                                                      
      IOUTTI(NB,2)=INFE(2)                                                      
      IOUTTI(NB,3)=INFE(3)                                                      
      GO TO 1291                                                                
 1293 IS(-NB)=-1                                                                
      ISTOTI(-NB,1)=INFE(1)                                                     
      ISTOTI(-NB,2)=INFE(2)                                                     
      ISTOTI(-NB,3)=INFE(3)                                                     
 1291 CONTINUE                                                                  
      GO TO 1340                                                                
 1292 CONTINUE                                                                  
      TMN = TSEC                                                        MAI 6560
      TMB = DELT                                                        MAI 6570
      VECT(MIRR+1) = 0.0                                                MAI 6580
      IF (INDEX.NE.2 .AND. NXSEC.GT.0) CALL FUDGE(NFUDFP,               MAI 6590
     *     NIRR(IFIND(L),1), IFD, XSTORE, LX, XNEW, MX, TOCAP, A, FISS, MAI 6600
     *     MIRR, DELT)                                                  MAI 6610
      MIRR = MIRR + 1                                                   MAI 6620
      M = MIRR                                                          MAI 6630
      CALL FLUXO(DELT, OLDFLX, M, TSEC, A, LOC, B, D, NON0, KD, LONG,   MAI 6640
     *     NUCL, TOCAP, XPAR, XTEMP, AP, IOUT, FISS, DIS, NFUDFP, IFD,  MAI 6650
     *     MIRR, FPRATN, IJKO, IFIRST)                                  MAI 6660
      CALL DECAY(DELT, A, D, LOC, NON0, KD, LONG, XP, XPAR, XTEMP, NQ,  MAI 6670
     *     XNEW, MX, IOUT, NUCL)                                        MAI 6680
      CALL TERM(TMB, A, B, D, LOC, NON0, KD, NONP, LOCP, NQ, LONG,      MAI 6690
     *     CIMN, CSUM, XPAR, XP, XTEMP, AP, XNEW, MX, IOUT, LOCA,       MAI 6700
     *     IAPMAX, NUCL)                                                MAI 6710
      CALL EQUIL(A, LOC, B, D, NON0, KD, LONG, XP, XPAR, XTEMP, XNEW,   MAI 6720
     *     MX, IOUT, NUCL)                                              MAI 6730
      IFIRST = 1                                                        MAI 6740
      OUTFPB(13,2) = FLUX                                               MAI 6750
      IF (NEGF.LT.0 .OR. NJUMP.NE.0) GO TO 1300                         MAI 6760
      FLSTOR(MIRR) = FLUX                                               MAI 6770
      POSTOR(MIRR) = POWER                                              MAI 6780
 1300 CONTINUE                                                          MAI 6790
      OUTFPB(13,3) = POWER                                              MAI 6800
      OUTFPB(13,4) = POWER*DELT/86400.0                                 MAI 6810
      FPBST(MIRR,1) = DELT                                              MAI 6820
      FPBST(MIRR,2) = FLUX                                              MAI 6830
      FPBST(MIRR,3) = POWER                                             MAI 6840
      DO 1310 I=1,IACT                                                  MAI 6850
         OUTFPB(13,6) = OUTFPB(13,6) + XNEW(13,I+ILITE)*GENNEU(I)*0.6023MAI 6860
 1310 CONTINUE                                                          MAI 6870
      DO 1320 I=1,ITOT                                                  MAI 6880
         OUTFPB(13,7) = OUTFPB(13,7) + XNEW(13,I)*TOCAP(I)*0.6023       MAI 6890
 1320 CONTINUE                                                          MAI 6900
      IF (OUTFPB(13,7).GT.ERR) OUTFPB(13,5) = OUTFPB(13,6)/OUTFPB(13,7) MAI 6910
      TMO = TMN                                                         MAI 6920
C    SBL ADDED XP TO ADDMOV LIST                                               
      CALL ADDMOV(2, 13, NIRR(IFIND(L),2), ITOT, 1.0, XSTORE, LX,       MAI 6930
     *     XTEMP, XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI, XP)             MAI 6940
      NTEST = NIRR(IFIND(L),2)                                          MAI 6950
      NVECST(MIRR) = NTEST                                              MAI 6960
      IF (NTEST.EQ.0) GO TO 1340                                        MAI 6970
      IF (NTEST.LT.0) GO TO 1330                                        MAI 6980
      ROUTTI(NTEST) = RIRR(IFIND(L),1)                                  MAI 6990
      IO(NTEST) = NIRR(IFIND(L),3)                                      MAI 7000
      GO TO 1340                                                        MAI 7010
 1330 RSTOTI(-NTEST) = RIRR(IFIND(L),1)                                 MAI 7020
      IS(-NTEST) = NIRR(IFIND(L),3)                                     MAI 7030
 1340 L = L + 1                                                         MAI 7040
      GO TO 50                                                          MAI 7050
C    PRO--EXECUTION                                                     MAI 7060
C    SBL ADDED XP TO ADDMOV LIST                                               
 1350 IF (NPRO(IFIND(L),1).NE.0) CALL ADDMOV(2, NPRO(IFIND(L),1), 0,    MAI 7070
     *     ITOT, 1.0, XSTORE, LX, XTEMP, XNEW, MX, STTFPB, ISTOTI, IS,  MAI 7080
     *     RSTOTI,XP)                                                   MAI 7090
      IF (NPRO(IFIND(L),4).LT.0) GO TO 1470                             MAI 7100
C    THIS SECTION IS USED IF INDIVIDUAL FRACTIONAL RECOVERIES ARE USED  MAI 7110
      DO 1360 I=1,ITOT                                                  MAI 7120
         NELEM = (NUCL(I)-MOD(NUCL(I),10000))/10000                     MAI 7130
C    FREC(SET #, ELEMENT #)                                             MAI 7140
         XNEW(13,I) = XTEMP(I)*FREC(NPRO(IFIND(L),4),NELEM)             MAI 7150
 1360 CONTINUE                                                          MAI 7160
C    SBL ADDED XP TO ADDMOV LIST                                               
 1370 CALL ADDMOV(2, 13, NPRO(IFIND(L),2), ITOT, 1.0, XSTORE, LX,       MAI 7170
     *     XTEMP, XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI, XP)             MAI 7180
      DO 1380 I=1,ITOT                                                  MAI 7190
         XNEW(13,I) = XTEMP(I) - XNEW(13,I)                             MAI 7200
 1380 CONTINUE                                                          MAI 7210
C   SBL ADDED XP TO ADDMOV LIST                                                
      CALL ADDMOV(2, 13, NPRO(IFIND(L),3), ITOT, 1.0, XSTORE, LX,       MAI 7220
     *     XTEMP, XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI, XP)             MAI 7230
      DO 1460 J=2,3                                                     MAI 7240
         IF (NPRO(IFIND(L),J).LT.0) GO TO 1420                          MAI 7250
         IF (NPRO(IFIND(L),1).LT.0) GO TO 1400                          MAI 7260
         DO 1390 JJ=8,10                                                MAI 7270
            OUTFPB(NPRO(IFIND(L),J),JJ) = OUTFPB(NPRO(IFIND(L),1),JJ)   MAI 7280
 1390    CONTINUE                                                       MAI 7290
         GO TO 1460                                                     MAI 7300
 1400    DO 1410 JJ=8,10                                                MAI 7310
            OUTFPB(NPRO(IFIND(L),J),JJ) = STTFPB(-NPRO(IFIND(L),1),JJ)  MAI 7320
 1410    CONTINUE                                                       MAI 7330
         GO TO 1460                                                     MAI 7340
 1420    IF (NPRO(IFIND(L),1).LT.0) GO TO 1440                          MAI 7350
         DO 1430 JJ=8,10                                                MAI 7360
            STTFPB(-NPRO(IFIND(L),J),JJ) = OUTFPB(NPRO(IFIND(L),1),JJ)  MAI 7370
 1430    CONTINUE                                                       MAI 7380
         GO TO 1460                                                     MAI 7390
 1440    DO 1450 JJ=8,10                                                MAI 7400
            STTFPB(-NPRO(IFIND(L),J),JJ) = STTFPB(-NPRO(IFIND(L),1),JJ) MAI 7410
 1450    CONTINUE                                                       MAI 7420
 1460 CONTINUE                                                          MAI 7430
      GO TO 1490                                                        MAI 7440
C    THIS SECTION IS USED IF GROUP FRACTIONAL RECOVERIES ARE USED       MAI 7450
 1470 DO 1480 I=1,ITOT                                                  MAI 7460
         NELEM = (NUCL(I)-MOD(NUCL(I),10000))/10000                     MAI 7470
C    FRECG(SET #,GROUP #)                                               MAI 7480
         XNEW(13,I) = XTEMP(I)*FRECG(-NPRO(IFIND(L),4),NGROUP(NELEM))   MAI 7490
 1480 CONTINUE                                                          MAI 7500
      GO TO 1370                                                        MAI 7510
 1490 L = L + 1                                                         MAI 7520
      GO TO 50                                                          MAI 7530
C    OPT--EXECUTION                                                     MAI 7540
 1500 N1 = 1                                                            MAI 7550
      N2 = NTAB                                                         MAI 7560
      GO TO 1530                                                        MAI 7570
 1510 N1 = NTAB + 1                                                     MAI 7580
      N2 = 2*NTAB                                                       MAI 7590
      GO TO 1530                                                        MAI 7600
 1520 N1 = 2*NTAB + 1                                                   MAI 7610
      N2 = 3*NTAB                                                       MAI 7620
 1530 DO 1540 J=N1,N2                                                   MAI 7630
         NTO(J) = NOPT(IFIND(L),J)                                      MAI 7640
 1540 CONTINUE                                                          MAI 7650
      L = L + 1                                                         MAI 7660
      GO TO 50                                                          MAI 7670
C    CON--EXECUTION                                                     MAI 7680
 1550 CONTINUE                                                          MAI 7690
      IF (NJUMP.EQ.0) GO TO 1570                                        MAI 7700
      IF (L.EQ.NCEND) GO TO 1560                                        MAI 7710
      L = L + 1                                                         MAI 7720
      GO TO 50                                                          MAI 7730
 1560 L = NJUMP + 1                                                     MAI 7740
      NJUMP = 0                                                         MAI 7750
      RPF = 1.0                                                         MAI 7760
      GO TO 50                                                          MAI 7770
 1570 CONTINUE                                                          MAI 7780
      NF = NCON(IFIND(L),1)                                             MAI 7790
      NFF = IFIND(NF)                                                   MAI 7800
      L = L + 1                                                         MAI 7810
      IF (NINSTR(NF).NE.11 .OR. NDOL(NFF,4).NE.(L-1)) GO TO 50          MAI 7820
      NDOL(NFF,3) = NDOL(NFF,3) + 1                                     MAI 7830
      IF (NDOL(NFF,3).LT.NDOL(NFF,2)) GO TO 1580                        MAI 7840
      NDOL(NFF,3) = 0                                                   MAI 7850
      GO TO 50                                                          MAI 7860
 1580 L = NCON(IFIND(L-1),1)                                            MAI 7870
      GO TO 50                                                          MAI 7880
C    STP--EXECUTION                                                     MAI 7890
 1590 L = 1                                                             MAI 7900
      RETURN                                                            MAI 7910
C    GTO--EXECUTION                                                     MAI 7920
 1600 CONTINUE                                                          MAI 7930
      IF (NGTO(IFIND(L),1).LT.0) GO TO 1650                             MAI 7940
      NCBEG = 0                                                         MAI 7950
      NCEND = 0                                                         MAI 7960
      DO 1610 J=1,300                                                   MAI 7970
         IF (NINSTR(J).EQ.29 .OR. NINSTR(J).EQ.30) GO TO 1620           MAI 7980
         IF (NINSTR(J).NE.28) GO TO 1610                                MAI 7990
         IF (NCON(IFIND(J),1).EQ.NGTO(IFIND(L),1)) NCBEG = J            MAI 8000
         IF (NCON(IFIND(J),1).EQ.NGTO(IFIND(L),2)) NCEND = J            MAI 8010
         IF (NCBEG.NE.0 .AND. NCEND.NE.0) GO TO 1630                    MAI 8020
 1610 CONTINUE                                                          MAI 8030
 1620 WRITE (15,99987) (NGTO(IF,K),K=1,2), RGTO(IF), NCBEG, NCEND       MAI 8040
      GO TO 1640                                                        MAI 8050
 1630 NJUMP = L                                                         MAI 8060
      RPF = RGTO(IFIND(L))                                              MAI 8070
      L = NCBEG                                                         MAI 8080
      GO TO 50                                                          MAI 8090
 1640 L = L + 1                                                         MAI 8100
      GO TO 50                                                          MAI 8110
C    TOFF--EXECUTION                                                            
 1641 IF=IFIND(L)                                                               
      TOFF=ROFF(IF )*RTIME(NOFF(IF ))                                           
      L=L+1                                                                     
      GO TO 50                                                                  
C    END--EXECUTION                                                     MAI 8120
 1650 NSTP = 4                                                          MAI 8130
      RETURN                                                            MAI 8140
99999 FORMAT (10H **FAC LD=, I3, 2X, 8HFAC(LD)=, 1PE11.4)               MAI 8150
99998 FORMAT (19H0FRACTION NEGATIVE=, 1PE10.3, 10X, 3(4HIMF=, 1PE10.3,  MAI 8160
     *     12H FOR VECTOR=, I3, 1X))                                    MAI 8170
99997 FORMAT (6H **KEQ, 1X, 6HNPROA=, 1PE10.3, 1X, 6HNDESA=, 1PE10.3,   MAI 8180
     *     1X, 5HIMFA=, 1PE10.3, 1X, 6HNPROB=, 1PE10.3, 1X, 6HNDESB=,   MAI 8190
     *     1PE10.3, 1X, 5HIMFB=, 1PE10.3, /, 7H NPROC=, 1PE10.3, 1X,    MAI 8200
     *     6HNDESC=, 1PE10.3, 1X, 5HIMFC=, 1PE10.3, 1X, 12HFRC=FRD BEFO,MAI 8210
     *     18HRE N DEST SCALING=, 1PE10.3, 1X, 21HFRD=FRAC OF VECT C IN,MAI 8220
     *     7HCLUDED=, 1PE10.3)                                          MAI 8230
99996 FORMAT (2X, I2, 6X, 1PE12.5, 2X, 1PE12.5, 2X, 1PE9.2)             MAI 8240
99995 FORMAT (45H **WARNING 9600** NUMBER OF MATRIX ELEMENTS =, I5,     MAI 8250
     *     28H IS GREATER THAN ALLOWABLE =, I5)                         MAI 8260
99994 FORMAT (49H **WARNING 9601** NUMBER OF ACTIVATION PRODUCTS =, I5, MAI 8270
     *     28H IS GREATER THAN ALLOWABLE =, I5)                         MAI 8280
99993 FORMAT (39H **WARNING 9602** NUMBER OF ACTINIDES =, I5, 7H IS GRE,MAI 8290
     *     21HATER THAN ALLOWABLE =, I5)                                MAI 8300
99992 FORMAT (46H **WARNING 9603** NUMBER OF FISSION PRODUCTS =, I5,    MAI 8310
     *     28H IS GREATER THAN ALLOWABLE =, I5)                         MAI 8320
99991 FORMAT (44H **WARNING 9604** TOTAL NUMBER OF NUCLIDES =, I5,      MAI 8330
     *     28H IS GREATER THAN ALLOWABLE =, I5)                         MAI 8340
99990 FORMAT (48H **WARNING 9607** NUMBER OF NATURAL ABUNDANCES =, I5,  MAI 8350
     *     28H IS GREATER THAN ALLOWABLE =, I5)                         MAI 8360
99989 FORMAT (52H **WARNING 9608** NUMBER OF FISSION PRODUCT YIELDS =,  MAI 8370
     *     I5, 28H IS GREATER THAN ALLOWABLE =, I5)                     MAI 8380
99988 FORMAT (43H **WARNING 9605** NUMBER OF PHOTON YIELDS =, I5,       MAI 8390
     *     28H IS GREATER THAN ALLOWABLE =, I5)                         MAI 8400
99987 FORMAT (1H , 31HCON COMMANDS FOR GTO NOT FOUND , 3I5, 1PE10.2,    MAI 8410
     *     2I5)                                                         MAI 8420
99986 FORMAT (1H1)                                                      MAI 8430
99985 FORMAT (48H0NO CON CARD FOUND FOR DO LOOP IN INSTRUCTION # , I5)  MAI 8440
99984 FORMAT (20H NUMBER OF COMMAND= , I4, 23H ; THIS IS INSTRUCTION ,  MAI 8450
     *     I3, 19H OUT OF A TOTAL OF , I4, 1X, 1H*, A4, 12H* INSTRUCTIO,MAI 8460
     *     11HNS. UNIT=15)                                              MAI 8470
99983 FORMAT (2X, I2, 4(1X, I6, 1X, 1PE11.4))                           MAI 8480
      END                                                               MAI 8490
      SUBROUTINE OUTPUT(L, NUCL, Q, FG, AMPC, WMPC, TOCAP, FISS, DIS,   OUT   10
     *     XNEW, MX, XN, XTOT, XA, XSUM, XSAV, FX, FFA, XCUM)           OUT   20
      INTEGER ELSF(2), ELE, STA                                         OUT   30
      INTEGER BASIS                                                     OUT   40
      DIMENSION XCUM(12,3,11)                                           OUT   50
      DIMENSION NUCL(1), Q(1), FG(1), AMPC(1), WMPC(1), TOCAP(1),       OUT   60
     *     FISS(1), DIS(1), XNEW(MX,1), FFA(1)                          OUT   70
      DIMENSION XN(100,13), XTOT(14), XA(14), XNT(14), XSUM(14),        OUT   80
     *     XSAV(14), FX(14)                                             OUT   90
      DIMENSION OUTTYP(6,3)                                             OUT  100
      COMMON /ASSN/ IL(1), LA, LB, LC, LD, LE, LF, LG, LH, LI, LJ, LK,  OUT  110
     *     LL, LM, LN, LO, LP, LQ, LR, LS, LT, LU, LV, LW, LX           OUT  120
      COMMON /INST/ INSTR(300), NINSTR(300), IFIND(300), NNTYP(32),     OUT  130
     *     LAB(300)                                                     OUT  140
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  OUT  150
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           OUT  160
     *     TFPBAV(4), IPHMAX                                            OUT  170
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 OUT  180
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    OUT  190
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  OUT  200
     *     ITITLE(18)                                                   OUT  210
      COMMON /OUTOPT/ NOUT(20,4), NTO(80),MCOUT                         OUT  220
      COMMON /SZ/ NSIZE(25)                                             OUT  230
      COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)   OUT  240
C      LEVEL 2,NUCL,Q,FG,AMPC,WMPC,TOCAP,FISS,DIS,XNEW,XN,XTOT,XA,XSUM, OUT  250
C    1       XSAV,FX,XCUM                                               OUT  260
      DATA OUTTYP /4H*ACT,4HIVAT,4HION ,4HPROD,4HUCTS,4H****,4H*ACT,    OUT  270
     *     4HINID,4HES +,4H DAU,4HGHTE,4HRS**,4H*FIS,4HSION,4H PRO,     OUT  280
     *     4HDUCT,4HS***,4H****/                                        OUT  290
      DATA ELSF /2H S,2HSF/                                             OUT  300
      IF (NPAGE.EQ.0) NPAGE = IPAGE(1)                                  OUT  310
      DO 30 M=1,MOUT                                                    OUT  320
         DO 20 JTAB=1,11                                                OUT  330
            DO 10 KIND=1,3                                              OUT  340
               XCUM(M,KIND,JTAB) = 0.0                                  OUT  350
   10       CONTINUE                                                    OUT  360
   20    CONTINUE                                                       OUT  370
   30 CONTINUE                                                          OUT  380
      NTO(11) = 8                                                       OUT  390
      NTO(12) = 8                                                       OUT  400
      NTO(35) = 8                                                       OUT  410
      NTO(36) = 8                                                       OUT  420
      NTO(59) = 8                                                       OUT  430
      NTO(60) = 8                                                       OUT  440
      JO = IOUT                                                         OUT  450
      IDM(1) = 0                                                        OUT  460
      DO 40 J=2,4                                                       OUT  470
         TFPBAV(J) = 1.0                                                OUT  480
   40 CONTINUE                                                          OUT  490
      DO 50 M=1,MOUT                                                    OUT  500
         IF (OUTFPB(M,8).LT.ERR .AND. OUTFPB(M,9).LT.ERR .AND.          OUT  510
     *        OUTFPB(M,10).LT.ERR) GO TO 50                             OUT  520
         TFPBAV(2) = OUTFPB(M,9)                                        OUT  530
         TFPBAV(3) = OUTFPB(M,10)                                       OUT  540
         TFPBAV(4) = OUTFPB(M,8)                                        OUT  550
         GO TO 60                                                       OUT  560
   50 CONTINUE                                                          OUT  570
   60 CONTINUE                                                          OUT  580
      CALL TOC(NPAGE, 04, 01, 00, 00, ITITLE, IOUT, NTOC, NREC)         OUT  590
      WRITE (IOUT,99994) IOUT, NPAGE                                    OUT  600
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 04, 00, 01, 08, ITITLE, IOUT, NTOC, NREC)         OUT  610
      NPAGE = IPAGE(1)                                                  OUT  620
      WRITE (JO,99992) ITITLE                                           OUT  630
      IF (NOUT(IFIND(L),3).GT.0) WRITE (JO,99979) NREC                  OUT  640
      WRITE (JO,99993) TFPBAV(3), TFPBAV(4), TFPBAV(2)                  OUT  650
      WRITE (JO,99990) BASIS                                            OUT  660
      CALL HEAD(MOUT, 1, IOUT)                                          OUT  670
      WRITE (JO,99989) (OUTFPB(K,1),K=1,MOUT)                           OUT  680
      WRITE (JO,99988) (OUTFPB(K,2),K=1,MOUT)                           OUT  690
      WRITE (JO,99987) (OUTFPB(K,3),K=1,MOUT)                           OUT  700
      WRITE (JO,99986) (OUTFPB(K,4),K=1,MOUT)                           OUT  710
      WRITE (JO,99985) (OUTFPB(K,5),K=1,MOUT)                           OUT  720
      WRITE (JO,99984) (OUTFPB(K,6),K=1,MOUT)                           OUT  730
      WRITE (JO,99983) (OUTFPB(K,7),K=1,MOUT)                           OUT  740
      WRITE (JO,99982) (OUTFPB(K,08),K=1,MOUT)                          OUT  750
      WRITE (JO,99981) (OUTFPB(K,09),K=1,MOUT)                          OUT  760
      WRITE (JO,99980) (OUTFPB(K,10),K=1,MOUT)                          OUT  770
      WRITE (IOUT,99991) (M,NSIZE(M),M=1,12), (NSIZE(N),N=13,21)        OUT  780
C                                                                       OUT  800
C    DO FOR EACH TYPE OF NUCLIDE                                        OUT  810
C                                                                       OUT  820
      DO 100 KIND=1,3                                                   OUT  830
         ELE(16) = ELSF(1)                                              OUT  840
         IF (KIND.EQ.2) ELE(16) = ELSF(2)                               OUT  850
         ILOW = 0                                                       OUT  860
         DO 70 J=1,KIND                                                 OUT  870
            ILOW = ILOW + IDM(J)                                        OUT  880
   70    CONTINUE                                                       OUT  890
         IHI = ILOW + IDM(KIND+1)                                       OUT  900
         IF (ILOW.EQ.IHI) GO TO 100                                     OUT  910
         ILOW = ILOW + 1                                                OUT  920
         WRITE (NTOC,99995) ((OUTTYP(J,KIND),J=1,6),K=1,4)              OUT  930
C                                                                       OUT  940
C    DO FOR EACH TABLE                                                  OUT  950
C                                                                       OUT  960
         DO 90 JTAB=1,NTAB                                              OUT  970
            JJTAB = NTAB*(KIND-1) + JTAB                                OUT  980
            NNN = NTO(JJTAB)                                            OUT  990
            IF (NNN.EQ.8) GO TO 90                                      OUT 1000
            IF (JTAB.GT.20 .AND. JTAB.LT.25 .AND. KIND.NE.2) GO TO 90   OUT 1010
            CALL OUT1(NUCL, LNOW, XNEW, MX, DIS, Q, FG, AMPC, WMPC,     OUT 1020
     *           TOCAP, FISS, FX, XSAV, XN, XA, XTOT, XSUM, J1, J2, J3, OUT 1030
     *           J4, J5, J6, JTAB, ILOW, IHI, KIND, L, FFA)             OUT 1040
            CALL OUT2(NUCL, LNOW, XNEW, MX, DIS, Q, FG, AMPC, WMPC,     OUT 1050
     *           TOCAP, FISS, FX, XSAV, XN, XA, XTOT, XSUM, J1, J2, J3, OUT 1060
     *           J4, J5, J6, JTAB, ILOW, IHI, KIND, L, FFA)             OUT 1070
            IF (JTAB.LT.2 .OR. MOD(JTAB,2).EQ.0) GO TO 90               OUT 1080
            LTAB = (JTAB-2)/2 + 1                                       OUT 1090
            DO 80 M=1,MOUT                                              OUT 1100
               IF (KIND.NE.2) XCUM(M,1,LTAB) = XTOT(M) + XCUM(M,1,LTAB) OUT 1110
               IF (KIND.NE.1) XCUM(M,2,LTAB) = XTOT(M) + XCUM(M,2,LTAB) OUT 1120
               XCUM(M,3,LTAB) = XTOT(M) + XCUM(M,3,LTAB)                OUT 1130
   80       CONTINUE                                                    OUT 1140
            WRITE (JO,99999)                                            OUT 1150
            WRITE (JO,99998) (XCUM(M,1,LTAB),M=1,MOUT)                  OUT 1160
            WRITE (JO,99997) (XCUM(M,2,LTAB),M=1,MOUT)                  OUT 1170
            WRITE (JO,99996) (XCUM(M,3,LTAB),M=1,MOUT)                  OUT 1180
   90    CONTINUE                                                       OUT 1190
  100 CONTINUE                                                          OUT 1200
      RETURN                                                            OUT 1210
99999 FORMAT (1H0, 40X, 23HCUMULATIVE TABLE TOTALS, /)                  OUT 1220
99998 FORMAT (1H , 5HAP+FP, 5X,1P,12E10.3)                              OUT 1230
99997 FORMAT (1H , 6HACT+FP, 4X,1P,12E10.3)                             OUT 1240
99996 FORMAT (1H , 9HAP+ACT+FP, 1X,1P,12E10.3)                          OUT 1250
99995 FORMAT (1H , T20, 24A4)                                           OUT 1260
99994 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        OUT 1270
99993 FORMAT (1H , 4X, 6HPOWER=, 1PE12.5, 13H MW,  BURNUP=, 1PE12.5,    OUT 1280
     *     12H MWD,  FLUX=, 1PE9.2, 12H N/CM**2-SEC)                    OUT 1290
99992 FORMAT (1H , /, 5X, 1H*, 18A4)                                    OUT 1300
99991 FORMAT (/, /, /, 17H SIZE OF MMAX(I):, 2(6(6H MMAX=, I2, 2X, 2H#=,OUT 1310
     *     I4, 3X), /), /, 35H THE NUMBER OF NON-ZERO TERMS IN A=, I4,  OUT 1320
     *     /, 47H THE NUMBER OF NON-ZERO FISSION PRODUCT YIELDS=, I4,   OUT 1330
     *     /, 7H ILITE=, I4, 5X, 5HIACT=, I4, 5X, 4HIFP=, I4, 5X,       OUT 1340
     *     5HITOT=, I4, /, 41H THE NUMBER OF NON-ZERO NATURAL ABUNDANCE,OUT 1350
     *     2HS=, I4, /, 38H THE NUMBER OF NON-ZERO PHOTON YIELDS=, I5,  OUT 1360
     *     /, 35H THE MAXIMUM NUMBER OF TERMS IN AP=, I5)               OUT 1370
99990 FORMAT (1H0, 45X, 26HREACTIVITY AND BURNUP DATA, /, 1H , 25X,     OUT 1380
     *     6HBASIS=, 18A4, /)                                           OUT 1390
99989 FORMAT (11H TIME, SEC , 1X, 12(1PE9.2, 1X))                       OUT 1400
99988 FORMAT (11H NEUT. FLUX, 1X, 12(1PE9.2, 1X))                       OUT 1410
99987 FORMAT (11H SP POW,MW , 1X, 12(1PE9.2, 1X))                       OUT 1420
99986 FORMAT (11H BURNUP,MWD, 1X, 12(1PE9.2, 1X))                       OUT 1430
99985 FORMAT (11H K INFINITY, 2X, 12(F8.5, 2X))                         OUT 1440
99984 FORMAT (11H NEUT PRODN, 1X, 12(1PE9.2, 1X))                       OUT 1450
99983 FORMAT (11H NEUT DESTN, 1X, 12(1PE9.2, 1X))                       OUT 1460
99982 FORMAT (11H TOT BURNUP, 1X, 12(1PE9.2, 1X))                       OUT 1470
99981 FORMAT (11H AVG N FLUX, 1X, 12(1PE9.2, 1X))                       OUT 1480
99980 FORMAT (11H AVG SP POW, 1X, 12(1PE9.2, 1X))                       OUT 1490
99979 FORMAT (1H+, 80X, 9HRECYCLE #, I3)                                OUT 1500
      END                                                               OUT 1510
      SUBROUTINE OUT1(NUCL, LNOW, XNEW, MX, DIS, Q, FG, AMPC, WMPC,     OUT   10
     *     TOCAP, FISS, FX, XSAV, XN, XA, XTOT, XSUM, J1, J2, J3, J4,   OUT   20
     *     J5, J6, JTAB, ILOW, IHI, KIND, L, FFA)                       OUT   30
      DIMENSION NUCL(1), XNEW(MX,1), DIS(1), Q(1), FG(1), AMPC(1),      OUT   40
     *     WMPC(1), TOCAP(1), FISS(1), XN(100,13), FFA(1)               OUT   50
      DIMENSION XSAV(14), XA(14), XTOT(14), XSUM(14), FX(14)            OUT   60
      COMMON /ASSN/ IL(1), LA, LB, LC, LD, LE, LF, LG, LH, LI, LJ, LK,  OUT   70
     *     LL, LM, LN, LO, LP, LQ, LR, LS, LT, LU, LV, LW, LX           OUT   80
      COMMON /BDAT/ ZAP(130), KZAP(12), JZAP(6,26), NTB(3,4),           OUT   90
     *     INSTR0(32), RRCUT(30), NDI(11,7), NDL(6,5), NDO(9,5)         OUT  100
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  OUT  110
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           OUT  120
     *     TFPBAV(4), IPHMAX                                            OUT  130
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    OUT  140
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  OUT  150
     *     ITITLE(18)                                                   OUT  160
      COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)   OUT  170
C      LEVEL 2,NUCL,XNEW,DIS,Q,FG,AMPC,WMPC,TOCAP,FISS,FX,XSAV,XN,XA,   OUT  180
C           XTOT,XSUM                                                   OUT  190
      ASSIGN 40 TO LA                                                   OUT  200
      ASSIGN 50 TO LB                                                   OUT  210
      ASSIGN 60 TO LG                                                   OUT  220
      ASSIGN 70 TO LI                                                   OUT  230
      ASSIGN 80 TO LK                                                   OUT  240
      ASSIGN 90 TO LM                                                   OUT  250
      ASSIGN 100 TO LO                                                  OUT  260
      ASSIGN 110 TO LQ                                                  OUT  270
      ASSIGN 120 TO LS                                                  OUT  280
      ASSIGN 130 TO LU                                                  OUT  290
      ASSIGN 140 TO LW                                                  OUT  300
      LC = LA                                                           OUT  310
      LD = LA                                                           OUT  320
      LE = LB                                                           OUT  330
      LF = LB                                                           OUT  340
      LH = LG                                                           OUT  350
      LJ = LI                                                           OUT  360
      LL = LK                                                           OUT  370
      LN = LM                                                           OUT  380
      LP = LO                                                           OUT  390
      LR = LQ                                                           OUT  400
      LT = LS                                                           OUT  410
      LV = LU                                                           OUT  420
      LX = LW                                                           OUT  430
C     LNOW = IL(JTAB+1)                                                 OUT  440
      LNOW = JTAB
      DO 20 M=1,MOUT                                                    OUT  450
         XSAV(M) = 1.0E-35                                              OUT  460
         XTOT(M) = 0.0                                                  OUT  470
         DO 10 J=1,100                                                  OUT  480
            XN(J,M) = 0.0                                               OUT  490
   10    CONTINUE                                                       OUT  500
   20 CONTINUE                                                          OUT  510
      DO 30 M=1,MOUT                                                    OUT  520
         FX(M) = OUTFPB(M,2)                                            OUT  530
         IF (FX(M).LT.ERR) FX(M) = TFPBAV(2)                            OUT  540
         IF (FX(M).LT.ERR) FX(M) = 1.0                                  OUT  550
   30 CONTINUE                                                          OUT  560
      DO 170 I=ILOW,IHI                                                 OUT  570
         JAM = NUCL(I)/10                                               OUT  580
         JAN = JAM/1000                                                 OUT  590
         WT = JAM - (1000*JAN)                                          OUT  600
         DO 160 M=1,MOUT                                                OUT  610
            XA(M) = 0.0                                                 OUT  620
            GO TO       (40, 50, 40, 40, 50, 50, 60, 60, 70, 70, 80,    OUT  630
     *           80, 90, 90, 100, 100, 110, 110, 120, 120, 130, 130,    OUT  640
     *           140, 140), JTAB                                        OUT  650
   40       XA(M) = XNEW(M,I)                                           OUT  660
            GO TO 150                                                   OUT  670
   50       XA(M) = XNEW(M,I)*WT                                        OUT  680
            GO TO 150                                                   OUT  690
   60       XA(M) = XNEW(M,I)*DIS(I)*1.6283E+13                         OUT  700
            GO TO 150                                                   OUT  710
   70       XA(M) = XNEW(M,I)*Q(I)                                      OUT  720
            GO TO 150                                                   OUT  730
   80       XA(M) = XNEW(M,I)*FG(I)                                     OUT  740
            GO TO 150                                                   OUT  750
   90       XA(M) = XNEW(M,I)*AMPC(I)                                   OUT  760
            GO TO 150                                                   OUT  770
  100       XA(M) = XNEW(M,I)*WMPC(I)                                   OUT  780
            GO TO 150                                                   OUT  790
  110       XA(M) = XNEW(M,I)*WT*CEMTOX(JAN)                            OUT  800
            GO TO 150                                                   OUT  810
  120       XA(M) = XNEW(M,I)*TOCAP(I)*0.6023*FX(M)                     OUT  820
            GO TO 150                                                   OUT  830
  130       XA(M) = XNEW(M,I)*FISS(I-ILITE)*0.6022169*FX(M)             OUT  840
            GO TO 150                                                   OUT  850
  140       IF (DIS(I).LE.1.0000E+08) XA(M) = XNEW(M,I)*DIS(I)*         OUT  860
     *           1.6283E+13*FFA(I-ILITE)                                OUT  870
  150       XSAV(M) = XSAV(M) + XA(M)                                   OUT  880
            XN(JAN,M) = XN(JAN,M) + XA(M)                               OUT  890
  160    CONTINUE                                                       OUT  900
  170 CONTINUE                                                          OUT  910
      IF (JTAB.LT.3) GO TO 210                                          OUT  920
      DO 200 M=1,MOUT                                                   OUT  930
         IF (XSAV(M).LT.ERR) GO TO 180                                  OUT  940
         XTOT(M) = XSAV(M)                                              OUT  950
         XSAV(M) = 1.0/XSAV(M)                                          OUT  960
         XN(100,M) = XSAV(M)                                            OUT  970
  180    IF (MOD(JTAB,2).NE.0) XSAV(M) = 1.0                            OUT  980
         DO 190 LEL=1,99                                                OUT  990
            XN(LEL,M) = XN(LEL,M)*XSAV(M)                               OUT 1000
  190    CONTINUE                                                       OUT 1010
  200 CONTINUE                                                          OUT 1020
  210 CONTINUE                                                          OUT 1030
      J1 = JZAP(1,JTAB)                                                 OUT 1040
      J2 = JZAP(2,JTAB)                                                 OUT 1050
      J3 = JZAP(3,JTAB)                                                 OUT 1060
      J4 = JZAP(4,JTAB)                                                 OUT 1070
      J5 = JZAP(5,JTAB)                                                 OUT 1080
      J6 = JZAP(6,JTAB)                                                 OUT 1090
      RETURN                                                            OUT 1100
      END                                                               OUT 1110
      SUBROUTINE OUT2(NUCL, LNOW, XNEW, MX, DIS, Q, FG, AMPC, WMPC,     OUT   10
     *     TOCAP, FISS, FX, XSAV, XN, XA, XTOT, XSUM, J1, J2, J3, J4,   OUT   20
     *     J5, J6, JTAB, ILOW, IHI, KIND, L, FFA)                       OUT   30
      LOGICAL P(3), ESKIP(100)                                          OUT   40
      INTEGER ELE, STA                                                  OUT   50
      INTEGER BASIS                                                     OUT   60
      DIMENSION TATY(5,3)                                               OUT   70
      DIMENSION NUCL(1), XNEW(MX,1), DIS(1), Q(1), FG(1), AMPC(1),      OUT   80
     *     WMPC(1), TOCAP(1), FISS(1), XN(100,13), FFA(1)               OUT   90
      DIMENSION FX(1), XSAV(1), XA(1), XTOT(1), XSUM(1)                 OUT  100
      COMMON /ASSN/ IL(1), LA, LB, LC, LD, LE, LF, LG, LH, LI, LJ, LK,  OUT  110
     *     LL, LM, LN, LO, LP, LQ, LR, LS, LT, LU, LV, LW, LX           OUT  120
      COMMON /BDAT/ ZAP(130), KZAP(12), JZAP(6,26), NTB(3,4),           OUT  130
     *     INSTR0(32), RRCUT(30), NDI(11,7), NDL(6,5), NDO(9,5)         OUT  140
      COMMON /INST/ INSTR(300), NINSTR(300), IFIND(300), NNTYP(32),     OUT  150
     *     LAB(300)                                                     OUT  160
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  OUT  170
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           OUT  180
     *     TFPBAV(4), IPHMAX                                            OUT  190
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 OUT  200
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    OUT  210
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  OUT  220
     *     ITITLE(18)                                                   OUT  230
      COMMON /OUTOPT/ NOUT(20,4), NTO(80),MCOUT                         OUT  240
C      LEVEL 2,NUCL,XNEW,DIS,Q,FG,AMPC,WMPC,TOCAP,FISS,FX,XSAV,XN,XA,   OUT  250
C    1  XTOT,XSUM                                                       OUT  260
      DATA TATY /4HACTI,4HVATI,4HON P,4HRODU,4HCTS ,4HACTI,4HNIDE,      OUT  280
     *     4HS+DA,4HUGHT,4HERS ,4HFISS,4HION ,4HPROD,4HUCTS,4H    /     OUT  290
      MSTAR = NOUT(IFIND(L),4)                                          OUT  300
      JO = IOUT                                                         OUT  310
      ASSIGN 50 TO LA                                                   OUT  320
      ASSIGN 60 TO LC                                                   OUT  330
      ASSIGN 70 TO LE                                                   OUT  340
      ASSIGN 80 TO LG                                                   OUT  350
      ASSIGN 90 TO LI                                                   OUT  360
      ASSIGN 100 TO LK                                                  OUT  370
      ASSIGN 110 TO LM                                                  OUT  380
      ASSIGN 120 TO LO                                                  OUT  390
      ASSIGN 130 TO LQ                                                  OUT  400
      ASSIGN 140 TO LS                                                  OUT  410
      ASSIGN 150 TO LU                                                  OUT  420
      ASSIGN 160 TO LW                                                  OUT  430
      LB = LA                                                           OUT  440
      LD = LC                                                           OUT  450
      LF = LE                                                           OUT  460
      LH = LG                                                           OUT  470
      LJ = LI                                                           OUT  480
      LL = LK                                                           OUT  490
      LN = LM                                                           OUT  500
      LP = LO                                                           OUT  510
      LR = LQ                                                           OUT  520
      LT = LS                                                           OUT  530
      LV = LU                                                           OUT  540
      LX = LW                                                           OUT  550
C
      DO 1 IE=1,99                                                      OUT  561
      ESKIP(IE)=.TRUE.                                                  OUT  562
    1 CONTINUE                                                          OUT  563
      DO 2 I=ILOW,IHI                                                   OUT  564
      JAN=NUCL(I)/10000                                                 OUT  565
      ESKIP(JAN)=.FALSE.                                                OUT  566
    2 CONTINUE                                                          OUT  567
C
C     LNOW = IL(JTAB+1)                                                 OUT  560
      LNOW = JTAB
      DO 360 KTAB=1,3                                                   OUT  570
         IF (JTAB.LT.3 .AND. KTAB.NE.1) GO TO 360                       OUT  580
         JJTAB = NTAB*(KIND-1) + JTAB                                   OUT  590
         NNN = NTO(JJTAB)                                               OUT  600
         IF (NNN.EQ.8) GO TO 360                                        OUT  610
         DO 10 MM=1,3                                                   OUT  620
            P(MM) = .TRUE.                                              OUT  630
   10    CONTINUE                                                       OUT  640
         DO 20 MM=1,4                                                   OUT  650
            IF (NNN.EQ.NTB(KTAB,MM)) P(KTAB) = .FALSE.                  OUT  660
   20    CONTINUE                                                       OUT  670
         IF (P(KTAB)) GO TO 360                                         OUT  680
         CALL TOC(NPAGE, 04, 01, KTAB, JTAB, ITITLE, IOUT, NTOC, NREC)  OUT  690
         K1 = KZAP(4*(KTAB-1)+1)                                        OUT  700
         K2 = K1 + 3                                                    OUT  710
         DO 30 M=1,MOUT                                                 OUT  720
            XSUM(M) = 0.0                                               OUT  730
            XTOT(M) = 0.0                                               OUT  740
   30    CONTINUE                                                       OUT  750
         IF (.NOT.P(2)) GO TO 320                                       OUT  760
         KOUNT = 0                                                      OUT  770
         DO 250 I=ILOW,IHI                                              OUT  780
            DO 40 M=1,MOUT                                              OUT  790
               XA(M) = 0.0                                              OUT  800
   40       CONTINUE                                                    OUT  810
            CALL NOAH(NUCL(I))                                          OUT  820
            JAM = NUCL(I)/10                                            OUT  830
            JAN = JAM/1000                                              OUT  840
            WT = JAM - (1000*JAN)                                       OUT  860
C                                                                       OUT  870
C    DO FOR EACH TIME                                                   OUT  880
C                                                                       OUT  890
            DO 180 M=1,MOUT                                             OUT  900
               GO TO       (50, 50, 60, 60, 70, 70, 80, 80, 90, 90,     OUT  910
     *              100, 100, 110, 110, 120, 120, 130, 130, 140, 140,   OUT  920
     *              150, 150, 160, 160), JTAB                           OUT  930
   50          IF (XN(JAN,M).LT.ERR) GO TO 170                          OUT  940
               XA(M) = XNEW(M,I)/XN(JAN,M)                              OUT  950
               IF (JTAB.GT.1) XA(M) = XA(M)*WT                          OUT  960
               GO TO 170                                                OUT  970
   60          XA(M) = XNEW(M,I)*XSAV(M)                                OUT  980
               GO TO 170                                                OUT  990
   70          XA(M) = XNEW(M,I)*WT*XSAV(M)                             OUT 1000
               GO TO 170                                                OUT 1010
   80          XA(M) = XNEW(M,I)*DIS(I)*1.6283E+13*XSAV(M)              OUT 1020
               GO TO 170                                                OUT 1030
   90          XA(M) = XNEW(M,I)*Q(I)*XSAV(M)                           OUT 1040
               GO TO 170                                                OUT 1050
  100          XA(M) = XNEW(M,I)*FG(I)*XSAV(M)                          OUT 1060
               GO TO 170                                                OUT 1070
  110          XA(M) = XNEW(M,I)*AMPC(I)*XSAV(M)                        OUT 1080
               GO TO 170                                                OUT 1090
  120          XA(M) = XNEW(M,I)*WMPC(I)*XSAV(M)                        OUT 1100
               GO TO 170                                                OUT 1110
  130          XA(M) = XNEW(M,I)*WT*CEMTOX(JAN)*XSAV(M)                 OUT 1120
               GO TO 170                                                OUT 1130
  140          XA(M) = XNEW(M,I)*TOCAP(I)*0.6023*FX(M)*XSAV(M)          OUT 1140
               GO TO 170                                                OUT 1150
  150          XA(M) = XNEW(M,I)*FISS(I-ILITE)*6.022169E-01*FX(M)*      OUT 1160
     *              XSAV(M)                                             OUT 1170
               GO TO 170                                                OUT 1180
  160          IF (DIS(I).LE.1.0000E+08) XA(M) = XNEW(M,I)*DIS(I)*      OUT 1190
     *              1.6283E+13*FFA(I-ILITE)*XSAV(M)                     OUT 1200
  170          XTOT(M) = XTOT(M) + XA(M)                                OUT 1210
  180       CONTINUE                                                    OUT 1220
            IF (P(1)) GO TO 200                                         OUT 1230
            KOUNT = KOUNT + 1                                           OUT 1240
            IF (MOD(KOUNT-1,51).NE.0) GO TO 190   !REDUCED 52 TO 51     OUT 1250
            WRITE (IOUT,99998) IOUT, NPAGE                              OUT 1260
      CALL IDENTIFY (IOUT, 0)                                           PC
            NPAGE = IPAGE(1)                                            OUT 1270
            WRITE (JO,99995) ITITLE                                     OUT 1280
            IF (NOUT(IFIND(L),3).GT.0) WRITE (JO,99988) NREC            OUT 1290
            WRITE (JO,99987) (TATY(N,KIND),N=1,5)                       OUT 1300
            WRITE (JO,99997) TFPBAV(3), TFPBAV(4), TFPBAV(2)            OUT 1310
            WRITE (IOUT,99996) JTAB, (ZAP(N),N=K1,K2),                  OUT 1320
     *           (ZAP(N),N=J1,J2), (ZAP(N),N=J3,J4), (ZAP(N),N=J5,J6)   OUT 1330
            WRITE (IOUT,99999) BASIS                                    OUT 1340
            CALL HEAD(MOUT, 1, IOUT)                                    OUT 1350
  190       WRITE (IOUT,99993) NAME, (XA(M),M=1,MOUT)                   OUT 1360
  200       IF (P(3)) GO TO 250                                         OUT 1370
            DO 210 M=1,MOUT                                             OUT 1380
               IF (MSTAR.GT.0 .AND. MSTAR.NE.M .OR. MSTAR.LT.0 .AND.    OUT 1390
     *              (-MSTAR).EQ.M) GO TO 210                            OUT 1400
               IF (MOD(JTAB,2).NE.0 .AND. (XA(M)*XN(100,M))             OUT 1410
     *              .GE.CUTOFF(JTAB)) GO TO 220                         OUT 1420
               IF (MOD(JTAB,2).EQ.0 .AND. XA(M).GE.CUTOFF(JTAB)) GO TO  OUT 1430
     *              220                                                 OUT 1440
  210       CONTINUE                                                    OUT 1450
            GO TO 250                                                   OUT 1460
  220       KOUNT = KOUNT + 1                                           OUT 1470
            IF (MOD(KOUNT-1,51).NE.0) GO TO 230  !REDUCED 52 TO 51      OUT 1480
            WRITE (IOUT,99998) IOUT, NPAGE                              OUT 1490
      CALL IDENTIFY (IOUT, 0)                                           PC
            NPAGE = IPAGE(1)                                            OUT 1500
            WRITE (JO,99995) ITITLE                                     OUT 1510
            IF (NOUT(IFIND(L),3).GT.0) WRITE (JO,99988) NREC            OUT 1520
            WRITE (JO,99987) (TATY(N,KIND),N=1,5)                       OUT 1530
            WRITE (JO,99997) TFPBAV(3), TFPBAV(4), TFPBAV(2)            OUT 1540
            WRITE (IOUT,99996) JTAB, (ZAP(N),N=K1,K2),                  OUT 1550
     *           (ZAP(N),N=J1,J2), (ZAP(N),N=J3,J4), (ZAP(N),N=J5,J6)   OUT 1560
            WRITE (IOUT,99999) BASIS                                    OUT 1570
            CALL HEAD(MOUT, 1, IOUT)                                    OUT 1580
  230       WRITE (IOUT,99993) NAME, (XA(M),M=1,MOUT)                   OUT 1590
            DO 240 M=1,MOUT                                             OUT 1600
               XSUM(M) = XSUM(M) + XA(M)                                OUT 1610
  240       CONTINUE                                                    OUT 1620
  250    CONTINUE                                                       OUT 1630
C                                                                       OUT 1640
C    WRITE OUT TOTALS FOR EACH TABLE                                    OUT 1650
C                                                                       OUT 1660
         IF (P(3)) GO TO 310                                            OUT 1670
         WRITE (JO,99994) (XSUM(M),M=1,MOUT)                            OUT 1680
         WRITE (JO,99990) (XTOT(M),M=1,MOUT)                            OUT 1690
         DO 260 M=1,MOUT                                                OUT 1700
            XSUM(M) = 0.0                                               OUT 1710
  260    CONTINUE                                                       OUT 1720
C    ELEMENT SUMMARY TABLES                                             OUT 1730
         WRITE (IOUT,99998) IOUT, NPAGE                                 OUT 1740
      CALL IDENTIFY (IOUT, 0)                                           PC
         NPAGE = IPAGE(1)                                               OUT 1750
         WRITE (JO,99995) ITITLE                                        OUT 1760
         IF (NOUT(IFIND(L),3).GT.0) WRITE (JO,99988) NREC               OUT 1770
         WRITE (JO,99987) (TATY(N,KIND),N=1,5)                          OUT 1780
         WRITE (JO,99997) TFPBAV(3), TFPBAV(4), TFPBAV(2)               OUT 1790
         WRITE (IOUT,99996) JTAB, (ZAP(N),N=K1,K2), (ZAP(N),N=J1,J2),   OUT 1800
     *        (ZAP(N),N=J3,J4), (ZAP(N),N=J5,J6)                        OUT 1810
         WRITE (IOUT,99999) BASIS                                       OUT 1820
         CALL HEAD(MOUT, 1, IOUT)                                       OUT 1840
         DO 300 J=1,99                                                  OUT 1850
            DO 270 MM=1,MOUT                                            OUT 1860
               IF (MSTAR.GT.0 .AND. MSTAR.NE.MM .OR. MSTAR.LT.0 .AND.   OUT 1870
     *              (-MSTAR).EQ.MM) GO TO 270                           OUT 1880
               IF (MOD(JTAB,2).NE.0 .AND. (XN(J,MM)*XN(100,MM))         OUT 1890
     *              .GE.CUTOFF(JTAB)) GO TO 280                         OUT 1900
               IF (MOD(JTAB,2).EQ.0 .AND. XN(J,MM).GE.CUTOFF(JTAB)) GO  OUT 1910
     *              TO 280                                              OUT 1920
  270       CONTINUE                                                    OUT 1930
            GO TO 300                                                   OUT 1940
  280       WRITE (IOUT,99989) ELE(J), (XN(J,M),M=1,MOUT)               OUT 1950
            DO 290 M=1,MOUT                                             OUT 1960
               XSUM(M) = XSUM(M) + XN(J,M)                              OUT 1970
  290       CONTINUE                                                    OUT 1980
  300    CONTINUE                                                       OUT 1990
         WRITE (JO,99994) (XSUM(M),M=1,MOUT)                            OUT 2000
         WRITE (JO,99990) (XTOT(M),M=1,MOUT)                            OUT 2010
         GO TO 360                                                      OUT 2020
  310    IF (P(1)) GO TO 320                                            OUT 2030
         WRITE (JO,99992) (XTOT(M),M=1,MOUT)                            OUT 2040
         GO TO 360                                                      OUT 2050
C    ****PRINT ELEMENT TABLES  -- COMPLETE                              OUT 2060
  320    IF (P(2)) GO TO 360                                            OUT 2070
         IF (JTAB.LT.3) GO TO 360                                       OUT 2080
         IPR = -1                                                       OUT 2090
         DO 350 J=1,99                                                  OUT 2100
            IF (ESKIP(J)) GO TO 350                                     OUT 2110
            IPR = IPR + 1                                               OUT 2120
            IF (MOD(IPR,50).NE.0) GO TO 330                             OUT 2130
            WRITE (IOUT,99998) IOUT, NPAGE                              OUT 2140
      CALL IDENTIFY (IOUT, 0)                                           PC
            NPAGE = IPAGE(1)                                            OUT 2150
            WRITE (JO,99995) ITITLE                                     OUT 2160
            IF (NOUT(IFIND(L),3).GT.0) WRITE (JO,99988) NREC            OUT 2170
            WRITE (JO,99987) (TATY(N,KIND),N=1,5)                       OUT 2180
            WRITE (JO,99997) TFPBAV(3), TFPBAV(4), TFPBAV(2)            OUT 2190
            WRITE (IOUT,99996) JTAB, (ZAP(N),N=K1,K2),                  OUT 2200
     *           (ZAP(N),N=J1,J2), (ZAP(N),N=J3,J4), (ZAP(N),N=J5,J6)   OUT 2210
            WRITE (IOUT,99999) BASIS                                    OUT 2220
            CALL HEAD(MOUT, 1, IOUT)                                    OUT 2230
  330       CONTINUE                                                    OUT 2240
            DO 340 M=1,MOUT                                             OUT 2250
               XTOT(M) = XTOT(M) + XN(J,M)                              OUT 2260
  340       CONTINUE                                                    OUT 2270
            WRITE (JO,99991) ELE(J), (XN(J,M),M=1,MOUT)                 OUT 2280
  350    CONTINUE                                                       OUT 2290
         WRITE (JO,99992) (XTOT(M),M=1,MOUT)                            OUT 2300
  360 CONTINUE                                                          OUT 2310
C                                                                       OUT 2320
C    FORMATS     FORMATS     FORMATS     FORMATS     FORMATS            OUT 2330
C                                                                       OUT 2340
      RETURN                                                            OUT 2350
99999 FORMAT (1H , 28X, 18A4)                                           OUT 2360
99998 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        OUT 2370
99997 FORMAT (1H , 4X, 6HPOWER=, 1PE12.5, 13H MW,  BURNUP=, 1PE12.5,    OUT 2380
     *     12H MWD,  FLUX=, 1PE9.2, 12H N/CM**2-SEC)                    OUT 2390
99996 FORMAT (1H0, 20X, I2, 1X, 25A4)                                   OUT 2400
99995 FORMAT (1H , /, 5X, 1H*, 18A4)                                    OUT 2410
99994 FORMAT (7H SUMTOT, 4X,1P,12E10.3)                                 OUT 2420
99993 FORMAT (1H , A2, I3, A1, 4X,1P,12E10.3)                           OUT 2430
99992 FORMAT (6H TOTAL, 5X,1P,12E10.3)                                  OUT 2440
99991 FORMAT (1H , A2, 8X,1P,12E10.3)                                   OUT 2450
99990 FORMAT (6H0TOTAL, 5X,1P,12E10.3)                                  OUT 2460
99989 FORMAT (1H , A2, 4X, 4X,1P,12E10.3)                               OUT 2470
99988 FORMAT (1H+, 80X, 9HRECYCLE #, I3)                                OUT 2480
99987 FORMAT (1H+, T100, 5A4)                                           OUT 2490
      END                                                               OUT 2500
      SUBROUTINE NUDAT1(NUCL, Q, FG, GENNEU, AMPC, WMPC, TOCAP, FISS,   NUD   10
     *     DIS, YIELD, COEFF, LC, NPROD, FFSF, FFA, NY, YY, NUCAB,      NUD   20
     *     NYIELD, LOC, NON0, KD, MMAX, KAP, AP, ABUND, NNSORS, LPUCNT) NUD   30
      INTEGER*2       NYIELD,LOC,NON0,KD,MMAX,KAP                       NUD   40
      DOUBLE PRECISION DNUM, UNIT, TENS                                 NUD   50
      INTEGER ELE, STA                                                  NUD   60
      INTEGER AREAD, BASIS, TITLE                                       NUD   70
      DIMENSION TYLD(11), NEWCX(11), NSRS(11), CAPT(11), TITLE(18),     NUD   80
     *     YZ(8)                                                        NUD   90
      DIMENSION NUCL(1), Q(1), FG(1), GENNEU(1), AMPC(1), WMPC(1),      NUD  100
     *     TOCAP(1), A(1), FISS(1), DIS(1), YIELD(1), COEFF(LC,1),      NUD  110
     *     NPROD(LC,1), FFSF(1), FFA(1), NY(1), YY(1), NUCAB(1),        NUD  120
     *     NYIELD(1), LOC(1), NON0(1), KD(1), MMAX(1), KAP(1), AP(1),   NUD  130
     *     ABUND(1)                                                     NUD  140
      COMMON /DECR/ DUM(1), DLAM, FB1, FP, FP1, FA, FT, FSF, LIBE, IU,  NUD  150
     *     ABUN, FN                                                     NUD  160
      COMMON /INST/ INSTR(300), NINSTR(300), IFIND(300), NNTYP(32),     NUD  170
     *     LAB(300)                                                     NUD  180
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  NUD  190
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           NUD  200
     *     TFPBAV(4), IPHMAX                                            NUD  210
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 NUD  220
      COMMON /MANE4I/ NINP(15,06), NIRR(150,04), NPRO(20,05),           NUD  230
     *     NCON(20,01), NDOL(02,04), NMOV(99,03), NADD(30,03),          NUD  240
     *     NHED(50,01), NWAC(09,02), NLPU(09,100), NPCH(54,03),         NUD  250
     *     NBUP(20), INUCL(5), NKEQ(03,05), NFAC(20,4), NGTO(10,2),     NUD  260
     *     NOFF(5)                                                              
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            NUD  270
     *     INCYLD(8,4), MYLD(4)                                         NUD  280
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    NUD  290
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  NUD  300
     *     ITITLE(18)                                                   NUD  310
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      NUD  320
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           NUD  330
      COMMON /SZ/ NSIZE(25)                                             NUD  340
C      LEVEL 2,NUCL,Q,FG,GENNEU,AMPC,WMPC,TOCAP,FISS,DIS,YIELD,COEFF,   NUD  350
C    1  NPROD,FFSF,NY,YY,NUCAB,NYIELD,LOC,NON0,KD,MMAX,KAP,AP,ABUND     NUD  360
      INPT = MLIB(8)                                                    NUD  370
      LPUCNT = 0                                                        NUD  380
      IAB = 1                                                           NUD  390
      ABUND(1) = 0.0                                                    NUD  400
      LOUT = IOUT                                                       NUD  410
      IF (MLIP(1).GT.0 .OR. MLIP(2).GT.0 .OR. MLIP(3).GT.0 .OR.         NUD  420
     *     MLIB(1).GT.0) CALL TOC(NPAGE, 3, 0, 0, 0, TITLE, IOUT, NTOC, NUD  430
     *     NREC)                                                        NUD  440
C    READ DECAY DATA FOR ALL NUCLIDES                                   NUD  450
      I = 0                                                             NUD  460
      DO 230 LITYP=2,4                                                  NUD  470
         LPUN = 1                                                       NUD  480
         L = 0                                                          NUD  490
         IF (MLIB(LITYP).EQ.0) GO TO 220                                NUD  500
         IF (MLIB(LITYP).GT.0) GO TO 30                                 NUD  510
         LPUCNT = 0                                                     NUD  520
         DO 10 MLVAR=2,LITYP                                            NUD  530
            IF (MLIB(MLVAR).LT.0) LPUCNT = LPUCNT + 1                   NUD  540
   10    CONTINUE                                                       NUD  550
         DO 20 J=1,100                                                  NUD  560
            MLPU(J) = NLPU(LPUCNT,J)                                    NUD  570
   20    CONTINUE                                                       NUD  580
   30    CONTINUE                                                       NUD  590
C    FIND CORRECT LIBRARY                                               NUD  600
         LIHUNT = IABS(MLIB(LITYP))                                     NUD  610
         DO 40 J=1,5000                                                 NUD  620
            LIBE = IREAD(1)                                             NUD  630
            IF (LIBE.EQ.LIHUNT) GO TO 50                                NUD  640
   40    CONTINUE                                                       NUD  650
         GO TO 260                                                      NUD  660
   50    CONTINUE                                                       NUD  670
         ICHR = 8                                                       NUD  680
         DO 60 ID=1,18                                                  NUD  690
            TITLE(ID) = AREAD(0)                                        NUD  700
   60    CONTINUE                                                       NUD  710
         IF (MLIP(1).GT.0) CALL TOC(NPAGE, 03, 01, 01, 00, TITLE, IOUT, NUD  720
     *        NTOC, NREC)                                               NUD  730
   70    CALL DECRED(I, NUCL, Q, FG, AMPC, WMPC)                        NUD  740
         IF (LIBE.LT.1) GO TO 220                                       NUD  750
         L = L + 1                                                      NUD  760
         IF (MLIB(LITYP).GT.0) GO TO 80                                 NUD  770
         IF (MLPU(LPUN).NE.NUCL(I)) GO TO 80                            NUD  780
         INPT = MLIB(9)                                                 NUD  790
         I = I - 1                                                      NUD  800
         CALL DECRED(I, NUCL, Q, FG, AMPC, WMPC)                        NUD  810
         LPUN = LPUN + 1                                                NUD  820
         INPT = MLIB(8)                                                 NUD  830
   80    M = 0                                                          NUD  840
         NUCAB(I) = 1                                                   NUD  850
         IF (ABUN.LT.ERR) GO TO 90                                      NUD  860
         IAB = IAB + 1                                                  NUD  870
         NUCAB(I) = IAB                                                 NUD  880
         ABUND(IAB) = ABUN                                              NUD  890
   90    CONTINUE                                                       NUD  900
         IF (DLAM.GT.ERR) DLAM = 0.6931472/(DLAM*RTIME(IU))             NUD  910
         DIS(I) = DLAM                                                  NUD  920
         IF (LITYP.NE.3) GO TO 110                                      NUD  930
         FFA(I-ILITE) = FA                                              NUD  940
         FFSF(I-ILITE) = FSF                                            NUD  950
         IF (NUCL(I).LT.992541) GO TO 110                               NUD  960
C    SET ES-254M/FM-254 DECAY                                           NUD  970
         IF (NUCL(I).EQ.992550) GO TO 100                               NUD  980
         NPROD(1,I) = 982540 + 1                                        NUD  990
         COEFF(1,I) = FP*DLAM                                           NUD 1000
         NPROD(2,I) = 982500 + 1                                        NUD 1010
         COEFF(2,I) = FA*DLAM                                           NUD 1020
         NPROD(3,I) = 20040 + 1                                         NUD 1030
         COEFF(3,I) = FA*DLAM                                           NUD 1040
         M = 3                                                          NUD 1050
         GO TO 110                                                      NUD 1060
C    SET ES-255/FM-255 DECAY                                            NUD 1070
  100    NPROD(1,I) = 982510                                            NUD 1080
         COEFF(1,I) = FA*DLAM                                           NUD 1090
         NPROD(2,I) = 20040                                             NUD 1100
         COEFF(2,I) = FA*DLAM                                           NUD 1110
         M = 2                                                          NUD 1120
         GO TO 110                                                      NUD 1130
  110    CONTINUE                                                       NUD 1140
         IF (MLIP(1).LT.1) GO TO 130                                    NUD 1150
         IF (MOD((L-1),51).NE.0) GO TO 120    !REDUCED 52 TO 51         NUD 1160
         WRITE (IOUT,99999) IOUT, NPAGE                                 NUD 1170
      CALL IDENTIFY (IOUT, 0)                                           PC
         WRITE (LOUT,99998) TITLE                                       NUD 1180
         NPAGE = IPAGE(1)                                               NUD 1190
         WRITE (LOUT,99997)                                             NUD 1200
  120    CALL NOAH(NUCL(I))                                             NUD 1210
         NQZ = NUCAB(I)                                                 NUD 1220
         WRITE (LOUT,99996) NAME, DLAM, FB1, FP, FP1, FT, FA, FSF, FN,  NUD 1230
     *        Q(I), ABUND(NQZ), AMPC(I), WMPC(I)                        NUD 1240
  130    CONTINUE                                                       NUD 1250
         IF (NUCL(I).GT.992540) GO TO 210                               NUD 1260
C       PUT DECAY DATA IN COEFF(M,I)                                    NUD 1270
C    TEST RADIOACTIVITY                                                 NUD 1280
         IF (DLAM.LT.ERR) GO TO 210                                     NUD 1290
         ABETA = 1.0                                                    NUD 1300
C    TEST FOR SPONTANEOUS FISSION                                       NUD 1310
         IF (FSF.LT.ERR) GO TO 140                                      NUD 1320
         M = M + 1                                                      NUD 1330
         COEFF(M,I) = FSF*DLAM                                          NUD 1340
         NPROD(M,I) = 162500                                            NUD 1350
         IF(MOD(NUCL(I),10).NE.0)NPROD(M,I)=NPROD(M,I)+1                NUD 1351
C    *** NUCLID=162500 IS A DUMMY ACTINIDE NUCLID AND IS A              NUD 1360
C    PSUEDO-DAUGHTER                                                    NUD 1370
C        FOR ALL NUCLIDES WHICH DECAY VIA SF.  IT IS STABLE WITH A ZERO NUD 1380
C        CROSS SECTION                                                  NUD 1390
         ABETA = ABETA - FSF                                            NUD 1400
C    TEST POSITRON EMISSION                                             NUD 1410
  140    IF (FP.LT.ERR) GO TO 150                                       NUD 1420
         M = M + 1                                                      NUD 1450
         COEFF(M,I) = FP*DLAM                                           NUD 1460
         NPROD(M,I) = NUCL(I) - 10000                                   NUD 1470
         ABETA = ABETA - FP                                             NUD 1480
C    TEST POSITRON EMISSION TO EXCITED STATE OF PRODUCT NUCLIDE         NUD 1490
  150    IF (FP1.LT.ERR) GO TO 160                                      NUD 1500
         M = M + 1                                                      NUD 1510
         NPROD(M,I) = NUCL(I) - 10000 + 1                               NUD 1520
         COEFF(M,I) = FP1*DLAM                                          NUD 1530
         ABETA = ABETA - FP1                                            NUD 1535
C    TEST ISOMERIC TRANSITION                                           NUD 1540
  160    IF (FT.LT.ERR) GO TO 170                                       NUD 1550
         M = M + 1                                                      NUD 1560
         COEFF(M,I) = FT*DLAM                                           NUD 1570
         NPROD(M,I) = NUCL(I)                                           NUD 1580
         ABETA = ABETA - FT                                             NUD 1590
C    TEST FOR (BETA+NEUTRON) DECAY                                      NUD 1600
  170    IF (FN.LT.ERR) GO TO 180                                       NUD 1610
         M = M + 1                                                      NUD 1620
         COEFF(M,I) = FN*DLAM                                           NUD 1630
         NPROD(M,I) = NUCL(I) + 9990                                    NUD 1640
         ABETA = ABETA - FN                                             NUD 1650
C    TEST ALPHA EMISSION                                                NUD 1660
  180    IF (FA.LT.ERR) GO TO 190                                       NUD 1670
         M = M + 1                                                      NUD 1680
         COEFF(M,I) = FA*DLAM                                           NUD 1690
         NPROD(M,I) = NUCL(I) - 20040                                   NUD 1700
         M = M + 1                                                      NUD 1710
         COEFF(M,I) = COEFF(M-1,I)                                      NUD 1720
         NPROD(M,I) = 20040                                             NUD 1730
         IF (MOD(NUCL(I),10).GT.0) NPROD(M,I) = 20041                   NUD 1740
         ABETA = ABETA - FA                                             NUD 1750
C    TEST NEGATRON EMISSION TO EXCITED STATE OF PRODUCT NUCLIDE         NUD 1760
  190    IF (FB1.LT.ERR) GO TO 200                                      NUD 1770
         M = M + 1                                                      NUD 1780
         COEFF(M,I) = FB1*DLAM                                          NUD 1790
         NPROD(M,I) = NUCL(I) + 10000 + 1                               NUD 1800
         ABETA = ABETA - FB1                                            NUD 1810
C    TEST NEGATRON EMISSION                                             NUD 1820
  200    IF (ABETA.LT.1.E-4) GO TO 210                                  NUD 1830
         M = M + 1                                                      NUD 1840
         COEFF(M,I) = ABETA*DLAM                                        NUD 1850
         NPROD(M,I) = NUCL(I) + 10000                                   NUD 1860
  210    CONTINUE                                                       NUD 1870
         KAP(I) = M                                                     NUD 1880
         GO TO 70                                                       NUD 1890
  220    IF (LITYP.EQ.2) ILITE = I                                      NUD 1900
         IF (LITYP.EQ.3) IACT = I - ILITE                               NUD 1910
         IF (LITYP.EQ.4) IFP = I - IACT - ILITE                         NUD 1920
  230 CONTINUE                                                          NUD 1930
      ITOT = ILITE + IACT + IFP                                         NUD 1950
      DO 240 I=1,ITOT                                                   NUD 1960
         Q(I) = DIS(I)*Q(I)*9.652E+10                                   NUD 1970
         IF (AMPC(I).GT.ERR) AMPC(I) = DIS(I)*1.6283E+13/AMPC(I)        NUD 1980
         IF (WMPC(I).GT.ERR) WMPC(I) = DIS(I)*1.6283E+13/WMPC(I)        NUD 1990
  240 CONTINUE                                                          NUD 2000
      DO 250 IE=1,100                                                   NUD 2010
         IF (CEMTOX(IE).GT.ERR) CEMTOX(IE) = 1.0/CEMTOX(IE)             NUD 2020
  250 CONTINUE                                                          NUD 2030
      NSIZE(19) = IAB                                                   NUD 2040
      RETURN                                                            NUD 2050
  260 WRITE (15,99995) LIHUNT                                           NUD 2060
      RETURN                                                            NUD 2070
99999 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        NUD 2080
99998 FORMAT (1H , 20X, 18A4)                                           NUD 2090
99997 FORMAT (1H , T2, 7HNUCLIDE, T10, 9HDECAY CON, T20, 9HF BETA EX,   NUD 2100
     *     T30, 7HF POSIT, T40, 8HF POS EX, T50, 8HF INT TR, T60,       NUD 2110
     *     7HF ALPHA, T70, 9HF SP FISS, T80, 8HF N+BETA, T90, 7HQ VALUE,NUD 2120
     *     T100, 9HNAT ABUND, T110, 7HAIR RCG, T120, 9HWATER RCG)       NUD 2130
99996 FORMAT (1H , A2, I3, A1, T10, 12(1PE9.3, 1X))                     NUD 2140
99995 FORMAT (41H CORRECT LIBRARY CANNOT BE FOUND; LIHUNT=, I5)         NUD 2150
      END                                                               NUD 2160
      SUBROUTINE NUDAT2(NUCL, Q, FG, GENNEU, AMPC, WMPC, TOCAP, FISS,   NUD   10
     *     DIS, YIELD, COEFF, LC, NPROD, FFSF, FFA, NY, YY, NUCAB,      NUD   20
     *     NYIELD, LOC, NON0, KD, MMAX, KAP, AP, ABUND, NNSORS, LPUCNT, NUD   30
     *     ICNMAX)                                                      NUD   40
      INTEGER*2       NYIELD,LOC,NON0,KD,MMAX,KAP                       NUD   50
      DOUBLE PRECISION DNUM, UNIT, TENS                                 NUD   60
      INTEGER ELE, STA                                                  NUD   70
      INTEGER AREAD, BASIS, TITLE                                       NUD   80
      DIMENSION NUCAL(10), NKFUD(10)                                    NUD   90
      DIMENSION TYLD(11), NEWCX(11), NSRS(11), CAPT(11), TITLE(18)      NUD  100
      DIMENSION NUCL(1), Q(1), FG(1), GENNEU(1), AMPC(1), WMPC(1),      NUD  110
     *     TOCAP(1), A(1), FISS(1), DIS(1), YIELD(1), COEFF(LC,1),      NUD  120
     *     NPROD(LC,1), FFSF(1), FFA(1), NY(1), YY(1), NUCAB(1),        NUD  130
     *     NYIELD(1), LOC(1), NON0(1), KD(1), MMAX(1), KAP(1), AP(1),   NUD  140
     *     ABUND(1)                                                     NUD  150
      COMMON /FUDG/ NFUD(7,20), FUDSIG(32,20), CBUP(32), NXSEC, NSM,    NUD  160
     *     NPT, NCO(2,20), FPBST(050,3), VECT(050), FOLD(20), FNEW(20)  NUD  170
      COMMON /INST/ INSTR(300), NINSTR(300), IFIND(300), NNTYP(32),     NUD  180
     *     LAB(300)                                                     NUD  190
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  NUD  200
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           NUD  210
     *     TFPBAV(4), IPHMAX                                            NUD  220
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 NUD  230
      COMMON /MANE4I/ NINP(15,06), NIRR(150,04), NPRO(20,05),           NUD  240
     *     NCON(20,01), NDOL(02,04), NMOV(99,03), NADD(30,03),          NUD  250
     *     NHED(50,01), NWAC(09,02), NLPU(09,100), NPCH(54,03),         NUD  260
     *     NBUP(20), INUCL(5), NKEQ(03,05), NFAC(20,4), NGTO(10,2),     NUD  270
     *     NOFF(5)                                                              
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            NUD  280
     *     INCYLD(8,4), MYLD(4)                                         NUD  290
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    NUD  300
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  NUD  310
     *     ITITLE(18)                                                   NUD  320
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      NUD  330
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           NUD  340
      COMMON /SIGR/ DDUM(1), SIGNG, SIGN2N, SIGN3N, SIGF, SNG1, SN2N1,  NUD  350
     *     Y(08), JIBE, NYLD, NCL, NCSK                                 NUD  360
      COMMON /SZ/ NSIZE(25)                                             NUD  370
C      LEVEL 2,NUCL,Q,FG,GENNEU,AMPC,WMPC,TOCAP,FISS,DIS,YIELD,COEFF,   NUD  380
C    1  NPROD,FFSF,NY,YY,NUCAB,NYIELD,LOC,NON0,KD,MMAX,KAP,AP,ABUND     NUD  390
      DATA NUCAL /-9,-20,10,11,-10,-20030,-10000,3*0/                   NUD  400
      DATA NKFUD /6,3,1,2,5,7,4,3*0/                                    NUD  410
C    1=N,G   2=N,G EX   3=N,3N    4=N,FISS   5=N,2N   6=N,2N EX         NUD  420
C    7=N,ALPHA                                                          NUD  430
      DATA NSRS /902320,922330,922350,922380,942390,942410,962450,      NUD  440
     *     982520,3*0/                                                  NUD  450
      NSM = 0                                                           NUD  460
C    FIND ACTINIDES WHICH HAVE FISSION PRODUCT YIELDS                   NUD  470
      IF (NXSEC.LE.0 .OR. NNSORS.LE.0) GO TO 30                         NUD  480
      DO 20 KK=1,NXSEC                                                  NUD  490
         DO 10 K=1,NNSORS                                               NUD  500
            J = INCYLD(K,MLIB(11))                                      NUD  510
            IF (NSRS(J).NE.NFUD(1,KK) .OR. NFUD(2,KK).NE.4) GO TO 10    NUD  520
            NSM = NSM + 1                                               NUD  530
            NFUD(6,KK) = NSM                                            NUD  540
   10    CONTINUE                                                       NUD  550
   20 CONTINUE                                                          NUD  560
   30 CONTINUE                                                          NUD  570
      IF (NSM.GT.ICNMAX) WRITE (15,99999) NSM, ICNMAX                   NUD  580
      IDM(1) = 0                                                        NUD  590
      DO 40 K=1,7                                                       NUD  600
         CAPT(K) = 0.0                                                  NUD  610
   40 CONTINUE                                                          NUD  620
      NYTOT = 2                                                         NUD  630
      YIELD(1) = 0.0                                                    NUD  640
      NYLD = 8                                                          NUD  650
      DO 50 I=1,IFP                                                     NUD  660
         NYIELD(I) = 1                                                  NUD  670
   50 CONTINUE                                                          NUD  680
      DO 60 I=1,IACT                                                    NUD  690
         GENNEU(I) = 0.0                                                NUD  700
         FISS(I) = 0.0                                                  NUD  710
   60 CONTINUE                                                          NUD  720
      DO 70 I=1,ITOT                                                    NUD  730
         TOCAP(I) = 0.0                                                 NUD  740
         MMAX(I) = KAP(I)                                               NUD  750
   70 CONTINUE                                                          NUD  760
      DO 80 I=1,NYLD                                                    NUD  770
         TYLD(I) = 0.0                                                  NUD  780
   80 CONTINUE                                                          NUD  790
      LPUN = 1                                                          NUD  800
      LOUT = IOUT                                                       NUD  810
C    READ CROSS-SECTION DATA FOR ALL NUCLIDES                           NUD  820
      DO 360 LITYP=5,7                                                  NUD  830
         IF (MLIB(LITYP).EQ.0) GO TO 360                                NUD  840
         IF (MLIB(LITYP).GT.0) GO TO 110                                NUD  850
         LPUN = 1                                                       NUD  860
         LPUCNT = 0                                                     NUD  870
         DO 90 MLVAR=2,LITYP                                            NUD  880
            IF (MLIB(MLVAR).LT.0) LPUCNT = LPUCNT + 1                   NUD  890
   90    CONTINUE                                                       NUD  900
         DO 100 J=1,100                                                 NUD  910
            MLPU(J) = NLPU(LPUCNT,J)                                    NUD  920
  100    CONTINUE                                                       NUD  930
  110    CONTINUE                                                       NUD  940
C    FIND CORRECT LIBRARY                                               NUD  950
         LIHUNT = IABS(MLIB(LITYP))                                     NUD  960
         DO 120 J=1,50000                                               NUD  970
            LIBE = IREAD(1)                                             NUD  980
            IF (LIBE.EQ.LIHUNT) GO TO 130                               NUD  990
  120    CONTINUE                                                       NUD 1000
         GO TO 430                                                      NUD 1010
  130    CONTINUE                                                       NUD 1020
         ICHR = 8                                                       NUD 1030
         DO 140 ID=1,18                                                 NUD 1040
            TITLE(ID) = AREAD(0)                                        NUD 1050
  140    CONTINUE                                                       NUD 1060
         IF (MLIP(2).GT.0) CALL TOC(NPAGE, 03, 02, 01, 00, TITLE, 72,   NUD 1070
     *        NTOC, NREC)                                               NUD 1080
         I0 = 0                                                         NUD 1090
         IG = IDM(LITYP-3)                                              NUD 1100
         KMZ = LITYP - 4                                                NUD 1110
         DO 150 K=1,KMZ                                                 NUD 1120
            I0 = I0 + IDM(K)                                            NUD 1130
  150    CONTINUE                                                       NUD 1140
         JJ = I0                                                        NUD 1150
         I1 = I0 + IG                                                   NUD 1160
         I0 = I0 + 1                                                    NUD 1170
         DO 350 IP=1,IG                                                 NUD 1180
            JJ = JJ + 1                                                 NUD 1190
  160       CALL SIGRED(JJ, I, I0, I1, LITYP, NUCL)                     NUD 1200
            IF (I.EQ.-1) GO TO 360                                      NUD 1210
            IF (I.NE.-2) GO TO 170                                      NUD 1220
            WRITE (IOUT,99982) JJ                                       NUD 1230
            GO TO 160                                                   NUD 1240
  170       CONTINUE                                                    NUD 1250
            IF (MLIB(LITYP).GT.0) GO TO 180                             NUD 1260
            IF (MLPU(LPUN).NE.NUCL(I)) GO TO 180                        NUD 1270
            INPT = MLIB(9)                                              NUD 1280
            CALL SIGRED(JJ, I, I0, I1, LITYP, NUCL)                     NUD 1290
            IF (I.EQ.-2) GO TO 440                                      NUD 1300
            LPUN = LPUN + 1                                             NUD 1310
            INPT = MLIB(8)                                              NUD 1320
  180       IF (MLIP(2).LT.1) GO TO 210                                 NUD 1330
            IF (MOD((IP-1),51).NE.0) GO TO 190   !REDUCED 52 TO 51      NUD 1340
            WRITE (IOUT,99998) IOUT, NPAGE                              NUD 1350
      CALL IDENTIFY (IOUT, 0)                                           PC
            WRITE (LOUT,99997) TITLE                                    NUD 1360
            NPAGE = IPAGE(1)                                            NUD 1370
            IF (LITYP.EQ.5) WRITE (LOUT,99994)                          NUD 1380
            IF (LITYP.EQ.6) WRITE (LOUT,99996)                          NUD 1390
            IF (LITYP.EQ.7) WRITE (LOUT,99995)                          NUD 1400
            IF (LITYP.EQ.7 .AND. MLIB(11).NE.4) WRITE (IOUT,99993)      NUD 1410
            IF (LITYP.EQ.7 .AND. MLIB(11).EQ.4) WRITE (IOUT,99992)      NUD 1420
  190       CALL NOAH(NUCL(I))                                          NUD 1430
            WRITE (IOUT,99987) NAME, SIGNG, SIGN2N, SIGN3N, SIGF, SNG1, NUD 1440
     *           SN2N1                                                  NUD 1450
            IF (LITYP.NE.7) GO TO 200                                   NUD 1460
            IF (NNSORS.LE.6) WRITE (IOUT,99991) (Y(J),J=1,6)            NUD 1470
            IF (NNSORS.EQ.8) WRITE (IOUT,99990) (Y(J),J=1,8)            NUD 1480
  200       CONTINUE                                                    NUD 1490
  210       CONTINUE                                                    NUD 1500
            M = KAP(I)                                                  NUD 1510
            IF (LITYP.NE.6) CAPT(6) = SIGN3N                            NUD 1520
            IF (LITYP.NE.6) CAPT(7) = SIGF                              NUD 1530
            CAPT(4) = SNG1                                              NUD 1540
            CAPT(3) = SIGNG                                             NUD 1550
            CAPT(1) = SN2N1                                             NUD 1560
            CAPT(5) = SIGN2N                                            NUD 1570
            IF (LITYP.EQ.6) CAPT(2) = SIGN3N                            NUD 1580
            IF (LITYP.EQ.6) FISS(I-ILITE) = SIGF                        NUD 1590
C    PUT CROSS-SECTIONS IN COEFF(M,I) ARRAY                             NUD 1600
            DO 260 K=1,7                                                NUD 1610
               IF (CAPT(K).LT.ERR) GO TO 250                            NUD 1620
               M = M + 1                                                NUD 1630
               NPROD(M,I) = NUCL(I) + NUCAL(K)                          NUD 1640
               COEFF(M,I) = CAPT(K)                                     NUD 1650
               TOCAP(I) = TOCAP(I) + CAPT(K)                            NUD 1660
C    STORE INITIAL XSEC AND LOCATION IN COEFF FOR EACH VARIALBE XSEC;   NUD 1670
C    FISSION XSEC NOT STORED IN COEFF, ONLY IN TOCAP AND FISSION        NUD 1680
C    PRODUCT YIELDS                                                     NUD 1690
               IF (NXSEC.EQ.0) GO TO 230                                NUD 1700
               DO 220 NN=1,NXSEC                                        NUD 1710
                  IF (NUCL(I).NE.NFUD(1,NN)) GO TO 220                  NUD 1720
                  NFUD(4,NN) = I                                        NUD 1730
                  IF (NKFUD(K).EQ.NFUD(2,NN)) NFUD(7,NN) = M            NUD 1740
                  IF (NFUD(2,NN).EQ.4) FNEW(NN) = SIGF                  NUD 1750
                  IF (NFUD(2,NN).EQ.2 .AND. K.EQ.4) FNEW(NN) = SNG1     NUD 1760
                  IF (NFUD(2,NN).EQ.1 .AND. K.EQ.3) FNEW(NN) = SIGNG    NUD 1770
  220          CONTINUE                                                 NUD 1780
  230          CONTINUE                                                 NUD 1790
               IF (K.NE.6 .OR. M.GE.LC) GO TO 240                       NUD 1800
               IF (M .EQ.(LC-1) .AND. CAPT(7).GE.ERR) GO TO 250         NUD 1810
               M = M + 1                                                NUD 1820
               NPROD(M,I) = 20040                                       NUD 1830
               COEFF(M,I) = COEFF(M-1,I)                                NUD 1840
  240          IF (K.NE.7 .OR. M.GE.LC) GO TO 250                       NUD 1850
               M = M + 1                                                NUD 1860
               COEFF(M,I) = COEFF(M-1,I)                                NUD 1870
               NPROD(M,I) = 10010                                       NUD 1880
  250          CAPT(K) = 0.0                                            NUD 1890
  260       CONTINUE                                                    NUD 1900
            MMAX(I) = M                                                 NUD 1910
            IF (LITYP.EQ.6) TOCAP(I) = TOCAP(I) + SIGF                  NUD 1920
            IF (M.GT.LC) WRITE (15,99985) M, LC, NUCL(I),               NUD 1930
     *           (NPROD(MM,I),MM=1,M)                                   NUD 1940
            IF (M.GT.LC) MMAX(I) = LC                                   NUD 1950
C    CALCULATE NEUTRON PRODUCTION RATE                                  NUD 1960
C       MATCH INPUT ACTINIDE NEUTRON YIELDS WITH NUCLIDE ID             NUD 1970
            IF (LITYP.NE.6) GO TO 290                                   NUD 1980
            YNEU = 0.0                                                  NUD 1990
            DO 270 II=1,IAMAX                                           NUD 2000
               IF (NY(II).LT.1) GO TO 290                               NUD 2010
               IF (NY(II).EQ.NUCL(I)) GO TO 280                         NUD 2020
  270       CONTINUE                                                    NUD 2030
  280       YNEU = YY(II)                                               NUD 2040
            GENNEU(I-ILITE) = SIGF*YNEU + SIGN2N*2.0 + SIGN3N*3.0       NUD 2050
  290       CONTINUE                                                    NUD 2060
C       WRITE(LOUT,9026) NAME,GENNEU(IP),YNEU                           NUD 2070
C    PLACE FISSION PRODUCT YIELDS IN YIELD(J,IL)                        NUD 2080
            IF (LITYP.NE.7) GO TO 340                                   NUD 2090
            IF (MLIB(11).EQ.0) GO TO 340                                NUD 2100
            YTOT = 0.0                                                  NUD 2110
            IE = IACT + ILITE                                           NUD 2120
            DO 300 KJ=1,NNSORS                                          NUD 2130
               J = INCYLD(KJ,MLIB(11))                                  NUD 2140
               YTOT = YTOT + Y(J)                                       NUD 2150
  300       CONTINUE                                                    NUD 2160
            NYIELD(I-IE) = 1                                            NUD 2170
            IF (YTOT.LT.ERR) GO TO 320                                  NUD 2180
            NYIELD(I-IE) = NYTOT                                        NUD 2190
            DO 310 KJ=1,NNSORS                                          NUD 2200
               J = INCYLD(KJ,MLIB(11))                                  NUD 2210
               YJ = Y(J)*0.01                                           NUD 2220
               Y(J) = 0.0                                               NUD 2230
               TYLD(J) = TYLD(J) + YJ                                   NUD 2240
               YIELD(NYTOT) = YJ                                        NUD 2250
               NYTOT = NYTOT + 1                                        NUD 2260
  310       CONTINUE                                                    NUD 2270
  320       DO 330 KJ=1,NYLD                                            NUD 2280
               Y(KJ) = 0.0                                              NUD 2290
  330       CONTINUE                                                    NUD 2300
  340       CONTINUE                                                    NUD 2310
  350    CONTINUE                                                       NUD 2320
  360 CONTINUE                                                          NUD 2330
      IF (IFP.LT.1) GO TO 380                                           NUD 2340
      IF (NNSORS.EQ.8) WRITE (IOUT,99989) (TYLD(J),J=1,8)               NUD 2350
      IF (NNSORS.LT.8) WRITE (IOUT,99988) (TYLD(J),J=1,6)               NUD 2360
      ML = MLIB(11)                                                     NUD 2370
      DO 370 J=1,NNSORS                                                 NUD 2380
         INC = INCYLD(J,ML)                                             NUD 2390
         NEWCX(J) = NSRS(INC)                                           NUD 2400
  370 CONTINUE                                                          NUD 2410
      WRITE (IOUT,99986) (NEWCX(J),J=1,NNSORS)                          NUD 2420
  380 CONTINUE                                                          NUD 2430
      NSIZE(14) = NYTOT                                                 NUD 2440
      IF (MLIB(10).EQ.0) GO TO 420                                      NUD 2450
C    ** READ ODD REACTIONS**                                            NUD 2460
      INPT = MLIB(8)                                                    NUD 2470
      IF (MLIB(10).LT.0) INPT = MLIB(9)                                 NUD 2480
      KOUN = IABS(MLIB(10))                                             NUD 2490
      DO 410 K=1,KOUN                                                   NUD 2500
         NINT = IREAD(1)                                                NUD 2510
         NFIN = IREAD(0)                                                NUD 2520
         SIG = READ(0)                                                  NUD 2530
         DO 400 II=1,ITOT                                               NUD 2540
            IF (NINT.NE.NUCL(II)) GO TO 400                             NUD 2550
            WRITE (LOUT,99981) NINT, NFIN, SIG                          NUD 2560
            MMAX(II) = MMAX(II) + 1                                     NUD 2570
            LM = MMAX(II)                                               NUD 2580
            IF (LM.LE.LC) GO TO 390                                     NUD 2590
            WRITE (15,99985) LM, LC, NUCL(II), (NPROD(NM,II),NM=1,LM)   NUD 2600
            MMAX(II) = LC                                               NUD 2610
  390       CONTINUE                                                    NUD 2620
            COEFF(MMAX(II),II) = SIG                                    NUD 2630
            NPROD(MMAX(II),II) = NFIN                                   NUD 2640
      TOCAP(II)=TOCAP(II)+SIG                                                   
  400    CONTINUE                                                       NUD 2660
  410 CONTINUE                                                          NUD 2670
      INPT = MLIB(8)                                                    NUD 2680
  420 CONTINUE                                                          NUD 2690
      IF (IFP.GT.IFMAX) WRITE (15,99984) IFP, IFMAX                     NUD 2700
      IF (ITOT.GT.ITMAX) WRITE (15,99983) ITOT, ITMAX                   NUD 2710
      RETURN                                                            NUD 2720
  430 WRITE (15,99980) LIHUNT                                           NUD 2730
      XDIV = 0.0                                                        NUD 2740
      XDIV = 1.0/XDIV                                                   NUD 2750
      RETURN                                                            NUD 2760
  440 WRITE (15,99982) JJ                                               NUD 2770
      RETURN                                                            NUD 2780
99999 FORMAT (54H **WARNING 9609** NUMBER OF ACTINIDES DIRECTLY CONNECT,NUD 2790
     *     35HED TO FISS PROD AND WITH VAR XSEC =, I5, 13H IS GREATER T,NUD 2800
     *     15HHAN ALLOWABLE =, I5)                                      NUD 2810
99998 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        NUD 2820
99997 FORMAT (1H , 20X, 18A4)                                           NUD 2830
99996 FORMAT (1H , T2, 7HNUCLIDE, T10, 9H(N,GAMMA), T20, 7H (N,2N),     NUD 2840
     *     T30, 7H (N,3N), T40, 8H(N,FISS), T49, 8H (N,G)EX, T57,       NUD 2850
     *     9H (N,2N)EX)                                                 NUD 2860
99995 FORMAT (1H , T2, 7HNUCLIDE, T10, 9H(N,GAMMA), T20, 7H (N,2N),     NUD 2870
     *     T30, 9H(N,ALPHA), T40, 8H(N,PROT), T49, 8H (N,G)EX, T57,     NUD 2880
     *     8H(N,2N)EX, T66, 18(1H*), 22HFISSION PRODUCT YIELDS, 18(1H*))NUD 2890
99994 FORMAT (1H , T2, 7HNUCLIDE, T10, 9H(N,GAMMA), T20, 7H (N,2N),     NUD 2900
     *     T30, 9H(N,ALPHA), T40, 8H(N,PROT), T49, 8H (N,G)EX, T57,     NUD 2910
     *     8H(N,2N)EX)                                                  NUD 2920
99993 FORMAT (1H , T70, 5HTH232, T81, 5H U233, T92, 5H U235, T103,      NUD 2930
     *     5H U238, T114, 5HPU239, T125, 5HPU241)                       NUD 2940
99992 FORMAT (1H , T66, 5HTH232, T74, 5H U233, T82, 5H U235, T90,       NUD 2950
     *     5H U238, T98, 5HPU239, T106, 5HPU241, T114, 5HCM245, T122,   NUD 2960
     *     5HCF249)                                                     NUD 2970
99991 FORMAT (1H+, 65X,1P,6E11.4)                                       NUD 2980
99990 FORMAT (1H+, 65X, 8F8.5)                                          NUD 2990
99989 FORMAT (40H0SUM OF YIELDS FOR EACH FISSILE ISOTOPE:, T66,         NUD 3000
     *     8(2PF7.3, 1X))                                               NUD 3010
99988 FORMAT (40H0SUM OF YIELDS FOR EACH FISSILE ISOTOPE:, T67,         NUD 3020
     *    1P,6E11.4)                                                    NUD 3030
99987 FORMAT (1H , A2, I3, A1, 1X,1P,4E10.3,1P,2E8.1,0P,8F8.5)          NUD 3040
99986 FORMAT (54H0ISOTOPES FOR WHICH YIELDS ARE EXPLICITLY ACCOUNTED FO,NUD 3050
     *     2HR:, 8(1X, I6))                                             NUD 3060
99985 FORMAT (41H WARNING: M IS OUT OF RANGE IN NUDAT2. M=, I2, 3HLC=,  NUD 3070
     *     I2, 10I8)                                                    NUD 3080
99984 FORMAT (19H WARNING: # OF FP =, I4, 22H AND EXCEEDS LIMIT OF , I4)NUD 3090
99983 FORMAT (25H WARNING: # OF NUCLIDES =, I4, 20H AND EXCEEDS LIMIT O,NUD 3100
     *     2HF , I4)                                                    NUD 3110
99982 FORMAT (43H NUCLIDE ID NOT FOUND IN DECAY LIBRARY. JJ=, I8)       NUD 3120
99981 FORMAT (18H ODD RXN; PARENT =, I8, 3X, 10HDAUGHTER =, I8, 3X,     NUD 3130
     *     14HCROSS SECTION=, 1PE12.4, 5HBARNS)                         NUD 3140
99980 FORMAT (41H CORRECT LIBRARY CANNOT BE FOUND; LIHUNT=, I5)         NUD 3150
      END                                                               NUD 3160
      SUBROUTINE NUDAT3(NUCL, Q, FG, GENNEU, AMPC, WMPC, TOCAP, A,      NUD   10
     *     FISS, DIS, YIELD, COEFF, LC, NPROD, FFSF, FFA, NY, YY,       NUD   20
     *     NUCAB, NYIELD, LOC, NON0, KD, MMAX, KAP, AP, ABUND, NNSORS,  NUD   30
     *     NFUDFP, IFD)                                                 NUD   40
      INTEGER*2       NYIELD,LOC,NON0,KD,MMAX,KAP ,NFUDFP(IFD,1)        NUD   50
      DOUBLE PRECISION DNUM, UNIT, TENS                                 NUD   60
      INTEGER ELE, STA                                                  NUD   70
      INTEGER TITLE, TDUM                                               NUD   80
      DIMENSION NSORS(08),IDUM(1)                                       NUD   90
      DIMENSION TYLD(11), NEWCX(11), NSRS(11), CAPT(11), TITLE(18)      NUD  100
      DIMENSION NUCL(1), Q(1), FG(1), GENNEU(1), AMPC(1), WMPC(1),      NUD  110
     *     TOCAP(1), A(1), FISS(1), DIS(1), YIELD(1), COEFF(LC,1),      NUD  120
     *     NPROD(LC,1), FFSF(1), FFA(1), NY(1), YY(1), NUCAB(1),        NUD  130
     *     NYIELD(1), LOC(1), NON0(1), KD(1), MMAX(1), KAP(1), AP(1),   NUD  140
     *     ABUND(1)                                                     NUD  150
      COMMON /FUDG/ NFUD(7,20), FUDSIG(32,20), CBUP(32), NXSEC, NSM,    NUD  160
     *     NPT, NCO(2,20), FPBST(050,3), VECT(050), FOLD(20), FNEW(20)  NUD  170
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  NUD  180
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           NUD  190
     *     TFPBAV(4), IPHMAX                                            NUD  200
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            NUD  210
     *     INCYLD(8,4), MYLD(4)                                         NUD  220
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    NUD  230
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  NUD  240
     *     ITITLE(18)                                                   NUD  250
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      NUD  260
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           NUD  270
      COMMON /SZ/ NSIZE(25)                                             NUD  280
      EQUIVALENCE (IDUM(1),TDUM)                                        NUD  281
C      LEVEL 2,NUCL,Q,FG,GENNEU,AMPC,WMPC,TOCAP,FISS,DIS,YIELD,COEFF,   NUD  290
C    1  NPROD,FFSF,NY,YY,NUCAB,NYIELD,LOC,NON0,KD,MMAX,KAP,AP,ABUND     NUD  300
C    2  ,A,NFUDFP                                                       NUD  310
      DATA IDUM /4H    /                                                NUD  320
      DATA NSORS /902320,922330,922350,922380,942390,942410,962450,     NUD  330
     *     982520/                                                      NUD  340
      IF (NXSEC.LT.1) NSM = 3                                           NUD  350
      LOUT = IOUT                                                       NUD  360
C                                                                       NUD  370
C    ALL DATA ON NUCLIDES HAS BEEN READ, BEGIN TO COMPUTE MATRIX COEFF  NUD  380
C                                                                       NUD  390
C                                                                       NUD  400
C    FIND PRODUCT NUCLIDES FOR REACTIONS OF LIGHT ELEMENTS(已经读取了所有核素的数据，开始计算矩阵系数。接下来的部分是查找轻元素反应产生的产物核素。)              NUD  410
C                                                                       NUD  420
      NON = 0                                                           NUD  430
      DO 30 I=1,ITOT                                                    NUD  440
         IF (MOD(NUCL(I),10).EQ.0) GO TO 20                             NUD  450
         M = MMAX(I)                                                    NUD  460
         DO 10 K=1,M                                                    NUD  470
            NPROD(K,I) = NPROD(K,I) - 1                                 NUD  480
   10    CONTINUE                                                       NUD  490
   20    NON0(I) = 0                                                    NUD  500
   30 CONTINUE                                                          NUD  510
      IF (ILITE.LT.1) GO TO 90                                          NUD  520
      DO 80 I=1,ILITE                                                   NUD  530
         NUCLI = NUCL(I)                                                NUD  540
         DO 50 J=1,ILITE                                                NUD  550
            KMAX = KAP(J)                                               NUD  560
            IF (KMAX.LT.1) GO TO 50                                     NUD  570
            DO 40 M=1,KMAX                                              NUD  580
               IF (NUCLI.NE.NPROD(M,J)) GO TO 40                        NUD  590
               NON0(I) = NON0(I) + 1                                    NUD  600
               NON = NON + 1                                            NUD  610
               IF (NON.GT.IZMAX) GO TO 410                              NUD  620
               A(NON) = COEFF(M,J)                                      NUD  630
               JT = J                                                   NUD  640
               LOC(NON) = JT                                            NUD  650
   40       CONTINUE                                                    NUD  660
   50    CONTINUE                                                       NUD  670
         KD(I) = NON0(I)                                                NUD  680
         DO 70 J=1,ILITE                                                NUD  690
            K1 = KAP(J) + 1                                             NUD  700
            KMAX = MMAX(J)                                              NUD  710
            IF (KMAX.LT.K1) GO TO 70                                    NUD  720
            DO 60 M=K1,KMAX                                             NUD  730
               IF (NUCLI.NE.NPROD(M,J)) GO TO 60                        NUD  740
               NON0(I) = NON0(I) + 1                                    NUD  750
               NON = NON + 1                                            NUD  760
               IF (NON.GT.IZMAX) GO TO 410                              NUD  770
               A(NON) = COEFF(M,J)                                      NUD  780
               JT = J                                                   NUD  790
               LOC(NON) = JT                                            NUD  800
   60       CONTINUE                                                    NUD  810
   70    CONTINUE                                                       NUD  820
   80 CONTINUE                                                          NUD  830
C                                                                       NUD  840
C    NON ZERO MATRIX ELEMENTS FOR THE ACTINIDES                         NUD  850
C                                                                       NUD  860
   90 IF (IACT.LT.1) GO TO 160                                          NUD  870
      I0 = ILITE + 1                                                    NUD  880
      I1 = ILITE + IACT                                                 NUD  890
      IF (I1.LT.I0) GO TO 160                                           NUD  900
      DO 150 I=I0,I1                                                    NUD  910
         NUCLI = NUCL(I)                                                NUD  920
         DO 110 J=I0,I1                                                 NUD  930
            MAX = KAP(J)                                                NUD  940
            IF (MAX.LT.1) GO TO 110                                     NUD  950
            DO 100 M=1,MAX                                              NUD  960
               IF (NUCLI.NE.NPROD(M,J)) GO TO 100                       NUD  970
               NON0(I) = NON0(I) + 1                                    NUD  980
               NON = NON + 1                                            NUD  990
               IF (NON.GT.IZMAX) GO TO 410                              NUD 1000
               A(NON) = COEFF(M,J)                                      NUD 1010
               JT = J                                                   NUD 1020
               LOC(NON) = JT                                            NUD 1030
  100       CONTINUE                                                    NUD 1040
  110    CONTINUE                                                       NUD 1050
         KD(I) = NON0(I)                                                NUD 1060
         DO 140 J=I0,I1                                                 NUD 1070
            M1 = KAP(J) + 1                                             NUD 1080
            M2 = MMAX(J)                                                NUD 1090
            IF (M2.LT.M1) GO TO 140                                     NUD 1100
            DO 130 M=M1,M2                                              NUD 1110
               IF (NUCLI.NE.NPROD(M,J)) GO TO 130                       NUD 1120
               NON0(I) = NON0(I) + 1                                    NUD 1130
               NON = NON + 1                                            NUD 1140
               IF (NON.GT.IZMAX) GO TO 410                              NUD 1150
C    STORE LOCATION OF VARIABLE XSEC IN ARRAY A; NOT APPLICABLE TO      NUD 1160
C    FISSION CROSS SECTIONS                                             NUD 1170
               DO 120 KK=1,NXSEC                                        NUD 1180
                  IF (M.NE.NFUD(7,KK) .OR. J.NE.NFUD(4,KK)) GO TO 120   NUD 1190
                  IF (NFUD(2,KK).EQ.0) GO TO 120                        NUD 1200
                  NFUD(5,KK) = NON                                      NUD 1210
  120          CONTINUE                                                 NUD 1220
               A(NON) = COEFF(M,J)                                      NUD 1230
               JT = J                                                   NUD 1240
               LOC(NON) = JT                                            NUD 1250
  130       CONTINUE                                                    NUD 1260
  140    CONTINUE                                                       NUD 1270
  150 CONTINUE                                                          NUD 1280
C                                                                       NUD 1290
C    MATRIX ELEMENTS FOR FISSION PRODUCTS                               NUD 1300
C                                                                       NUD 1310
  160 IF (IFP.LT.1) GO TO 300                                           NUD 1320
C    DETERMINE IF ACTINIDE (I) HAS A FISSION YIELD                      NUD 1330
      IF (IACT.LT.1) GO TO 190                                          NUD 1340
      IF (MLIB(11).EQ.0) GO TO 190                                      NUD 1350
      DO 180 I=I0,I1                                                    NUD 1360
         DO 170 K=1,NNSORS                                              NUD 1370
            J = INCYLD(K,MLIB(11))                                      NUD 1380
            IF (NUCL(I).EQ.NSORS(J)) NSORS(J) = I                       NUD 1390
  170    CONTINUE                                                       NUD 1400
  180 CONTINUE                                                          NUD 1410
  190 CONTINUE                                                          NUD 1420
      IM = ILITE + IACT                                                 NUD 1430
      I0 = IM + 1                                                       NUD 1440
      IF (ITOT.LT.I0) GO TO 300                                         NUD 1450
      DO 290 I=I0,ITOT                                                  NUD 1460
         NUCLI = NUCL(I)                                                NUD 1470
         DO 210 J=I0,ITOT                                               NUD 1480
            KMAX = KAP(J)                                               NUD 1490
            IF (KMAX.LT.1) GO TO 210                                    NUD 1500
            DO 200 M=1,KMAX                                             NUD 1510
               IF (NUCLI.NE.NPROD(M,J)) GO TO 200                       NUD 1520
               NON0(I) = NON0(I) + 1                                    NUD 1530
               NON = NON + 1                                            NUD 1540
               IF (NON.GT.IZMAX) GO TO 410                              NUD 1550
               A(NON) = COEFF(M,J)                                      NUD 1560
               JT = J                                                   NUD 1570
               LOC(NON) = JT                                            NUD 1580
  200       CONTINUE                                                    NUD 1590
  210    CONTINUE                                                       NUD 1600
         KD(I) = NON0(I)                                                NUD 1610
         DO 230 J=I0,ITOT                                               NUD 1620
            K1 = KAP(J) + 1                                             NUD 1630
            KMAX = MMAX(J)                                              NUD 1640
            IF (KMAX.LT.K1) GO TO 230                                   NUD 1650
            DO 220 M=K1,KMAX                                            NUD 1660
               IF (NUCLI.NE.NPROD(M,J)) GO TO 220                       NUD 1670
               NON0(I) = NON0(I) + 1                                    NUD 1680
               NON = NON + 1                                            NUD 1690
               IF (NON.GT.IZMAX) GO TO 410                              NUD 1700
               A(NON) = COEFF(M,J)                                      NUD 1710
               JT = J                                                   NUD 1720
               LOC(NON) = JT                                            NUD 1730
  220       CONTINUE                                                    NUD 1740
  230    CONTINUE                                                       NUD 1750
         IF (IACT.LT.1) GO TO 290                                       NUD 1760
         IF (MLIB(11).EQ.0) GO TO 290                                   NUD 1770
         IL = I - IM                                                    NUD 1780
         DO 240 KK=1,NSM                                                NUD 1790
            NFUDFP(IL,KK) = 0                                           NUD 1800
  240    CONTINUE                                                       NUD 1810
         NYNOW = NYIELD(IL) - 1                                         NUD 1820
         IF (NYNOW.LT.1) GO TO 290                                      NUD 1830
         DO 280 K=1,NNSORS                                              NUD 1840
            IF (YIELD(NYNOW+K).LT.ERR) GO TO 280                        NUD 1850
            J = INCYLD(K,MLIB(11))                                      NUD 1860
            NON = NON + 1                                               NUD 1870
            IF (NON.GT.IZMAX) GO TO 410                                 NUD 1880
C    STORE LOCATIONS OF FP YIELDS IN ARRAY A WHICH ARE CONNECTED TO     NUD 1890
C    VARIABLE FISSION XSECS                                             NUD 1900
C    SET DIRECT FP YIELDS WHEN NO VARIABLE XSEC ARE USED                NUD 1910
            IF (NXSEC.GT.0) GO TO 250                                   NUD 1920
            MUCL = NUCL(NSORS(J))                                       NUD 1930
            IF (MUCL.EQ.922330) NFUDFP(IL,1) = NON                      NUD 1940
            IF (MUCL.EQ.922350) NFUDFP(IL,2) = NON                      NUD 1950
            IF (MUCL.EQ.942390) NFUDFP(IL,3) = NON                      NUD 1960
            GO TO 270                                                   NUD 1970
  250       CONTINUE                                                    NUD 1980
            DO 260 KK=1,NXSEC                                           NUD 1990
               IF (NUCL(NSORS(J)).EQ.NFUD(1,KK) .AND. NFUD(2,KK).EQ.4)  NUD 2000
     *              NFUDFP(IL,NFUD(6,KK)) = NON                         NUD 2010
  260       CONTINUE                                                    NUD 2020
  270       CONTINUE                                                    NUD 2030
            NON0(I) = NON0(I) + 1                                       NUD 2040
            KK = NSORS(J)                                               NUD 2050
            LOC(NON) = KK                                               NUD 2060
            KF = KK - ILITE                                             NUD 2070
            A(NON) = YIELD(NYNOW+K)*FISS(KF)                            NUD 2080
  280    CONTINUE                                                       NUD 2090
  290 CONTINUE                                                          NUD 2100
      IF (IFP.LE.0) GO TO 300                                           NUD 2110
C                                                                       NUD 2120
C    ALL MATRIX ELEMENTS ARE NOW COMPUTED                               NUD 2130
C    BEGIN TRANSIENT SOLUTION(所有矩阵元素已经计算完成，现在开始进行矩阵的瞬态解决方案)                                           NUD 2140
C                                                                       NUD 2150
C                                                                       NUD 2160
C    TEMPORARILY WRITE OUT MATRIX ELEMENTS                              NUD 2170
C                                                                       NUD 2180
  300 CONTINUE                                                          NUD 2190
      DO 310 IV=1,18                                                    NUD 2200
         TITLE(IV) = TDUM                                               NUD 2210
  310 CONTINUE                                                          NUD 2220
      INPT = IUNIT                                                      NUD 2230
      IOL = 0                                                           NUD 2240
      DO 320 M=1,12                                                     NUD 2250
         NSIZE(M) = 0                                                   NUD 2260
  320 CONTINUE                                                          NUD 2270
      DO 330 I=1,ITOT                                                   NUD 2280
         MM = MMAX(I)                                                   NUD 2290
         NSIZE(MM) = NSIZE(MM) + 1                                      NUD 2300
         IF (MM.GT.LC) WRITE (15,99999) I, NUCL(I), MM,                 NUD 2310
     *        (NPROD(MQ,I),MQ=1,MM)                                     NUD 2320
  330 CONTINUE                                                          NUD 2330
      NSIZE(13) = NON                                                   NUD 2340
      DO 340 JQ=2,5                                                     NUD 2350
         NSIZE(13+JQ) = IDM(JQ)                                         NUD 2360
  340 CONTINUE                                                          NUD 2370
      IF (MLIB(1).LE.0) RETURN                                          NUD 2380
      CALL TOC(NPAGE, 03, 05, 01, 00, TITLE, 00, NTOC, NREC)            NUD 2390
      N = 0                                                             NUD 2400
      DO 370 I=1,ITOT                                                   NUD 2410
         NUM = NON0(I)                                                  NUD 2420
         NAZ = NUCAB(I)                                                 NUD 2430
         IF (NUM.GT.0 .OR. ABUND(NAZ).GT.ERR) GO TO 350                 NUD 2440
         WRITE (LOUT,99997) I, NUCL(I), DIS(I), TOCAP(I)                NUD 2450
         GO TO 370                                                      NUD 2460
  350    N1 = N + NUM                                                   NUD 2470
         N = N + 1                                                      NUD 2480
         MUM = (NUM/5)*5 + 5                                            NUD 2490
         IOL = IOL + MUM                                                NUD 2500
         IF (IOL.LT.225 .AND. N.NE.1) GO TO 360                         NUD 2510
         IOL = 0                                                        NUD 2520
         WRITE (IOUT,99996) IOUT, NPAGE                                 NUD 2530
      CALL IDENTIFY (IOUT, 0)                                           PC
         NPAGE = IPAGE(1)                                               NUD 2540
         WRITE (LOUT,99994)                                             NUD 2550
  360    WRITE (LOUT,99995) I, NUCL(I), DIS(I), TOCAP(I),               NUD 2560
     *        (A(K),LOC(K),K=N,N1)                                      NUD 2570
         N = N1                                                         NUD 2580
  370 CONTINUE                                                          NUD 2590
      DO 400 I=1,ITOT                                                   NUD 2600
         MM = MMAX(I)                                                   NUD 2610
         IF (MM.LT.1) GO TO 400                                         NUD 2620
         DO 390 M=1,MM                                                  NUD 2630
            NP = NPROD(M,I)                                             NUD 2640
            IM=0                                                        NUD 2640
            DO 385 J=1,3                                                NUD 2641
            IL=IM+1                                                     NUD 2642
            IF(J.EQ.1)IM=IM+ILITE                                       NUD 2643
            IF(J.EQ.2)IM=IM+IACT                                        NUD 2644
            IF(J.EQ.3)IM=IM+IFP                                         NUD 2645
            IF(IM.LT.IL)GO TO 385                                       NUD 2646
            DO 380 II=IL,IM                                             NUD 2650
               IF (NP.EQ.NUCL(II)) GO TO 390                            NUD 2660
  380       CONTINUE                                                    NUD 2670
  385       CONTINUE                                                    NUD 2671
            WRITE (LOUT,99998) I, NUCL(I), M, NP, MM                    NUD 2680
  390    CONTINUE                                                       NUD 2690
  400 CONTINUE                                                          NUD 2700
      RETURN                                                            NUD 2720
  410 WRITE (15,99993) IZMAX, NON, NUCL(I), I                           NUD 2730
      STOP 'ORIGEN2 V2.1 - Error, execution terminated'                 NUD 2740
99999 FORMAT (7H0NUC #=, I3, 2X, 4HPAR=, I6, 2X, 6H# RXN=, I3, 2X,      NUD 2750
     *     6HDAUGH=, 10(I6, 1X))                                        NUD 2760
99998 FORMAT (42H NPROD(M,I) NOT FOUND IN NUCL(I) ARRAY. I=, I6, 5X,    NUD 2770
     *     7HPARENT=, I6, 5X, 2HM=, I6, 5X, 11HNPROD(M,I)=, I7, 5X,     NUD 2780
     *     8HMMAX(I)=, I6)                                              NUD 2790
99997 FORMAT (1H , I4, T7, I6, T15, 2(1PE10.3, 1X), 10(1H*), 8HTHIS NUC,NUD 2800
     *     21HLIDE IS NOT PRODUCED , 10(1H*))                           NUD 2810
99996 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        NUD 2820
99995 FORMAT (1H , I4, T7, I6, T15, 2(1PE10.3, 1X), T40, 5(1PE10.3, 1X, NUD 2830
     *     I4, 3X), 50(/, 1H , T40, 5(1PE10.3, 1X, I4, 3X)))            NUD 2840
99994 FORMAT (54H NON-ZERO MATRIX ELEMENTS; J=# OF PARENT,I=# OF DAUGHT,NUD 2850
     *     11HER, A= RATE, /, 1H , T4, 1HI, T7, 7HNUCL(I), T17,         NUD 2860
     *     6HDIS(I), T27, 8HTOCAP(I), T40, 5(2X, 6HA(I,J), 4X, 1HJ, 5X))NUD 2870
99993 FORMAT (17H0NON HAS EXCEEDED, I5, 12H  AND EQUALS, 4I6)           NUD 2880
      END                                                               NUD 2890
      SUBROUTINE NUDOC(LOUT)                                            NUD   10
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    NUD   20
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  NUD   30
     *     ITITLE(18)                                                   NUD   40
      WRITE (LOUT,99999) LOUT, NPAGE                                    NUD   50
      CALL IDENTIFY (LOUT, 1)                                           PC
      CALL TOC(NPAGE, 6, 0, 0, 0, ITITLE, IOUT, NTOC, NREC)             NUD   60
      NPAGE = IPAGE(1)                                                  NUD   70
      WRITE (LOUT,99998)                                                NUD   80
      WRITE (LOUT,99997)                                                NUD   90
      WRITE (LOUT,99996)                                                NUD  100
      WRITE (LOUT,99995)                                                NUD  110
      WRITE (LOUT,99994)                                                NUD  120
      WRITE (LOUT,99993)                                                NUD  130
      WRITE (LOUT,99992)                                                NUD  140
      WRITE (LOUT,99991)                                                NUD  150
      WRITE (LOUT,99990)                                                NUD  160
      WRITE (LOUT,99999) LOUT, NPAGE                                    NUD  170
      CALL IDENTIFY (LOUT, 0)                                           PC
      NPAGE = IPAGE(1)                                                  NUD  180
      WRITE (LOUT,99998)                                                NUD  190
      WRITE (LOUT,99989)                                                NUD  200
      WRITE (LOUT,99988)                                                NUD  210
      WRITE (LOUT,99987)                                                NUD  220
      WRITE (LOUT,99986)                                                NUD  230
      WRITE (LOUT,99985)                                                NUD  240
      RETURN                                                            NUD  250
99999 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        NUD  260
99998 FORMAT (1H0, 130(1H*), /, 1H , 130(1H*), /, 1H , 30(1H*),         NUD  270
     *     57H ORIGEN2:  A REVISED AND UPDATED VERSION OF THE ORIGEN CO,NUD  280
     *     12HMPUTER CODE , 31(1H*), /, 2(1H , 130(1H*), /))            NUD  290
99997 FORMAT (1H )                                                      NUD  300
99996 FORMAT (13H INTRODUCTION, /, 1H , T5, 24HTHIS TEXT IS INTENDED TO,NUD  310
     *     57H BE A BRIEF OUTLINE OF THE ORIGEN2 COMPUTER CODE, WHICH I,NUD  320
     *     44HS A REVISED AND UPDATED VERSION OF THE      , /, 1H , T5, NUD  330
     *     57HORIGEN DOCUMENTED IN REPORT ORNL-4628 (MAY 1973).  INCLUD,NUD  340
     *     57HED HERE ARE A BRIEF DESCRIPTION OF THE FUNCTIONS OF ORIGE,NUD  350
     *     11HN2, A      , /, 1H , T5, 29HLISTING OF THE MAJOR DATA SOU,NUD  360
     *     57HRCES, A LISTING OF THE PUBLISHED DOCUMENTATION CONCERNING,NUD  370
     *     39H ORIGEN2, AND AN OUTLINE OF THE ORIGEN2, /, 1H , T5,      NUD  380
     *     57HOUTPUT ORGANIZATION.  ORIGEN2 IS AVAILABLE FROM THE ORNL ,NUD  390
     *     57HRADIATION SHIELDING INFORMATION CENTER (RSIC) AT THE FOLL,NUD  400
     *     11HOWING      , /, 1H , T5, 29HADDRESS:                     ,NUD  410
     *     31H                               , /, 1H , T30, 9HCODES COO,NUD  420
     *     31HRDINATOR                       , /, 1H , T30, 9HRADIATION,NUD  430
     *     31H SHIELDING INFORMATION CENTER  , /, 1H , T30, 9HBLDG. 602,NUD  440
     *     31H5                              , /, 1H , T30, 9HOAK RIDGE,NUD  450
     *     31H NATIONAL LABORATORY           )                          NUD  460
99995 FORMAT (1H , T30, 40HOAK RIDGE, TENNESSEE  37830             , /, NUD  470
     *     1H , T30, 40HPHONE: (615) 574-6176                   , /,    NUD  480
     *     1H , T5, 48HQUESTIONS CONCERNING ORIGEN2 SHOULD BE ADDRESSED,NUD  490
     *     12H TO RSIC.   )                                             NUD  500
99994 FORMAT (20H-DESCRIPTION        , /, 1H , T5, 17HORIGEN2 IS A REVI,NUD  510
     *     57HSION AND UPDATE OF THE ORIGEN COMPUTER CODE. SPECIFICALLY,NUD  520
     *     51H, THE INPUT, OUTPUT, CONTROL, AND DATA BASE ASPECTS, /,   NUD  530
     *     1H , T5, 48H OF ORIGEN HAVE BEEN SIGNIFICANTLY REVISED AND U,NUD  540
     *     57HPDATED TO REFLECT CURRENT INFORMATION AND NEEDS.  IT SHOU,NUD  550
     *     20HLD BE NOTED THAT THE, /, 1H , T5, 20HMATHEMATICAL METHODS,NUD  560
     *     57H USED TO SOLVE THE NUCLIDE BUILDUP, DEPLETION, AND DECAY ,NUD  570
     *     48HEQUATIONS ARE ESSENTIALLY UNCHANGED FROM THAT IN, /, 1H , NUD  580
     *     T5, 53HORIGEN.  ORIGEN2 IS A COMPUTER CODE DESIGNED TO CALCU,NUD  590
     *     57HLATE THE COMPOSITION AND CHARACTERISTICS OF NUCLEAR MATER,NUD  600
     *     15HIALS AS A      , /, 1H , T5, 25HFUNCTION OF DECAY TIME AN,NUD  610
     *     57HD THE CHANGES THE MATERIALS UNDERGO DURING VARIOUS FUEL C,NUD  620
     *     43HYCLE OPERATIONS.  INPUT AND OUTPUT FEATURES, /, 1H , T5,  NUD  630
     *     57HHAVE BEEN DESIGNED TO FACILITATE FLEXIBILITY IN THE TYPE ,NUD  640
     *     57HOF CASES THAT CAN BE CONSIDERED AND IN CONTROLLING THE DE,NUD  650
     *     11HTAIL OF THE)                                              NUD  660
99993 FORMAT (1H , T5, 45HOUTPUT.  FOR FURTHER INFORMATION, THE USER IS,NUD  670
     *     45H REFERRED TO THE DOCUMENTATION LISTED BELOW. )            NUD  680
99992 FORMAT (30H-MAJOR DATA SOURCES           , /, 1H , T5, 8HVIRTUALL,NUD  690
     *     57HY ALL ASPECTS OF THE DATA INPUT TO ORIGEN2 HAVE BEEN UPDA,NUD  700
     *     57HTED OR REVISED TO REFLECT CURRENT INFORMATION.  THE      ,NUD  710
     *     3H   , /, 1H , T5, 38HPRINCIPAL SOURCES OF CROSS SECTION DAT,NUD  720
     *     57HA WERE THE ENDF/B-IV, ENDF/B-V, AND LENDL COMPILATIONS.  ,NUD  730
     *     30HDECAY AND PHOTON INFORMATION  , /, 1H , T5, 10HWERE PRIMA,NUD  740
     *     57HRILY BASED ON THE EVALUATED NUCLEAR STRUCTURE DATA FILE (,NUD  750
     *     57HENSDF) AT ORNL AND ENDF/B-IV.  DATA CONCERNING REACTOR   ,NUD  760
     *     1H , /, 1H , T5, 40HAND FUEL CHARACTERISTICS WERE OBTAINED F,NUD  770
     *     57HROM REFERENCE SAFETY ANALYSIS REPORTS AND, WHERE POSSIBLE,NUD  780
     *     28H, THE COMMERCIAL REACTOR    , /, 1H , T5, 8HVENDORS.)     NUD  790
99991 FORMAT (30H-DOCUMENTATION                , /, 1H , T5, 8HTHE FOLL,NUD  800
     *     57HOWING ITEMS CONSTITUTE THE ORIGEN2 DOCUMENTATION PUBLISHE,NUD  810
     *     57HD AS OF THE DATE OF THIS CODE PACKAGE:                   ,NUD  820
     *     3H   , /, 1H , T10, 37HA.G. CROFF, "ORIGEN2 - A REVISED AND ,NUD  830
     *     57HUPDATED VERSION OF THE OAK RIDGE ISOTOPE GENERATION AND D,NUD  840
     *     26HEPLETION CODE", ORNL-5621 , /, 1H , T10, 13H(JULY 1980). ,NUD  850
     *     7H       , /, 1H , T10, 33HA.G. CROFF, "A USER'S MANUAL FOR ,NUD  860
     *     57HTHE ORIGEN2 COMPUTER CODE",  ORNL/TM-7175 (JULY 1980).   ,NUD  870
     *     30H                              , /, 1H , T10, 9HA.G. CROF, NUD  880
     *     57HF, M.A. BJERKE, G.W.MORRISON, AND L.M. PETRIE, "REVISED U,NUD  890
     *     54HRANIUM-PLUTONIUM CYCLE PWR AND BWR MODELS FOR THE     ,   NUD  900
     *     /, 1H , T10, 44HORIGEN COMPUTER CODE", ORNL/TM-6051 (SEPTEMB,NUD  910
     *     16HER 1978).       , /, 1H , T10, 23HA.G. CROFF AND M.A. BJE,NUD  920
     *     57HRKE, "ALTERNATIVE FUEL CYCLE PWR MODELS FOR THE ORIGEN CO,NUD  930
     *     40HMPUTER CODE", ORNL/TM-7005 (FEB 1980).  )                 NUD  940
99990 FORMAT (1H , T10, 44HA.G. CROFF, R.L. HAESE, AND N.B. GOVE, "UPDA,NUD  950
     *     57HTED DECAY AND PHOTON LIBRARIES FOR THE ORIGEN CODE", ORNL,NUD  960
     *     19H/TM-6055 (FEB 1979), /, 1H , T10, 20HA.G. CROFF, "ORIGEN2,NUD  970
     *     57H: A REVISED AND UPDATED VERSION OF ORIGEN", TRANS. AM. NU,NUD  980
     *     43HCL. SOC., VOL. 34, P. 349-50 (JUNE 1980).  )              NUD  990
99989 FORMAT (35H ORGANIZATION OF ORIGEN2 OUTPUT    , /, 1H , T5,       NUD 1000
     *     57HPAST EXPERIENCE HAS INDICATED THAT MANY USERS ENCOUNTER C,NUD 1010
     *     57HONSIDERABLE DIFFICULTY IN FINDING THE DESIRED INFORMATION,NUD 1020
     *     11H IN A      , /, 1H , T5, 29HORIGEN2 OUT PUT WHICH IS SOME,NUD 1030
     *     57HTIMES RATHER MASSIVE.  THIS SECTION IS INTENDED AS A BRIE,NUD 1040
     *     39HF OUTLINE OF THE ORGANIZATION OF       , /, 1H , T5,      NUD 1050
     *     57HORIGEN2 OUTPUT.  FOR DETAILS REFER TO THE USER'S MANUAL (,NUD 1060
     *     57HORNL/TM-7175, SECT. 8.2).  THE ORIGEN2 OUTPUT IS EXTREMEL,NUD 1070
     *     11HY          , /, 1H , T5, 29HHIERARCHICAL, AND IS ORGANIZE,NUD 1080
     *     16HD AS FOLLOWS:   , /, 1H0, T10, 16HCARD INPUT ECHO , /,    NUD 1090
     *     1H , T10, 47HMISCELLANEOUS INPUT DATA (NEUTRON YIELDS, REPRO,NUD 1100
     *     57HCESSING LOSSES,ELEMENT CHEMICAL TOXICITIES               ,NUD 1110
     *     16H                , /, 1H , T10, 23HLISTING OF ORIGEN2 COMM,NUD 1120
     *     29HANDS CURRENTLY BEING EXECUTED, /, 1H , T10, 10HLISTING OF,NUD 1130
     *     40H ORIGEN2 DATA LIBRARIES (IF SPECIFIED)  )                 NUD 1140
99988 FORMAT (1H , T12, 13HDECAY LIBRARY, /, 1H , T14, 13HACTIVATION PR,NUD 1150
     *     6HODUCTS, /, 1H , T14, 9HACTINIDES, /, 1H , T14, 9HFISSION P,NUD 1160
     *     7HRODUCTS, /, 1H , T12, 33HCROSS SECTION/FISSION PRODUCT YIE,NUD 1170
     *     11HLD LIBRARY , /, 1H , T14, 28HACTIVATION PRODUCTS, ACTINID,NUD 1180
     *     25HES, AND FISSION PRODUCTS , /, 1H , T12, 14HPHOTON LIBRARY,NUD 1190
     *     /, 1H , T14, 44HACTIVATION PRODUCTS, ACTINIDES, AND FISSION ,NUD 1200
     *     9HPRODUCTS , /, 1H , T10, 9HOUTPUT 1 , /, 1H , T12, 6HREACTI,NUD 1210
     *     20HVITY AND BURNUP DATA, /, 1H , T12, 19HACTIVATION PRODUCT ,NUD 1220
     *     6HTABLES, /, 1H , T14, 34HGRAM TABLES (NUCLIDE, ELEMENT, NUC,NUD 1230
     *     31HLIDE SUMMARY, ELEMENT SUMMARY) , /, 1H , T14, 9HCURIE TAB,NUD 1240
     *     57HLES (NUCLIDE, ELEMENT, NUCLIDE SUMMARY, ELEMENT SUMMARY) ,NUD 1250
     *     1H , /, 1H , T14, 39HETC. (DEPENDING ON THE OUTPUT OPTIONS S,NUD 1260
     *     46HPECIFIED, MANY OF THESE TABLES MAY BE OMITTED , /, 1H ,   NUD 1270
     *     T15, 1H., /, 1H , T15, 1H.)                                  NUD 1280
99987 FORMAT (1H , T12, 15HACTINIDE TABLES, /, 1H , T14, 11HSAME SUBHEA,NUD 1290
     *     57HDINGS POSSIBLE AS UNDER ACTIVATION PRODUCT TABLES        ,NUD 1300
     *     31H                               , /, 1H , T12, 9HFISSION P,NUD 1310
     *     14HRODUCT TABLES , /, 1H , T14, 25HSAME SUBHEADINGS POSSIBLE,NUD 1320
     *     57H AS UNDER ACTIVATION PRODUCT TABLES                      ,NUD 1330
     *     17H                 , /, 1H , T12, 22HNEUTRON PRODUCTION RAT,NUD 1340
     *     57HE TABLES: (ALPHA,N) AND SPONTANEOUS FISSION              ,NUD 1350
     *     20H                    , /, 1H , T12, 13HPHOTON TABLES, /,   NUD 1360
     *     1H , T14, 47HACTIVATION PRODUCTS (SUMMATION AND PRINCIPAL CO,NUD 1370
     *     14HNTRIBUTORS)   , /, 1H , T14, 25HACTINIDES (SUMMATION AND ,NUD 1380
     *     36HPRINCIPAL CONTRIBUTORS)             , /, 1H , T14,        NUD 1390
     *     57HFISSION PRODUCTS (SUMMATION AND PRINCIPAL CONTRIBUTORS)  ,NUD 1400
     *     4H    , /, 1H , T10, 12HOUTPUT 2    , /, 1H , T12, 7HSAME GE,NUD 1410
     *     43HNERAL CONTENT AND ORDER AS OUTPUT 1        , /, 1H , T11, NUD 1420
     *     1H., /, 1H , T11, 1H., /, 1H , T11, 1H., /, 1H , T10,        NUD 1430
     *     10HOUTPUT N  , /, 1H , T12, 29HSAME GENERAL CONTENT AND ORDE,NUD 1440
     *     21HR AS OUTPUT 1        )                                    NUD 1450
99986 FORMAT (1H , T10, 44HTABLE OF CONTENTS (UNIT 12) FOR THE ABOVE (U,NUD 1460
     *     16HNIT 6) OUTPUT   , /, 1H , T10, 23HVARIABLE CROSS SECTION ,NUD 1470
     *     37HINFORMATION OUTPUT (UNIT 16)         , /, 1H , T10,       NUD 1480
     *     57HDEBUGGING AND OTHER INTERNAL INFORMATION OUTPUT (UNIT 15),NUD 1490
     *     3H   , /, 1H0, T5, 38HTHE SCENARIO LISTED ABOVE CONSTITUTES ,NUD 1500
     *     57HA TYPICAL ORIGEN2 OUTPUT FOR MANY CASES.  ONE POSSIBLE MO,NUD 1510
     *     30HDIFICATION IS THE USE OF AN   , /, 1H , T5, 10HSTP COMMAN,NUD 1520
     *     57HD TO EXECUTE AN ADDITIONAL SET OF INSTRUCTIONS AFTER THE ,NUD 1530
     *     57HFIRST SET HAS BEEN EXECUTED.  IF THIS IS DONE, THE OUTPUT,NUD 1540
     *     1H , /, 1H , T5, 40HWILL BEGIN WITH "MISCELLANEOUS INPUT DAT,NUD 1550
     *     57HA" IF NSTP=1, "ORIGEN2COMMANDS CURRENTLY BEING EXECUTED" ,NUD 1560
     *     28HIF NSTP=2, OR "OUTPUT 1"    , /, 1H , T5, 12HFOR NSTP=3. ,NUD 1570
     *     57H ANOTHER OFTEN-USED OPTION IS TO EMPLOY BOTH THE PRIMARY ,NUD 1580
     *     56H(UNIT 6) AND ALTERNATE (UNIT 11) OUTPUT UNITS.  IF BOTH ) NUD 1590
99985 FORMAT (1H , T5, 45HARE ROUTED TO PAPER, THE TABLE OF CONTENTS FO,NUD 1600
     *     57HR UNIT 11, WHICH IS ON UNIT 13, WILL IMMEDIATELY FOLLOW T,NUD 1610
     *     23HHE DEBUGGING AND       , /, 1H , T5, 17HINTERNAL INFORMAT,NUD 1620
     *     57HION (UNIT 15) OUTPUT.  THE "OUTPUT N" TABLES FOR UNIT 11 ,NUD 1630
     *     51HWILL BE PRINTED FOLLOWING THE TABLE OF CONTENTS.   )      NUD 1640
      END                                                               NUD 1650
      SUBROUTINE DECRED(I, NUCL, Q, FG, AMPC, WMPC)                     DEC   10
      DOUBLE PRECISION DNUM, UNIT, TENS                                 DEC   20
      DIMENSION NUCL(1), Q(1), FG(1), AMPC(1), WMPC(1)                  DEC   30
      COMMON /DECR/ DUM(1), DLAM, FB1, FP, FP1, FA, FT, FSF, LIBE, IU,  DEC   40
     *     ABUN, FN                                                     DEC   50
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      DEC   60
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           DEC   70
C      LEVEL 2,NUCL,Q,FG,AMPC,WMPC                                      DEC   80
      LIBE = IREAD(1)                                                   DEC   90
      IF (LIBE.LT.0) RETURN                                             DEC  100
      I = I + 1                                                         DEC  110
      NUCL(I) = IREAD(0)                                                DEC  120
      IU = IREAD(0)                                                     DEC  130
      DO 10 J=2,7                                                       DEC  140
         DUM(J) = READ(0)                                               DEC  150
   10 CONTINUE                                                          DEC  160
      LIBEE = IREAD(0)                                                  DEC  170
      FSF = READ(0)                                                     DEC  180
      FN = READ(0)                                                      DEC  190
      Q(I) = READ(0)                                                    DEC  200
      ABUN = READ(0)                                                    DEC  210
      AMPC(I) = READ(0)                                                 DEC  220
      WMPC(I) = READ(0)                                                 DEC  230
      RETURN                                                            DEC  240
      END                                                               DEC  250
      SUBROUTINE SIGRED(JJ, I, I0, I1, LITYP, NUCL)                     SIG   10
      DOUBLE PRECISION DNUM, UNIT, TENS                                 SIG   20
      DIMENSION NUCL(1)                                                 SIG   30
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      SIG   40
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           SIG   50
      COMMON /SIGR/ DDUM(1), SIGNG, SIGN2N, SIGN3N, SIGF, SNG1, SN2N1,  SIG   60
     *     Y(08), JIBE, NYLD, NCL, NCSK                                 SIG   70
C      LEVEL 2,NUCL                                                     SIG   80
      JIBE = IREAD(1)                                                   SIG   90
      IF (JIBE.LT.1) GO TO 60                                           SIG  100
      I = JJ                                                            SIG  110
      NCL = IREAD(0)                                                    SIG  120
      IF (NCL.EQ.NUCL(JJ)) GO TO 30                                     SIG  130
      DO 10 II=I0,I1                                                    SIG  140
         IF (NCL.EQ.NUCL(II)) GO TO 20                                  SIG  150
   10 CONTINUE                                                          SIG  160
      I = -2                                                            SIG  170
      RETURN                                                            SIG  180
   20 I = II                                                            SIG  190
   30 DO 40 J=2,7                                                       SIG  200
         DDUM(J) = READ(0)                                              SIG  210
   40 CONTINUE                                                          SIG  220
      IF (LITYP.NE.7) RETURN                                            SIG  230
      YESNO = READ(0)                                                   SIG  240
      IF (YESNO.LT.0.0) RETURN                                          SIG  250
      JIBE2 = IREAD(0)                                                  SIG  260
      DO 50 J=1,NYLD                                                    SIG  270
         Y(J) = READ(0)                                                 SIG  280
   50 CONTINUE                                                          SIG  290
      RETURN                                                            SIG  300
   60 I = -1                                                            SIG  310
      RETURN                                                            SIG  320
      END                                                               SIG  330
      SUBROUTINE ANSF(NUCL, Q, ALPHAN,        SFNU, ALPHN, NUCAN,       ANS   10
     *     NUCSFU, FFSF, FFA, DIS, ANMUL, ANEXP)                        ANS   20
      DIMENSION NUCL(1), Q(1), ALPHAN(1),           SFNU(1), ALPHN(1),  ANS   30
     *     NUCAN(1), NUCSFU(1), FFSF(1), FFA(1), DIS(1)                 ANS   40
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  ANS   50
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           ANS   60
     *     TFPBAV(4), IPHMAX                                            ANS   70
C      LEVEL 2,NUCL,Q,ALPHAN,      SFNU,ALPHN,NUCAN,NUCSFU,FFSF,DIS     ANS   80
      IF (IACT.EQ.0) RETURN                                             ANS   90
      BNMUL = ANMUL*ANEXP  *6.022169E+23                                ANS  100
C   ANMUL=MULTIPLIER FOR IMPURITIES (ALPHA,N ONLY)                              
C   ANEXP=MULTIPLIES FOR OTHER FACTORS (FOR ALPHA,N AND SF)                     
C     SEE INP COMMAND IN MAIN3 OR BLOCK DATA FOR DETERMINATION OF THESE         
      DO 20 I=1,IACT                                                    ANS  130
         ALPHAN(I) = 0.0                                                ANS  140
         J = I + ILITE                                                  ANS  150
         DO 10 II=1,IACT                                                ANS  160
            IF (NUCL(J).NE.NUCAN(II) .OR. NUCAN(II).LT.1) GO TO 10      ANS  170
            ALPHAN(I) = ALPHN(II)                                       ANS  180
            GO TO 20                                                    ANS  190
   10    CONTINUE                                                       ANS  200
         QJ = 0.0                                                       ANS  210
         QK = 0.0                                                       ANS  220
         IF (DIS(J).GT.ERR) QK = Q(J)/(DIS(J)*9.652E+10)                ANS  230
         IF (FFA(I).GT.ERR) QJ = (QK-FFSF(I)*200.0)/FFA(I)              ANS  240
         IF (QJ.LE.0.0 .OR. FFA(I).LT.0.05) QJ = 0.0                    ANS  250
      PNPERA=0.0                                                        ANS  251
      IF(QJ.GT.0.0)                                                     ANS  252
     $PNPERA=EXP(-30.654+11.701*ALOG(QJ)-2.382*(ALOG(QJ))**2)           ANS  255
         ALPHAN(I) = FFA(I)*DIS(J)*BNMUL*PNPERA                         ANS  260
   20 CONTINUE                                                          ANS  270
      DO 50 I=1,IACT                                                    ANS  280
         J = I + ILITE                                                  ANS  290
         DO 30 II=1,IACT                                                ANS  300
            IF (NUCL(J).NE.NUCSFU(II) .OR. NUCSFU(II).LT.1) GO TO 30    ANS  310
            ALPHN(I) = FFSF(I)*DIS(J)*SFNU(II)*6.022169E+23             ANS  320
            GO TO 40                                                    ANS  330
   30    CONTINUE                                                       ANS  340
         WT = NUCL(J)/10 - (NUCL(J)/10000)*1000                         ANS  350
         ALPHN(I) = FFSF(I)*DIS(J)*6.022169E+23*(2.84+0.1225*(WT-244))  ANS  360
      ALPHN(I)=ALPHN(I)*ANEXP                                           ANS  361
   40 CONTINUE                                                          ANS  370
   50 CONTINUE                                                          ANS  380
      DO 60 I=1,IACT                                                    ANS  390
         SFNU(I) = ALPHN(I)                                             ANS  400
   60 CONTINUE                                                          ANS  410
      RETURN                                                            ANS  420
      END                                                               ANS  430
      SUBROUTINE PHOLIB(NUCL, NGF, NGN, NGR, GGR, GFRAC, TITLE)         PHO   10
      INTEGER*2 NGF,NGN,NGR                                             PHO   20
      DOUBLE PRECISION DNUM, UNIT, TENS                                 PHO   30
      INTEGER ELE, STA                                                  PHO   40
      INTEGER AREAD, BASIS, TITLE                                       PHO   50
      DIMENSION NUCL(1), NGF(1), NGN(1), NGR(1), GGR(1)                 PHO   60
      DIMENSION EACTGP(18)                                              PHO   70
      DIMENSION GFRAC(19), TITLE(18)                                    PHO   80
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  PHO   90
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           PHO  100
     *     TFPBAV(4), IPHMAX                                            PHO  110
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 PHO  120
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            PHO  130
     *     INCYLD(8,4), MYLD(4)                                         PHO  140
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    PHO  150
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  PHO  160
     *     ITITLE(18)                                                   PHO  170
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      PHO  180
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           PHO  190
      COMMON /SZ/ NSIZE(25)                                             PHO  200
C      LEVEL 2,NUCL,NGF,NGN,NGR,GGR,GFRAC,TITLE                         PHO  210
      DATA EACTGP /0.010,0.025,0.0375,0.0575,0.085,0.125,0.225,0.375,   PHO  220
     *     0.575,0.850,1.25,1.75,2.25,2.75,3.5,5.0,7.0, 9.5/            PHO  230
      IDM(1) = 0                                                        PHO  250
      IGPMX = 18                                                        PHO  260
      INPT = MPHO(4)                                                    PHO  270
      DO 10 I=1,IPHMAX                                                  PHO  280
         GGR(I) = 0.0                                                   PHO  290
   10 CONTINUE                                                          PHO  300
      DO 20 I=1,ITOT                                                    PHO  310
         NGN(I) = 0                                                     PHO  320
         NGF(I) = 0                                                     PHO  330
   20 CONTINUE                                                          PHO  340
      NKN = 0                                                           PHO  350
      DO 160 IOV=1,3                                                    PHO  360
         IP = 0                                                         PHO  400
         IF (MPHO(IOV).LT.1 .OR. IDM(IOV+1).LT.1) GO TO 160             PHO  410
         J1 = 0                                                         PHO  420
         DO 30 J=1,IOV                                                  PHO  430
            J1 = J1 + IDM(J)                                            PHO  440
   30    CONTINUE                                                       PHO  450
         J2 = J1 + IDM(IOV+1)                                           PHO  460
         J1 = J1 + 1                                                    PHO  470
         LIHUNT = MPHO(IOV)                                             PHO  480
         DO 40 I=1,5000                                                 PHO  490
            LIBE = IREAD(1)                                             PHO  500
            IF (LIBE.EQ.LIHUNT) GO TO 50                                PHO  510
   40    CONTINUE                                                       PHO  520
         GO TO 170                                                      PHO  530
   50    CONTINUE                                                       PHO  540
         ICHR = 8                                                       PHO  550
         DO 60 JN=1,IGPMX                                               PHO  560
            TITLE(JN) = AREAD(0)                                        PHO  570
   60    CONTINUE                                                       PHO  580
         IF (MLIP(3).GT.0) CALL TOC(NPAGE, 03, 03, 01, 01, TITLE, IOUT, PHO  590
     *        NTOC, NREC)                                               PHO  600
         DO 150 J=J1,J2                                                 PHO  610
            LIBE = IREAD(1)                                             PHO  620
            IF (LIBE.LT.0) GO TO 160                                    PHO  630
            NUCLID = IREAD(0)                                           PHO  640
C    THE PROCEDURE ON THE NEXT FEW CARDS SATISFIES CDC LCM REQUIREMENTS PHO  650
            NSAV = NUCL(1)                                              PHO  660
            NUCL(1) = NUCLID                                            PHO  670
            CALL NOAH(NUCL(1))                                          PHO  680
            NUCL(1) = NSAV                                              PHO  690
            DO 70 NPH=1,IGPMX                                           PHO  700
               GFRAC(NPH) = 0.0                                         PHO  710
   70       CONTINUE                                                    PHO  720
            DO 80 NPH=1,IGPMX                                           PHO  730
               NDUM = IREAD(0)                                          PHO  740
               IF (NDUM.LT.1) GO TO 90                                  PHO  750
               GFRAC(NDUM) = READ(0)                                    PHO  760
   80       CONTINUE                                                    PHO  770
   90       IF (MLIP(3).LT.1) GO TO 110                                 PHO  780
            IF (MOD(IP,25).NE.0) GO TO 100                              PHO  790
            WRITE (IOUT,99997) IOUT, NPAGE                              PHO  800
      CALL IDENTIFY (IOUT, 0)                                           PC
            NPAGE = IPAGE(1)                                            PHO  810
            WRITE (JOUT,99996) TITLE                                    PHO  820
            WRITE (JOUT,99998)  (EACTGP(N),N= 1,IGPMX)                  PHO  830
  100       IP = IP + 1                                                 PHO  840
            WRITE (JOUT,99999) NAME, (GFRAC(N),N=1,IGPMX)               PHO  850
C    NKN=CUMULATIVE NUMBER OF NON-ZERO TERMS                            PHO  860
C    NKN=PRESENT NUMBER OF NON-ZERO TERMS                               PHO  870
C    NGF(J)=NKN OF FIRST NON-ZERO TERM FOR NUCLIDE J                    PHO  880
C    NGN(J)=NUMBER OF NON-ZERO TERMS FOR NUCLIDE J                      PHO  890
C    NGR(NKN)=PHOTON GROUP # FOR TERM NKN                               PHO  900
C    GGR(NKN)=PHOTON INTENSITY FOR TERM NKN                             PHO  910
  110       CONTINUE                                                    PHO  920
            DO 120 LL=J1,J2                                             PHO  930
               L = LL                                                   PHO  940
               IF (NUCLID.EQ.NUCL(LL)) GO TO 130                        PHO  950
  120       CONTINUE                                                    PHO  960
  130       CONTINUE                                                    PHO  970
            NGF(L) = NKN + 1                                            PHO  980
            NGN(L) = 0                                                  PHO  990
            DO 140 N=1,IGPMX                                            PHO 1000
               IF (GFRAC(N).LT.ERR) GO TO 140                           PHO 1010
               NKN = NKN + 1                                            PHO 1020
               NGN(L) = NGN(L) + 1                                      PHO 1030
               NGR(NKN) = N                                             PHO 1040
               GGR(NKN) = GFRAC(N) + GGR(NKN)                           PHO 1050
  140       CONTINUE                                                    PHO 1060
  150    CONTINUE                                                       PHO 1070
  160 CONTINUE                                                          PHO 1080
      NSIZE(20) = NKN                                                   PHO 1090
      RETURN                                                            PHO 1100
  170 WRITE (15,99995) LIHUNT                                           PHO 1110
      RETURN                                                            PHO 1120
99999 FORMAT (1H , A2, I3, A1,1P,12E10.2, /, 7X,1P,6E10.2)              PHO 1130
99998 FORMAT (8H NUCLIDE, T50, 22HMEAN GROUP ENERGY, MEV, /, 1H , 3X,   PHO 1140
     *     12F10.4, /, 1H , 3X, 6F10.4)                                 PHO 1150
99997 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        PHO 1160
99996 FORMAT (1H , 30X, 18A4)                                           PHO 1170
99995 FORMAT (48H CORRECT LIBRARY NUMBER CANNOT BE FOUND; LIHUNT=, I5)  PHO 1180
      END                                                               PHO 1190
      SUBROUTINE NUTRON(XNEW, MX, NUCL, ALPHAN,        SFNU, SFRATE,    NUT   10
     *     ALFAN, SPONT, ALN, SPN, RALFAN, RSPONT)                      NUT   20
      INTEGER ELE, STA                                                  NUT   30
      INTEGER BASIS                                                     NUT   40
      DIMENSION XNEW(MX,1), NUCL(1), ALPHAN(1),           SFNU(1)       NUT   50
      DIMENSION SFRATE(14), ALFAN(14), SPONT(14), ALN(14), SPN(14),     NUT   60
     *     RALFAN(14), RSPONT(14)                                       NUT   70
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  NUT   80
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           NUT   90
     *     TFPBAV(4), IPHMAX                                            NUT  100
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 NUT  110
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    NUT  120
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  NUT  130
     *     ITITLE(18)                                                   NUT  140
C      LEVEL 2,XNEW,NUCL,ALPHAN,      SFNU,SFRATE,ALFAN,SPONT,ALN,      NUT  150
C           SPN,RALFAN,RSPONT                                           NUT  160
      IF (IACT.LT.1) RETURN                                             NUT  170
      DO 10 M=1,MOUT                                                    NUT  180
         RSPONT(M) = 0.0                                                NUT  190
         RALFAN(M) = 0.0                                                NUT  200
         ALFAN(M) = 0.0                                                 NUT  210
         SFRATE(M) = 0.0                                                NUT  220
         SPONT(M) = 0.0                                                 NUT  230
   10 CONTINUE                                                          NUT  240
      JO = JOUT                                                         NUT  250
      ACUT = CUTOFF(NTAB+1)                                             NUT  260
      IF (ACUT.LT.0.0) GO TO 70                                         NUT  270
      DO 30 I=1,IACT                                                    NUT  280
         II = I + ILITE                                                 NUT  290
         DO 20 M=1,MOUT                                                 NUT  300
            RALFAN(M) = RALFAN(M) + XNEW(M,II)*ALPHAN(I)                NUT  310
   20    CONTINUE                                                       NUT  320
   30 CONTINUE                                                          NUT  330
      WRITE (IOUT,99999) IOUT, NPAGE                                    NUT  340
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 04, 02, 02, 00, ITITLE, IOUT, NTOC, NREC)         NUT  350
      NPAGE = IPAGE(1)                                                  NUT  360
      WRITE (JO,99995) ITITLE                                           NUT  370
C    ***          COMPUTE ALPHA-NEUTRON NEUTRON SOURCE              **  NUT  380
      WRITE (JO,99993) BASIS                                            NUT  390
      CALL HEAD(MOUT, 1, IOUT)                                          NUT  400
      DO 60 I=1,IACT                                                    NUT  410
         I1 = I + ILITE                                                 NUT  420
         NUCLI = NUCL(I1)                                               NUT  430
         CALL NOAH(NUCL(I1))                                            NUT  440
         DO 40 M=1,MOUT                                                 NUT  450
            ALN(M) = XNEW(M,I1)*ALPHAN(I)                               NUT  460
            ALFAN(M) = ALFAN(M) + ALN(M)                                NUT  470
   40    CONTINUE                                                       NUT  480
         DO 50 M=1,MOUT                                                 NUT  490
            IF (RALFAN(M).LT.ERR) GO TO 50                              NUT  500
            IF ((ALN(M)/RALFAN(M)).LT.ACUT) GO TO 50                    NUT  510
            WRITE (IOUT,99991) NAME, (ALN(MM),MM=1,MOUT)                NUT  520
            GO TO 60                                                    NUT  530
   50    CONTINUE                                                       NUT  540
   60 CONTINUE                                                          NUT  550
      WRITE (IOUT,99997)                                                NUT  560
      WRITE (IOUT,99998) (ALFAN(M),M=1,MOUT)                            NUT  570
      WRITE (IOUT,99994) (RALFAN(M),M=1,MOUT)                           NUT  580
   70 SCUT = CUTOFF(NTAB+2)                                             NUT  590
      IF (SCUT.LT.0.0) GO TO 140                                        NUT  600
      DO 90 I=1,IACT                                                    NUT  610
         II = I + ILITE                                                 NUT  620
         DO 80 M=1,MOUT                                                 NUT  630
            RSPONT(M) = RSPONT(M) + XNEW(M,II)*SFNU(I)                  NUT  640
   80    CONTINUE                                                       NUT  650
   90 CONTINUE                                                          NUT  660
      WRITE (IOUT,99999) IOUT, NPAGE                                    NUT  670
      CALL IDENTIFY (IOUT, 0)                                           PC
      CALL TOC(NPAGE, 04, 02, 01, 00, ITITLE, IOUT, NTOC, NREC)         NUT  680
      NPAGE = IPAGE(1)                                                  NUT  690
      WRITE (JO,99995) ITITLE                                           NUT  700
C    ***          COMPUTE SPONTANEOUS FISSION NEUTRON SOURCE        **  NUT  710
      WRITE (JO,99992) BASIS                                            NUT  720
      CALL HEAD(MOUT, 1, IOUT)                                          NUT  730
      DO 120 I=1,IACT                                                   NUT  740
         I1 = I + ILITE                                                 NUT  750
         NUCLI = NUCL(I1)                                               NUT  760
         CALL NOAH(NUCL(I1))                                            NUT  770
         DO 100 M=1,MOUT                                                NUT  780
            SPN(M) = XNEW(M,I1)*SFNU(I)                                 NUT  790
            SPONT(M) = SPONT(M) + SPN(M)                                NUT  800
  100    CONTINUE                                                       NUT  820
         DO 110 M=1,MOUT                                                NUT  830
            IF (RSPONT(M).LT.ERR) GO TO 110                             NUT  840
            IF ((SPN(M)/RSPONT(M)).LT.SCUT) GO TO 110                   NUT  850
            WRITE (IOUT,99991) NAME, (SPN(J),J=1,MOUT)                  NUT  860
            GO TO 120                                                   NUT  870
  110    CONTINUE                                                       NUT  880
  120 CONTINUE                                                          NUT  890
      WRITE (IOUT,99997)                                                NUT  900
      WRITE (IOUT,99998) (SPONT(M),M=1,MOUT)                            NUT  910
      WRITE (IOUT,99994) (RSPONT(M),M=1,MOUT)                           NUT  920
      DO 130 M=1,MOUT                                                   NUT  930
         RALFAN(M) = RALFAN(M) + RSPONT(M)                              NUT  940
         ALFAN(M) = ALFAN(M) + SPONT(M)                                 NUT  950
  130 CONTINUE                                                          NUT  960
      WRITE (IOUT,99996)                                                NUT  970
      WRITE (IOUT,99998) (ALFAN(M),M=1,MOUT)                            NUT  980
      WRITE (IOUT,99994) (RALFAN(M),M=1,MOUT)                           NUT  990
  140 RETURN                                                            NUT 1000
99999 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        NUT 1010
99998 FORMAT (7H TOTALS, /, 6H TABLE, 5X,1P,12E10.3)                    NUT 1020
99997 FORMAT (1H0, 130(1H-))                                            NUT 1030
99996 FORMAT (1H0, 130(1H-), /, 8H OVERALL)                             NUT 1040
99995 FORMAT (1H , 30X, 18A4)                                           NUT 1050
99994 FORMAT (7H ACTUAL, 4X,1P,12E10.3)                                 NUT 1060
99993 FORMAT (1H , T30, 38H(ALPHA,N) NEUTRON SOURCE, NEUTRONS/SEC, /,   NUT 1070
     *     1H , 6HBASIS=, 18A4)                                         NUT 1080
99992 FORMAT (1H , T30, 44HSPONTANEOUS FISSION NEUTRON SOURCE, NEUTRONS,NUT 1090
     *     4H/SEC, /, 1H , T28, 6HBASIS=, 18A4)                         NUT 1100
99991 FORMAT (1H , A2, I3, A1, 4X,1P,12E10.3)                           NUT 1110
      END                                                               NUT 1120
      SUBROUTINE GAMMA(NUCL, DIS, NGF, NGN, NGR, GGR, XNEW, MX, SFRATE, GAM   10
     *     GSUM, ENERGY, GTOTAL, PSTAR, GTOTA, MSTAR)                   GAM   20
      INTEGER*2 NGF,NGN,NGR                                             GAM   30
      INTEGER ELE, STA                                                  GAM   40
      INTEGER BASIS                                                     GAM   50
      DIMENSION IGP(3), JGP(3), EGROUP(18,2), GSUM(18,13), SFRATE(14),  GAM   60
     *     ENERGY(14), GTOTAL(14), PSTAR(14), GTOTA(14)                 GAM   70
      DIMENSION NUCL(1), DIS(1), NGF(1), NGN(1), NGR(1), GGR(1),        GAM   80
     *     XNEW(MX,1)                                                   GAM   90
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  GAM  100
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           GAM  110
     *     TFPBAV(4), IPHMAX                                            GAM  120
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 GAM  130
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            GAM  140
     *     INCYLD(8,4), MYLD(4)                                         GAM  150
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    GAM  160
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  GAM  170
     *     ITITLE(18)                                                   GAM  180
C       LEVEL 2,NUCL,DIS,NGF,NGN,NGR,GGR,XNEW,SFRATE,GSUM,ENERGY,GTOTAL,GAM  190
C           PSTAR,GTOTA                                                 GAM  200
C                                                                       GAM  210
C    ENERGY GROUP STRUCTURES                                            GAM  220
C                                                                       GAM  230
C     IGP = # OF ENERGY GROUPS FOR ACTIV PROD, ACTINIDES, AND FISS PROD GAM  240
C     JGP = WHICH SET OF ENERGY GROUP MIDPOINTS (1 OR 2) IS USED FOR    GAM  250
C           ACTIV PROD, ACTINIDES, AND FISS PROD                        GAM  260
C     EGROUP = ENERGY GROUP MIDPOINTS IN TWO SETS (MAX)                 GAM  270
C                                                                       GAM  280
      DATA IGP, JGP /18,18,18,1,1,1/                                    GAM  290
      DATA EGROUP /0.010,0.025,0.0375,0.0575,0.085,0.125,0.225,0.375,   GAM  300
     *     0.575,0.850,1.25,1.75,2.25,2.75,3.5,5.0,7.0, 9.5,0.150,0.375,GAM  310
     *     0.575,0.850,1.25,1.75,2.25,2.75,3.5,5.0,7.0,11.0,6*0.0/      GAM  320
      IDM(1) = 0                                                        GAM  330
      I = 0                                                             GAM  340
      DO 250 JTYP=1,3                                                   GAM  350
         IF (IDM(JTYP+1).LT.1 .OR. MPHO(JTYP).LT.1) GO TO 250           GAM  360
         IGPMX = IGP(JTYP)                                              GAM  370
         JGPN = JGP(JTYP)                                               GAM  380
         I0 = 0                                                         GAM  390
         DO 10 J=1,JTYP                                                 GAM  400
            I0 = I0 + IDM(J)                                            GAM  410
   10    CONTINUE                                                       GAM  420
         I1 = I0 + IDM(JTYP+1)                                          GAM  430
         I0 = I0 + 1                                                    GAM  440
         WRITE (IOUT,99989) IOUT, NPAGE                                 GAM  450
      CALL IDENTIFY (IOUT, 0)                                           PC
         CALL TOC(NPAGE, 04, 02, JTYP+2, 00, ITITLE, IOUT, NTOC, NREC)  GAM  460
         NPAGE = IPAGE(1)                                               GAM  470
         IF (JTYP.EQ.1) WRITE (IOUT,99986)                              GAM  480
         IF (JTYP.EQ.2) WRITE (IOUT,99984)                              GAM  490
         IF (JTYP.EQ.3) WRITE (IOUT,99985)                              GAM  500
         WRITE (IOUT,99987) ITITLE                                      GAM  510
         DO 30 M=1,MOUT                                                 GAM  520
            ENERGY(M) = 0.0                                             GAM  530
            GTOTAL(M) = 0.0                                             GAM  540
            GTOTA(M) = 0.0                                              GAM  550
            DO 20 N=1,IGPMX                                             GAM  560
               GSUM(N,M) = 0.0                                          GAM  570
   20       CONTINUE                                                    GAM  580
   30    CONTINUE                                                       GAM  590
C    *** PHOTONS/SEC TABLE                                              GAM  600
         DO 60 I=I0,I1                                                  GAM  610
            IF (NGN(I).LT.1) GO TO 60                                   GAM  620
            DISI = DIS(I)*6.023E+23                                     GAM  630
            DO 50 M=1,MOUT                                              GAM  640
               N1 = NGF(I)                                              GAM  650
               N2 = NGF(I) + NGN(I) - 1                                 GAM  660
               DO 40 K=N1,N2                                            GAM  670
                  N = NGR(K)                                            GAM  680
                  GSUM(N,M) = GSUM(N,M) + XNEW(M,I)*DISI*GGR(K)         GAM  690
   40          CONTINUE                                                 GAM  700
   50       CONTINUE                                                    GAM  710
   60    CONTINUE                                                       GAM  720
         WRITE (IOUT,99988) TFPBAV(3), TFPBAV(4), TFPBAV(2)             GAM  730
         WRITE (IOUT,99998) IGPMX, BASIS                                GAM  740
         WRITE (IOUT,99993)                                             GAM  750
         CALL HEAD(MOUT, 1, IOUT)                                       GAM  760
         DO 80 N=1,IGPMX                                                GAM  770
            WRITE (IOUT,99999) EGROUP(N,JGPN), (GSUM(N,M),M=1,MOUT)     GAM  780
            DO 70 M=1,MOUT                                              GAM  790
               GTOTAL(M) = GTOTAL(M) + GSUM(N,M)                        GAM  800
               ENERGY(M) = ENERGY(M) + GSUM(N,M)*EGROUP(N,JGPN)         GAM  810
   70       CONTINUE                                                    GAM  820
   80    CONTINUE                                                       GAM  830
         WRITE (IOUT,99994) (GTOTAL(M),M=1,MOUT)                        GAM  840
         WRITE (IOUT,99996) (ENERGY(M),M=1,MOUT)                        GAM  850
C    *** MEV/WATT-SEC TABLE                                             GAM  860
         IF (TFPBAV(3).LT.ERR) GO TO 140                                GAM  870
         CONST = 1.0E-06/TFPBAV(3)                                      GAM  880
         DO 90 M=1,MOUT                                                 GAM  890
            IF (GTOTAL(M).GT.ERR) GTOTA(M) = 1.0/GTOTAL(M)              GAM  900
            GTOTAL(M) = ENERGY(M)*CONST                                 GAM  910
   90    CONTINUE                                                       GAM  920
         WRITE (IOUT,99997) IGPMX, BASIS                                GAM  930
         WRITE (IOUT,99993)                                             GAM  940
         CALL HEAD(MOUT, 1, IOUT)                                       GAM  950
         DO 110 N=1,IGPMX                                               GAM  960
            DO 100 M=1,MOUT                                             GAM  970
               GSUM(N,M) = GSUM(N,M)*EGROUP(N,JGPN)*CONST               GAM  980
  100       CONTINUE                                                    GAM  990
            WRITE (IOUT,99999) EGROUP(N,JGPN), (GSUM(N,M),M=1,MOUT)     GAM 1000
  110    CONTINUE                                                       GAM 1010
         DO 130 M=1,MOUT                                                GAM 1020
            DO 120 N=1,IGPMX                                            GAM 1030
               GSUM(N,M) = GSUM(N,M)/(EGROUP(N,JGPN)*CONST)             GAM 1040
  120       CONTINUE                                                    GAM 1050
            ENERGY(M) = ENERGY(M)*1.603E-13                             GAM 1060
  130    CONTINUE                                                       GAM 1070
         WRITE (IOUT,99994) (GTOTAL(M),M=1,MOUT)                        GAM 1080
         WRITE (IOUT,99995) (ENERGY(M),M=1,MOUT)                        GAM 1090
  140    CONTINUE                                                       GAM 1100
C    *** PRINCIPAL PHOTON SOURCES                                       GAM 1110
         WRITE (IOUT,99989) IOUT, NPAGE                                 GAM 1120
      CALL IDENTIFY (IOUT, 0)                                           PC
         NPAGE = IPAGE(1)                                               GAM 1130
         PHCUT = CUTOFF(NTAB+3)                                         GAM 1140
         LINE = 50                                                      GAM 1150
         DO 240 N=1,IGPMX                                               GAM 1160
            DO 150 M=1,MOUT                                             GAM 1170
C    ELIMINATE LOW-INTENSITY GROUPS                                     GAM 1180
               IF (MSTAR.GT.0 .AND. MSTAR.NE.M .OR. MSTAR.LT.0 .AND.    GAM 1190
     *              (-MSTAR).EQ.M) GO TO 150                            GAM 1200
               IF ((GSUM(N,M)*GTOTA(M)).GT.ERR) GO TO 160               GAM 1210
  150       CONTINUE                                                    GAM 1220
            GO TO 240                                                   GAM 1230
  160       IF (LINE.GT.10) GO TO 170                                   GAM 1240
            WRITE (IOUT,99989) IOUT, NPAGE                              GAM 1250
      CALL IDENTIFY (IOUT, 0)                                           PC
            NPAGE = IPAGE(1)                                            GAM 1260
            LINE = 50                                                   GAM 1270
  170       WRITE (IOUT,99992) N, EGROUP(N,JGPN)                        GAM 1280
            WRITE (IOUT,99990)                                          GAM 1290
            CALL HEAD(MOUT, 1, IOUT)                                    GAM 1300
            LINE = LINE - 5                                             GAM 1310
            DO 230 I=I0,I1                                              GAM 1320
               IF (NGN(I).LT.1) GO TO 230                               GAM 1330
               N1 = NGF(I)                                              GAM 1340
               N2 = N1 + NGN(I) - 1                                     GAM 1350
               DO 180 K=N1,N2                                           GAM 1360
                  IF (N.NE.NGR(K)) GO TO 180                            GAM 1370
                  DISI = DIS(I)*GGR(K)*6.023E+23                        GAM 1380
                  GO TO 190                                             GAM 1390
  180          CONTINUE                                                 GAM 1400
               GO TO 230                                                GAM 1410
  190          DO 200 M=1,MOUT                                          GAM 1420
C    ELIMINATE LOW-INTENSITY CONTRIBUTORS IN HIGH-INTENSITY GROUPS      GAM 1430
                  IF (MSTAR.GT.0 .AND. MSTAR.NE.M .OR. MSTAR.LT.0 .AND. GAM 1440
     *                 (-MSTAR).EQ.M) GO TO 200                         GAM 1450
                  IF (GSUM(N,M).LT.ERR) GO TO 200                       GAM 1460
                  IF (((XNEW(M,I)*DISI)/GSUM(N,M)).GT.PHCUT) GO TO 210  GAM 1470
  200          CONTINUE                                                 GAM 1480
               GO TO 230                                                GAM 1490
  210          CALL NOAH(NUCL(I))                                       GAM 1500
               DO 220 M=1,MOUT                                          GAM 1510
                  PSTAR(M) = DISI*XNEW(M,I)                             GAM 1520
  220          CONTINUE                                                 GAM 1530
               WRITE (IOUT,99991) NAME, (PSTAR(M),M=1,MOUT)             GAM 1540
               LINE = LINE - 1                                          GAM 1550
               IF (LINE.GT.0) GO TO 230                                 GAM 1560
               WRITE (IOUT,99989) IOUT, NPAGE                           GAM 1570
      CALL IDENTIFY (IOUT, 0)                                           PC
               NPAGE = IPAGE(1)                                         GAM 1580
               LINE = 50                                                GAM 1590
  230       CONTINUE                                                    GAM 1600
  240    CONTINUE                                                       GAM 1610
  250 CONTINUE                                                          GAM 1620
      RETURN                                                            GAM 1630
99999 FORMAT (1PE10.3, 1X,1P,12E10.3)                                   GAM 1640
99998 FORMAT (1H0, T35, I2, 40H GROUP PHOTON RELEASE RATES, PHOTONS/SEC,GAM 1650
     *     3HOND, /, 1H , T35, 6HBASIS=, 18A4)                          GAM 1660
99997 FORMAT (1H0, T35, I2, 40H GROUP SPECIFIC ENERGY RELEASE RATES, ME,GAM 1670
     *     10HV/WATT-SEC, /, 1H , T35, 6HBASIS=, 18A4)                  GAM 1680
99996 FORMAT (8H0MEV/SEC, 3X,1P,12E10.3)                                GAM 1690
99995 FORMAT (9H0 GAM POW, 2X,1P,12E10.3)                               GAM 1700
99994 FORMAT (8H0  TOTAL, 3X,1P,12E10.3)                                GAM 1710
99993 FORMAT (7H  EMEAN)                                                GAM 1720
99992 FORMAT (1H0, 20X, 33HPRINCIPAL PHOTON SOURCES IN GROUP, I3,       GAM 1730
     *     12H,PHOTONS/SEC, /, 1H , 21X, 12HMEAN ENERGY=, F6.3, 3HMEV,  GAM 1740
     *     /)                                                           GAM 1750
99991 FORMAT (1H , A2, I3, A1, 4X,1P,12E10.3)                           GAM 1760
99990 FORMAT (8H NUCLIDE)                                               GAM 1770
99989 FORMAT (1H1, T80, 14HOUTPUT UNIT = , I2, T120, 4HPAGE, I5)        GAM 1780
99988 FORMAT (1H , 35X, 6HPOWER=, F7.2, 12H MW, BURNUP=, F8.0, 6H MWD, ,GAM 1790
     *     5HFLUX=, 1PE9.2, 12H N/CM**2-SEC)                            GAM 1800
99987 FORMAT (1H0, 30X, 18A4)                                           GAM 1810
99986 FORMAT (1H , 15X, 39HPHOTON SPECTRUM FOR ACTIVATION PRODUCTS)     GAM 1820
99985 FORMAT (1H , 15X, 36HPHOTON SPECTRUM FOR FISSION PRODUCTS)        GAM 1830
99984 FORMAT (1H , 15X, 41HPHOTON SPECTRUM FOR ACTINIDES + DAUGHTERS)   GAM 1840
      END                                                               GAM 1850
C
C  pc vers of HEAD - SEEMS TO BE MORE COMPATIBLE WITH OTHER HARDWARE
C
      SUBROUTINE HEAD(MMN, K, IOUT)
      INTEGER ELE,STA
      CHARACTER*132 BUFF
      CHARACTER*4 LABTIM, IOUTTI
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)
      COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)
C    THIS SUBROUTINE PRINTS VECTOR HEADINGS
      WRITE (BUFF,9000)
      IBUF = 13
      DO 10 M=1,MMN
          IOM = IO(M)
          IF(IOM .LT. 0) THEN
              WRITE (BUFF(IBUF:IBUF+9),9010) (IOUTTI(M,J), J=1,3)
          ELSE IF(ROUTTI(M) .GE. 10000.0) THEN
              WRITE (BUFF(IBUF:IBUF+9),9020) ROUTTI(M), LABTIM(IOM)
          ELSE
              WRITE (BUFF(IBUF:IBUF+9),9030) ROUTTI(M), LABTIM(IOM)
          ENDIF
          IBUF = IBUF + 10
   10 CONTINUE
      WRITE (IOUT,9040) BUFF
      RETURN
 9000     FORMAT(12X)
 9010     FORMAT(2A4,A2)
 9020     FORMAT(1PE7.1,A2,1X)
 9030     FORMAT(0PF7.1,A2,1X)
 9040     FORMAT(A132)
      END
C
C   OLD MAINFRAME VERSION OF HEAD
C     SUBROUTINE HEAD(MMN, K, IOUT)                                     HEA   10
C     INTEGER ELE,STA                                                   HEA   20
C     REAL LABTIM,IOUTTI                                                HEA   30
C     DIMENSION FOR(54),AL(4),RIMP(4),SPA(2),PAREN(3),SIMP(4),HEADO(36) HEA   40
C     COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 HEA   50
C     COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)   HEA   60
C    THIS SUBROUTINE PRINTS VECTOR HEADINGS                             HEA   70
C     DATA SIMP /4H1PE7, 4H.1, , 4HA2,1, 4HX,  /                        HEA   80
C     DATA RIMP /4H0PF7, 4H.1, , 4HA2,1, 4HX,  /                        HEA   90
C     DATA AL   /4H2A4,, 4HA2, , 4H    , 4H    /                        HEA  100
C     DATA SPA  /4H1H ,, 4H11X,/                                        HEA  110
C     DATA PAREN /1H(   , 1H)   , 1H/   /                               HEA  120
C     FOR(1)=PAREN(1)                                                   HEA  130
C     FOR(2)=SPA(1)                                                     HEA  140
C     FOR(3)=SPA(2)                                                     HEA  150
C     JJ=3                                                              HEA  160
C     JK=0                                                              HEA  170
C     DO 10 M=1,MMN                                                     HEA  180
C     DO 20 J=1,4                                                       HEA  190
C     FOR(JJ+J)=RIMP(J)                                                 HEA  200
C     IF(ROUTTI(M).GE.10000.0) FOR(JJ+J)=SIMP(J)                        HEA  210
C     IF(IO(M).LT.0)FOR(JJ+J)=AL(J)                                     HEA  220
C  20 CONTINUE                                                          HEA  230
C     JJ=JJ+4                                                           HEA  240
C     IF(IO(M).LT.0)GO TO 30                                            HEA  250
C     HEADO(JK+1)=ROUTTI(M)                                             HEA  260
C     HEADO(JK+2)=LABTIM(IO(M))                                         HEA  270
C     JK=JK+2                                                           HEA  280
C     GO TO 10                                                          HEA  290
C  30 DO 40 JL=1,3                                                      HEA  300
C     HEADO(JK+JL)=IOUTTI(M,JL)                                         HEA  310
C  40 CONTINUE                                                          HEA  320
C     JK=JK+3                                                           HEA  330
C  10 CONTINUE                                                          HEA  340
C     FOR(JJ+1)=PAREN(3)                                                HEA  350
C     FOR(JJ+2)=PAREN(2)                                                HEA  360
C     WRITE(IOUT,FOR) (HEADO(J),J=1,JK)                                 HEA  370
C     RETURN                                                            HEA  380
C     END                                                               HEA  390
      SUBROUTINE TOC(NPAG, NA, NB, MC, ND, TITLE, NT, MB, MV)           TOC   10
      INTEGER TITLE                                                     TOC   20
      DIMENSION TITLE(18)                                               TOC   30
      COMMON /BDAT/ ZAP(130), KZAP(12), JZAP(6,26), NTB(3,4),           TOC   40
     *     INSTR0(32), RRCUT(30), NDI(11,7), NDL(6,5), NDO(9,5)         TOC   50
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    TOC   60
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  TOC   70
     *     ITITLE(18)                                                   TOC   80
      IF (NPAGE.NE.1) GO TO 10                                          TOC   90
      IF (ND.EQ.08) GO TO 10                                            TOC  100
      IF (NTOC.EQ.NTOCA) IPA = 0                                        TOC  110
      IF (NTOC.EQ.NTOCB) IPB = 0                                        TOC  120
   10 CONTINUE                                                          TOC  130
      IF (NTOC.EQ.NTOCA) IP = IPA                                       TOC  140
      IF (NTOC.EQ.NTOCB) IP = IPB                                       TOC  150
      IF (IP.LT.50 .AND. IP.GT.0) GO TO 20                              TOC  160
      WRITE (MB,99995) MB, NT                                           TOC  170
      CALL IDENTIFY (MB, 0)                                             PC
      IP = 1                                                            TOC  180
   20 CONTINUE                                                          TOC  190
      GO TO (30, 50, 60, 80, 120, 130), NA                              TOC  200
   30 IF (NPAGE.NE.1) GO TO 40                                          TOC  210
      WRITE (MB,99998) NPAGE                                            TOC  220
      IP = IP + 3                                                       TOC  230
      GO TO 140                                                         TOC  240
   40 WRITE (MB,99997) NPAGE, (NDI(J,NB),J=1,11)                        TOC  250
      IP = IP + 1                                                       TOC  260
      GO TO 140                                                         TOC  270
   50 WRITE (MB,99994) NPAGE                                            TOC  280
      IP = IP + 3                                                       TOC  290
      GO TO 140                                                         TOC  300
   60 IF (MC.NE.0) GO TO 70                                             TOC  310
      WRITE (MB,99993) NPAGE                                            TOC  320
      IP = IP + 3                                                       TOC  330
      GO TO 140                                                         TOC  340
   70 WRITE (MB,99992) NPAGE, (NDL(J,NB),J=1,6), TITLE                  TOC  350
      IP = IP + 1                                                       TOC  360
      GO TO 140                                                         TOC  370
   80 IF (MC.NE.0) GO TO 90                                             TOC  380
      WRITE (MB,99991) NPAGE, TITLE, NREC                               TOC  390
      IP = IP + 3                                                       TOC  400
      GO TO 140                                                         TOC  410
   90 IF (NB.NE.0) GO TO 100                                            TOC  420
      WRITE (MB,99990) NPAGE                                            TOC  430
      IP = IP + 1                                                       TOC  440
      GO TO 140                                                         TOC  450
  100 IF (NB.NE.1) GO TO 110                                            TOC  460
      K1 = KZAP(4*(MC-1)+1)                                             TOC  470
      K2 = K1 + 3                                                       TOC  480
      J1 = JZAP(1,ND)                                                   TOC  490
      J2 = JZAP(2,ND)                                                   TOC  500
      J3 = JZAP(3,ND)                                                   TOC  510
      J4 = JZAP(4,ND)                                                   TOC  520
      J5 = JZAP(5,ND)                                                   TOC  530
      J6 = JZAP(6,ND)                                                   TOC  540
      WRITE (MB,99997) NPAGE, (ZAP(N),N=J1,J2), (ZAP(N),N=J3,J4),       TOC  550
     *     (ZAP(N),N=J5,J6), (ZAP(N),N=K1,K2)                           TOC  560
      IP = IP + 1                                                       TOC  570
      GO TO 140                                                         TOC  580
  110 WRITE (MB,99997) NPAGE, (NDO(J,MC),J=1,9)                         TOC  590
      IP = IP + 1                                                       TOC  600
      GO TO 140                                                         TOC  610
  120 WRITE (MB,99996) NPAGE, NB, MC, ND                                TOC  620
      GO TO 140                                                         TOC  630
  130 WRITE (MB,99999) NPAGE                                            TOC  640
  140 IF (NTOC.EQ.NTOCA) IPA = IP                                       TOC  650
      IF (NTOC.EQ.NTOCB) IPB = IP                                       TOC  660
      RETURN                                                            TOC  670
99999 FORMAT (1H , I5, T20, 35HORIGEN2 OVERVIEW AND DOCUMENTATION )     TOC  680
99998 FORMAT (1H0, I5, T10, 19HDATA INITIALIZATION, /)                  TOC  690
99997 FORMAT (1H , I5, T20, 25A4)                                       TOC  700
99996 FORMAT (1H , I5, T20, 20HINPUT ECHO; READ ON , I2, 9H LIST ON ,   TOC  710
     *     I2, 9H COPY TO , I2)                                         TOC  720
99995 FORMAT (1H1, 5H PAGE, T30, 28HTABLE OF CONTENTS ON UNIT = , I2,   TOC  730
     *     19H FOR OUTPUT UNIT = , I2)                                  TOC  740
99994 FORMAT (1H0, I5, T10, 33HORIGEN INSTRUCTIONS FOR THIS CASE, /)    TOC  750
99993 FORMAT (1H0, I5, T10, 22HNUCLIDE DATA LIBRARIES, /)               TOC  760
99992 FORMAT (1H , I5, T20, 6A4, 2H--, 18A4)                            TOC  770
99991 FORMAT (1H0, I5, T10, 21HOUTPUT TABLES--TITLE=, 18A4, 5X,         TOC  780
     *     11HRECYCLE # =, I3)                                          TOC  790
99990 FORMAT (1H , I5, T20, 26HREACTIVITY AND BURNUP DATA)              TOC  800
      END                                                               TOC  810
      SUBROUTINE FLUXO(DELT, OLDFLX, MM, TSEC, A, LOC, B, D, NON0, KD,  FLU   10
     *     LONG, NUCL, TOCAP, XDOT, XTEMP, XTEM, IOUT, FISS, DIS,       FLU   20
     *     NFUDFP, IFD, MIRR, FPRATN, IJKO, IFIRST)                     FLU   30
      LOGICAL LONG                                                      FLU   40
      INTEGER*2 NFUDFP                                                  FLU   50
      INTEGER*2 LOC,NON0,KD                                             FLU   60
      DIMENSION NDEF(3), LDEF(3), NSORS(8), NFUDFP(IFD,1)               FLU   70
      DIMENSION A(1), LOC(1), B(1), D(1), NON0(1), KD(1), LONG(1),      FLU   80
     *     NUCL(1), TOCAP(1), XDOT(1), XTEMP(1), XTEM(1), FISS(1),      FLU   90
     *     DIS(1)                                                       FLU  100
      COMMON /FUDG/ NFUD(7,20), FUDSIG(32,20), CBUP(32), NXSEC, NSM,    FLU  110
     *     NPT, NCO(2,20), FPBST(050,3), VECT(050), FOLD(20), FNEW(20)  FLU  120
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  FLU  130
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           FLU  140
     *     TFPBAV(4), IPHMAX                                            FLU  150
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            FLU  160
     *     INCYLD(8,4), MYLD(4)                                         FLU  170
      COMMON /PROCSS/ MPROS, PRATE(8), NOPROS(8), NZPROS(8,20)          FLU  180
      COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)   FLU  190
C      LEVEL 2,LOC,B,D,NON0,KD,NUCL,TOCAP,XDOT,XTEMP,XTEM,FISS,DIS      FLU  200
C    1  ,A,NFUDFP                                                       FLU  210
      DATA FDOT, FDDOT, FISSR /3*0.0/                                   FLU  220
      DATA EPF1, EPF2, EPF3, EPFAVG /4*1.0/                             FLU  230
      DATA NSORS /902320,922330,922350,922380,942390,942410,962450,     FLU  240
     *     982520/                                                      FLU  250
      DATA NDEF, LDEF /922330,922350,942390,3*0/                        FLU  260
C                                                                       FLU  270
C    FLUX0 COMPUTES THE AVERAGE FLUX OVER A TIME INTERVAL USING A       FLU  280
C    TAYLOR SERIES EXPANSION ABOUT THE BEGINNING OF THE INTERVAL        FLU  290
C                                                                       FLU  300
      FISSEN(AWT,Z) = 1.29927E-03*(Z**2*SQRT(AWT)) + 33.12              FLU  310
      KOUNT = 0                                                         FLU  320
      IF (INDEX.EQ.2) GO TO 380                                         FLU  330
      IF (IFP.LT.1) GO TO 150                                           FLU  340
      NNSORS = MYLD(MLIB(11))                                           FLU  350
      IF (MLIB(11).EQ.0) GO TO 150                                      FLU  360
      IDIF = 100000                                                     FLU  370
      IHI = 0                                                           FLU  380
      XSHI = 0.0                                                        FLU  390
      XSTOT = 0.0                                                       FLU  400
C    FIND LARGEST UNCONNECTED ACTINIDE FISSION SOURCE                   FLU  410
      DO 50 I=1,IACT                                                    FLU  420
         L = I + ILITE                                                  FLU  430
C    STORE LOCATION OF DEFAULT CONNECTED NUCLIDES                       FLU  440
         DO 20 K=1,NNSORS                                               FLU  450
            J = INCYLD(K,MLIB(11))                                      FLU  460
            DO 10 M=1,3                                                 FLU  470
               IF (NDEF(M).EQ.NSORS(J) .AND. NDEF(M).EQ.NUCL(L))        FLU  480
     *              LDEF(M) = L                                         FLU  490
   10       CONTINUE                                                    FLU  500
   20    CONTINUE                                                       FLU  510
         PR = XTEMP(L)*FISS(I)                                          FLU  520
         IF (PR.LE.XSHI) GO TO 40                                       FLU  530
C    SEE IF NUCLIDE IS CONNECTED; IF SO, DON'T STORE DATA               FLU  540
         DO 30 K=1,NNSORS                                               FLU  550
            J = INCYLD(K,MLIB(11))                                      FLU  560
            IF (NUCL(L).EQ.NSORS(J)) GO TO 40                           FLU  570
   30    CONTINUE                                                       FLU  580
C    STORE LOCATION OF LARGEST UNCONNECTED ACTINIDE                     FLU  590
         IHI = L                                                        FLU  600
         XSHI = PR                                                      FLU  610
C    TOTAL XTEMP*FISS                                                   FLU  620
         XSTOT = XSTOT + PR                                             FLU  630
   40    CONTINUE                                                       FLU  640
   50 CONTINUE                                                          FLU  650
      IF (IHI.EQ.0) GO TO 150                                           FLU  660
C    FIND CLOSEST CONNECTED ACT TO THE LARGEST UNCONNECTED CONTRIBUTOR  FLU  670
C        FIRST CASE: WITH VARIABLE XSEC                                 FLU  680
      IF (NXSEC.LT.1) GO TO 70                                          FLU  690
      DO 60 L=1,NXSEC                                                   FLU  700
         IF (NFUD(6,L).LT.1) GO TO 60                                   FLU  710
         NDIF = IABS(NUCL(IHI)-NFUD(1,L))                               FLU  720
         IF (NDIF.GE.IDIF) GO TO 60                                     FLU  730
         ICLOS = L                                                      FLU  740
         IDIF = NDIF                                                    FLU  750
   60 CONTINUE                                                          FLU  760
      LL = NFUD(4,ICLOS)                                                FLU  770
      IJK = NFUD(6,ICLOS)                                               FLU  780
      GO TO 90                                                          FLU  790
C    SECOND CASE: NO VARIABLE XSEC                                      FLU  800
   70 DO 80 M=1,3                                                       FLU  810
         IF (LDEF(M).LT.1) GO TO 80                                     FLU  820
         NDIF = IABS(NUCL(IHI)-NDEF(M))                                 FLU  830
         IF (NDIF.GE.IDIF) GO TO 80                                     FLU  840
         IJK = M                                                        FLU  850
         IDIF = NDIF                                                    FLU  860
   80 CONTINUE                                                          FLU  870
      LL = LDEF(IJK)                                                    FLU  880
C    GENERATE MULTIPLIER FPRATN AS RATIO OF (TOTAL UNCONNECTED FISSION  FLU  890
C    RATE)+(FISSION RATE OF CONNECTED ACTINIDE CLOSEST TO THE           FLU  900
C    UNCONNECTED ACTINIDE WITH THE LARGEST FISSION RATE) TO THE         FLU  910
C    PRECEDING QUANTITY IN THE SECOND ().  RETURN FISSION PRODUCT       FLU  920
C    YIELDS TO THEIR ORIGINAL VALUES IF THE IDENTITY OF THE CLOSEST     FLU  930
C    CONNECTED ACTINIDE CHANGES.                                        FLU  940
   90 PRCON = XTEMP(LL)*FISS(LL-ILITE)                                  FLU  950
      FPRATO = FPRATN                                                   FLU  960
      IF (PRCON.LT.ERR) GO TO 110                                       FLU  970
      FPRATN = (XSTOT+PRCON)/PRCON                                      FLU  980
      IF (FPRATN.LT.100000.0) GO TO 100                                 FLU  990
      WRITE (15,99999) NUCL(LL), FPRATN, FPRATO, XTEMP(LL)              FLU 1000
      FPRATN = 1.0                                                      FLU 1010
  100 CONTINUE                                                          FLU 1020
C    RETURN FP YIELDS TO ORIGINAL VALUES IF THE ID OF THE CLOSEST       FLU 1030
C    CONNECTE ACTINIDE HAS CHANGED                                      FLU 1040
  110 CONTINUE                                                          FLU 1050
      IF (IJKO.EQ.IJK .OR. IJKO.EQ.0) GO TO 130                         FLU 1060
      WRITE (16,99998) IJKO, FPRATO                                     FLU 1070
      DO 120 I=1,IFP                                                    FLU 1080
         LOCAT = NFUDFP(I,IJKO)                                         FLU 1090
         IF (LOCAT.GT.0) A(LOCAT) = A(LOCAT)/FPRATO                     FLU 1100
  120 CONTINUE                                                          FLU 1110
      FPRATO = 1.0                                                      FLU 1120
  130 CONTINUE                                                          FLU 1130
      DO 140 I=1,IFP                                                    FLU 1140
         LOCAT = NFUDFP(I,IJK)                                          FLU 1150
         IF (LOCAT.GT.0) A(LOCAT) = A(LOCAT)*FPRATN/FPRATO              FLU 1160
  140 CONTINUE                                                          FLU 1170
      LOCAT = NFUDFP(IFP/3,IJK)                                         FLU 1180
      WRITE (16,99997) NUCL(LL), NUCL(IHI), FPRATN, FPRATO, A(LOCAT)    FLU 1190
      IJKO = IJK                                                        FLU 1200
  150 CONTINUE                                                          FLU 1210
C                                                                       FLU 1220
C    COMPUTE FLUX BASED ON CONDITIONS AT THE START OF THE               FLU 1230
C    INTERVAL                                                           FLU 1240
      DO 160 J=1,ITOT                                                   FLU 1250
         XTEM(J) = XTEMP(J)                                             FLU 1260
  160 CONTINUE                                                          FLU 1270
      XFSM = 0.0                                                        FLU 1280
      FISSR = 0.0                                                       FLU 1290
      IF (IACT.LT.1) GO TO 180                                          FLU 1300
      DO 170 I=1,IACT                                                   FLU 1310
         L = I + ILITE                                                  FLU 1320
         XFSM = XFSM + XTEMP(L)*FISS(I)                                 FLU 1330
         JAM = NUCL(L)/10                                               FLU 1340
         AWT = JAM - (1000*(JAM/1000))                                  FLU 1350
         Z = NUCL(L)/10000                                              FLU 1360
         FISSR = FISSR + XTEMP(L)*FISS(I)*FISSEN(AWT,Z)                 FLU 1370
  170 CONTINUE                                                          FLU 1380
  180 CONTINUE                                                          FLU 1390
      IF (XFSM.GT.ERR) EPF1 = FISSR/XFSM                                FLU 1400
      IF (INDEX.EQ.0) FLUX = POWER*1.0365E+19/FISSR                     FLU 1410
      IF (OLDFLX.GT.0.0) GO TO 190                                      FLU 1420
      FLX = FLUX*1.E-24                                                 FLU 1430
      GO TO 200                                                         FLU 1440
  190 FLX = FLUX/OLDFLX                                                 FLU 1450
C                                                                       FLU 1460
C    MULTIPLY CROSS SECTIONS BY FLUX                                    FLU 1470
C                                                                       FLU 1480
  200 N = 0                                                             FLU 1490
      DO 220 I=1,ITOT                                                   FLU 1500
         D(I) = -(DIS(I)+TOCAP(I)*FLUX*1.0E-24)                         FLU 1510
         N = KD(I)                                                      FLU 1520
         NDIFF = NON0(I) - KD(I)                                        FLU 1530
         IF (NDIFF.EQ.0) GO TO 220                                      FLU 1540
         DO 210 L=1,NDIFF                                               FLU 1550
            N = N + 1                                                   FLU 1560
            A(N) = A(N)*FLX                                             FLU 1570
  210    CONTINUE                                                       FLU 1580
  220 CONTINUE                                                          FLU 1590
C                                                                       FLU 1600
C    ADD CHEMICAL PROCESSING TERM TO DIAGONAL ELEMENT                   FLU 1610
C                                                                       FLU 1620
      IF (MPROS.LE.0) GO TO 260                                         FLU 1630
      DO 250 MP=1,MPROS                                                 FLU 1640
         NOPR = NOPROS(MP)                                              FLU 1650
         DO 240 N=1,NOPR                                                FLU 1660
            NZ = NZPROS(MP,N)                                           FLU 1670
            DO 230 I=1,ITOT                                             FLU 1680
               NZI = NUCL(I)/10000                                      FLU 1690
               IF (NZ.NE.NZI) GO TO 230                                 FLU 1700
               D(I) = D(I) - PRATE(MP)                                  FLU 1710
  230       CONTINUE                                                    FLU 1720
  240    CONTINUE                                                       FLU 1730
  250 CONTINUE                                                          FLU 1740
C                                                                       FLU 1750
C    COMPUTE CONCENTRATIONS AT THE BEGINNING OF THE INTERVAL            FLU 1760
C                                                                       FLU 1770
  260 IF (IACT.LT.1) GO TO 360                                          FLU 1780
      IF (KOUNT.GT.0) GO TO 370                                         FLU 1790
      NB = 0                                                            FLU 1800
      BIG = 0.0                                                         FLU 1810
      IA = ILITE                                                        FLU 1820
      DO 300 I=1,IACT                                                   FLU 1830
         IA = IA + 1                                                    FLU 1840
         DI = -D(IA)                                                    FLU 1850
         IF (DI*DELT.LE.AXN) GO TO 290                                  FLU 1860
         NZ = NON0(IA) - NB                                             FLU 1870
         XNW = B(IA)                                                    FLU 1880
         IF (NZ.LT.1) GO TO 290                                         FLU 1890
         DO 270 K=1,NZ                                                  FLU 1900
            NB = NB + 1                                                 FLU 1910
            J = LOC(NB)                                                 FLU 1920
            XNW = XNW + A(NB)*XTEM(J)                                   FLU 1930
  270    CONTINUE                                                       FLU 1940
         XNW = XNW/DI                                                   FLU 1950
         IF (XNW.LT.1.0E-35) GO TO 280                                  FLU 1960
         ARG = ABS((XNW-XTEM(IA))/XNW)                                  FLU 1970
         IF (ARG.GT.BIG) BIG = ARG                                      FLU 1980
  280    XTEM(IA) = XNW                                                 FLU 1990
  290    NB = NON0(IA)                                                  FLU 2000
  300 CONTINUE                                                          FLU 2010
      IF (BIG.GT.QXN) GO TO 260                                         FLU 2020
C                                                                       FLU 2030
C    COMPUTE FIRST DERIVATIVE FOR TAYLOR SERIES EXPANSION               FLU 2040
C                                                                       FLU 2050
      FDOT = 0.0                                                        FLU 2060
      XFSM = 0.0                                                        FLU 2070
      L = ILITE                                                         FLU 2080
      NUM = 0                                                           FLU 2090
      IF (ILITE.GT.0) NUM = NON0(ILITE)                                 FLU 2100
      DO 330 I=1,IACT                                                   FLU 2110
         L = L + 1                                                      FLU 2120
         NX = NON0(L) - NUM                                             FLU 2130
         XDOTI = D(L)*XTEM(L) + B(L)                                    FLU 2140
         IF (NX.LT.1) GO TO 320                                         FLU 2150
         DO 310 N=1,NX                                                  FLU 2160
            NUM = NUM + 1                                               FLU 2170
            J = LOC(NUM)                                                FLU 2180
            XDOTI = XDOTI + A(NUM)*XTEM(J)                              FLU 2190
  310    CONTINUE                                                       FLU 2200
  320    XDOT(L) = XDOTI                                                FLU 2210
         XFSM = XFSM + XDOTI*FISS(I)                                    FLU 2220
         JAM = NUCL(L)/10                                               FLU 2230
         AWT = JAM - (1000*(JAM/1000))                                  FLU 2240
         Z = NUCL(L)/10000                                              FLU 2250
         FDOT = FDOT + XDOTI*FISS(I)*FISSEN(AWT,Z)                      FLU 2260
  330 CONTINUE                                                          FLU 2270
      IF (ABS(XFSM).GT.ERR) EPF2 = FDOT/XFSM                            FLU 2280
C    COMPUTE SECOND DERIVATIVE FOR TAYLOR SERIES EXPANSION              FLU 2290
      XFSM = 0.0                                                        FLU 2300
      FDDOT = 0.0                                                       FLU 2310
      L = ILITE                                                         FLU 2320
      NUM = 0                                                           FLU 2330
      IF (ILITE.GT.0) NUM = NON0(ILITE)                                 FLU 2340
      DO 350 I=1,IACT                                                   FLU 2350
         L = L + 1                                                      FLU 2360
         NX = NON0(L) - NUM                                             FLU 2370
         XDOTI = D(L)*XDOT(L)                                           FLU 2380
         DO 340 N=1,NX                                                  FLU 2390
            NUM = NUM + 1                                               FLU 2400
            J = LOC(NUM)                                                FLU 2410
            XDOTI = XDOTI + A(NUM)*XDOT(J)                              FLU 2420
  340    CONTINUE                                                       FLU 2430
         XFSM = XFSM + XDOTI*FISS(I)                                    FLU 2440
         JAM = NUCL(L)/10                                               FLU 2450
         AWT = JAM - (1000*(JAM/1000))                                  FLU 2460
         Z = NUCL(L)/10000                                              FLU 2470
         FDDOT = FDDOT + XDOTI*FISS(I)*FISSEN(AWT,Z)                    FLU 2480
  350 CONTINUE                                                          FLU 2490
      IF (ABS(XFSM).GT.ERR) EPF3 = FDDOT/XFSM                           FLU 2500
      IF (INDEX.EQ.1) GO TO 360                                         FLU 2510
C                                                                       FLU 2520
C    CORRECT FLUX USING FIRST TWO TERMS OF SERIES                       FLU 2530
      AVG = 1.0 - 0.50*FDOT*DELT/FISSR                                  FLU 2540
      T3 = (FDDOT-2.*FDOT**2/FISSR)/FISSR*DELT**2*0.16667               FLU 2550
      T1 = 1.0/FISSR                                                    FLU 2560
      T2 = 0.5*DELT*FDOT/FISSR**2                                       FLU 2570
      T4 = (DELT**2/6.0)*(2.0*FDOT**2-FISSR*FDDOT)/FISSR**3             FLU 2580
      T1M = T1*EPF1                                                     FLU 2590
      T2M = T2*EPF1**2/EPF2                                             FLU 2600
      T3M = (DELT**2/6.0)*(2.0*(FDOT/EPF2)**2-(FISSR*FDDOT)/(EPF1*EPF3))FLU 2610
     *     /(FISSR/EPF1)**3                                             FLU 2620
      EPFAVG = (T1M-T2M+T3M)/(T1-T2+T4)                                 FLU 2630
      AVG = AVG - T3                                                    FLU 2640
      DELF = ABS(1.0-AVG)                                               FLU 2650
      NSTEPS = 1. + DELF/0.2                                            FLU 2660
      FLUX = FLUX*AVG                                                   FLU 2670
      WRITE (15,99996) TSEC, DELT, T1, EPF1, FDOT, EPF2, FDDOT, EPF3,   FLU 2680
     *     T2, T4, T1M, T2M, T3M, EPFAVG, FLUX, POWER,T3,AVG,FISSR      FLU 2690
      IF (DELF.GT.0.2) WRITE (15,99995) TSEC, NSTEPS, FLUX, POWER       FLU 2700
      FLX = AVG                                                         FLU 2710
      KOUNT = 1                                                         FLU 2720
C                                                                       FLU 2730
C    GO BACK AND CORRECT CROSS SECTIONS                                 FLU 2740
C                                                                       FLU 2750
      GO TO 200                                                         FLU 2760
  360 T2 = 0.5*DELT*FDOT                                                FLU 2770
      T1 = FISSR                                                        FLU 2780
      T2M = T2/EPF2                                                     FLU 2790
      T1M = FISSR/EPF1                                                  FLU 2800
      T3 = DELT**2*FDDOT/6.0                                            FLU 2810
      T3M = T3/EPF3                                                     FLU 2820
      POWER = FLUX*(FISSR+T2+T3)*9.664E-20                              FLU 2830
      IF (ABS(T1M+T2M+T3M).GT.ERR) EPFAVG = (FISSR+T2+T3)/(T1M+T2M+T3M) FLU 2840
      WRITE (15,99994) TSEC, DELT, T1, EPF1, FDOT, EPF2, FDDOT, EPF3,   FLU 2850
     *     T2, T3, T1M, T2M, T3M, EPFAVG, FLUX, POWER                   FLU 2860
  370 OLDFLX = FLUX                                                     FLU 2870
      OUTFPB(13,2) = FLUX                                               FLU 2880
      OUTFPB(13,3) = POWER                                              FLU 2890
      RETURN                                                            FLU 2900
  380 DO 390 I=1,ITOT                                                   FLU 2910
         D(I) = -DIS(I)                                                 FLU 2920
  390 CONTINUE                                                          FLU 2930
      IF(INDEX.EQ.2)WRITE(15,99993)TSEC,DELT                                    
C                                                                       FLU 2940
C    ADD CHEMICAL PROCESSING TERM TO DIAGONAL ELEMENT                   FLU 2950
C                                                                       FLU 2960
      IF (MPROS.LE.0) GO TO 430                                         FLU 2970
      DO 420 MP=1,MPROS                                                 FLU 2980
         NOPR = NOPROS(MP)                                              FLU 2990
         DO 410 N=1,NOPR                                                FLU 3000
            NZ = NZPROS(MP,N)                                           FLU 3010
            DO 400 I=1,ITOT                                             FLU 3020
               NZI = NUCL(I)/10000                                      FLU 3030
               IF (NZ.NE.NZI) GO TO 400                                 FLU 3040
               D(I) = D(I) - PRATE(MP)                                  FLU 3050
  400       CONTINUE                                                    FLU 3060
  410    CONTINUE                                                       FLU 3070
  420 CONTINUE                                                          FLU 3080
  430 CONTINUE                                                          FLU 3090
      RETURN                                                            FLU 3100
99999 FORMAT (31H FP YIELD ADJUSTMENT RATIO FOR , I7, 3H = , 1PE12.4,   FLU 3110
     *     12H;OLD RATIO =, 1PE12.4, 17H;MACRO FISS XSEC=, 1PE12.4,     FLU 3120
     *     6H;MASS=, 1PE12.4)                                           FLU 3130
99998 FORMAT (21H0THE FP YIELDS IN THE, I2, 24H COLUMN OF NFUDFP HAVE B,FLU 3140
     *     56HEEN RETURNED TO THEIR ORIGINAL VALUES USING A DIVISOR OF, FLU 3150
     *     F10.7)                                                       FLU 3160
99997 FORMAT (54H0FP YIELD ADJUSTMENT FOR UNCONNECTED ACTINIDES: CONNEC,FLU 3170
     *     8HTED ACT=, I6, /, 25H LARGEST UNCONNECTED ACT=, I6, /,      FLU 3180
     *     12H NEW RATIO =, F10.7, /, 14H   OLD RATIO =, F10.7, /,      FLU 3190
     *     18H   FP YIELD IN A =, 1PE12.5)                              FLU 3200
99996 FORMAT (6H0TSEC=, 1PE10.3, 6H DELT=, 1PE10.3, 4H T1=, 1PE10.3,    FLU 3210
     *     6H EPF1=, 1PE10.3, 6H FDOT=, 1PE10.3, 6H EPF2=, 1PE10.3,     FLU 3220
     *     7H FDDOT=, 1PE10.3, 6H EPF3=, 1PE10.3, /, 4H T2=, 1PE10.3,   FLU 3230
     *     4H T4=, 1PE10.3, 5H T1M=, 1PE10.3, 5H T2M=, 1PE10.3, 5H T3M=,FLU 3240
     *     1PE10.3, 8H EPFAVG=, 1PE10.3, 6H FLUX=, 1PE10.3, 7H POWER=,  FLU 3250
     *     1PE10.3,/,4H T3=,1PE10.3,4HAVG=,1PE10.3,6HFISSR=,1PE10.3)    FLU 3260
99995 FORMAT (26H0**** TIME STEP ENDING AT , 1PE9.2, 15H S IS TOO LARGE,FLU 3270
     *     14H.  DIVIDE INTO, I5, 19H INCREMENTS.  FLUX=, 1PE9.2,       FLU 3280
     *     9H; SP POW=, 1PE9.2)                                         FLU 3290
99994 FORMAT (6H0TSEC=, 1PE10.3, 6H DELT=, 1PE10.3, 4H T1=, 1PE10.3,    FLU 3300
     *     6H EPF1=, 1PE10.3, 6H FDOT=, 1PE10.3, 6H EPF2=, 1PE10.3,     FLU 3310
     *     7H FDDOT=, 1PE10.3, 6H EPF3=, 1PE10.3, /, 4H T2=, 1PE10.3,   FLU 3320
     *     4H T3=, 1PE10.3, 5H T1M=, 1PE10.3, 5H T2M=, 1PE10.3, 5H T3M=,FLU 3330
     *     1PE10.3, 8H EPFAVG=, 1PE10.3, 6H FLUX=, 1PE10.3, 7H POWER=,  FLU 3340
     *     1PE10.3)                                                     FLU 3350
99993 FORMAT(6H0TSEC=,1PE10.3,6H DELT=,1PE10.3)                                 
      END                                                               FLU 3360
      SUBROUTINE FUDGE(NFUDFP, IVEC, IFD, XSTORE, LX, XNEW, MX, TOCAP,  FUD   10
     *     A, FISS, MIRR, DELT)                                         FUD   20
      INTEGER*2 NFUDFP                                                  FUD   30
      DIMENSION NFUDFP(IFD,1), XSTORE(LX,1), XNEW(MX,1), TOCAP(1),      FUD   40
     *     A(1), FISS(1), F(32)                                         FUD   50
      COMMON /FUDG/ NFUD(7,20), FUDSIG(32,20), CBUP(32), NXSEC, NSM,    FUD   60
     *     NPT, NCO(2,20), FPBST(050,3), VECT(050), FOLD(20), FNEW(20)  FUD   70
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  FUD   80
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           FUD   90
     *     TFPBAV(4), IPHMAX                                            FUD  100
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    FUD  110
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  FUD  120
     *     ITITLE(18)                                                   FUD  130
C      LEVEL 2,XSTORE,XNEW,TOCAP,FISS,A,NFUDFP                          FUD  140
      ILOW = ILITE + 1                                                  FUD  150
      VECTOT = 0.0                                                      FUD  160
      DO 10 N=1,NPT                                                     FUD  170
         F(N) = 0.0                                                     FUD  180
   10 CONTINUE                                                          FUD  190
      IF (MIRR.LT.1) F(1) = 1.0                                         FUD  200
      XF = 0.0                                                          FUD  210
      CBUP(NPT+1) = CBUP(NPT)                                           FUD  220
      TIRR = 0.0                                                        FUD  230
      ANTIC = 0.0                                                       FUD  240
      LOCAT = 1                                                         FUD  250
      LOCFIS = 1                                                        FUD  260
      KIRR = MIRR + 1                                                   FUD  270
      X = 0.0                                                           FUD  280
      IHI = ILITE + IACT                                                FUD  290
C    COMPUTE INITIAL HM MASS(=ACT)    FROM PREVIOUS VECTOR              FUD  300
      IF (IVEC.LT.0) GO TO 30                                           FUD  310
      DO 20 I=ILOW,IHI                                                  FUD  320
         VECTOT = VECTOT + XNEW(IVEC,I)                                 FUD  330
   20 CONTINUE                                                          FUD  340
      GO TO 50                                                          FUD  350
   30 DO 40 I=ILOW,IHI                                                  FUD  360
         VECTOT = VECTOT + XSTORE(-IVEC,I)                              FUD  370
   40 CONTINUE                                                          FUD  380
   50 CONTINUE                                                          FUD  390
      VECT(KIRR) = VECTOT                                               FUD  400
C    DETERMINE BURNUP IN MWD/G-A ACTUAL HEAVY METAL                     FUD  410
      IF (MIRR.EQ.0) GO TO 130                                          FUD  420
      DO 60 M=1,MIRR                                                    FUD  430
         IF (FPBST(M,2).GT.ERR) TIRR = TIRR + FPBST(M,1)                FUD  440
         IF (VECT(M).GT.ERR) X = X + FPBST(M,3)*FPBST(M,1)/(86400.0*    FUD  450
     *        VECT(M))                                                  FUD  460
   60 CONTINUE                                                          FUD  470
C    CALCULATE BURNUP ANTICIPATOR FROM SPECIFIED POWER                  FUD  480
      IF (POWER.LT.ERR) GO TO 70                                        FUD  490
      IF (VECT(KIRR).GT.ERR) XAN = POWER*DELT/(86400.0*VECT(KIRR))      FUD  500
      IF (X.GT.ERR) ANTIC = (XAN+X)/X                                   FUD  510
      GO TO 80                                                          FUD  520
C    CALCULATE BURNUP ANTICIPATOR FROM TIME                             FUD  530
   70 IF (TIRR.GT.ERR) ANTIC = (DELT+TIRR)/TIRR                         FUD  540
   80 CONTINUE                                                          FUD  550
      XF = X*ANTIC                                                      FUD  560
      DX = XF - X                                                       FUD  570
      DXX = DX                                                          FUD  580
C    CALCULATE XSEC WEIGHTING FACTORS                                   FUD  590
      NZZ = 0                                                           FUD  600
      DO 100 N=1,NPT                                                    FUD  610
         IF (X.GT.CBUP(N+1)) GO TO 100                                  FUD  620
         IF (N.EQ.NPT) GO TO 110                                        FUD  630
C    ANTICIPATED BURNUP .GE. WIDTH OF BURNUP GROUP                      FUD  640
         IF (DXX.LT.1.0E-05) GO TO 120                                  FUD  650
         DCB = CBUP(N+1) - CBUP(N)                                      FUD  660
         NZZ = NZZ + 1                                                  FUD  670
         IF (NZZ.EQ.1) DCB = CBUP(N+1) - X                              FUD  680
         IF (DCB.LT.ERR) GO TO 100                                      FUD  690
         IF (DXX.LE.DCB) GO TO 90                                       FUD  700
         F(N) = DCB/DX                                                  FUD  710
         DXX = DXX - DCB                                                FUD  720
         GO TO 100                                                      FUD  730
   90    CONTINUE                                                       FUD  740
C    ANTICIPATED BURNUP .LT. WIDTH OF BURNUP GROUP                      FUD  750
         F(N) = DXX/DX                                                  FUD  760
         DXX = 0.0                                                      FUD  770
  100 CONTINUE                                                          FUD  780
  110 CONTINUE                                                          FUD  790
C    LAST XSEC FACTOR                                                   FUD  800
      IF (DXX.GT.1.0E-05) F(NPT) = DXX/DX                               FUD  810
  120 CONTINUE                                                          FUD  820
  130 CONTINUE                                                          FUD  830
      WRITE (16,99998)                                                  FUD  840
      WRITE (16,99997) ITITLE, NREC                                     FUD  850
      WRITE (16,99996) BASIS                                            FUD  860
      WRITE (16,99995) IVEC, VECTOT, X, XF, ANTIC                       FUD  870
      WRITE (16,99999) (N,F(N),N=1,NPT)                                 FUD  880
      WRITE (16,99994)                                                  FUD  890
      X = X*ANTIC                                                       FUD  900
C    BEGIN LOOP OVER THE XSEC TO BE ALTERED                             FUD  910
      DO 180 L=1,NXSEC                                                  FUD  920
         I = NFUD(4,L)                                                  FUD  930
C    CALCULATE WEIGHTED XSEC                                            FUD  940
         Y = 0.0                                                        FUD  950
         DO 140 N=1,NPT                                                 FUD  960
            Y = Y + FUDSIG(N,L)*F(N)                                    FUD  970
  140    CONTINUE                                                       FUD  980
         IF (Y.EQ.0.0) GO TO 180                                        FUD  990
C    STORE OLD XSEC                                                     FUD 1000
         FOLD(L) = FNEW(L)                                              FUD 1010
C    STORE NEWLY CALCULATED XSEC                                        FUD 1020
         FNEW(L) = Y                                                    FUD 1030
C    CORRECT TOCAP(I)                                                   FUD 1040
         TOCAP(I) = TOCAP(I) - FOLD(L) + FNEW(L)                        FUD 1050
         LOCFIS = I - ILITE                                             FUD 1060
C    CORRECT VARIABLE XSEC IN ARRAY A IF VAR XSEC IS NOT FISSION        FUD 1070
         IF (NFUD(2,L).NE.4) A(NFUD(5,L)) = A(NFUD(5,L))*FNEW(L)/FOLD(L)FUD 1080
C    CORRECT FISSION PRODUCT YIELDS IF THIS IS A FISSION XSEC AND IF    FUD 1090
C    THIS ACTINIDE HAS YIELDS.  ALSO CORRECT FISS(I) FOR ALL FISSION    FUD 1100
C    XSEC                                                               FUD 1110
         IF (NFUD(2,L).NE.4) GO TO 170                                  FUD 1120
         IF (NFUD(6,L).LE.0) GO TO 160                                  FUD 1130
         DO 150 N=1,IFP                                                 FUD 1140
            MOCAT = NFUDFP(N,NFUD(6,L))                                 FUD 1150
            IF (MOCAT.LT.1) GO TO 150                                   FUD 1160
            LOCAT = MOCAT                                               FUD 1170
            A(LOCAT) = A(LOCAT)*FNEW(L)/FOLD(L)                         FUD 1180
  150    CONTINUE                                                       FUD 1190
  160    CONTINUE                                                       FUD 1200
         FISS(LOCFIS) = FNEW(L)                                         FUD 1210
  170    CONTINUE                                                       FUD 1220
C    PRINT NEW CROSS-SECTIONS CALCULATED FROM CURVE-FIT                 FUD 1230
         NFU5 = NFUD(5,L)                                               FUD 1240
         WRITE (16,99993) L, NFUD(1,L), NFUD(2,L), (NFUD(J,L),J=4,6),   FUD 1250
     *        LOCFIS, A(NFU5), TOCAP(I), A(LOCAT), FISS(LOCFIS),        FUD 1260
     *        FOLD(L), FNEW(L)                                          FUD 1270
  180 CONTINUE                                                          FUD 1280
      RETURN                                                            FUD 1290
99999 FORMAT ((1H , 7(2HN=, I2, 1X, 2HF=, 1PE10.3, 1X), /))             FUD 1300
99998 FORMAT (1H-)                                                      FUD 1310
99997 FORMAT (1H0, 10X, 18A4, 10X, 12H RECYCLE # =, I3, 5X, 8HUNIT=16 ) FUD 1320
99996 FORMAT (1H0, 15X, 7HBASIS =, 18A4)                                FUD 1330
99995 FORMAT (17H0INITIAL VECTOR =, I4, 5X, 16HTOTAL ACT G-A = ,        FUD 1340
     *     1PE10.3, 5X, 10HMWD/G-A = , 1PE10.3, 1X, 3HTO , 1PE10.3, 2X, FUD 1350
     *     21HANTICIPATION FACTOR= , 1PE10.3, /)                        FUD 1360
99994 FORMAT (54H0L  NUCLID  XSEC       TOCAP(I)  A(N)  FP YIELD   FISS,FUD 1370
     *     57H(J)   A(N)     TOCAP(I)      A(N)       FISS(J)      OLD ,FUD 1380
     *     10H       NEW, /, 1H , 11X, 4HTYPE, 10X, 2HI=, 6X, 2HN=, 3X, FUD 1390
     *     9HINDIC ARR, 33X, 8HFP YIELD, 18X, 4HXSEC, 07X, 4HXSEC)      FUD 1400
99993 FORMAT (1H , I2, 1X, I6, 3X, I2, 10X, I4, 4X, I4, 4X, I3, 4X, I4, FUD 1410
     *     2X, 6(2X, 1PE10.3))                                          FUD 1420
      END                                                               FUD 1430
      SUBROUTINE DECAY(T, A, D, LOC, NON0, KD, LONG, XP, XPAR, XTEMP,   DEC   10
     *     NQ, XNEW, MX, IOUT, NUCL)                                    DEC   20
      LOGICAL LONG                                                      DEC   30
      INTEGER*2 LOC,NON0,KD,NQ                                          DEC   40
      DOUBLE PRECISION BATE                                             DEC   50
      DIMENSION A(1), D(1), LOC(1), NON0(1), KD(1), LONG(1), XP(1),     DEC   60
     *     XPAR(1), XTEMP(1), NQ(1), XNEW(MX,1), NUCL(1)                DEC   70
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  DEC   80
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           DEC   90
     *     TFPBAV(4), IPHMAX                                            DEC  100
      COMMON /TERMD/ QUB(50), DD(100), DXP(100), QUEUE(50), NQUEUE(50), DEC  110
     *     NQU(50)                                                      DEC  120
C        LEVEL 2,D,LOC,NON0,KD,XP,XPAR,XTEMP,NQ,XNEW,A,NUCL             DEC  130
C                                                                       DEC  140
C       DECAY TREATS SHORT-LIVED NUCLIDES AT BEGINNING OF CHAINS USING  DEC  150
C       BATEMAN EQUATIONS                                               DEC  160
C                                                                       DEC  170
C      SCHEME WORKS BY SAVING ALL SECOND PARENTS OF THE FIRST           DEC  180
C      FIRST-PARENT WHILE CALCULATING THE CONTRIBUTION OF THE SECOND    DEC  190
C      PARENTS.  IT THEN EXAMINES THE FIRST SECOND PARENT (WHICH APPEARSDEC  200
C      IN THE ROUTINE AS A FIRST PARENT WITH INDEX J) TO SEE IF IT, IN  DEC  210
C      TURN, HAS PARENTS.  THIS SEQUENCE IS PURSUED AWAY FROM THE       DEC  220
C      ULTIMATE DAUGHTER L UNTIL THERE ARE NO MORE PARENTS.  THE SCHEME DEC  230
C      THEN RETURNS AND PURSUES THE SECOND SECOND-PARENT IN A SIMILAR   DEC  240
C      MANNER, THEN THE THIRD SECOND-PARENT, FOURTH, ETC.  AFTER ALL    DEC  250
C      SECOND PARENTS ARE EXHAUSTED, NSAVE IS BACK TO ZERO, AND THE     DEC  260
C      SECOND FIRST-PARENT (INDEX N) IS IDENTIFIED BY INCREMENTING      DEC  270
C      N AND THEN APPLYING THE SCHEME.                                  DEC  280
C                                                                       DEC  290
C      L= NUCLIDE # OF ULTIMATE DAUGHTER                                DEC  300
C      NUL= # OF START OF REACTION RATES OF FIRST PARENTS OF L IN       DEC  310
C           MATRIX A                                                    DEC  320
C      NUM= # OF END OF REACTION RATES OF FIRST PARENTS OF L IN         DEC  330
C           MATRIX A                                                    DEC  340
C      N= RUNNING INDEX OVER THE REACTION RATES OF THE FIRST PARENTS    DEC  350
C         OF L IN A                                                     DEC  360
C      J= INITIALLY, THE NUCLIDE # OF A PARTICULAR FIRST PARENT OF L    DEC  370
C         VIA REACTION N; THROUGH THE ACTIONS OF THE QUEUE, J BECOMES   DEC  380
C         THE NUCLIDE # OF FURTHER-REMOVED PARENTS OF L UNTIL THE END   DEC  390
C         OF THE CHAIN IS REACHED; J IS ALWAYS THE NEXT-TO-THE          DEC  400
C         FARTHEST NUCLIDE FROM L; EVENTUALLY, J IS BACK TO ITS         DEC  410
C         ORIGINAL VALUE, NSAVE IS THEN ZERO AGAIN, AND REACTION N      DEC  420
C         IS FINISHED;  N IS THEN INCREMENTED TO THE NEXT REACTION.     DEC  430
C      K= RUNNING INDEX OVER THE FARTHEST-REMOVED REACTION RATES        DEC  440
C         (AND THEREBY PARENTS) OF NUCLIDE L                            DEC  450
C      J1= NUCLIDE # OF THE FARTHEST-REMOVED PARENT OF L WHICH IS       DEC  460
C          BEING CONSIDERED AT THE PRESENT TIME                         DEC  470
C      NUF= # OF THE START OF RECTION RATES OF THE SECOND PARENTS       DEC  480
C           OF L (N-PRODUCING) IN MATRIX A                              DEC  490
C      NUX= # OF THE END OF THE REACTION RATES OF THE SECOND PARENTS    DEC  500
C           OF L (N-PRODUCING) IN MATRIX A                              DEC  510
C      NQUEUE, QUEUE, NQU= IDENTITIES OF LAST (FARTHEST REMOVED FROM    DEC  520
C           L) NUCLIDE, NEXT-TO-LAST NUCLIDE, AND FRACTIONAL            DEC  530
C           CONTRIBUTION OF LAST NUCLIDE TO L; THESE ARE THE YET-TO-BE  DEC  540
C           PURSUED CHAINS                                              DEC  550
C    D(I)= TOTAL REMOVAL RATE OF NUCLIDE I                              DEC  560
C    T= LENGTH OF TIME STEP IN SECONDS                                  DEC  570
C    XP(I)= EXP(D(I)*T)                                                 DEC  580
C    LONG(I)= TRUE IF D(I)*T.LE.6.9, WHICH IS THE SAME AS HAVING THE    DEC  590
C       "REMOVAL HALF-LIFE".LE.0.1*T.                                   DEC  600
C    NON0(L)= # OF LAST REACTION IN A THAT IS A PARENT OF NUCLIDE L     DEC  610
C    KD(L)= SAME AS NON0(L) EXCEPT FLUX-DEPENDENT REACTIONS ARE         DEC  620
C           EXCLUDED                                                    DEC  630
C    LOC(M)= NUCLIDE # OF PARENT OF NUCLIDE L FOR M=NON0(L-1)+1 TO      DEC  640
C    NON0(L)                                                            DEC  650
C      AKDJQ= RATIO OF PRODUCTION OF NUCLIDE NB FROM NUCLIDE NA DIVIDED DEC  660
C             BY THE TOTAL REMOVAL RATE OF NUCLIDE NA                   DEC  670
C    A(K)= MATRIX OF REACTION RATES; THOSE PRODUCING NUCLIDE L ARE      DEC  680
C       K=NON0(L-1)+1 TO NON0(L)                                        DEC  690
C                                                                       DEC  700
C    INITIALIZE ARRAYS                                                  DEC  710
C      TEST D(I)*T AND SET UP ARRAY LONG(I) FOR USE IN DECAY, TERM,     DEC  720
C      AND EQUIL                                                        DEC  730
C    STORE EXP(D(I)*T) IN XP(I)                                         DEC  740
      DO 20 I=1,ITOT                                                    DEC  750
         XPAR(I) = 0.0                                                  DEC  760
         LONG(I) = .FALSE.                                              DEC  770
         XPI = 0.0                                                      DEC  780
         DT = D(I)*T                                                    DEC  790
C    SET EXP(D(I)*T)=0.0 IF REMOVAL LIFE IS LESS THAN 1.4% OF TIME STEP DEC  800
         IF (DT.LT.-50.) GO TO 10                                       DEC  810
         IF (ABS(DT).LE.AXN) LONG(I) = .TRUE.                           DEC  820
         XPI = EXP(DT)                                                  DEC  830
   10    XP(I) = XPI                                                    DEC  840
   20 CONTINUE                                                          DEC  850
      NUL = 1                                                           DEC  860
C    LOOP OVER THE DAUGHTER NUCLIDE TO BE CONSIDERED                    DEC  870
      DO 160 L=1,ITOT                                                   DEC  880
         XTEM = 0.0                                                     DEC  890
         DL = -D(L)                                                     DEC  900
C    USE ALL REACTIONS PRODUCING L                                      DEC  910
         NUM = NON0(L)                                                  DEC  920
C    USE ONLY FLUX-INDEPENDENT REACTIONS PRODUCING L                    DEC  930
         IF (INDEX.EQ.2) NUM = KD(L)                                    DEC  940
C    SKIP IF NO REACTIONS PRODUCING L                                   DEC  950
         IF (NUM.LT.NUL) GO TO 150                                      DEC  960
C    LOOP OVER THE FIRST PARENTS OF L                                   DEC  970
         DO 140 N=NUL,NUM                                               DEC  980
            J = LOC(N)                                                  DEC  990
            DJ = -D(J)                                                  DEC 1000
C    SKIP PARENT LOC(N) IF IT IS LONG-LIVED                             DEC 1010
            IF (LONG(J)) GO TO 140                                      DEC 1020
C    CALCULATE CONTRIBUTION OF FIRST PARENT LOC(N) TO NUCLIDE L         DEC 1030
C    USE THIS FORM IF HALF-LIVES ARE NEARLY EQUAL                       DEC 1040
            IF (ABS(DL/DJ-1.0).LE.1.0E-5) XTEM = XTEM +                 DEC 1050
     *           XTEMP(J)*A(N)*XP(J)*T                                  DEC 1060
C    USE THIS FORM IS HALF-LIVES ARE NOT NEARLY EQUAL                   DEC 1070
            IF (ABS(DL/DJ-1.0).GT.1.0E-5) XTEM = XTEM +                 DEC 1080
     *           XTEMP(J)*A(N)*(XP(J)-XP(L))/(DL-DJ)                    DEC 1090
C    FRACTION OF LOC(N) GOING TO L                                      DEC 1100
            QUE = A(N)/DJ                                               DEC 1110
C    NQ(L)=0 INDICATES ULTIMATE DAUGHTER (BEGINNING OF CHAIN)           DEC 1120
            NQ(L) = 0                                                   DEC 1130
C    NQ(J=LOC(N))=L CONNECTS CURRENT FIRST PARENT TO DAUGHTER           DEC 1140
            NQ(J) = L                                                   DEC 1150
            NSAVE = 0                                                   DEC 1160
C      LABEL 20 IS THE BEGINNING OF A LOOP ENDING JUST ABOVE LABEL 140. DEC 1170
C      THIS LOOP SHIFTS THE IDENTITY OF THE FARTHEST-REMOVED PARENT AND DEC 1180
C      NEXT-TO-THE-FARTHEST-REMOVED PARENT (RELATIVE TO DAUGHTER L AND  DEC 1190
C      ITS PARENT LOC(N)) VIA CHANGING THE REACTION INDICES J AND J1.   DEC 1200
C      THIS RESULTS IN A CHAIN BEING PRUSUED UNTIL EITHER (L) IT HAS    DEC 1210
C      NO PRECURSORS OR (2) A LONG-LIVED PRECURSOR IS ENCOUNTERED, IN   DEC 1220
C      WHICH CASE THE CHAIN IS SKIPPED.  THE CHAINS THAT NEED TO BE     DEC 1230
C      PURSUED ARE STORED FOR FUTURE ACCESS.  NOTE THAT THIS METHOD     DEC 1240
C      RESULTS IN THE IDENTITY OF NUCLIDES J AND J1 (AND THEIR DISTANCE DEC 1250
C      FROM L) VARYING WITH EACH PASS THROUGH LABEL 20.                 DEC 1260
C                                                                       DEC 1270
C    SET UP REACTION RATE INDICES FOR PARENTS OF CURRENT NUCLIDE J      DEC 1280
   30       NUX = NON0(J)                                               DEC 1290
            IF (INDEX.EQ.2) NUX = KD(J)                                 DEC 1300
            NUF = 1                                                     DEC 1310
            IF (J.GT.1) NUF = NON0(J-1) + 1                             DEC 1320
C    SKIP IF NUCLIDE HAS NO PRECURSORS AND "WORK OFF" ANOTHER SAVED     DEC 1330
C    CHAIN                                                              DEC 1340
            IF (NUF.GT.NUX) GO TO 130                                   DEC 1350
C    LOOP OVER REACTIONS PRODUCING CURRENT NUCLIDE J                    DEC 1360
            DO 120 K=NUF,NUX                                            DEC 1370
               J1 = LOC(K)                                              DEC 1380
C    SKIP THIS CHAIN IF CURRENT FURTHEST PARENT NUCLIDE J1 IS LONG-LIVEDDEC 1390
               IF (LONG(J1)) GO TO 120                                  DEC 1400
               KP = J                                                   DEC 1410
C    SKIP THIS CHAIN IF NUCLIDE J1 IS ALREADY IN THE CHAIN; I.E.,       DEC 1420
C    IF THE NUCLIDE PRODUCES ITSELF.  THE ENTIRE CHAIN (OR QUEUE)       DEC 1430
C    AS STORED IN NQ IS SEARCHED BY THE FOLLOWING 3 CARDS               DEC 1440
   40          IF (J1.EQ.NQ(KP)) GO TO 120                              DEC 1450
               KP = NQ(KP)                                              DEC 1460
               IF (KP.NE.0) GO TO 40                                    DEC 1470
               DJ = -D(J1)                                              DEC 1480
C      CALCULATE THE FRACTION OF NUCLIDE J1 REMOVALS THAT               DEC 1490
C      PRODUCE NUCLIDE L; END CHAIN IS .LT. 1.0E-06                     DEC 1500
C    CHAIN IF .LT. 1 PPM                                                DEC 1510
               AKDJQ = A(K)/DJ*QUE                                      DEC 1520
               IF (AKDJQ.LE.1.0E-06) GO TO 120                          DEC 1530
C    THE FOLLOWING 10 CARDS STORE THE REMVAL RATES AND EXPONENTIALS     DEC 1540
C    THEREOF FOR THE CURRENT CHAIN IN PACKED ORDER IN ARRAYS DD(I)      DEC 1550
C    AND DXP(I), RESPECTIVELY.  IF THE CHAIN IS LONGER THAN 100         DEC 1560
C    NUCLIDES A MESSAGE IS PRINTED AND THE CHAIN IS SKIPPED.            DEC 1570
C                                                                       DEC 1580
C    SUCCESSFUL COMPLETION BRANCHES TO STATEMENT LABEL 50               DEC 1590
               NQ(J1) = J                                               DEC 1600
               I = 1                                                    DEC 1610
               KP = J1                                                  DEC 1620
   50          DD(I) = -D(KP)                                           DEC 1630
               DXP(I) = XP(KP)                                          DEC 1640
               KP = NQ(KP)                                              DEC 1650
               IF (KP.EQ.0) GO TO 60                                    DEC 1660
               I = I + 1                                                DEC 1670
               IF (I.LE.100) GO TO 50                                   DEC 1680
               WRITE (15,99999) L, J1, J, AKDJQ                         DEC 1690
               GO TO 130                                                DEC 1700
   60          BATE = 0.                                                DEC 1710
               I1 = I - 1                                               DEC 1720
               XPL = XP(L)                                              DEC 1730
C      D R VONDY FORM OF BATEMAN EQUATIONS -- ORNL-TM-361               DEC 1740
C                                                                       DEC 1750
C    THE LOOP DOWN TO LABEL 100 CALCULATES THE CONTRIBUTION OF          DEC 1760
C    NUCLIDE J1 TO L.  FOR MOST NUCLIDES, EQN. 9 IN ORNL-4628           DEC 1770
C    (BELL), WHICH IS THE SAME AS EQN. 3.4 IN ORNL/TM-361, IS           DEC 1780
C    USED.  IF THE "REMOVAL LIVES" OF NUCLIDE J1 AND ANY OTHER          DEC 1790
C    NUCLIDE ARE CLOSE (0.01%), A SPECIAL EQUATION IS USED              DEC 1800
C    (EQN. 10 IN ORNL-4628 EXCEPT AS NOTED).                            DEC 1810
               DO 110 KB=1,I1                                           DEC 1820
                  XPJ = DXP(KB)                                         DEC 1830
C    SKIP THIS NUCLIDE IF ITS REMOVAL RATE PLUS THAT OF THE             DEC 1840
C    ULTIMATE DAUGHTER ARE SMALL; I.E., IF BOTH THE NUCLIDE             DEC 1850
C    AND ITS DAUGHTER HAVE REMOVAL LIVES LESS THAN 1.4% OF              DEC 1860
C    THE TIME STEP LENGTH                                               DEC 1870
                  IF (XPL+XPJ.LT.ERR) GO TO 110                         DEC 1880
                  DK = DD(KB)                                           DEC 1890
                  PROD = (DL/DK-1.0)                                    DEC 1900
C    SAVE ORIGINAL PROD AS DKR                                          DEC 1910
                  DKR = PROD                                            DEC 1920
C    FOR THE CASE WHEN THE REMOVAL LIVES OF THE CURRENT NUCLIDE         DEC 1930
C    AND THE ULTIMATE DAUGHTER ARE VERY NEARLY EQUAL, USE               DEC 1940
C    EQN. 4.2 IN ORNL/TM-361 TRUNCATED TO TWO TERMS OF THE              DEC 1950
C    SERIES                                                             DEC 1960
                  IF (ABS(PROD).GT.1.E-4) GO TO 70                      DEC 1970
                  PROD = T*DK*XPJ*(1.0-0.5*(DL-DK)*T)                   DEC 1980
                  GO TO 80                                              DEC 1990
   70             PROD = (XPJ-XPL)/PROD                                 DEC 2000
                  PRO1 = XPJ/DKR                                        DEC 2010
   80             PI = 1.0                                              DEC 2020
                  S1 = 2./(DK*T)                                        DEC 2030
C    THE LOOP DOWN TO LABEL 90 CALCULATES THE MULTIPLE PRODUCT          DEC 2040
C    WITHIN THE BRACKET IN EQN. 9 OF ORNL-4628.  THE RESULT             DEC 2050
C    IS STORED AS PI.  THE CURRENT TERM OF THE PRODUCT IS S.            DEC 2060
                  DO 100 JK=1,I1                                        DEC 2070
                     IF (JK.EQ.KB) GO TO 100                            DEC 2080
                     S = 1.0 - DK/DD(JK)                                DEC 2090
C    IF S IS LT 0.01%, THE TWO CURRENT NUCLIDES HAVE                    DEC 2100
C    CLOSE REMOVAL LIVES AND A SPECIAL FORM OF THE                      DEC 2110
C    BATEMAN EQUATION IS USED                                           DEC 2120
                     IF (ABS(S).GT.1.E-4) GO TO 90                      DEC 2130
C    PRO1 IS SUBSTITUTED FOR PROD IF THE TWO CURRENT NUCLIDES AND       DEC 2140
C    THE ULTIMATE DAUGHTER ALL HAVE REMOVAL LIVES WITHIN 0.01%          DEC 2150
                     IF (ABS(DKR).GT.1.0E-4) PROD = PRO1                DEC 2160
                     S = S1                                             DEC 2170
   90                PI = PI*S                                          DEC 2180
                     IF (ABS(PI).GT.1.E25) GO TO 110                    DEC 2190
  100             CONTINUE                                              DEC 2200
C    BATE ACCUMULATES THE TERMS OF THE BRACKET SUMMATION OF EQN. 9.     DEC 2210
                  BATE = BATE + PROD/PI                                 DEC 2220
  110          CONTINUE                                                 DEC 2230
               IF (BATE.LT.0.0) WRITE (15,99998) L, I, BATE, XTEM,      DEC 2240
     *              XTEMP(J1), AKDJQ                                    DEC 2250
               IF (BATE.LT.0.) BATE = 0.                                DEC 2260
C    XTEM ACCUMULATES THE TERMS OF THE OUTERMOST SUMMATION IN EQN. 9.   DEC 2270
               XTEM = XTEM + XTEMP(J1)*AKDJQ*BATE                       DEC 2280
C    STORE A MAX OF 50 TO-BE-PURSUED EQUATIONS.                         DEC 2290
               IF (NSAVE.GE.50) GO TO 120                               DEC 2300
C    THE FOLLOWING 4 CARDS STORE THE CURRENT FARTHEST-REMOVED           DEC 2310
C    PARENT NUCLIDE, ITS IMMEDIATE DAUGHTER J, AND ITS CONTRI           DEC 2320
C    BUTION TO THE ULTIMATE DAUGHTER L SO THAT J1 CAN BE                DEC 2330
C    CHECKED FOR PRECURSORS AND POSSIBLE FURTHER PURSUIT AT A           DEC 2340
C    LATER TIME                                                         DEC 2350
               NSAVE = NSAVE + 1                                        DEC 2360
               NQUEUE(NSAVE) = J1                                       DEC 2370
               QUEUE(NSAVE) = AKDJQ                                     DEC 2380
               NQU(NSAVE) = J                                           DEC 2390
  120       CONTINUE                                                    DEC 2400
  130       IF (NSAVE.LE.0) GO TO 140                                   DEC 2410
C    IF NSAVE.GT.0, THEN NUCLIDES REMAIN TO BE CHECKED FOR              DEC 2420
C    PRECURSORS AND THEIR CONTRIBUTIONS TO THE DAUGHTER                 DEC 2430
C    INCLUDED IF THEY HAVE ONE.  THE FOLLOWING 4 CARDS                  DEC 2440
C    RETRIEVE THE NEEDED DATA ON THE LAST-SAVED NUCLIDE.                DEC 2450
C    NOTE THAT THE NUCLIDE J WAS J1 WHEN IT WAS ORIGINALLY              DEC 2460
C    STORED IN THE ARRAYS.                                              DEC 2470
            J = NQUEUE(NSAVE)                                           DEC 2480
            QUE = QUEUE(NSAVE)                                          DEC 2490
            NQ(J) = NQU(NSAVE)                                          DEC 2500
            NSAVE = NSAVE - 1                                           DEC 2510
C    BRANCH TO LABEL 20 TO PURSUE NUCLIDE J                             DEC 2520
            GO TO 30                                                    DEC 2530
  140    CONTINUE                                                       DEC 2540
C    IF NUCLIDE L IS LONG-LIVED, DIVIDE BY XP(L).  THIS WILL BE         DEC 2550
C    REVERSED IN EQUIL, WHICH WILL HANDLE THE LONG-LIVED NUCLIDES       DEC 2560
C    AND THEIR DAUGHTERS.                                               DEC 2570
         IF (LONG(L)) XPAR(L) = XTEM/XP(L)                              DEC 2580
  150    NUL = NON0(L) + 1                                              DEC 2590
C    FOR SHORT-LIVED NUCLIDE L, ADD AMOUNT OF NUCLIDE L PRESENT         DEC 2600
C    AT BEGINNING OF TIME STEP THAT HAS NOT DECAYES AWAY                DEC 2610
C    (=XTEMP*XP) TO CONTRIBUTIONS FROM PRECURSORS AND STORE             DEC 2620
C    PERMANENTLY IN VECTOR 13.                                          DEC 2630
         IF (.NOT.LONG(L)) XNEW(13,L) = XTEM + XTEMP(L)*XP(L)           DEC 2640
  160 CONTINUE                                                          DEC 2650
      DO 170 I=1,ITOT                                                   DEC 2660
C    ADD THE CURRENT CONTRIBUTION TO THE LONG-LIVED NUCLIDES TO THE     DEC 2670
C    AMOUNT AT THE BEGINNING OF THE STEP IN XTEMP.  ZERO OUT THE        DEC 2680
C    SHORT-LIVED NUCLIDES SINCE THEY WILL HAVE GONE THROUGH AT LEAST    DEC 2690
C    10 HALF-LIVES DURING THE TIME STEP AND WILL BE INSIGNIFICANT.      DEC 2700
         IF (LONG(I)) XTEMP(I) = XTEMP(I) + XPAR(I)                     DEC 2710
         IF (.NOT.LONG(I)) XTEMP(I) = 0.0                               DEC 2720
  170 CONTINUE                                                          DEC 2730
      RETURN                                                            DEC 2740
99999 FORMAT (1H1, 3I5, 1PE12.5)                                        DEC 2750
99998 FORMAT (32H L,I,BATE,XTEM,XTEMP(J1),AKDJQ =, 2I5,1P,4E12.5)       DEC 2760
      END                                                               DEC 2770
      SUBROUTINE TERM(T, A, B, D, LOC, NON0, KD, NONP, LOCP, NQ, LONG,  TER   10
     *     CIMN, CSUM, CIMB, XP, XTEMP, AP, XNEW, MX, IOUT, LOCA,       TER   20
     *     IAPMAX, NUCL)                                                TER   30
      LOGICAL LONG                                                      TER   40
      INTEGER*2 LOC,NON0,KD,NONP,LOCP,NQ,LOCA                           TER   50
      DOUBLE PRECISION BATE, BATM, CIMN, CSUM                           TER   60
      DIMENSION LOCA(1)                                                 TER   70
      DIMENSION A(1), B(1), D(1), LOC(1), NON0(1), KD(1), NONP(1),      TER   80
     *     LOCP(1), NQ(1), LONG(1), CIMN(1), CSUM(1), CIMB(1), XP(1),   TER   90
     *     XTEMP(1), AP(1), XNEW(MX,1), NUCL(1)                         TER  100
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  TER  110
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           TER  120
     *     TFPBAV(4), IPHMAX                                            TER  130
      COMMON /SZ/ NSIZE(25)                                             TER  140
      COMMON /TERMD/ QUB(50), DD(100), DXP(100), QUEUE(50), NQUEUE(50), TER  150
     *     NQU(50)                                                      TER  160
C      LEVEL 2,B,D,LOC,NON0,KD,NONP,LOCP,NQ,CIMN,CSUM,CIMB,XP,XTEMP,AP, TER  170
C    1  XNEW,A,LOCA,NUCL                                                TER  180
      DATA NPRNT /0/                                                    TER  190
C                                                                       TER  200
C     TERM ADDS ONE TERM TO EACH ELEMENT OF THE SOLUTION VECTOR         TER  210
C     CSUM(J) IS THE CURRENT APPROXIMATION TO XNEW(M,J)                 TER  220
C     CIMB(J) IS THE VECTOR CONTAINING THE LAST TERM ADDED TO EACH      TER  230
C     ELEMENT OF CSUM(J); SAME AS CIM0 IN MATREX                        TER  240
C     CIMN(J) IS THE VECTOR CONTAINING 1/TON TIMES THE NEW TERM TO BE   TER  250
C     ADDED TO  CSUM(J)                                                 TER  260
C     CIMN(J) IS GENERATED FROM CIM0(J) BY A RECURSION RELATION:        TER  270
C          CIMN(J)= SUM OVER L OF (AP(J,L)*CIM0(L))                     TER  280
C     AP(I,J) IS THE REDUCED TRANSITION MATRIX FOR THE LONG-LIVED       TER  290
C     NUCLIDES                                                          TER  300
C                                                                       TER  310
C    MUCH OF TERM IS THE SAME AS DECAY AND, AS A RESULT, DETAILED       TER  320
C    COMMENT CARDS WILL NOT BE REPEATED.  ONLY THOSE UNIQUE ASPECTS     TER  330
C    OF TERM WILL BE HIGHLIGHTED.                                       TER  340
C                                                                       TER  350
      NUL = 0                                                           TER  360
      NN = 0                                                            TER  370
C     FIRST CONSTRUCT REDUCED TRANSITION MATRIX FOR LONG-LIVED ISOTOPES TER  380
C    CONSTRUCT REDUCED TRANSITION MATRIX FROM HERE THROUGH LABEL 220    TER  390
      DO 250 L=1,ITOT                                                   TER  400
C    SKIP NUCLIDE IF IT IS SHORT-LIVED OR IT HAS NO REACTIONS           TER  410
         IF (.NOT.LONG(L)) GO TO 240                                    TER  420
         NUM = NON0(L)                                                  TER  430
         IF (INDEX.EQ.2) NUM = KD(L)                                    TER  440
C    CIMB, WHICH CONTAINS THE NONHOMOGENEOUS TERMS, IS SET UP           TER  450
C    THROUGHOUT TERM FOR USE IN SUBROUTINE MATREX                       TER  460
         CIMB(L) = B(L)                                                 TER  470
         IF (NUM.LE.NUL) GO TO 240                                      TER  480
C    NS= # OF FIRST REACTION FOR NUCLIDE L                              TER  490
         NS = NN + 1                                                    TER  500
         N = NUL                                                        TER  510
         NL = NUM - NUL                                                 TER  520
         IF (NL.LT.1) GO TO 240                                         TER  530
C    LOOP OVER FIRST PARENTS OF DAUGHTER NUCLIDE L                      TER  540
         DO 230 N1=1,NL                                                 TER  550
            N = N + 1                                                   TER  560
            J = LOC(N)                                                  TER  570
            DJ = -D(J)                                                  TER  580
C                                                                       TER  590
C      THIS IS A TEST TO SEE IF ONE OF THE ASSYMPTOTIC SOLUTIONS APPLIESTER  600
C                                                                       TER  610
C     THE CARDS THROUGH "GO TO 200" PUT THE RECTION RATES OF THE        TER  620
C     LONG-LIVED NUCLIDES DIRECTLY IN THE REDUCED MATRIX AP             TER  630
C    NUCLIDES DIRECTLY IN AP                                            TER  640
            IF (.NOT.LONG(J)) GO TO 10                                  TER  650
            NN = NN + 1                                                 TER  660
            AP(NN) = A(N)                                               TER  670
            LOCP(NN) = J                                                TER  680
            GO TO 230                                                   TER  690
C                                                                       TER  700
C    GOING BACK UP THE CHAIN TO FIND A PARENT WHICH IS NOT IN           TER  710
C    EQUILIBRIUM                                                        TER  720
C                                                                       TER  730
   10       NSAVE = 0                                                   TER  740
            QUE = A(N)/DJ                                               TER  750
C    DRB= (1-(1/(DJ*T))) IS A MEASURE OF THE CLOSENESS OF THE CURRENT   TER  760
C    FARTHEST PARENT TO PERFECT LONG-LIVED BEHAVIOR.  IT IS A           TER  770
C    CUMULTIVE FUNCTION FOR A SPECIFIC CHAIN.                           TER  780
            DRB = 1.0                                                   TER  790
            CIMB(L) = CIMB(L) + QUE*B(J)                                TER  800
            NQ(L) = 0                                                   TER  810
            NQ(J) = L                                                   TER  820
   20       NUX = NON0(J)                                               TER  830
            IF (INDEX.EQ.2) NUX = KD(J)                                 TER  840
            NUF = 0                                                     TER  850
            IF (J.GT.1) NUF = NON0(J-1)                                 TER  860
            NX = NUX - NUF                                              TER  870
            IF (NX.LT.1) GO TO 220                                      TER  880
            K = NUF                                                     TER  890
            IF (NX.LT.1) GO TO 220                                      TER  900
            DO 210 KK=1,NX                                              TER  910
               K = K + 1                                                TER  920
               J1 = LOC(K)                                              TER  930
               DJ = -D(J1)                                              TER  940
               KP = J                                                   TER  950
               IF (J1.EQ.L) GO TO 40                                    TER  960
   30          IF (J1.EQ.NQ(KP)) GO TO 210                              TER  970
               KP = NQ(KP)                                              TER  980
               IF (KP.NE.0) GO TO 30                                    TER  990
   40          CONTINUE                                                 TER 1000
               AKDJQ = QUE*A(K)/DJ                                      TER 1010
C    SKIP CALCULATION UNTIL A LONG-LIVED J1 IS FOUND                    TER 1020
               IF (.NOT.LONG(J1)) GO TO 200                             TER 1030
C    TRM=FRACTION OF LONG-LIVED NUCLIDE REMAINING AT END OF TIME STEP   TER 1040
               TRM = 1.0 - XP(J1)                                       TER 1050
               IF (TRM.LT.1.0E-6) GO TO 160                             TER 1060
               IF (J1.NE.L) NQ(J1) = J                                  TER 1070
               DD(1) = -D(J1)                                           TER 1080
               DXP(1) = XP(J1)                                          TER 1090
               I = 2                                                    TER 1100
               KP = J                                                   TER 1110
   50          DD(I) = -D(KP)                                           TER 1120
               DXP(I) = XP(KP)                                          TER 1130
               KP = NQ(KP)                                              TER 1140
               IF (KP.EQ.0) GO TO 60                                    TER 1150
               I = I + 1                                                TER 1160
               IF (I.LE.100) GO TO 50                                   TER 1170
C    IF QUEUE OF SHORT-LIVED NUCLIDES EXCEEDS 100 ISOTOPES, TERMINATE   TER 1180
C    CHAIN AND WRITE MESSAGE                                            TER 1190
               WRITE (15,99999) L, J1, J, AKDJQ                         TER 1200
               GO TO 220                                                TER 1210
   60          BATM = 0.0                                               TER 1220
               IM = I - 1                                               TER 1230
               IF (IM.LT.2) GO TO 150                                   TER 1240
               DO 140 I=2,IM                                            TER 1250
                  DL = DD(I)                                            TER 1260
                  XPL = DXP(I)                                          TER 1270
                  BATE = 0.0                                            TER 1280
                  I1 = I - 1                                            TER 1290
C    D R VONDY FORM OF BATEMAN EQUATIONS -- ORNL-TM-361                 TER 1300
                  IF (I1.LT.1) GO TO 120                                TER 1310
                  DO 110 KB=1,I1                                        TER 1320
                     XPJ = DXP(KB)                                      TER 1330
                     IF (XPL+XPJ.LT.ERR) GO TO 110                      TER 1340
                     DK = DD(KB)                                        TER 1350
                     PROD = (DL/DK-1.0)                                 TER 1360
                     DKR = PROD                                         TER 1370
                     IF (ABS(PROD).GT.1.E-4) GO TO 70                   TER 1380
C    USE THIS FORM FOR TWO NEARLY EQUAL HALF-LIVES                      TER 1390
                     PROD = T*DK*XPJ*(1.0-0.5*(DL-DK)*T)                TER 1400
                     GO TO 80                                           TER 1410
   70                PROD = (XPJ-XPL)/PROD                              TER 1420
                     PRO1 = XPJ/DKR                                     TER 1430
   80                PI = 1.0                                           TER 1440
                     S1 = 2./(DK*T)                                     TER 1450
                     DO 100 JK=1,I1                                     TER 1460
                        IF (JK.EQ.KB) GO TO 100                         TER 1470
                        S = 1.0 - DK/DD(JK)                             TER 1480
                        IF (ABS(S).GT.1.E-4) GO TO 90                   TER 1490
                        IF (ABS(DKR).GT.1.0E-4) PROD = PRO1             TER 1500
                        S = S1                                          TER 1510
   90                   PI = PI*S                                       TER 1520
                        IF (ABS(PI).GT.1.E25) GO TO 110                 TER 1530
  100                CONTINUE                                           TER 1540
                     BATE = BATE + PROD/PI                              TER 1550
  110             CONTINUE                                              TER 1560
  120             CONTINUE                                              TER 1570
C    IF SUMMATION IS NEGATIVE, SET EQUAL TO ZERO AND PRINT MESSAGE      TER 1580
                  IF (BATE.GE.0.) GO TO 130                             TER 1590
                  NPRNT = NPRNT + 1                                     TER 1600
                  IF (NPRNT.LE.10) WRITE (15,99998) L, IM, BATE, BATM   TER 1610
                  BATE = 0.                                             TER 1620
  130             CONTINUE                                              TER 1630
C    ACCUMULATE THE EFFECT OF ALL SHORT-LIVED PRECURSORS ON THE         TER 1640
C    LONG-LIVED DAUGHTER IN PARAMETER BATM                              TER 1650
                  BATM = BATM + BATE                                    TER 1660
  140          CONTINUE                                                 TER 1670
  150          CONTINUE                                                 TER 1680
C    CALCULATE DRA, THE EFFECTIVE PRODUCTION RATE OF LONG-LIVED         TER 1690
C    NUCLIDE L FROM LONG-LIVED NUCLIDE J1 VIA THE CURRENT CHAIN.        TER 1700
C    USE FIRST EQUATION FOR A 2-MEMBER CHAIN.  APPROACH IS TO           TER 1710
C    MODIFY PRODUCTION RATE AKDJQ*DJ BY TERMS ACCOUNTING FOR            TER 1720
C    THE EFFECTS OF SHORT-LIVED INTERMEDIATES.                          TER 1730
               DRA = AKDJQ*DJ*(TRM-BATM)/TRM                            TER 1740
               GO TO 170                                                TER 1750
  160          DRA = AKDJQ*AMAX1(DRB,0.0)*DJ                            TER 1760
  170          IF (NS.GT.NN) GO TO 190                                  TER 1770
C    ACCUMULATE PRODUCTION OF L FROM J1 IN AP FOR ALL CHAINS PRODUCING LTER 1780
               DO 180 LJ=NS,NN                                          TER 1790
                  IF (LOCP(LJ).NE.J1) GO TO 180                         TER 1800
C    ACCUMULATE IF AP AND LOCP HAVE ALREADY BEEN INITIALIZED            TER 1810
                  AP(LJ) = AP(LJ) + DRA                                 TER 1820
                  GO TO 210                                             TER 1830
  180          CONTINUE                                                 TER 1840
  190          NN = NN + 1                                              TER 1850
C    INITIALIZE AP WITH DRA AND LOCP WITH J1, THE LONG-LIVED PARENT     TER 1860
               AP(NN) = DRA                                             TER 1870
               LOCP(NN) = J1                                            TER 1880
               GO TO 210                                                TER 1890
C    ONLY SAVE-FOR-FURTHER-PURSUIT THOSE NUCLIDES THAT ARE SHORT-LIVED  TER 1900
C    AND WHERE ITS LONG-LIVED PARENT HAS NOT YET BEEN REACHED.          TER 1910
  200          IF (AKDJQ.LE.1.0E-06) GO TO 210                          TER 1920
               IF (NSAVE.GE.50) GO TO 210                               TER 1930
               NSAVE = NSAVE + 1                                        TER 1940
               NQUEUE(NSAVE) = J1                                       TER 1950
               QUEUE(NSAVE) = AKDJQ                                     TER 1960
               NQU(NSAVE) = J                                           TER 1970
C    QUB(NSAVE)= ACCUMULATED DRB                                        TER 1980
               QUB(NSAVE) = DRB - 1./(DJ*T)                             TER 1990
  210       CONTINUE                                                    TER 2000
  220       IF (NSAVE.LE.0) GO TO 230                                   TER 2010
            J = NQUEUE(NSAVE)                                           TER 2020
            QUE = QUEUE(NSAVE)                                          TER 2030
            NQ(J) = NQU(NSAVE)                                          TER 2040
            DRB = QUB(NSAVE)                                            TER 2050
            CIMB(L) = CIMB(L) + QUE*B(J)*AMAX1(DRB,0.0)                 TER 2060
            NSAVE = NSAVE - 1                                           TER 2070
            GO TO 20                                                    TER 2080
  230    CONTINUE                                                       TER 2090
  240    NUL = NON0(L)                                                  TER 2100
C    SET NONP(L) AS THE NUMBER (NN) OF THE LAST REACTION IN AP FOR      TER 2110
C    NUCLIDE                                                            TER 2120
C      SET NONP(L) AS THE NUMBER (NN) OF THE LASTE RECTION              TER 2130
C      IN AP FOR NUCLIDE L                                              TER 2140
         NONP(L) = NN                                                   TER 2150
  250 CONTINUE                                                          TER 2160
C    CHECK TO SEE IF AP OVERFLOWS.                                      TER 2170
      IF (NN.GT.IAPMAX) WRITE (15,99997) NN, IAPMAX                     TER 2180
C     FIND NORM OF MATRIX AND ESTIMATE ERROR AS DESCRIBED IN LAPIDUS    TER 2190
C     AND LUUS, OPTIMAL CONTROL OF ENGINEERING PROCESSES BLAISDELL 1967 TER 2200
C     FIND THE MINIMUM OF THE MAXIMUM ROW SUM AND THE MAXIMUM COLUMN SUMTER 2210
C                                                                       TER 2220
C    CALCULATE THE NUMBER OF TERMS (=NLARGE) NEEDED IN                  TER 2230
C    THE MATRIX EXPONENTIAL                                             TER 2240
C                                                                       TER 2250
      ASUM = 0.0                                                        TER 2260
      ASUMJ = 0.0                                                       TER 2270
      NUL = 1                                                           TER 2280
      DO 310 I=1,ITOT                                                   TER 2290
         IF (.NOT.LONG(I)) GO TO 300                                    TER 2300
         DI = -D(I)*T                                                   TER 2310
         AJ = DI                                                        TER 2320
         NUM = NONP(I)                                                  TER 2330
         IF (NUL.GT.NUM) GO TO 270                                      TER 2340
         DO 260 N=NUL,NUM                                               TER 2350
C    AJ= ROW SUM                                                        TER 2360
C    AI= COLUMN SUM                                                     TER 2370
C    ASUM= MAX AI OVER ALL COLUMNS                                      TER 2380
C    ASUMJ= MAX AJ OVER ALL ROWS                                        TER 2390
            AJ = AJ + AP(N)*T                                           TER 2400
  260    CONTINUE                                                       TER 2410
  270    AI = DI + DI                                                   TER 2420
         IF (AI.LE.ASUM) GO TO 280                                      TER 2430
         ASUM = AI                                                      TER 2440
         II = I                                                         TER 2450
  280    IF (AJ.LE.ASUMJ) GO TO 290                                     TER 2460
         ASUMJ = AJ                                                     TER 2470
         IJ = I                                                         TER 2480
  290    CONTINUE                                                       TER 2490
  300    NUL = NONP(I) + 1                                              TER 2500
  310 CONTINUE                                                          TER 2510
      BSUM = ASUM                                                       TER 2520
C    ASUM= MIN(ASUM,ASUMJ)                                              TER 2530
      IF (ASUMJ.LT.ASUM) ASUM = ASUMJ                                   TER 2540
      IK = II                                                           TER 2550
      IF (ASUMJ.LT.ASUM) IK = IJ                                        TER 2560
C    USE ASUM TO DECIDE HOW MANY TERMS ARE REQUIRED AND ESTIMATE ERROR  TER 2570
      NLARGE = 3.5*ASUM + 6.0                                           TER 2580
      XLARGE = NLARGE                                                   TER 2590
      IF (NN.GT.NSIZE(21)) NSIZE(21) = NN                               TER 2600
      ERR1 = EXP(ASUM)*(ASUM*2.71828/XLARGE)**NLARGE/SQRT(6.2832*XLARGE)TER 2610
      WRITE (15,99996) BSUM, ASUMJ, XLARGE, ERR1, NLARGE, NUCL(II),     TER 2620
     *     NUCL(IJ)                                                     TER 2630
      IF (ERR1.GT.1.0E-03) WRITE (15,99995) ERR1, ASUM, NUCL(IK), NLARGETER 2640
C       NEXT GENERATE MATRIX EXPONENTIAL SOLUTION                       TER 2650
C    SET ARRAYS CSUM,CIMN EQUAL TO XTEMP                                TER 2660
      DO 320 I=1,ITOT                                                   TER 2670
         CSUM(I) = XTEMP(I)                                             TER 2680
         CIMN(I) = XTEMP(I)                                             TER 2690
  320 CONTINUE                                                          TER 2700
      CALL MATREX(NLARGE, ITOT, T, ERR, CIMB, CIMN, LONG, NONP, LOCP,   TER 2710
     *     AP, D, CSUM, LOCA)                                           TER 2720
C    STORE CALCULATED AMOUNTS FOR LONG-LIVED NUCLIDES                   TER 2730
      DO 330 I=1,ITOT                                                   TER 2740
         IF (LONG(I)) XNEW(13,I) = CSUM(I)                              TER 2750
  330 CONTINUE                                                          TER 2760
      RETURN                                                            TER 2770
99999 FORMAT (41H1TOO LONG A QUEUE HAS BEEN FORMED IN TERM, 3I5, E12.5) TER 2780
99998 FORMAT (54H1BATE IS NEGATIVE IN TERM. THERE ARE MORE THAN TWO SHO,TER 2790
     *     57HRT-LIVED NUCLIDES IN A CHAIN WITH NEARLY EQUAL DIAGONAL E,TER 2800
     *     7HLEMENTS, /, 18H L,IM,BATE,BATM = , 2I5,1P,2E12.5)          TER 2810
99997 FORMAT (53H **WARNING 9614** NUMBER OF REDUCED MATRIX ELEMENTS =, TER 2820
     *     I5, 28H IS GREATER THAN ALLOWABLE =, I5)                     TER 2830
99996 FORMAT (1H ,1P,4E12.5, 3I7)                                       TER 2840
99995 FORMAT (26H0MAXIMUM ERROR GT 0.001, =, 1PE12.3, 7H,ASUM= ,        TER 2850
     *     0PF12.3, 10H,NUCLIDE= , I7, 9H,NLARGE= , I6)                 TER 2860
      END                                                               TER 2870
      SUBROUTINE MATREX(NLARGE, ITOT, T, ERR, CIM0, CIMN, LONG, NONP,   MAT   10
     *     LOCP, AP, D, CSUM, LOCA)                                     MAT   20
      LOGICAL LONG                                                      MAT   30
      INTEGER*2 LOCP,NONP,LOCA                                          MAT   40
      DOUBLE PRECISION CIMN, CSUM                                       MAT   50
C    CIMNI SHOULD REMAIN DOUBLE PRECISION FOR ALL COMPUTERS             MAT   60
      DOUBLE PRECISION CIMNI                                            MAT   70
      DIMENSION CIM0(1), CIMN(1), LONG(1), NONP(1), LOCA(1), LOCP(1),   MAT   80
     *     AP(1), D(1), CSUM(1)                                         MAT   90
C      LEVEL 2,CIM0,CIMN,NONP,LOCP,AP,D,CSUM,LOCA                       MAT  100
      NUL = 1                                                           MAT  110
      DO 30 IA=1,ITOT                                                   MAT  120
         CIMB = CIM0(IA)                                                MAT  130
         CIM0(IA) = CIMN(IA)                                            MAT  140
         CIMN(IA) = CIMB                                                MAT  150
C    NOTE THAT CIM0 (CIMB IN TERM) IS DESTROYED IN THIS OPERATION       MAT  160
         NUM = NONP(IA)                                                 MAT  170
C    WE WILL SKIP THE TEST ON LONG                                      MAT  180
C    WE WILL ASSUME AS PER EXAMINATION OF                               MAT  190
C    TERM THAT THIS TEST IS UNNECESSARY.                                MAT  200
         IF (NUL.GT.NUM) GO TO 20                                       MAT  210
         DO 10 NA=NUL,NUM                                               MAT  220
            LOCA(NA) = IA                                               MAT  230
   10    CONTINUE                                                       MAT  240
   20    NUL = NONP(IA) + 1                                             MAT  250
   30 CONTINUE                                                          MAT  260
      NONT = NONP(ITOT)                                                 MAT  270
C    THIS IS PART OF THE ASSUMPTION THAT THE VALUES                     MAT  280
C    OF NUL AND NUM ABOVE ARE EQUAL WHEN LONG(IA) IS FALSE              MAT  290
C    THIS COMPLETES THE SETTING UP OF LOCA                              MAT  300
      DO 60 NT=1,NLARGE                                                 MAT  310
         TON = T/FLOAT(NT)                                              MAT  320
C    THE FOLLOWING DO-LOOP HAS BEEN GREATLY SHORTENED                   MAT  330
         DO 40 N=1,NONT                                                 MAT  340
            I = LOCA(N)                                                 MAT  350
            J = LOCP(N)                                                 MAT  360
            CIMN(I) = CIMN(I) + AP(N)*CIM0(J)                           MAT  370
   40    CONTINUE                                                       MAT  380
         DO 50 IB=1,ITOT                                                MAT  390
            IF (.NOT.LONG(IB)) GO TO 50                                 MAT  400
            CIMNI = (CIMN(IB)+D(IB)*CIM0(IB))*TON                       MAT  410
            CIMN(IB) = 0.0                                              MAT  420
            IF (DABS(CIMNI).LE.ERR) CIMNI = 0.0                         MAT  430
            CIM0(IB) = CIMNI                                            MAT  440
            CSUM(IB) = CIMNI + CSUM(IB)                                 MAT  450
   50    CONTINUE                                                       MAT  460
   60 CONTINUE                                                          MAT  470
      RETURN                                                            MAT  480
      END                                                               MAT  490
      SUBROUTINE EQUIL(A, LOC, B, D, NON0, KD, LONG, XP, XPAR, XTEMP,   EQU   10
     *     XNEW, MX, IOUT, NUCL)                                        EQU   20
C    EQUIL PUTS SHORT-LIVED PROGENY IN EQUILIBRIUM WITH PARENTS;        EQU   30
C    EQUIL USES GAUSS-SEIDEL ITERATION TO GENERATE STEADY-STATE         EQU   40
C    CONCENTRATIONS                                                     EQU   50
      LOGICAL LONG                                                      EQU   60
      INTEGER*2 LOC,NON0,KD                                             EQU   70
      DIMENSION A(1), LOC(1), B(1), D(1), NON0(1), KD(1), LONG(1),      EQU   80
     *     XP(1), XPAR(1), XTEMP(1), XNEW(MX,1), NUCL(1)                EQU   90
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  EQU  100
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           EQU  110
     *     TFPBAV(4), IPHMAX                                            EQU  120
C      LEVEL 2,LOC,B,D,NON0,KD,XP,XPAR,XTEMP,XNEW,A,NUCL                EQU  130
C                                                                       EQU  140
      DO 10 I=1,ITOT                                                    EQU  150
         XPAR(I) = 0.0                                                  EQU  160
         IF (.NOT.LONG(I)) GO TO 10                                     EQU  170
         XTEMP(I) = XTEMP(I)*XP(I)                                      EQU  180
         XPAR(I) = AMAX1(XNEW(13,I)-XTEMP(I),0.0)                       EQU  190
   10 CONTINUE                                                          EQU  200
      ITER = 1                                                          EQU  210
   20 N = 0                                                             EQU  220
      BIG = 0.0                                                         EQU  230
      DO 70 I=1,ITOT                                                    EQU  240
         NUM = NON0(I) - N                                              EQU  250
         DI = -D(I)                                                     EQU  260
         IF (LONG(I)) GO TO 60                                          EQU  270
         XNW = B(I)                                                     EQU  280
         IF (INDEX.EQ.2) NUM = KD(I) - N                                EQU  290
         IF (NUM.EQ.0) GO TO 40                                         EQU  300
         DO 30 K=1,NUM                                                  EQU  310
            N = N + 1                                                   EQU  320
            J = LOC(N)                                                  EQU  330
            DJ = -D(J)                                                  EQU  340
            XJ = XPAR(J)                                                EQU  350
            IF (LONG(J)) XJ = XJ + XTEMP(J)/(1.0-DJ/DI)                 EQU  360
            XNW = XNW + A(N)*XJ                                         EQU  370
   30    CONTINUE                                                       EQU  380
   40    XNW = XNW/DI                                                   EQU  390
         IF (XNW.LT.1.0E-35) GO TO 50                                   EQU  400
         ARG = ABS((XNW-XPAR(I))/XNW)                                   EQU  410
         IF (ARG.GT.BIG) BIG = ARG                                      EQU  420
         IF (ARG.GT.BIG) ISV = I                                        EQU  430
         IF (ARG.GT.QXN .AND. ITER.EQ.99) WRITE (IOUT,99998) BIG, ARG,  EQU  440
     *        QXN, I, ITER                                              EQU  450
   50    XPAR(I) = XNW                                                  EQU  460
   60    N = NON0(I)                                                    EQU  470
   70 CONTINUE                                                          EQU  480
      IF (BIG.LT.QXN) GO TO 80                                          EQU  490
      ITER = ITER + 1                                                   EQU  500
      IF (ITER.LT.100) GO TO 20                                         EQU  510
      WRITE (15,99998) BIG, ARG, QXN, ISV, ITER                         EQU  520
   80 DO 90 I=1,ITOT                                                    EQU  530
         IF (.NOT.LONG(I)) XNEW(13,I) = XNEW(13,I) + XPAR(I)            EQU  540
C    ZERO OUT NGEATIVE CONCENTRATIONS                                   EQU  550
C    PRINT MESSAGE IF THEY ARE BECOMING TOO "LARGE"                     EQU  560
         IF (-XNEW(13,I).GT.1.0E-15) WRITE (15,99999) I, XNEW(13,I)     EQU  570
         IF (XNEW(13,I).LT.0.0) XNEW(13,I) = 0.0                        EQU  580
   90 CONTINUE                                                          EQU  590
      RETURN                                                            EQU  600
99999 FORMAT (29H XNEW.LT.-1.0E-15 IN EQUIL;I=, I6, 2X,8HXNEW(I)=,      EQU  610
     *     1PE12.4)                                                     EQU  620
99998 FORMAT (49H GAUSS-SEIDEL ITERATION DID NOT CONVERGE IN EQUIL,     EQU  630
     *    1P,3E10.2, 2I8)                                               EQU  640
      END                                                               EQU  650
      SUBROUTINE NOAH(NUCLI)                                            NOA   10
C    SUBROUTINE NOAH CONVERTS SIX DIGIT IDENTIFIER TO ALPHAMERIC SYMBOL NOA   20
      INTEGER ELE, STA                                                  NOA   30
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 NOA   40
C      LEVEL 2,NUCLI                                                    NOA   50
      IS = MOD(NUCLI,10) + 1                                            NOA   60
      NZ = NUCLI/10000                                                  NOA   70
      MW = NUCLI/10 - NZ*1000                                           NOA   80
      NAME(1) = ELE(NZ)                                                 NOA   90
      NAME(2) = MW                                                      NOA  100
      NAME(3) = STA(IS)                                                 NOA  110
      RETURN                                                            NOA  120
      END                                                               NOA  130
      FUNCTION RMASS(NA, I)                                             RMA   10
C    THIS FUNCTION CONVERTS MASS UNITS                                  RMA   20
      RMASS = 1.                                                        RMA   30
      IF (I.EQ.1) RMASS = (MOD(NA,10000)/10)                            RMA   40
      RETURN                                                            RMA   50
      END                                                               RMA   60
      FUNCTION RTIME(I)                                                 RTI   10
C    THIS SUBROUTINE SUPPLIES TIME CONVERSION FACTORS                   RTI   20
      GO TO (10, 20, 30, 40, 50, 10, 60, 70, 80), I                     RTI   30
   10 RTIME = 1.                                                        RTI   40
      RETURN                                                            RTI   50
   20 RTIME = 60.                                                       RTI   60
      RETURN                                                            RTI   70
   30 RTIME = 3600.                                                     RTI   80
      RETURN                                                            RTI   90
   40 RTIME = 86400.                                                    RTI  100
      RETURN                                                            RTI  110
   50 RTIME = 3.15576E+07                                               RTI  120
      RETURN                                                            RTI  130
   60 RTIME = 3.15576E+10                                               RTI  140
      RETURN                                                            RTI  150
   70 RTIME = 3.15576E+13                                               RTI  160
      RETURN                                                            RTI  170
   80 RTIME = 3.15576E+16                                               RTI  180
      RETURN                                                            RTI  190
      END                                                               RTI  200
      FUNCTION IPAGE(J)                                                 IPA   10
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    IPA   20
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  IPA   30
     *     ITITLE(18)                                                   IPA   40
      IF (J.NE.2) GO TO 10                                              IPA   50
      IF (NTOC.EQ.NTOCA) IPAGE = IPA                                    IPA   60
      IF (NTOC.EQ.NTOCB) IPAGE = IPB                                    IPA   70
      RETURN                                                            IPA   80
   10 CONTINUE                                                          IPA   90
      IF (J.NE.0) GO TO 20                                              IPA  100
      IPA = 1                                                           IPA  110
      IPB = 1                                                           IPA  120
      NTOC = NTOCA                                                      IPA  130
      IPAGE = IPA                                                       IPA  140
      GO TO 40                                                          IPA  150
   20 IF (NTOC.NE.NTOCA) GO TO 30                                       IPA  160
      IPA = IPA + 1                                                     IPA  170
      IPAGE = IPA                                                       IPA  180
      GO TO 40                                                          IPA  190
   30 IPB = IPB + 1                                                     IPA  200
      IPAGE = IPB                                                       IPA  210
   40 RETURN                                                            IPA  220
      END                                                               IPA  230
      SUBROUTINE LISTIT(LOOK, LIST, LWRITE)                             LIS   10
      DIMENSION CARD(21)                                                LIS   20
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    LIS   30
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  LIS   40
     *     ITITLE(18)                                                   LIS   50
      DATA LT /4H$   /                                                  LIS   60
      I = 0                                                             LIS   70
      CALL TOC(NPAGE, 5, LOOK, LIST, LWRITE, ITITLE, IOUT, NTOC, NREC)  LIS   80
      WRITE (LIST,99998) NPAGE                                          LIS   90
      CALL IDENTIFY (IOUT, 0)                                           PC
      NPAGE = IPAGE(1)                                                  LIS  100
      WRITE (LIST,99997) LOOK                                           LIS  110
C   FOR THE CDC VERSION, DELETE THE FOLLOWING CARD AND INSERT THE TWO   LIS  111
C   SUCCEEDING COMMENT CARDS                                            LIS  112
   40 READ(LOOK,1000,END=30) (CARD(J),J=1,21)                           LIS  120
C  40 READ(LOOK,1000) (CARD(J),J=1,21)                                  LIS  121
C     IF(EOF(LOOK).NE.0)GO TO 30                                        LIS  122
 1000 FORMAT(A1,A3,19A4)                                                LIS  130
      I = I + 1                                                         LIS  140
      IF (MOD(I,50).NE.0) GO TO 10                                      LIS  150
      WRITE (LIST,99998) NPAGE                                          LIS  160
      CALL IDENTIFY (IOUT, 0)                                           PC
      WRITE (LIST,99997) LOOK                                           LIS  170
      NPAGE = IPAGE(1)                                                  LIS  180
   10 CONTINUE                                                          LIS  190
      IF (LT.EQ.CARD(1)) GO TO 20                                       LIS  200
      WRITE (LWRITE, 1000) CARD                                         LIS  210
   20 WRITE (LIST,99999) LOOK, LWRITE, I, CARD                          LIS  220
      GO TO    40                                                       LIS  230
   30 ENDFILE LWRITE                                                    LIS  240
      RETURN                                                            LIS  250
99999 FORMAT (1H , 2X, I2, T12, I2, T22, I4, T30, A1, A3, 19A4)         LIS  260
99998 FORMAT (1H1, T120, 4HPAGE, I5)                                    LIS  270
99997 FORMAT (32H LISTING OF INPUT DATA ON UNIT =, I2, /, /, 6H INPUT,  LIS  280
     *     T10, 5HWRITE, T20, 4HCARD, T60, 10HCARD IMAGE, /, 5H UNIT,   LIS  290
     *     T10, 4HUNIT, T20, 6HNUMBER)                                  LIS  300
      END                                                               LIS  310
      SUBROUTINE ADDMOV(IOPN, IFROM, ITO, ITOT, RM, XSTORE, LX, XTEMP,  ADD   10
     *     XNEW, MX, STTFPB, ISTOTI, IS, RSTOTI, XP)                    ADD   20
      DIMENSION STTFPB(LX,10), ISTOTI(LX,3), IS(LX), RSTOTI(LX)         ADD   30
      DIMENSION MYN(10), XP(1)                                          ADD   40
      DIMENSION XSTORE(LX,1), XTEMP(1), XNEW(MX,1)                      ADD   50
      COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)   ADD   60
C      LEVEL 2,XSTORE,XTEMP,XNEW                                        ADD   70
      DATA MYN /0,1,1,1,0,1,1,1,1,1/                                    ADD   80
C       THIS SUBROUTINE MANIPULATES THE OUTPUT,STORAGE, AND TEMPORARY   ADD   90
C       VECTORS.  MOV COPIES ONE VECTOR OVER ANOTHER.  ADD ADDS OR      ADD  100
C       SUBTRACTS TWO VECTORS.  TITLES, FLUXES, SPECIFIC POWERS,        ADD  110
C       BURNUPS, AND K-INFINITY DATA ARE ALSO MANIPULATED.              ADD  120
C    IF A=1.0,ADD IS BEING EXECUTED; IF A=0.0,MOV IS BEING EXECUTED     ADD  130
      (个子程序主要用于处理输出、存储和临时向量。它包括以下功能：

MOV：将一个向量复制到另一个向量。
ADD：对两个向量进行加法或减法运算。
处理标题、通量、特定功率、燃耗以及 K 无穷大数据等的操作。
根据参数 A 的取值，确定执行 ADD 还是 MOV 操作。当 A=1.0 时执行 ADD 操作，当 A=0.0 时执行 MOV 操作。)
      L9=0                                                                      
      IF(RM.EQ.XP(1).AND.XP(1).EQ.XP(ITOT+1))L9=1                               
      MFROM = -IFROM                                                    ADD  140
      MTO = -ITO                                                        ADD  150
      IF (IFROM.GT.0 .AND. IFROM.GT.MX) WRITE (15,99999) IFROM, MX      ADD  160
      IF (MFROM.GT.0 .AND. MFROM.GT.LX) WRITE (15,99998) MFROM, LX      ADD  170
      IF (ITO.GT.0 .AND. ITO.GT.MX) WRITE (15,99999) ITO, MX            ADD  180
      IF (MTO.GT.0 .AND. MTO.GT.LX) WRITE (15,99998) MTO, LX            ADD  190
      A = 0.0                                                           ADD  200
      IOU = 10                                                          ADD  210
      IF (IOPN.LT.2) A = 1.0                                            ADD  220
      IF (IFROM.LT.0) GO TO 120                                         ADD  230
      IF (ITO.LE.0) GO TO 50                                            ADD  240
C    OUTPUT VECTOR TO OUTPUT VECTOR                                     ADD  250
      DO 10 I=1,ITOT                                                    ADD  260
      IF(L9.EQ.1)RM=XP(I)                                                       
         XNEW(ITO,I) = XNEW(ITO,I)*A + XNEW(IFROM,I)*RM                 ADD  270
   10 CONTINUE                                                          ADD  280
      DO 20 J=1,IOU                                                     ADD  290
         OUTFPB(ITO,J) = 0.0                                            ADD  300
         IF (MYN(J).GT.0) OUTFPB(ITO,J) = OUTFPB(IFROM,J)*RM +          ADD  310
     *        OUTFPB(ITO,J)                                             ADD  320
         IF (MYN(J).LT.1) OUTFPB(ITO,J) = OUTFPB(IFROM,J) +             ADD  330
     *        OUTFPB(ITO,J)                                             ADD  340
   20 CONTINUE                                                          ADD  350
      IF (IO(IFROM).LT.0) GO TO 30                                      ADD  360
      ROUTTI(ITO) = ROUTTI(IFROM)                                       ADD  370
      IO(ITO) = IO(IFROM)                                               ADD  380
      RETURN                                                            ADD  390
   30 DO 40 J=1,3                                                       ADD  400
         IOUTTI(ITO,J) = IOUTTI(IFROM,J)                                ADD  410
   40 CONTINUE                                                          ADD  420
      IO(ITO) = -1                                                      ADD  430
      RETURN                                                            ADD  440
C    FROM OUTPUT VECTOR TO TEMPORARY VECTOR                             ADD  450
   50 IF (ITO.LT.0) GO TO 70                                            ADD  460
      DO 60 I=1,ITOT                                                    ADD  470
      IF(L9.EQ.1)RM=XP(I)                                                       
         XTEMP(I) = XTEMP(I)*A + XNEW(IFROM,I)*RM                       ADD  480
   60 CONTINUE                                                          ADD  490
      RETURN                                                            ADD  500
C    FROM OUTPUT VECTOR TO STORAGE VECTOR                               ADD  510
   70 DO 80 I=1,ITOT                                                    ADD  520
      IF(L9.EQ.1)RM=XP(I)                                                       
         XSTORE(-ITO,I) = XSTORE(-ITO,I)*A + XNEW(IFROM,I)*RM           ADD  530
   80 CONTINUE                                                          ADD  540
      DO 90 J=1,IOU                                                     ADD  550
         STTFPB(-ITO,J) = 0.0                                           ADD  560
         IF (MYN(J).LT.1) STTFPB(-ITO,J) = OUTFPB(IFROM,J) +            ADD  570
     *        STTFPB(-ITO,J)                                            ADD  580
         IF (MYN(J).GT.0) STTFPB(-ITO,J) = OUTFPB(IFROM,J)*RM +         ADD  590
     *        STTFPB(-ITO,J)                                            ADD  600
   90 CONTINUE                                                          ADD  610
      IF (IO(IFROM).LT.0) GO TO 100                                     ADD  620
      RSTOTI(-ITO) = ROUTTI(IFROM)                                      ADD  630
      IS(-ITO) = IO(IFROM)                                              ADD  640
      RETURN                                                            ADD  650
  100 DO 110 J=1,3                                                      ADD  660
         ISTOTI(-ITO,J) = IOUTTI(IFROM,J)                               ADD  670
  110 CONTINUE                                                          ADD  680
      IS(-ITO) = -1                                                     ADD  690
      RETURN                                                            ADD  700
C    FROM STORAGE VECTOR TO OUTPUT VECTOR                               ADD  710
  120 IF (ITO.LE.0) GO TO 170                                           ADD  720
      DO 130 I=1,ITOT                                                   ADD  730
      IF(L9.EQ.1)RM=XP(I)                                                       
         XNEW(ITO,I) = XNEW(ITO,I)*A + XSTORE(-IFROM,I)*RM              ADD  740
  130 CONTINUE                                                          ADD  750
      DO 140 J=1,IOU                                                    ADD  760
         OUTFPB(ITO,J) = 0.0                                            ADD  770
         IF (MYN(J).GT.0) OUTFPB(ITO,J) = STTFPB(-IFROM,J)*RM +         ADD  780
     *        OUTFPB(ITO,J)                                             ADD  790
         IF (MYN(J).LT.1) OUTFPB(ITO,J) = STTFPB(-IFROM,J) +            ADD  800
     *        OUTFPB(ITO,J)                                             ADD  810
  140 CONTINUE                                                          ADD  820
      IF (IS(-IFROM).LT.0) GO TO 150                                    ADD  830
      ROUTTI(ITO) = RSTOTI(-IFROM)                                      ADD  840
      IO(ITO) = IS(-IFROM)                                              ADD  850
      RETURN                                                            ADD  860
  150 DO 160 J=1,3                                                      ADD  870
         IOUTTI(ITO,J) = ISTOTI(-IFROM,J)                               ADD  880
  160 CONTINUE                                                          ADD  890
      IO(ITO) = -1                                                      ADD  900
      RETURN                                                            ADD  910
C    FROM STORAGE VECTOR TO TEMPORARY VECTOR                            ADD  920
  170 IF (ITO.LT.0) GO TO 190                                           ADD  930
      DO 180 I=1,ITOT                                                   ADD  940
      IF(L9.EQ.1)RM=XP(I)                                                       
         XTEMP(I) = XTEMP(I)*A + XSTORE(-IFROM,I)*RM                    ADD  950
  180 CONTINUE                                                          ADD  960
      RETURN                                                            ADD  970
C    FROM STORAGE VECTOR TO STORAGE VECTOR                              ADD  980
  190 DO 200 I=1,ITOT                                                   ADD  990
      IF(L9.EQ.1)RM=XP(I)                                                       
         XSTORE(-ITO,I) = XSTORE(-ITO,I)*A + XSTORE(-IFROM,I)*RM        ADD 1000
  200 CONTINUE                                                          ADD 1010
      DO 210 J=1,IOU                                                    ADD 1020
         STTFPB(-ITO,J) = 0.0                                           ADD 1030
         IF (MYN(J).LT.1) STTFPB(-ITO,J) = STTFPB(-IFROM,J) +           ADD 1040
     *        STTFPB(-ITO,J)                                            ADD 1050
         IF (MYN(J).GT.0) STTFPB(-ITO,J) = STTFPB(-IFROM,J)*RM +        ADD 1060
     *        STTFPB(-ITO,J)                                            ADD 1070
  210 CONTINUE                                                          ADD 1080
      IF (IS(-IFROM).LT.0) GO TO 220                                    ADD 1090
      RSTOTI(-ITO) = RSTOTI(-IFROM)                                     ADD 1100
      IS(-ITO) = IS(-IFROM)                                             ADD 1110
      RETURN                                                            ADD 1120
  220 DO 230 J=1,3                                                      ADD 1130
         ISTOTI(-ITO,J) = ISTOTI(-IFROM,J)                              ADD 1140
  230 CONTINUE                                                          ADD 1150
      IS(-ITO) = -1                                                     ADD 1160
      RETURN                                                            ADD 1170
99999 FORMAT (49H **WARNING 9610** NUMBER OF OUTPUT VECTORS USED =, I5, ADD 1180
     *     28H IS GREATER THAN ALLOWABLE =, I5)                         ADD 1190
99998 FORMAT (50H **WARNING 9610** NUMBER OF STORAGE VECTORS USED =,    ADD 1200
     *     I5, 28H IS GREATER THAN ALLOWABLE =, I5)                     ADD 1210
      END                                                               ADD 1220
      BLOCK DATA                                                        BLK   10
      LOGICAL WRDFLG, LSCAN, LALLC, LNEG, LBIN                          BLK   20
      INTEGER ELE, STA                                                  BLK   30
      DOUBLE PRECISION UNIT, TENS, DNUM                                 BLK   40
      DIMENSION FA(100), FB(100), FC(100), FD(100), FE(100), FF(100),   BLK   50
     *     FG(100), FH(100), FI(100), FJ(100)                           BLK   60
      COMMON /BDAT/ ZAP(130), KZAP(12), JZAP(6,26), NTB(3,4),           BLK   70
     *     INSTR0(32), RRCUT(30), NDI(11,7), NDL(6,5), NDO(9,5)         BLK   80
      COMMON /BDCOM/ NUCLA(100), TNPNIF(100), NUCLB(020), ALN(020),     BLK   90
     *     NUCLC(100), SFNDAT(100), CNMUL, CNEXP, FNPNIF(100)           BLK  100
      COMMON /FRECOV/ FRECG(10,20), NGROUP(100), FREC(10,100)           BLK  110
      COMMON /JUNK/ ERR, IDM(1), ILITE, IACT, IFP, ITOT, ILMAX, IAMAX,  BLK  120
     *     IFMAX, ITMAX, IZMAX, AXN, QXN, FLUX, POWER, INDEX,           BLK  130
     *     TFPBAV(4), IPHMAX                                            BLK  140
      COMMON /LABL/ ELE(99), STA(2), LABTIM(9), NAME(3)                 BLK  150
      COMMON /NUPHMN/ MLIB(11), MLPU(100), MLIP(3), MPHO(4),            BLK  160
     *     INCYLD(8,4), MYLD(4)                                         BLK  170
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    BLK  180
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  BLK  190
     *     ITITLE(18)                                                   BLK  200
      COMMON /OUTOPT/ NOUT(20,4), NTO(80),MCOUT                         BLK  210
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      BLK  220
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           BLK  230
      COMMON /VECTIT/ OUTFPB(13,10), IO(13), ROUTTI(13), IOUTTI(13,3)   BLK  240
      EQUIVALENCE (FREC(1,1),FA(1)), (FREC(1,11),FB(1)),                BLK  250
     *     (FREC(1,21),FC(1)), (FREC(1,31),FD(1)), (FREC(1,41),FE(1)),  BLK  260
     *     (FREC(1,51),FF(1)), (FREC(1,61),FG(1)), (FREC(1,71),FH(1)),  BLK  270
     *     (FREC(1,81),FI(1)), (FREC(1,91),FJ(1))                       BLK  280
C        LEVEL 2,FRECG,FA,FB,FC,FD,FE,FF,FG,FH,FI,FJ                    BLK  290
C    DATA FOR FREE FORMAT READING ROUTINES                              BLK  300
      DATA LSCAN /.FALSE./, LALLC /.FALSE./, LNEG /.FALSE./             BLK  310
      DATA CHAR /1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1HA,1HB,1HC,   BLK  320
     *     1HD,1HE,1HF,1H ,1H,,1HR,1H*,1H$,1H&,1H+,1H-,1HZ,1H.,1HO,1HP/ BLK  330
C    SET NREC NEGATIVE TO INITIALIZE IS AND STTFPB IN SUB MAIN3         BLK  340
      DATA NREC /-1/                                                    BLK  350
      DATA UNIT /1.,10.,100.,1.E3,1.E4,1.E5,1.E6,1.E7,1.E8,1.E9/        BLK  360
      DATA TENS /1.,1.E10,1.E20,1.E30,1.E30,1.E30,1.E30,1.E30/          BLK  370
      DATA IRPT /0/, ICHR /80/, LCHR /80/                               BLK  380
C    ZERO OUTFPB INITIALLY                                              BLK  390
      DATA OUTFPB /130*0.0/                                             BLK  400
C    IO DETERMINES THE TYPE OF VECTOR HEADINGS; 1=TIME+UNITS            BLK  410
      DATA IO /13*1/                                                    BLK  420
C    ZERO MPHO INITIALLY                                                BLK  430
      DATA MPHO /4*0/                                                   BLK  440
C                                                                       BLK  450
C    SPECIFICATION OF WHICH ACTINIDES HAVE DIRECT FP YIELDS             BLK  460
      DATA MYLD /4,4,6,8/                                               BLK  470
      DATA INCYLD /3,4,5,6,4*0,1,2,3,5,4*0,1,2,3,4,5,6,2*0,1,2,3,4,5,6, BLK  480
     *     7,8/                                                         BLK  490
C                                                                       BLK  500
C    ZERO NUMBER OF NUCLIDES IN CASE                                    BLK  510
      DATA ILITE, IACT, IFP, ITOT /4*0/                                 BLK  520
C    ERR USED INSTEAD OF ZERO                                           BLK  530
      DATA ERR /1.0E-25/                                                BLK  540
C    NTB CONTROLS THE COMBINATION OF NUCL,ELEM, AND SUMM TABLES PRINTED BLK  550
      DATA NTB /1,1,1,2,2,3,3,4,4,5,6,7/                                BLK  560
C    CUTOFF VALUES FOR SUMMARY TABLES                                   BLK  570
      DATA RRCUT /26*0.001,0.01,1.0E-25,2*0.0/                          BLK  580
      DATA CUTOFF /26*0.001,0.01,1.0E-25,2*0.0/                         BLK  581
C    DEFALT VALUES FOR TABLE PRINT OPTIONS                              BLK  590
      DATA NTO /4*8,1,19*8,4*8,1,19*8,4*8,1,19*8,8*8/                   BLK  600
C    OUTPUT UNITS FOR TABLE OF CONTENTS                                 BLK  610
      DATA NTOCA, NTOCB /12,13/                                         BLK  620
      DATA NPAGE /1/                                                    BLK  630
C    INSTRUCTION INPUT UNIT, 3 OUTPUT UNITS (ALL THE SAME),             BLK  640
C    AND THE NUMBER OF TABLE TYPES                                      BLK  650
      DATA IUNIT, IOUT, JOUT, KOUT, NTAB /50,3*6,24/                    BLK  660
C    TITLES FOR PRINTING DATA FIXED IN MAIN1                            BLK  670
      DATA NDI /4HNEUT,4HRON ,4HYIEL,4HD PE,4HR NE,4HUTRO,4HN-IN,4HDUCE,BLK  680
     *     4HD FI,4HSSIO,4HN   ,4H(ALP,4HHA,N,4H) NE,4HUTRO,4HN YI,     BLK  690
     *     4HELD ,4HPER ,4HFISS,4HION ,4H    ,4H    ,4HSPON,4HTANE,     BLK  700
     *     4HOUS ,4HFISS,4HION ,4HNEUT,4HRON ,4HYIEL,4HD PE,4HR FI,     BLK  710
     *     4HSSN ,4HINDI,4HVIDU,4HAL E,4HLEME,4HNT F,4HRACT,4HIONA,     BLK  720
     *     4HL RE,4HCOVE,4HRIES,4H    ,4HGROU,4HP EL,4HEMEN,4HTAL ,     BLK  730
     *     4HFRAC,4HTION,4HAL R,4HECOV,4HERIE,4HS   ,4H    ,4HELEM,     BLK  740
     *     4HENTA,4HL AS,4HSIGN,4HNMEN,4HT TO,4HFRAC,4H REC,4HOVER,     BLK  750
     *     4HY GR,4HOUP ,4HELEM,4HENTA,4HL CH,4HEMIC,4HAL T,4HOXIC,     BLK  760
     *     4HITIE,4HS   ,4H    ,4H    ,4H    /                          BLK  770
C    DATA USED TO PRINT THE TABLE TITLES IN THE TABLE OF                BLK  780
C    CONTENTS AND ABOVE EACH TABLE PAGE                                 BLK  790
      DATA NDL /4HDECA,4HY DA,4HTA L,4HIBRA,4HRY--,4H----,4HCROS,4HS SE,BLK  800
     *     4HCTIO,4HN LI,4HBRAR,4HY---,4HPHOT,4HON L,4HIBRA,4HRY--,     BLK  810
     *     4H----,4H----,4HDATA,4HLIBR,4HARY ,4HNOME,4HNCLA,4HTURE,     BLK  820
     *     4HNON-,4HZERO,4H MAT,4HRIX ,4HELEM,4HENTS/                   BLK  830
      DATA NDO /4HSPON,4HTANE,4HOUS ,4HFISS,4HION ,4HNEUT,4HRON ,4HSOUR,BLK  840
     *     4HCE  ,4H(ALP,4HHA,N,4H) NE,4HUTRO,4HN SO,4HURCE,4H    ,     BLK  850
     *     4H    ,4H    ,4HACTI,4HVATI,4HON P,4HRODU,4HCT P,4HHOTO,     BLK  860
     *     4HN TA,4HBLE ,4H    ,4HACTI,4HNIDE,4H NUC,4HLIDE,4H PHO,     BLK  870
     *     4HTON ,4HTABL,4HE   ,4H    ,4HFISS,4HION ,4HPROD,4HUCT ,     BLK  880
     *     4HNUCL,4HIDE ,4HPHOT,4HON T,4HABLE/                          BLK  890
      DATA ZAP /4HNUCL,4HIDE ,4HTABL,4HE:  ,4HELEM,4HENT ,4HTABL,4HE:  ,BLK  900
     *     4HSUMM,4HARY ,4HTABL,4HE:  ,4HISOT,4HOPIC,4H COM,4HPOSI,     BLK  910
     *     4HTION,4H OF ,4HEACH,4H ELE,4HMENT,4H,   ,4HCONC,4HENTR,     BLK  920
     *     4HATIO,4HNS, ,4HRADI,4HOACT,4HIVIT,4HY,  ,4HTHER,4HMAL ,     BLK  930
     *     4HPOWE,4HR,  ,4HGAMM,4HA PO,4HWER,,4HRADI,4HOACT,4HIVE ,     BLK  940
     *     4HINHA,4HLATI,4HON H,4HAZAR,4HD,  ,4HRADI,4HOACT,4HIVE ,     BLK  950
     *     4HINGE,4HSTIO,4HN HA,4HZARD,4H,   ,4HCHEM,4HICAL,4H ING,     BLK  960
     *     4HESTI,4HON H,4HAZAR,4HD,  ,4HNEUT,4HRON ,4HABSO,4HRPTI,     BLK  970
     *     4HON R,4HATE,,4HFISS,4HION ,4HRATE,4H,   ,4HATOM,4H FRA,     BLK  980
     *     4HCTIO,4HN   ,4HWEIG,4HHT F,4HRACT,4HION ,4HGRAM,4H-ATO,     BLK  990
     *     4HMS  ,4HGRAM,4HS   ,4HCURI,4HES  ,4HFRAC,4HTION,4HAL  ,     BLK 1000
     *     4HWATT,4HS   ,4HM**3,4H AIR,4H AT ,4HRCG ,4HM**3,4H WAT,     BLK 1010
     *     4HER A,4HT RC,4HG   ,4HM**3,4H WAT,4HER  ,4HNEUT,4HRONS,     BLK 1020
     *     4H PER,4H SEC,4HOND ,4HFISS,4HIONS,4H PER,4H SEC,4HOND ,     BLK 1030
     *     4H    ,4HALPH,4HA RA,4HDIOA,4HCTIV,4HITY ,12*0.0/            BLK 1040
      DATA KZAP /01,02,03,04,05,06,07,08,09,10,11,12/                   BLK 1050
      DATA JZAP /013,022,071,074,113,113,013,022,075,078,113,113,023,   BLK 1060
     *     026,079,081,113,113,023,026,086,088,079,081,023,026,082,083, BLK 1070
     *     113,113,023,026,086,088,082,083,027,030,084,085,113,113,027, BLK 1080
     *     030,086,088,084,085,031,034,089,090,113,113,031,034,086,088, BLK 1090
     *     089,090,035,037,089,090,113,113,035,037,086,088,089,090,038, BLK 1100
     *     045,091,094,113,113,038,045,086,088,091,094,046,053,095,099, BLK 1110
     *     113,113,046,053,086,088,095,099,054,060,100,102,113,113,054, BLK 1120
     *     060,086,088,100,102,061,066,103,107,113,113,061,066,086,088, BLK 1130
     *     103,107,067,070,108,112,113,113,067,070,086,088,108,112,114, BLK 1140
     *     118,084,085,113,113,114,118,086,088,084,085,12*0/            BLK 1150
C    THESE ARE THE ORIGEN COMMANDS                                      BLK 1160
      DATA INSTR0 /4HRDA ,4HTIT ,4HBAS ,4HFAC ,4HOUT ,4HINP ,4HHED ,    BLK 1170
     *     4HREC ,4HCUT ,4HKEQ ,4HDOL ,4HMOV ,4HADD ,4HBUP ,4HPCH ,     BLK 1180
     *     4HLIP ,4HWAC ,4HLIB ,4HPHO ,4HLPU ,4HIRF ,4HIRP ,4HDEC ,     BLK 1190
     *     4HPRO ,4HOPTL,4HOPTA,4HOPTF,4HCON ,4HSTP ,4HEND ,4HGTO ,     BLK 1200
     *     4HTOF  /                                                             
C    ALPHANUMERIC ELEMENT SYMBOLS                                       BLK 1210
      DATA ELE /2H H,2HHE,2HLI,2HBE,2H B,2H C,2H N,2H O,2H F,2HNE,2HNA, BLK 1220
     *     2HMG,2HAL,2HSI,2H P,2H S,2HCL,2HAR,2H K,2HCA,2HSC,2HTI,2H V, BLK 1230
     *     2HCR,2HMN,2HFE,2HCO,2HNI,2HCU,2HZN,2HGA,2HGE,2HAS,2HSE,2HBR, BLK 1240
     *     2HKR,2HRB,2HSR,2H Y,2HZR,2HNB,2HMO,2HTC,2HRU,2HRH,2HPD,2HAG, BLK 1250
     *     2HCD,2HIN,2HSN,2HSB,2HTE,2H I,2HXE,2HCS,2HBA,2HLA,2HCE,2HPR, BLK 1260
     *     2HND,2HPM,2HSM,2HEU,2HGD,2HTB,2HDY,2HHO,2HER,2HTM,2HYB,2HLU, BLK 1270
     *     2HHF,2HTA,2H W,2HRE,2HOS,2HIR,2HPT,2HAU,2HHG,2HTL,2HPB,2HBI, BLK 1280
     *     2HPO,2HAT,2HRN,2HFR,2HRA,2HAC,2HTH,2HPA,2H U,2HNP,2HPU,2HAM, BLK 1290
     *     2HCM,2HBK,2HCF,2HES/                                         BLK 1300
C    GROUND AND EXCITED STATE INDICATOR                                 BLK 1310
      DATA STA /2H  ,2HM /                                              BLK 1320
C    ALPHANUMERIC TIME UNITS                                            BLK 1330
      DATA LABTIM /2HS ,2HMI,2HHR,2HD ,2HYR,2H  ,2HKY,2HMY,2HGY/        BLK 1340
C    NEUTRON YIELDS PER NEUTRON-INDUCED FISSION (THERMAL)               BLK 1350
C    NUCLIDE ID LIST FOR NEUTRON YIELDS PER NEUTRON-INDUCED FISSION     BLK 1360
      DATA NUCLA /832090,882230,882260,882280,892270,902270,902280,     BLK 1370
     *     902290,902300,902320,902330,902340,912310,912320,912330,     BLK 1380
     *     912341,912340,922300,922310,922320,922330,922340,922350,     BLK 1390
     *     922360,922370,922380,922390,922400,932360,932370,932380,     BLK 1400
     *     932390,942360,942370,942380,942390,942400,942410,942420,     BLK 1410
     *     942430,952410,952421,952420,952430,952441,952440,962420,     BLK 1420
     *     962430,962440,962450,962460,962470,962480,962490,972490,     BLK 1430
     *     982490,982500,982510,982520,982530,982540,992530,992541,     BLK 1440
     *     992540,-1,35*0/                                              BLK 1450
      DATA TNPNIF /7*0.0,2.049,0.0,2.418,4*0.0,2.663,3*0.0,0.0,0.000,   BLK 1460
     *     2.499,2.631,2.421,2.734,0.0,2.801,0.0,0.0,0.0,3.005,0.0,0.0, BLK 1470
     *     2.870,0.0,2.833,2.875,3.135,2.934,3.280,0.0,3.277,3.162,     BLK 1480
     *     3.360,3.732,0.0,0.0,3.746,3.434,3.725,3.832,3.858,3.592,     BLK 1490
     *     3.796,0.0,3.760,4.062,3.970,4.140,4.126,4.150,3*0.0,2*0.0,35*BLK 1500
     *     0.0/                                                         BLK 1510
C    NEUTRON YIELDS PER NEUTRON-INDUCED FISSION (FAST)                  BLK 1520
      DATA FNPNIF /7*0.0,2.049,0.0,2.396,4*0.0,2.631,3*0.0,0.0,0.000,   BLK 1530
     *     2.520,2.555,2.468,2.614,0.0,2.776,0.0,0.0,0.0,2.935,0.0,0.0, BLK 1540
     *     2.946,0.0,3.009,2.946,3.024,2.978,3.075,0.0,3.402,3.311,     BLK 1550
     *     3.361,3.653,0.0,0.0,3.868,3.496,3.743,3.898,3.870,3.680,     BLK 1560
     *     3.866,0.0,3.671,4.130,3.813,4.227,4.364,4.151,3*0.0,2*0.0,35*BLK 1570
     *     0.0/                                                         BLK 1580
C    MEASURED (ALPHA,N) NEUTRON YIELDS-- NEUTRONS/(SEC-G ATOM)          BLK 1590
C    FOR OXIDES OF ALPHA EMITTERS                                       BLK 1600
      DATA NUCLB /-1,942380,942390,942400,942420,952410,962420,962440,-1BLK 1610
     *    ,11*0/                                                        BLK 1620
      DATA ALN /4.76E+06,10755.0,40800.0,653.0,9.64E+05,6.46E+09,       BLK 1630
     *     1.40E+08,13*0.0/                                             BLK 1640
C    DEFAULT PARAMETERS FOR CALCULATING (ALPHA,N) DATA IN SUBR ANSF     BLK 1650
      DATA CNMUL, CNEXP /1.1,1.0        /                               BLK 1660
C    MEASURED SPONTANEOUS FISSION NEUTRON                               BLK 1670
C    NUCLIDE IDS FOR MEASURED SF DATA                                   BLK 1680
      DATA NUCLC /922350,922360,922370,922380,922390,932361,932360,     BLK 1690
     *     932370,932380,932390,942360,942370,942380,942390,942400,     BLK 1700
     *     942410,942420,942430,942440,952400,952410,952420,952421,     BLK 1710
     *     952430,952440,952441,962410,962420,962430,962440,962450,     BLK 1720
     *     962460,962480,962500,972490,982490,982500,982520,982540,     BLK 1730
     *     992530,992540,992550,-1,57*0/                                BLK 1740
      DATA SFNDAT /1.695,1.650,1.872,2.000,2.048,1.790,1.783,1.873,     BLK 1750
     *     1.963,2.053,2.220,1.886,2.280,2.240,2.160,2.250,2.150,2.430, BLK 1760
     *     2.300,2.290,2.383,2.475,2.590,2.520,2.657,2.665,2.500,2.590, BLK 1770
     *     2.687,2.760,2.872,3.000,3.320,3.560,3.720,3.440,3.560,3.725, BLK 1780
     *     3.900,3.920,4.040,4.160,0.0,57*0.0/                          BLK 1790
C    INDIVIDUAL ELEMENT FRACTIONAL RECOVERIES                           BLK 1800
      DATA FA /2*0.0,0.0005,0.0,1.0,1.0,0.0,0.0,0.9,1.0,2*0.0,0.0005,   BLK 1810
     *     0.0,1.0,1.0,3*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,BLK 1820
     *     1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*   BLK 1830
     *     0.0,1.0,2*0.0,0.0005,0.0,2*1.0,3*0.0,1.0,2*0.0,0.0005,0.0,2* BLK 1840
     *     1.0,3*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,0.001,  BLK 1850
     *     0.0005,0.0,3*1.0,2*0.0,1.0,2*0.0,0.0005,0.0,2*1.0,3*0.0,1.0/ BLK 1860
      DATA FB /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 1870
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 1880
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 1890
     *     1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,0.001,0.0005,0.0, BLK 1900
     *     3*1.0,2*0.0,1.0,2*0.0,0.0005,0.0,2*1.0,3*0.0,1.0,0.0,1.0,    BLK 1910
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 1920
     *     1.0/                                                         BLK 1930
      DATA FC /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 1940
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 1950
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 1960
     *     1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 1970
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 1980
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 1990
     *     1.0/                                                         BLK 2000
      DATA FD /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2010
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2020
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,0.001,0.0005,0.0,3*1.0,2*0.0, BLK 2030
     *     1.0,0.0,0.0,0.0005,0.0,2*1.0,3*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2040
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2050
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2060
     *     1.0/                                                         BLK 2070
      DATA FE /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2080
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2090
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2100
     *     1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2110
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2120
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2130
     *     1.0/                                                         BLK 2140
      DATA FF /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2150
     *     1.0,2*0.0,1.0,0.0,0.001,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,0.0,  BLK 2160
     *     0.0005,0.0,2*1.0,3*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2170
     *     1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2180
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2190
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2200
     *     1.0/                                                         BLK 2210
      DATA FG /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2220
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2230
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2240
     *     1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2250
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2260
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2270
     *     1.0/                                                         BLK 2280
      DATA FH /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2290
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2300
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2310
     *     1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2320
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2330
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2340
     *     1.0/                                                         BLK 2350
      DATA FI /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2360
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2370
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2380
     *     1.0,0.0,0.0,0.0005,0.0,2*1.0,3*0.0,1.0,0.0,1.0,0.0005,0.0,3* BLK 2390
     *     1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,    BLK 2400
     *     0.0005,0.0,3*1.0,2*0.0,1.0,0.0,1.0,0.0005,0.0,3*1.0,2*0.0,   BLK 2410
     *     1.0/                                                         BLK 2420
      DATA FJ /0.0,1.0,0.0005,0.0,3*1.0,2*0.0,1.0,0.995,1.0,0.0005,     BLK 2430
     *     0.999,1.0,0.2,0.6,1.0,0.0,1.0,0.0,1.0,0.0005,0.0,0.05,0.05,3*BLK 2440
     *     0.0,1.0,0.995,1.0,0.0005,0.9999,1.0,0.02,3*0.0,1.0,0.0,1.0,  BLK 2450
     *     0.0005,0.0,0.001,0.001,3*0.0,1.0,0.0,1.0,0.0005,0.0,0.001,   BLK 2460
     *     0.001,3*0.0,1.0,0.0,1.0,0.0005,0.0,0.001,0.001,3*0.0,1.0,0.0,BLK 2470
     *     1.0,0.0005,0.0,0.001,0.001,3*0.0,1.0,0.0,1.0,0.0005,0.0,     BLK 2480
     *     0.001,0.001,3*0.0,1.0,10*0.0/                                BLK 2490
C    ELEMENT GROUP FRACTIONAL RECOVERIES                                BLK 2500
      DATA FRECG /0.0,1.0,0.0005,0.0,1.0,1.0,1.0,3*0.0,0.0,1.0,0.0005,  BLK 2510
     *     0.0,1.0,1.0,4*0.0,0.0,1.0,0.0005,0.0,1.0,1.0,4*0.0,0.995,1.0,BLK 2520
     *     0.0005,0.999,1.0,0.2,0.6,1.0,0.0,0.0,0.0,1.0,0.0005,0.0,0.05,BLK 2530
     *     0.05,4*0.0,0.995,1.0,0.0005,0.9999,1.0,0.02,4*0.0,0.0,1.0,   BLK 2540
     *     0.0005,0.0,0.001,0.001,4*0.0,0.0,1.0,0.0005,0.0,0.001,0.001, BLK 2550
     *     4*0.0,0.0,1.0,0.0005,0.0,0.001,0.001,4*0.0,0.0,1.0,0.0005,   BLK 2560
     *     0.0,0.001,0.001,4*0.0,0.0,1.0,0.0005,0.0,0.001,0.001,4*0.0,  BLK 2570
     *     0.0,0.001,0.0005,0.0,3*1.0,3*0.0,0.0,0.0,0.0005,0.0,1.0,1.0, BLK 2580
     *     4*0.0,0.0,0.0,0.0005,0.0,1.0,1.0,0.0,0.0,0.9,0.0,60*0.0/     BLK 2590
C    ASSIGNMENT OF ELEMENTS TO GROUPS                                   BLK 2600
      DATA NGROUP /14,13,01,01,01,13,13,01,12,13,6*1,12,13,1,1,10*1,4*1,BLK 2610
     *     12,13,4*1,10*1,1,1,12,13,6*1,10*1,10*1,5*1,13,1,1,1,2,3,4,5, BLK 2620
     *     6,7,8,9,10,11,1/                                             BLK 2630
C    ELEMENTAL CHEMICAL TOXICITIES                                      BLK 2640
      DATA CEMTOX /3.50E+03,2.00E-01,5.00E+00,1.00E+00,1.00E+00,        BLK 2650
     *     4.00E+02,0.01,9.45E+05,1.0,1.0,1000.0,10.0,0.01,5.0,0.01,    BLK 2660
     *     50.0,0.15,10.0,1000.0,30.0,0.5,0.1,0.1,0.02,0.01,0.05,0.05,  BLK 2670
     *     0.05,0.01,0.05,0.2,0.5,0.01,0.01,3.0,40.0,50.0,10.0,0.001,   BLK 2680
     *     1.0,0.02,0.5,100.0,1.0,0.05,0.05,0.001,0.01,0.02,0.05,0.05,  BLK 2690
     *     0.2,10.0,150.0,5.0,0.5,1.0,2.0,1.0,0.2,1.0,0.2,0.2,0.2,0.5,  BLK 2700
     *     1.0,1.0,0.1,0.2,0.1,0.1,0.05,1.0,100.0,10.0,1.0,0.8,0.3,0.02,BLK 2710
     *     0.002,0.005,0.01,0.1,0.2,10.0,500.0,5.0,0.001,0.02,0.0005,   BLK 2720
     *     0.005,0.5,0.008,0.0008,0.04,0.5,0.005,0.01,0.01,0.0/         BLK 2730
      END                                                               BLK 2740
      INTEGER FUNCTION AREAD(IBS)                                       ARE   10
C                                                                       ARE   20
C    TYPE DEFINITIONS                                                   ARE   30
C                                                                       ARE   40
      INTEGER CARD, BLANK, QUOTE, QQPACK, OUTPT                         ARE   50
      DOUBLE PRECISION DNUM, UNIT, TENS                                 ARE   60
      LOGICAL WRDFLG, LSCAN, LALLC, LNEG, LBIN                          ARE   70
      LOGICAL LSKIP, CRDFLG, LASCN                                      ARE   80
C                                                                       ARE   90
C    ENTRY TO RETURN HOLLERITH DATA IN FOUR BYTE WORDS                  ARE  100
C    IBR = 0  RETURN NEXT 4 CHARACTERS ON CARD                          ARE  110
C    IBR = 1  READ A NEW CARD AND RETURN FIRST 4 CHARACTERS BEGINNING   ARE  120
C             WITH FIRST NON-BLANK CHARACTER; SAME AS IBR=4             ARE  130
C    IBR = 2  RETURN 4 CHARACTERS STARTING WITH THE FIRST NON BLANK     ARE  140
C             CHARACTER AND STOPPING WITH 2 CONSECUTIVE BLANKS          ARE  150
C    IBR = 3  RETURN NEXT 4 CHARACTERS, BUT STOP AT A DOUBLE BLANK      ARE  160
C    IBR = 4  SAME AS IBR=1                                             ARE  170
C    IBR =-1  RETURNS LAST CHARACTER SCANNED                            ARE  180
C                                                                       ARE  190
C    COMMONS                                                            ARE  200
C                                                                       ARE  210
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    ARE  220
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  ARE  230
     *     ITITLE(18)                                                   ARE  240
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      ARE  250
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           ARE  260
      DIMENSION RUOTE(1)                                                ARE  261
C                                                                       ARE  270
C                                                                       ARE  280
C    EQUIVALENCES                                                       ARE  290
C                                                                       ARE  300
      EQUIVALENCE (BLANK,CHAR(17)),(RUOTE(1),QUOTE)                     ARE  310
      EQUIVALENCE (E,CHAR(15)), (D,CHAR(14))                            ARE  320
C                                                                       ARE  330
C    DATA STATEMENTS                                                    ARE  340
C                                                                       ARE  350
      DATA RUOTE /1H*/, IBLNK /1/                                       ARE  360
C                                                                       ARE  370
C    END SPECIFICATION STATEMENTS                                       ARE  380
C                                                                       ARE  390
      OUTPT = IOUT                                                      ARE  400
      IBR = IBS                                                         ARE  410
      IF (IBS.EQ.1) IBR = 4                                             ARE  420
      IF (IBR.EQ.-1) GO TO 50                                           ARE  430
      IBLNK = 0                                                         ARE  440
      LSKIP = IBR.EQ.4                                                  ARE  450
      IF (IBR.NE.4) GO TO 20                                            ARE  460
   10 CONTINUE                                                          ARE  470
      ICHR = 0                                                          ARE  480
      READ (INPT,99999) CARD                                            ARE  490
      CRDFLG = .TRUE.                                                   ARE  500
      IF (CARD(1).NE.QUOTE) GO TO 20                                    ARE  510
      WRITE (OUTPT,99998) CARD                                          ARE  520
      GO TO 10                                                          ARE  530
C    FIND FIRST NON-BLANK CHARACTER                                     ARE  540
   20 CONTINUE                                                          ARE  550
      IF (.NOT.LSKIP) GO TO 30                                          ARE  560
      ICHR = ICHR + 1                                                   ARE  570
      IF (ICHR.GT.LCHR) GO TO 10                                        ARE  580
      IF (CARD(ICHR).EQ.BLANK) GO TO 20                                 ARE  590
      ICHR = ICHR - 1                                                   ARE  600
   30 CONTINUE                                                          ARE  610
      NCHR = 0                                                          ARE  620
      AREAD = BLANK                                                     ARE  630
      IRPT = 0                                                          ARE  640
      IF (ICHR.GE.LCHR) GO TO 10                                        ARE  650
      IF (IBR.EQ.3) GO TO 40                                            ARE  660
      WRDFLG = .FALSE.                                                  ARE  670
      INFG = 0                                                          ARE  680
   40 CONTINUE                                                          ARE  690
      NCHR = NCHR + 1                                                   ARE  700
      IF (WRDFLG .OR. NCHR.GT.4) GO TO 60                               ARE  710
      ICHR = ICHR + 1                                                   ARE  720
      IF (ICHR.GT.LCHR) RETURN                                          ARE  730
      INFG = INFG + 1                                                   ARE  740
      IF (CARD(ICHR).NE.BLANK) INFG = 0                                 ARE  750
      AREAD = QQPACK(AREAD,CARD(ICHR),NCHR)                             ARE  760
      WRDFLG = (INFG.GT.IBLNK .AND. IBR.GT.0)                           ARE  770
      GO TO 40                                                          ARE  780
C                                                                       ARE  790
C    RETURN SINGLE CHARACTER                                            ARE  800
C                                                                       ARE  810
   50 CONTINUE                                                          ARE  820
      AREAD = BLANK                                                     ARE  830
      IF (ICHR.GT.0 .AND. ICHR.LE.LCHR) AREAD = CARD(ICHR)              ARE  840
      RETURN                                                            ARE  850
C                                                                       ARE  860
   60 CONTINUE                                                          ARE  870
      RETURN                                                            ARE  880
C                                                                       ARE  890
99999 FORMAT (80A1)                                                     ARE  900
99998 FORMAT (1H0, 80A1)                                                ARE  910
      END                                                               ARE  920
      FUNCTION DREAD(IBR)                                               DRE   10
C                                                                       DRE   20
C    DREAD RETURNS A DOUBLE PRECISION FLOATING POINT NUMBER             DRE   30
C    SEE COMMENTS FOR AREAD                                             DRE   40
C                                                                       DRE   50
C    COMMONS                                                            DRE   60
C                                                                       DRE   70
      COMMON /OUT/ CUTOFF(30), CEMTOX(100), IUNIT, IOUT, JOUT, KOUT,    DRE   80
     *     MOUT, NREC, NTAB, NC, NPAGE, NTOC, NTOCA, NTOCB, BASIS(18),  DRE   90
     *     ITITLE(18)                                                   DRE  100
      COMMON /QRDBUF/ TENS(8), UNIT(10), DNUM, CARD(80), CHAR(28),      DRE  110
     *     ICHR, LCHR, IRPT, WRDFLG, LSCAN, LALLC, LNEG, INPT           DRE  120
      DIMENSION DHN(1),RUOTE(1)                                         DRE  121
C                                                                       DRE  130
C    TYPE DEFINITIONS                                                   DRE  140
C                                                                       DRE  150
      INTEGER OUTPT                                                     DRE  160
      DOUBLE PRECISION DREAD, DEXPN, TEN, DI, DNUM, UNIT, TENS, RLG     DRE  170
      LOGICAL WRDFLG, LSCAN, LALLC, LNEG, LBIN                          DRE  180
      LOGICAL ERRFLG, EXPFLG, NUMSGN, EXPSGN, PERFLG, ZFLG, CRDFLG      DRE  190
C                                                                       DRE  200
C    EQUIVALENCES                                                       DRE  210
C                                                                       DRE  220
      EQUIVALENCE (BLANK,CHAR(17)),(DHN(1),CHN),(RUOTE(1),QUOTE)        DRE  230
      EQUIVALENCE (TEN,UNIT(2)), (E,CHAR(15)), (D,CHAR(14))             DRE  240
C                                                                       DRE  250
C    DATA STATEMENTS                                                    DRE  260
C                                                                       DRE  270
      DATA DHN /1HN/, NEXP /3/, MANT /17/                               DRE  280
      DATA RUOTE /1H*/                                                  DRE  290
C                                                                       DRE  300
C    END SPECIFICATION STATEMENTS                                       DRE  310
C                                                                       DRE  320
      OUTPT = IOUT                                                      DRE  330
      INFG = 0                                                          DRE  340
      NSIG = MANT                                                       DRE  350
      EXPFLG = .FALSE.                                                  DRE  360
      ERRFLG = .FALSE.                                                  DRE  370
      WRDFLG = .FALSE.                                                  DRE  380
      ZFLG = .FALSE.                                                    DRE  390
      IF (IBR.NE.1) GO TO 20                                            DRE  400
   10 CONTINUE                                                          DRE  410
      IF (INFG.GT.0) GO TO 90                                           DRE  420
      ICHR = 0                                                          DRE  430
      IRPT = 0                                                          DRE  440
      READ (INPT,99999) CARD                                            DRE  450
      CRDFLG = .TRUE.                                                   DRE  460
      IF (CARD(1).NE.QUOTE) GO TO 20                                    DRE  470
C    SKIP CARDS WITH QUOTES IN COLUMN ONE                               DRE  480
      WRITE (OUTPT,99996) CARD                                          DRE  490
      GO TO 10                                                          DRE  500
   20 CONTINUE                                                          DRE  510
      IF (IRPT.GT.0) GO TO 200                                          DRE  520
      DNUM = 0                                                          DRE  530
      IPER = 0                                                          DRE  540
      NUM = 0                                                           DRE  550
      NUMSGN = .FALSE.                                                  DRE  560
      EXPSGN = .FALSE.                                                  DRE  570
      PERFLG = .FALSE.                                                  DRE  580
      LNEG = .FALSE.                                                    DRE  590
   30 CONTINUE                                                          DRE  600
      ICHR = ICHR + 1                                                   DRE  610
      IF (ICHR.GT.LCHR) GO TO 10                                        DRE  620
      DO 40 I=1,28                                                      DRE  640
         IF (CHAR(I).EQ.CARD(ICHR)) GO TO (50, 50, 50, 50, 50, 50, 50,  DRE  650
     *        50, 50, 50, 60, 60, 60, 70, 70, 60, 80, 90, 150, 150,     DRE  660
     *        150, 110, 110, 130, 170, 180, 90, 160), I                 DRE  670
   40 CONTINUE                                                          DRE  680
      GO TO 60                                                          DRE  690
C                                                                       DRE  700
C    NUMERIC CHARACTER                                                  DRE  710
C                                                                       DRE  720
   50 CONTINUE                                                          DRE  730
      INFG = INFG + 1                                                   DRE  740
      IF (EXPFLG) GO TO 190                                             DRE  750
      IF (INFG.GT.NSIG) GO TO 30                                        DRE  760
      DI = I - 1                                                        DRE  770
      DNUM = TEN*DNUM + DI                                              DRE  780
      IF (PERFLG) IPER = IPER + 1                                       DRE  790
      GO TO 30                                                          DRE  800
C                                                                       DRE  810
C    ILLEGAL CHARACTER                                                  DRE  820
C                                                                       DRE  830
   60 CONTINUE                                                          DRE  840
      IF (LALLC) GO TO 90                                               DRE  850
      IF (ERRFLG) GO TO 30                                              DRE  860
      ERRFLG = .TRUE.                                                   DRE  870
      IF (CRDFLG) WRITE (OUTPT,99998) CARD                              DRE  880
      CRDFLG = .FALSE.                                                  DRE  890
      WRITE (OUTPT,99997) ICHR, CARD(ICHR)                              DRE  900
      GO TO 30                                                          DRE  910
C                                                                       DRE  920
C    EXPONENT FLAG - D OR E                                             DRE  930
C                                                                       DRE  940
   70 CONTINUE                                                          DRE  950
      EXPFLG = .TRUE.                                                   DRE  960
      NSIG = NEXP                                                       DRE  970
      INFG = 1                                                          DRE  980
      IF (ICHR+1.LE.LCHR .AND. CARD(ICHR+1).EQ.BLANK) ICHR = ICHR + 1   DRE  990
      GO TO 30                                                          DRE 1000
C                                                                       DRE 1010
C    BLANK                                                              DRE 1020
C                                                                       DRE 1030
   80 CONTINUE                                                          DRE 1040
      IF (INFG.EQ.0) GO TO 30                                           DRE 1050
C                                                                       DRE 1060
C    COMMA - FIELD TERMINATOR                                           DRE 1070
C                                                                       DRE 1080
   90 CONTINUE                                                          DRE 1090
      IF (NUMSGN) DNUM = -DNUM                                          DRE 1100
      IF (EXPSGN) NUM = -NUM                                            DRE 1110
      I = NUM - IPER                                                    DRE 1120
      IU = IABS(I)                                                      DRE 1130
      RLG = 1.0                                                         DRE 1140
      IF (IU.LT.36) GO TO 100                                           DRE 1150
      RLG = 10.0**(IU-35)                                               DRE 1160
      IU = 35                                                           DRE 1170
  100 CONTINUE                                                          DRE 1180
      IT = IU/10                                                        DRE 1190
      IU = IU - 10*IT                                                   DRE 1200
      DEXPN = UNIT(IU+1)*TENS(IT+1)                                     DRE 1210
      EXPSGN = I.LT.0                                                   DRE 1220
      IF (EXPSGN) DNUM = (DNUM/DEXPN)/RLG                               DRE 1230
      IF (.NOT.EXPSGN) DNUM = DNUM*DEXPN*RLG                            DRE 1240
      GO TO 200                                                         DRE 1250
C                                                                       DRE 1260
C    PLUS SIGN                                                          DRE 1270
C                                                                       DRE 1280
  110 CONTINUE                                                          DRE 1290
      IF (INFG.GT.0) GO TO 120                                          DRE 1300
      INFG = INFG + 1                                                   DRE 1310
      GO TO 30                                                          DRE 1320
  120 CONTINUE                                                          DRE 1330
      EXPFLG = .TRUE.                                                   DRE 1340
      NSIG = NEXP                                                       DRE 1350
      INFG = 1                                                          DRE 1360
      GO TO 30                                                          DRE 1370
C                                                                       DRE 1380
C    MINUS SIGN                                                         DRE 1390
C                                                                       DRE 1400
  130 CONTINUE                                                          DRE 1410
      IF (INFG.GT.0) GO TO 140                                          DRE 1420
      NUMSGN = .TRUE.                                                   DRE 1430
      INFG = INFG + 1                                                   DRE 1440
      GO TO 30                                                          DRE 1450
  140 CONTINUE                                                          DRE 1460
      EXPFLG = .TRUE.                                                   DRE 1470
      EXPSGN = .TRUE.                                                   DRE 1480
      NSIG = NEXP                                                       DRE 1490
      INFG = 1                                                          DRE 1500
      GO TO 30                                                          DRE 1510
C                                                                       DRE 1520
C    REPEAT FLAG - R,*, OR $                                            DRE 1530
C                                                                       DRE 1540
  150 CONTINUE                                                          DRE 1550
      IRPT = DNUM                                                       DRE 1560
      DNUM = 0                                                          DRE 1570
      NUMSGN = .FALSE.                                                  DRE 1580
      EXPSGN = .FALSE.                                                  DRE 1590
      EXPFLG = .FALSE.                                                  DRE 1600
      PERFLG = .FALSE.                                                  DRE 1610
      IPER = 0                                                          DRE 1620
      IF (ZFLG) GO TO 30                                                DRE 1630
      INFG = 0                                                          DRE 1640
      GO TO 30                                                          DRE 1650
C                                                                       DRE 1660
C    P - ALTERNATING SIGN REPEAT                                        DRE 1670
C                                                                       DRE 1680
  160 CONTINUE                                                          DRE 1690
      LNEG = .TRUE.                                                     DRE 1700
      GO TO 150                                                         DRE 1710
C                                                                       DRE 1720
C    Z - ZERO FIELD                                                     DRE 1730
C                                                                       DRE 1740
  170 CONTINUE                                                          DRE 1750
      ZFLG = .TRUE.                                                     DRE 1760
      GO TO 150                                                         DRE 1770
C                                                                       DRE 1780
C    DECIMAL POINT                                                      DRE 1790
C                                                                       DRE 1800
  180 CONTINUE                                                          DRE 1810
      PERFLG = .TRUE.                                                   DRE 1820
      GO TO 30                                                          DRE 1830
  190 CONTINUE                                                          DRE 1840
      NUM = 10*NUM + I - 1                                              DRE 1850
      IF(INFG.GT.NSIG)WRITE(OUTPT,99995)NUM,ICHR,CARD                           
      GO TO 30                                                          DRE 1870
  200 CONTINUE                                                          DRE 1880
      DREAD = DNUM                                                      DRE 1890
      IRPT = IRPT - 1                                                   DRE 1900
      IF (LNEG) DNUM = -DNUM                                            DRE 1910
      GO TO 210                                                         DRE 1920
  210 CONTINUE                                                          DRE 1930
      IF (IRPT.GT.0) RETURN                                             DRE 1940
      IF (.NOT.LSCAN) RETURN                                            DRE 1950
      IF (ICHR.LT.LCHR) GO TO 230                                       DRE 1960
  220 CONTINUE                                                          DRE 1970
      READ (INPT,99999) CARD                                            DRE 1980
      ICHR = 0                                                          DRE 1990
      CRDFLG = .TRUE.                                                   DRE 2000
      IF (CARD(1).NE.QUOTE) GO TO 230                                   DRE 2010
      WRITE (OUTPT,99996) CARD                                          DRE 2020
      GO TO 220                                                         DRE 2030
  230 CONTINUE                                                          DRE 2040
      IE = ICHR + 1                                                     DRE 2050
      IF (CARD(IE).EQ.BLANK) GO TO 240                                  DRE 2060
      IN = IE + 1                                                       DRE 2070
      ID = IN + 1                                                       DRE 2080
      IF (ID.GT.LCHR) RETURN                                            DRE 2090
      IF (CARD(IE).NE.E .OR. CARD(IN).NE.CHN .OR. CARD(ID).NE.D) RETURN DRE 2100
      RETURN                                                            DRE 2110
  240 CONTINUE                                                          DRE 2120
      ICHR = ICHR + 1                                                   DRE 2130
      GO TO 210                                                         DRE 2140
C                                                                       DRE 2150
C    FORMATS                                                            DRE 2160
C                                                                       DRE 2170
99995 FORMAT(1H0,29H***ERROR IN INPUT: EXPONENT= ,I6,                           
     *     20H IS TOO LONG. ICHR= ,I6,/,1H ,80A1)                               
99999 FORMAT (80A1)                                                     DRE 2180
99998 FORMAT (54H0*****ERROR IN INPUT. CARD IMAGE PRINTED ON NEXT LINE ,DRE 2190
     *     6H***** /10X, 80A1)                                          DRE 2200
99997 FORMAT (37H ON THE ABOVE CARD, CHARACTER NUMBER , I2, 8H (IMAGE=, DRE 2210
     *     A1, 15H) IS NOT VALID.)                                      DRE 2220
99996 FORMAT (1H0, 80A1)                                                DRE 2230
      END                                                               DRE 2240
      INTEGER FUNCTION QQPACK(A, B, I)                                  QQP   10
      INTEGER A, B, C, D                                                QQP   20
      LOGICAL*1 L(4),L1                                                 QQP   30
      EQUIVALENCE (C,L1), (D,L(1))                                      QQP   40
      C = B                                                             QQP   50
      D = A                                                             QQP   60
      L(I) = L1                                                         QQP   70
      QQPACK = D                                                        QQP   80
C    ON A CDC COMPUTER SUBSTITUTE THE FOLLOWING SIX  CARDS FOR THE      QQP   90
C    PRECEDING SIX CARDS                                                QQP  100
C      II=6*(I-1)                                                       QQP  110
C      IMSK=MASK(II).OR..NOT.MASK(II+6)                                 QQP  111
C      C= .NOT.IMSK.AND.SHIFT(B,-II)                                    QQP  112
C      D=IMSK.AND.A                                                     QQP  120
C      A=C.OR.D                                                         QQP  130
C      QQPACK=A                                                         QQP  140
      RETURN                                                            QQP  150
      END                                                               QQP  160
      FUNCTION READ(IBR)                                                REA   10
      DOUBLE PRECISION DREAD                                            REA   20
      READ = DREAD(IBR)                                                 REA   30
      RETURN                                                            REA   40
      END                                                               REA   50
      FUNCTION IREAD(IBR)                                               IRE   10
      DOUBLE PRECISION DREAD                                            IRE   20
      IREAD = DREAD(IBR)                                                IRE   30
      RETURN                                                            IRE   40
      END                                                               IRE   50
      SUBROUTINE XSEC01(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PWR-U  ;XSDRNPM/CITATION,XSEC FROM CITATION                        XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,18/                                           XSE  110
      DATA DBUP /0.0,5.708E-01,1.151E00,1.738E00,2.329E00,2.923E00,     XSE  120
     *     2.923E00,3.461E00,3.995E00,4.528E00,5.061E00,5.593E00,       XSE  130
     *     5.593E00,6.077E00,6.558E00,7.038E00,7.518E00,8.000E00,0.0,13*XSE  140
     *     0.0/                                                         XSE  150
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  160
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  170
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  180
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  190
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  200
      DATA X01 /1.935E01,1.935E01,1.935E01,1.937E01,1.939E01,1.941E01,  XSE  210
     *     1.926E01,1.926E01,1.926E01,1.928E01,1.930E01,1.932E01,       XSE  220
     *     1.933E01,1.932E01,1.933E01,1.934E01,1.936E01,1.939E01,0.0,13*XSE  230
     *     0.0/                                                         XSE  240
      DATA X02 /1.027E01,1.023E01,1.021E01,1.020E01,1.020E01,1.021E01,  XSE  250
     *     1.042E01,1.037E01,1.035E01,1.034E01,1.034E01,1.035E01,       XSE  260
     *     1.065E01,1.060E01,1.058E01,1.057E01,1.057E01,1.058E01,0.0,13*XSE  270
     *     0.0/                                                         XSE  280
      DATA X03 /4.578E01,4.551E01,4.534E01,4.526E01,4.525E01,4.531E01,  XSE  290
     *     4.660E01,4.632E01,4.615E01,4.607E01,4.606E01,4.612E01,       XSE  300
     *     4.781E01,4.752E01,4.735E01,4.726E01,4.726E01,4.732E01,0.0,13*XSE  310
     *     0.0/                                                         XSE  320
      DATA X04 /7.691E00,7.712E00,7.730E00,7.746E00,7.758E00,7.768E00,  XSE  330
     *     7.573E00,7.593E00,7.611E00,7.626E00,7.638E00,7.648E00,       XSE  340
     *     7.512E00,7.532E00,7.550E00,7.564E00,7.576E00,7.586E00,0.0,13*XSE  350
     *     0.0/                                                         XSE  360
      DATA X05 /8.878E-01,8.883E-01,8.889E-01,8.896E-01,8.903E-01,      XSE  370
     *     8.911E-01,8.940E-01,8.944E-01,8.950E-01,8.957E-01,8.964E-01, XSE  380
     *     8.973E-01,9.027E-01,9.031E-01,9.036E-01,9.043E-01,9.051E-01, XSE  390
     *     9.059E-01,0.0,13*0.0/                                        XSE  400
      DATA X06 /3.285E01,3.272E01,3.263E01,3.257E01,3.255E01,3.255E01,  XSE  410
     *     3.325E01,3.312E01,3.302E01,3.296E01,3.294E01,3.294E01,       XSE  420
     *     3.333E01,3.320E01,3.310E01,3.304E01,3.301E01,3.301E01,0.0,13*XSE  430
     *     0.0/                                                         XSE  440
      DATA X07 /3.314E01,3.293E01,3.281E01,3.275E01,3.275E01,3.280E01,  XSE  450
     *     3.418E01,3.396E01,3.383E01,3.377E01,3.377E01,3.382E01,       XSE  460
     *     3.539E01,3.516E01,3.503E01,3.497E01,3.497E01,3.502E01,0.0,13*XSE  470
     *     0.0/                                                         XSE  480
      DATA X08 /2.400E00,2.396E00,2.394E00,2.393E00,2.394E00,2.395E00,  XSE  490
     *     2.434E00,2.430E00,2.427E00,2.427E00,2.427E00,2.429E00,       XSE  500
     *     2.472E00,2.468E00,2.465E00,2.464E00,2.465E00,2.466E00,0.0,13*XSE  510
     *     0.0/                                                         XSE  520
      DATA X09 /6.564E01,6.520E01,6.493E01,6.481E01,6.480E01,6.491E01,  XSE  530
     *     6.211E01,6.169E01,6.144E01,6.133E01,6.132E01,6.142E01,       XSE  540
     *     6.105E01,6.065E01,6.040E01,6.029E01,6.028E01,6.038E01,0.0,13*XSE  550
     *     0.0/                                                         XSE  560
      DATA X10 /1.154E02,1.146E02,1.141E02,1.139E02,1.139E02,1.141E02,  XSE  570
     *     1.110E02,1.103E02,1.098E02,1.096E02,1.096E02,1.098E02,       XSE  580
     *     1.103E02,1.095E02,1.091E02,1.089E02,1.089E02,1.090E02,0.0,13*XSE  590
     *     0.0/                                                         XSE  600
      DATA X11 /1.850E02,1.840E02,1.828E02,1.814E02,1.800E02,1.787E02,  XSE  610
     *     1.300E02,1.294E02,1.285E02,1.276E02,1.267E02,1.259E02,       XSE  620
     *     1.073E02,1.068E02,1.061E02,1.054E02,1.047E02,1.041E02,0.0,13*XSE  630
     *     0.0/                                                         XSE  640
      DATA X12 /4.027E01,4.001E01,3.985E01,3.977E01,3.977E01,3.983E01,  XSE  650
     *     3.960E01,3.934E01,3.918E01,3.911E01,3.911E01,3.917E01,       XSE  660
     *     3.987E01,3.960E01,3.945E01,3.937E01,3.937E01,3.944E01,0.0,13*XSE  670
     *     0.0/                                                         XSE  680
      DATA X13 /1.213E02,1.205E02,1.201E02,1.199E02,1.199E02,1.201E02,  XSE  690
     *     1.202E02,1.194E02,1.190E02,1.188E02,1.188E02,1.190E02,       XSE  700
     *     1.214E02,1.206E02,1.202E02,1.200E02,1.200E02,1.201E02,0.0,13*XSE  710
     *     0.0/                                                         XSE  720
      DATA X14 /3.177E01,3.186E01,3.195E01,3.201E01,3.207E01,3.211E01,  XSE  730
     *     3.167E01,3.176E01,3.184E01,3.191E01,3.197E01,3.201E01,       XSE  740
     *     3.121E01,3.130E01,3.138E01,3.145E01,3.150E01,3.155E01,0.0,13*XSE  750
     *     0.0/                                                         XSE  760
      DATA X15 /1.162E02,1.155E02,1.151E02,1.148E02,1.147E02,1.148E02,  XSE  770
     *     1.121E02,1.115E02,1.110E02,1.108E02,1.107E02,1.107E02,       XSE  780
     *     1.105E02,1.098E02,1.094E02,1.091E02,1.090E02,1.091E02,0.0,13*XSE  790
     *     0.0/                                                         XSE  800
      DATA X16 /1.436E01,1.428E01,1.422E01,1.419E01,1.418E01,1.418E01,  XSE  810
     *     1.386E01,1.378E01,1.373E01,1.369E01,1.368E01,1.368E01,       XSE  820
     *     1.365E01,1.357E01,1.352E01,1.349E01,1.348E01,1.348E01,0.0,13*XSE  830
     *     0.0/                                                         XSE  840
      DATA X17 /2.508E00,2.498E00,2.486E00,2.472E00,2.459E00,2.446E00,  XSE  850
     *     2.549E00,2.540E00,2.527E00,2.513E00,2.499E00,2.487E00,       XSE  860
     *     2.565E00,2.556E00,2.543E00,2.529E00,2.515E00,2.503E00,0.0,13*XSE  870
     *     0.0/                                                         XSE  880
      DATA X18 /4.765E01,4.747E01,4.723E01,4.697E01,4.671E01,4.648E01,  XSE  890
     *     4.844E01,4.825E01,4.801E01,4.775E01,4.749E01,4.725E01,       XSE  900
     *     4.874E01,4.856E01,4.831E01,4.805E01,4.779E01,4.755E01,0.0,13*XSE  910
     *     0.0/                                                         XSE  920
      DATA X19 /5.422E00,5.426E00,5.432E00,5.438E00,5.444E00,5.450E00,  XSE  930
     *     5.450E00,5.454E00,5.460E00,5.466E00,5.472E00,5.479E00,       XSE  940
     *     5.494E00,5.498E00,5.503E00,5.509E00,5.516E00,5.522E00,0.0,13*XSE  950
     *     0.0/                                                         XSE  960
      DATA X20 /1.367E01,1.367E01,1.368E01,1.369E01,1.371E01,1.372E01,  XSE  970
     *     1.375E01,1.376E01,1.376E01,1.378E01,1.379E01,1.380E01,       XSE  980
     *     1.387E01,1.387E01,1.388E01,1.389E01,1.390E01,1.392E01,0.0,13*XSE  990
     *     0.0/                                                         XSE 1000
      NPT = MPT                                                         XSE 1010
      NXSEC = MXSEC                                                     XSE 1020
      DO 10 I=1,30                                                      XSE 1030
         F(I,01) = X01(I)                                               XSE 1040
         F(I,02) = X02(I)                                               XSE 1050
         F(I,03) = X03(I)                                               XSE 1060
         F(I,04) = X04(I)                                               XSE 1070
         F(I,05) = X05(I)                                               XSE 1080
         F(I,06) = X06(I)                                               XSE 1090
         F(I,07) = X07(I)                                               XSE 1100
         F(I,08) = X08(I)                                               XSE 1110
         F(I,09) = X09(I)                                               XSE 1120
         F(I,10) = X10(I)                                               XSE 1130
         F(I,11) = X11(I)                                               XSE 1140
         F(I,12) = X12(I)                                               XSE 1150
         F(I,13) = X13(I)                                               XSE 1160
         F(I,14) = X14(I)                                               XSE 1170
         F(I,15) = X15(I)                                               XSE 1180
         F(I,16) = X16(I)                                               XSE 1190
         F(I,17) = X17(I)                                               XSE 1200
         F(I,18) = X18(I)                                               XSE 1210
         F(I,19) = X19(I)                                               XSE 1220
         F(I,20) = X20(I)                                               XSE 1230
         CBUP(I) = DBUP(I)                                              XSE 1240
   10 CONTINUE                                                          XSE 1250
      DO 30 I=1,20                                                      XSE 1260
         DO 20 J=1,7                                                    XSE 1270
            NFUD(J,I) = MFUD(J,I)                                       XSE 1280
   20    CONTINUE                                                       XSE 1290
   30 CONTINUE                                                          XSE 1300
      RETURN                                                            XSE 1310
      END                                                               XSE 1320
      SUBROUTINE XSEC02(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PWR-PUU ;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,18/                                           XSE  110
      DATA DBUP /0.0,8.542E-01,1.662E00,2.420E00,3.136E00,3.821E00,     XSE  120
     *     3.821E00,4.299E00,4.786E00,5.281E00,5.784E00,6.290E00,       XSE  130
     *     6.290E00,6.556E00,6.855E00,7.189E00,7.554E00,7.943E00,0.0,13*XSE  140
     *     0.0/                                                         XSE  150
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  160
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  170
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  180
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  190
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  200
      DATA X01 /2.025E01,2.012E01,2.007E01,2.006E01,2.007E01,2.009E01,  XSE  210
     *     1.995E01,1.997E01,1.998E01,2.000E01,2.003E01,2.005E01,       XSE  220
     *     2.015E01,2.016E01,2.017E01,2.017E01,2.018E01,2.020E01,0.0,13*XSE  230
     *     0.0/                                                         XSE  240
      DATA X02 /1.049E01,1.023E01,1.007E01,9.987E00,9.942E00,9.927E00,  XSE  250
     *     9.987E00,9.983E00,9.981E00,9.986E00,1.000E01,1.002E01,       XSE  260
     *     1.022E01,1.023E01,1.024E01,1.025E01,1.026E01,1.029E01,0.0,13*XSE  270
     *     0.0/                                                         XSE  280
      DATA X03 /4.659E01,4.498E01,4.402E01,4.345E01,4.315E01,4.303E01,  XSE  290
     *     4.329E01,4.324E01,4.321E01,4.323E01,4.331E01,4.343E01,       XSE  300
     *     4.435E01,4.442E01,4.446E01,4.453E01,4.463E01,4.476E01,0.0,13*XSE  310
     *     0.0/                                                         XSE  320
      DATA X04 /8.164E00,8.229E00,8.284E00,8.327E00,8.359E00,8.383E00,  XSE  330
     *     8.284E00,8.300E00,8.311E00,8.319E00,8.326E00,8.331E00,       XSE  340
     *     8.315E00,8.315E00,8.312E00,8.309E00,8.308E00,8.308E00,0.0,13*XSE  350
     *     0.0/                                                         XSE  360
      DATA X05 /9.316E-01,9.304E-01,9.308E-01,9.316E-01,9.327E-01,      XSE  370
     *     9.339E-01,9.379E-01,9.388E-01,9.394E-01,9.400E-01,9.407E-01, XSE  380
     *     9.415E-01,9.512E-01,9.514E-01,9.513E-01,9.514E-01,9.516E-01, XSE  390
     *     9.521E-01,0.0,13*0.0/                                        XSE  400
      DATA X06 /3.351E01,3.281E01,3.235E01,3.205E01,3.188E01,3.179E01,  XSE  410
     *     3.205E01,3.201E01,3.197E01,3.195E01,3.196E01,3.199E01,       XSE  420
     *     3.228E01,3.229E01,3.228E01,3.229E01,3.231E01,3.235E01,0.0,13*XSE  430
     *     0.0/                                                         XSE  440
      DATA X07 /3.388E01,3.263E01,3.189E01,3.146E01,3.124E01,3.115E01,  XSE  450
     *     3.159E01,3.156E01,3.154E01,3.157E01,3.163E01,3.173E01,       XSE  460
     *     3.259E01,3.265E01,3.269E01,3.274E01,3.282E01,3.293E01,0.0,13*XSE  470
     *     0.0/                                                         XSE  480
      DATA X08 /2.423E00,2.397E00,2.382E00,2.374E00,2.370E00,2.369E00,  XSE  490
     *     2.385E00,2.385E00,2.385E00,2.386E00,2.387E00,2.390E00,       XSE  500
     *     2.414E00,2.416E00,2.417E00,2.418E00,2.420E00,2.422E00,0.0,13*XSE  510
     *     0.0/                                                         XSE  520
      DATA X09 /6.614E01,6.357E01,6.205E01,6.116E01,6.069E01,6.051E01,  XSE  530
     *     5.747E01,5.740E01,5.737E01,5.740E01,5.752E01,5.771E01,       XSE  540
     *     5.706E01,5.716E01,5.722E01,5.732E01,5.747E01,5.767E01,0.0,13*XSE  550
     *     0.0/                                                         XSE  560
      DATA X10 /1.165E02,1.119E02,1.092E02,1.076E02,1.068E02,1.064E02,  XSE  570
     *     1.024E02,1.023E02,1.022E02,1.023E02,1.025E02,1.028E02,       XSE  580
     *     1.024E02,1.026E02,1.027E02,1.029E02,1.032E02,1.035E02,0.0,13*XSE  590
     *     0.0/                                                         XSE  600
      DATA X11 /1.948E02,1.936E02,1.906E02,1.874E02,1.847E02,1.825E02,  XSE  610
     *     1.371E02,1.362E02,1.352E02,1.343E02,1.335E02,1.328E02,       XSE  620
     *     1.141E02,1.137E02,1.133E02,1.129E02,1.126E02,1.123E02,0.0,13*XSE  630
     *     0.0/                                                         XSE  640
      DATA X12 /4.067E01,3.910E01,3.818E01,3.764E01,3.736E01,3.725E01,  XSE  650
     *     3.644E01,3.640E01,3.638E01,3.640E01,3.648E01,3.660E01,       XSE  660
     *     3.682E01,3.688E01,3.693E01,3.699E01,3.709E01,3.721E01,0.0,13*XSE  670
     *     0.0/                                                         XSE  680
      DATA X13 /1.223E02,1.179E02,1.153E02,1.138E02,1.130E02,1.127E02,  XSE  690
     *     1.109E02,1.108E02,1.108E02,1.108E02,1.110E02,1.114E02,       XSE  700
     *     1.123E02,1.125E02,1.127E02,1.128E02,1.131E02,1.135E02,0.0,13*XSE  710
     *     0.0/                                                         XSE  720
      DATA X14 /3.144E01,3.170E01,3.192E01,3.209E01,3.222E01,3.231E01,  XSE  730
     *     3.223E01,3.229E01,3.234E01,3.237E01,3.240E01,3.242E01,       XSE  740
     *     3.222E01,3.222E01,3.221E01,3.220E01,3.219E01,3.219E01,0.0,13*XSE  750
     *     0.0/                                                         XSE  760
      DATA X15 /1.059E02,1.027E02,1.006E02,9.927E01,9.851E01,9.814E01,  XSE  770
     *     9.526E01,9.509E01,9.496E01,9.492E01,9.498E01,9.515E01,       XSE  780
     *     9.465E01,9.473E01,9.476E01,9.482E01,9.496E01,9.515E01,0.0,13*XSE  790
     *     0.0/                                                         XSE  800
      DATA X16 /1.309E01,1.269E01,1.243E01,1.227E01,1.218E01,1.213E01,  XSE  810
     *     1.177E01,1.175E01,1.174E01,1.173E01,1.174E01,1.176E01,       XSE  820
     *     1.170E01,1.171E01,1.171E01,1.172E01,1.174E01,1.176E01,0.0,13*XSE  830
     *     0.0/                                                         XSE  840
      DATA X17 /2.559E00,2.544E00,2.512E00,2.482E00,2.456E00,2.435E00,  XSE  850
     *     2.520E00,2.507E00,2.493E00,2.480E00,2.469E00,2.460E00,       XSE  860
     *     2.521E00,2.515E00,2.507E00,2.500E00,2.493E00,2.487E00,0.0,13*XSE  870
     *     0.0/                                                         XSE  880
      DATA X18 /4.862E01,4.833E01,4.774E01,4.715E01,4.666E01,4.626E01,  XSE  890
     *     4.788E01,4.762E01,4.736E01,4.712E01,4.691E01,4.673E01,       XSE  900
     *     4.790E01,4.778E01,4.764E01,4.750E01,4.737E01,4.725E01,0.0,13*XSE  910
     *     0.0/                                                         XSE  920
      DATA X19 /5.558E00,5.551E00,5.556E00,5.563E00,5.572E00,5.581E00,  XSE  930
     *     5.591E00,5.598E00,5.603E00,5.607E00,5.613E00,5.618E00,       XSE  940
     *     5.667E00,5.668E00,5.668E00,5.668E00,5.670E00,5.673E00,0.0,13*XSE  950
     *     0.0/                                                         XSE  960
      DATA X20 /1.400E01,1.397E01,1.396E01,1.397E01,1.399E01,1.401E01,  XSE  970
     *     1.404E01,1.406E01,1.407E01,1.408E01,1.409E01,1.410E01,       XSE  980
     *     1.422E01,1.423E01,1.423E01,1.423E01,1.423E01,1.424E01,0.0,13*XSE  990
     *     0.0/                                                         XSE 1000
      NPT = MPT                                                         XSE 1010
      NXSEC = MXSEC                                                     XSE 1020
      DO 10 I=1,30                                                      XSE 1030
         F(I,01) = X01(I)                                               XSE 1040
         F(I,02) = X02(I)                                               XSE 1050
         F(I,03) = X03(I)                                               XSE 1060
         F(I,04) = X04(I)                                               XSE 1070
         F(I,05) = X05(I)                                               XSE 1080
         F(I,06) = X06(I)                                               XSE 1090
         F(I,07) = X07(I)                                               XSE 1100
         F(I,08) = X08(I)                                               XSE 1110
         F(I,09) = X09(I)                                               XSE 1120
         F(I,10) = X10(I)                                               XSE 1130
         F(I,11) = X11(I)                                               XSE 1140
         F(I,12) = X12(I)                                               XSE 1150
         F(I,13) = X13(I)                                               XSE 1160
         F(I,14) = X14(I)                                               XSE 1170
         F(I,15) = X15(I)                                               XSE 1180
         F(I,16) = X16(I)                                               XSE 1190
         F(I,17) = X17(I)                                               XSE 1200
         F(I,18) = X18(I)                                               XSE 1210
         F(I,19) = X19(I)                                               XSE 1220
         F(I,20) = X20(I)                                               XSE 1230
         CBUP(I) = DBUP(I)                                              XSE 1240
   10 CONTINUE                                                          XSE 1250
      DO 30 I=1,20                                                      XSE 1260
         DO 20 J=1,7                                                    XSE 1270
            NFUD(J,I) = MFUD(J,I)                                       XSE 1280
   20    CONTINUE                                                       XSE 1290
   30 CONTINUE                                                          XSE 1300
      RETURN                                                            XSE 1310
      END                                                               XSE 1320
      SUBROUTINE XSEC03(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PWR-PUPU;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,18/                                           XSE  110
      DATA DBUP /0.0,3.065E-01,6.525E-01,1.042E00,1.469E00,1.929E00,    XSE  120
     *     1.929E00,2.455E00,2.997E00,3.554E00,4.123E00,4.701E00,       XSE  130
     *     4.701E00,5.501E00,6.267E00,6.999E00,7.700E00,8.377E00,0.0,13*XSE  140
     *     0.0/                                                         XSE  150
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  160
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  170
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  180
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  190
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  200
      DATA X01 /1.657E01,1.658E01,1.659E01,1.661E01,1.664E01,1.667E01,  XSE  210
     *     1.720E01,1.726E01,1.731E01,1.736E01,1.742E01,1.748E01,       XSE  220
     *     1.764E01,1.775E01,1.785E01,1.795E01,1.805E01,1.814E01,0.0,13*XSE  230
     *     0.0/                                                         XSE  240
      DATA X02 /5.946E00,5.970E00,5.993E00,6.021E00,6.056E00,6.097E00,  XSE  250
     *     6.436E00,6.496E00,6.555E00,6.618E00,6.687E00,6.760E00,       XSE  260
     *     6.944E00,7.033E00,7.129E00,7.229E00,7.332E00,7.437E00,0.0,13*XSE  270
     *     0.0/                                                         XSE  280
      DATA X03 /2.027E01,2.042E01,2.056E01,2.073E01,2.094E01,2.118E01,  XSE  290
     *     2.274E01,2.308E01,2.342E01,2.379E01,2.419E01,2.462E01,       XSE  300
     *     2.556E01,2.605E01,2.659E01,2.715E01,2.774E01,2.835E01,0.0,13*XSE  310
     *     0.0/                                                         XSE  320
      DATA X04 /8.353E00,8.350E00,8.344E00,8.338E00,8.334E00,8.332E00,  XSE  330
     *     8.500E00,8.502E00,8.500E00,8.498E00,8.496E00,8.494E00,       XSE  340
     *     8.510E00,8.527E00,8.536E00,8.542E00,8.546E00,8.548E00,0.0,13*XSE  350
     *     0.0/                                                         XSE  360
      DATA X05 /8.567E-01,8.570E-01,8.570E-01,8.571E-01,8.575E-01,      XSE  370
     *     8.581E-01,8.777E-01,8.790E-01,8.799E-01,8.809E-01,8.820E-01, XSE  380
     *     8.833E-01,8.913E-01,8.939E-01,8.962E-01,8.984E-01,9.005E-01, XSE  390
     *     9.025E-01,0.0,13*0.0/                                        XSE  400
      DATA X06 /2.239E01,2.248E01,2.255E01,2.265E01,2.277E01,2.292E01,  XSE  410
     *     2.349E01,2.369E01,2.389E01,2.410E01,2.433E01,2.457E01,       XSE  420
     *     2.459E01,2.489E01,2.520E01,2.552E01,2.585E01,2.619E01,0.0,13*XSE  430
     *     0.0/                                                         XSE  440
      DATA X07 /1.342E01,1.354E01,1.365E01,1.379E01,1.395E01,1.413E01,  XSE  450
     *     1.533E01,1.560E01,1.587E01,1.615E01,1.647E01,1.680E01,       XSE  460
     *     1.756E01,1.794E01,1.836E01,1.880E01,1.926E01,1.973E01,0.0,13*XSE  470
     *     0.0/                                                         XSE  480
      DATA X08 /2.000E00,2.002E00,2.004E00,2.007E00,2.010E00,2.014E00,  XSE  490
     *     2.041E00,2.047E00,2.052E00,2.057E00,2.064E00,2.070E00,       XSE  500
     *     2.091E00,2.098E00,2.107E00,2.115E00,2.125E00,2.134E00,0.0,13*XSE  510
     *     0.0/                                                         XSE  520
      DATA X09 /2.267E01,2.288E01,2.309E01,2.334E01,2.364E01,2.399E01,  XSE  530
     *     2.622E01,2.672E01,2.721E01,2.774E01,2.832E01,2.893E01,       XSE  540
     *     3.017E01,3.086E01,3.162E01,3.243E01,3.327E01,3.413E01,0.0,13*XSE  550
     *     0.0/                                                         XSE  560
      DATA X10 /4.034E01,4.072E01,4.111E01,4.156E01,4.211E01,4.274E01,  XSE  570
     *     4.666E01,4.756E01,4.845E01,4.941E01,5.046E01,5.158E01,       XSE  580
     *     5.383E01,5.507E01,5.645E01,5.791E01,5.943E01,6.100E01,0.0,13*XSE  590
     *     0.0/                                                         XSE  600
      DATA X11 /4.012E01,4.020E01,4.026E01,4.035E01,4.048E01,4.065E01,  XSE  610
     *     4.337E01,4.362E01,4.386E01,4.413E01,4.443E01,4.477E01,       XSE  620
     *     4.705E01,4.756E01,4.804E01,4.853E01,4.904E01,4.956E01,0.0,13*XSE  630
     *     0.0/                                                         XSE  640
      DATA X12 /1.460E01,1.473E01,1.486E01,1.502E01,1.521E01,1.542E01,  XSE  650
     *     1.683E01,1.714E01,1.745E01,1.778E01,1.814E01,1.853E01,       XSE  660
     *     1.938E01,1.982E01,2.030E01,2.081E01,2.134E01,2.189E01,0.0,13*XSE  670
     *     0.0/                                                         XSE  680
      DATA X13 /4.828E01,4.866E01,4.904E01,4.949E01,5.003E01,5.066E01,  XSE  690
     *     5.481E01,5.570E01,5.660E01,5.756E01,5.861E01,5.973E01,       XSE  700
     *     6.220E01,6.347E01,6.487E01,6.635E01,6.788E01,6.947E01,0.0,13*XSE  710
     *     0.0/                                                         XSE  720
      DATA X14 /2.502E01,2.501E01,2.499E01,2.498E01,2.497E01,2.497E01,  XSE  730
     *     2.515E01,2.517E01,2.517E01,2.517E01,2.518E01,2.518E01,       XSE  740
     *     2.484E01,2.491E01,2.495E01,2.499E01,2.501E01,2.503E01,0.0,13*XSE  750
     *     0.0/                                                         XSE  760
      DATA X15 /5.133E01,5.164E01,5.194E01,5.231E01,5.276E01,5.329E01,  XSE  770
     *     5.595E01,5.668E01,5.742E01,5.820E01,5.906E01,5.999E01,       XSE  780
     *     6.096E01,6.204E01,6.318E01,6.439E01,6.565E01,6.694E01,0.0,13*XSE  790
     *     0.0/                                                         XSE  800
      DATA X16 /6.344E00,6.382E00,6.420E00,6.465E00,6.521E00,6.586E00,  XSE  810
     *     6.915E00,7.006E00,7.096E00,7.194E00,7.300E00,7.415E00,       XSE  820
     *     7.535E00,7.667E00,7.809E00,7.959E00,8.114E00,8.274E00,0.0,13*XSE  830
     *     0.0/                                                         XSE  840
      DATA X17 /1.970E00,1.970E00,1.969E00,1.969E00,1.970E00,1.972E00,  XSE  850
     *     2.038E00,2.042E00,2.046E00,2.050E00,2.056E00,2.062E00,       XSE  860
     *     2.068E00,2.082E00,2.094E00,2.104E00,2.116E00,2.127E00,0.0,13*XSE  870
     *     0.0/                                                         XSE  880
      DATA X18 /3.742E01,3.742E01,3.741E01,3.741E01,3.743E01,3.748E01,  XSE  890
     *     3.872E01,3.880E01,3.887E01,3.895E01,3.906E01,3.918E01,       XSE  900
     *     3.930E01,3.956E01,3.978E01,3.999E01,4.020E01,4.042E01,0.0,13*XSE  910
     *     0.0/                                                         XSE  920
      DATA X19 /4.933E00,4.935E00,4.935E00,4.936E00,4.939E00,4.944E00,  XSE  930
     *     5.083E00,5.092E00,5.099E00,5.106E00,5.115E00,5.124E00,       XSE  940
     *     5.167E00,5.187E00,5.204E00,5.220E00,5.236E00,5.251E00,0.0,13*XSE  950
     *     0.0/                                                         XSE  960
      DATA X20 /1.357E01,1.357E01,1.357E01,1.358E01,1.359E01,1.360E01,  XSE  970
     *     1.391E01,1.394E01,1.396E01,1.398E01,1.400E01,1.402E01,       XSE  980
     *     1.413E01,1.418E01,1.422E01,1.426E01,1.430E01,1.435E01,0.0,13*XSE  990
     *     0.0/                                                         XSE 1000
      NPT = MPT                                                         XSE 1010
      NXSEC = MXSEC                                                     XSE 1020
      DO 10 I=1,30                                                      XSE 1030
         F(I,01) = X01(I)                                               XSE 1040
         F(I,02) = X02(I)                                               XSE 1050
         F(I,03) = X03(I)                                               XSE 1060
         F(I,04) = X04(I)                                               XSE 1070
         F(I,05) = X05(I)                                               XSE 1080
         F(I,06) = X06(I)                                               XSE 1090
         F(I,07) = X07(I)                                               XSE 1100
         F(I,08) = X08(I)                                               XSE 1110
         F(I,09) = X09(I)                                               XSE 1120
         F(I,10) = X10(I)                                               XSE 1130
         F(I,11) = X11(I)                                               XSE 1140
         F(I,12) = X12(I)                                               XSE 1150
         F(I,13) = X13(I)                                               XSE 1160
         F(I,14) = X14(I)                                               XSE 1170
         F(I,15) = X15(I)                                               XSE 1180
         F(I,16) = X16(I)                                               XSE 1190
         F(I,17) = X17(I)                                               XSE 1200
         F(I,18) = X18(I)                                               XSE 1210
         F(I,19) = X19(I)                                               XSE 1220
         F(I,20) = X20(I)                                               XSE 1230
         CBUP(I) = DBUP(I)                                              XSE 1240
   10 CONTINUE                                                          XSE 1250
      DO 30 I=1,20                                                      XSE 1260
         DO 20 J=1,7                                                    XSE 1270
            NFUD(J,I) = MFUD(J,I)                                       XSE 1280
   20    CONTINUE                                                       XSE 1290
   30 CONTINUE                                                          XSE 1300
      RETURN                                                            XSE 1310
      END                                                               XSE 1320
      SUBROUTINE XSEC04(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    *BWR-U* ;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,3.479E-01,7.019E-01,1.060E00,1.422E00,1.787E00,    XSE  120
     *     1.787E00,2.132E00,2.477E00,2.821E00,3.166E00,3.510E00,       XSE  130
     *     3.510E00,3.836E00,4.160E00,4.483E00,4.805E00,5.127E00,       XSE  140
     *     5.127E00,5.432E00,5.735E00,6.037E00,6.339E00,6.640E00,0.0,7* XSE  150
     *     0.0/                                                         XSE  160
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  170
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  180
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  190
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  200
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  210
      DATA X01 /2.092E01,2.091E01,2.090E01,2.089E01,2.089E01,2.090E01,  XSE  220
     *     2.083E01,2.082E01,2.081E01,2.080E01,2.080E01,2.081E01,       XSE  230
     *     2.083E01,2.081E01,2.080E01,2.079E01,2.079E01,2.079E01,       XSE  240
     *     2.088E01,2.086E01,2.085E01,2.084E01,2.084E01,2.084E01,0.0,7* XSE  250
     *     0.0/                                                         XSE  260
      DATA X02 /1.099E01,1.094E01,1.090E01,1.087E01,1.085E01,1.083E01,  XSE  270
     *     1.106E01,1.101E01,1.097E01,1.094E01,1.091E01,1.090E01,       XSE  280
     *     1.118E01,1.113E01,1.109E01,1.106E01,1.103E01,1.102E01,       XSE  290
     *     1.133E01,1.128E01,1.124E01,1.121E01,1.118E01,1.117E01,0.0,7* XSE  300
     *     0.0/                                                         XSE  310
      DATA X03 /4.907E01,4.877E01,4.852E01,4.833E01,4.818E01,4.807E01,  XSE  320
     *     4.945E01,4.914E01,4.889E01,4.870E01,4.855E01,4.844E01,       XSE  330
     *     5.006E01,4.974E01,4.949E01,4.929E01,4.914E01,4.904E01,       XSE  340
     *     5.083E01,5.051E01,5.026E01,5.006E01,4.990E01,4.979E01,0.0,7* XSE  350
     *     0.0/                                                         XSE  360
      DATA X04 /8.334E00,8.352E00,8.368E00,8.382E00,8.395E00,8.406E00,  XSE  370
     *     8.252E00,8.270E00,8.286E00,8.300E00,8.313E00,8.324E00,       XSE  380
     *     8.198E00,8.216E00,8.232E00,8.246E00,8.258E00,8.269E00,       XSE  390
     *     8.164E00,8.181E00,8.197E00,8.211E00,8.223E00,8.234E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /9.005E-01,9.005E-01,9.006E-01,9.008E-01,9.011E-01,      XSE  420
     *     9.014E-01,9.033E-01,9.033E-01,9.034E-01,9.036E-01,9.039E-01, XSE  430
     *     9.042E-01,9.076E-01,9.076E-01,9.077E-01,9.079E-01,9.081E-01, XSE  440
     *     9.084E-01,9.132E-01,9.132E-01,9.133E-01,9.135E-01,9.137E-01, XSE  450
     *     9.140E-01,0.0,7*0.0/                                         XSE  460
      DATA X06 /3.395E01,3.382E01,3.370E01,3.360E01,3.352E01,3.345E01,  XSE  470
     *     3.437E01,3.423E01,3.411E01,3.401E01,3.392E01,3.386E01,       XSE  480
     *     3.451E01,3.437E01,3.425E01,3.415E01,3.406E01,3.400E01,       XSE  490
     *     3.453E01,3.439E01,3.427E01,3.416E01,3.408E01,3.401E01,0.0,7* XSE  500
     *     0.0/                                                         XSE  510
      DATA X07 /3.592E01,3.569E01,3.550E01,3.535E01,3.524E01,3.516E01,  XSE  520
     *     3.647E01,3.623E01,3.604E01,3.589E01,3.578E01,3.570E01,       XSE  530
     *     3.713E01,3.689E01,3.669E01,3.654E01,3.642E01,3.634E01,       XSE  540
     *     3.790E01,3.765E01,3.745E01,3.729E01,3.718E01,3.709E01,0.0,7* XSE  550
     *     0.0/                                                         XSE  560
      DATA X08 /2.437E00,2.432E00,2.428E00,2.426E00,2.424E00,2.422E00,  XSE  570
     *     2.455E00,2.450E00,2.446E00,2.443E00,2.441E00,2.440E00,       XSE  580
     *     2.476E00,2.471E00,2.467E00,2.464E00,2.462E00,2.461E00,       XSE  590
     *     2.500E00,2.495E00,2.491E00,2.488E00,2.486E00,2.484E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /6.910E01,6.862E01,6.824E01,6.794E01,6.772E01,6.756E01,  XSE  620
     *     6.622E01,6.576E01,6.539E01,6.511E01,6.489E01,6.474E01,       XSE  630
     *     6.490E01,6.446E01,6.410E01,6.382E01,6.361E01,6.346E01,       XSE  640
     *     6.419E01,6.375E01,6.340E01,6.312E01,6.291E01,6.277E01,0.0,7* XSE  650
     *     0.0/                                                         XSE  660
      DATA X10 /1.219E02,1.211E02,1.204E02,1.199E02,1.195E02,1.192E02,  XSE  670
     *     1.182E02,1.174E02,1.167E02,1.162E02,1.158E02,1.155E02,       XSE  680
     *     1.168E02,1.159E02,1.153E02,1.148E02,1.144E02,1.141E02,       XSE  690
     *     1.162E02,1.154E02,1.148E02,1.142E02,1.139E02,1.136E02,0.0,7* XSE  700
     *     0.0/                                                         XSE  710
      DATA X11 /2.060E02,2.051E02,2.039E02,2.026E02,2.014E02,2.001E02,  XSE  720
     *     1.577E02,1.570E02,1.561E02,1.551E02,1.542E02,1.533E02,       XSE  730
     *     1.309E02,1.303E02,1.295E02,1.288E02,1.280E02,1.273E02,       XSE  740
     *     1.150E02,1.144E02,1.138E02,1.131E02,1.125E02,1.118E02,0.0,7* XSE  750
     *     0.0/                                                         XSE  760
      DATA X12 /4.290E01,4.261E01,4.238E01,4.219E01,4.206E01,4.196E01,  XSE  770
     *     4.219E01,4.190E01,4.168E01,4.150E01,4.136E01,4.127E01,       XSE  780
     *     4.208E01,4.179E01,4.156E01,4.138E01,4.125E01,4.116E01,       XSE  790
     *     4.224E01,4.195E01,4.172E01,4.154E01,4.141E01,4.131E01,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.291E02,1.283E02,1.276E02,1.271E02,1.267E02,1.265E02,  XSE  820
     *     1.277E02,1.268E02,1.262E02,1.257E02,1.253E02,1.250E02,       XSE  830
     *     1.277E02,1.269E02,1.262E02,1.257E02,1.253E02,1.250E02,       XSE  840
     *     1.284E02,1.276E02,1.269E02,1.264E02,1.260E02,1.258E02,0.0,7* XSE  850
     *     0.0/                                                         XSE  860
      DATA X14 /3.322E01,3.329E01,3.336E01,3.342E01,3.347E01,3.352E01,  XSE  870
     *     3.321E01,3.329E01,3.335E01,3.341E01,3.347E01,3.351E01,       XSE  880
     *     3.306E01,3.313E01,3.320E01,3.326E01,3.331E01,3.336E01,       XSE  890
     *     3.269E01,3.276E01,3.283E01,3.289E01,3.294E01,3.298E01,0.0,7* XSE  900
     *     0.0/                                                         XSE  910
      DATA X15 /1.158E02,1.151E02,1.146E02,1.141E02,1.137E02,1.134E02,  XSE  920
     *     1.130E02,1.123E02,1.117E02,1.113E02,1.109E02,1.106E02,       XSE  930
     *     1.114E02,1.108E02,1.102E02,1.098E02,1.094E02,1.091E02,       XSE  940
     *     1.103E02,1.097E02,1.091E02,1.087E02,1.083E02,1.080E02,0.0,7* XSE  950
     *     0.0/                                                         XSE  960
      DATA X16 /1.431E01,1.423E01,1.416E01,1.410E01,1.405E01,1.402E01,  XSE  970
     *     1.396E01,1.388E01,1.381E01,1.375E01,1.371E01,1.367E01,       XSE  980
     *     1.377E01,1.369E01,1.362E01,1.357E01,1.352E01,1.349E01,       XSE  990
     *     1.363E01,1.355E01,1.349E01,1.343E01,1.339E01,1.335E01,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /2.581E00,2.572E00,2.561E00,2.549E00,2.537E00,2.526E00,  XSE 1020
     *     2.621E00,2.612E00,2.600E00,2.589E00,2.577E00,2.565E00,       XSE 1030
     *     2.641E00,2.631E00,2.620E00,2.608E00,2.596E00,2.584E00,       XSE 1040
     *     2.652E00,2.642E00,2.631E00,2.619E00,2.607E00,2.595E00,0.0,7* XSE 1050
     *     0.0/                                                         XSE 1060
      DATA X18 /4.904E01,4.886E01,4.866E01,4.843E01,4.821E01,4.800E01,  XSE 1070
     *     4.980E01,4.962E01,4.941E01,4.918E01,4.895E01,4.873E01,       XSE 1080
     *     5.018E01,5.000E01,4.978E01,4.956E01,4.933E01,4.910E01,       XSE 1090
     *     5.039E01,5.020E01,4.999E01,4.976E01,4.953E01,4.930E01,0.0,7* XSE 1100
     *     0.0/                                                         XSE 1110
      DATA X19 /5.733E00,5.735E00,5.737E00,5.740E00,5.743E00,5.746E00,  XSE 1120
     *     5.744E00,5.746E00,5.748E00,5.751E00,5.753E00,5.757E00,       XSE 1130
     *     5.765E00,5.766E00,5.768E00,5.771E00,5.774E00,5.777E00,       XSE 1140
     *     5.792E00,5.794E00,5.796E00,5.798E00,5.801E00,5.804E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /1.420E01,1.420E01,1.420E01,1.421E01,1.421E01,1.422E01,  XSE 1170
     *     1.424E01,1.424E01,1.424E01,1.424E01,1.425E01,1.425E01,       XSE 1180
     *     1.430E01,1.430E01,1.430E01,1.430E01,1.430E01,1.431E01,       XSE 1190
     *     1.437E01,1.437E01,1.437E01,1.437E01,1.437E01,1.438E01,0.0,7* XSE 1200
     *     0.0/                                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC05(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    BWR-PUU ;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,4.067E-01,8.123E-01,1.214E00,1.610E00,2.001E00,    XSE  120
     *     2.001E00,2.343E00,2.683E00,3.023E00,3.362E00,3.700E00,       XSE  130
     *     3.700E00,3.972E00,4.248E00,4.527E00,4.809E00,5.096E00,       XSE  140
     *     5.096E00,5.330E00,5.568E00,5.812E00,6.062E00,6.319E00,0.0,7* XSE  150
     *     0.0/                                                         XSE  160
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  170
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  180
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  190
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  200
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  210
      DATA X01 /2.067E01,2.058E01,2.052E01,2.049E01,2.047E01,2.047E01,  XSE  220
     *     2.023E01,2.023E01,2.024E01,2.024E01,2.025E01,2.027E01,       XSE  230
     *     1.998E01,2.000E01,2.001E01,2.002E01,2.003E01,2.004E01,       XSE  240
     *     1.997E01,1.998E01,1.999E01,2.000E01,2.001E01,2.002E01,0.0,7* XSE  250
     *     0.0/                                                         XSE  260
      DATA X02 /1.070E01,1.054E01,1.043E01,1.035E01,1.029E01,1.025E01,  XSE  270
     *     1.013E01,1.010E01,1.008E01,1.007E01,1.007E01,1.007E01,       XSE  280
     *     9.855E00,9.856E00,9.858E00,9.863E00,9.870E00,9.881E00,       XSE  290
     *     9.867E00,9.876E00,9.882E00,9.891E00,9.901E00,9.913E00,0.0,7* XSE  300
     *     0.0/                                                         XSE  310
      DATA X03 /4.741E01,4.648E01,4.580E01,4.529E01,4.493E01,4.467E01,  XSE  320
     *     4.392E01,4.373E01,4.360E01,4.352E01,4.347E01,4.347E01,       XSE  330
     *     4.211E01,4.210E01,4.211E01,4.213E01,4.217E01,4.222E01,       XSE  340
     *     4.196E01,4.200E01,4.204E01,4.209E01,4.215E01,4.222E01,0.0,7* XSE  350
     *     0.0/                                                         XSE  360
      DATA X04 /8.341E00,8.370E00,8.398E00,8.423E00,8.445E00,8.465E00,  XSE  370
     *     8.400E00,8.416E00,8.430E00,8.441E00,8.450E00,8.457E00,       XSE  380
     *     8.425E00,8.431E00,8.436E00,8.440E00,8.443E00,8.445E00,       XSE  390
     *     8.438E00,8.441E00,8.442E00,8.443E00,8.443E00,8.442E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /8.996E-01,8.982E-01,8.975E-01,8.973E-01,8.973E-01,      XSE  420
     *     8.975E-01,8.979E-01,8.982E-01,8.986E-01,8.989E-01,8.993E-01, XSE  430
     *     8.997E-01,8.993E-01,8.997E-01,9.000E-01,9.003E-01,9.006E-01, XSE  440
     *     9.009E-01,9.053E-01,9.056E-01,9.058E-01,9.060E-01,9.062E-01, XSE  450
     *     9.064E-01,0.0,7*0.0/                                         XSE  460
      DATA X06 /3.398E01,3.358E01,3.327E01,3.303E01,3.284E01,3.270E01,  XSE  470
     *     3.255E01,3.245E01,3.237E01,3.231E01,3.227E01,3.225E01,       XSE  480
     *     3.162E01,3.160E01,3.159E01,3.158E01,3.158E01,3.159E01,       XSE  490
     *     3.138E01,3.139E01,3.139E01,3.140E01,3.141E01,3.142E01,0.0,7* XSE  500
     *     0.0/                                                         XSE  510
      DATA X07 /3.451E01,3.378E01,3.325E01,3.287E01,3.259E01,3.239E01,  XSE  520
     *     3.202E01,3.188E01,3.178E01,3.172E01,3.169E01,3.169E01,       XSE  530
     *     3.074E01,3.074E01,3.074E01,3.076E01,3.079E01,3.084E01,       XSE  540
     *     3.071E01,3.075E01,3.078E01,3.082E01,3.087E01,3.093E01,0.0,7* XSE  550
     *     0.0/                                                         XSE  560
      DATA X08 /2.405E00,2.390E00,2.380E00,2.372E00,2.367E00,2.363E00,  XSE  570
     *     2.360E00,2.357E00,2.356E00,2.355E00,2.354E00,2.355E00,       XSE  580
     *     2.339E00,2.339E00,2.339E00,2.340E00,2.341E00,2.342E00,       XSE  590
     *     2.343E00,2.344E00,2.345E00,2.346E00,2.347E00,2.349E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /6.731E01,6.583E01,6.474E01,6.395E01,6.338E01,6.297E01,  XSE  620
     *     5.874E01,5.847E01,5.827E01,5.815E01,5.809E01,5.810E01,       XSE  630
     *     5.435E01,5.435E01,5.436E01,5.439E01,5.445E01,5.454E01,       XSE  640
     *     5.289E01,5.295E01,5.301E01,5.308E01,5.316E01,5.327E01,0.0,7* XSE  650
     *     0.0/                                                         XSE  660
      DATA X10 /1.185E02,1.158E02,1.139E02,1.125E02,1.115E02,1.108E02,  XSE  670
     *     1.044E02,1.039E02,1.036E02,1.033E02,1.033E02,1.033E02,       XSE  680
     *     9.720E01,9.719E01,9.721E01,9.726E01,9.737E01,9.753E01,       XSE  690
     *     9.505E01,9.516E01,9.526E01,9.539E01,9.554E01,9.573E01,0.0,7* XSE  700
     *     0.0/                                                         XSE  710
      DATA X11 /2.053E02,2.049E02,2.036E02,2.018E02,1.999E02,1.981E02,  XSE  720
     *     1.534E02,1.525E02,1.515E02,1.506E02,1.497E02,1.488E02,       XSE  730
     *     1.237E02,1.233E02,1.228E02,1.223E02,1.218E02,1.214E02,       XSE  740
     *     1.071E02,1.068E02,1.065E02,1.062E02,1.060E02,1.057E02,0.0,7* XSE  750
     *     0.0/                                                         XSE  760
      DATA X12 /4.142E01,4.052E01,3.986E01,3.938E01,3.903E01,3.878E01,  XSE  770
     *     3.713E01,3.696E01,3.684E01,3.677E01,3.673E01,3.673E01,       XSE  780
     *     3.491E01,3.491E01,3.492E01,3.494E01,3.498E01,3.503E01,       XSE  790
     *     3.438E01,3.442E01,3.446E01,3.451E01,3.456E01,3.463E01,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.247E02,1.221E02,1.202E02,1.189E02,1.179E02,1.172E02,  XSE  820
     *     1.130E02,1.125E02,1.122E02,1.120E02,1.119E02,1.119E02,       XSE  830
     *     1.069E02,1.069E02,1.069E02,1.070E02,1.071E02,1.073E02,       XSE  840
     *     1.056E02,1.057E02,1.058E02,1.059E02,1.061E02,1.063E02,0.0,7* XSE  850
     *     0.0/                                                         XSE  860
      DATA X14 /3.276E01,3.287E01,3.298E01,3.308E01,3.318E01,3.325E01,  XSE  870
     *     3.328E01,3.335E01,3.341E01,3.345E01,3.349E01,3.352E01,       XSE  880
     *     3.342E01,3.345E01,3.347E01,3.349E01,3.350E01,3.351E01,       XSE  890
     *     3.326E01,3.328E01,3.328E01,3.328E01,3.328E01,3.328E01,0.0,7* XSE  900
     *     0.0/                                                         XSE  910
      DATA X15 /1.108E02,1.089E02,1.074E02,1.063E02,1.054E02,1.048E02,  XSE  920
     *     1.001E02,9.969E01,9.936E01,9.911E01,9.897E01,9.890E01,       XSE  930
     *     9.428E01,9.422E01,9.418E01,9.417E01,9.419E01,9.425E01,       XSE  940
     *     9.222E01,9.226E01,9.229E01,9.233E01,9.240E01,9.250E01,0.0,7* XSE  950
     *     0.0/                                                         XSE  960
      DATA X16 /1.370E01,1.346E01,1.328E01,1.314E01,1.303E01,1.295E01,  XSE  970
     *     1.237E01,1.232E01,1.228E01,1.225E01,1.223E01,1.222E01,       XSE  980
     *     1.165E01,1.165E01,1.164E01,1.164E01,1.164E01,1.165E01,       XSE  990
     *     1.140E01,1.140E01,1.141E01,1.141E01,1.142E01,1.143E01,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /2.641E00,2.634E00,2.620E00,2.603E00,2.584E00,2.566E00,  XSE 1020
     *     2.610E00,2.598E00,2.586E00,2.574E00,2.562E00,2.551E00,       XSE 1030
     *     2.570E00,2.562E00,2.554E00,2.546E00,2.538E00,2.531E00,       XSE 1040
     *     2.554E00,2.548E00,2.543E00,2.537E00,2.531E00,2.526E00,0.0,7* XSE 1050
     *     0.0/                                                         XSE 1060
      DATA X18 /5.018E01,5.005E01,4.978E01,4.945E01,4.910E01,4.876E01,  XSE 1070
     *     4.959E01,4.936E01,4.913E01,4.890E01,4.868E01,4.848E01,       XSE 1080
     *     4.883E01,4.868E01,4.852E01,4.837E01,4.823E01,4.809E01,       XSE 1090
     *     4.852E01,4.842E01,4.831E01,4.820E01,4.810E01,4.800E01,0.0,7* XSE 1100
     *     0.0/                                                         XSE 1110
      DATA X19 /5.700E00,5.692E00,5.689E00,5.690E00,5.692E00,5.695E00,  XSE 1120
     *     5.688E00,5.692E00,5.695E00,5.699E00,5.703E00,5.707E00,       XSE 1130
     *     5.693E00,5.696E00,5.699E00,5.701E00,5.704E00,5.706E00,       XSE 1140
     *     5.724E00,5.727E00,5.729E00,5.730E00,5.732E00,5.733E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /1.445E01,1.443E01,1.441E01,1.441E01,1.441E01,1.441E01,  XSE 1170
     *     1.440E01,1.441E01,1.441E01,1.442E01,1.443E01,1.443E01,       XSE 1180
     *     1.441E01,1.441E01,1.442E01,1.442E01,1.443E01,1.443E01,       XSE 1190
     *     1.448E01,1.449E01,1.449E01,1.449E01,1.450E01,1.450E01,0.0,7* XSE 1200
     *     0.0/                                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC06(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    BWR-PUPU;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,3.166E-01,6.384E-01,9.664E-01,1.301E00,1.642E00,   XSE  120
     *     1.642E00,1.989E00,2.336E00,2.687E00,3.040E00,3.397E00,       XSE  130
     *     3.397E00,3.794E00,4.187E00,4.576E00,4.962E00,5.345E00,       XSE  140
     *     5.345E00,5.774E00,6.194E00,6.604E00,7.005E00,7.396E00,0.0,7* XSE  150
     *     0.0/                                                         XSE  160
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  170
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  180
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  190
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  200
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  210
      DATA X01 /1.754E01,1.759E01,1.762E01,1.766E01,1.770E01,1.774E01,  XSE  220
     *     1.797E01,1.802E01,1.807E01,1.812E01,1.816E01,1.821E01,       XSE  230
     *     1.859E01,1.866E01,1.872E01,1.878E01,1.884E01,1.890E01,       XSE  240
     *     1.906E01,1.911E01,1.917E01,1.923E01,1.929E01,1.935E01,0.0,7* XSE  250
     *     0.0/                                                         XSE  260
      DATA X02 /6.881E00,6.924E00,6.965E00,7.008E00,7.052E00,7.097E00,  XSE  270
     *     7.328E00,7.381E00,7.431E00,7.483E00,7.535E00,7.589E00,       XSE  280
     *     7.996E00,8.058E00,8.120E00,8.181E00,8.244E00,8.308E00,       XSE  290
     *     8.633E00,8.669E00,8.713E00,8.761E00,8.813E00,8.867E00,0.0,7* XSE  300
     *     0.0/                                                         XSE  310
      DATA X03 /2.543E01,2.567E01,2.591E01,2.616E01,2.642E01,2.668E01,  XSE  320
     *     2.781E01,2.810E01,2.840E01,2.869E01,2.900E01,2.931E01,       XSE  330
     *     3.150E01,3.185E01,3.219E01,3.255E01,3.290E01,3.327E01,       XSE  340
     *     3.516E01,3.534E01,3.557E01,3.583E01,3.612E01,3.642E01,0.0,7* XSE  350
     *     0.0/                                                         XSE  360
      DATA X04 /8.452E00,8.454E00,8.454E00,8.453E00,8.452E00,8.450E00,  XSE  370
     *     8.494E00,8.497E00,8.498E00,8.499E00,8.498E00,8.497E00,       XSE  380
     *     8.526E00,8.533E00,8.538E00,8.541E00,8.543E00,8.544E00,       XSE  390
     *     8.478E00,8.493E00,8.507E00,8.518E00,8.527E00,8.534E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /8.296E-01,8.305E-01,8.313E-01,8.320E-01,8.328E-01,      XSE  420
     *     8.336E-01,8.435E-01,8.447E-01,8.457E-01,8.467E-01,8.476E-01, XSE  430
     *     8.486E-01,8.599E-01,8.615E-01,8.629E-01,8.643E-01,8.655E-01, XSE  440
     *     8.668E-01,8.726E-01,8.742E-01,8.758E-01,8.773E-01,8.787E-01, XSE  450
     *     8.801E-01,0.0,7*0.0/                                         XSE  460
      DATA X06 /2.521E01,2.535E01,2.548E01,2.561E01,2.575E01,2.590E01,  XSE  470
     *     2.626E01,2.643E01,2.658E01,2.674E01,2.691E01,2.708E01,       XSE  480
     *     2.798E01,2.817E01,2.836E01,2.855E01,2.874E01,2.893E01,       XSE  490
     *     2.944E01,2.956E01,2.970E01,2.984E01,3.000E01,3.016E01,0.0,7* XSE  500
     *     0.0/                                                         XSE  510
      DATA X07 /1.752E01,1.771E01,1.790E01,1.810E01,1.830E01,1.851E01,  XSE  520
     *     1.937E01,1.960E01,1.983E01,2.006E01,2.030E01,2.055E01,       XSE  530
     *     2.225E01,2.253E01,2.280E01,2.308E01,2.336E01,2.364E01,       XSE  540
     *     2.514E01,2.529E01,2.547E01,2.567E01,2.589E01,2.613E01,0.0,7* XSE  550
     *     0.0/                                                         XSE  560
      DATA X08 /2.044E00,2.048E00,2.052E00,2.056E00,2.060E00,2.064E00,  XSE  570
     *     2.084E00,2.089E00,2.093E00,2.098E00,2.103E00,2.108E00,       XSE  580
     *     2.145E00,2.151E00,2.156E00,2.162E00,2.167E00,2.173E00,       XSE  590
     *     2.208E00,2.211E00,2.215E00,2.219E00,2.224E00,2.229E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /2.894E01,2.928E01,2.961E01,2.995E01,3.031E01,3.068E01,  XSE  620
     *     3.251E01,3.292E01,3.333E01,3.374E01,3.417E01,3.461E01,       XSE  630
     *     3.780E01,3.829E01,3.877E01,3.926E01,3.976E01,4.028E01,       XSE  640
     *     4.299E01,4.324E01,4.356E01,4.392E01,4.431E01,4.473E01,0.0,7* XSE  650
     *     0.0/                                                         XSE  660
      DATA X10 /5.209E01,5.271E01,5.332E01,5.395E01,5.460E01,5.527E01,  XSE  670
     *     5.843E01,5.918E01,5.992E01,6.068E01,6.146E01,6.225E01,       XSE  680
     *     6.798E01,6.885E01,6.974E01,7.063E01,7.155E01,7.248E01,       XSE  690
     *     7.738E01,7.783E01,7.841E01,7.906E01,7.978E01,8.054E01,0.0,7* XSE  700
     *     0.0/                                                         XSE  710
      DATA X11 /4.418E01,4.432E01,4.445E01,4.459E01,4.474E01,4.490E01,  XSE  720
     *     4.651E01,4.669E01,4.687E01,4.705E01,4.724E01,4.744E01,       XSE  730
     *     5.055E01,5.079E01,5.102E01,5.125E01,5.149E01,5.174E01,       XSE  740
     *     5.465E01,5.489E01,5.511E01,5.532E01,5.552E01,5.574E01,0.0,7* XSE  750
     *     0.0/                                                         XSE  760
      DATA X12 /1.899E01,1.921E01,1.943E01,1.965E01,1.988E01,2.012E01,  XSE  770
     *     2.119E01,2.145E01,2.172E01,2.198E01,2.226E01,2.254E01,       XSE  780
     *     2.455E01,2.485E01,2.517E01,2.548E01,2.581E01,2.613E01,       XSE  790
     *     2.788E01,2.804E01,2.824E01,2.848E01,2.873E01,2.900E01,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /6.136E01,6.200E01,6.262E01,6.327E01,6.393E01,6.463E01,  XSE  820
     *     6.770E01,6.847E01,6.923E01,7.000E01,7.080E01,7.161E01,       XSE  830
     *     7.741E01,7.831E01,7.922E01,8.014E01,8.107E01,8.203E01,       XSE  840
     *     8.699E01,8.746E01,8.806E01,8.875E01,8.949E01,9.027E01,0.0,7* XSE  850
     *     0.0/                                                         XSE  860
      DATA X14 /2.630E01,2.631E01,2.632E01,2.632E01,2.632E01,2.633E01,  XSE  870
     *     2.626E01,2.628E01,2.629E01,2.630E01,2.630E01,2.631E01,       XSE  880
     *     2.620E01,2.624E01,2.626E01,2.628E01,2.630E01,2.631E01,       XSE  890
     *     2.581E01,2.587E01,2.592E01,2.596E01,2.600E01,2.603E01,0.0,7* XSE  900
     *     0.0/                                                         XSE  910
      DATA X15 /5.555E01,5.600E01,5.644E01,5.689E01,5.736E01,5.786E01,  XSE  920
     *     5.971E01,6.026E01,6.079E01,6.134E01,6.190E01,6.248E01,       XSE  930
     *     6.625E01,6.689E01,6.754E01,6.818E01,6.884E01,6.952E01,       XSE  940
     *     7.235E01,7.272E01,7.317E01,7.365E01,7.418E01,7.473E01,0.0,7* XSE  950
     *     0.0/                                                         XSE  960
      DATA X16 /6.865E00,6.921E00,6.975E00,7.032E00,7.090E00,7.151E00,  XSE  970
     *     7.380E00,7.447E00,7.514E00,7.581E00,7.650E00,7.722E00,       XSE  980
     *     8.189E00,8.268E00,8.347E00,8.427E00,8.509E00,8.592E00,       XSE  990
     *     8.942E00,8.988E00,9.043E00,9.103E00,9.168E00,9.237E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /2.189E00,2.190E00,2.191E00,2.192E00,2.193E00,2.195E00,  XSE 1020
     *     2.229E00,2.231E00,2.233E00,2.235E00,2.238E00,2.240E00,       XSE 1030
     *     2.297E00,2.301E00,2.305E00,2.309E00,2.313E00,2.317E00,       XSE 1040
     *     2.333E00,2.341E00,2.347E00,2.352E00,2.357E00,2.361E00,0.0,7* XSE 1050
     *     0.0/                                                         XSE 1060
      DATA X18 /4.160E01,4.162E01,4.163E01,4.165E01,4.168E01,4.171E01,  XSE 1070
     *     4.234E01,4.239E01,4.243E01,4.247E01,4.251E01,4.257E01,       XSE 1080
     *     4.364E01,4.372E01,4.380E01,4.387E01,4.394E01,4.402E01,       XSE 1090
     *     4.433E01,4.448E01,4.460E01,4.470E01,4.478E01,4.486E01,0.0,7* XSE 1100
     *     0.0/                                                         XSE 1110
      DATA X19 /5.134E00,5.141E00,5.147E00,5.152E00,5.158E00,5.164E00,  XSE 1120
     *     5.238E00,5.246E00,5.254E00,5.261E00,5.268E00,5.275E00,       XSE 1130
     *     5.358E00,5.370E00,5.381E00,5.391E00,5.401E00,5.410E00,       XSE 1140
     *     5.435E00,5.447E00,5.459E00,5.471E00,5.482E00,5.493E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /1.379E01,1.381E01,1.383E01,1.384E01,1.386E01,1.387E01,  XSE 1170
     *     1.404E01,1.406E01,1.408E01,1.410E01,1.412E01,1.414E01,       XSE 1180
     *     1.433E01,1.436E01,1.439E01,1.442E01,1.444E01,1.447E01,       XSE 1190
     *     1.454E01,1.457E01,1.460E01,1.463E01,1.465E01,1.468E01,0.0,7* XSE 1200
     *     0.0/                                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC07(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PTD3   ;XSDRNPM/CITATION,XSEC FROM CITATION                        XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION MFUD(7,21), DBUP(30)                                    XSE   60
      DIMENSION X01(30), X02(30), X03(30), X04(30), X05(30), X06(30),   XSE   70
     *     X07(30), X08(30), X09(30), X10(30), X11(30), X12(30),        XSE   80
     *     X13(30), X14(30), X15(30), X16(30), X17(30), X18(30),        XSE   90
     *     X19(30), X20(30)                                             XSE  100
      DATA MXSEC, MPT /19,18/                                           XSE  110
      DATA DBUP /0.0,5.566E-01,1.113E00,1.671E00,2.229E00,2.787E00,     XSE  120
     *     2.787E00,3.305E00,3.825E00,4.346E00,4.867E00,5.390E00,       XSE  130
     *     5.390E00,5.876E00,6.364E00,6.854E00,7.347E00,7.842E00,0.0,11*XSE  140
     *     0.0/                                                         XSE  150
      DATA MFUD /902320,1,5*0,912330,1,5*0,922330,1,5*0,922330,4,5*0,   XSE  160
     *     922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,922380,1,XSE  170
     *     5*0,932370,1,5*0,942380,1,5*0,942390,1,5*0,942390,4,5*0,     XSE  180
     *     942400,1,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,952410,1,XSE  190
     *     5*0,952410,2,5*0,7*0,7*0/                                    XSE  200
      DATA X01 /1.129E00,1.131E00,1.133E00,1.134E00,1.137E00,1.139E00,  XSE  210
     *     1.133E00,1.135E00,1.137E00,1.139E00,1.141E00,1.143E00,       XSE  220
     *     1.137E00,1.139E00,1.141E00,1.143E00,1.145E00,1.148E00,0.0,11*XSE  230
     *     0.0/                                                         XSE  240
      DATA X02 /1.099E01,1.101E01,1.103E01,1.104E01,1.105E01,1.106E01,  XSE  250
     *     1.108E01,1.110E01,1.111E01,1.112E01,1.113E01,1.114E01,       XSE  260
     *     1.116E01,1.118E01,1.119E01,1.120E01,1.121E01,1.122E01,0.0,11*XSE  270
     *     0.0/                                                         XSE  280
      DATA X03 /6.241E00,6.260E00,6.275E00,6.292E00,6.309E00,6.328E00,  XSE  290
     *     6.296E00,6.315E00,6.331E00,6.347E00,6.365E00,6.384E00,       XSE  300
     *     6.353E00,6.372E00,6.388E00,6.404E00,6.422E00,6.441E00,0.0,11*XSE  310
     *     0.0/                                                         XSE  320
      DATA X04 /5.220E01,5.237E01,5.251E01,5.266E01,5.282E01,5.300E01,  XSE  330
     *     5.251E01,5.268E01,5.282E01,5.297E01,5.314E01,5.332E01,       XSE  340
     *     5.283E01,5.300E01,5.315E01,5.330E01,5.346E01,5.365E01,0.0,11*XSE  350
     *     0.0/                                                         XSE  360
      DATA X05 /1.421E01,1.425E01,1.429E01,1.432E01,1.435E01,1.439E01,  XSE  370
     *     1.421E01,1.425E01,1.428E01,1.432E01,1.435E01,1.439E01,       XSE  380
     *     1.418E01,1.422E01,1.425E01,1.428E01,1.432E01,1.436E01,0.0,11*XSE  390
     *     0.0/                                                         XSE  400
      DATA X06 /9.983E00,1.001E01,1.004E01,1.007E01,1.010E01,1.013E01,  XSE  410
     *     1.001E01,1.004E01,1.007E01,1.009E01,1.013E01,1.016E01,       XSE  420
     *     1.004E01,1.007E01,1.009E01,1.012E01,1.015E01,1.019E01,0.0,11*XSE  430
     *     0.0/                                                         XSE  440
      DATA X07 /4.329E01,4.345E01,4.358E01,4.372E01,4.388E01,4.407E01,  XSE  450
     *     4.348E01,4.364E01,4.377E01,4.392E01,4.408E01,4.426E01,       XSE  460
     *     4.368E01,4.384E01,4.397E01,4.412E01,4.428E01,4.446E01,0.0,11*XSE  470
     *     0.0/                                                         XSE  480
      DATA X08 /7.996E00,8.009E00,8.021E00,8.033E00,8.043E00,8.051E00,  XSE  490
     *     7.991E00,8.004E00,8.017E00,8.028E00,8.038E00,8.047E00,       XSE  500
     *     7.991E00,8.004E00,8.016E00,8.027E00,8.038E00,8.046E00,0.0,11*XSE  510
     *     0.0/                                                         XSE  520
      DATA X09 /1.516E00,1.517E00,1.519E00,1.521E00,1.522E00,1.524E00,  XSE  530
     *     1.523E00,1.525E00,1.527E00,1.529E00,1.530E00,1.532E00,       XSE  540
     *     1.531E00,1.533E00,1.535E00,1.537E00,1.538E00,1.540E00,0.0,11*XSE  550
     *     0.0/                                                         XSE  560
      DATA X10 /3.106E01,3.115E01,3.122E01,3.130E01,3.138E01,3.147E01,  XSE  570
     *     3.117E01,3.126E01,3.133E01,3.140E01,3.149E01,3.158E01,       XSE  580
     *     3.126E01,3.135E01,3.142E01,3.150E01,3.158E01,3.167E01,0.0,11*XSE  590
     *     0.0/                                                         XSE  600
      DATA X11 /3.381E01,3.395E01,3.405E01,3.417E01,3.431E01,3.446E01,  XSE  610
     *     3.401E01,3.414E01,3.425E01,3.437E01,3.450E01,3.466E01,       XSE  620
     *     3.421E01,3.434E01,3.445E01,3.457E01,3.471E01,3.486E01,0.0,11*XSE  630
     *     0.0/                                                         XSE  640
      DATA X12 /6.113E01,6.138E01,6.157E01,6.179E01,6.205E01,6.233E01,  XSE  650
     *     6.038E01,6.063E01,6.082E01,6.104E01,6.129E01,6.157E01,       XSE  660
     *     5.967E01,5.991E01,6.010E01,6.032E01,6.056E01,6.084E01,0.0,11*XSE  670
     *     0.0/                                                         XSE  680
      DATA X13 /1.074E02,1.078E02,1.082E02,1.086E02,1.090E02,1.095E02,  XSE  690
     *     1.065E02,1.069E02,1.073E02,1.077E02,1.081E02,1.086E02,       XSE  700
     *     1.056E02,1.061E02,1.064E02,1.068E02,1.072E02,1.077E02,0.0,11*XSE  710
     *     0.0/                                                         XSE  720
      DATA X14 /1.766E02,1.769E02,1.771E02,1.772E02,1.772E02,1.773E02,  XSE  730
     *     1.632E02,1.635E02,1.637E02,1.638E02,1.638E02,1.639E02,       XSE  740
     *     1.525E02,1.528E02,1.529E02,1.530E02,1.531E02,1.531E02,0.0,11*XSE  750
     *     0.0/                                                         XSE  760
      DATA X15 /3.752E01,3.767E01,3.780E01,3.793E01,3.809E01,3.826E01,  XSE  770
     *     3.740E01,3.755E01,3.767E01,3.781E01,3.796E01,3.814E01,       XSE  780
     *     3.729E01,3.744E01,3.756E01,3.770E01,3.785E01,3.802E01,0.0,11*XSE  790
     *     0.0/                                                         XSE  800
      DATA X16 /1.130E02,1.134E02,1.138E02,1.142E02,1.146E02,1.151E02,  XSE  810
     *     1.128E02,1.133E02,1.136E02,1.140E02,1.145E02,1.150E02,       XSE  820
     *     1.127E02,1.131E02,1.135E02,1.139E02,1.143E02,1.148E02,0.0,11*XSE  830
     *     0.0/                                                         XSE  840
      DATA X17 /3.175E01,3.182E01,3.189E01,3.195E01,3.200E01,3.205E01,  XSE  850
     *     3.174E01,3.181E01,3.188E01,3.194E01,3.200E01,3.205E01,       XSE  860
     *     3.174E01,3.181E01,3.188E01,3.194E01,3.199E01,3.204E01,0.0,11*XSE  870
     *     0.0/                                                         XSE  880
      DATA X18 /9.802E01,9.838E01,9.867E01,9.899E01,9.933E01,9.972E01,  XSE  890
     *     9.750E01,9.786E01,9.815E01,9.846E01,9.881E01,9.920E01,       XSE  900
     *     9.699E01,9.735E01,9.763E01,9.794E01,9.829E01,9.867E01,0.0,11*XSE  910
     *     0.0/                                                         XSE  920
      DATA X19 /1.211E01,1.216E01,1.220E01,1.223E01,1.228E01,1.233E01,  XSE  930
     *     1.205E01,1.210E01,1.213E01,1.217E01,1.221E01,1.226E01,       XSE  940
     *     1.199E01,1.203E01,1.207E01,1.211E01,1.215E01,1.220E01,0.0,11*XSE  950
     *     0.0/                                                         XSE  960
      DATA X20 /30*0.0/                                                 PC
      NPT = MPT                                                         XSE  970
      NXSEC = MXSEC                                                     XSE  980
      DO 10 I=1,30                                                      XSE  990
         F(I,01) = X01(I)                                               XSE 1000
         F(I,02) = X02(I)                                               XSE 1010
         F(I,03) = X03(I)                                               XSE 1020
         F(I,04) = X04(I)                                               XSE 1030
         F(I,05) = X05(I)                                               XSE 1040
         F(I,06) = X06(I)                                               XSE 1050
         F(I,07) = X07(I)                                               XSE 1060
         F(I,08) = X08(I)                                               XSE 1070
         F(I,09) = X09(I)                                               XSE 1080
         F(I,10) = X10(I)                                               XSE 1090
         F(I,11) = X11(I)                                               XSE 1100
         F(I,12) = X12(I)                                               XSE 1110
         F(I,13) = X13(I)                                               XSE 1120
         F(I,14) = X14(I)                                               XSE 1130
         F(I,15) = X15(I)                                               XSE 1140
         F(I,16) = X16(I)                                               XSE 1150
         F(I,17) = X17(I)                                               XSE 1160
         F(I,18) = X18(I)                                               XSE 1170
         F(I,19) = X19(I)                                               XSE 1180
         F(I,20) = X20(I)                                               XSE 1190
         CBUP(I) = DBUP(I)                                              XSE 1200
   10 CONTINUE                                                          XSE 1210
      DO 30 I=1,20                                                      XSE 1220
         DO 20 J=1,7                                                    XSE 1230
            NFUD(J,I) = MFUD(J,I)                                       XSE 1240
   20    CONTINUE                                                       XSE 1250
   30 CONTINUE                                                          XSE 1260
      RETURN                                                            XSE 1270
      END                                                               XSE 1280
      SUBROUTINE XSEC08(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PTPU   ;XSDRNPM/CITATION,XSEC FROM CITATION                        XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION MFUD(7,21), DBUP(30)                                    XSE   60
      DIMENSION X01(30), X02(30), X03(30), X04(30), X05(30), X06(30),   XSE   70
     *     X07(30), X08(30), X09(30), X10(30), X11(30), X12(30),        XSE   80
     *     X13(30), X14(30), X15(30), X16(30), X17(30), X18(30),        XSE   90
     *     X19(30), X20(30)                                             XSE  100
      DATA MXSEC, MPT /19,18/                                           XSE  110
      DATA DBUP /0.0,5.494E-01,1.099E00,1.650E00,2.201E00,2.753E00,     XSE  120
     *     2.753E00,3.271E00,3.789E00,4.309E00,4.830E00,5.352E00,       XSE  130
     *     5.352E00,5.842E00,6.334E00,6.827E00,7.321E00,7.817E00,0.0,11*XSE  140
     *     0.0/                                                         XSE  150
      DATA MFUD /902320,1,5*0,912330,1,5*0,922330,1,5*0,922330,4,5*0,   XSE  160
     *     922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,922380,1,XSE  170
     *     5*0,932370,1,5*0,942380,1,5*0,942390,1,5*0,942390,4,5*0,     XSE  180
     *     942400,1,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,952410,1,XSE  190
     *     5*0,952410,2,5*0,7*0,7*0/                                    XSE  200
      DATA X01 /7.376E-01,7.395E-01,7.415E-01,7.435E-01,7.457E-01,      XSE  210
     *     7.478E-01,7.421E-01,7.441E-01,7.462E-01,7.483E-01,7.505E-01, XSE  220
     *     7.528E-01,7.469E-01,7.490E-01,7.511E-01,7.533E-01,7.556E-01, XSE  230
     *     7.580E-01,0.0,11*0.0/                                        XSE  240
      DATA X02 /1.002E01,1.002E01,1.003E01,1.003E01,1.004E01,1.005E01,  XSE  250
     *     9.983E00,9.985E00,9.989E00,9.993E00,9.999E00,1.001E01,       XSE  260
     *     9.945E00,9.947E00,9.950E00,9.955E00,9.961E00,9.969E00,0.0,11*XSE  270
     *     0.0/                                                         XSE  280
      DATA X03 /3.634E00,3.648E00,3.662E00,3.677E00,3.693E00,3.710E00,  XSE  290
     *     3.617E00,3.632E00,3.647E00,3.662E00,3.679E00,3.696E00,       XSE  300
     *     3.604E00,3.620E00,3.635E00,3.652E00,3.669E00,3.687E00,0.0,11*XSE  310
     *     0.0/                                                         XSE  320
      DATA X04 /2.338E01,2.353E01,2.369E01,2.385E01,2.402E01,2.420E01,  XSE  330
     *     2.348E01,2.364E01,2.380E01,2.397E01,2.415E01,2.433E01,       XSE  340
     *     2.361E01,2.377E01,2.394E01,2.412E01,2.430E01,2.450E01,0.0,11*XSE  350
     *     0.0/                                                         XSE  360
      DATA X05 /1.704E01,1.706E01,1.708E01,1.709E01,1.711E01,1.713E01,  XSE  370
     *     1.702E01,1.704E01,1.706E01,1.708E01,1.710E01,1.712E01,       XSE  380
     *     1.701E01,1.703E01,1.705E01,1.707E01,1.709E01,1.711E01,0.0,11*XSE  390
     *     0.0/                                                         XSE  400
      DATA X06 /4.768E00,4.794E00,4.821E00,4.849E00,4.878E00,4.908E00,  XSE  410
     *     4.821E00,4.849E00,4.877E00,4.906E00,4.936E00,4.967E00,       XSE  420
     *     4.880E00,4.908E00,4.938E00,4.968E00,4.999E00,5.032E00,0.0,11*XSE  430
     *     0.0/                                                         XSE  440
      DATA X07 /1.336E01,1.351E01,1.367E01,1.383E01,1.400E01,1.418E01,  XSE  450
     *     1.361E01,1.376E01,1.393E01,1.410E01,1.428E01,1.446E01,       XSE  460
     *     1.388E01,1.404E01,1.421E01,1.439E01,1.457E01,1.476E01,0.0,11*XSE  470
     *     0.0/                                                         XSE  480
      DATA X08 /9.447E00,9.445E00,9.443E00,9.441E00,9.438E00,9.436E00,  XSE  490
     *     9.430E00,9.429E00,9.427E00,9.424E00,9.422E00,9.419E00,       XSE  500
     *     9.413E00,9.412E00,9.410E00,9.407E00,9.405E00,9.402E00,0.0,11*XSE  510
     *     0.0/                                                         XSE  520
      DATA X09 /7.636E00,7.635E00,7.634E00,7.632E00,7.630E00,7.628E00,  XSE  530
     *     7.658E00,7.657E00,7.656E00,7.654E00,7.652E00,7.650E00,       XSE  540
     *     7.681E00,7.680E00,7.679E00,7.677E00,7.675E00,7.673E00,0.0,11*XSE  550
     *     0.0/                                                         XSE  560
      DATA X10 /2.053E01,2.068E01,2.084E01,2.101E01,2.118E01,2.136E01,  XSE  570
     *     2.043E01,2.058E01,2.074E01,2.090E01,2.107E01,2.125E01,       XSE  580
     *     2.031E01,2.046E01,2.061E01,2.077E01,2.094E01,2.112E01,0.0,11*XSE  590
     *     0.0/                                                         XSE  600
      DATA X11 /8.225E00,8.340E00,8.458E00,8.581E00,8.708E00,8.841E00,  XSE  610
     *     8.407E00,8.527E00,8.650E00,8.778E00,8.910E00,9.049E00,       XSE  620
     *     8.605E00,8.730E00,8.858E00,8.991E00,9.129E00,9.273E00,0.0,11*XSE  630
     *     0.0/                                                         XSE  640
      DATA X12 /1.333E01,1.357E01,1.381E01,1.406E01,1.432E01,1.459E01,  XSE  650
     *     1.391E01,1.415E01,1.441E01,1.467E01,1.494E01,1.523E01,       XSE  660
     *     1.455E01,1.481E01,1.508E01,1.536E01,1.565E01,1.595E01,0.0,11*XSE  670
     *     0.0/                                                         XSE  680
      DATA X13 /2.330E01,2.371E01,2.415E01,2.459E01,2.506E01,2.554E01,  XSE  690
     *     2.420E01,2.464E01,2.509E01,2.556E01,2.605E01,2.656E01,       XSE  700
     *     2.521E01,2.567E01,2.615E01,2.664E01,2.716E01,2.769E01,0.0,11*XSE  710
     *     0.0/                                                         XSE  720
      DATA X14 /1.999E01,2.010E01,2.023E01,2.036E01,2.051E01,2.066E01,  XSE  730
     *     2.049E01,2.061E01,2.073E01,2.087E01,2.102E01,2.117E01,       XSE  740
     *     2.102E01,2.115E01,2.128E01,2.142E01,2.157E01,2.173E01,0.0,11*XSE  750
     *     0.0/                                                         XSE  760
      DATA X15 /8.402E00,8.540E00,8.681E00,8.828E00,8.981E00,9.140E00,  XSE  770
     *     8.674E00,8.818E00,8.967E00,9.121E00,9.281E00,9.448E00,       XSE  780
     *     8.976E00,9.128E00,9.283E00,9.445E00,9.614E00,9.789E00,0.0,11*XSE  790
     *     0.0/                                                         XSE  800
      DATA X16 /3.042E01,3.082E01,3.122E01,3.164E01,3.208E01,3.253E01,  XSE  810
     *     3.119E01,3.160E01,3.203E01,3.247E01,3.292E01,3.340E01,       XSE  820
     *     3.204E01,3.247E01,3.291E01,3.337E01,3.385E01,3.435E01,0.0,11*XSE  830
     *     0.0/                                                         XSE  840
      DATA X17 /9.640E00,9.642E00,9.643E00,9.645E00,9.646E00,9.648E00,  XSE  850
     *     9.677E00,9.680E00,9.681E00,9.683E00,9.685E00,9.686E00,       XSE  860
     *     9.715E00,9.718E00,9.720E00,9.722E00,9.724E00,9.726E00,0.0,11*XSE  870
     *     0.0/                                                         XSE  880
      DATA X18 /3.303E01,3.344E01,3.387E01,3.431E01,3.478E01,3.526E01,  XSE  890
     *     3.336E01,3.378E01,3.422E01,3.467E01,3.514E01,3.564E01,       XSE  900
     *     3.374E01,3.417E01,3.461E01,3.507E01,3.556E01,3.606E01,0.0,11*XSE  910
     *     0.0/                                                         XSE  920
      DATA X19 /4.082E00,4.133E00,4.186E00,4.241E00,4.298E00,4.358E00,  XSE  930
     *     4.123E00,4.175E00,4.229E00,4.285E00,4.344E00,4.405E00,       XSE  940
     *     4.170E00,4.223E00,4.278E00,4.335E00,4.395E00,4.457E00,0.0,11*XSE  950
     *     0.0/                                                         XSE  960
      DATA X20  /30*0.0/                                                PC
      NPT = MPT                                                         XSE  970
      NXSEC = MXSEC                                                     XSE  980
      DO 10 I=1,30                                                      XSE  990
         F(I,01) = X01(I)                                               XSE 1000
         F(I,02) = X02(I)                                               XSE 1010
         F(I,03) = X03(I)                                               XSE 1020
         F(I,04) = X04(I)                                               XSE 1030
         F(I,05) = X05(I)                                               XSE 1040
         F(I,06) = X06(I)                                               XSE 1050
         F(I,07) = X07(I)                                               XSE 1060
         F(I,08) = X08(I)                                               XSE 1070
         F(I,09) = X09(I)                                               XSE 1080
         F(I,10) = X10(I)                                               XSE 1090
         F(I,11) = X11(I)                                               XSE 1100
         F(I,12) = X12(I)                                               XSE 1110
         F(I,13) = X13(I)                                               XSE 1120
         F(I,14) = X14(I)                                               XSE 1130
         F(I,15) = X15(I)                                               XSE 1140
         F(I,16) = X16(I)                                               XSE 1150
         F(I,17) = X17(I)                                               XSE 1160
         F(I,18) = X18(I)                                               XSE 1170
         F(I,19) = X19(I)                                               XSE 1180
         F(I,20) = X20(I)                                               XSE 1190
         CBUP(I) = DBUP(I)                                              XSE 1200
   10 CONTINUE                                                          XSE 1210
      DO 30 I=1,20                                                      XSE 1220
         DO 20 J=1,7                                                    XSE 1230
            NFUD(J,I) = MFUD(J,I)                                       XSE 1240
   20    CONTINUE                                                       XSE 1250
   30 CONTINUE                                                          XSE 1260
      RETURN                                                            XSE 1270
      END                                                               XSE 1280
      SUBROUTINE XSEC09(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PUD50  ;XSDRNPM/CITATION,XSEC FROM CITATION                        XSE   30
C                                                                       XSE   40
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   50
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   60
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   70
     *     X19(32), X20(32)                                             XSE   80
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE   90
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE  100
      DATA MXSEC, MPT /20,30/                                           XSE  110
      DATA DBUP /0.0,5.786E-01,1.165E00,1.757E00,2.353E00,2.952E00,     XSE  120
     *     2.952E00,3.500E00,4.048E00,4.595E00,5.141E00,5.686E00,       XSE  130
     *     5.686E00,6.184E00,6.680E00,7.175E00,7.670E00,8.164E00,       XSE  140
     *     8.164E00,8.615E00,9.065E00,9.516E00,9.967E00,1.042E01,       XSE  150
     *     1.042E01,1.083E01,1.125E01,1.166E01,1.208E01,1.250E01,2*0.0/ XSE  160
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  170
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  180
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  190
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  200
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  210
      DATA X01 /2.066E01,2.068E01,2.070E01,2.073E01,2.075E01,2.078E01,  XSE  220
     *     2.065E01,2.067E01,2.069E01,2.071E01,2.074E01,2.077E01,       XSE  230
     *     2.066E01,2.068E01,2.070E01,2.073E01,2.076E01,2.079E01,       XSE  240
     *     2.070E01,2.071E01,2.074E01,2.076E01,2.079E01,2.082E01,       XSE  250
     *     2.075E01,2.077E01,2.079E01,2.082E01,2.085E01,2.088E01,2*0.0/ XSE  260
      DATA X02 /9.374E00,9.372E00,9.377E00,9.389E00,9.408E00,9.432E00,  XSE  270
     *     9.489E00,9.486E00,9.492E00,9.504E00,9.523E00,9.547E00,       XSE  280
     *     9.619E00,9.616E00,9.622E00,9.634E00,9.653E00,9.678E00,       XSE  290
     *     9.764E00,9.762E00,9.767E00,9.780E00,9.799E00,9.825E00,       XSE  300
     *     9.927E00,9.925E00,9.930E00,9.944E00,9.963E00,9.989E00,2*0.0/ XSE  310
      DATA X03 /4.022E01,4.018E01,4.020E01,4.025E01,4.035E01,4.047E01,  XSE  320
     *     4.079E01,4.076E01,4.078E01,4.083E01,4.093E01,4.106E01,       XSE  330
     *     4.145E01,4.142E01,4.143E01,4.149E01,4.159E01,4.172E01,       XSE  340
     *     4.218E01,4.215E01,4.217E01,4.223E01,4.232E01,4.246E01,       XSE  350
     *     4.300E01,4.297E01,4.298E01,4.304E01,4.315E01,4.328E01,2*0.0/ XSE  360
      DATA X04 /8.107E00,8.119E00,8.129E00,8.138E00,8.145E00,8.150E00,  XSE  370
     *     7.668E00,7.679E00,7.688E00,7.697E00,7.703E00,7.708E00,       XSE  380
     *     7.302E00,7.313E00,7.322E00,7.329E00,7.336E00,7.340E00,       XSE  390
     *     6.992E00,7.002E00,7.011E00,7.018E00,7.024E00,7.028E00,       XSE  400
     *     6.725E00,6.734E00,6.742E00,6.749E00,6.755E00,6.759E00,2*0.0/ XSE  410
      DATA X05 /9.196E-01,9.203E-01,9.211E-01,9.219E-01,9.227E-01,      XSE  420
     *     9.235E-01,9.249E-01,9.256E-01,9.264E-01,9.272E-01,9.280E-01, XSE  430
     *     9.289E-01,9.308E-01,9.315E-01,9.323E-01,9.331E-01,9.339E-01, XSE  440
     *     9.348E-01,9.371E-01,9.378E-01,9.386E-01,9.394E-01,9.403E-01, XSE  450
     *     9.411E-01,9.437E-01,9.444E-01,9.452E-01,9.461E-01,9.469E-01, XSE  460
     *     9.478E-01,2*0.0/                                             XSE  470
      DATA X06 /3.075E01,3.073E01,3.072E01,3.073E01,3.075E01,3.079E01,  XSE  480
     *     3.093E01,3.090E01,3.090E01,3.090E01,3.093E01,3.097E01,       XSE  490
     *     3.104E01,3.101E01,3.101E01,3.101E01,3.104E01,3.108E01,       XSE  500
     *     3.111E01,3.109E01,3.108E01,3.109E01,3.111E01,3.115E01,       XSE  510
     *     3.116E01,3.113E01,3.112E01,3.113E01,3.115E01,3.119E01,2*0.0/ XSE  520
      DATA X07 /2.899E01,2.896E01,2.898E01,2.903E01,2.910E01,2.920E01,  XSE  530
     *     2.962E01,2.960E01,2.961E01,2.966E01,2.974E01,2.984E01,       XSE  540
     *     3.031E01,3.029E01,3.030E01,3.035E01,3.043E01,3.054E01,       XSE  550
     *     3.106E01,3.104E01,3.106E01,3.111E01,3.119E01,3.130E01,       XSE  560
     *     3.189E01,3.186E01,3.188E01,3.193E01,3.202E01,3.213E01,2*0.0/ XSE  570
      DATA X08 /2.177E00,2.176E00,2.177E00,2.178E00,2.180E00,2.183E00,  XSE  580
     *     2.197E00,2.197E00,2.197E00,2.199E00,2.201E00,2.203E00,       XSE  590
     *     2.220E00,2.219E00,2.220E00,2.221E00,2.223E00,2.226E00,       XSE  600
     *     2.244E00,2.243E00,2.244E00,2.245E00,2.248E00,2.250E00,       XSE  610
     *     2.270E00,2.270E00,2.270E00,2.272E00,2.274E00,2.277E00,2*0.0/ XSE  620
      DATA X09 /5.893E01,5.887E01,5.890E01,5.900E01,5.916E01,5.938E01,  XSE  630
     *     5.740E01,5.735E01,5.738E01,5.747E01,5.763E01,5.784E01,       XSE  640
     *     5.599E01,5.594E01,5.596E01,5.606E01,5.621E01,5.642E01,       XSE  650
     *     5.468E01,5.464E01,5.466E01,5.475E01,5.490E01,5.510E01,       XSE  660
     *     5.347E01,5.342E01,5.345E01,5.354E01,5.368E01,5.388E01,2*0.0/ XSE  670
      DATA X10 /1.026E02,1.025E02,1.025E02,1.027E02,1.030E02,1.034E02,  XSE  680
     *     1.008E02,1.007E02,1.008E02,1.009E02,1.012E02,1.016E02,       XSE  690
     *     9.925E01,9.916E01,9.921E01,9.937E01,9.964E01,1.000E02,       XSE  700
     *     9.788E01,9.778E01,9.783E01,9.799E01,9.825E01,9.861E01,       XSE  710
     *     9.667E01,9.658E01,9.662E01,9.678E01,9.704E01,9.739E01,2*0.0/ XSE  720
      DATA X11 /1.767E02,1.759E02,1.750E02,1.741E02,1.731E02,1.721E02,  XSE  730
     *     1.433E02,1.427E02,1.420E02,1.413E02,1.405E02,1.398E02,       XSE  740
     *     1.242E02,1.237E02,1.232E02,1.225E02,1.219E02,1.214E02,       XSE  750
     *     1.118E02,1.114E02,1.109E02,1.103E02,1.098E02,1.093E02,       XSE  760
     *     1.030E02,1.026E02,1.022E02,1.017E02,1.013E02,1.009E02,2*0.0/ XSE  770
      DATA X12 /3.533E01,3.530E01,3.532E01,3.538E01,3.547E01,3.560E01,  XSE  780
     *     3.516E01,3.513E01,3.515E01,3.521E01,3.530E01,3.543E01,       XSE  790
     *     3.505E01,3.503E01,3.504E01,3.510E01,3.520E01,3.533E01,       XSE  800
     *     3.501E01,3.498E01,3.500E01,3.506E01,3.516E01,3.529E01,       XSE  810
     *     3.504E01,3.501E01,3.502E01,3.508E01,3.518E01,3.531E01,2*0.0/ XSE  820
      DATA X13 /1.066E02,1.065E02,1.065E02,1.067E02,1.070E02,1.074E02,  XSE  830
     *     1.065E02,1.064E02,1.065E02,1.066E02,1.069E02,1.073E02,       XSE  840
     *     1.066E02,1.065E02,1.066E02,1.067E02,1.070E02,1.074E02,       XSE  850
     *     1.068E02,1.068E02,1.068E02,1.070E02,1.073E02,1.076E02,       XSE  860
     *     1.073E02,1.072E02,1.073E02,1.074E02,1.077E02,1.081E02,2*0.0/ XSE  870
      DATA X14 /3.298E01,3.303E01,3.307E01,3.311E01,3.314E01,3.316E01,  XSE  880
     *     3.182E01,3.187E01,3.192E01,3.195E01,3.198E01,3.200E01,       XSE  890
     *     3.078E01,3.082E01,3.087E01,3.090E01,3.093E01,3.095E01,       XSE  900
     *     2.982E01,2.987E01,2.991E01,2.994E01,2.997E01,2.999E01,       XSE  910
     *     2.894E01,2.899E01,2.903E01,2.906E01,2.908E01,2.910E01,2*0.0/ XSE  920
      DATA X15 /9.367E01,9.357E01,9.356E01,9.363E01,9.378E01,9.400E01,  XSE  930
     *     9.220E01,9.210E01,9.208E01,9.214E01,9.229E01,9.252E01,       XSE  940
     *     9.075E01,9.066E01,9.064E01,9.070E01,9.084E01,9.106E01,       XSE  950
     *     8.938E01,8.928E01,8.926E01,8.932E01,8.946E01,8.967E01,       XSE  960
     *     8.806E01,8.796E01,8.794E01,8.800E01,8.814E01,8.834E01,2*0.0/ XSE  970
      DATA X16 /1.158E01,1.156E01,1.156E01,1.157E01,1.159E01,1.162E01,  XSE  980
     *     1.139E01,1.138E01,1.138E01,1.139E01,1.141E01,1.143E01,       XSE  990
     *     1.122E01,1.120E01,1.120E01,1.121E01,1.123E01,1.125E01,       XSE 1000
     *     1.105E01,1.103E01,1.103E01,1.104E01,1.106E01,1.108E01,       XSE 1010
     *     1.088E01,1.087E01,1.087E01,1.088E01,1.089E01,1.092E01,2*0.0/ XSE 1020
      DATA X17 /2.454E00,2.447E00,2.438E00,2.428E00,2.418E00,2.408E00,  XSE 1030
     *     2.488E00,2.480E00,2.471E00,2.461E00,2.451E00,2.441E00,       XSE 1040
     *     2.506E00,2.499E00,2.490E00,2.480E00,2.469E00,2.459E00,       XSE 1050
     *     2.518E00,2.511E00,2.501E00,2.491E00,2.481E00,2.471E00,       XSE 1060
     *     2.526E00,2.519E00,2.510E00,2.499E00,2.489E00,2.479E00,2*0.0/ XSE 1070
      DATA X18 /4.662E01,4.649E01,4.632E01,4.613E01,4.594E01,4.576E01,  XSE 1080
     *     4.727E01,4.713E01,4.695E01,4.677E01,4.657E01,4.638E01,       XSE 1090
     *     4.762E01,4.748E01,4.731E01,4.712E01,4.692E01,4.673E01,       XSE 1100
     *     4.784E01,4.770E01,4.753E01,4.733E01,4.714E01,4.695E01,       XSE 1110
     *     4.800E01,4.786E01,4.768E01,4.749E01,4.729E01,4.710E01,2*0.0/ XSE 1120
      DATA X19 /5.669E00,5.674E00,5.680E00,5.686E00,5.693E00,5.699E00,  XSE 1130
     *     5.686E00,5.692E00,5.698E00,5.704E00,5.710E00,5.717E00,       XSE 1140
     *     5.704E00,5.709E00,5.715E00,5.722E00,5.728E00,5.735E00,       XSE 1150
     *     5.722E00,5.727E00,5.733E00,5.740E00,5.746E00,5.753E00,       XSE 1160
     *     5.741E00,5.746E00,5.752E00,5.758E00,5.765E00,5.772E00,2*0.0/ XSE 1170
      DATA X20 /1.384E01,1.386E01,1.388E01,1.389E01,1.390E01,1.391E01,  XSE 1180
     *     1.386E01,1.388E01,1.390E01,1.392E01,1.393E01,1.394E01,       XSE 1190
     *     1.389E01,1.391E01,1.392E01,1.394E01,1.395E01,1.396E01,       XSE 1200
     *     1.391E01,1.393E01,1.394E01,1.396E01,1.397E01,1.398E01,       XSE 1210
     *     1.392E01,1.395E01,1.396E01,1.398E01,1.399E01,1.400E01,2*0.0/ XSE 1220
      NPT = MPT                                                         XSE 1230
      NXSEC = MXSEC                                                     XSE 1240
      DO 10 I=1,30                                                      XSE 1250
         F(I,01) = X01(I)                                               XSE 1260
         F(I,02) = X02(I)                                               XSE 1270
         F(I,03) = X03(I)                                               XSE 1280
         F(I,04) = X04(I)                                               XSE 1290
         F(I,05) = X05(I)                                               XSE 1300
         F(I,06) = X06(I)                                               XSE 1310
         F(I,07) = X07(I)                                               XSE 1320
         F(I,08) = X08(I)                                               XSE 1330
         F(I,09) = X09(I)                                               XSE 1340
         F(I,10) = X10(I)                                               XSE 1350
         F(I,11) = X11(I)                                               XSE 1360
         F(I,12) = X12(I)                                               XSE 1370
         F(I,13) = X13(I)                                               XSE 1380
         F(I,14) = X14(I)                                               XSE 1390
         F(I,15) = X15(I)                                               XSE 1400
         F(I,16) = X16(I)                                               XSE 1410
         F(I,17) = X17(I)                                               XSE 1420
         F(I,18) = X18(I)                                               XSE 1430
         F(I,19) = X19(I)                                               XSE 1440
         F(I,20) = X20(I)                                               XSE 1450
         CBUP(I) = DBUP(I)                                              XSE 1460
   10 CONTINUE                                                          XSE 1470
      DO 30 I=1,20                                                      XSE 1480
         DO 20 J=1,7                                                    XSE 1490
            NFUD(J,I) = MFUD(J,I)                                       XSE 1500
   20    CONTINUE                                                       XSE 1510
   30 CONTINUE                                                          XSE 1520
      RETURN                                                            XSE 1530
      END                                                               XSE 1540
      SUBROUTINE XSEC10(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PTD5D35 ;XSDRNPM/CITATION, XSEC FROM CITAT                         XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /19,18/                                           XSE  110
      DATA DBUP /0.0,3.059E-01,6.421E-01,1.005E00,1.393E00,1.804E00,    XSE  120
     *     1.804E00,2.298E00,2.802E00,3.314E00,3.833E00,4.358E00,       XSE  130
     *     4.358E00,5.094E00,5.803E00,6.489E00,7.156E00,7.803E00,0.0,13*XSE  140
     *     0.0/                                                         XSE  150
      DATA MFUD /902320,1,5*0,912330,1,5*0,922330,1,5*0,922330,4,5*0,   XSE  160
     *     922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,922380,1,XSE  170
     *     5*0,932370,1,5*0,942380,1,5*0,942390,1,5*0,942390,4,5*0,     XSE  180
     *     942400,1,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,952410,1,XSE  190
     *     5*0,952410,2,5*0,14*0/                                       XSE  200
      DATA X01 /1.047E00,1.047E00,1.047E00,1.048E00,1.049E00,1.050E00,  XSE  210
     *     1.056E00,1.056E00,1.058E00,1.060E00,1.063E00,1.065E00,       XSE  220
     *     1.076E00,1.081E00,1.085E00,1.089E00,1.094E00,1.098E00,0.0,13*XSE  230
     *     0.0/                                                         XSE  240
      DATA X02 /1.191E01,1.189E01,1.186E01,1.183E01,1.180E01,1.177E01,  XSE  250
     *     1.177E01,1.173E01,1.170E01,1.167E01,1.164E01,1.161E01,       XSE  260
     *     1.160E01,1.158E01,1.156E01,1.155E01,1.154E01,1.153E01,0.0,13*XSE  270
     *     0.0/                                                         XSE  280
      DATA X03 /6.304E00,6.299E00,6.298E00,6.298E00,6.299E00,6.303E00,  XSE  290
     *     6.246E00,6.244E00,6.251E00,6.262E00,6.275E00,6.291E00,       XSE  300
     *     6.176E00,6.211E00,6.238E00,6.266E00,6.295E00,6.325E00,0.0,13*XSE  310
     *     0.0/                                                         XSE  320
      DATA X04 /4.971E01,4.970E01,4.971E01,4.974E01,4.978E01,4.985E01,  XSE  330
     *     4.955E01,4.956E01,4.967E01,4.981E01,4.997E01,5.017E01,       XSE  340
     *     4.977E01,5.014E01,5.044E01,5.073E01,5.105E01,5.137E01,0.0,13*XSE  350
     *     0.0/                                                         XSE  360
      DATA X05 /1.994E01,1.993E01,1.992E01,1.992E01,1.992E01,1.992E01,  XSE  370
     *     1.849E01,1.849E01,1.851E01,1.853E01,1.855E01,1.858E01,       XSE  380
     *     1.452E01,1.459E01,1.464E01,1.470E01,1.475E01,1.481E01,0.0,13*XSE  390
     *     0.0/                                                         XSE  400
      DATA X06 /9.025E00,9.025E00,9.032E00,9.041E00,9.054E00,9.071E00,  XSE  410
     *     9.095E00,9.103E00,9.128E00,9.159E00,9.195E00,9.235E00,       XSE  420
     *     9.282E00,9.355E00,9.413E00,9.471E00,9.532E00,9.594E00,0.0,13*XSE  430
     *     0.0/                                                         XSE  440
      DATA X07 /3.832E01,3.833E01,3.838E01,3.843E01,3.851E01,3.862E01,  XSE  450
     *     3.857E01,3.863E01,3.877E01,3.895E01,3.915E01,3.938E01,       XSE  460
     *     3.954E01,3.995E01,4.027E01,4.060E01,4.094E01,4.129E01,0.0,13*XSE  470
     *     0.0/                                                         XSE  480
      DATA X08 /8.520E00,8.510E00,8.503E00,8.497E00,8.491E00,8.485E00,  XSE  490
     *     7.496E00,7.490E00,7.488E00,7.486E00,7.483E00,7.480E00,       XSE  500
     *     5.234E00,5.238E00,5.241E00,5.244E00,5.246E00,5.247E00,0.0,13*XSE  510
     *     0.0/                                                         XSE  520
      DATA X09 /1.766E00,1.765E00,1.764E00,1.763E00,1.763E00,1.762E00,  XSE  530
     *     1.787E00,1.786E00,1.786E00,1.786E00,1.787E00,1.787E00,       XSE  540
     *     1.731E00,1.732E00,1.734E00,1.735E00,1.737E00,1.737E00,0.0,13*XSE  550
     *     0.0/                                                         XSE  560
      DATA X10 /2.992E01,2.990E01,2.990E01,2.991E01,2.992E01,2.994E01,  XSE  570
     *     2.958E01,2.957E01,2.961E01,2.967E01,2.974E01,2.982E01,       XSE  580
     *     2.916E01,2.933E01,2.946E01,2.959E01,2.973E01,2.988E01,0.0,13*XSE  590
     *     0.0/                                                         XSE  600
      DATA X11 /2.739E01,2.740E01,2.744E01,2.749E01,2.755E01,2.764E01,  XSE  610
     *     2.772E01,2.777E01,2.789E01,2.803E01,2.819E01,2.837E01,       XSE  620
     *     2.865E01,2.897E01,2.923E01,2.948E01,2.975E01,3.003E01,0.0,13*XSE  630
     *     0.0/                                                         XSE  640
      DATA X12 /5.689E01,5.692E01,5.701E01,5.712E01,5.726E01,5.745E01,  XSE  650
     *     5.565E01,5.575E01,5.599E01,5.630E01,5.664E01,5.703E01,       XSE  660
     *     5.506E01,5.571E01,5.623E01,5.675E01,5.729E01,5.785E01,0.0,13*XSE  670
     *     0.0/                                                         XSE  680
      DATA X13 /9.865E01,9.870E01,9.885E01,9.904E01,9.929E01,9.962E01,  XSE  690
     *     9.702E01,9.720E01,9.763E01,9.816E01,9.877E01,9.945E01,       XSE  700
     *     9.679E01,9.794E01,9.885E01,9.977E01,1.007E02,1.017E02,0.0,13*XSE  710
     *     0.0/                                                         XSE  720
      DATA X14 /1.992E02,1.985E02,1.978E02,1.970E02,1.961E02,1.951E02,  XSE  730
     *     1.795E02,1.784E02,1.774E02,1.765E02,1.756E02,1.749E02,       XSE  740
     *     1.618E02,1.611E02,1.606E02,1.600E02,1.596E02,1.593E02,0.0,13*XSE  750
     *     0.0/                                                         XSE  760
      DATA X15 /3.376E01,3.378E01,3.383E01,3.389E01,3.398E01,3.409E01,  XSE  770
     *     3.352E01,3.358E01,3.372E01,3.391E01,3.411E01,3.434E01,       XSE  780
     *     3.384E01,3.424E01,3.456E01,3.488E01,3.521E01,3.555E01,0.0,13*XSE  790
     *     0.0/                                                         XSE  800
      DATA X16 /1.021E02,1.021E02,1.022E02,1.024E02,1.026E02,1.029E02,  XSE  810
     *     1.015E02,1.017E02,1.021E02,1.026E02,1.032E02,1.039E02,       XSE  820
     *     1.027E02,1.038E02,1.047E02,1.056E02,1.066E02,1.076E02,0.0,13*XSE  830
     *     0.0/                                                         XSE  840
      DATA X17 /3.305E01,3.300E01,3.297E01,3.294E01,3.291E01,3.288E01,  XSE  850
     *     3.283E01,3.280E01,3.278E01,3.276E01,3.274E01,3.273E01,       XSE  860
     *     3.258E01,3.260E01,3.262E01,3.263E01,3.264E01,3.264E01,0.0,13*XSE  870
     *     0.0/                                                         XSE  880
      DATA X18 /9.112E01,9.112E01,9.118E01,9.127E01,9.140E01,9.159E01,  XSE  890
     *     8.912E01,8.917E01,8.943E01,8.977E01,9.017E01,9.061E01,       XSE  900
     *     8.785E01,8.865E01,8.928E01,8.993E01,9.060E01,9.131E01,0.0,13*XSE  910
     *     0.0/                                                         XSE  920
      DATA X19 /1.126E01,1.126E01,1.127E01,1.128E01,1.130E01,1.132E01,  XSE  930
     *     1.101E01,1.102E01,1.105E01,1.109E01,1.114E01,1.120E01,       XSE  940
     *     1.086E01,1.096E01,1.104E01,1.111E01,1.120E01,1.128E01,0.0,13*XSE  950
     *     0.0/                                                         XSE  960
      DATA X20 /32*0.0/                                                 PC
      NPT = MPT                                                         XSE  970
      NXSEC = MXSEC                                                     XSE  980
      DO 10 I=1,30                                                      XSE  990
         F(I,01) = X01(I)                                               XSE 1000
         F(I,02) = X02(I)                                               XSE 1010
         F(I,03) = X03(I)                                               XSE 1020
         F(I,04) = X04(I)                                               XSE 1030
         F(I,05) = X05(I)                                               XSE 1040
         F(I,06) = X06(I)                                               XSE 1050
         F(I,07) = X07(I)                                               XSE 1060
         F(I,08) = X08(I)                                               XSE 1070
         F(I,09) = X09(I)                                               XSE 1080
         F(I,10) = X10(I)                                               XSE 1090
         F(I,11) = X11(I)                                               XSE 1100
         F(I,12) = X12(I)                                               XSE 1110
         F(I,13) = X13(I)                                               XSE 1120
         F(I,14) = X14(I)                                               XSE 1130
         F(I,15) = X15(I)                                               XSE 1140
         F(I,16) = X16(I)                                               XSE 1150
         F(I,17) = X17(I)                                               XSE 1160
         F(I,18) = X18(I)                                               XSE 1170
         F(I,19) = X19(I)                                               XSE 1180
         F(I,20) = X20(I)                                               XSE 1190
         CBUP(I) = DBUP(I)                                              XSE 1200
   10 CONTINUE                                                          XSE 1210
      DO 30 I=1,20                                                      XSE 1220
         DO 20 J=1,7                                                    XSE 1230
            NFUD(J,I) = MFUD(J,I)                                       XSE 1240
   20    CONTINUE                                                       XSE 1250
   30 CONTINUE                                                          XSE 1260
      RETURN                                                            XSE 1270
      END                                                               XSE 1280
      SUBROUTINE XSEC11(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PTD5D33 ;XSDRNPM/CITATION, XSEC FROM CITAT                         XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /19,18/                                           XSE  110
      DATA DBUP /0.0,8.195E-01,1.606E00,2.365E00,3.099E00,3.811E00,     XSE  120
     *     3.811E00,4.298E00,4.794E00,5.298E00,5.807E00,6.322E00,       XSE  130
     *     6.322E00,6.603E00,6.912E00,7.245E00,7.601E00,7.977E00,0.0,13*XSE  140
     *     0.0/                                                         XSE  150
      DATA MFUD /902320,1,5*0,912330,1,5*0,922330,1,5*0,922330,4,5*0,   XSE  160
     *     922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,922380,1,XSE  170
     *     5*0,932370,1,5*0,942380,1,5*0,942390,1,5*0,942390,4,5*0,     XSE  180
     *     942400,1,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,952410,1,XSE  190
     *     5*0,952410,2,5*0,14*0/                                       XSE  200
      DATA X01 /1.028E00,1.031E00,1.033E00,1.034E00,1.036E00,1.038E00,  XSE  210
     *     1.049E00,1.050E00,1.051E00,1.053E00,1.055E00,1.057E00,       XSE  220
     *     1.073E00,1.072E00,1.073E00,1.074E00,1.076E00,1.077E00,0.0,13*XSE  230
     *     0.0/                                                         XSE  240
      DATA X02 /1.089E01,1.091E01,1.090E01,1.089E01,1.087E01,1.085E01,  XSE  250
     *     1.102E01,1.100E01,1.099E01,1.097E01,1.096E01,1.095E01,       XSE  260
     *     1.114E01,1.112E01,1.111E01,1.110E01,1.109E01,1.108E01,0.0,13*XSE  270
     *     0.0/                                                         XSE  280
      DATA X03 /5.735E00,5.756E00,5.762E00,5.768E00,5.778E00,5.792E00,  XSE  290
     *     5.866E00,5.866E00,5.875E00,5.886E00,5.899E00,5.915E00,       XSE  300
     *     6.029E00,6.021E00,6.026E00,6.034E00,6.044E00,6.056E00,0.0,13*XSE  310
     *     0.0/                                                         XSE  320
      DATA X04 /4.618E01,4.636E01,4.640E01,4.645E01,4.655E01,4.669E01,  XSE  330
     *     4.730E01,4.731E01,4.741E01,4.753E01,4.767E01,4.783E01,       XSE  340
     *     4.886E01,4.880E01,4.886E01,4.895E01,4.906E01,4.919E01,0.0,13*XSE  350
     *     0.0/                                                         XSE  360
      DATA X05 /1.317E01,1.322E01,1.325E01,1.328E01,1.331E01,1.335E01,  XSE  370
     *     1.347E01,1.347E01,1.350E01,1.352E01,1.355E01,1.358E01,       XSE  380
     *     1.377E01,1.375E01,1.376E01,1.378E01,1.380E01,1.382E01,0.0,13*XSE  390
     *     0.0/                                                         XSE  400
      DATA X06 /8.619E00,8.656E00,8.669E00,8.686E00,8.710E00,8.742E00,  XSE  410
     *     8.866E00,8.873E00,8.894E00,8.919E00,8.948E00,8.980E00,       XSE  420
     *     9.181E00,9.171E00,9.184E00,9.202E00,9.224E00,9.249E00,0.0,13*XSE  430
     *     0.0/                                                         XSE  440
      DATA X07 /3.639E01,3.656E01,3.660E01,3.667E01,3.679E01,3.694E01,  XSE  450
     *     3.750E01,3.754E01,3.765E01,3.779E01,3.795E01,3.813E01,       XSE  460
     *     3.911E01,3.906E01,3.914E01,3.925E01,3.937E01,3.952E01,0.0,13*XSE  470
     *     0.0/                                                         XSE  480
      DATA X08 /4.868E00,4.882E00,4.894E00,4.904E00,4.913E00,4.921E00,  XSE  490
     *     4.910E00,4.910E00,4.913E00,4.915E00,4.918E00,4.920E00,       XSE  500
     *     4.913E00,4.909E00,4.908E00,4.907E00,4.906E00,4.906E00,0.0,13*XSE  510
     *     0.0/                                                         XSE  520
      DATA X09 /1.651E00,1.656E00,1.659E00,1.662E00,1.665E00,1.668E00,  XSE  530
     *     1.689E00,1.689E00,1.690E00,1.691E00,1.692E00,1.693E00,       XSE  540
     *     1.717E00,1.715E00,1.715E00,1.715E00,1.715E00,1.715E00,0.0,13*XSE  550
     *     0.0/                                                         XSE  560
      DATA X10 /2.780E01,2.789E01,2.792E01,2.794E01,2.799E01,2.806E01,  XSE  570
     *     2.819E01,2.819E01,2.823E01,2.829E01,2.835E01,2.843E01,       XSE  580
     *     2.874E01,2.870E01,2.872E01,2.876E01,2.881E01,2.887E01,0.0,13*XSE  590
     *     0.0/                                                         XSE  600
      DATA X11 /2.598E01,2.611E01,2.614E01,2.620E01,2.629E01,2.642E01,  XSE  610
     *     2.695E01,2.698E01,2.707E01,2.718E01,2.730E01,2.744E01,       XSE  620
     *     2.831E01,2.827E01,2.834E01,2.842E01,2.852E01,2.863E01,0.0,13*XSE  630
     *     0.0/                                                         XSE  640
      DATA X12 /5.361E01,5.389E01,5.396E01,5.408E01,5.427E01,5.454E01,  XSE  650
     *     5.390E01,5.397E01,5.415E01,5.438E01,5.465E01,5.494E01,       XSE  660
     *     5.490E01,5.484E01,5.497E01,5.514E01,5.535E01,5.558E01,0.0,13*XSE  670
     *     0.0/                                                         XSE  680
      DATA X13 /9.298E01,9.346E01,9.358E01,9.378E01,9.412E01,9.458E01,  XSE  690
     *     9.399E01,9.411E01,9.443E01,9.484E01,9.530E01,9.582E01,       XSE  700
     *     9.631E01,9.621E01,9.643E01,9.674E01,9.710E01,9.752E01,0.0,13*XSE  710
     *     0.0/                                                         XSE  720
      DATA X14 /1.823E02,1.823E02,1.815E02,1.805E02,1.795E02,1.786E02,  XSE  730
     *     1.656E02,1.649E02,1.643E02,1.638E02,1.634E02,1.630E02,       XSE  740
     *     1.539E02,1.535E02,1.532E02,1.529E02,1.527E02,1.525E02,0.0,13*XSE  750
     *     0.0/                                                         XSE  760
      DATA X15 /3.184E01,3.200E01,3.205E01,3.212E01,3.223E01,3.239E01,  XSE  770
     *     3.248E01,3.252E01,3.263E01,3.277E01,3.293E01,3.311E01,       XSE  780
     *     3.357E01,3.353E01,3.361E01,3.372E01,3.384E01,3.398E01,0.0,13*XSE  790
     *     0.0/                                                         XSE  800
      DATA X16 /9.611E01,9.660E01,9.674E01,9.695E01,9.729E01,9.774E01,  XSE  810
     *     9.827E01,9.838E01,9.870E01,9.910E01,9.955E01,1.000E02,       XSE  820
     *     1.017E02,1.016E02,1.018E02,1.021E02,1.024E02,1.028E02,0.0,13*XSE  830
     *     0.0/                                                         XSE  840
      DATA X17 /3.081E01,3.093E01,3.103E01,3.112E01,3.119E01,3.126E01,  XSE  850
     *     3.142E01,3.143E01,3.145E01,3.147E01,3.149E01,3.150E01,       XSE  860
     *     3.170E01,3.167E01,3.165E01,3.165E01,3.164E01,3.164E01,0.0,13*XSE  870
     *     0.0/                                                         XSE  880
      DATA X18 /8.488E01,8.526E01,8.532E01,8.543E01,8.564E01,8.595E01,  XSE  890
     *     8.545E01,8.550E01,8.571E01,8.598E01,8.630E01,8.667E01,       XSE  900
     *     8.698E01,8.686E01,8.700E01,8.721E01,8.746E01,8.775E01,0.0,13*XSE  910
     *     0.0/                                                         XSE  920
      DATA X19 /1.049E01,1.054E01,1.055E01,1.056E01,1.059E01,1.062E01,  XSE  930
     *     1.056E01,1.057E01,1.059E01,1.063E01,1.067E01,1.071E01,       XSE  940
     *     1.075E01,1.074E01,1.075E01,1.078E01,1.081E01,1.085E01,0.0,13*XSE  950
     *     0.0/                                                         XSE  960
      DATA X20 /32*0.0/                                                 PC
      NPT = MPT                                                         XSE  970
      NXSEC = MXSEC                                                     XSE  980
      DO 10 I=1,30                                                      XSE  990
         F(I,01) = X01(I)                                               XSE 1000
         F(I,02) = X02(I)                                               XSE 1010
         F(I,03) = X03(I)                                               XSE 1020
         F(I,04) = X04(I)                                               XSE 1030
         F(I,05) = X05(I)                                               XSE 1040
         F(I,06) = X06(I)                                               XSE 1050
         F(I,07) = X07(I)                                               XSE 1060
         F(I,08) = X08(I)                                               XSE 1070
         F(I,09) = X09(I)                                               XSE 1080
         F(I,10) = X10(I)                                               XSE 1090
         F(I,11) = X11(I)                                               XSE 1100
         F(I,12) = X12(I)                                               XSE 1110
         F(I,13) = X13(I)                                               XSE 1120
         F(I,14) = X14(I)                                               XSE 1130
         F(I,15) = X15(I)                                               XSE 1140
         F(I,16) = X16(I)                                               XSE 1150
         F(I,17) = X17(I)                                               XSE 1160
         F(I,18) = X18(I)                                               XSE 1170
         F(I,19) = X19(I)                                               XSE 1180
         F(I,20) = X20(I)                                               XSE 1190
         CBUP(I) = DBUP(I)                                              XSE 1200
   10 CONTINUE                                                          XSE 1210
      DO 30 I=1,20                                                      XSE 1220
         DO 20 J=1,7                                                    XSE 1230
            NFUD(J,I) = MFUD(J,I)                                       XSE 1240
   20    CONTINUE                                                       XSE 1250
   30 CONTINUE                                                          XSE 1260
      RETURN                                                            XSE 1270
      END                                                               XSE 1280
      SUBROUTINE XSEC12(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPUUUC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.120E00,6.144E00,9.080E00,9.080E00,1.095E01,      XSE  120
     *     1.279E01,1.458E01,1.458E01,1.748E01,2.027E01,2.297E01,0.0,19*XSE  130
     *     0.0/                                                         XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /3.888E-01,3.891E-01,3.894E-01,3.896E-01,3.888E-01,      XSE  200
     *     3.891E-01,3.894E-01,3.896E-01,3.888E-01,3.891E-01,3.894E-01, XSE  210
     *     3.896E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /1.067E-02,1.064E-02,1.061E-02,1.058E-02,1.067E-02,      XSE  230
     *     1.064E-02,1.061E-02,1.058E-02,1.067E-02,1.064E-02,1.061E-02, XSE  240
     *     1.058E-02,0.0,19*0.0/                                        XSE  250
      DATA X03 /2.578E-01,2.581E-01,2.583E-01,2.584E-01,2.578E-01,      XSE  260
     *     2.581E-01,2.583E-01,2.584E-01,2.578E-01,2.581E-01,2.583E-01, XSE  270
     *     2.584E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.743E00,2.745E00,2.745E00,2.746E00,2.743E00,2.745E00,  XSE  290
     *     2.745E00,2.746E00,2.743E00,2.745E00,2.745E00,2.746E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /5.681E-01,5.686E-01,5.689E-01,5.692E-01,5.681E-01,      XSE  320
     *     5.686E-01,5.689E-01,5.692E-01,5.681E-01,5.686E-01,5.689E-01, XSE  330
     *     5.692E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /3.190E-01,3.182E-01,3.175E-01,3.169E-01,3.190E-01,      XSE  350
     *     3.182E-01,3.175E-01,3.169E-01,3.190E-01,3.182E-01,3.175E-01, XSE  360
     *     3.169E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /5.561E-01,5.567E-01,5.571E-01,5.576E-01,5.561E-01,      XSE  380
     *     5.567E-01,5.571E-01,5.576E-01,5.561E-01,5.567E-01,5.571E-01, XSE  390
     *     5.576E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.895E00,1.896E00,1.897E00,1.897E00,1.895E00,1.896E00,  XSE  410
     *     1.897E00,1.897E00,1.895E00,1.896E00,1.897E00,1.897E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /2.869E-01,2.872E-01,2.875E-01,2.877E-01,2.869E-01,      XSE  440
     *     2.872E-01,2.875E-01,2.877E-01,2.869E-01,2.872E-01,2.875E-01, XSE  450
     *     2.877E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /4.681E-02,4.666E-02,4.654E-02,4.642E-02,4.681E-02,      XSE  470
     *     4.666E-02,4.654E-02,4.642E-02,4.681E-02,4.666E-02,4.654E-02, XSE  480
     *     4.642E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /4.922E-01,4.927E-01,4.932E-01,4.936E-01,4.922E-01,      XSE  500
     *     4.927E-01,4.932E-01,4.936E-01,4.922E-01,4.927E-01,4.932E-01, XSE  510
     *     4.936E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.791E00,1.791E00,1.791E00,1.791E00,1.791E00,1.791E00,  XSE  530
     *     1.791E00,1.791E00,1.791E00,1.791E00,1.791E00,1.791E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /4.997E-01,5.003E-01,5.007E-01,5.011E-01,4.997E-01,      XSE  560
     *     5.003E-01,5.007E-01,5.011E-01,4.997E-01,5.003E-01,5.007E-01, XSE  570
     *     5.011E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /3.855E-01,3.847E-01,3.840E-01,3.833E-01,3.855E-01,      XSE  590
     *     3.847E-01,3.840E-01,3.833E-01,3.855E-01,3.847E-01,3.840E-01, XSE  600
     *     3.833E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /4.611E-01,4.616E-01,4.619E-01,4.622E-01,4.611E-01,      XSE  620
     *     4.616E-01,4.619E-01,4.622E-01,4.611E-01,4.616E-01,4.619E-01, XSE  630
     *     4.622E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /2.462E00,2.463E00,2.464E00,2.464E00,2.462E00,2.463E00,  XSE  650
     *     2.464E00,2.464E00,2.462E00,2.463E00,2.464E00,2.464E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /5.495E-01,5.502E-01,5.506E-01,5.511E-01,5.495E-01,      XSE  680
     *     5.502E-01,5.506E-01,5.511E-01,5.495E-01,5.502E-01,5.506E-01, XSE  690
     *     5.511E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /2.987E-01,2.979E-01,2.973E-01,2.967E-01,2.987E-01,      XSE  710
     *     2.979E-01,2.973E-01,2.967E-01,2.987E-01,2.979E-01,2.973E-01, XSE  720
     *     2.967E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.366E00,1.367E00,1.368E00,1.369E00,1.366E00,1.367E00,  XSE  740
     *     1.368E00,1.369E00,1.366E00,1.367E00,1.368E00,1.369E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /2.484E-01,2.477E-01,2.471E-01,2.466E-01,2.484E-01,      XSE  770
     *     2.477E-01,2.471E-01,2.466E-01,2.484E-01,2.477E-01,2.471E-01, XSE  780
     *     2.466E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC13(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPUUUA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,1.478E-01,3.159E-01,5.043E-01,5.043E-01,5.972E-01, XSE  120
     *     6.981E-01,8.072E-01,8.072E-01,8.800E-01,9.683E-01,1.072E00,  XSE  130
     *     0.0,19*0.0/                                                  XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /6.435E-01,6.348E-01,6.265E-01,6.185E-01,6.435E-01,      XSE  200
     *     6.348E-01,6.265E-01,6.185E-01,6.435E-01,6.348E-01,6.265E-01, XSE  210
     *     6.185E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /5.122E-03,5.283E-03,5.446E-03,5.609E-03,5.122E-03,      XSE  230
     *     5.283E-03,5.446E-03,5.609E-03,5.122E-03,5.283E-03,5.446E-03, XSE  240
     *     5.609E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /3.639E-01,3.604E-01,3.570E-01,3.538E-01,3.639E-01,      XSE  260
     *     3.604E-01,3.570E-01,3.538E-01,3.639E-01,3.604E-01,3.570E-01, XSE  270
     *     3.538E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /3.483E00,3.458E00,3.434E00,3.411E00,3.483E00,3.458E00,  XSE  290
     *     3.434E00,3.411E00,3.483E00,3.458E00,3.434E00,3.411E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /9.274E-01,9.148E-01,9.026E-01,8.910E-01,9.274E-01,      XSE  320
     *     9.148E-01,9.026E-01,8.910E-01,9.274E-01,9.148E-01,9.026E-01, XSE  330
     *     8.910E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /1.904E-01,1.941E-01,1.978E-01,2.015E-01,1.904E-01,      XSE  350
     *     1.941E-01,1.978E-01,2.015E-01,1.904E-01,1.941E-01,1.978E-01, XSE  360
     *     2.015E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /8.727E-01,8.623E-01,8.521E-01,8.424E-01,8.727E-01,      XSE  380
     *     8.623E-01,8.521E-01,8.424E-01,8.727E-01,8.623E-01,8.521E-01, XSE  390
     *     8.424E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /2.508E00,2.487E00,2.468E00,2.449E00,2.508E00,2.487E00,  XSE  410
     *     2.468E00,2.449E00,2.508E00,2.487E00,2.468E00,2.449E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /3.853E-01,3.822E-01,3.791E-01,3.761E-01,3.853E-01,      XSE  440
     *     3.822E-01,3.791E-01,3.761E-01,3.853E-01,3.822E-01,3.791E-01, XSE  450
     *     3.761E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /2.269E-02,2.339E-02,2.411E-02,2.482E-02,2.269E-02,      XSE  470
     *     2.339E-02,2.411E-02,2.482E-02,2.269E-02,2.339E-02,2.411E-02, XSE  480
     *     2.482E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /9.203E-01,9.057E-01,8.916E-01,8.782E-01,9.203E-01,      XSE  500
     *     9.057E-01,8.916E-01,8.782E-01,9.203E-01,9.057E-01,8.916E-01, XSE  510
     *     8.782E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /2.133E00,2.120E00,2.108E00,2.096E00,2.133E00,2.120E00,  XSE  530
     *     2.108E00,2.096E00,2.133E00,2.120E00,2.108E00,2.096E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /9.409E-01,9.245E-01,9.088E-01,8.940E-01,9.409E-01,      XSE  560
     *     9.245E-01,9.088E-01,8.940E-01,9.409E-01,9.245E-01,9.088E-01, XSE  570
     *     8.940E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /2.611E-01,2.647E-01,2.683E-01,2.719E-01,2.611E-01,      XSE  590
     *     2.647E-01,2.683E-01,2.719E-01,2.611E-01,2.647E-01,2.683E-01, XSE  600
     *     2.719E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /7.389E-01,7.299E-01,7.211E-01,7.126E-01,7.389E-01,      XSE  620
     *     7.299E-01,7.211E-01,7.126E-01,7.389E-01,7.299E-01,7.211E-01, XSE  630
     *     7.126E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /3.264E00,3.237E00,3.210E00,3.185E00,3.264E00,3.237E00,  XSE  650
     *     3.210E00,3.185E00,3.264E00,3.237E00,3.210E00,3.185E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /1.043E00,1.026E00,1.009E00,9.937E-01,1.043E00,1.026E00, XSE  680
     *     1.009E00,9.937E-01,1.043E00,1.026E00,1.009E00,9.937E-01,0.0, XSE  690
     *     19*0.0/                                                      XSE  700
      DATA X18 /1.764E-01,1.800E-01,1.835E-01,1.871E-01,1.764E-01,      XSE  710
     *     1.800E-01,1.835E-01,1.871E-01,1.764E-01,1.800E-01,1.835E-01, XSE  720
     *     1.871E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /2.123E00,2.098E00,2.074E00,2.051E00,2.123E00,2.098E00,  XSE  740
     *     2.074E00,2.051E00,2.123E00,2.098E00,2.074E00,2.051E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /1.436E-01,1.467E-01,1.498E-01,1.528E-01,1.436E-01,      XSE  770
     *     1.467E-01,1.498E-01,1.528E-01,1.436E-01,1.467E-01,1.498E-01, XSE  780
     *     1.528E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC14(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPUUUR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,5.304E-03,1.143E-02,1.844E-02,1.844E-02,2.937E-02, XSE  120
     *     4.181E-02,5.583E-02,5.583E-02,7.964E-02,1.064E-01,1.363E-01, XSE  130
     *     1.363E-01,1.880E-01,2.454E-01,3.089E-01,3.089E-01,4.191E-01, XSE  140
     *     5.402E-01,6.720E-01,6.720E-01,9.014E-01,1.148E00,1.412E00,   XSE  150
     *     0.0,7*0.0/                                                   XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /1.047E00,1.032E00,1.016E00,1.002E00,8.590E-01,8.455E-01,XSE  220
     *     8.325E-01,8.201E-01,7.653E-01,7.531E-01,7.416E-01,7.306E-01, XSE  230
     *     6.905E-01,6.800E-01,6.701E-01,6.606E-01,6.158E-01,6.078E-01, XSE  240
     *     6.002E-01,5.929E-01,5.364E-01,5.321E-01,5.279E-01,5.240E-01, XSE  250
     *     0.0,7*0.0/                                                   XSE  260
      DATA X02 /2.156E-03,2.320E-03,2.481E-03,2.639E-03,2.682E-03,      XSE  270
     *     2.869E-03,3.051E-03,3.229E-03,3.360E-03,3.565E-03,3.764E-03, XSE  280
     *     3.959E-03,4.249E-03,4.462E-03,4.670E-03,4.872E-03,5.484E-03, XSE  290
     *     5.680E-03,5.870E-03,6.056E-03,7.378E-03,7.498E-03,7.616E-03, XSE  300
     *     7.733E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /5.178E-01,5.117E-01,5.057E-01,5.000E-01,4.509E-01,      XSE  320
     *     4.454E-01,4.402E-01,4.353E-01,4.155E-01,4.105E-01,4.058E-01, XSE  330
     *     4.013E-01,3.857E-01,3.814E-01,3.773E-01,3.733E-01,3.548E-01, XSE  340
     *     3.514E-01,3.482E-01,3.452E-01,3.203E-01,3.185E-01,3.168E-01, XSE  350
     *     3.151E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /4.574E00,4.530E00,4.488E00,4.447E00,4.080E00,4.042E00,  XSE  370
     *     4.005E00,3.970E00,3.822E00,3.787E00,3.754E00,3.722E00,       XSE  380
     *     3.607E00,3.577E00,3.549E00,3.521E00,3.388E00,3.365E00,       XSE  390
     *     3.343E00,3.322E00,3.151E00,3.138E00,3.126E00,3.115E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /1.538E00,1.514E00,1.491E00,1.469E00,1.247E00,1.227E00,  XSE  420
     *     1.207E00,1.189E00,1.106E00,1.088E00,1.071E00,1.055E00,       XSE  430
     *     9.960E-01,9.807E-01,9.662E-01,9.524E-01,8.875E-01,8.758E-01, XSE  440
     *     8.648E-01,8.543E-01,7.732E-01,7.669E-01,7.608E-01,7.551E-01, XSE  450
     *     0.0,7*0.0/                                                   XSE  460
      DATA X06 /1.066E-01,1.110E-01,1.153E-01,1.195E-01,1.224E-01,      XSE  470
     *     1.272E-01,1.318E-01,1.363E-01,1.415E-01,1.465E-01,1.513E-01, XSE  480
     *     1.560E-01,1.654E-01,1.703E-01,1.750E-01,1.796E-01,1.967E-01, XSE  490
     *     2.009E-01,2.049E-01,2.089E-01,2.405E-01,2.429E-01,2.453E-01, XSE  500
     *     2.476E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /1.328E00,1.310E00,1.292E00,1.275E00,1.123E00,1.107E00,  XSE  520
     *     1.092E00,1.077E00,1.016E00,1.001E00,9.871E-01,9.737E-01,     XSE  530
     *     9.263E-01,9.136E-01,9.015E-01,8.899E-01,8.347E-01,8.249E-01, XSE  540
     *     8.156E-01,8.067E-01,7.346E-01,7.294E-01,7.243E-01,7.194E-01, XSE  550
     *     0.0,7*0.0/                                                   XSE  560
      DATA X08 /3.397E00,3.361E00,3.327E00,3.294E00,2.992E00,2.961E00,  XSE  570
     *     2.931E00,2.902E00,2.780E00,2.752E00,2.725E00,2.699E00,       XSE  580
     *     2.605E00,2.581E00,2.557E00,2.535E00,2.427E00,2.408E00,       XSE  590
     *     2.390E00,2.373E00,2.236E00,2.226E00,2.216E00,2.207E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /4.999E-01,4.950E-01,4.902E-01,4.856E-01,4.573E-01,      XSE  620
     *     4.526E-01,4.482E-01,4.439E-01,4.295E-01,4.252E-01,4.210E-01, XSE  630
     *     4.170E-01,4.030E-01,3.992E-01,3.955E-01,3.919E-01,3.738E-01, XSE  640
     *     3.708E-01,3.680E-01,3.652E-01,3.397E-01,3.382E-01,3.366E-01, XSE  650
     *     3.351E-01,0.0,7*0.0/                                         XSE  660
      DATA X10 /9.638E-03,1.036E-02,1.106E-02,1.175E-02,1.195E-02,      XSE  670
     *     1.276E-02,1.356E-02,1.434E-02,1.492E-02,1.582E-02,1.669E-02, XSE  680
     *     1.754E-02,1.882E-02,1.975E-02,2.066E-02,2.154E-02,2.424E-02, XSE  690
     *     2.509E-02,2.591E-02,2.672E-02,3.252E-02,3.304E-02,3.356E-02, XSE  700
     *     3.406E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /1.556E00,1.530E00,1.506E00,1.482E00,1.257E00,1.235E00,  XSE  720
     *     1.214E00,1.194E00,1.106E00,1.086E00,1.067E00,1.049E00,       XSE  730
     *     9.840E-01,9.669E-01,9.507E-01,9.353E-01,8.614E-01,8.483E-01, XSE  740
     *     8.359E-01,8.241E-01,7.307E-01,7.236E-01,7.168E-01,7.104E-01, XSE  750
     *     0.0,7*0.0/                                                   XSE  760
      DATA X12 /2.778E00,2.754E00,2.730E00,2.707E00,2.441E00,2.420E00,  XSE  770
     *     2.401E00,2.383E00,2.289E00,2.271E00,2.254E00,2.239E00,       XSE  780
     *     2.178E00,2.163E00,2.149E00,2.136E00,2.075E00,2.063E00,       XSE  790
     *     2.053E00,2.043E00,1.973E00,1.967E00,1.962E00,1.956E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.835E00,1.801E00,1.768E00,1.737E00,1.387E00,1.359E00,  XSE  820
     *     1.333E00,1.307E00,1.186E00,1.162E00,1.139E00,1.118E00,       XSE  830
     *     1.039E00,1.019E00,9.999E-01,9.817E-01,8.995E-01,8.842E-01,   XSE  840
     *     8.698E-01,8.563E-01,7.556E-01,7.474E-01,7.396E-01,7.322E-01, XSE  850
     *     0.0,7*0.0/                                                   XSE  860
      DATA X14 /1.881E-01,1.921E-01,1.960E-01,1.998E-01,1.973E-01,      XSE  870
     *     2.017E-01,2.060E-01,2.102E-01,2.131E-01,2.178E-01,2.224E-01, XSE  880
     *     2.269E-01,2.346E-01,2.394E-01,2.440E-01,2.485E-01,2.640E-01, XSE  890
     *     2.682E-01,2.723E-01,2.762E-01,3.071E-01,3.096E-01,3.119E-01, XSE  900
     *     3.143E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /1.126E00,1.110E00,1.096E00,1.081E00,9.546E-01,9.408E-01,XSE  920
     *     9.276E-01,9.150E-01,8.627E-01,8.501E-01,8.381E-01,8.267E-01, XSE  930
     *     7.852E-01,7.743E-01,7.639E-01,7.539E-01,7.056E-01,6.971E-01, XSE  940
     *     6.891E-01,6.814E-01,6.189E-01,6.144E-01,6.100E-01,6.057E-01, XSE  950
     *     0.0,7*0.0/                                                   XSE  960
      DATA X16 /4.493E00,4.444E00,4.398E00,4.353E00,3.924E00,3.881E00,  XSE  970
     *     3.841E00,3.803E00,3.634E00,3.596E00,3.560E00,3.526E00,       XSE  980
     *     3.400E00,3.367E00,3.336E00,3.307E00,3.165E00,3.140E00,       XSE  990
     *     3.116E00,3.094E00,2.915E00,2.901E00,2.888E00,2.876E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /1.797E00,1.767E00,1.738E00,1.711E00,1.446E00,1.420E00,  XSE 1020
     *     1.395E00,1.372E00,1.268E00,1.245E00,1.222E00,1.201E00,       XSE 1030
     *     1.123E00,1.103E00,1.084E00,1.066E00,9.773E-01,9.620E-01,     XSE 1040
     *     9.474E-01,9.336E-01,8.215E-01,8.133E-01,8.054E-01,7.978E-01, XSE 1050
     *     0.0,7*0.0/                                                   XSE 1060
      DATA X18 /9.958E-02,1.036E-01,1.076E-01,1.114E-01,1.139E-01,      XSE 1070
     *     1.183E-01,1.226E-01,1.268E-01,1.312E-01,1.358E-01,1.403E-01, XSE 1080
     *     1.447E-01,1.530E-01,1.576E-01,1.621E-01,1.665E-01,1.818E-01, XSE 1090
     *     1.859E-01,1.898E-01,1.937E-01,2.233E-01,2.256E-01,2.279E-01, XSE 1100
     *     2.302E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /3.217E00,3.173E00,3.131E00,3.090E00,2.727E00,2.688E00,  XSE 1120
     *     2.651E00,2.616E00,2.470E00,2.434E00,2.401E00,2.369E00,       XSE 1130
     *     2.255E00,2.225E00,2.196E00,2.168E00,2.035E00,2.012E00,       XSE 1140
     *     1.989E00,1.968E00,1.794E00,1.782E00,1.769E00,1.758E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /8.483E-02,8.802E-02,9.114E-02,9.419E-02,9.299E-02,      XSE 1170
     *     9.658E-02,1.001E-01,1.035E-01,1.058E-01,1.097E-01,1.135E-01, XSE 1180
     *     1.172E-01,1.231E-01,1.271E-01,1.310E-01,1.347E-01,1.470E-01, XSE 1190
     *     1.506E-01,1.541E-01,1.575E-01,1.829E-01,1.851E-01,1.872E-01, XSE 1200
     *     1.892E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC15(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMORUUUC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.403E00,6.681E00,9.847E00,9.847E00,1.176E01,      XSE  120
     *     1.363E01,1.546E01,1.546E01,1.797E01,2.040E01,2.276E01,0.0,19*XSE  130
     *     0.0/                                                         XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /3.842E-01,3.842E-01,3.842E-01,3.843E-01,3.842E-01,      XSE  200
     *     3.842E-01,3.842E-01,3.843E-01,3.842E-01,3.842E-01,3.842E-01, XSE  210
     *     3.843E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /1.065E-02,1.062E-02,1.060E-02,1.058E-02,1.065E-02,      XSE  230
     *     1.062E-02,1.060E-02,1.058E-02,1.065E-02,1.062E-02,1.060E-02, XSE  240
     *     1.058E-02,0.0,19*0.0/                                        XSE  250
      DATA X03 /2.562E-01,2.564E-01,2.564E-01,2.565E-01,2.562E-01,      XSE  260
     *     2.564E-01,2.564E-01,2.565E-01,2.562E-01,2.564E-01,2.564E-01, XSE  270
     *     2.565E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.728E00,2.729E00,2.729E00,2.729E00,2.728E00,2.729E00,  XSE  290
     *     2.729E00,2.729E00,2.728E00,2.729E00,2.729E00,2.729E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /5.615E-01,5.616E-01,5.616E-01,5.616E-01,5.615E-01,      XSE  320
     *     5.616E-01,5.616E-01,5.616E-01,5.615E-01,5.616E-01,5.616E-01, XSE  330
     *     5.616E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /3.196E-01,3.190E-01,3.184E-01,3.179E-01,3.196E-01,      XSE  350
     *     3.190E-01,3.184E-01,3.179E-01,3.196E-01,3.190E-01,3.184E-01, XSE  360
     *     3.179E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /5.503E-01,5.505E-01,5.507E-01,5.508E-01,5.503E-01,      XSE  380
     *     5.505E-01,5.507E-01,5.508E-01,5.503E-01,5.505E-01,5.507E-01, XSE  390
     *     5.508E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.882E00,1.882E00,1.882E00,1.882E00,1.882E00,1.882E00,  XSE  410
     *     1.882E00,1.882E00,1.882E00,1.882E00,1.882E00,1.882E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /2.853E-01,2.855E-01,2.856E-01,2.857E-01,2.853E-01,      XSE  440
     *     2.855E-01,2.856E-01,2.857E-01,2.853E-01,2.855E-01,2.856E-01, XSE  450
     *     2.857E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /4.670E-02,4.659E-02,4.649E-02,4.640E-02,4.670E-02,      XSE  470
     *     4.659E-02,4.649E-02,4.640E-02,4.670E-02,4.659E-02,4.649E-02, XSE  480
     *     4.640E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /4.839E-01,4.840E-01,4.841E-01,4.841E-01,4.839E-01,      XSE  500
     *     4.840E-01,4.841E-01,4.841E-01,4.839E-01,4.840E-01,4.841E-01, XSE  510
     *     4.841E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.783E00,1.783E00,1.782E00,1.782E00,1.783E00,1.783E00,  XSE  530
     *     1.782E00,1.782E00,1.783E00,1.783E00,1.782E00,1.782E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /4.932E-01,4.933E-01,4.934E-01,4.935E-01,4.932E-01,      XSE  560
     *     4.933E-01,4.934E-01,4.935E-01,4.932E-01,4.933E-01,4.934E-01, XSE  570
     *     4.935E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /3.857E-01,3.850E-01,3.844E-01,3.839E-01,3.857E-01,      XSE  590
     *     3.850E-01,3.844E-01,3.839E-01,3.857E-01,3.850E-01,3.844E-01, XSE  600
     *     3.839E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /4.555E-01,4.556E-01,4.557E-01,4.558E-01,4.555E-01,      XSE  620
     *     4.556E-01,4.557E-01,4.558E-01,4.555E-01,4.556E-01,4.557E-01, XSE  630
     *     4.558E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /2.444E00,2.445E00,2.445E00,2.445E00,2.444E00,2.445E00,  XSE  650
     *     2.445E00,2.445E00,2.444E00,2.445E00,2.445E00,2.445E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /5.384E-01,5.384E-01,5.384E-01,5.384E-01,5.384E-01,      XSE  680
     *     5.384E-01,5.384E-01,5.384E-01,5.384E-01,5.384E-01,5.384E-01, XSE  690
     *     5.384E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /2.991E-01,2.985E-01,2.980E-01,2.975E-01,2.991E-01,      XSE  710
     *     2.985E-01,2.980E-01,2.975E-01,2.991E-01,2.985E-01,2.980E-01, XSE  720
     *     2.975E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.351E00,1.352E00,1.352E00,1.352E00,1.351E00,1.352E00,  XSE  740
     *     1.352E00,1.352E00,1.351E00,1.352E00,1.352E00,1.352E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /2.482E-01,2.477E-01,2.473E-01,2.468E-01,2.482E-01,      XSE  770
     *     2.477E-01,2.473E-01,2.468E-01,2.482E-01,2.477E-01,2.473E-01, XSE  780
     *     2.468E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC16(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMORUUUA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,9.301E-02,2.092E-01,3.481E-01,3.481E-01,4.203E-01, XSE  120
     *     5.007E-01,5.894E-01,5.894E-01,7.419E-01,9.074E-01,1.086E00,  XSE  130
     *     0.0,19*0.0/                                                  XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /6.418E-01,6.328E-01,6.242E-01,6.161E-01,6.418E-01,      XSE  200
     *     6.328E-01,6.242E-01,6.161E-01,6.418E-01,6.328E-01,6.242E-01, XSE  210
     *     6.161E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /5.115E-03,5.281E-03,5.447E-03,5.612E-03,5.115E-03,      XSE  230
     *     5.281E-03,5.447E-03,5.612E-03,5.115E-03,5.281E-03,5.447E-03, XSE  240
     *     5.612E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /3.633E-01,3.597E-01,3.563E-01,3.530E-01,3.633E-01,      XSE  260
     *     3.597E-01,3.563E-01,3.530E-01,3.633E-01,3.597E-01,3.563E-01, XSE  270
     *     3.530E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /3.478E00,3.453E00,3.428E00,3.404E00,3.478E00,3.453E00,  XSE  290
     *     3.428E00,3.404E00,3.478E00,3.453E00,3.428E00,3.404E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /9.250E-01,9.119E-01,8.994E-01,8.876E-01,9.250E-01,      XSE  320
     *     9.119E-01,8.994E-01,8.876E-01,9.250E-01,9.119E-01,8.994E-01, XSE  330
     *     8.876E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /1.906E-01,1.944E-01,1.982E-01,2.020E-01,1.906E-01,      XSE  350
     *     1.944E-01,1.982E-01,2.020E-01,1.906E-01,1.944E-01,1.982E-01, XSE  360
     *     2.020E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /8.708E-01,8.600E-01,8.496E-01,8.396E-01,8.708E-01,      XSE  380
     *     8.600E-01,8.496E-01,8.396E-01,8.708E-01,8.600E-01,8.496E-01, XSE  390
     *     8.396E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /2.503E00,2.482E00,2.462E00,2.443E00,2.503E00,2.482E00,  XSE  410
     *     2.462E00,2.443E00,2.503E00,2.482E00,2.462E00,2.443E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /3.849E-01,3.817E-01,3.785E-01,3.755E-01,3.849E-01,      XSE  440
     *     3.817E-01,3.785E-01,3.755E-01,3.849E-01,3.817E-01,3.785E-01, XSE  450
     *     3.755E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /2.266E-02,2.338E-02,2.411E-02,2.483E-02,2.266E-02,      XSE  470
     *     2.338E-02,2.411E-02,2.483E-02,2.266E-02,2.338E-02,2.411E-02, XSE  480
     *     2.483E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /9.173E-01,9.023E-01,8.878E-01,8.741E-01,9.173E-01,      XSE  500
     *     9.023E-01,8.878E-01,8.741E-01,9.173E-01,9.023E-01,8.878E-01, XSE  510
     *     8.741E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /2.130E00,2.117E00,2.104E00,2.092E00,2.130E00,2.117E00,  XSE  530
     *     2.104E00,2.092E00,2.130E00,2.117E00,2.104E00,2.092E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /9.380E-01,9.211E-01,9.049E-01,8.899E-01,9.380E-01,      XSE  560
     *     9.211E-01,9.049E-01,8.899E-01,9.380E-01,9.211E-01,9.049E-01, XSE  570
     *     8.899E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /2.612E-01,2.649E-01,2.686E-01,2.722E-01,2.612E-01,      XSE  590
     *     2.649E-01,2.686E-01,2.722E-01,2.612E-01,2.649E-01,2.686E-01, XSE  600
     *     2.722E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /7.371E-01,7.277E-01,7.187E-01,7.101E-01,7.371E-01,      XSE  620
     *     7.277E-01,7.187E-01,7.101E-01,7.371E-01,7.277E-01,7.187E-01, XSE  630
     *     7.101E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /3.258E00,3.230E00,3.203E00,3.177E00,3.258E00,3.230E00,  XSE  650
     *     3.203E00,3.177E00,3.258E00,3.230E00,3.203E00,3.177E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /1.039E00,1.022E00,1.005E00,9.890E-01,1.039E00,1.022E00, XSE  680
     *     1.005E00,9.890E-01,1.039E00,1.022E00,1.005E00,9.890E-01,0.0, XSE  690
     *     19*0.0/                                                      XSE  700
      DATA X18 /1.766E-01,1.802E-01,1.839E-01,1.875E-01,1.766E-01,      XSE  710
     *     1.802E-01,1.839E-01,1.875E-01,1.766E-01,1.802E-01,1.839E-01, XSE  720
     *     1.875E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /2.119E00,2.093E00,2.068E00,2.044E00,2.119E00,2.093E00,  XSE  740
     *     2.068E00,2.044E00,2.119E00,2.093E00,2.068E00,2.044E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /1.436E-01,1.467E-01,1.499E-01,1.530E-01,1.436E-01,      XSE  770
     *     1.467E-01,1.499E-01,1.530E-01,1.436E-01,1.467E-01,1.499E-01, XSE  780
     *     1.530E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC17(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMORUUUR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,5.246E-03,1.127E-02,1.812E-02,1.812E-02,2.907E-02, XSE  120
     *     4.144E-02,5.531E-02,5.531E-02,7.892E-02,1.053E-01,1.347E-01, XSE  130
     *     1.347E-01,1.858E-01,2.424E-01,3.046E-01,3.046E-01,4.141E-01, XSE  140
     *     5.338E-01,6.635E-01,6.635E-01,8.906E-01,1.134E00,1.394E00,   XSE  150
     *     0.0,7*0.0/                                                   XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /1.047E00,1.032E00,1.017E00,1.003E00,8.596E-01,8.461E-01,XSE  220
     *     8.334E-01,8.213E-01,7.666E-01,7.545E-01,7.431E-01,7.323E-01, XSE  230
     *     6.925E-01,6.820E-01,6.721E-01,6.627E-01,6.186E-01,6.104E-01, XSE  240
     *     6.028E-01,5.955E-01,5.402E-01,5.356E-01,5.313E-01,5.271E-01, XSE  250
     *     0.0,7*0.0/                                                   XSE  260
      DATA X02 /2.169E-03,2.328E-03,2.483E-03,2.635E-03,2.698E-03,      XSE  270
     *     2.879E-03,3.055E-03,3.226E-03,3.375E-03,3.574E-03,3.766E-03, XSE  280
     *     3.953E-03,4.264E-03,4.471E-03,4.671E-03,4.866E-03,5.498E-03, XSE  290
     *     5.689E-03,5.873E-03,6.052E-03,7.380E-03,7.500E-03,7.616E-03, XSE  300
     *     7.729E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /5.178E-01,5.117E-01,5.059E-01,5.004E-01,4.510E-01,      XSE  320
     *     4.456E-01,4.405E-01,4.357E-01,4.159E-01,4.110E-01,4.063E-01, XSE  330
     *     4.019E-01,3.864E-01,3.821E-01,3.780E-01,3.741E-01,3.558E-01, XSE  340
     *     3.524E-01,3.492E-01,3.461E-01,3.217E-01,3.198E-01,3.180E-01, XSE  350
     *     3.163E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /4.574E00,4.530E00,4.489E00,4.450E00,4.082E00,4.044E00,  XSE  370
     *     4.008E00,3.973E00,3.825E00,3.791E00,3.758E00,3.727E00,       XSE  380
     *     3.613E00,3.583E00,3.554E00,3.527E00,3.396E00,3.373E00,       XSE  390
     *     3.350E00,3.329E00,3.162E00,3.149E00,3.136E00,3.124E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /1.538E00,1.514E00,1.492E00,1.471E00,1.247E00,1.228E00,  XSE  420
     *     1.209E00,1.191E00,1.108E00,1.090E00,1.073E00,1.057E00,       XSE  430
     *     9.991E-01,9.837E-01,9.692E-01,9.555E-01,8.916E-01,8.797E-01, XSE  440
     *     8.686E-01,8.581E-01,7.787E-01,7.720E-01,7.657E-01,7.597E-01, XSE  450
     *     0.0,7*0.0/                                                   XSE  460
      DATA X06 /1.069E-01,1.112E-01,1.154E-01,1.194E-01,1.228E-01,      XSE  470
     *     1.274E-01,1.319E-01,1.362E-01,1.418E-01,1.467E-01,1.513E-01, XSE  480
     *     1.558E-01,1.657E-01,1.704E-01,1.750E-01,1.794E-01,1.968E-01, XSE  490
     *     2.009E-01,2.049E-01,2.087E-01,2.402E-01,2.426E-01,2.450E-01, XSE  500
     *     2.472E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /1.328E00,1.310E00,1.293E00,1.276E00,1.124E00,1.108E00,  XSE  520
     *     1.093E00,1.078E00,1.017E00,1.002E00,9.887E-01,9.757E-01,     XSE  530
     *     9.285E-01,9.158E-01,9.038E-01,8.924E-01,8.378E-01,8.279E-01, XSE  540
     *     8.185E-01,8.096E-01,7.391E-01,7.335E-01,7.282E-01,7.231E-01, XSE  550
     *     0.0,7*0.0/                                                   XSE  560
      DATA X08 /3.397E00,3.362E00,3.329E00,3.296E00,2.993E00,2.962E00,  XSE  570
     *     2.933E00,2.905E00,2.783E00,2.755E00,2.728E00,2.703E00,       XSE  580
     *     2.610E00,2.585E00,2.562E00,2.540E00,2.434E00,2.415E00,       XSE  590
     *     2.397E00,2.379E00,2.245E00,2.234E00,2.224E00,2.214E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /4.998E-01,4.950E-01,4.904E-01,4.859E-01,4.573E-01,      XSE  620
     *     4.527E-01,4.484E-01,4.442E-01,4.298E-01,4.255E-01,4.214E-01, XSE  630
     *     4.175E-01,4.036E-01,3.997E-01,3.961E-01,3.926E-01,3.747E-01, XSE  640
     *     3.716E-01,3.688E-01,3.660E-01,3.411E-01,3.394E-01,3.378E-01, XSE  650
     *     3.362E-01,0.0,7*0.0/                                         XSE  660
      DATA X10 /9.695E-03,1.039E-02,1.107E-02,1.174E-02,1.202E-02,      XSE  670
     *     1.281E-02,1.358E-02,1.433E-02,1.499E-02,1.586E-02,1.670E-02, XSE  680
     *     1.751E-02,1.889E-02,1.979E-02,2.066E-02,2.151E-02,2.429E-02, XSE  690
     *     2.513E-02,2.593E-02,2.671E-02,3.252E-02,3.305E-02,3.355E-02, XSE  700
     *     3.404E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /1.556E00,1.531E00,1.507E00,1.484E00,1.258E00,1.236E00,  XSE  720
     *     1.215E00,1.196E00,1.108E00,1.088E00,1.070E00,1.052E00,       XSE  730
     *     9.873E-01,9.702E-01,9.540E-01,9.387E-01,8.659E-01,8.526E-01, XSE  740
     *     8.401E-01,8.283E-01,7.369E-01,7.294E-01,7.223E-01,7.156E-01, XSE  750
     *     0.0,7*0.0/                                                   XSE  760
      DATA X12 /2.779E00,2.754E00,2.731E00,2.708E00,2.442E00,2.422E00,  XSE  770
     *     2.403E00,2.384E00,2.291E00,2.273E00,2.257E00,2.241E00,       XSE  780
     *     2.181E00,2.167E00,2.153E00,2.140E00,2.079E00,2.068E00,       XSE  790
     *     2.057E00,2.047E00,1.979E00,1.973E00,1.967E00,1.961E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.835E00,1.802E00,1.770E00,1.739E00,1.388E00,1.361E00,  XSE  820
     *     1.335E00,1.310E00,1.189E00,1.165E00,1.143E00,1.121E00,       XSE  830
     *     1.044E00,1.023E00,1.004E00,9.860E-01,9.049E-01,8.893E-01,    XSE  840
     *     8.748E-01,8.611E-01,7.625E-01,7.538E-01,7.456E-01,7.380E-01, XSE  850
     *     0.0,7*0.0/                                                   XSE  860
      DATA X14 /1.885E-01,1.923E-01,1.961E-01,1.997E-01,1.977E-01,      XSE  870
     *     2.020E-01,2.061E-01,2.102E-01,2.135E-01,2.181E-01,2.225E-01, XSE  880
     *     2.268E-01,2.350E-01,2.396E-01,2.440E-01,2.484E-01,2.643E-01, XSE  890
     *     2.684E-01,2.723E-01,2.762E-01,3.071E-01,3.095E-01,3.118E-01, XSE  900
     *     3.142E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /1.126E00,1.111E00,1.096E00,1.082E00,9.551E-01,9.414E-01,XSE  920
     *     9.285E-01,9.161E-01,8.639E-01,8.514E-01,8.396E-01,8.284E-01, XSE  930
     *     7.872E-01,7.763E-01,7.659E-01,7.561E-01,7.084E-01,6.998E-01, XSE  940
     *     6.917E-01,6.841E-01,6.229E-01,6.181E-01,6.135E-01,6.091E-01, XSE  950
     *     0.0,7*0.0/                                                   XSE  960
      DATA X16 /4.493E00,4.445E00,4.400E00,4.356E00,3.925E00,3.883E00,  XSE  970
     *     3.844E00,3.806E00,3.638E00,3.601E00,3.565E00,3.531E00,       XSE  980
     *     3.407E00,3.374E00,3.343E00,3.314E00,3.174E00,3.148E00,       XSE  990
     *     3.124E00,3.102E00,2.927E00,2.912E00,2.899E00,2.886E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /1.797E00,1.768E00,1.740E00,1.713E00,1.447E00,1.422E00,  XSE 1020
     *     1.397E00,1.374E00,1.270E00,1.247E00,1.225E00,1.204E00,       XSE 1030
     *     1.127E00,1.107E00,1.088E00,1.070E00,9.827E-01,9.671E-01,     XSE 1040
     *     9.525E-01,9.386E-01,8.290E-01,8.202E-01,8.120E-01,8.041E-01, XSE 1050
     *     0.0,7*0.0/                                                   XSE 1060
      DATA X18 /9.988E-02,1.038E-01,1.076E-01,1.113E-01,1.143E-01,      XSE 1070
     *     1.185E-01,1.227E-01,1.267E-01,1.315E-01,1.360E-01,1.403E-01, XSE 1080
     *     1.446E-01,1.532E-01,1.577E-01,1.621E-01,1.663E-01,1.820E-01, XSE 1090
     *     1.860E-01,1.898E-01,1.935E-01,2.231E-01,2.254E-01,2.277E-01, XSE 1100
     *     2.300E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /3.217E00,3.174E00,3.132E00,3.093E00,2.728E00,2.690E00,  XSE 1120
     *     2.654E00,2.619E00,2.473E00,2.438E00,2.405E00,2.374E00,       XSE 1130
     *     2.261E00,2.230E00,2.201E00,2.174E00,2.043E00,2.019E00,       XSE 1140
     *     1.996E00,1.975E00,1.805E00,1.792E00,1.779E00,1.766E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /8.510E-02,8.819E-02,9.119E-02,9.412E-02,9.332E-02,      XSE 1170
     *     9.680E-02,1.002E-01,1.035E-01,1.061E-01,1.099E-01,1.135E-01, XSE 1180
     *     1.171E-01,1.234E-01,1.273E-01,1.310E-01,1.346E-01,1.473E-01, XSE 1190
     *     1.508E-01,1.541E-01,1.574E-01,1.829E-01,1.851E-01,1.871E-01, XSE 1200
     *     1.892E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC18(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    EMOPUUUC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,8/                                            XSE  110
      DATA DBUP /0.0,2.120E00,4.162E00,6.136E00,6.136E00,9.042E00,      XSE  120
     *     1.185E01,1.457E01,23*0.0,0.0/                                XSE  130
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  140
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  150
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  160
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  170
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  180
      DATA X01 /3.988E-01,3.986E-01,3.982E-01,3.980E-01,3.988E-01,      XSE  190
     *     3.986E-01,3.982E-01,3.980E-01,23*0.0,0.0/                    XSE  200
      DATA X02 /1.035E-02,1.033E-02,1.032E-02,1.030E-02,1.035E-02,      XSE  210
     *     1.033E-02,1.032E-02,1.030E-02,23*0.0,0.0/                    XSE  220
      DATA X03 /2.626E-01,2.625E-01,2.625E-01,2.624E-01,2.626E-01,      XSE  230
     *     2.625E-01,2.625E-01,2.624E-01,23*0.0,0.0/                    XSE  240
      DATA X04 /2.777E00,2.777E00,2.776E00,2.775E00,2.777E00,2.777E00,  XSE  250
     *     2.776E00,2.775E00,23*0.0,0.0/                                XSE  260
      DATA X05 /5.824E-01,5.820E-01,5.816E-01,5.812E-01,5.824E-01,      XSE  270
     *     5.820E-01,5.816E-01,5.812E-01,23*0.0,0.0/                    XSE  280
      DATA X06 /3.117E-01,3.113E-01,3.109E-01,3.106E-01,3.117E-01,      XSE  290
     *     3.113E-01,3.109E-01,3.106E-01,23*0.0,0.0/                    XSE  300
      DATA X07 /5.701E-01,5.699E-01,5.696E-01,5.694E-01,5.701E-01,      XSE  310
     *     5.699E-01,5.696E-01,5.694E-01,23*0.0,0.0/                    XSE  320
      DATA X08 /1.922E00,1.921E00,1.921E00,1.920E00,1.922E00,1.921E00,  XSE  330
     *     1.921E00,1.920E00,23*0.0,0.0/                                XSE  340
      DATA X09 /2.929E-01,2.929E-01,2.928E-01,2.928E-01,2.929E-01,      XSE  350
     *     2.929E-01,2.928E-01,2.928E-01,23*0.0,0.0/                    XSE  360
      DATA X10 /4.542E-02,4.534E-02,4.527E-02,4.521E-02,4.542E-02,      XSE  370
     *     4.534E-02,4.527E-02,4.521E-02,23*0.0,0.0/                    XSE  380
      DATA X11 /5.093E-01,5.089E-01,5.084E-01,5.080E-01,5.093E-01,      XSE  390
     *     5.089E-01,5.084E-01,5.080E-01,23*0.0,0.0/                    XSE  400
      DATA X12 /1.802E00,1.801E00,1.801E00,1.800E00,1.802E00,1.801E00,  XSE  410
     *     1.801E00,1.800E00,23*0.0,0.0/                                XSE  420
      DATA X13 /5.151E-01,5.148E-01,5.144E-01,5.140E-01,5.151E-01,      XSE  430
     *     5.148E-01,5.144E-01,5.140E-01,23*0.0,0.0/                    XSE  440
      DATA X14 /3.783E-01,3.779E-01,3.775E-01,3.771E-01,3.783E-01,      XSE  450
     *     3.779E-01,3.775E-01,3.771E-01,23*0.0,0.0/                    XSE  460
      DATA X15 /4.734E-01,4.731E-01,4.728E-01,4.725E-01,4.734E-01,      XSE  470
     *     4.731E-01,4.728E-01,4.725E-01,23*0.0,0.0/                    XSE  480
      DATA X16 /2.496E00,2.495E00,2.494E00,2.493E00,2.496E00,2.495E00,  XSE  490
     *     2.494E00,2.493E00,23*0.0,0.0/                                XSE  500
      DATA X17 /5.720E-01,5.714E-01,5.706E-01,5.699E-01,5.720E-01,      XSE  510
     *     5.714E-01,5.706E-01,5.699E-01,23*0.0,0.0/                    XSE  520
      DATA X18 /2.917E-01,2.913E-01,2.910E-01,2.907E-01,2.917E-01,      XSE  530
     *     2.913E-01,2.910E-01,2.907E-01,23*0.0,0.0/                    XSE  540
      DATA X19 /1.400E00,1.399E00,1.399E00,1.398E00,1.400E00,1.399E00,  XSE  550
     *     1.399E00,1.398E00,23*0.0,0.0/                                XSE  560
      DATA X20 /2.423E-01,2.419E-01,2.416E-01,2.413E-01,2.423E-01,      XSE  570
     *     2.419E-01,2.416E-01,2.413E-01,23*0.0,0.0/                    XSE  580
      NPT = MPT                                                         XSE  590
      NXSEC = MXSEC                                                     XSE  600
      DO 10 I=1,30                                                      XSE  610
         F(I,01) = X01(I)                                               XSE  620
         F(I,02) = X02(I)                                               XSE  630
         F(I,03) = X03(I)                                               XSE  640
         F(I,04) = X04(I)                                               XSE  650
         F(I,05) = X05(I)                                               XSE  660
         F(I,06) = X06(I)                                               XSE  670
         F(I,07) = X07(I)                                               XSE  680
         F(I,08) = X08(I)                                               XSE  690
         F(I,09) = X09(I)                                               XSE  700
         F(I,10) = X10(I)                                               XSE  710
         F(I,11) = X11(I)                                               XSE  720
         F(I,12) = X12(I)                                               XSE  730
         F(I,13) = X13(I)                                               XSE  740
         F(I,14) = X14(I)                                               XSE  750
         F(I,15) = X15(I)                                               XSE  760
         F(I,16) = X16(I)                                               XSE  770
         F(I,17) = X17(I)                                               XSE  780
         F(I,18) = X18(I)                                               XSE  790
         F(I,19) = X19(I)                                               XSE  800
         F(I,20) = X20(I)                                               XSE  810
         CBUP(I) = DBUP(I)                                              XSE  820
   10 CONTINUE                                                          XSE  830
      DO 30 I=1,20                                                      XSE  840
         DO 20 J=1,7                                                    XSE  850
            NFUD(J,I) = MFUD(J,I)                                       XSE  860
   20    CONTINUE                                                       XSE  870
   30 CONTINUE                                                          XSE  880
      RETURN                                                            XSE  890
      END                                                               XSE  900
      SUBROUTINE XSEC19(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    EMOPUUUA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,8/                                            XSE  110
      DATA DBUP /0.0,5.324E-02,1.161E-01,1.885E-01,1.885E-01,3.242E-01, XSE  120
     *     4.775E-01,6.483E-01,23*0.0,0.0/                              XSE  130
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  140
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  150
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  160
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  170
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  180
      DATA X01 /6.584E-01,6.483E-01,6.386E-01,6.296E-01,6.584E-01,      XSE  190
     *     6.483E-01,6.386E-01,6.296E-01,23*0.0,0.0/                    XSE  200
      DATA X02 /4.869E-03,5.049E-03,5.228E-03,5.403E-03,4.869E-03,      XSE  210
     *     5.049E-03,5.228E-03,5.403E-03,23*0.0,0.0/                    XSE  220
      DATA X03 /3.699E-01,3.659E-01,3.620E-01,3.584E-01,3.699E-01,      XSE  230
     *     3.659E-01,3.620E-01,3.584E-01,23*0.0,0.0/                    XSE  240
      DATA X04 /3.526E00,3.498E00,3.470E00,3.444E00,3.526E00,3.498E00,  XSE  250
     *     3.470E00,3.444E00,23*0.0,0.0/                                XSE  260
      DATA X05 /9.490E-01,9.343E-01,9.202E-01,9.070E-01,9.490E-01,      XSE  270
     *     9.343E-01,9.202E-01,9.070E-01,23*0.0,0.0/                    XSE  280
      DATA X06 /1.842E-01,1.883E-01,1.924E-01,1.965E-01,1.842E-01,      XSE  290
     *     1.883E-01,1.924E-01,1.965E-01,23*0.0,0.0/                    XSE  300
      DATA X07 /8.910E-01,8.789E-01,8.672E-01,8.561E-01,8.910E-01,      XSE  310
     *     8.789E-01,8.672E-01,8.561E-01,23*0.0,0.0/                    XSE  320
      DATA X08 /2.543E00,2.520E00,2.497E00,2.475E00,2.543E00,2.520E00,  XSE  330
     *     2.497E00,2.475E00,23*0.0,0.0/                                XSE  340
      DATA X09 /3.909E-01,3.874E-01,3.838E-01,3.805E-01,3.909E-01,      XSE  350
     *     3.874E-01,3.838E-01,3.805E-01,23*0.0,0.0/                    XSE  360
      DATA X10 /2.158E-02,2.237E-02,2.315E-02,2.392E-02,2.158E-02,      XSE  370
     *     2.237E-02,2.315E-02,2.392E-02,23*0.0,0.0/                    XSE  380
      DATA X11 /9.454E-01,9.284E-01,9.122E-01,8.969E-01,9.454E-01,      XSE  390
     *     9.284E-01,9.122E-01,8.969E-01,23*0.0,0.0/                    XSE  400
      DATA X12 /2.155E00,2.140E00,2.125E00,2.112E00,2.155E00,2.140E00,  XSE  410
     *     2.125E00,2.112E00,23*0.0,0.0/                                XSE  420
      DATA X13 /9.683E-01,9.491E-01,9.309E-01,9.140E-01,9.683E-01,      XSE  430
     *     9.491E-01,9.309E-01,9.140E-01,23*0.0,0.0/                    XSE  440
      DATA X14 /2.552E-01,2.592E-01,2.632E-01,2.671E-01,2.552E-01,      XSE  450
     *     2.592E-01,2.632E-01,2.671E-01,23*0.0,0.0/                    XSE  460
      DATA X15 /7.548E-01,7.443E-01,7.342E-01,7.246E-01,7.548E-01,      XSE  470
     *     7.443E-01,7.342E-01,7.246E-01,23*0.0,0.0/                    XSE  480
      DATA X16 /3.311E00,3.279E00,3.249E00,3.220E00,3.311E00,3.279E00,  XSE  490
     *     3.249E00,3.220E00,23*0.0,0.0/                                XSE  500
      DATA X17 /1.072E00,1.052E00,1.033E00,1.016E00,1.072E00,1.052E00,  XSE  510
     *     1.033E00,1.016E00,23*0.0,0.0/                                XSE  520
      DATA X18 /1.705E-01,1.745E-01,1.785E-01,1.824E-01,1.705E-01,      XSE  530
     *     1.745E-01,1.785E-01,1.824E-01,23*0.0,0.0/                    XSE  540
      DATA X19 /2.167E00,2.138E00,2.110E00,2.084E00,2.167E00,2.138E00,  XSE  550
     *     2.110E00,2.084E00,23*0.0,0.0/                                XSE  560
      DATA X20 /1.387E-01,1.421E-01,1.455E-01,1.488E-01,1.387E-01,      XSE  570
     *     1.421E-01,1.455E-01,1.488E-01,23*0.0,0.0/                    XSE  580
      NPT = MPT                                                         XSE  590
      NXSEC = MXSEC                                                     XSE  600
      DO 10 I=1,30                                                      XSE  610
         F(I,01) = X01(I)                                               XSE  620
         F(I,02) = X02(I)                                               XSE  630
         F(I,03) = X03(I)                                               XSE  640
         F(I,04) = X04(I)                                               XSE  650
         F(I,05) = X05(I)                                               XSE  660
         F(I,06) = X06(I)                                               XSE  670
         F(I,07) = X07(I)                                               XSE  680
         F(I,08) = X08(I)                                               XSE  690
         F(I,09) = X09(I)                                               XSE  700
         F(I,10) = X10(I)                                               XSE  710
         F(I,11) = X11(I)                                               XSE  720
         F(I,12) = X12(I)                                               XSE  730
         F(I,13) = X13(I)                                               XSE  740
         F(I,14) = X14(I)                                               XSE  750
         F(I,15) = X15(I)                                               XSE  760
         F(I,16) = X16(I)                                               XSE  770
         F(I,17) = X17(I)                                               XSE  780
         F(I,18) = X18(I)                                               XSE  790
         F(I,19) = X19(I)                                               XSE  800
         F(I,20) = X20(I)                                               XSE  810
         CBUP(I) = DBUP(I)                                              XSE  820
   10 CONTINUE                                                          XSE  830
      DO 30 I=1,20                                                      XSE  840
         DO 20 J=1,7                                                    XSE  850
            NFUD(J,I) = MFUD(J,I)                                       XSE  860
   20    CONTINUE                                                       XSE  870
   30 CONTINUE                                                          XSE  880
      RETURN                                                            XSE  890
      END                                                               XSE  900
      SUBROUTINE XSEC20(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    EMOPUUUR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,3.084E-03,6.396E-03,9.945E-03,9.945E-03,1.867E-02, XSE  120
     *     2.819E-02,3.852E-02,3.852E-02,5.949E-02,8.248E-02,1.075E-01, XSE  130
     *     1.075E-01,1.547E-01,2.062E-01,2.620E-01,2.620E-01,3.629E-01, XSE  140
     *     4.719E-01,5.889E-01,5.889E-01,8.062E-01,1.037E00,1.281E00,   XSE  150
     *     0.0,7*0.0/                                                   XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /7.776E-01,7.668E-01,7.567E-01,7.471E-01,7.777E-01,      XSE  220
     *     7.669E-01,7.567E-01,7.471E-01,7.356E-01,7.252E-01,7.155E-01, XSE  230
     *     7.062E-01,6.744E-01,6.653E-01,6.567E-01,6.486E-01,6.025E-01, XSE  240
     *     5.956E-01,5.891E-01,5.829E-01,5.325E-01,5.288E-01,5.251E-01, XSE  250
     *     5.216E-01,0.0,7*0.0/                                         XSE  260
      DATA X02 /3.151E-03,3.334E-03,3.509E-03,3.676E-03,2.979E-03,      XSE  270
     *     3.153E-03,3.319E-03,3.479E-03,3.449E-03,3.637E-03,3.816E-03, XSE  280
     *     3.989E-03,4.255E-03,4.450E-03,4.637E-03,4.817E-03,5.481E-03, XSE  290
     *     5.656E-03,5.825E-03,5.988E-03,7.409E-03,7.503E-03,7.598E-03, XSE  300
     *     7.691E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /4.204E-01,4.159E-01,4.117E-01,4.077E-01,4.213E-01,      XSE  320
     *     4.168E-01,4.127E-01,4.087E-01,4.047E-01,4.004E-01,3.963E-01, XSE  330
     *     3.925E-01,3.800E-01,3.762E-01,3.726E-01,3.692E-01,3.501E-01, XSE  340
     *     3.472E-01,3.444E-01,3.418E-01,3.188E-01,3.173E-01,3.158E-01, XSE  350
     *     3.144E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /3.855E00,3.823E00,3.794E00,3.766E00,3.860E00,3.829E00,  XSE  370
     *     3.800E00,3.772E00,3.740E00,3.710E00,3.682E00,3.655E00,       XSE  380
     *     3.563E00,3.537E00,3.512E00,3.488E00,3.352E00,3.332E00,       XSE  390
     *     3.313E00,3.295E00,3.140E00,3.129E00,3.119E00,3.109E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /1.125E00,1.109E00,1.094E00,1.080E00,1.124E00,1.108E00,  XSE  420
     *     1.093E00,1.079E00,1.062E00,1.046E00,1.032E00,1.019E00,       XSE  430
     *     9.721E-01,9.588E-01,9.463E-01,9.345E-01,8.677E-01,8.577E-01, XSE  440
     *     8.483E-01,8.394E-01,7.674E-01,7.619E-01,7.567E-01,7.516E-01, XSE  450
     *     0.0,7*0.0/                                                   XSE  460
      DATA X06 /1.363E-01,1.410E-01,1.454E-01,1.496E-01,1.322E-01,      XSE  470
     *     1.366E-01,1.407E-01,1.448E-01,1.449E-01,1.495E-01,1.538E-01, XSE  480
     *     1.579E-01,1.664E-01,1.708E-01,1.750E-01,1.791E-01,1.972E-01, XSE  490
     *     2.009E-01,2.045E-01,2.080E-01,2.413E-01,2.431E-01,2.449E-01, XSE  500
     *     2.468E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /1.030E00,1.017E00,1.004E00,9.926E-01,1.032E00,1.019E00, XSE  520
     *     1.006E00,9.948E-01,9.820E-01,9.694E-01,9.575E-01,9.463E-01,  XSE  530
     *     9.082E-01,8.971E-01,8.866E-01,8.766E-01,8.198E-01,8.114E-01, XSE  540
     *     8.033E-01,7.957E-01,7.301E-01,7.256E-01,7.212E-01,7.170E-01, XSE  550
     *     0.0,7*0.0/                                                   XSE  560
      DATA X08 /2.807E00,2.781E00,2.757E00,2.735E00,2.811E00,2.785E00,  XSE  570
     *     2.761E00,2.739E00,2.713E00,2.688E00,2.665E00,2.644E00,       XSE  580
     *     2.568E00,2.547E00,2.527E00,2.507E00,2.397E00,2.381E00,       XSE  590
     *     2.365E00,2.350E00,2.227E00,2.218E00,2.210E00,2.202E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /4.327E-01,4.287E-01,4.249E-01,4.213E-01,4.355E-01,      XSE  620
     *     4.315E-01,4.278E-01,4.242E-01,4.210E-01,4.172E-01,4.135E-01, XSE  630
     *     4.101E-01,3.986E-01,3.951E-01,3.919E-01,3.887E-01,3.702E-01, XSE  640
     *     3.676E-01,3.651E-01,3.627E-01,3.385E-01,3.372E-01,3.359E-01, XSE  650
     *     3.347E-01,0.0,7*0.0/                                         XSE  660
      DATA X10 /1.400E-02,1.480E-02,1.556E-02,1.629E-02,1.325E-02,      XSE  670
     *     1.401E-02,1.474E-02,1.544E-02,1.532E-02,1.614E-02,1.692E-02, XSE  680
     *     1.768E-02,1.885E-02,1.970E-02,2.052E-02,2.130E-02,2.422E-02, XSE  690
     *     2.499E-02,2.572E-02,2.643E-02,3.265E-02,3.306E-02,3.347E-02, XSE  700
     *     3.388E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /1.126E00,1.108E00,1.091E00,1.076E00,1.127E00,1.109E00,  XSE  720
     *     1.092E00,1.077E00,1.058E00,1.041E00,1.025E00,1.010E00,       XSE  730
     *     9.580E-01,9.431E-01,9.290E-01,9.158E-01,8.399E-01,8.286E-01, XSE  740
     *     8.180E-01,8.079E-01,7.243E-01,7.182E-01,7.123E-01,7.066E-01, XSE  750
     *     0.0,7*0.0/                                                   XSE  760
      DATA X12 /2.309E00,2.294E00,2.279E00,2.266E00,2.303E00,2.287E00,  XSE  770
     *     2.273E00,2.259E00,2.240E00,2.225E00,2.211E00,2.198E00,       XSE  780
     *     2.151E00,2.138E00,2.127E00,2.115E00,2.052E00,2.043E00,       XSE  790
     *     2.034E00,2.026E00,1.967E00,1.962E00,1.957E00,1.952E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.216E00,1.195E00,1.175E00,1.156E00,1.208E00,1.187E00,  XSE  820
     *     1.167E00,1.148E00,1.124E00,1.104E00,1.085E00,1.067E00,       XSE  830
     *     1.007E00,9.891E-01,9.728E-01,9.575E-01,8.722E-01,8.593E-01,  XSE  840
     *     8.472E-01,8.358E-01,7.479E-01,7.408E-01,7.340E-01,7.275E-01, XSE  850
     *     0.0,7*0.0/                                                   XSE  860
      DATA X14 /2.079E-01,2.123E-01,2.165E-01,2.205E-01,2.037E-01,      XSE  870
     *     2.079E-01,2.118E-01,2.156E-01,2.151E-01,2.195E-01,2.236E-01, XSE  880
     *     2.276E-01,2.347E-01,2.390E-01,2.432E-01,2.472E-01,2.638E-01, XSE  890
     *     2.676E-01,2.712E-01,2.747E-01,3.077E-01,3.096E-01,3.114E-01, XSE  900
     *     3.133E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /8.742E-01,8.629E-01,8.522E-01,8.421E-01,8.767E-01,      XSE  920
     *     8.654E-01,8.548E-01,8.447E-01,8.335E-01,8.226E-01,8.124E-01, XSE  930
     *     8.027E-01,7.693E-01,7.597E-01,7.506E-01,7.420E-01,6.925E-01, XSE  940
     *     6.852E-01,6.782E-01,6.716E-01,6.149E-01,6.110E-01,6.072E-01, XSE  950
     *     6.035E-01,0.0,7*0.0/                                         XSE  960
      DATA X16 /3.672E00,3.638E00,3.606E00,3.576E00,3.674E00,3.640E00,  XSE  970
     *     3.608E00,3.578E00,3.542E00,3.510E00,3.479E00,3.450E00,       XSE  980
     *     3.350E00,3.322E00,3.295E00,3.269E00,3.123E00,3.102E00,       XSE  990
     *     3.081E00,3.062E00,2.902E00,2.891E00,2.879E00,2.868E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /1.290E00,1.270E00,1.250E00,1.231E00,1.293E00,1.272E00,  XSE 1020
     *     1.252E00,1.233E00,1.211E00,1.191E00,1.172E00,1.155E00,       XSE 1030
     *     1.092E00,1.075E00,1.058E00,1.043E00,9.520E-01,9.388E-01,     XSE 1040
     *     9.262E-01,9.144E-01,8.140E-01,8.069E-01,8.000E-01,7.934E-01, XSE 1050
     *     0.0,7*0.0/                                                   XSE 1060
      DATA X18 /1.262E-01,1.305E-01,1.346E-01,1.385E-01,1.224E-01,      XSE 1070
     *     1.265E-01,1.304E-01,1.341E-01,1.341E-01,1.383E-01,1.424E-01, XSE 1080
     *     1.463E-01,1.536E-01,1.579E-01,1.619E-01,1.658E-01,1.822E-01, XSE 1090
     *     1.858E-01,1.893E-01,1.927E-01,2.239E-01,2.258E-01,2.276E-01, XSE 1100
     *     2.294E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /2.504E00,2.472E00,2.442E00,2.414E00,2.509E00,2.477E00,  XSE 1120
     *     2.448E00,2.420E00,2.389E00,2.359E00,2.330E00,2.303E00,       XSE 1130
     *     2.212E00,2.185E00,2.160E00,2.136E00,2.000E00,1.979E00,       XSE 1140
     *     1.960E00,1.941E00,1.784E00,1.773E00,1.762E00,1.752E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /1.016E-01,1.052E-01,1.086E-01,1.118E-01,9.826E-02,      XSE 1170
     *     1.016E-01,1.048E-01,1.079E-01,1.074E-01,1.110E-01,1.144E-01, XSE 1180
     *     1.177E-01,1.232E-01,1.268E-01,1.303E-01,1.337E-01,1.469E-01, XSE 1190
     *     1.501E-01,1.532E-01,1.562E-01,1.835E-01,1.851E-01,1.868E-01, XSE 1200
     *     1.884E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC21(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    CANDNATUXSDRNPM/CITATION,XSEC FROM CITATION                        XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,16/                                           XSE  110
      DATA DBUP /0.0,1.302E-01,2.605E-01,3.908E-01,5.213E-01,6.518E-01, XSE  120
     *     7.824E-01,9.130E-01,1.044E00,1.175E00,1.305E00,1.436E00,     XSE  130
     *     1.567E00,1.698E00,1.830E00,1.961E00,16*0.0/                  XSE  140
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  150
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  160
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  170
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  180
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  190
      DATA X01 /3.741E01,3.716E01,3.677E01,3.647E01,3.624E01,3.634E01,  XSE  200
     *     3.619E01,3.608E01,3.598E01,3.590E01,3.610E01,3.604E01,       XSE  210
     *     3.599E01,3.594E01,3.589E01,3.585E01,16*0.0/                  XSE  220
      DATA X02 /3.173E01,3.141E01,3.093E01,3.056E01,3.027E01,3.029E01,  XSE  230
     *     3.011E01,2.997E01,2.984E01,2.974E01,2.991E01,2.984E01,       XSE  240
     *     2.977E01,2.971E01,2.965E01,2.960E01,16*0.0/                  XSE  250
      DATA X03 /1.764E02,1.745E02,1.716E02,1.693E02,1.676E02,1.678E02,  XSE  260
     *     1.668E02,1.659E02,1.651E02,1.645E02,1.656E02,1.652E02,       XSE  270
     *     1.648E02,1.644E02,1.640E02,1.637E02,16*0.0/                  XSE  280
      DATA X04 /6.173E00,6.202E00,6.243E00,6.275E00,6.300E00,6.230E00,  XSE  290
     *     6.245E00,6.258E00,6.270E00,6.279E00,6.196E00,6.203E00,       XSE  300
     *     6.209E00,6.215E00,6.221E00,6.226E00,16*0.0/                  XSE  310
      DATA X05 /1.269E00,1.263E00,1.254E00,1.248E00,1.242E00,1.248E00,  XSE  320
     *     1.245E00,1.242E00,1.240E00,1.238E00,1.246E00,1.245E00,       XSE  330
     *     1.244E00,1.243E00,1.242E00,1.241E00,16*0.0/                  XSE  340
      DATA X06 /6.464E01,6.403E01,6.311E01,6.239E01,6.184E01,6.200E01,  XSE  350
     *     6.166E01,6.138E01,6.114E01,6.094E01,6.137E01,6.122E01,       XSE  360
     *     6.109E01,6.097E01,6.086E01,6.076E01,16*0.0/                  XSE  370
      DATA X07 /1.571E02,1.553E02,1.527E02,1.506E02,1.490E02,1.498E02,  XSE  380
     *     1.488E02,1.480E02,1.473E02,1.468E02,1.482E02,1.478E02,       XSE  390
     *     1.474E02,1.471E02,1.468E02,1.464E02,16*0.0/                  XSE  400
      DATA X08 /5.568E00,5.524E00,5.460E00,5.410E00,5.371E00,5.402E00,  XSE  410
     *     5.378E00,5.358E00,5.341E00,5.327E00,5.374E00,5.364E00,       XSE  420
     *     5.354E00,5.346E00,5.338E00,5.330E00,16*0.0/                  XSE  430
      DATA X09 /1.449E02,1.434E02,1.410E02,1.391E02,1.377E02,1.334E02,  XSE  440
     *     1.326E02,1.319E02,1.313E02,1.308E02,1.281E02,1.278E02,       XSE  450
     *     1.275E02,1.272E02,1.269E02,1.266E02,16*0.0/                  XSE  460
      DATA X10 /3.078E02,3.044E02,2.993E02,2.953E02,2.922E02,2.866E02,  XSE  470
     *     2.847E02,2.832E02,2.819E02,2.808E02,2.780E02,2.772E02,       XSE  480
     *     2.765E02,2.758E02,2.753E02,2.747E02,16*0.0/                  XSE  490
      DATA X11 /2.038E02,2.015E02,1.981E02,1.955E02,1.934E02,1.661E02,  XSE  500
     *     1.651E02,1.642E02,1.634E02,1.628E02,1.504E02,1.499E02,       XSE  510
     *     1.496E02,1.492E02,1.489E02,1.486E02,16*0.0/                  XSE  520
      DATA X12 /1.307E02,1.292E02,1.270E02,1.253E02,1.240E02,1.229E02,  XSE  530
     *     1.221E02,1.214E02,1.209E02,1.204E02,1.202E02,1.199E02,       XSE  540
     *     1.195E02,1.193E02,1.190E02,1.188E02,16*0.0/                  XSE  550
      DATA X13 /3.831E02,3.789E02,3.726E02,3.678E02,3.640E02,3.609E02,  XSE  560
     *     3.586E02,3.567E02,3.551E02,3.538E02,3.533E02,3.523E02,       XSE  570
     *     3.514E02,3.506E02,3.498E02,3.491E02,16*0.0/                  XSE  580
      DATA X14 /2.538E01,2.511E01,2.472E01,2.441E01,2.417E01,2.475E01,  XSE  590
     *     2.460E01,2.448E01,2.438E01,2.429E01,2.484E01,2.477E01,       XSE  600
     *     2.471E01,2.466E01,2.461E01,2.456E01,16*0.0/                  XSE  610
      DATA X15 /2.378E02,2.352E02,2.313E02,2.283E02,2.260E02,2.240E02,  XSE  620
     *     2.226E02,2.214E02,2.204E02,2.196E02,2.191E02,2.185E02,       XSE  630
     *     2.180E02,2.175E02,2.170E02,2.166E02,16*0.0/                  XSE  640
      DATA X16 /2.939E01,2.907E01,2.859E01,2.822E01,2.793E01,2.768E01,  XSE  650
     *     2.751E01,2.737E01,2.724E01,2.714E01,2.708E01,2.701E01,       XSE  660
     *     2.694E01,2.688E01,2.682E01,2.677E01,16*0.0/                  XSE  670
      DATA X17 /3.754E00,3.716E00,3.660E00,3.616E00,3.581E00,3.552E00,  XSE  680
     *     3.531E00,3.514E00,3.499E00,3.487E00,3.480E00,3.471E00,       XSE  690
     *     3.463E00,3.456E00,3.449E00,3.443E00,16*0.0/                  XSE  700
      DATA X18 /7.133E01,7.061E01,6.953E01,6.870E01,6.804E01,6.748E01,  XSE  710
     *     6.709E01,6.676E01,6.649E01,6.626E01,6.612E01,6.595E01,       XSE  720
     *     6.580E01,6.566E01,6.553E01,6.542E01,16*0.0/                  XSE  730
      DATA X19 /1.393E01,1.383E01,1.367E01,1.355E01,1.346E01,1.337E01,  XSE  740
     *     1.331E01,1.327E01,1.323E01,1.319E01,1.317E01,1.315E01,       XSE  750
     *     1.313E01,1.311E01,1.309E01,1.307E01,16*0.0/                  XSE  760
      DATA X20 /1.683E01,1.682E01,1.680E01,1.679E01,1.678E01,1.675E01,  XSE  770
     *     1.674E01,1.674E01,1.673E01,1.673E01,1.671E01,1.671E01,       XSE  780
     *     1.671E01,1.671E01,1.671E01,1.671E01,16*0.0/                  XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC22(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    CANDSEU XSDRNPM/CITATION,XSEC FROM CITATION                        XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,16/                                           XSE  110
      DATA DBUP /0.0,3.655E-01,7.316E-01,1.098E00,1.466E00,1.833E00,    XSE  120
     *     2.202E00,2.571E00,2.940E00,3.310E00,3.681E00,4.052E00,       XSE  130
     *     4.424E00,4.796E00,5.169E00,5.543E00,16*0.0/                  XSE  140
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  150
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  160
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  170
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  180
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  190
      DATA X01 /3.279E01,3.272E01,3.267E01,3.270E01,3.279E01,3.351E01,  XSE  200
     *     3.368E01,3.383E01,3.398E01,3.410E01,3.486E01,3.495E01,       XSE  210
     *     3.501E01,3.504E01,3.505E01,3.505E01,16*0.0/                  XSE  220
      DATA X02 /2.634E01,2.622E01,2.614E01,2.617E01,2.626E01,2.695E01,  XSE  230
     *     2.714E01,2.733E01,2.749E01,2.763E01,2.839E01,2.849E01,       XSE  240
     *     2.856E01,2.859E01,2.860E01,2.859E01,16*0.0/                  XSE  250
      DATA X03 /1.440E02,1.433E02,1.428E02,1.429E02,1.435E02,1.477E02,  XSE  260
     *     1.489E02,1.499E02,1.509E02,1.518E02,1.564E02,1.570E02,       XSE  270
     *     1.574E02,1.576E02,1.576E02,1.576E02,16*0.0/                  XSE  280
      DATA X04 /6.355E00,6.376E00,6.392E00,6.400E00,6.402E00,6.232E00,  XSE  290
     *     6.225E00,6.219E00,6.212E00,6.207E00,6.058E00,6.055E00,       XSE  300
     *     6.054E00,6.054E00,6.056E00,6.059E00,16*0.0/                  XSE  310
      DATA X05 /1.159E00,1.158E00,1.157E00,1.158E00,1.160E00,1.180E00,  XSE  320
     *     1.184E00,1.188E00,1.191E00,1.194E00,1.215E00,1.217E00,       XSE  330
     *     1.218E00,1.219E00,1.219E00,1.219E00,16*0.0/                  XSE  340
      DATA X06 /5.608E01,5.585E01,5.569E01,5.575E01,5.593E01,5.669E01,  XSE  350
     *     5.707E01,5.743E01,5.776E01,5.804E01,5.887E01,5.907E01,       XSE  360
     *     5.920E01,5.927E01,5.928E01,5.926E01,16*0.0/                  XSE  370
      DATA X07 /1.271E02,1.264E02,1.260E02,1.261E02,1.266E02,1.312E02,  XSE  380
     *     1.322E02,1.332E02,1.341E02,1.348E02,1.397E02,1.403E02,       XSE  390
     *     1.406E02,1.408E02,1.408E02,1.408E02,16*0.0/                  XSE  400
      DATA X08 /4.795E00,4.779E00,4.769E00,4.772E00,4.783E00,4.921E00,  XSE  410
     *     4.946E00,4.970E00,4.992E00,5.011E00,5.155E00,5.169E00,       XSE  420
     *     5.177E00,5.182E00,5.183E00,5.181E00,16*0.0/                  XSE  430
      DATA X09 /1.232E02,1.225E02,1.221E02,1.222E02,1.227E02,1.203E02,  XSE  440
     *     1.212E02,1.221E02,1.229E02,1.236E02,1.220E02,1.225E02,       XSE  450
     *     1.228E02,1.229E02,1.230E02,1.229E02,16*0.0/                  XSE  460
      DATA X10 /2.573E02,2.560E02,2.551E02,2.553E02,2.563E02,2.558E02,  XSE  470
     *     2.578E02,2.597E02,2.614E02,2.629E02,2.636E02,2.647E02,       XSE  480
     *     2.654E02,2.657E02,2.658E02,2.657E02,16*0.0/                  XSE  490
      DATA X11 /1.675E02,1.666E02,1.660E02,1.662E02,1.668E02,1.339E02,  XSE  500
     *     1.349E02,1.359E02,1.368E02,1.376E02,1.255E02,1.260E02,       XSE  510
     *     1.263E02,1.264E02,1.265E02,1.264E02,16*0.0/                  XSE  520
      DATA X12 /1.075E02,1.070E02,1.066E02,1.067E02,1.071E02,1.086E02,  XSE  530
     *     1.095E02,1.103E02,1.110E02,1.117E02,1.136E02,1.141E02,       XSE  540
     *     1.143E02,1.145E02,1.145E02,1.145E02,16*0.0/                  XSE  550
      DATA X13 /3.166E02,3.151E02,3.140E02,3.143E02,3.154E02,3.201E02,  XSE  560
     *     3.226E02,3.249E02,3.271E02,3.289E02,3.344E02,3.358E02,       XSE  570
     *     3.366E02,3.370E02,3.371E02,3.370E02,16*0.0/                  XSE  580
      DATA X14 /2.596E01,2.584E01,2.576E01,2.578E01,2.587E01,2.562E01,  XSE  590
     *     2.581E01,2.599E01,2.615E01,2.629E01,2.571E01,2.581E01,       XSE  600
     *     2.587E01,2.590E01,2.591E01,2.590E01,16*0.0/                  XSE  610
      DATA X15 /2.032E02,2.022E02,2.015E02,2.017E02,2.024E02,2.029E02,  XSE  620
     *     2.045E02,2.059E02,2.073E02,2.084E02,2.092E02,2.100E02,       XSE  630
     *     2.106E02,2.108E02,2.109E02,2.108E02,16*0.0/                  XSE  640
      DATA X16 /2.511E01,2.499E01,2.491E01,2.493E01,2.502E01,2.508E01,  XSE  650
     *     2.527E01,2.545E01,2.562E01,2.576E01,2.586E01,2.596E01,       XSE  660
     *     2.602E01,2.606E01,2.606E01,2.605E01,16*0.0/                  XSE  670
      DATA X17 /3.244E00,3.229E00,3.219E00,3.223E00,3.234E00,3.241E00,  XSE  680
     *     3.264E00,3.286E00,3.306E00,3.323E00,3.334E00,3.347E00,       XSE  690
     *     3.354E00,3.358E00,3.359E00,3.358E00,16*0.0/                  XSE  700
      DATA X18 /6.163E01,6.136E01,6.117E01,6.124E01,6.144E01,6.158E01,  XSE  710
     *     6.201E01,6.243E01,6.281E01,6.313E01,6.336E01,6.359E01,       XSE  720
     *     6.373E01,6.381E01,6.383E01,6.380E01,16*0.0/                  XSE  730
      DATA X19 /1.250E01,1.246E01,1.244E01,1.245E01,1.248E01,1.250E01,  XSE  740
     *     1.257E01,1.263E01,1.269E01,1.274E01,1.277E01,1.281E01,       XSE  750
     *     1.283E01,1.284E01,1.285E01,1.284E01,16*0.0/                  XSE  760
      DATA X20 /1.650E01,1.651E01,1.652E01,1.654E01,1.656E01,1.658E01,  XSE  770
     *     1.660E01,1.663E01,1.665E01,1.667E01,1.668E01,1.669E01,       XSE  780
     *     1.670E01,1.671E01,1.671E01,1.672E01,16*0.0/                  XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC23(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO1TTTC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.399E00,6.760E00,1.006E01,1.006E01,1.180E01,      XSE  120
     *     1.360E01,1.543E01,1.543E01,1.754E01,1.968E01,2.185E01,0.0,19*XSE  130
     *     0.0/                                                         XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /3.337E-01,3.355E-01,3.373E-01,3.389E-01,3.337E-01,      XSE  200
     *     3.355E-01,3.373E-01,3.389E-01,3.337E-01,3.355E-01,3.373E-01, XSE  210
     *     3.389E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /1.090E-02,1.080E-02,1.069E-02,1.060E-02,1.090E-02,      XSE  230
     *     1.080E-02,1.069E-02,1.060E-02,1.090E-02,1.080E-02,1.069E-02, XSE  240
     *     1.060E-02,0.0,19*0.0/                                        XSE  250
      DATA X03 /2.351E-01,2.362E-01,2.372E-01,2.382E-01,2.351E-01,      XSE  260
     *     2.362E-01,2.372E-01,2.382E-01,2.351E-01,2.362E-01,2.372E-01, XSE  270
     *     2.382E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.556E00,2.563E00,2.569E00,2.575E00,2.556E00,2.563E00,  XSE  290
     *     2.569E00,2.575E00,2.556E00,2.563E00,2.569E00,2.575E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /4.832E-01,4.857E-01,4.880E-01,4.902E-01,4.832E-01,      XSE  320
     *     4.857E-01,4.880E-01,4.902E-01,4.832E-01,4.857E-01,4.880E-01, XSE  330
     *     4.902E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /3.238E-01,3.212E-01,3.187E-01,3.163E-01,3.238E-01,      XSE  350
     *     3.212E-01,3.187E-01,3.163E-01,3.238E-01,3.212E-01,3.187E-01, XSE  360
     *     3.163E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /4.822E-01,4.851E-01,4.878E-01,4.903E-01,4.822E-01,      XSE  380
     *     4.851E-01,4.878E-01,4.903E-01,4.822E-01,4.851E-01,4.878E-01, XSE  390
     *     4.903E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.738E00,1.743E00,1.748E00,1.752E00,1.738E00,1.743E00,  XSE  410
     *     1.748E00,1.752E00,1.738E00,1.743E00,1.748E00,1.752E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /2.563E-01,2.577E-01,2.591E-01,2.604E-01,2.563E-01,      XSE  440
     *     2.577E-01,2.591E-01,2.604E-01,2.563E-01,2.577E-01,2.591E-01, XSE  450
     *     2.604E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /4.790E-02,4.744E-02,4.699E-02,4.656E-02,4.790E-02,      XSE  470
     *     4.744E-02,4.699E-02,4.656E-02,4.790E-02,4.744E-02,4.699E-02, XSE  480
     *     4.656E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /3.951E-01,3.983E-01,4.012E-01,4.039E-01,3.951E-01,      XSE  500
     *     3.983E-01,4.012E-01,4.039E-01,3.951E-01,3.983E-01,4.012E-01, XSE  510
     *     4.039E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.712E00,1.713E00,1.714E00,1.715E00,1.712E00,1.713E00,  XSE  530
     *     1.714E00,1.715E00,1.712E00,1.713E00,1.714E00,1.715E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /4.204E-01,4.234E-01,4.263E-01,4.289E-01,4.204E-01,      XSE  560
     *     4.234E-01,4.263E-01,4.289E-01,4.204E-01,4.234E-01,4.263E-01, XSE  570
     *     4.289E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /3.856E-01,3.829E-01,3.803E-01,3.779E-01,3.856E-01,      XSE  590
     *     3.829E-01,3.803E-01,3.779E-01,3.856E-01,3.829E-01,3.803E-01, XSE  600
     *     3.779E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /3.919E-01,3.942E-01,3.964E-01,3.984E-01,3.919E-01,      XSE  620
     *     3.942E-01,3.964E-01,3.984E-01,3.919E-01,3.942E-01,3.964E-01, XSE  630
     *     3.984E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /2.258E00,2.264E00,2.270E00,2.275E00,2.258E00,2.264E00,  XSE  650
     *     2.270E00,2.275E00,2.258E00,2.264E00,2.270E00,2.275E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /4.149E-01,4.188E-01,4.224E-01,4.257E-01,4.149E-01,      XSE  680
     *     4.188E-01,4.224E-01,4.257E-01,4.149E-01,4.188E-01,4.224E-01, XSE  690
     *     4.257E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /3.025E-01,3.001E-01,2.977E-01,2.955E-01,3.025E-01,      XSE  710
     *     3.001E-01,2.977E-01,2.955E-01,3.025E-01,3.001E-01,2.977E-01, XSE  720
     *     2.955E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.187E00,1.193E00,1.200E00,1.206E00,1.187E00,1.193E00,  XSE  740
     *     1.200E00,1.206E00,1.187E00,1.193E00,1.200E00,1.206E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /2.484E-01,2.462E-01,2.442E-01,2.422E-01,2.484E-01,      XSE  770
     *     2.462E-01,2.442E-01,2.422E-01,2.484E-01,2.462E-01,2.442E-01, XSE  780
     *     2.422E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC24(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO1TTTA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,6.728E-02,1.526E-01,2.574E-01,2.574E-01,3.066E-01, XSE  120
     *     3.641E-01,4.306E-01,4.306E-01,5.354E-01,6.549E-01,7.899E-01, XSE  130
     *     0.0,19*0.0/                                                  XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /6.317E-01,6.259E-01,6.199E-01,6.139E-01,6.317E-01,      XSE  200
     *     6.259E-01,6.199E-01,6.139E-01,6.317E-01,6.259E-01,6.199E-01, XSE  210
     *     6.139E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /4.937E-03,5.054E-03,5.177E-03,5.309E-03,4.937E-03,      XSE  230
     *     5.054E-03,5.177E-03,5.309E-03,4.937E-03,5.054E-03,5.177E-03, XSE  240
     *     5.309E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /3.602E-01,3.579E-01,3.555E-01,3.531E-01,3.602E-01,      XSE  260
     *     3.579E-01,3.555E-01,3.531E-01,3.602E-01,3.579E-01,3.555E-01, XSE  270
     *     3.531E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /3.450E00,3.433E00,3.417E00,3.400E00,3.450E00,3.433E00,  XSE  290
     *     3.417E00,3.400E00,3.450E00,3.433E00,3.417E00,3.400E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /9.104E-01,9.018E-01,8.931E-01,8.844E-01,9.104E-01,      XSE  320
     *     9.018E-01,8.931E-01,8.844E-01,9.104E-01,9.018E-01,8.931E-01, XSE  330
     *     8.844E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /1.877E-01,1.903E-01,1.930E-01,1.959E-01,1.877E-01,      XSE  350
     *     1.903E-01,1.930E-01,1.959E-01,1.877E-01,1.903E-01,1.930E-01, XSE  360
     *     1.959E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /8.597E-01,8.527E-01,8.456E-01,8.384E-01,8.597E-01,      XSE  380
     *     8.527E-01,8.456E-01,8.384E-01,8.597E-01,8.527E-01,8.456E-01, XSE  390
     *     8.384E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /2.479E00,2.465E00,2.451E00,2.437E00,2.479E00,2.465E00,  XSE  410
     *     2.451E00,2.437E00,2.479E00,2.465E00,2.451E00,2.437E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /3.822E-01,3.803E-01,3.782E-01,3.761E-01,3.822E-01,      XSE  440
     *     3.803E-01,3.782E-01,3.761E-01,3.822E-01,3.803E-01,3.782E-01, XSE  450
     *     3.761E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /2.200E-02,2.252E-02,2.306E-02,2.363E-02,2.200E-02,      XSE  470
     *     2.252E-02,2.306E-02,2.363E-02,2.200E-02,2.252E-02,2.306E-02, XSE  480
     *     2.363E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /9.003E-01,8.904E-01,8.804E-01,8.703E-01,9.003E-01,      XSE  500
     *     8.904E-01,8.804E-01,8.703E-01,9.003E-01,8.904E-01,8.804E-01, XSE  510
     *     8.703E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /2.112E00,2.103E00,2.094E00,2.085E00,2.112E00,2.103E00,  XSE  530
     *     2.094E00,2.085E00,2.112E00,2.103E00,2.094E00,2.085E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /9.207E-01,9.092E-01,8.979E-01,8.866E-01,9.207E-01,      XSE  560
     *     9.092E-01,8.979E-01,8.866E-01,9.207E-01,9.092E-01,8.979E-01, XSE  570
     *     8.866E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /2.575E-01,2.600E-01,2.627E-01,2.656E-01,2.575E-01,      XSE  590
     *     2.600E-01,2.627E-01,2.656E-01,2.575E-01,2.600E-01,2.627E-01, XSE  600
     *     2.656E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /7.266E-01,7.206E-01,7.144E-01,7.082E-01,7.266E-01,      XSE  620
     *     7.206E-01,7.144E-01,7.082E-01,7.266E-01,7.206E-01,7.144E-01, XSE  630
     *     7.082E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /3.225E00,3.207E00,3.188E00,3.170E00,3.225E00,3.207E00,  XSE  650
     *     3.188E00,3.170E00,3.225E00,3.207E00,3.188E00,3.170E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /1.018E00,1.007E00,9.957E-01,9.841E-01,1.018E00,1.007E00,XSE  680
     *     9.957E-01,9.841E-01,1.018E00,1.007E00,9.957E-01,9.841E-01,   XSE  690
     *     0.0,19*0.0/                                                  XSE  700
      DATA X18 /1.738E-01,1.763E-01,1.789E-01,1.818E-01,1.738E-01,      XSE  710
     *     1.763E-01,1.789E-01,1.818E-01,1.738E-01,1.763E-01,1.789E-01, XSE  720
     *     1.818E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /2.092E00,2.075E00,2.058E00,2.041E00,2.092E00,2.075E00,  XSE  740
     *     2.058E00,2.041E00,2.092E00,2.075E00,2.058E00,2.041E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /1.408E-01,1.430E-01,1.454E-01,1.479E-01,1.408E-01,      XSE  770
     *     1.430E-01,1.454E-01,1.479E-01,1.408E-01,1.430E-01,1.454E-01, XSE  780
     *     1.479E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC25(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO1TTTR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,5.703E-04,1.637E-03,3.332E-03,3.332E-03,6.580E-03, XSE  120
     *     1.098E-02,1.674E-02,1.674E-02,2.639E-02,3.870E-02,5.404E-02, XSE  130
     *     5.404E-02,7.902E-02,1.100E-01,1.476E-01,1.476E-01,2.087E-01, XSE  140
     *     2.826E-01,3.704E-01,3.704E-01,5.061E-01,6.667E-01,8.538E-01, XSE  150
     *     0.0,7*0.0/                                                   XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /5.245E-01,5.202E-01,5.161E-01,5.118E-01,4.858E-01,      XSE  220
     *     4.816E-01,4.777E-01,4.736E-01,4.591E-01,4.551E-01,4.514E-01, XSE  230
     *     4.476E-01,4.336E-01,4.300E-01,4.268E-01,4.234E-01,4.052E-01, XSE  240
     *     4.024E-01,3.999E-01,3.974E-01,3.711E-01,3.697E-01,3.685E-01, XSE  250
     *     3.672E-01,0.0,7*0.0/                                         XSE  260
      DATA X02 /1.586E-03,1.737E-03,1.889E-03,2.050E-03,2.129E-03,      XSE  270
     *     2.315E-03,2.491E-03,2.677E-03,2.852E-03,3.065E-03,3.261E-03, XSE  280
     *     3.467E-03,3.802E-03,4.027E-03,4.234E-03,4.450E-03,5.110E-03, XSE  290
     *     5.323E-03,5.516E-03,5.716E-03,7.071E-03,7.215E-03,7.338E-03, XSE  300
     *     7.464E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /4.810E-01,4.760E-01,4.711E-01,4.660E-01,4.156E-01,      XSE  320
     *     4.111E-01,4.069E-01,4.026E-01,3.810E-01,3.770E-01,3.732E-01, XSE  330
     *     3.694E-01,3.540E-01,3.505E-01,3.473E-01,3.440E-01,3.266E-01, XSE  340
     *     3.240E-01,3.216E-01,3.191E-01,2.954E-01,2.940E-01,2.927E-01, XSE  350
     *     2.915E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /4.297E00,4.261E00,4.226E00,4.191E00,3.812E00,3.781E00,  XSE  370
     *     3.752E00,3.722E00,3.560E00,3.533E00,3.507E00,3.481E00,       XSE  380
     *     3.367E00,3.343E00,3.322E00,3.300E00,3.174E00,3.156E00,       XSE  390
     *     3.140E00,3.124E00,2.958E00,2.949E00,2.941E00,2.933E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /1.405E00,1.386E00,1.367E00,1.348E00,1.101E00,1.085E00,  XSE  420
     *     1.070E00,1.055E00,9.579E-01,9.444E-01,9.319E-01,9.192E-01,   XSE  430
     *     8.597E-01,8.486E-01,8.383E-01,8.279E-01,7.683E-01,7.601E-01, XSE  440
     *     7.526E-01,7.450E-01,6.697E-01,6.654E-01,6.617E-01,6.580E-01, XSE  450
     *     0.0,7*0.0/                                                   XSE  460
      DATA X06 /9.523E-02,9.967E-02,1.040E-01,1.086E-01,1.125E-01,      XSE  470
     *     1.176E-01,1.223E-01,1.273E-01,1.337E-01,1.392E-01,1.441E-01, XSE  480
     *     1.493E-01,1.598E-01,1.653E-01,1.702E-01,1.754E-01,1.932E-01, XSE  490
     *     1.980E-01,2.023E-01,2.067E-01,2.381E-01,2.411E-01,2.436E-01, XSE  500
     *     2.462E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /1.209E00,1.195E00,1.180E00,1.165E00,1.010E00,9.967E-01, XSE  520
     *     9.845E-01,9.719E-01,9.056E-01,8.940E-01,8.832E-01,8.723E-01, XSE  530
     *     8.255E-01,8.156E-01,8.065E-01,7.972E-01,7.454E-01,7.379E-01, XSE  540
     *     7.311E-01,7.242E-01,6.550E-01,6.511E-01,6.476E-01,6.442E-01, XSE  550
     *     0.0,7*0.0/                                                   XSE  560
      DATA X08 /3.139E00,3.110E00,3.083E00,3.054E00,2.746E00,2.721E00,  XSE  570
     *     2.698E00,2.674E00,2.544E00,2.522E00,2.501E00,2.481E00,       XSE  580
     *     2.389E00,2.370E00,2.353E00,2.336E00,2.235E00,2.221E00,       XSE  590
     *     2.209E00,2.196E00,2.064E00,2.057E00,2.051E00,2.044E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /9.899E-01,9.746E-01,9.596E-01,9.443E-01,7.343E-01,      XSE  620
     *     7.223E-01,7.108E-01,6.992E-01,6.192E-01,6.092E-01,5.999E-01, XSE  630
     *     5.904E-01,5.449E-01,5.367E-01,5.292E-01,5.215E-01,4.784E-01, XSE  640
     *     4.724E-01,4.670E-01,4.615E-01,4.088E-01,4.057E-01,4.029E-01, XSE  650
     *     4.002E-01,0.0,7*0.0/                                         XSE  660
      DATA X10 /7.178E-03,7.845E-03,8.513E-03,9.222E-03,9.575E-03,      XSE  670
     *     1.039E-02,1.117E-02,1.199E-02,1.277E-02,1.370E-02,1.456E-02, XSE  680
     *     1.547E-02,1.695E-02,1.794E-02,1.885E-02,1.979E-02,2.270E-02, XSE  690
     *     2.363E-02,2.448E-02,2.535E-02,3.130E-02,3.192E-02,3.246E-02, XSE  700
     *     3.301E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /1.445E00,1.423E00,1.402E00,1.380E00,1.117E00,1.098E00,  XSE  720
     *     1.081E00,1.064E00,9.577E-01,9.419E-01,9.273E-01,9.124E-01,   XSE  730
     *     8.442E-01,8.310E-01,8.189E-01,8.065E-01,7.359E-01,7.261E-01, XSE  740
     *     7.172E-01,7.082E-01,6.175E-01,6.124E-01,6.079E-01,6.035E-01, XSE  750
     *     0.0,7*0.0/                                                   XSE  760
      DATA X12 /2.621E00,2.601E00,2.582E00,2.563E00,2.270E00,2.255E00,  XSE  770
     *     2.241E00,2.227E00,2.117E00,2.105E00,2.094E00,2.083E00,       XSE  780
     *     2.022E00,2.013E00,2.004E00,1.995E00,1.941E00,1.934E00,       XSE  790
     *     1.928E00,1.922E00,1.859E00,1.856E00,1.853E00,1.850E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.717E00,1.689E00,1.661E00,1.633E00,1.229E00,1.208E00,  XSE  820
     *     1.187E00,1.165E00,1.014E00,9.963E-01,9.796E-01,9.626E-01,    XSE  830
     *     8.796E-01,8.651E-01,8.517E-01,8.382E-01,7.621E-01,7.516E-01, XSE  840
     *     7.420E-01,7.324E-01,6.411E-01,6.356E-01,6.308E-01,6.261E-01, XSE  850
     *     0.0,7*0.0/                                                   XSE  860
      DATA X14 /1.726E-01,1.767E-01,1.807E-01,1.849E-01,1.831E-01,      XSE  870
     *     1.879E-01,1.924E-01,1.972E-01,2.011E-01,2.064E-01,2.113E-01, XSE  880
     *     2.164E-01,2.253E-01,2.307E-01,2.357E-01,2.408E-01,2.573E-01, XSE  890
     *     2.622E-01,2.667E-01,2.712E-01,3.020E-01,3.052E-01,3.079E-01, XSE  900
     *     3.106E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /1.016E00,1.004E00,9.915E-01,9.789E-01,8.482E-01,        XSE  920
     *     8.371E-01,8.267E-01,8.160E-01,7.597E-01,7.498E-01,7.406E-01, XSE  930
     *     7.312E-01,6.906E-01,6.821E-01,6.743E-01,6.664E-01,6.211E-01, XSE  940
     *     6.147E-01,6.089E-01,6.030E-01,5.427E-01,5.393E-01,5.365E-01, XSE  950
     *     5.336E-01,0.0,7*0.0/                                         XSE  960
      DATA X16 /4.163E00,4.123E00,4.086E00,4.047E00,3.593E00,3.560E00,  XSE  970
     *     3.529E00,3.497E00,3.311E00,3.282E00,3.255E00,3.227E00,       XSE  980
     *     3.103E00,3.079E00,3.056E00,3.034E00,2.902E00,2.884E00,       XSE  990
     *     2.868E00,2.851E00,2.681E00,2.672E00,2.664E00,2.655E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /1.660E00,1.634E00,1.609E00,1.584E00,1.283E00,1.261E00,  XSE 1020
     *     1.241E00,1.220E00,1.097E00,1.078E00,1.061E00,1.043E00,       XSE 1030
     *     9.610E-01,9.452E-01,9.307E-01,9.160E-01,8.295E-01,8.179E-01, XSE 1040
     *     8.073E-01,7.966E-01,6.850E-01,6.790E-01,6.738E-01,6.686E-01, XSE 1050
     *     0.0,7*0.0/                                                   XSE 1060
      DATA X18 /9.056E-02,9.453E-02,9.845E-02,1.026E-01,1.061E-01,      XSE 1070
     *     1.107E-01,1.150E-01,1.196E-01,1.251E-01,1.301E-01,1.347E-01, XSE 1080
     *     1.396E-01,1.488E-01,1.539E-01,1.585E-01,1.634E-01,1.795E-01, XSE 1090
     *     1.841E-01,1.883E-01,1.926E-01,2.222E-01,2.252E-01,2.277E-01, XSE 1100
     *     2.303E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /2.932E00,2.897E00,2.862E00,2.827E00,2.454E00,2.422E00,  XSE 1120
     *     2.393E00,2.363E00,2.204E00,2.176E00,2.150E00,2.124E00,       XSE 1130
     *     2.012E00,1.988E00,1.966E00,1.944E00,1.819E00,1.801E00,       XSE 1140
     *     1.785E00,1.768E00,1.601E00,1.592E00,1.583E00,1.575E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /7.221E-02,7.537E-02,7.850E-02,8.182E-02,8.134E-02,      XSE 1170
     *     8.518E-02,8.878E-02,9.257E-02,9.579E-02,1.001E-01,1.040E-01, XSE 1180
     *     1.082E-01,1.151E-01,1.196E-01,1.236E-01,1.279E-01,1.412E-01, XSE 1190
     *     1.453E-01,1.491E-01,1.529E-01,1.790E-01,1.817E-01,1.840E-01, XSE 1200
     *     1.864E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC26(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO2TTTC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.570E00,7.046E00,1.043E01,1.043E01,1.236E01,      XSE  120
     *     1.430E01,1.624E01,1.624E01,1.858E01,2.091E01,2.321E01,0.0,19*XSE  130
     *     0.0/                                                         XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /3.184E-01,3.199E-01,3.214E-01,3.228E-01,3.184E-01,      XSE  200
     *     3.199E-01,3.214E-01,3.228E-01,3.184E-01,3.199E-01,3.214E-01, XSE  210
     *     3.228E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /1.182E-02,1.171E-02,1.161E-02,1.151E-02,1.182E-02,      XSE  230
     *     1.171E-02,1.161E-02,1.151E-02,1.182E-02,1.171E-02,1.161E-02, XSE  240
     *     1.151E-02,0.0,19*0.0/                                        XSE  250
      DATA X03 /2.311E-01,2.322E-01,2.332E-01,2.341E-01,2.311E-01,      XSE  260
     *     2.322E-01,2.332E-01,2.341E-01,2.311E-01,2.322E-01,2.332E-01, XSE  270
     *     2.341E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.558E00,2.564E00,2.570E00,2.576E00,2.558E00,2.564E00,  XSE  290
     *     2.570E00,2.576E00,2.558E00,2.564E00,2.570E00,2.576E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /4.864E-01,4.889E-01,4.912E-01,4.933E-01,4.864E-01,      XSE  320
     *     4.889E-01,4.912E-01,4.933E-01,4.864E-01,4.889E-01,4.912E-01, XSE  330
     *     4.933E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /3.614E-01,3.587E-01,3.562E-01,3.538E-01,3.614E-01,      XSE  350
     *     3.587E-01,3.562E-01,3.538E-01,3.614E-01,3.587E-01,3.562E-01, XSE  360
     *     3.538E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /4.782E-01,4.811E-01,4.837E-01,4.861E-01,4.782E-01,      XSE  380
     *     4.811E-01,4.837E-01,4.861E-01,4.782E-01,4.811E-01,4.837E-01, XSE  390
     *     4.861E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.746E00,1.751E00,1.755E00,1.760E00,1.746E00,1.751E00,  XSE  410
     *     1.755E00,1.760E00,1.746E00,1.751E00,1.755E00,1.760E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /2.577E-01,2.592E-01,2.604E-01,2.617E-01,2.577E-01,      XSE  440
     *     2.592E-01,2.604E-01,2.617E-01,2.577E-01,2.592E-01,2.604E-01, XSE  450
     *     2.617E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /5.226E-02,5.179E-02,5.134E-02,5.092E-02,5.226E-02,      XSE  470
     *     5.179E-02,5.134E-02,5.092E-02,5.226E-02,5.179E-02,5.134E-02, XSE  480
     *     5.092E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /4.012E-01,4.044E-01,4.073E-01,4.101E-01,4.012E-01,      XSE  500
     *     4.044E-01,4.073E-01,4.101E-01,4.012E-01,4.044E-01,4.073E-01, XSE  510
     *     4.101E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.735E00,1.736E00,1.738E00,1.739E00,1.735E00,1.736E00,  XSE  530
     *     1.738E00,1.739E00,1.735E00,1.736E00,1.738E00,1.739E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /4.208E-01,4.237E-01,4.265E-01,4.290E-01,4.208E-01,      XSE  560
     *     4.237E-01,4.265E-01,4.290E-01,4.208E-01,4.237E-01,4.265E-01, XSE  570
     *     4.290E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /4.274E-01,4.246E-01,4.220E-01,4.196E-01,4.274E-01,      XSE  590
     *     4.246E-01,4.220E-01,4.196E-01,4.274E-01,4.246E-01,4.220E-01, XSE  600
     *     4.196E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /3.931E-01,3.954E-01,3.976E-01,3.996E-01,3.931E-01,      XSE  620
     *     3.954E-01,3.976E-01,3.996E-01,3.931E-01,3.954E-01,3.976E-01, XSE  630
     *     3.996E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /2.273E00,2.279E00,2.285E00,2.290E00,2.273E00,2.279E00,  XSE  650
     *     2.285E00,2.290E00,2.273E00,2.279E00,2.285E00,2.290E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /4.375E-01,4.415E-01,4.452E-01,4.487E-01,4.375E-01,      XSE  680
     *     4.415E-01,4.452E-01,4.487E-01,4.375E-01,4.415E-01,4.452E-01, XSE  690
     *     4.487E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /3.368E-01,3.342E-01,3.319E-01,3.296E-01,3.368E-01,      XSE  710
     *     3.342E-01,3.319E-01,3.296E-01,3.368E-01,3.342E-01,3.319E-01, XSE  720
     *     3.296E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.178E00,1.184E00,1.191E00,1.197E00,1.178E00,1.184E00,  XSE  740
     *     1.191E00,1.197E00,1.178E00,1.184E00,1.191E00,1.197E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /2.802E-01,2.780E-01,2.758E-01,2.739E-01,2.802E-01,      XSE  770
     *     2.780E-01,2.758E-01,2.739E-01,2.802E-01,2.780E-01,2.758E-01, XSE  780
     *     2.739E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC27(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO2TTTA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,2.100E-02,6.160E-02,1.238E-01,1.238E-01,1.548E-01, XSE  120
     *     1.935E-01,2.407E-01,2.407E-01,3.310E-01,4.354E-01,5.550E-01, XSE  130
     *     0.0,19*0.0/                                                  XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /4.329E-01,4.307E-01,4.283E-01,4.257E-01,4.329E-01,      XSE  200
     *     4.307E-01,4.283E-01,4.257E-01,4.329E-01,4.307E-01,4.283E-01, XSE  210
     *     4.257E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /5.259E-03,5.380E-03,5.521E-03,5.670E-03,5.259E-03,      XSE  230
     *     5.380E-03,5.521E-03,5.670E-03,5.259E-03,5.380E-03,5.521E-03, XSE  240
     *     5.670E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /3.378E-01,3.358E-01,3.335E-01,3.312E-01,3.378E-01,      XSE  260
     *     3.358E-01,3.335E-01,3.312E-01,3.378E-01,3.358E-01,3.335E-01, XSE  270
     *     3.312E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /3.289E00,3.274E00,3.259E00,3.243E00,3.289E00,3.274E00,  XSE  290
     *     3.259E00,3.243E00,3.289E00,3.274E00,3.259E00,3.243E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /8.358E-01,8.285E-01,8.208E-01,8.130E-01,8.358E-01,      XSE  320
     *     8.285E-01,8.208E-01,8.130E-01,8.358E-01,8.285E-01,8.208E-01, XSE  330
     *     8.130E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /2.097E-01,2.124E-01,2.156E-01,2.190E-01,2.097E-01,      XSE  350
     *     2.124E-01,2.156E-01,2.190E-01,2.097E-01,2.124E-01,2.156E-01, XSE  360
     *     2.190E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /7.905E-01,7.845E-01,7.779E-01,7.711E-01,7.905E-01,      XSE  380
     *     7.845E-01,7.779E-01,7.711E-01,7.905E-01,7.845E-01,7.779E-01, XSE  390
     *     7.711E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /2.342E00,2.331E00,2.318E00,2.305E00,2.342E00,2.331E00,  XSE  410
     *     2.318E00,2.305E00,2.342E00,2.331E00,2.318E00,2.305E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /5.255E-01,5.201E-01,5.144E-01,5.086E-01,5.255E-01,      XSE  440
     *     5.201E-01,5.144E-01,5.086E-01,5.255E-01,5.201E-01,5.144E-01, XSE  450
     *     5.086E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /2.354E-02,2.407E-02,2.469E-02,2.534E-02,2.354E-02,      XSE  470
     *     2.407E-02,2.469E-02,2.534E-02,2.354E-02,2.407E-02,2.469E-02, XSE  480
     *     2.534E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /8.151E-01,8.065E-01,7.974E-01,7.881E-01,8.151E-01,      XSE  500
     *     8.065E-01,7.974E-01,7.881E-01,8.151E-01,8.065E-01,7.974E-01, XSE  510
     *     7.881E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /2.039E00,2.032E00,2.025E00,2.017E00,2.039E00,2.032E00,  XSE  530
     *     2.025E00,2.017E00,2.039E00,2.032E00,2.025E00,2.017E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /8.385E-01,8.291E-01,8.191E-01,8.092E-01,8.385E-01,      XSE  560
     *     8.291E-01,8.191E-01,8.092E-01,8.385E-01,8.291E-01,8.191E-01, XSE  570
     *     8.092E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /2.772E-01,2.800E-01,2.832E-01,2.867E-01,2.772E-01,      XSE  590
     *     2.800E-01,2.832E-01,2.867E-01,2.772E-01,2.800E-01,2.832E-01, XSE  600
     *     2.867E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /6.657E-01,6.605E-01,6.548E-01,6.490E-01,6.657E-01,      XSE  620
     *     6.605E-01,6.548E-01,6.490E-01,6.657E-01,6.605E-01,6.548E-01, XSE  630
     *     6.490E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /3.050E00,3.034E00,3.018E00,3.001E00,3.050E00,3.034E00,  XSE  650
     *     3.018E00,3.001E00,3.050E00,3.034E00,3.018E00,3.001E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /9.197E-01,9.098E-01,8.993E-01,8.885E-01,9.197E-01,      XSE  680
     *     9.098E-01,8.993E-01,8.885E-01,9.197E-01,9.098E-01,8.993E-01, XSE  690
     *     8.885E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /1.925E-01,1.952E-01,1.983E-01,2.016E-01,1.925E-01,      XSE  710
     *     1.952E-01,1.983E-01,2.016E-01,1.925E-01,1.952E-01,1.983E-01, XSE  720
     *     2.016E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.929E00,1.915E00,1.899E00,1.882E00,1.929E00,1.915E00,  XSE  740
     *     1.899E00,1.882E00,1.929E00,1.915E00,1.899E00,1.882E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /1.541E-01,1.565E-01,1.592E-01,1.622E-01,1.541E-01,      XSE  770
     *     1.565E-01,1.592E-01,1.622E-01,1.541E-01,1.565E-01,1.592E-01, XSE  780
     *     1.622E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC28(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO2TTTR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,5.513E-04,1.558E-03,3.109E-03,3.109E-03,6.210E-03, XSE  120
     *     1.031E-02,1.554E-02,1.554E-02,2.482E-02,3.633E-02,5.034E-02, XSE  130
     *     5.034E-02,7.445E-02,1.035E-01,1.379E-01,1.379E-01,1.966E-01, XSE  140
     *     2.654E-01,3.453E-01,3.453E-01,4.739E-01,6.217E-01,7.896E-01, XSE  150
     *     0.0,7*0.0/                                                   XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /5.406E-01,5.363E-01,5.321E-01,5.279E-01,5.005E-01,      XSE  220
     *     4.963E-01,4.924E-01,4.884E-01,4.714E-01,4.674E-01,4.637E-01, XSE  230
     *     4.601E-01,4.427E-01,4.392E-01,4.360E-01,4.328E-01,4.105E-01, XSE  240
     *     4.077E-01,4.054E-01,4.030E-01,3.714E-01,3.699E-01,3.689E-01, XSE  250
     *     3.678E-01,0.0,7*0.0/                                         XSE  260
      DATA X02 /1.575E-03,1.718E-03,1.860E-03,2.009E-03,2.126E-03,      XSE  270
     *     2.301E-03,2.465E-03,2.636E-03,2.868E-03,3.067E-03,3.248E-03, XSE  280
     *     3.435E-03,3.864E-03,4.073E-03,4.261E-03,4.454E-03,5.284E-03, XSE  290
     *     5.480E-03,5.649E-03,5.823E-03,7.523E-03,7.656E-03,7.754E-03, XSE  300
     *     7.854E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /5.054E-01,5.001E-01,4.950E-01,4.898E-01,4.363E-01,      XSE  320
     *     4.316E-01,4.272E-01,4.227E-01,3.978E-01,3.936E-01,3.897E-01, XSE  330
     *     3.858E-01,3.666E-01,3.630E-01,3.597E-01,3.564E-01,3.346E-01, XSE  340
     *     3.318E-01,3.294E-01,3.270E-01,2.979E-01,2.964E-01,2.953E-01, XSE  350
     *     2.942E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /4.490E00,4.452E00,4.415E00,4.378E00,3.980E00,3.946E00,  XSE  370
     *     3.915E00,3.884E00,3.699E00,3.669E00,3.642E00,3.615E00,       XSE  380
     *     3.475E00,3.450E00,3.427E00,3.405E00,3.249E00,3.230E00,       XSE  390
     *     3.214E00,3.198E00,2.998E00,2.988E00,2.981E00,2.974E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /1.510E00,1.489E00,1.469E00,1.449E00,1.192E00,1.174E00,  XSE  420
     *     1.158E00,1.142E00,1.033E00,1.018E00,1.005E00,9.915E-01,      XSE  430
     *     9.192E-01,9.070E-01,8.961E-01,8.851E-01,8.112E-01,8.023E-01, XSE  440
     *     7.945E-01,7.867E-01,6.958E-01,6.911E-01,6.875E-01,6.839E-01, XSE  450
     *     0.0,7*0.0/                                                   XSE  460
      DATA X06 /9.842E-02,1.027E-01,1.069E-01,1.112E-01,1.171E-01,      XSE  470
     *     1.220E-01,1.265E-01,1.312E-01,1.403E-01,1.455E-01,1.501E-01, XSE  480
     *     1.550E-01,1.693E-01,1.745E-01,1.790E-01,1.837E-01,2.077E-01, XSE  490
     *     2.122E-01,2.160E-01,2.199E-01,2.622E-01,2.650E-01,2.669E-01, XSE  500
     *     2.690E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /1.291E00,1.275E00,1.260E00,1.244E00,1.079E00,1.065E00,  XSE  520
     *     1.052E00,1.039E00,9.621E-01,9.497E-01,9.384E-01,9.270E-01,   XSE  530
     *     8.691E-01,8.585E-01,8.490E-01,8.395E-01,7.750E-01,7.671E-01, XSE  540
     *     7.602E-01,7.533E-01,6.694E-01,6.652E-01,6.621E-01,6.588E-01, XSE  550
     *     0.0,7*0.0/                                                   XSE  560
      DATA X08 /3.314E00,3.283E00,3.254E00,3.224E00,2.896E00,2.870E00,  XSE  570
     *     2.845E00,2.819E00,2.668E00,2.645E00,2.623E00,2.601E00,       XSE  580
     *     2.487E00,2.467E00,2.449E00,2.431E00,2.306E00,2.291E00,       XSE  590
     *     2.278E00,2.265E00,2.106E00,2.099E00,2.093E00,2.087E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /1.070E00,1.053E00,1.037E00,1.021E00,8.011E-01,7.878E-01,XSE  620
     *     7.754E-01,7.629E-01,6.736E-01,6.625E-01,6.523E-01,6.422E-01, XSE  630
     *     5.870E-01,5.780E-01,5.699E-01,5.618E-01,5.081E-01,5.016E-01, XSE  640
     *     4.959E-01,4.902E-01,4.258E-01,4.224E-01,4.198E-01,4.171E-01, XSE  650
     *     0.0,7*0.0/                                                   XSE  660
      DATA X10 /7.175E-03,7.806E-03,8.437E-03,9.094E-03,9.620E-03,      XSE  670
     *     1.039E-02,1.112E-02,1.187E-02,1.291E-02,1.379E-02,1.459E-02, XSE  680
     *     1.541E-02,1.732E-02,1.825E-02,1.907E-02,1.992E-02,2.360E-02, XSE  690
     *     2.446E-02,2.520E-02,2.596E-02,3.347E-02,3.405E-02,3.447E-02, XSE  700
     *     3.491E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /1.583E00,1.559E00,1.536E00,1.513E00,1.231E00,1.211E00,  XSE  720
     *     1.192E00,1.173E00,1.051E00,1.034E00,1.018E00,1.002E00,       XSE  730
     *     9.177E-01,9.032E-01,8.902E-01,8.772E-01,7.887E-01,7.781E-01, XSE  740
     *     7.688E-01,7.595E-01,6.495E-01,6.440E-01,6.397E-01,6.354E-01, XSE  750
     *     0.0,7*0.0/                                                   XSE  760
      DATA X12 /2.768E00,2.746E00,2.724E00,2.703E00,2.389E00,2.372E00,  XSE  770
     *     2.355E00,2.339E00,2.214E00,2.200E00,2.187E00,2.174E00,       XSE  780
     *     2.100E00,2.089E00,2.079E00,2.069E00,2.000E00,1.992E00,       XSE  790
     *     1.986E00,1.979E00,1.903E00,1.899E00,1.896E00,1.893E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.855E00,1.824E00,1.795E00,1.765E00,1.343E00,1.319E00,  XSE  820
     *     1.296E00,1.274E00,1.106E00,1.086E00,1.068E00,1.050E00,       XSE  830
     *     9.497E-01,9.337E-01,9.193E-01,9.049E-01,8.105E-01,7.991E-01, XSE  840
     *     7.890E-01,7.790E-01,6.674E-01,6.616E-01,6.569E-01,6.522E-01, XSE  850
     *     0.0,7*0.0/                                                   XSE  860
      DATA X14 /1.791E-01,1.830E-01,1.868E-01,1.907E-01,1.906E-01,      XSE  870
     *     1.953E-01,1.995E-01,2.040E-01,2.103E-01,2.153E-01,2.199E-01, XSE  880
     *     2.246E-01,2.370E-01,2.422E-01,2.467E-01,2.514E-01,2.739E-01, XSE  890
     *     2.785E-01,2.824E-01,2.864E-01,3.285E-01,3.314E-01,3.334E-01, XSE  900
     *     3.356E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /1.091E00,1.078E00,1.065E00,1.052E00,9.134E-01,9.014E-01,XSE  920
     *     8.903E-01,8.791E-01,8.141E-01,8.034E-01,7.937E-01,7.838E-01, XSE  930
     *     7.335E-01,7.244E-01,7.163E-01,7.081E-01,6.517E-01,6.449E-01, XSE  940
     *     6.390E-01,6.330E-01,5.603E-01,5.567E-01,5.540E-01,5.513E-01, XSE  950
     *     0.0,7*0.0/                                                   XSE  960
      DATA X16 /4.396E00,4.354E00,4.314E00,4.273E00,3.796E00,3.760E00,  XSE  970
     *     3.726E00,3.692E00,3.480E00,3.448E00,3.419E00,3.390E00,       XSE  980
     *     3.238E00,3.211E00,3.187E00,3.164E00,3.000E00,2.980E00,       XSE  990
     *     2.963E00,2.946E00,2.742E00,2.732E00,2.724E00,2.716E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /1.792E00,1.765E00,1.739E00,1.713E00,1.406E00,1.382E00,  XSE 1020
     *     1.361E00,1.339E00,1.202E00,1.182E00,1.163E00,1.145E00,       XSE 1030
     *     1.046E00,1.029E00,1.014E00,9.988E-01,8.931E-01,8.806E-01,    XSE 1040
     *     8.698E-01,8.589E-01,7.265E-01,7.201E-01,7.151E-01,7.101E-01, XSE 1050
     *     0.0,7*0.0/                                                   XSE 1060
      DATA X18 /9.239E-02,9.625E-02,1.000E-01,1.039E-01,1.090E-01,      XSE 1070
     *     1.135E-01,1.176E-01,1.219E-01,1.296E-01,1.345E-01,1.388E-01, XSE 1080
     *     1.433E-01,1.559E-01,1.607E-01,1.650E-01,1.695E-01,1.911E-01, XSE 1090
     *     1.954E-01,1.991E-01,2.029E-01,2.425E-01,2.453E-01,2.472E-01, XSE 1100
     *     2.492E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /3.126E00,3.088E00,3.051E00,3.014E00,2.620E00,2.586E00,  XSE 1120
     *     2.555E00,2.524E00,2.341E00,2.311E00,2.284E00,2.257E00,       XSE 1130
     *     2.118E00,2.093E00,2.070E00,2.047E00,1.892E00,1.873E00,       XSE 1140
     *     1.856E00,1.839E00,1.637E00,1.627E00,1.619E00,1.612E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /7.586E-02,7.890E-02,8.190E-02,8.503E-02,8.566E-02,      XSE 1170
     *     8.934E-02,9.277E-02,9.632E-02,1.011E-01,1.052E-01,1.090E-01, XSE 1180
     *     1.128E-01,1.223E-01,1.265E-01,1.303E-01,1.341E-01,1.519E-01, XSE 1190
     *     1.558E-01,1.591E-01,1.625E-01,1.971E-01,1.996E-01,2.015E-01, XSE 1200
     *     2.034E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC29(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPTTTC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.169E00,6.224E00,9.183E00,9.183E00,1.110E01,      XSE  120
     *     1.296E01,1.477E01,1.477E01,1.721E01,1.956E01,2.184E01,0.0,19*XSE  130
     *     0.0/                                                         XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /2.149E-01,2.150E-01,2.150E-01,2.150E-01,2.149E-01,      XSE  200
     *     2.150E-01,2.150E-01,2.150E-01,2.149E-01,2.150E-01,2.150E-01, XSE  210
     *     2.150E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /1.003E-02,9.981E-03,9.947E-03,9.912E-03,1.003E-02,      XSE  230
     *     9.981E-03,9.947E-03,9.912E-03,1.003E-02,9.981E-03,9.947E-03, XSE  240
     *     9.912E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /1.768E-01,1.770E-01,1.770E-01,1.771E-01,1.768E-01,      XSE  260
     *     1.770E-01,1.770E-01,1.771E-01,1.768E-01,1.770E-01,1.770E-01, XSE  270
     *     1.771E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.159E00,2.159E00,2.159E00,2.159E00,2.159E00,2.159E00,  XSE  290
     *     2.159E00,2.159E00,2.159E00,2.159E00,2.159E00,2.159E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /3.299E-01,3.300E-01,3.300E-01,3.300E-01,3.299E-01,      XSE  320
     *     3.300E-01,3.300E-01,3.300E-01,3.299E-01,3.300E-01,3.300E-01, XSE  330
     *     3.300E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /3.774E-01,3.762E-01,3.756E-01,3.749E-01,3.774E-01,      XSE  350
     *     3.762E-01,3.756E-01,3.749E-01,3.774E-01,3.762E-01,3.756E-01, XSE  360
     *     3.749E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /2.955E-01,2.957E-01,2.957E-01,2.958E-01,2.955E-01,      XSE  380
     *     2.957E-01,2.957E-01,2.958E-01,2.955E-01,2.957E-01,2.957E-01, XSE  390
     *     2.958E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.400E00,1.400E00,1.400E00,1.400E00,1.400E00,1.400E00,  XSE  410
     *     1.400E00,1.400E00,1.400E00,1.400E00,1.400E00,1.400E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /1.626E-01,1.627E-01,1.626E-01,1.626E-01,1.626E-01,      XSE  440
     *     1.627E-01,1.626E-01,1.626E-01,1.626E-01,1.627E-01,1.626E-01, XSE  450
     *     1.626E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /4.409E-02,4.387E-02,4.373E-02,4.358E-02,4.409E-02,      XSE  470
     *     4.387E-02,4.373E-02,4.358E-02,4.409E-02,4.387E-02,4.373E-02, XSE  480
     *     4.358E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /2.041E-01,2.043E-01,2.043E-01,2.043E-01,2.041E-01,      XSE  500
     *     2.043E-01,2.043E-01,2.043E-01,2.041E-01,2.043E-01,2.043E-01, XSE  510
     *     2.043E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.611E00,1.611E00,1.610E00,1.610E00,1.611E00,1.611E00,  XSE  530
     *     1.610E00,1.610E00,1.611E00,1.611E00,1.610E00,1.610E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /2.320E-01,2.322E-01,2.321E-01,2.321E-01,2.320E-01,      XSE  560
     *     2.322E-01,2.321E-01,2.321E-01,2.320E-01,2.322E-01,2.321E-01, XSE  570
     *     2.321E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /4.083E-01,4.070E-01,4.063E-01,4.055E-01,4.083E-01,      XSE  590
     *     4.070E-01,4.063E-01,4.055E-01,4.083E-01,4.070E-01,4.063E-01, XSE  600
     *     4.055E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /2.370E-01,2.371E-01,2.371E-01,2.371E-01,2.370E-01,      XSE  620
     *     2.371E-01,2.371E-01,2.371E-01,2.370E-01,2.371E-01,2.371E-01, XSE  630
     *     2.371E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /1.827E00,1.828E00,1.827E00,1.827E00,1.827E00,1.828E00,  XSE  650
     *     1.827E00,1.827E00,1.827E00,1.828E00,1.827E00,1.827E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /1.762E-01,1.764E-01,1.763E-01,1.762E-01,1.762E-01,      XSE  680
     *     1.764E-01,1.763E-01,1.762E-01,1.762E-01,1.764E-01,1.763E-01, XSE  690
     *     1.762E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /3.243E-01,3.231E-01,3.224E-01,3.217E-01,3.243E-01,      XSE  710
     *     3.231E-01,3.224E-01,3.217E-01,3.243E-01,3.231E-01,3.224E-01, XSE  720
     *     3.217E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /7.617E-01,7.624E-01,7.626E-01,7.627E-01,7.617E-01,      XSE  740
     *     7.624E-01,7.626E-01,7.627E-01,7.617E-01,7.624E-01,7.626E-01, XSE  750
     *     7.627E-01,0.0,19*0.0/                                        XSE  760
      DATA X20 /2.493E-01,2.483E-01,2.476E-01,2.469E-01,2.493E-01,      XSE  770
     *     2.483E-01,2.476E-01,2.469E-01,2.493E-01,2.483E-01,2.476E-01, XSE  780
     *     2.469E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC30(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPTTTA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.650E-02,1.047E-01,2.046E-01,2.046E-01,2.616E-01, XSE  120
     *     3.309E-01,4.125E-01,4.125E-01,5.634E-01,7.321E-01,9.181E-01, XSE  130
     *     0.0,19*0.0/                                                  XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /2.420E-01,2.411E-01,2.400E-01,2.389E-01,2.420E-01,      XSE  200
     *     2.411E-01,2.400E-01,2.389E-01,2.420E-01,2.411E-01,2.400E-01, XSE  210
     *     2.389E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /4.166E-03,4.340E-03,4.536E-03,4.728E-03,4.166E-03,      XSE  230
     *     4.340E-03,4.536E-03,4.728E-03,4.166E-03,4.340E-03,4.536E-03, XSE  240
     *     4.728E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /2.022E-01,2.014E-01,2.004E-01,1.995E-01,2.022E-01,      XSE  260
     *     2.014E-01,2.004E-01,1.995E-01,2.022E-01,2.014E-01,2.004E-01, XSE  270
     *     1.995E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.269E00,2.265E00,2.260E00,2.256E00,2.269E00,2.265E00,  XSE  290
     *     2.260E00,2.256E00,2.269E00,2.265E00,2.260E00,2.256E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /3.695E-01,3.681E-01,3.664E-01,3.649E-01,3.695E-01,      XSE  320
     *     3.681E-01,3.664E-01,3.649E-01,3.695E-01,3.681E-01,3.664E-01, XSE  330
     *     3.649E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /2.548E-01,2.585E-01,2.628E-01,2.670E-01,2.548E-01,      XSE  350
     *     2.585E-01,2.628E-01,2.670E-01,2.548E-01,2.585E-01,2.628E-01, XSE  360
     *     2.670E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /3.496E-01,3.477E-01,3.456E-01,3.435E-01,3.496E-01,      XSE  380
     *     3.477E-01,3.456E-01,3.435E-01,3.496E-01,3.477E-01,3.456E-01, XSE  390
     *     3.435E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.474E00,1.471E00,1.468E00,1.465E00,1.474E00,1.471E00,  XSE  410
     *     1.468E00,1.465E00,1.474E00,1.471E00,1.468E00,1.465E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /1.910E-01,1.900E-01,1.888E-01,1.876E-01,1.910E-01,      XSE  440
     *     1.900E-01,1.888E-01,1.876E-01,1.910E-01,1.900E-01,1.888E-01, XSE  450
     *     1.876E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /1.858E-02,1.934E-02,2.019E-02,2.103E-02,1.858E-02,      XSE  470
     *     1.934E-02,2.019E-02,2.103E-02,1.858E-02,1.934E-02,2.019E-02, XSE  480
     *     2.103E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /2.562E-01,2.544E-01,2.522E-01,2.502E-01,2.562E-01,      XSE  500
     *     2.544E-01,2.522E-01,2.502E-01,2.562E-01,2.544E-01,2.522E-01, XSE  510
     *     2.502E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.603E00,1.603E00,1.603E00,1.603E00,1.603E00,1.603E00,  XSE  530
     *     1.603E00,1.603E00,1.603E00,1.603E00,1.603E00,1.603E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /2.829E-01,2.810E-01,2.790E-01,2.770E-01,2.829E-01,      XSE  560
     *     2.810E-01,2.790E-01,2.770E-01,2.829E-01,2.810E-01,2.790E-01, XSE  570
     *     2.770E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /2.717E-01,2.759E-01,2.806E-01,2.853E-01,2.717E-01,      XSE  590
     *     2.759E-01,2.806E-01,2.853E-01,2.717E-01,2.759E-01,2.806E-01, XSE  600
     *     2.853E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /2.770E-01,2.756E-01,2.740E-01,2.724E-01,2.770E-01,      XSE  620
     *     2.756E-01,2.740E-01,2.724E-01,2.770E-01,2.756E-01,2.740E-01, XSE  630
     *     2.724E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /1.918E00,1.915E00,1.911E00,1.908E00,1.918E00,1.915E00,  XSE  650
     *     1.911E00,1.908E00,1.918E00,1.915E00,1.911E00,1.908E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /2.298E-01,2.279E-01,2.256E-01,2.235E-01,2.298E-01,      XSE  680
     *     2.279E-01,2.256E-01,2.235E-01,2.298E-01,2.279E-01,2.256E-01, XSE  690
     *     2.235E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /1.983E-01,2.021E-01,2.064E-01,2.107E-01,1.983E-01,      XSE  710
     *     2.021E-01,2.064E-01,2.107E-01,1.983E-01,2.021E-01,2.064E-01, XSE  720
     *     2.107E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /9.012E-01,8.965E-01,8.910E-01,8.857E-01,9.012E-01,      XSE  740
     *     8.965E-01,8.910E-01,8.857E-01,9.012E-01,8.965E-01,8.910E-01, XSE  750
     *     8.857E-01,0.0,19*0.0/                                        XSE  760
      DATA X20 /1.337E-01,1.372E-01,1.411E-01,1.450E-01,1.337E-01,      XSE  770
     *     1.372E-01,1.411E-01,1.450E-01,1.337E-01,1.372E-01,1.411E-01, XSE  780
     *     1.450E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC31(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPTTTR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,1.584E-03,4.489E-03,8.795E-03,8.795E-03,1.748E-02, XSE  120
     *     2.868E-02,4.245E-02,4.245E-02,6.667E-02,9.546E-02,1.289E-01, XSE  130
     *     1.289E-01,1.850E-01,2.492E-01,3.215E-01,3.215E-01,4.404E-01, XSE  140
     *     5.726E-01,7.182E-01,7.182E-01,9.521E-01,1.206E00,1.480E00,   XSE  150
     *     0.0,7*0.0/                                                   XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /2.642E-01,2.629E-01,2.617E-01,2.605E-01,2.548E-01,      XSE  220
     *     2.535E-01,2.523E-01,2.512E-01,2.483E-01,2.470E-01,2.459E-01, XSE  230
     *     2.448E-01,2.420E-01,2.408E-01,2.398E-01,2.388E-01,2.352E-01, XSE  240
     *     2.341E-01,2.333E-01,2.325E-01,2.273E-01,2.266E-01,2.261E-01, XSE  250
     *     2.255E-01,0.0,7*0.0/                                         XSE  260
      DATA X02 /1.107E-03,1.256E-03,1.407E-03,1.557E-03,1.561E-03,      XSE  270
     *     1.745E-03,1.920E-03,2.093E-03,2.234E-03,2.446E-03,2.640E-03, XSE  280
     *     2.830E-03,3.149E-03,3.375E-03,3.576E-03,3.774E-03,4.412E-03, XSE  290
     *     4.629E-03,4.813E-03,4.993E-03,6.286E-03,6.444E-03,6.560E-03, XSE  300
     *     6.675E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /2.250E-01,2.239E-01,2.228E-01,2.217E-01,2.168E-01,      XSE  320
     *     2.156E-01,2.146E-01,2.135E-01,2.110E-01,2.099E-01,2.088E-01, XSE  330
     *     2.078E-01,2.052E-01,2.041E-01,2.031E-01,2.021E-01,1.986E-01, XSE  340
     *     1.976E-01,1.968E-01,1.960E-01,1.908E-01,1.901E-01,1.896E-01, XSE  350
     *     1.891E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /2.391E00,2.385E00,2.379E00,2.374E00,2.335E00,2.329E00,  XSE  370
     *     2.324E00,2.319E00,2.303E00,2.298E00,2.293E00,2.288E00,       XSE  380
     *     2.276E00,2.271E00,2.266E00,2.262E00,2.246E00,2.242E00,       XSE  390
     *     2.238E00,2.235E00,2.214E00,2.211E00,2.209E00,2.207E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /4.158E-01,4.136E-01,4.116E-01,4.095E-01,3.929E-01,      XSE  420
     *     3.909E-01,3.890E-01,3.872E-01,3.811E-01,3.792E-01,3.774E-01, XSE  430
     *     3.757E-01,3.713E-01,3.695E-01,3.679E-01,3.663E-01,3.611E-01, XSE  440
     *     3.595E-01,3.582E-01,3.570E-01,3.498E-01,3.487E-01,3.479E-01, XSE  450
     *     3.471E-01,0.0,7*0.0/                                         XSE  460
      DATA X06 /1.753E-01,1.796E-01,1.838E-01,1.879E-01,1.898E-01,      XSE  470
     *     1.947E-01,1.992E-01,2.036E-01,2.091E-01,2.142E-01,2.189E-01, XSE  480
     *     2.235E-01,2.333E-01,2.385E-01,2.430E-01,2.475E-01,2.640E-01, XSE  490
     *     2.687E-01,2.725E-01,2.763E-01,3.040E-01,3.072E-01,3.095E-01, XSE  500
     *     3.117E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /4.037E-01,4.009E-01,3.983E-01,3.958E-01,3.809E-01,      XSE  520
     *     3.783E-01,3.758E-01,3.734E-01,3.669E-01,3.643E-01,3.619E-01, XSE  530
     *     3.597E-01,3.538E-01,3.514E-01,3.492E-01,3.471E-01,3.396E-01, XSE  540
     *     3.376E-01,3.358E-01,3.341E-01,3.236E-01,3.222E-01,3.212E-01, XSE  550
     *     3.201E-01,0.0,7*0.0/                                         XSE  560
      DATA X08 /1.561E00,1.557E00,1.553E00,1.549E00,1.518E00,1.514E00,  XSE  570
     *     1.510E00,1.507E00,1.495E00,1.491E00,1.488E00,1.485E00,       XSE  580
     *     1.476E00,1.472E00,1.469E00,1.466E00,1.456E00,1.453E00,       XSE  590
     *     1.451E00,1.448E00,1.435E00,1.433E00,1.432E00,1.430E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /2.254E-01,2.237E-01,2.221E-01,2.206E-01,2.084E-01,      XSE  620
     *     2.069E-01,2.055E-01,2.042E-01,1.996E-01,1.982E-01,1.969E-01, XSE  630
     *     1.957E-01,1.923E-01,1.910E-01,1.899E-01,1.887E-01,1.848E-01, XSE  640
     *     1.837E-01,1.828E-01,1.819E-01,1.767E-01,1.760E-01,1.754E-01, XSE  650
     *     1.748E-01,0.0,7*0.0/                                         XSE  660
      DATA X10 /5.094E-03,5.752E-03,6.416E-03,7.078E-03,7.104E-03,      XSE  670
     *     7.914E-03,8.681E-03,9.439E-03,1.007E-02,1.100E-02,1.185E-02, XSE  680
     *     1.268E-02,1.410E-02,1.509E-02,1.596E-02,1.682E-02,1.964E-02, XSE  690
     *     2.058E-02,2.138E-02,2.216E-02,2.781E-02,2.850E-02,2.900E-02, XSE  700
     *     2.950E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /3.160E-01,3.132E-01,3.105E-01,3.079E-01,2.868E-01,      XSE  720
     *     2.842E-01,2.818E-01,2.794E-01,2.715E-01,2.690E-01,2.667E-01, XSE  730
     *     2.645E-01,2.585E-01,2.562E-01,2.541E-01,2.521E-01,2.449E-01, XSE  740
     *     2.430E-01,2.413E-01,2.396E-01,2.300E-01,2.287E-01,2.276E-01, XSE  750
     *     2.266E-01,0.0,7*0.0/                                         XSE  760
      DATA X12 /1.619E00,1.618E00,1.618E00,1.617E00,1.601E00,1.601E00,  XSE  770
     *     1.601E00,1.601E00,1.597E00,1.597E00,1.597E00,1.597E00,       XSE  780
     *     1.597E00,1.597E00,1.597E00,1.597E00,1.598E00,1.599E00,       XSE  790
     *     1.599E00,1.599E00,1.602E00,1.603E00,1.603E00,1.603E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /3.459E-01,3.429E-01,3.401E-01,3.373E-01,3.147E-01,      XSE  820
     *     3.120E-01,3.095E-01,3.071E-01,2.989E-01,2.964E-01,2.941E-01, XSE  830
     *     2.919E-01,2.859E-01,2.836E-01,2.816E-01,2.796E-01,2.727E-01, XSE  840
     *     2.707E-01,2.691E-01,2.675E-01,2.584E-01,2.570E-01,2.560E-01, XSE  850
     *     2.550E-01,0.0,7*0.0/                                         XSE  860
      DATA X14 /1.850E-01,1.896E-01,1.941E-01,1.985E-01,1.993E-01,      XSE  870
     *     2.045E-01,2.095E-01,2.143E-01,2.198E-01,2.255E-01,2.306E-01, XSE  880
     *     2.356E-01,2.462E-01,2.520E-01,2.570E-01,2.619E-01,2.802E-01, XSE  890
     *     2.854E-01,2.897E-01,2.939E-01,3.252E-01,3.288E-01,3.314E-01, XSE  900
     *     3.339E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /3.199E-01,3.178E-01,3.158E-01,3.138E-01,3.002E-01,      XSE  920
     *     2.982E-01,2.963E-01,2.945E-01,2.890E-01,2.871E-01,2.853E-01, XSE  930
     *     2.836E-01,2.791E-01,2.773E-01,2.757E-01,2.742E-01,2.686E-01, XSE  940
     *     2.671E-01,2.658E-01,2.646E-01,2.571E-01,2.560E-01,2.552E-01, XSE  950
     *     2.545E-01,0.0,7*0.0/                                         XSE  960
      DATA X16 /2.027E00,2.022E00,2.017E00,2.012E00,1.971E00,1.967E00,  XSE  970
     *     1.962E00,1.958E00,1.943E00,1.938E00,1.934E00,1.930E00,       XSE  980
     *     1.919E00,1.915E00,1.911E00,1.908E00,1.895E00,1.891E00,       XSE  990
     *     1.889E00,1.886E00,1.870E00,1.867E00,1.866E00,1.864E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /3.010E-01,2.979E-01,2.949E-01,2.920E-01,2.605E-01,      XSE 1020
     *     2.578E-01,2.552E-01,2.528E-01,2.425E-01,2.399E-01,2.376E-01, XSE 1030
     *     2.353E-01,2.289E-01,2.266E-01,2.245E-01,2.225E-01,2.156E-01, XSE 1040
     *     2.137E-01,2.120E-01,2.104E-01,2.018E-01,2.005E-01,1.994E-01, XSE 1050
     *     1.984E-01,0.0,7*0.0/                                         XSE 1060
      DATA X18 /1.191E-01,1.232E-01,1.273E-01,1.313E-01,1.325E-01,      XSE 1070
     *     1.372E-01,1.417E-01,1.460E-01,1.510E-01,1.562E-01,1.608E-01, XSE 1080
     *     1.654E-01,1.749E-01,1.801E-01,1.847E-01,1.892E-01,2.057E-01, XSE 1090
     *     2.104E-01,2.144E-01,2.182E-01,2.469E-01,2.502E-01,2.525E-01, XSE 1100
     *     2.549E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /1.037E00,1.030E00,1.024E00,1.017E00,9.801E-01,9.734E-01,XSE 1120
     *     9.673E-01,9.614E-01,9.451E-01,9.386E-01,9.327E-01,9.270E-01, XSE 1130
     *     9.123E-01,9.061E-01,9.007E-01,8.954E-01,8.763E-01,8.710E-01, XSE 1140
     *     8.666E-01,8.622E-01,8.349E-01,8.313E-01,8.286E-01,8.259E-01, XSE 1150
     *     0.0,7*0.0/                                                   XSE 1160
      DATA X20 /6.430E-02,6.785E-02,7.136E-02,7.483E-02,7.520E-02,      XSE 1170
     *     7.936E-02,8.326E-02,8.711E-02,9.106E-02,9.566E-02,9.981E-02, XSE 1180
     *     1.039E-01,1.119E-01,1.166E-01,1.208E-01,1.249E-01,1.394E-01, XSE 1190
     *     1.437E-01,1.473E-01,1.509E-01,1.772E-01,1.802E-01,1.824E-01, XSE 1200
     *     1.846E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC32(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPUUTC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.222E00,6.324E00,9.317E00,9.317E00,1.132E01,      XSE  120
     *     1.327E01,1.518E01,1.518E01,1.765E01,2.005E01,2.238E01,0.0,19*XSE  130
     *     0.0/                                                         XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /3.884E-01,3.884E-01,3.883E-01,3.882E-01,3.884E-01,      XSE  200
     *     3.884E-01,3.883E-01,3.882E-01,3.884E-01,3.884E-01,3.883E-01, XSE  210
     *     3.882E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /1.032E-02,1.030E-02,1.028E-02,1.026E-02,1.032E-02,      XSE  230
     *     1.030E-02,1.028E-02,1.026E-02,1.032E-02,1.030E-02,1.028E-02, XSE  240
     *     1.026E-02,0.0,19*0.0/                                        XSE  250
      DATA X03 /2.583E-01,2.584E-01,2.584E-01,2.584E-01,2.583E-01,      XSE  260
     *     2.584E-01,2.584E-01,2.584E-01,2.583E-01,2.584E-01,2.584E-01, XSE  270
     *     2.584E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.741E00,2.741E00,2.741E00,2.741E00,2.741E00,2.741E00,  XSE  290
     *     2.741E00,2.741E00,2.741E00,2.741E00,2.741E00,2.741E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /5.672E-01,5.671E-01,5.670E-01,5.668E-01,5.672E-01,      XSE  320
     *     5.671E-01,5.670E-01,5.668E-01,5.672E-01,5.671E-01,5.670E-01, XSE  330
     *     5.668E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /3.118E-01,3.113E-01,3.108E-01,3.104E-01,3.118E-01,      XSE  350
     *     3.113E-01,3.108E-01,3.104E-01,3.118E-01,3.113E-01,3.108E-01, XSE  360
     *     3.104E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /5.562E-01,5.563E-01,5.563E-01,5.563E-01,5.562E-01,      XSE  380
     *     5.563E-01,5.563E-01,5.563E-01,5.562E-01,5.563E-01,5.563E-01, XSE  390
     *     5.563E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.893E00,1.893E00,1.893E00,1.893E00,1.893E00,1.893E00,  XSE  410
     *     1.893E00,1.893E00,1.893E00,1.893E00,1.893E00,1.893E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /2.867E-01,2.868E-01,2.868E-01,2.869E-01,2.867E-01,      XSE  440
     *     2.868E-01,2.868E-01,2.869E-01,2.867E-01,2.868E-01,2.868E-01, XSE  450
     *     2.869E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /4.529E-02,4.519E-02,4.511E-02,4.503E-02,4.529E-02,      XSE  470
     *     4.519E-02,4.511E-02,4.503E-02,4.529E-02,4.519E-02,4.511E-02, XSE  480
     *     4.503E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /4.913E-01,4.913E-01,4.911E-01,4.910E-01,4.913E-01,      XSE  500
     *     4.913E-01,4.911E-01,4.910E-01,4.913E-01,4.913E-01,4.911E-01, XSE  510
     *     4.910E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.788E00,1.787E00,1.787E00,1.787E00,1.788E00,1.787E00,  XSE  530
     *     1.787E00,1.787E00,1.788E00,1.787E00,1.787E00,1.787E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /4.986E-01,4.986E-01,4.985E-01,4.984E-01,4.986E-01,      XSE  560
     *     4.986E-01,4.985E-01,4.984E-01,4.986E-01,4.986E-01,4.985E-01, XSE  570
     *     4.984E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /3.774E-01,3.768E-01,3.763E-01,3.758E-01,3.774E-01,      XSE  590
     *     3.768E-01,3.763E-01,3.758E-01,3.774E-01,3.768E-01,3.763E-01, XSE  600
     *     3.758E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /4.608E-01,4.608E-01,4.607E-01,4.606E-01,4.608E-01,      XSE  620
     *     4.608E-01,4.607E-01,4.606E-01,4.608E-01,4.608E-01,4.607E-01, XSE  630
     *     4.606E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /2.460E00,2.460E00,2.459E00,2.459E00,2.460E00,2.460E00,  XSE  650
     *     2.459E00,2.459E00,2.460E00,2.460E00,2.459E00,2.459E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /5.468E-01,5.467E-01,5.464E-01,5.461E-01,5.468E-01,      XSE  680
     *     5.467E-01,5.464E-01,5.461E-01,5.468E-01,5.467E-01,5.464E-01, XSE  690
     *     5.461E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /2.914E-01,2.909E-01,2.905E-01,2.901E-01,2.914E-01,      XSE  710
     *     2.909E-01,2.905E-01,2.901E-01,2.914E-01,2.909E-01,2.905E-01, XSE  720
     *     2.901E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.367E00,1.367E00,1.367E00,1.367E00,1.367E00,1.367E00,  XSE  740
     *     1.367E00,1.367E00,1.367E00,1.367E00,1.367E00,1.367E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /2.416E-01,2.412E-01,2.408E-01,2.404E-01,2.416E-01,      XSE  770
     *     2.412E-01,2.408E-01,2.404E-01,2.416E-01,2.412E-01,2.408E-01, XSE  780
     *     2.404E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC33(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPUUTA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,8.804E-02,1.984E-01,3.307E-01,3.307E-01,4.166E-01, XSE  120
     *     5.134E-01,6.211E-01,6.211E-01,7.886E-01,9.710E-01,1.168E00,  XSE  130
     *     0.0,19*0.0/                                                  XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /6.390E-01,6.307E-01,6.227E-01,6.151E-01,6.390E-01,      XSE  200
     *     6.307E-01,6.227E-01,6.151E-01,6.390E-01,6.307E-01,6.227E-01, XSE  210
     *     6.151E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /4.913E-03,5.068E-03,5.224E-03,5.378E-03,4.913E-03,      XSE  230
     *     5.068E-03,5.224E-03,5.378E-03,4.913E-03,5.068E-03,5.224E-03, XSE  240
     *     5.378E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /3.626E-01,3.592E-01,3.560E-01,3.529E-01,3.626E-01,      XSE  260
     *     3.592E-01,3.560E-01,3.529E-01,3.626E-01,3.592E-01,3.560E-01, XSE  270
     *     3.529E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /3.458E00,3.435E00,3.412E00,3.390E00,3.458E00,3.435E00,  XSE  290
     *     3.412E00,3.390E00,3.458E00,3.435E00,3.412E00,3.390E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /9.211E-01,9.090E-01,8.973E-01,8.862E-01,9.211E-01,      XSE  320
     *     9.090E-01,8.973E-01,8.862E-01,9.211E-01,9.090E-01,8.973E-01, XSE  330
     *     8.862E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /1.857E-01,1.892E-01,1.928E-01,1.963E-01,1.857E-01,      XSE  350
     *     1.892E-01,1.928E-01,1.963E-01,1.857E-01,1.892E-01,1.928E-01, XSE  360
     *     1.963E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /8.646E-01,8.546E-01,8.449E-01,8.356E-01,8.646E-01,      XSE  380
     *     8.546E-01,8.449E-01,8.356E-01,8.646E-01,8.546E-01,8.449E-01, XSE  390
     *     8.356E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /2.486E00,2.467E00,2.448E00,2.430E00,2.486E00,2.467E00,  XSE  410
     *     2.448E00,2.430E00,2.486E00,2.467E00,2.448E00,2.430E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /3.809E-01,3.780E-01,3.751E-01,3.723E-01,3.809E-01,      XSE  440
     *     3.780E-01,3.751E-01,3.723E-01,3.809E-01,3.780E-01,3.751E-01, XSE  450
     *     3.723E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /2.177E-02,2.245E-02,2.313E-02,2.380E-02,2.177E-02,      XSE  470
     *     2.245E-02,2.313E-02,2.380E-02,2.177E-02,2.245E-02,2.313E-02, XSE  480
     *     2.380E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /9.069E-01,8.932E-01,8.799E-01,8.672E-01,9.069E-01,      XSE  500
     *     8.932E-01,8.799E-01,8.672E-01,9.069E-01,8.932E-01,8.799E-01, XSE  510
     *     8.672E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /2.120E00,2.108E00,2.096E00,2.085E00,2.120E00,2.108E00,  XSE  530
     *     2.096E00,2.085E00,2.120E00,2.108E00,2.096E00,2.085E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /9.420E-01,9.257E-01,9.102E-01,8.956E-01,9.420E-01,      XSE  560
     *     9.257E-01,9.102E-01,8.956E-01,9.420E-01,9.257E-01,9.102E-01, XSE  570
     *     8.956E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /2.547E-01,2.582E-01,2.617E-01,2.651E-01,2.547E-01,      XSE  590
     *     2.582E-01,2.617E-01,2.651E-01,2.547E-01,2.582E-01,2.617E-01, XSE  600
     *     2.651E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /7.300E-01,7.215E-01,7.132E-01,7.052E-01,7.300E-01,      XSE  620
     *     7.215E-01,7.132E-01,7.052E-01,7.300E-01,7.215E-01,7.132E-01, XSE  630
     *     7.052E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /3.240E00,3.214E00,3.189E00,3.165E00,3.240E00,3.214E00,  XSE  650
     *     3.189E00,3.165E00,3.240E00,3.214E00,3.189E00,3.165E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /1.019E00,1.003E00,9.879E-01,9.732E-01,1.019E00,1.003E00,XSE  680
     *     9.879E-01,9.732E-01,1.019E00,1.003E00,9.879E-01,9.732E-01,   XSE  690
     *     0.0,19*0.0/                                                  XSE  700
      DATA X18 /1.717E-01,1.751E-01,1.786E-01,1.820E-01,1.717E-01,      XSE  710
     *     1.751E-01,1.786E-01,1.820E-01,1.717E-01,1.751E-01,1.786E-01, XSE  720
     *     1.820E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /2.104E00,2.080E00,2.057E00,2.035E00,2.104E00,2.080E00,  XSE  740
     *     2.057E00,2.035E00,2.104E00,2.080E00,2.057E00,2.035E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /1.389E-01,1.419E-01,1.448E-01,1.478E-01,1.389E-01,      XSE  770
     *     1.419E-01,1.448E-01,1.478E-01,1.389E-01,1.419E-01,1.448E-01, XSE  780
     *     1.478E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC34(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMOPUUTR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,8.610E-03,1.897E-02,3.117E-02,3.117E-02,4.911E-02, XSE  120
     *     6.948E-02,9.243E-02,9.243E-02,1.274E-01,1.664E-01,2.094E-01, XSE  130
     *     2.094E-01,2.765E-01,3.499E-01,4.298E-01,4.298E-01,5.581E-01, XSE  140
     *     6.965E-01,8.449E-01,8.449E-01,1.091E00,1.354E00,1.632E00,0.0,XSE  150
     *     7*0.0/                                                       XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /2.879E-01,2.870E-01,2.860E-01,2.851E-01,2.616E-01,      XSE  220
     *     2.606E-01,2.597E-01,2.589E-01,2.537E-01,2.528E-01,2.519E-01, XSE  230
     *     2.510E-01,2.509E-01,2.500E-01,2.491E-01,2.482E-01,2.518E-01, XSE  240
     *     2.510E-01,2.501E-01,2.493E-01,2.619E-01,2.613E-01,2.607E-01, XSE  250
     *     2.601E-01,0.0,7*0.0/                                         XSE  260
      DATA X02 /1.266E-03,1.388E-03,1.510E-03,1.632E-03,1.748E-03,      XSE  270
     *     1.885E-03,2.021E-03,2.156E-03,2.421E-03,2.570E-03,2.716E-03, XSE  280
     *     2.861E-03,3.365E-03,3.520E-03,3.671E-03,3.820E-03,4.753E-03, XSE  290
     *     4.898E-03,5.038E-03,5.175E-03,6.918E-03,7.010E-03,7.097E-03, XSE  300
     *     7.184E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /2.423E-01,2.415E-01,2.406E-01,2.398E-01,2.206E-01,      XSE  320
     *     2.198E-01,2.189E-01,2.182E-01,2.139E-01,2.131E-01,2.123E-01, XSE  330
     *     2.115E-01,2.103E-01,2.095E-01,2.087E-01,2.079E-01,2.086E-01, XSE  340
     *     2.079E-01,2.072E-01,2.064E-01,2.130E-01,2.124E-01,2.119E-01, XSE  350
     *     2.114E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /2.523E00,2.518E00,2.514E00,2.510E00,2.359E00,2.355E00,  XSE  370
     *     2.351E00,2.348E00,2.320E00,2.316E00,2.312E00,2.309E00,       XSE  380
     *     2.308E00,2.304E00,2.301E00,2.297E00,2.315E00,2.312E00,       XSE  390
     *     2.308E00,2.305E00,2.379E00,2.376E00,2.373E00,2.370E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /4.703E-01,4.687E-01,4.672E-01,4.657E-01,4.011E-01,      XSE  420
     *     3.997E-01,3.984E-01,3.972E-01,3.863E-01,3.850E-01,3.837E-01, XSE  430
     *     3.824E-01,3.830E-01,3.816E-01,3.803E-01,3.790E-01,3.875E-01, XSE  440
     *     3.862E-01,3.849E-01,3.837E-01,4.149E-01,4.137E-01,4.126E-01, XSE  450
     *     4.115E-01,0.0,7*0.0/                                         XSE  460
      DATA X06 /1.681E-01,1.713E-01,1.746E-01,1.777E-01,1.870E-01,      XSE  470
     *     1.904E-01,1.938E-01,1.972E-01,2.064E-01,2.099E-01,2.134E-01, XSE  480
     *     2.168E-01,2.299E-01,2.333E-01,2.367E-01,2.400E-01,2.599E-01, XSE  490
     *     2.629E-01,2.658E-01,2.687E-01,2.982E-01,3.000E-01,3.017E-01, XSE  500
     *     3.034E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /4.568E-01,4.548E-01,4.529E-01,4.509E-01,3.922E-01,      XSE  520
     *     3.904E-01,3.886E-01,3.868E-01,3.752E-01,3.733E-01,3.715E-01, XSE  530
     *     3.698E-01,3.690E-01,3.672E-01,3.654E-01,3.637E-01,3.703E-01, XSE  540
     *     3.686E-01,3.670E-01,3.654E-01,3.930E-01,3.917E-01,3.904E-01, XSE  550
     *     3.892E-01,0.0,7*0.0/                                         XSE  560
      DATA X08 /1.665E00,1.662E00,1.659E00,1.657E00,1.536E00,1.533E00,  XSE  570
     *     1.531E00,1.528E00,1.507E00,1.504E00,1.501E00,1.499E00,       XSE  580
     *     1.500E00,1.497E00,1.495E00,1.492E00,1.509E00,1.506E00,       XSE  590
     *     1.504E00,1.502E00,1.565E00,1.563E00,1.560E00,1.558E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /2.652E-01,2.641E-01,2.630E-01,2.619E-01,2.147E-01,      XSE  620
     *     2.137E-01,2.127E-01,2.118E-01,2.043E-01,2.033E-01,2.023E-01, XSE  630
     *     2.014E-01,2.018E-01,2.008E-01,1.999E-01,1.989E-01,2.051E-01, XSE  640
     *     2.041E-01,2.032E-01,2.023E-01,2.242E-01,2.234E-01,2.227E-01, XSE  650
     *     2.219E-01,0.0,7*0.0/                                         XSE  660
      DATA X10 /5.743E-03,6.278E-03,6.814E-03,7.348E-03,7.865E-03,      XSE  670
     *     8.465E-03,9.060E-03,9.650E-03,1.082E-02,1.147E-02,1.211E-02, XSE  680
     *     1.274E-02,1.497E-02,1.564E-02,1.630E-02,1.695E-02,2.104E-02, XSE  690
     *     2.166E-02,2.227E-02,2.287E-02,3.044E-02,3.084E-02,3.122E-02, XSE  700
     *     3.159E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /3.828E-01,3.808E-01,3.788E-01,3.769E-01,2.957E-01,      XSE  720
     *     2.940E-01,2.923E-01,2.907E-01,2.763E-01,2.746E-01,2.729E-01, XSE  730
     *     2.713E-01,2.707E-01,2.690E-01,2.673E-01,2.657E-01,2.745E-01, XSE  740
     *     2.728E-01,2.712E-01,2.697E-01,3.073E-01,3.059E-01,3.045E-01, XSE  750
     *     3.032E-01,0.0,7*0.0/                                         XSE  760
      DATA X12 /1.663E00,1.663E00,1.663E00,1.663E00,1.603E00,1.603E00,  XSE  770
     *     1.603E00,1.603E00,1.596E00,1.596E00,1.597E00,1.597E00,       XSE  780
     *     1.600E00,1.600E00,1.600E00,1.600E00,1.612E00,1.612E00,       XSE  790
     *     1.612E00,1.612E00,1.647E00,1.647E00,1.646E00,1.646E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /4.171E-01,4.151E-01,4.131E-01,4.112E-01,3.243E-01,      XSE  820
     *     3.226E-01,3.209E-01,3.193E-01,3.062E-01,3.045E-01,3.028E-01, XSE  830
     *     3.012E-01,3.017E-01,3.000E-01,2.983E-01,2.967E-01,3.068E-01, XSE  840
     *     3.052E-01,3.036E-01,3.020E-01,3.393E-01,3.380E-01,3.366E-01, XSE  850
     *     3.353E-01,0.0,7*0.0/                                         XSE  860
      DATA X14 /1.835E-01,1.870E-01,1.904E-01,1.939E-01,1.985E-01,      XSE  870
     *     2.022E-01,2.060E-01,2.096E-01,2.188E-01,2.227E-01,2.265E-01, XSE  880
     *     2.303E-01,2.454E-01,2.492E-01,2.529E-01,2.566E-01,2.809E-01, XSE  890
     *     2.842E-01,2.874E-01,2.906E-01,3.293E-01,3.312E-01,3.331E-01, XSE  900
     *     3.349E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /3.666E-01,3.650E-01,3.635E-01,3.620E-01,3.087E-01,      XSE  920
     *     3.073E-01,3.060E-01,3.047E-01,2.946E-01,2.932E-01,2.919E-01, XSE  930
     *     2.906E-01,2.904E-01,2.890E-01,2.877E-01,2.865E-01,2.930E-01, XSE  940
     *     2.918E-01,2.905E-01,2.893E-01,3.157E-01,3.147E-01,3.137E-01, XSE  950
     *     3.126E-01,0.0,7*0.0/                                         XSE  960
      DATA X16 /2.160E00,2.157E00,2.153E00,2.149E00,1.991E00,1.988E00,  XSE  970
     *     1.985E00,1.982E00,1.954E00,1.951E00,1.948E00,1.945E00,       XSE  980
     *     1.946E00,1.943E00,1.940E00,1.937E00,1.957E00,1.954E00,       XSE  990
     *     1.952E00,1.949E00,2.028E00,2.026E00,2.023E00,2.021E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /3.947E-01,3.925E-01,3.904E-01,3.884E-01,2.683E-01,      XSE 1020
     *     2.667E-01,2.651E-01,2.635E-01,2.443E-01,2.426E-01,2.411E-01, XSE 1030
     *     2.396E-01,2.401E-01,2.385E-01,2.369E-01,2.354E-01,2.495E-01, XSE 1040
     *     2.478E-01,2.462E-01,2.446E-01,3.020E-01,3.003E-01,2.987E-01, XSE 1050
     *     2.971E-01,0.0,7*0.0/                                         XSE 1060
      DATA X18 /1.170E-01,1.202E-01,1.233E-01,1.264E-01,1.325E-01,      XSE 1070
     *     1.359E-01,1.393E-01,1.426E-01,1.511E-01,1.546E-01,1.580E-01, XSE 1080
     *     1.614E-01,1.751E-01,1.786E-01,1.819E-01,1.853E-01,2.073E-01, XSE 1090
     *     2.104E-01,2.133E-01,2.162E-01,2.512E-01,2.530E-01,2.547E-01, XSE 1100
     *     2.564E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /1.164E00,1.159E00,1.154E00,1.149E00,1.003E00,9.986E-01, XSE 1120
     *     9.942E-01,9.898E-01,9.609E-01,9.563E-01,9.519E-01,9.476E-01, XSE 1130
     *     9.436E-01,9.390E-01,9.347E-01,9.304E-01,9.426E-01,9.384E-01, XSE 1140
     *     9.343E-01,9.303E-01,9.918E-01,9.886E-01,9.854E-01,9.822E-01, XSE 1150
     *     0.0,7*0.0/                                                   XSE 1160
      DATA X20 /6.525E-02,6.800E-02,7.075E-02,7.347E-02,7.638E-02,      XSE 1170
     *     7.939E-02,8.236E-02,8.530E-02,9.212E-02,9.528E-02,9.838E-02, XSE 1180
     *     1.014E-01,1.134E-01,1.165E-01,1.196E-01,1.226E-01,1.428E-01, XSE 1190
     *     1.457E-01,1.484E-01,1.510E-01,1.851E-01,1.868E-01,1.884E-01, XSE 1200
     *     1.900E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC35(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO0TTTC;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,3.614E00,7.159E00,1.064E01,1.064E01,1.257E01,      XSE  120
     *     1.453E01,1.651E01,1.651E01,1.888E01,2.125E01,2.362E01,0.0,19*XSE  130
     *     0.0/                                                         XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /3.237E-01,3.251E-01,3.265E-01,3.279E-01,3.237E-01,      XSE  200
     *     3.251E-01,3.265E-01,3.279E-01,3.237E-01,3.251E-01,3.265E-01, XSE  210
     *     3.279E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /1.191E-02,1.182E-02,1.172E-02,1.162E-02,1.191E-02,      XSE  230
     *     1.182E-02,1.172E-02,1.162E-02,1.191E-02,1.182E-02,1.172E-02, XSE  240
     *     1.162E-02,0.0,19*0.0/                                        XSE  250
      DATA X03 /2.352E-01,2.362E-01,2.372E-01,2.381E-01,2.352E-01,      XSE  260
     *     2.362E-01,2.372E-01,2.381E-01,2.352E-01,2.362E-01,2.372E-01, XSE  270
     *     2.381E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /2.595E00,2.601E00,2.607E00,2.613E00,2.595E00,2.601E00,  XSE  290
     *     2.607E00,2.613E00,2.595E00,2.601E00,2.607E00,2.613E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /5.029E-01,5.053E-01,5.076E-01,5.098E-01,5.029E-01,      XSE  320
     *     5.053E-01,5.076E-01,5.098E-01,5.029E-01,5.053E-01,5.076E-01, XSE  330
     *     5.098E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /3.698E-01,3.674E-01,3.650E-01,3.627E-01,3.698E-01,      XSE  350
     *     3.674E-01,3.650E-01,3.627E-01,3.698E-01,3.674E-01,3.650E-01, XSE  360
     *     3.627E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /4.932E-01,4.959E-01,4.985E-01,5.010E-01,4.932E-01,      XSE  380
     *     4.959E-01,4.985E-01,5.010E-01,4.932E-01,4.959E-01,4.985E-01, XSE  390
     *     5.010E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /1.780E00,1.784E00,1.789E00,1.793E00,1.780E00,1.784E00,  XSE  410
     *     1.789E00,1.793E00,1.780E00,1.784E00,1.789E00,1.793E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /2.915E-01,2.932E-01,2.948E-01,2.963E-01,2.915E-01,      XSE  440
     *     2.932E-01,2.948E-01,2.963E-01,2.915E-01,2.932E-01,2.948E-01, XSE  450
     *     2.963E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /5.273E-02,5.233E-02,5.190E-02,5.149E-02,5.273E-02,      XSE  470
     *     5.233E-02,5.190E-02,5.149E-02,5.273E-02,5.233E-02,5.190E-02, XSE  480
     *     5.149E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /4.237E-01,4.268E-01,4.298E-01,4.327E-01,4.237E-01,      XSE  500
     *     4.268E-01,4.298E-01,4.327E-01,4.237E-01,4.268E-01,4.298E-01, XSE  510
     *     4.327E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /1.758E00,1.759E00,1.761E00,1.762E00,1.758E00,1.759E00,  XSE  530
     *     1.761E00,1.762E00,1.758E00,1.759E00,1.761E00,1.762E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /4.435E-01,4.464E-01,4.492E-01,4.518E-01,4.435E-01,      XSE  560
     *     4.464E-01,4.492E-01,4.518E-01,4.435E-01,4.464E-01,4.492E-01, XSE  570
     *     4.518E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /4.389E-01,4.365E-01,4.339E-01,4.315E-01,4.389E-01,      XSE  590
     *     4.365E-01,4.339E-01,4.315E-01,4.389E-01,4.365E-01,4.339E-01, XSE  600
     *     4.315E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /4.079E-01,4.102E-01,4.123E-01,4.144E-01,4.079E-01,      XSE  620
     *     4.102E-01,4.123E-01,4.144E-01,4.079E-01,4.102E-01,4.123E-01, XSE  630
     *     4.144E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /2.318E00,2.324E00,2.329E00,2.335E00,2.318E00,2.324E00,  XSE  650
     *     2.329E00,2.335E00,2.318E00,2.324E00,2.329E00,2.335E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /4.685E-01,4.725E-01,4.762E-01,4.798E-01,4.685E-01,      XSE  680
     *     4.725E-01,4.762E-01,4.798E-01,4.685E-01,4.725E-01,4.762E-01, XSE  690
     *     4.798E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /3.453E-01,3.431E-01,3.408E-01,3.386E-01,3.453E-01,      XSE  710
     *     3.431E-01,3.408E-01,3.386E-01,3.453E-01,3.431E-01,3.408E-01, XSE  720
     *     3.386E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.213E00,1.219E00,1.225E00,1.231E00,1.213E00,1.219E00,  XSE  740
     *     1.225E00,1.231E00,1.213E00,1.219E00,1.225E00,1.231E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /2.887E-01,2.867E-01,2.847E-01,2.827E-01,2.887E-01,      XSE  770
     *     2.867E-01,2.847E-01,2.827E-01,2.887E-01,2.867E-01,2.847E-01, XSE  780
     *     2.827E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC36(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO0TTTA;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,12/                                           XSE  110
      DATA DBUP /0.0,2.074E-02,6.106E-02,1.234E-01,1.234E-01,1.544E-01, XSE  120
     *     1.932E-01,2.408E-01,2.408E-01,3.311E-01,4.360E-01,5.566E-01, XSE  130
     *     0.0,19*0.0/                                                  XSE  140
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  150
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  160
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  170
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  180
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  190
      DATA X01 /4.293E-01,4.269E-01,4.244E-01,4.218E-01,4.293E-01,      XSE  200
     *     4.269E-01,4.244E-01,4.218E-01,4.293E-01,4.269E-01,4.244E-01, XSE  210
     *     4.218E-01,0.0,19*0.0/                                        XSE  220
      DATA X02 /5.314E-03,5.444E-03,5.585E-03,5.734E-03,5.314E-03,      XSE  230
     *     5.444E-03,5.585E-03,5.734E-03,5.314E-03,5.444E-03,5.585E-03, XSE  240
     *     5.734E-03,0.0,19*0.0/                                        XSE  250
      DATA X03 /3.334E-01,3.312E-01,3.290E-01,3.266E-01,3.334E-01,      XSE  260
     *     3.312E-01,3.290E-01,3.266E-01,3.334E-01,3.312E-01,3.290E-01, XSE  270
     *     3.266E-01,0.0,19*0.0/                                        XSE  280
      DATA X04 /3.260E00,3.244E00,3.229E00,3.213E00,3.260E00,3.244E00,  XSE  290
     *     3.229E00,3.213E00,3.260E00,3.244E00,3.229E00,3.213E00,0.0,19*XSE  300
     *     0.0/                                                         XSE  310
      DATA X05 /8.170E-01,8.095E-01,8.019E-01,7.942E-01,8.170E-01,      XSE  320
     *     8.095E-01,8.019E-01,7.942E-01,8.170E-01,8.095E-01,8.019E-01, XSE  330
     *     7.942E-01,0.0,19*0.0/                                        XSE  340
      DATA X06 /2.140E-01,2.170E-01,2.203E-01,2.237E-01,2.140E-01,      XSE  350
     *     2.170E-01,2.203E-01,2.237E-01,2.140E-01,2.170E-01,2.203E-01, XSE  360
     *     2.237E-01,0.0,19*0.0/                                        XSE  370
      DATA X07 /7.776E-01,7.712E-01,7.645E-01,7.578E-01,7.776E-01,      XSE  380
     *     7.712E-01,7.645E-01,7.578E-01,7.776E-01,7.712E-01,7.645E-01, XSE  390
     *     7.578E-01,0.0,19*0.0/                                        XSE  400
      DATA X08 /2.317E00,2.305E00,2.292E00,2.279E00,2.317E00,2.305E00,  XSE  410
     *     2.292E00,2.279E00,2.317E00,2.305E00,2.292E00,2.279E00,0.0,19*XSE  420
     *     0.0/                                                         XSE  430
      DATA X09 /5.104E-01,5.049E-01,4.993E-01,4.938E-01,5.104E-01,      XSE  440
     *     5.049E-01,4.993E-01,4.938E-01,5.104E-01,5.049E-01,4.993E-01, XSE  450
     *     4.938E-01,0.0,19*0.0/                                        XSE  460
      DATA X10 /2.379E-02,2.437E-02,2.499E-02,2.564E-02,2.379E-02,      XSE  470
     *     2.437E-02,2.499E-02,2.564E-02,2.379E-02,2.437E-02,2.499E-02, XSE  480
     *     2.564E-02,0.0,19*0.0/                                        XSE  490
      DATA X11 /7.966E-01,7.876E-01,7.785E-01,7.693E-01,7.966E-01,      XSE  500
     *     7.876E-01,7.785E-01,7.693E-01,7.966E-01,7.876E-01,7.785E-01, XSE  510
     *     7.693E-01,0.0,19*0.0/                                        XSE  520
      DATA X12 /2.023E00,2.015E00,2.008E00,2.001E00,2.023E00,2.015E00,  XSE  530
     *     2.008E00,2.001E00,2.023E00,2.015E00,2.008E00,2.001E00,0.0,19*XSE  540
     *     0.0/                                                         XSE  550
      DATA X13 /8.174E-01,8.075E-01,7.978E-01,7.880E-01,8.174E-01,      XSE  560
     *     8.075E-01,7.978E-01,7.880E-01,8.174E-01,8.075E-01,7.978E-01, XSE  570
     *     7.880E-01,0.0,19*0.0/                                        XSE  580
      DATA X14 /2.814E-01,2.844E-01,2.877E-01,2.912E-01,2.814E-01,      XSE  590
     *     2.844E-01,2.877E-01,2.912E-01,2.814E-01,2.844E-01,2.877E-01, XSE  600
     *     2.912E-01,0.0,19*0.0/                                        XSE  610
      DATA X15 /6.545E-01,6.489E-01,6.432E-01,6.374E-01,6.545E-01,      XSE  620
     *     6.489E-01,6.432E-01,6.374E-01,6.545E-01,6.489E-01,6.432E-01, XSE  630
     *     6.374E-01,0.0,19*0.0/                                        XSE  640
      DATA X16 /3.015E00,2.998E00,2.981E00,2.965E00,3.015E00,2.998E00,  XSE  650
     *     2.981E00,2.965E00,3.015E00,2.998E00,2.981E00,2.965E00,0.0,19*XSE  660
     *     0.0/                                                         XSE  670
      DATA X17 /9.032E-01,8.926E-01,8.819E-01,8.711E-01,9.032E-01,      XSE  680
     *     8.926E-01,8.819E-01,8.711E-01,9.032E-01,8.926E-01,8.819E-01, XSE  690
     *     8.711E-01,0.0,19*0.0/                                        XSE  700
      DATA X18 /1.962E-01,1.991E-01,2.022E-01,2.055E-01,1.962E-01,      XSE  710
     *     1.991E-01,2.022E-01,2.055E-01,1.962E-01,1.991E-01,2.022E-01, XSE  720
     *     2.055E-01,0.0,19*0.0/                                        XSE  730
      DATA X19 /1.899E00,1.883E00,1.867E00,1.851E00,1.899E00,1.883E00,  XSE  740
     *     1.867E00,1.851E00,1.899E00,1.883E00,1.867E00,1.851E00,0.0,19*XSE  750
     *     0.0/                                                         XSE  760
      DATA X20 /1.568E-01,1.594E-01,1.622E-01,1.652E-01,1.568E-01,      XSE  770
     *     1.594E-01,1.622E-01,1.652E-01,1.568E-01,1.594E-01,1.622E-01, XSE  780
     *     1.652E-01,0.0,19*0.0/                                        XSE  790
      NPT = MPT                                                         XSE  800
      NXSEC = MXSEC                                                     XSE  810
      DO 10 I=1,30                                                      XSE  820
         F(I,01) = X01(I)                                               XSE  830
         F(I,02) = X02(I)                                               XSE  840
         F(I,03) = X03(I)                                               XSE  850
         F(I,04) = X04(I)                                               XSE  860
         F(I,05) = X05(I)                                               XSE  870
         F(I,06) = X06(I)                                               XSE  880
         F(I,07) = X07(I)                                               XSE  890
         F(I,08) = X08(I)                                               XSE  900
         F(I,09) = X09(I)                                               XSE  910
         F(I,10) = X10(I)                                               XSE  920
         F(I,11) = X11(I)                                               XSE  930
         F(I,12) = X12(I)                                               XSE  940
         F(I,13) = X13(I)                                               XSE  950
         F(I,14) = X14(I)                                               XSE  960
         F(I,15) = X15(I)                                               XSE  970
         F(I,16) = X16(I)                                               XSE  980
         F(I,17) = X17(I)                                               XSE  990
         F(I,18) = X18(I)                                               XSE 1000
         F(I,19) = X19(I)                                               XSE 1010
         F(I,20) = X20(I)                                               XSE 1020
         CBUP(I) = DBUP(I)                                              XSE 1030
   10 CONTINUE                                                          XSE 1040
      DO 30 I=1,20                                                      XSE 1050
         DO 20 J=1,7                                                    XSE 1060
            NFUD(J,I) = MFUD(J,I)                                       XSE 1070
   20    CONTINUE                                                       XSE 1080
   30 CONTINUE                                                          XSE 1090
      RETURN                                                            XSE 1100
      END                                                               XSE 1110
      SUBROUTINE XSEC37(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    AMO0TTTR;XSDRNPM/CITATION,XSEC FROM CITATION                       XSE   30
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,24/                                           XSE  110
      DATA DBUP /0.0,5.443E-04,1.549E-03,3.119E-03,3.119E-03,6.254E-03, XSE  120
     *     1.041E-02,1.575E-02,1.575E-02,2.508E-02,3.670E-02,5.092E-02, XSE  130
     *     5.092E-02,7.497E-02,1.040E-01,1.387E-01,1.387E-01,1.969E-01, XSE  140
     *     2.655E-01,3.454E-01,3.454E-01,4.724E-01,6.188E-01,7.859E-01, XSE  150
     *     0.0,7*0.0/                                                   XSE  160
      DATA MFUD /902320,1,5*0,902320,4,5*0,922330,1,5*0,922330,4,5*0,   XSE  170
     *     922340,1,5*0,922340,4,5*0,922350,1,5*0,922350,4,5*0,922380,1,XSE  180
     *     5*0,922380,4,5*0,942390,1,5*0,942390,4,5*0,942400,1,5*0,     XSE  190
     *     942400,4,5*0,942410,1,5*0,942410,4,5*0,942420,1,5*0,942420,4,XSE  200
     *     5*0,952410,1,5*0,952410,2,5*0,7*0/                           XSE  210
      DATA X01 /5.296E-01,5.254E-01,5.214E-01,5.173E-01,4.906E-01,      XSE  220
     *     4.864E-01,4.827E-01,4.788E-01,4.622E-01,4.582E-01,4.547E-01, XSE  230
     *     4.512E-01,4.342E-01,4.306E-01,4.276E-01,4.245E-01,4.027E-01, XSE  240
     *     3.999E-01,3.976E-01,3.952E-01,3.646E-01,3.629E-01,3.618E-01, XSE  250
     *     3.607E-01,0.0,7*0.0/                                         XSE  260
      DATA X02 /1.535E-03,1.676E-03,1.815E-03,1.963E-03,2.079E-03,      XSE  270
     *     2.255E-03,2.415E-03,2.584E-03,2.818E-03,3.018E-03,3.195E-03, XSE  280
     *     3.381E-03,3.811E-03,4.024E-03,4.208E-03,4.400E-03,5.231E-03, XSE  290
     *     5.437E-03,5.604E-03,5.778E-03,7.479E-03,7.634E-03,7.733E-03, XSE  300
     *     7.837E-03,0.0,7*0.0/                                         XSE  310
      DATA X03 /4.982E-01,4.929E-01,4.879E-01,4.828E-01,4.301E-01,      XSE  320
     *     4.253E-01,4.211E-01,4.167E-01,3.920E-01,3.877E-01,3.840E-01, XSE  330
     *     3.801E-01,3.611E-01,3.574E-01,3.543E-01,3.510E-01,3.295E-01, XSE  340
     *     3.266E-01,3.243E-01,3.219E-01,2.933E-01,2.916E-01,2.904E-01, XSE  350
     *     2.893E-01,0.0,7*0.0/                                         XSE  360
      DATA X04 /4.428E00,4.390E00,4.354E00,4.318E00,3.926E00,3.893E00,  XSE  370
     *     3.863E00,3.832E00,3.650E00,3.621E00,3.595E00,3.568E00,       XSE  380
     *     3.430E00,3.405E00,3.383E00,3.361E00,3.209E00,3.190E00,       XSE  390
     *     3.174E00,3.158E00,2.964E00,2.953E00,2.945E00,2.938E00,0.0,7* XSE  400
     *     0.0/                                                         XSE  410
      DATA X05 /1.495E00,1.474E00,1.454E00,1.434E00,1.177E00,1.159E00,  XSE  420
     *     1.143E00,1.127E00,1.018E00,1.003E00,9.899E-01,9.764E-01,     XSE  430
     *     9.042E-01,8.919E-01,8.811E-01,8.701E-01,7.970E-01,7.878E-01, XSE  440
     *     7.801E-01,7.723E-01,6.831E-01,6.779E-01,6.742E-01,6.705E-01, XSE  450
     *     0.0,7*0.0/                                                   XSE  460
      DATA X06 /1.001E-01,1.044E-01,1.085E-01,1.129E-01,1.188E-01,      XSE  470
     *     1.238E-01,1.282E-01,1.329E-01,1.422E-01,1.475E-01,1.521E-01, XSE  480
     *     1.569E-01,1.716E-01,1.769E-01,1.814E-01,1.860E-01,2.104E-01, XSE  490
     *     2.152E-01,2.190E-01,2.229E-01,2.657E-01,2.690E-01,2.710E-01, XSE  500
     *     2.731E-01,0.0,7*0.0/                                         XSE  510
      DATA X07 /1.265E00,1.250E00,1.235E00,1.220E00,1.057E00,1.043E00,  XSE  520
     *     1.031E00,1.018E00,9.425E-01,9.301E-01,9.191E-01,9.079E-01,   XSE  530
     *     8.509E-01,8.402E-01,8.310E-01,8.217E-01,7.584E-01,7.502E-01, XSE  540
     *     7.435E-01,7.366E-01,6.548E-01,6.501E-01,6.469E-01,6.436E-01, XSE  550
     *     0.0,7*0.0/                                                   XSE  560
      DATA X08 /3.257E00,3.227E00,3.199E00,3.169E00,2.850E00,2.824E00,  XSE  570
     *     2.799E00,2.775E00,2.627E00,2.604E00,2.583E00,2.561E00,       XSE  580
     *     2.450E00,2.430E00,2.413E00,2.395E00,2.273E00,2.257E00,       XSE  590
     *     2.245E00,2.232E00,2.078E00,2.069E00,2.063E00,2.057E00,0.0,7* XSE  600
     *     0.0/                                                         XSE  610
      DATA X09 /1.068E00,1.051E00,1.035E00,1.019E00,7.963E-01,7.828E-01,XSE  620
     *     7.702E-01,7.575E-01,6.670E-01,6.557E-01,6.455E-01,6.352E-01, XSE  630
     *     5.794E-01,5.702E-01,5.621E-01,5.539E-01,5.001E-01,4.933E-01, XSE  640
     *     4.876E-01,4.818E-01,4.179E-01,4.142E-01,4.114E-01,4.087E-01, XSE  650
     *     0.0,7*0.0/                                                   XSE  660
      DATA X10 /6.998E-03,7.625E-03,8.240E-03,8.892E-03,9.418E-03,      XSE  670
     *     1.019E-02,1.090E-02,1.165E-02,1.269E-02,1.358E-02,1.436E-02, XSE  680
     *     1.517E-02,1.709E-02,1.803E-02,1.884E-02,1.969E-02,2.337E-02, XSE  690
     *     2.428E-02,2.501E-02,2.577E-02,3.328E-02,3.396E-02,3.439E-02, XSE  700
     *     3.485E-02,0.0,7*0.0/                                         XSE  710
      DATA X11 /1.542E00,1.519E00,1.496E00,1.474E00,1.201E00,1.181E00,  XSE  720
     *     1.163E00,1.145E00,1.026E00,1.008E00,9.930E-01,9.774E-01,     XSE  730
     *     8.947E-01,8.802E-01,8.676E-01,8.548E-01,7.683E-01,7.574E-01, XSE  740
     *     7.483E-01,7.391E-01,6.321E-01,6.260E-01,6.216E-01,6.172E-01, XSE  750
     *     0.0,7*0.0/                                                   XSE  760
      DATA X12 /2.726E00,2.705E00,2.684E00,2.663E00,2.363E00,2.346E00,  XSE  770
     *     2.331E00,2.315E00,2.194E00,2.181E00,2.168E00,2.156E00,       XSE  780
     *     2.083E00,2.072E00,2.062E00,2.053E00,1.986E00,1.978E00,       XSE  790
     *     1.971E00,1.965E00,1.892E00,1.888E00,1.884E00,1.881E00,0.0,7* XSE  800
     *     0.0/                                                         XSE  810
      DATA X13 /1.848E00,1.816E00,1.787E00,1.756E00,1.333E00,1.309E00,  XSE  820
     *     1.286E00,1.263E00,1.094E00,1.074E00,1.055E00,1.037E00,       XSE  830
     *     9.363E-01,9.200E-01,9.056E-01,8.911E-01,7.966E-01,7.846E-01, XSE  840
     *     7.746E-01,7.645E-01,6.538E-01,6.472E-01,6.425E-01,6.376E-01, XSE  850
     *     0.0,7*0.0/                                                   XSE  860
      DATA X14 /1.784E-01,1.823E-01,1.861E-01,1.901E-01,1.903E-01,      XSE  870
     *     1.950E-01,1.992E-01,2.037E-01,2.103E-01,2.154E-01,2.199E-01, XSE  880
     *     2.247E-01,2.374E-01,2.428E-01,2.473E-01,2.519E-01,2.749E-01, XSE  890
     *     2.799E-01,2.837E-01,2.878E-01,3.305E-01,3.340E-01,3.361E-01, XSE  900
     *     3.383E-01,0.0,7*0.0/                                         XSE  910
      DATA X15 /1.066E00,1.053E00,1.041E00,1.028E00,8.929E-01,8.811E-01,XSE  920
     *     8.703E-01,8.593E-01,7.957E-01,7.851E-01,7.757E-01,7.661E-01, XSE  930
     *     7.168E-01,7.077E-01,6.998E-01,6.918E-01,6.367E-01,6.297E-01, XSE  940
     *     6.240E-01,6.181E-01,5.474E-01,5.434E-01,5.407E-01,5.379E-01, XSE  950
     *     0.0,7*0.0/                                                   XSE  960
      DATA X16 /4.336E00,4.294E00,4.255E00,4.214E00,3.744E00,3.708E00,  XSE  970
     *     3.675E00,3.642E00,3.432E00,3.401E00,3.372E00,3.344E00,       XSE  980
     *     3.193E00,3.167E00,3.143E00,3.120E00,2.960E00,2.939E00,       XSE  990
     *     2.923E00,2.906E00,2.706E00,2.695E00,2.687E00,2.679E00,0.0,7* XSE 1000
     *     0.0/                                                         XSE 1010
      DATA X17 /1.752E00,1.725E00,1.700E00,1.674E00,1.370E00,1.347E00,  XSE 1020
     *     1.326E00,1.304E00,1.169E00,1.149E00,1.131E00,1.112E00,       XSE 1030
     *     1.015E00,9.984E-01,9.836E-01,9.686E-01,8.654E-01,8.526E-01,  XSE 1040
     *     8.420E-01,8.312E-01,7.029E-01,6.958E-01,6.908E-01,6.857E-01, XSE 1050
     *     0.0,7*0.0/                                                   XSE 1060
      DATA X18 /9.292E-02,9.680E-02,1.005E-01,1.044E-01,1.095E-01,      XSE 1070
     *     1.141E-01,1.181E-01,1.224E-01,1.303E-01,1.352E-01,1.395E-01, XSE 1080
     *     1.440E-01,1.567E-01,1.618E-01,1.660E-01,1.704E-01,1.924E-01, XSE 1090
     *     1.970E-01,2.006E-01,2.044E-01,2.444E-01,2.477E-01,2.497E-01, XSE 1100
     *     2.518E-01,0.0,7*0.0/                                         XSE 1110
      DATA X19 /3.066E00,3.028E00,2.993E00,2.956E00,2.569E00,2.536E00,  XSE 1120
     *     2.506E00,2.475E00,2.295E00,2.266E00,2.239E00,2.212E00,       XSE 1130
     *     2.076E00,2.050E00,2.028E00,2.005E00,1.853E00,1.833E00,       XSE 1140
     *     1.817E00,1.800E00,1.603E00,1.592E00,1.584E00,1.576E00,0.0,7* XSE 1150
     *     0.0/                                                         XSE 1160
      DATA X20 /7.517E-02,7.822E-02,8.118E-02,8.430E-02,8.506E-02,      XSE 1170
     *     8.877E-02,9.214E-02,9.568E-02,1.006E-01,1.048E-01,1.085E-01, XSE 1180
     *     1.123E-01,1.219E-01,1.263E-01,1.300E-01,1.339E-01,1.519E-01, XSE 1190
     *     1.560E-01,1.593E-01,1.627E-01,1.976E-01,2.006E-01,2.025E-01, XSE 1200
     *     2.045E-01,0.0,7*0.0/                                         XSE 1210
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC38(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PWRUS; DOCUMENTED IN ORNL/TM-11018                                    
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,05/                                                   
      DATA DBUP /0.0, 2.0, 4.0, 6.0, 8.0, 27*0.0 /                              
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  160
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  170
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  180
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  190
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  200
      DATA X01 /5*2.071E+01, 0.0, 26*0.0/                                       
      DATA X02 /5*1.068E+01, 0.0, 26*0.0/                                       
      DATA X03 /5*4.752E+01, 0.0, 26*0.0/                                       
      DATA X04 /5*8.348E+00, 0.0, 26*0.0/                                       
      DATA X05 /5*8.872E-01, 0.0, 26*0.0/                                       
      DATA X06 /5*3.328E+01, 0.0, 26*0.0/                                       
      DATA X07 /5*3.483E+01, 0.0, 26*0.0/                                       
      DATA X08 /5*2.308E+00, 0.0, 26*0.0/                                       
      DATA X09 /6.909E+01, 6.431E+01, 5.993E+01, 5.585E+01, 5.205E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X10 /1.211E+02, 1.137E+02, 1.067E+02, 1.002E+02, 9.407E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X11 /2.228E+02, 1.762E+02, 1.394E+02, 1.102E+02, 8.717E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X12 /4.202E+01, 4.023E+01, 3.851E+01, 3.687E+01, 3.529E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X13 /1.259E+02, 1.214E+02, 1.170E+02, 1.128E+02, 1.087E+02,          
     *     0.0, 26*0.0 /                                                        
      DATA X14 /3.320E+01, 3.240E+01, 3.163E+01, 3.087E+01, 3.013E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X15 /9.570E+01, 9.028E+01, 8.517E+01, 8.035E+01, 7.581E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X16 /5*2.392E+01, 0.0, 26*0.0 /                                      
      DATA X17 /5*2.518E+00, 0.0, 26*0.0 /                                      
      DATA X18 /5*4.784E+01, 0.0, 26*0.0 /                                      
      DATA X19 /5*5.677E+00, 0.0, 26*0.0 /                                      
      DATA X20 /5*1.281E+01, 0.0, 26*0.0 /                                      
C                                                                               
      NPT = MPT                                                         XSE 1010
      NXSEC = MXSEC                                                     XSE 1020
      DO 10 I=1,30                                                      XSE 1030
         F(I,01) = X01(I)                                               XSE 1040
         F(I,02) = X02(I)                                               XSE 1050
         F(I,03) = X03(I)                                               XSE 1060
         F(I,04) = X04(I)                                               XSE 1070
         F(I,05) = X05(I)                                               XSE 1080
         F(I,06) = X06(I)                                               XSE 1090
         F(I,07) = X07(I)                                               XSE 1100
         F(I,08) = X08(I)                                               XSE 1110
         F(I,09) = X09(I)                                               XSE 1120
         F(I,10) = X10(I)                                               XSE 1130
         F(I,11) = X11(I)                                               XSE 1140
         F(I,12) = X12(I)                                               XSE 1150
         F(I,13) = X13(I)                                               XSE 1160
         F(I,14) = X14(I)                                               XSE 1170
         F(I,15) = X15(I)                                               XSE 1180
         F(I,16) = X16(I)                                               XSE 1190
         F(I,17) = X17(I)                                               XSE 1200
         F(I,18) = X18(I)                                               XSE 1210
         F(I,19) = X19(I)                                               XSE 1220
         F(I,20) = X20(I)                                               XSE 1230
         CBUP(I) = DBUP(I)                                              XSE 1240
   10 CONTINUE                                                          XSE 1250
      DO 30 I=1,20                                                      XSE 1260
         DO 20 J=1,7                                                    XSE 1270
            NFUD(J,I) = MFUD(J,I)                                       XSE 1280
   20    CONTINUE                                                       XSE 1290
   30 CONTINUE                                                          XSE 1300
      RETURN                                                            XSE 1310
      END                                                               XSE 1320
      SUBROUTINE XSEC39(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    PWRUE; DOCUMENTED IN ORNL/TM-11018
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,08/                                                   
      DATA DBUP /0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 24*0.0 /            
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  160
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  170
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  180
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  190
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  200
      DATA X01 /8*1.911E+01, 0.0, 23*0.0/                                       
      DATA X02 /8*8.910E+00, 0.0, 23*0.0/                                       
      DATA X03 /8*3.775E+01, 0.0, 23*0.0/                                       
      DATA X04 /8*8.295E+00, 0.0, 23*0.0/                                       
      DATA X05 /8*8.546E-01, 0.0, 23*0.0/                                       
      DATA X06 /8*3.031E+01, 0.0, 23*0.0/                                       
      DATA X07 /8*2.690E+01, 0.0, 23*0.0/                                       
      DATA X08 /8*2.119E+00, 0.0, 23*0.0/                                       
      DATA X09 /5.753E+01, 5.361E+01, 4.996E+01, 4.656E+01, 4.339E+01,          
     *     4.044E+01, 3.768E+01, 3.512E+01, 0.0, 23*0.0 /                       
      DATA X10 /9.922E+01, 9.315E+01, 8.745E+01, 8.210E+01, 7.707E+01,          
     *     7.236E+01, 6.793E+01, 6.377E+01, 0.0, 23*0.0 /                       
      DATA X11 /2.135E+02, 1.689E+02, 1.335E+02, 1.056E+02, 8.353E+01,          
     *     6.606E+01, 5.225E+01, 4.132E+01, 0.0, 23*0.0 /                       
      DATA X12 /3.363E+01, 3.219E+01, 3.082E+01, 2.951E+01, 2.825E+01,          
     *     2.704E+01, 2.589E+01, 2.478E+01, 0.0, 23*0.0 /                       
      DATA X13 /1.013E+02, 9.765E+01, 9.412E+01, 9.073E+01, 8.745E+01,          
     *     8.430E+01, 8.126E+01, 7.832E+01, 0.0, 23*0.0 /                       
      DATA X14 /3.290E+01, 3.211E+01, 3.134E+01, 3.059E+01, 2.985E+01,          
     *     2.914E+01, 2.844E+01, 2.776E+01, 0.0, 23*0.0 /                       
      DATA X15 /8.303E+01, 7.833E+01, 7.390E+01, 6.972E+01, 6.577E+01,          
     *     6.205E+01, 5.854E+01, 5.522E+01, 0.0, 23*0.0 /                       
      DATA X16 /8*2.076E+01, 0.0, 23*0.0 /                                      
      DATA X17 /8*2.434E+00, 0.0, 23*0.0 /                                      
      DATA X18 /8*4.625E+01, 0.0, 23*0.0 /                                      
      DATA X19 /8*5.439E+00, 0.0, 23*0.0 /                                      
      DATA X20 /8*1.271E+01, 0.0, 23*0.0 /                                      
C                                                                               
      NPT = MPT                                                         XSE 1010
      NXSEC = MXSEC                                                     XSE 1020
      DO 10 I=1,30                                                      XSE 1030
         F(I,01) = X01(I)                                               XSE 1040
         F(I,02) = X02(I)                                               XSE 1050
         F(I,03) = X03(I)                                               XSE 1060
         F(I,04) = X04(I)                                               XSE 1070
         F(I,05) = X05(I)                                               XSE 1080
         F(I,06) = X06(I)                                               XSE 1090
         F(I,07) = X07(I)                                               XSE 1100
         F(I,08) = X08(I)                                               XSE 1110
         F(I,09) = X09(I)                                               XSE 1120
         F(I,10) = X10(I)                                               XSE 1130
         F(I,11) = X11(I)                                               XSE 1140
         F(I,12) = X12(I)                                               XSE 1150
         F(I,13) = X13(I)                                               XSE 1160
         F(I,14) = X14(I)                                               XSE 1170
         F(I,15) = X15(I)                                               XSE 1180
         F(I,16) = X16(I)                                               XSE 1190
         F(I,17) = X17(I)                                               XSE 1200
         F(I,18) = X18(I)                                               XSE 1210
         F(I,19) = X19(I)                                               XSE 1220
         F(I,20) = X20(I)                                               XSE 1230
         CBUP(I) = DBUP(I)                                              XSE 1240
   10 CONTINUE                                                          XSE 1250
      DO 30 I=1,20                                                      XSE 1260
         DO 20 J=1,7                                                    XSE 1270
            NFUD(J,I) = MFUD(J,I)                                       XSE 1280
   20    CONTINUE                                                       XSE 1290
   30 CONTINUE                                                          XSE 1300
      RETURN                                                            XSE 1310
      END                                                               XSE 1320
      SUBROUTINE XSEC40(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    *BWRUS; DOCUMENTED IN ORNL/TM-11018
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,05/                                                   
      DATA DBUP /0.0, 2.0, 4.0, 6.0, 8.0, 27*0.0 /                              
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  170
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  180
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  190
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  200
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  210
      DATA X01 /5*1.940E+01, 0.0, 26*0.0/                                       
      DATA X02 /5*9.742E+00, 0.0, 26*0.0/                                       
      DATA X03 /5*4.186E+01, 0.0, 26*0.0/                                       
      DATA X04 /5*8.132E+00, 0.0, 26*0.0/                                       
      DATA X05 /5*8.025E-01, 0.0, 26*0.0/                                       
      DATA X06 /5*3.260E+01, 0.0, 26*0.0/                                       
      DATA X07 /5*2.978E+01, 0.0, 26*0.0/                                       
      DATA X08 /5*2.143E+00, 0.0, 26*0.0/                                       
      DATA X09 /6.845E+01, 6.379E+01, 5.944E+01, 5.540E+01, 5.163E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X10 /1.168E+02, 1.097E+02, 1.029E+02, 9.664E+01, 9.073E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X11 /2.227E+02, 1.761E+02, 1.393E+02, 1.102E+02, 8.713E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X12 /3.890E+01, 3.724E+01, 3.565E+01, 3.413E+01, 3.267E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X13 /1.153E+02, 1.111E+02, 1.071E+02, 1.033E+02, 9.954E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X14 /3.402E+01, 3.320E+01, 3.241E+01, 3.163E+01, 3.087E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X15 /9.599E+01, 9.056E+01, 8.543E+01, 8.060E+01, 7.604E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X16 /5*2.400E+01, 0.0, 26*0.0 /                                      
      DATA X17 /5*2.545E+00, 0.0, 26*0.0 /                                      
      DATA X18 /5*4.835E+01, 0.0, 26*0.0 /                                      
      DATA X19 /5*5.709E+00, 0.0, 26*0.0 /                                      
      DATA X20 /5*1.266E+01, 0.0, 26*0.0 /                                      
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC41(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    *BWRUS0; DOCUMENTED IN ORNL/TM-11018
C      (Pu-239 (n,fiss) @ 0 burnup corrected from 1.1245E+2 to 1.245E+2  
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,05/                                                   
      DATA DBUP /0.0, 2.0, 4.0, 6.0, 8.0, 27*0.0 /                              
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  170
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  180
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  190
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  200
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  210
      DATA X01 /5*2.028E+01, 0.0, 26*0.0/                                       
      DATA X02 /5*1.038E+01, 0.0, 26*0.0/                                       
      DATA X03 /5*4.558E+01, 0.0, 26*0.0/                                       
      DATA X04 /5*8.281E+00, 0.0, 26*0.0/                                       
      DATA X05 /5*8.262E-01, 0.0, 26*0.0/                                       
      DATA X06 /5*3.364E+01, 0.0, 26*0.0/                                       
      DATA X07 /5*3.288E+01, 0.0, 26*0.0/                                       
      DATA X08 /5*2.225E+00, 0.0, 26*0.0/                                       
      DATA X09 /7.234E+01, 6.741E+01, 6.282E+01, 5.855E+01, 5.456E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X10 /1.245E+02, 1.169E+02, 1.097E+02, 1.030E+02, 9.671E+01,         
     *     0.0, 26*0.0 /                                                        
      DATA X11 /2.277E+02, 1.801E+02, 1.424E+02, 1.126E+02, 8.909E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X12 /4.202E+01, 4.023E+01, 3.851E+01, 3.687E+01, 3.529E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X13 /1.247E+02, 1.202E+02, 1.159E+02, 1.117E+02, 1.077E+02,          
     *     0.0, 26*0.0 /                                                        
      DATA X14 /3.412E+01, 3.330E+01, 3.250E+01, 3.172E+01, 3.096E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X15 /1.001E+02, 9.443E+01, 8.909E+01, 8.405E+01, 7.929E+01,          
     *     0.0, 26*0.0 /                                                        
      DATA X16 /5*2.503E+01, 0.0, 26*0.0 /                                      
      DATA X17 /5*2.581E+00, 0.0, 26*0.0 /                                      
      DATA X18 /5*4.904E+01, 0.0, 26*0.0 /                                      
      DATA X19 /5*5.736E+00, 0.0, 26*0.0 /                                      
      DATA X20 /5*1.285E+01, 0.0, 26*0.0 /                                      
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE XSEC42(NFUD, F, CBUP, NXSEC, NPT)                      XSE   10
C                                                                       XSE   20
C    *BWRUE; DOCUMENTED IN ORNL/TM-11018
C     (Pu-241 (n,g) and (n,fiss) @ 14 burnup added by curve fit
C                                                                       XSE   40
      DIMENSION NFUD(7,20), F(32,20), CBUP(32)                          XSE   50
      DIMENSION X01(32), X02(32), X03(32), X04(32), X05(32), X06(32),   XSE   60
     *     X07(32), X08(32), X09(32), X10(32), X11(32), X12(32),        XSE   70
     *     X13(32), X14(32), X15(32), X16(32), X17(32), X18(32),        XSE   80
     *     X19(32), X20(32)                                             XSE   90
      DIMENSION MFUD(7,21), DBUP(32)                                    XSE  100
      DATA MXSEC, MPT /20,08/                                                   
      DATA DBUP /0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 24*0.0 /            
      DATA MFUD /922340,1,5*0,922350,1,5*0,922350,4,5*0,922360,1,5*0,   XSE  170
     *     922380,1,5*0,932370,1,5*0,942380,1,5*0,942380,4,5*0,942390,1,XSE  180
     *     5*0,942390,4,5*0,942400,1,5*0,942410,1,5*0,942410,4,5*0,     XSE  190
     *     942420,1,5*0,952410,1,5*0,952410,2,5*0,952430,1,5*0,952430,2,XSE  200
     *     5*0,962420,1,5*0,962440,1,5*0,7*0/                           XSE  210
      DATA X01 /8*1.800E+01, 0.0, 23*0.0/                                       
      DATA X02 /8*8.120E+00, 0.0, 23*0.0/                                       
      DATA X03 /8*3.301E+01, 0.0, 23*0.0/                                       
      DATA X04 /8*8.084E+00, 0.0, 23*0.0/                                       
      DATA X05 /8*7.730E-01, 0.0, 23*0.0/                                       
      DATA X06 /8*2.992E+01, 0.0, 23*0.0/                                       
      DATA X07 /8*2.293E+01, 0.0, 23*0.0/                                       
      DATA X08 /8*1.981E+00, 0.0, 23*0.0/                                       
      DATA X09 /5.697E+01, 5.309E+01, 4.948E+01, 4.611E+01, 4.297E+01,          
     *     4.004E+01, 3.732E+01, 3.477E+01, 0.0, 23*0.0 /                       
      DATA X10 /9.581E+01, 8.995E+01, 8.444E+01, 7.928E+01, 7.443E+01,          
     *     6.987E+01, 6.560E+01, 6.158E+01, 0.0, 23*0.0 /                       
      DATA X11 /2.137E+02, 1.690E+02, 1.337E+02, 1.057E+02, 8.361E+01,          
     *     6.613E+01, 5.230E+01, 4.136E+01, 0.0, 23*0.0 /                       
      DATA X12 /3.113E+01, 2.980E+01, 2.853E+01, 2.731E+01, 2.615E+01,          
     *     2.503E+01, 2.396E+01, 2.25E+01, 0.0, 23*0.0 /  ! 22.5 from curve fit
      DATA X13 /9.282E+01, 8.947E+01, 8.624E+01, 8.313E+01, 8.013E+01,          
     *     7.724E+01, 7.445E+01, 7.14E+01, 0.0, 23*0.0 /  ! 71.4 from curve fit
      DATA X14 /3.363E+01, 3.282E+01, 3.204E+01, 3.127E+01, 3.052E+01,          
     *     2.978E+01, 2.907E+01, 2.837E+01, 0.0, 23*0.0 /                       
      DATA X15 /8.356E+01, 7.883E+01, 7.437E+01, 7.016E+01, 6.619E+01,          
     *     6.244E+01, 5.891E+01, 5.558E+01, 0.0, 23*0.0 /                       
      DATA X16 /8*2.089E+01, 0.0, 23*0.0 /                                      
      DATA X17 /8*2.470E+00, 0.0, 23*0.0 /                                      
      DATA X18 /8*4.694E+01, 0.0, 23*0.0 /                                      
      DATA X19 /8*5.499E+00, 0.0, 23*0.0 /                                      
      DATA X20 /8*1.258E+01, 0.0, 23*0.0 /                                      
      NPT = MPT                                                         XSE 1220
      NXSEC = MXSEC                                                     XSE 1230
      DO 10 I=1,30                                                      XSE 1240
         F(I,01) = X01(I)                                               XSE 1250
         F(I,02) = X02(I)                                               XSE 1260
         F(I,03) = X03(I)                                               XSE 1270
         F(I,04) = X04(I)                                               XSE 1280
         F(I,05) = X05(I)                                               XSE 1290
         F(I,06) = X06(I)                                               XSE 1300
         F(I,07) = X07(I)                                               XSE 1310
         F(I,08) = X08(I)                                               XSE 1320
         F(I,09) = X09(I)                                               XSE 1330
         F(I,10) = X10(I)                                               XSE 1340
         F(I,11) = X11(I)                                               XSE 1350
         F(I,12) = X12(I)                                               XSE 1360
         F(I,13) = X13(I)                                               XSE 1370
         F(I,14) = X14(I)                                               XSE 1380
         F(I,15) = X15(I)                                               XSE 1390
         F(I,16) = X16(I)                                               XSE 1400
         F(I,17) = X17(I)                                               XSE 1410
         F(I,18) = X18(I)                                               XSE 1420
         F(I,19) = X19(I)                                               XSE 1430
         F(I,20) = X20(I)                                               XSE 1440
         CBUP(I) = DBUP(I)                                              XSE 1450
   10 CONTINUE                                                          XSE 1460
      DO 30 I=1,20                                                      XSE 1470
         DO 20 J=1,7                                                    XSE 1480
            NFUD(J,I) = MFUD(J,I)                                       XSE 1490
   20    CONTINUE                                                       XSE 1500
   30 CONTINUE                                                          XSE 1510
      RETURN                                                            XSE 1520
      END                                                               XSE 1530
      SUBROUTINE ERRSET(I,J,K,L)
      RETURN
      END
      SUBROUTINE Q105F(N)                                                       
      RETURN                                                                    
      END                                                                       
      SUBROUTINE IDENTIFY (IOUT, J)
C
C     IDENTIFY ADDS A LINE TO EACH OUTPUT PAGE THAT INCLUDES THE CODE NAME
C     AND VERSION NUMBER PLUS THE DATE AND TIME OF EXECUTION.
C
      CHARACTER*11 MMDDYY
      CHARACTER*8 HHMMSS
C
C     J IS USED SO THAT THE DATE AND TIME ARE ONLY CHECKED ONE TIME
C
      IF (J .EQ. 1) THEN
C
C     MAINFRAME USERS WILL NEED TO DETERMINE IF THE DATE AND TIME CALLS IN
C     FORTRAN ARE COMPATIBLE WITH THOSE IN THIS CODE.
C
      CALL DATE(MMDDYY)
      CALL TIME(HHMMSS)
      ELSE
      END IF
      WRITE (IOUT,10) MMDDYY, HHMMSS
   10 FORMAT (' ORIGEN2 V2.1 (8-1-91), ',
     1  'Run on ',a11,' at ',a8)
      RETURN
      END
      SUBROUTINE QAPAGE (IOUT)
C
C     QAPAGE PRINTS A HEADER PAGE (BLOCK LETTERS) AT THE BEGINNING OF OUTPUT ON
C     UNIT 6 AND 11 OUTPUT - THIS IDENTIFIES THE CODE VERSION.
C
      include 'HEADER.O2'
      RETURN
      END
