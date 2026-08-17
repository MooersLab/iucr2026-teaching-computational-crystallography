
 ***** CORRECT ***** (VERSION Apr 16, 2026  BUILT=20260616)  16-Aug-2026

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=    1 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    44.917    53.271    63.081 114.722 102.238  90.160 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=   1048576     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.300
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/blaine/3285CCinJNiucr26/demos/mpro_xds/frames/l6p17_09_001_?????.cbf
 DATA_RANGE=       1      60
 ROTATION_AXIS=  1.000000  0.000000  0.000000
 OSCILLATION_RANGE=  0.200000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  1.033200
 INCIDENT_BEAM_DIRECTION=  0.001746  0.001574  0.999997
 FRACTION_OF_POLARIZATION=   0.990
 POLARIZATION_PLANE_NORMAL=  0.000000  1.000000  0.000000
 AIR=   0.00038
 DETECTOR=PILATUS
 SILICON=  4.345246 SENSOR_THICKNESS=  0.450000
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463  NY=  2527    QX=  0.172000  QY=  0.172000
 ORGX=   1259.00  ORGY=   1240.00
 DETECTOR_DISTANCE=   200.000
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 BEAM_DIVERGENCE_E.S.D.=   0.079
 REFLECTING_RANGE_E.S.D.=   0.253
 MINIMUM_ZETA=     0.050
 MAXIMUM_ERROR_OF_SPOT_POSITION=   3.0
 MAXIMUM_ERROR_OF_SPINDLE_POSITION=   2.0
 MINIMUM_I/SIGMA=   3.0
 REFLECTIONS/CORRECTION_FACTOR=      50
 SNRC=  0.10 I/sigma for use in scaling factor determination
 STRICT_ABSORPTION_CORRECTION=FALSE 
 CORRECTIONS= DECAY MODULATION ABSORPTION                      
 TEST_RESOLUTION_RANGE=    10.000     5.000     MIN_RFL_Rmeas=    50     MAX_FAC_Rmeas=    2.00
 REFERENCE_DATA_SET=                                                  
 FIT_B-FACTOR_TO_REFERENCE_DATA_SET=FALSE 

   13562 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...      60
       0 OVERLOADED REFLECTIONS (IGNORED)
       1 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
   13561 REFLECTIONS ACCEPTED



 ******************************************************************************
                      AUTOMATIC SPACE GROUP ASSIGNMENT
 ******************************************************************************

 XDS adopts the following approach.
 (1) it looks for possible symmetries of the crystal lattice
 (2) it computes a redundancy independent R-factor for all enantiomorphous
     point groups compatible with the observed lattice symmetry.
 (3) it selects the group which explains the intensity data at an acceptable,
     redundancy-independent R-factor (Rmeas, Rrim) using a minimum number of
     unique reflections.

 This approach does not test for the presence of screw axes. Consequently,
 orthorhombic cell axes will be specified in increasing length (following
 conventions), despite the possibility that different assignments for the
 cell axes could become necessary for space groups P222(1) and P2(1)2(1)2
 containing one or two screw axes, respectively.

 The user can always override the automatic decisions by specifying the
 correct space group number and unit cell constants in XDS.INP and repeating
 the CORRECT step of XDS. This provides a simple way to rename orthorhombic
 cell constants if screw axes are present.

 In addition, the user has the option to specify in XDS.INP
 (a) a reference data set or
 (b) a reindexing transformation or
 (c) the three basis vectors (if known from processing a previous data set
 taken at the same crystal orientation in a multi-wavelength experiment).
 These features of XDS are useful for resolving the issue of alternative
 settings of polar or rhombohedral cells (like P4, P6, R3).

 =====> Specifications provided by the user in XDS.INP:
 SPACE_GROUP_NUMBER=    0 (0: unknown symmetry and unit cell)
 UNIT_CELL_CONSTANTS=     0.00     0.00     0.00   0.000   0.000   0.000
 UNIT_CELL_A-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_B-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_C-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 REIDX=  0  0  0  0  0  0  0  0  0  0  0  0 (all 0 : not specified)
 TEST_RESOLUTION_RANGE=   10.00    5.00 (Angstrom) for space group determination
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas              230
 NUMBER OF ACCEPTED UNIQUE REFLECTIONS FROM REFERENCE DATA SET       0


 *********** DETERMINATION OF LATTICE CHARACTER AND BRAVAIS LATTICE ***********

 The CHARACTER OF A LATTICE is defined by the metrical parameters of its
 reduced cell as described in the INTERNATIONAL TABLES FOR CRYSTALLOGRAPHY
 Volume A, p. 746 (KLUWER ACADEMIC PUBLISHERS, DORDRECHT/BOSTON/LONDON, 1989).
 Note that more than one lattice character may have the same BRAVAIS LATTICE.

 A lattice character is marked "*" to indicate a lattice consistent with the
 observed locations of the diffraction spots. These marked lattices must have
 low values for the QUALITY OF FIT and their implicated UNIT CELL CONSTANTS
 should not violate the ideal values by more than
 MAX_CELL_AXIS_ERROR=  0.03
 MAX_CELL_ANGLE_ERROR=           2.0 (Degrees)

 The REINDEXING TRANSFORMATION REIDX() consists of 12 integers that relate
 the original indices H,K,L from file INTEGRATE.HKL to the indices H',K',L'
 with respect to the new cell.
    H'=(REIDX(1)*H+REIDX( 2)*K+REIDX( 3)*L)/IDXV+REIDX( 4)
    K'=(REIDX(5)*H+REIDX( 6)*K+REIDX( 7)*L)/IDXV+REIDX( 8)
    L'=(REIDX(9)*H+REIDX(10)*K+REIDX(11)*L)/IDXV+REIDX(12)
 The value of the integer IDXV depends on the lattice type used for
 specifying reflections on file INTEGRATE.HKL; IDXV is
 1 for a primitive,
 2 for a face or body centred,
 3 for a rhombohedral,
 4 for a lattice centred on all faces.
 IDXV is set by XDS and cannot be input by the user.

  LATTICE-  BRAVAIS-   QUALITY  UNIT CELL CONSTANTS (ANGSTROEM & DEGREES)    REINDEXING TRANSFORMATION
 CHARACTER  LATTICE     OF FIT      a      b      c   alpha  beta gamma

 *  44        aP          0.0      44.9   53.3   63.1 114.7 102.2  90.2    1  0  0  0  0  1  0  0  0  0  1  0
 *  31        aP          1.7      44.9   53.3   63.1  65.3  77.8  90.2   -1  0  0  0  0 -1  0  0  0  0  1  0
 *  41        mC          5.0     114.6   53.3   44.9  90.2 103.6  89.7    0 -1 -2  0  0 -1  0  0 -1  0  0  0
    37        mC        102.8     124.6   44.9   53.3  90.2 115.1  81.6    1  0  2  0  1  0  0  0  0  1  0  0
    42        oI        106.2      44.9   53.3  112.8  90.2  99.2  90.2   -1  0  0  0  0 -1  0  0  1  1  2  0
    35        mP        150.5      53.3   44.9   63.1 102.2 114.7  90.2    0 -1  0  0 -1  0  0  0  0  0 -1  0
    40        oC        153.9      53.3  114.6   44.9  76.4  90.2  90.3    0 -1  0  0  0  1  2  0 -1  0  0  0
    30        mC        154.7      53.3  114.6   44.9  76.4  89.8  89.7    0 -1  0  0  0 -1 -2  0  1  0  0  0
    33        mP        349.9      44.9   53.3   63.1 114.7 102.2  90.2    1  0  0  0  0  1  0  0  0  0  1  0
    27        mC        354.5     115.7   44.9   63.3 102.3 118.3  67.0    1 -2  0  0  1  0  0  0  0  1  1  0
    39        mC        397.2     115.5   44.9   63.1 102.2 117.9  67.3   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    14        mC        402.7      69.6   69.8   63.1 100.5 117.2  80.3    1  1  0  0 -1  1  0  0  0  0  1  0
    10        mC        404.4      69.6   69.8   63.1  79.5 117.2  99.7   -1 -1  0  0 -1  1  0  0  0  0 -1  0
    15        tI        404.4      44.9   53.3  112.8  89.8  80.8  90.2    1  0  0  0  0  1  0  0  1  1  2  0
    43        mI        406.4      69.8  112.8   53.3  89.8 139.9  84.2   -1  1  0  0 -1 -1 -2  0  0 -1  0  0
    25        mC        430.1      63.3   98.1   44.9  97.7 102.3  79.4    0  1  1  0  0 -1  1  0  1  0  0  0
    23        oC        433.4      63.3   98.1   44.9  97.7 102.3  79.4    0  1  1  0  0 -1  1  0  1  0  0  0
    20        mC        435.0      98.1   63.3   44.9  77.7  97.7 100.6    0 -1  1  0  0 -1 -1  0  1  0  0  0
    22        hP        436.7      53.3   63.1   44.9 102.2  90.2 114.7    0  1  0  0  0  0  1  0  1  0  0  0
    36        oC        451.1      44.9  124.6   53.3  64.9  90.2  98.4    1  0  0  0 -1  0 -2  0  0  1  0  0
    28        mC        451.9      44.9  124.6   53.3  64.9  89.8  81.6    1  0  0  0  1  0  2  0  0 -1  0  0
    34        mP        497.1      44.9   63.1   53.3 114.7  90.2 102.2   -1  0  0  0  0  0 -1  0  0 -1  0  0
    32        oP        498.8      44.9   53.3   63.1 114.7 102.2  90.2    1  0  0  0  0  1  0  0  0  0  1  0
    29        mC        527.2      44.9  115.7   63.1  72.4 102.2  67.0   -1  0  0  0 -1  2  0  0  0  0 -1  0
    26        oF        577.8      44.9  115.7  124.6  63.4  98.4 113.0   -1  0  0  0  1 -2  0  0  1  0  2  0
    19        oI        637.3      44.9   63.3  102.3  85.3  71.9 102.3    1  0  0  0  0  1  1  0  1 -1  1  0
    24        hR        653.6     102.5   98.1   44.9  97.7  72.2 111.6    1  2  1  0  0 -1  1  0  1  0  0  0
     9        hR        656.0      44.9   69.8  169.2  82.2  91.6 130.2   -1  0  0  0  1 -1  0  0  1  1  3  0
    13        oC        700.4      69.6   69.8   63.1 100.5 117.2  80.3    1  1  0  0 -1  1  0  0  0  0  1  0
    11        tP        702.1      44.9   53.3   63.1 114.7 102.2  90.2    1  0  0  0  0  1  0  0  0  0  1  0
    38        oC        745.5      44.9  115.5   63.1  62.1 102.2 112.7   -1  0  0  0  1  2  0  0  0  0 -1  0
    21        tP        781.7      53.3   63.1   44.9 102.2  90.2 114.7    0  1  0  0  0  0  1  0  1  0  0  0
    17        mC        809.1      69.8   69.6   69.2  89.9 125.8  80.3    1 -1  0  0 -1 -1  0  0 -1  0 -1  0
     4        hR        832.8      69.8   84.8   69.3  95.5  89.8 101.9    1 -1  0  0 -1  0  1  0 -1 -1 -1  0
     2        hR        837.8      69.8   69.2  113.3 115.5  78.1 125.8   -1  1  0  0  1  0  1  0  1  1 -1  0
    18        tI        860.6      63.3  102.3   44.9  71.9  77.7  94.7    0  1  1  0 -1  1 -1  0 -1  0  0  0
     1        cF        936.0      85.1   69.3  102.3 106.0  95.9 109.7   -1  1  1  0 -1 -1 -1  0  1 -1  1  0
    12        hP        948.8      44.9   53.3   63.1 114.7 102.2  90.2    1  0  0  0  0  1  0  0  0  0  1  0
     3        cP        999.0      44.9   53.3   63.1 114.7 102.2  90.2    1  0  0  0  0  1  0  0  0  0  1  0
     5        cI        999.0      69.2   69.6   63.3  79.2  63.3  89.9    1  0  1  0  1  1  0  0  0  1  1  0
     6        tI        999.0      63.3   69.2   69.6  89.9  79.2  63.3    0  1  1  0  1  0  1  0  1  1  0  0
     7        tI        999.0      69.2   69.6   63.3  79.2  63.3  89.9    1  0  1  0  1  1  0  0  0  1  1  0
     8        oI        999.0      69.6   69.2   63.3  63.3  79.2  89.9   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
    16        oF        999.0      69.6   69.8  112.8  84.2  96.1  80.3   -1 -1  0  0  1 -1  0  0  1  1  2  0

 ****** LATTICE SYMMETRY IMPLICATED BY SPACE GROUP SYMMETRY ******

BRAVAIS-            POSSIBLE SPACE-GROUPS FOR PROTEIN CRYSTALS
  TYPE                     [SPACE GROUP NUMBER,SYMBOL]
  aP      [1,P1]
  mP      [3,P2] [4,P2(1)]
 mC,mI    [5,C2]
  oP      [16,P222] [17,P222(1)] [18,P2(1)2(1)2] [19,P2(1)2(1)2(1)]
  oC      [21,C222] [20,C222(1)]
  oF      [22,F222]
  oI      [23,I222] [24,I2(1)2(1)2(1)]
  tP      [75,P4] [76,P4(1)] [77,P4(2)] [78,P4(3)] [89,P422] [90,P42(1)2]
          [91,P4(1)22] [92,P4(1)2(1)2] [93,P4(2)22] [94,P4(2)2(1)2]
          [95,P4(3)22] [96,P4(3)2(1)2]
  tI      [79,I4] [80,I4(1)] [97,I422] [98,I4(1)22]
  hP      [143,P3] [144,P3(1)] [145,P3(2)] [149,P312] [150,P321] [151,P3(1)12]
          [152,P3(1)21] [153,P3(2)12] [154,P3(2)21] [168,P6] [169,P6(1)]
          [170,P6(5)] [171,P6(2)] [172,P6(4)] [173,P6(3)] [177,P622]
          [178,P6(1)22] [179,P6(5)22] [180,P6(2)22] [181,P6(4)22] [182,P6(3)22]
  hR      [146,R3] [155,R32]
  cP      [195,P23] [198,P2(1)3] [207,P432] [208,P4(2)32] [212,P4(3)32]
          [213,P4(1)32]
  cF      [196,F23] [209,F432] [210,F4(1)32]
  cI      [197,I23] [199,I2(1)3] [211,I432] [214,I4(1)32]


 ********************* SYMMETRY OF REFLECTION INTENSITIES *********************

 Analysis of symmetry related reflections is based on strong data within a
 TEST_RESOLUTION_RANGE=   10.00    5.00 Angstrom (defined in XDS.INP).

 The symmetry test is carried out for the enantiomorphous space groups above
 (omitting those with screw axes) that conform with the observed lattice
 symmetry. The test is quantified by the redundancy independent R-factor,
 Rmeas, defined by Diederichs & Karplus, Nature Struct. Biol. 4, 269-275 (1997).

 For computation of Rmeas a minimum number of symmetry-related reflections is
 required and set by the user (in XDS.INP) as the parameter
                      MIN_RFL_Rmeas=   50

 Low values for Rmeas should be obtained if the crystal possesses the symmetry
 elements of the tested space group. The decision whether the crystal's data
 satisfy the requirements of the tested space group is based on a constant set
 by the user (in XDS.INP) as the parameter
                      MAX_FAC_Rmeas=      2.
 which provides an upper limit for acceptable Rmeas as a multiple of the
 lowest Rmeas found for all tested space groups. If Rmeas is larger than this
 upper limit, the corresponding space group symmetry is considered violated
 by the observed data.

 XDS will select (*) that space group that requires a minimum number of unique
 reflections to account for the observed intensities at an acceptable Rmeas.

 UNIQUE   = number of unique reflections
 Rmeas    = redundancy independent R-factor (intensities)
 COMPARED = number of reflections used for calculating Rmeas

 SPACE-GROUP         UNIT CELL CONSTANTS            UNIQUE   Rmeas  COMPARED  LATTICE-
   NUMBER      a      b      c   alpha beta gamma                            CHARACTER

       1      44.9   53.3   63.1  65.3  77.8  90.2     217     3.9       13    31 aP
   *   5     114.6   53.3   44.9  90.0 103.6  90.0     130     4.2      100    41 mC
       1      44.9   53.3   63.1 114.7 102.2  90.2     217     3.9       13    44 aP


 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=    5
 UNIT_CELL_CONSTANTS=   114.60    53.27    44.92  90.000 103.571  90.000
 UNIT_CELL_A-AXIS=    -0.077   -54.934   100.575
 UNIT_CELL_B-AXIS=   -53.229     1.755     1.187
 UNIT_CELL_C-AXIS=    -1.643   -33.239   -30.166
 REFLECTIONS H,K,L FROM "INTEGRATE.HKL" WILL BE REINDEXED BY
 H' = (  0*H + -1*K + -2*L)/ 1 +  0
 K' = (  0*H + -1*K +  0*L)/ 1 +  0
 L' = ( -1*H +  0*K +  0*L)/ 1 +  0

 *****  4 EQUIVALENT POSITIONS IN SPACE GROUP #  5 *****

    If x',y',z' is an equivalent position to x,y,z, then
        x'=x*ML(1)+y*ML( 2)+z*ML( 3)+ML( 4)/12.0
        y'=x*ML(5)+y*ML( 6)+z*ML( 7)+ML( 8)/12.0
        z'=x*ML(9)+y*ML(10)+z*ML(11)+ML(12)/12.0

    #    1  2  3  4    5  6  7  8    9 10 11 12 
    1    1  0  0  0    0  1  0  0    0  0  1  0
    2   -1  0  0  0    0  1  0  0    0  0 -1  0
    3    1  0  0  6    0  1  0  6    0  0  1  0
    4   -1  0  0  6    0  1  0  6    0  0 -1  0

 ******************************************************************************
      MEAN DISCREPANCIES BETWEEN OBSERVED AND CALCULATED SPOT LOCATIONS
 ******************************************************************************

 The discrepancies in X- and Y-coordinates of the spots are depicted in the
 two images DX-CORRECTIONS.cbf and DY-CORRECTIONS.cbf for inspection with
 the XDS-Viewer.



 ******************************************************************************
  REFINEMENT OF DIFFRACTION PARAMETERS USING ALL IMAGES
 ******************************************************************************


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM      1146 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.33
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      5
 UNIT CELL PARAMETERS    113.635    53.148    44.670  90.000 103.138  90.000
 E.S.D. OF CELL PARAMETERS  3.4E-01 6.7E-02 1.2E-01 0.0E+00 1.2E-01 0.0E+00
 REC. CELL PARAMETERS   0.009037  0.018815  0.022988  90.000  76.862  90.000
 COORDINATES OF UNIT CELL A-AXIS    -0.044   -54.943    99.469
 COORDINATES OF UNIT CELL B-AXIS   -53.111     1.738     0.937
 COORDINATES OF UNIT CELL C-AXIS    -1.612   -33.142   -29.907
 CRYSTAL MOSAICITY (DEGREES)     0.253
 LAB COORDINATES OF ROTATION AXIS  1.000000  0.000000  0.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.003396  0.001612  0.967860
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1260.80   1241.79
 DETECTOR ORIGIN (PIXELS) AT                     1256.73   1239.86
 CRYSTAL TO DETECTOR DISTANCE (mm)       199.54
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=    5
 UNIT_CELL_CONSTANTS=   113.63    53.15    44.67  90.000 103.138  90.000


 ******************************************************************************
        MEAN INTENSITY AS FUNCTION OF SPINDLE POSITION WITHIN DATA IMAGE
 ******************************************************************************

 This statistics could serve as a diagnostic tool for detecting
 shutter problems (suggested by Kay Diederichs).
 Data are corrected for this effect if PATCH_SHUTTER_PROBLEM=TRUE in XDS.INP.
 ===> Selected : PATCH_SHUTTER_PROBLEM=FALSE

 INTERVAL = Angular interval in units of the oscillation range
            covered by a data image.
 NUMBER   = Number of reflections in each interval
 INTENSITY= Mean intensity of the reflections in each interval
 FACTOR   = Correction factor applied to intensities;
            INTENSITY(overall)/INTENSITY(interval)

   INTERVAL     NUMBER      INTENSITY      FACTOR
 -0.05  0.05      1370         28.776       0.910
  0.05  0.15      1319         27.015       0.969
  0.15  0.25      1391         26.988       0.970
  0.25  0.35      1341         28.714       0.912
  0.35  0.45      1382         23.214       1.128
  0.45  0.55      1310         28.830       0.908
  0.55  0.65      1362         24.859       1.053
  0.65  0.75      1386         25.253       1.037
  0.75  0.85      1338         24.618       1.063
  0.85  0.95      1362         23.664       1.106

 -0.05  0.95     13561         26.176       1.000


 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       40
 DEGREES OF FREEDOM OF CHI^2 FIT              1904.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.780
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     0.3 XMAX=    59.6 NXBIN=    2
 YMIN= 0.00104 YMAX= 0.59622 NYBIN=   20
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS       3818


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       81
 DEGREES OF FREEDOM OF CHI^2 FIT              1895.8
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.772
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    40.5 XMAX=  2388.4 NXBIN=    9
 YMIN=     4.9 YMAX=  2481.2 NYBIN=    9
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS       3818


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       26
 DEGREES OF FREEDOM OF CHI^2 FIT              1904.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.766
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     0.3 XMAX=    59.6 NXBIN=    2
 DETECTOR_SURFACE_POSITION=    1214    1243
 DETECTOR_SURFACE_POSITION=    1613    1664
 DETECTOR_SURFACE_POSITION=     816    1664
 DETECTOR_SURFACE_POSITION=     816     822
 DETECTOR_SURFACE_POSITION=    1613     822
 DETECTOR_SURFACE_POSITION=    2117    1637
 DETECTOR_SURFACE_POSITION=    1588    2195
 DETECTOR_SURFACE_POSITION=     841    2195
 DETECTOR_SURFACE_POSITION=     312    1637
 DETECTOR_SURFACE_POSITION=     312     849
 DETECTOR_SURFACE_POSITION=     841     291
 DETECTOR_SURFACE_POSITION=    1588     291
 DETECTOR_SURFACE_POSITION=    2117     849
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS       3818


 ******************************************************************************
    CORRECTION PARAMETERS FOR THE STANDARD ERROR OF REFLECTION INTENSITIES
 ******************************************************************************

 The variance v0(I) of the intensity I obtained from counting statistics is
 replaced by v(I)=a*(v0(I)+b*I^2). The model parameters a, b are chosen to
 minimize the discrepancies between v(I) and the variance estimated from
 sample statistics of symmetry related reflections. This model implicates
 an asymptotic limit ISa=1/SQRT(a*b) for the highest I/Sigma(I) that the
 experimental setup can produce (Diederichs (2010) Acta Cryst D66, 733-740).

     a        b          ISa
 6.340E-01  3.763E-03   20.47


 ******************************************************************************
      STANDARD ERROR OF REFLECTION INTENSITIES AS FUNCTION OF RESOLUTION
      FOR DATA SET  XDS_ASCII.HKL                                     
 ******************************************************************************

  I/Sigma  = mean intensity/Sigma of a reflection in shell
  Chi^2    = goodness of fit between sample variances of
             symmetry-related intensities and their errors
             (Chi^2 = 1 for perfect agreement)
  R-FACTOR
  observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
  expected = expected R-FACTOR derived from Sigma(I)

   NUMBER  = number of reflections in resolution shell
             used for calculation of R-FACTOR
  ACCEPTED = number of accepted reflections
  REJECTED = number of rejected reflections (MISFITS),
             recognized by comparison with symmetry-related
             reflections.

 RESOLUTION RANGE  I/Sigma  Chi^2  R-FACTOR  R-FACTOR  NUMBER ACCEPTED REJECTED
                                   observed  expected

   30.950  10.037    18.07   1.85      4.15      3.24      30      34       2
   10.037   7.291    16.38   0.56      2.44      3.14      45      54       0
    7.291   6.009    14.72   1.31      3.62      3.39      56      68       0
    6.009   5.229    12.44   1.16      3.61      3.64      69      82       0
    5.229   4.690    15.59   0.83      2.69      3.08      52      75       0
    4.690   4.290    14.20   1.02      3.06      3.28      64      86       0
    4.290   3.977    14.41   0.79      3.17      3.29      82     101       0
    3.977   3.724    12.83   0.92      3.06      3.67      84     115       0
    3.724   3.514    11.24   0.84      3.62      3.73      86     113       0
    3.514   3.336    10.35   1.07      3.59      4.14      94     121       6
    3.336   3.182     8.43   1.26      6.18      5.39     110     138       2
    3.182   3.048     7.09   1.37      7.25      6.28      96     130       0
    3.048   2.929     6.62   1.38      7.69      6.57      94     130       0
    2.929   2.824     5.96   1.07      9.14      8.11      92     136       2
    2.824   2.729     4.63   1.22     12.51     11.32      92     137       0
    2.729   2.643     5.08   0.92      8.82      9.04     106     141       2
    2.643   2.564     3.34   1.07     16.19     15.26     124     159       2
    2.564   2.493     3.35   1.42     20.34     16.22     124     169       6
    2.493   2.427     2.96   1.11     18.99     17.30     132     173       2
    2.427   2.366     3.42   1.21     15.87     13.74     130     167       6
    2.366   2.309     2.76   1.03     22.80     20.10     120     170       2
    2.309   2.256     2.05   1.33     27.35     24.64     152     195       2
    2.256   2.207     2.12   0.83     21.57     23.81     168     203       0
    2.207   2.160     1.98   0.93     29.72     29.02     132     170       2
    2.160   2.117     1.40   1.15     58.86     51.39     132     185       8
    2.117   2.076     1.41   0.97     49.18     47.57     138     188       2
    2.076   2.037     1.61   1.17     35.11     31.70     132     183       4
    2.037   2.001     1.37   0.81     31.51     35.36     130     185       2
    2.001   1.966     1.12   1.17     57.52     51.82     136     192       4
    1.966   1.933     0.98   0.90     47.16     51.64     158     214       4
    1.933   1.902     0.69   0.86     67.98     77.27     144     202       6
    1.902   1.872     0.68   1.05     88.26     89.59     168     213       0
    1.872   1.844     0.58   0.81     73.70     83.25     152     211       4
    1.844   1.816     0.65   1.00     99.62     91.14     156     203       6
    1.816   1.790     0.49   1.04    104.39     99.68     162     222       4
    1.790   1.765     0.39   0.93    111.96    115.70     178     228       2
    1.765   1.742     0.27   0.96    189.36    188.43     178     226       0
    1.742   1.719     0.31   1.22    211.92    196.99     190     239       2
    1.719   1.696     0.31   1.02    161.55    155.14     186     236       2
    1.696   1.675     0.33   0.93    211.50    212.13     166     220       0
    1.675   1.655     0.07   1.06    681.60    659.07     164     226       4
    1.655   1.635     0.10   1.00    726.49    741.23     148     218       2
    1.635   1.616     0.06   1.16   1504.80   1329.58     154     217       4
    1.616   1.597     0.25   0.93    181.47    186.30     170     232       4
    1.597   1.580     0.11   0.98    341.22    353.38     184     249       2
    1.580   1.562     0.13   0.68    338.22    420.77     180     245       0
    1.562   1.546     0.15   1.03    302.37    291.53     190     257       2
    1.546   1.530     0.07   1.03   1102.74   1089.03     184     243       2
    1.530   1.514    -0.03   1.03    -99.90    -99.90     186     255       4
    1.514   1.499     0.01   1.12    914.09    852.10     198     265       2
    1.499   1.484     0.10   1.07    375.06    353.11     220     283       2
    1.484   1.470     0.08   0.88   1045.77   1026.27     170     244       2
    1.470   1.456     0.09   1.00    814.95    787.98     220     274       4
    1.456   1.442    -0.04   1.03    -99.90    -99.90     202     266       4
    1.442   1.429     0.15   0.92    284.63    280.41     182     251       0
    1.429   1.416     0.09   1.01    433.47    421.50     210     274       4
    1.416   1.404     0.01   0.84    -99.90    -99.90     180     241       4
    1.404   1.392    -0.15   0.99    -99.90    -99.90     198     261       4
    1.392   1.380     0.05   0.74   2124.80   2441.49     196     273       2
    1.380   1.368    -0.08   0.89    -99.90    -99.90     198     269       0
    1.368   1.357    -0.04   0.70    -99.90    -99.90     202     271       4
    1.357   1.346    -0.07   0.89    -99.90    -99.90     214     289       2
    1.346   1.335     0.05   1.00   1004.97    966.60     204     278       0
    1.335   1.325     0.01   1.10   8284.24   7699.93     222     289       2
    1.325   1.315     0.04   1.19  10241.76   9091.77     198     272       4
    1.315   1.305     0.04   0.88    657.43    673.47     206     280       2
    1.305   1.295     0.05   1.17   1828.67   1628.47     206     273       0
   --------------------------------------------------------------------------
   30.950   1.295     1.72   1.01     19.26     19.00   10026   13409     152



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR VARIOUS SUBSETS OF INCLUDED DATA IMAGES
 ******************************************************************************

 Data set statistics is reported below several times, each with a different
 upper limit on the number of images included. This provides the user with
 the information for deciding which data images should be excluded from the
 final data set because of radiation damage or other defects. If the user
 decides for a subset of "good" images that differs from the specification
 as given by the input parameter (see DATA_RANGE= in XDS.INP) the CORRECT-
 step of XDS must be repeated with the new parameter values for DATA_RANGE=.


 ******************************************************************************
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1      60)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.86         548     323      2532       12.8%       3.1%      3.3%      430   18.72      4.3%    99.6*   -50    0.767       8
     2.74         958     604      4434       13.6%       5.4%      5.2%      708    9.87      7.7%    99.5*     0    0.000       0
     2.24        1238     770      5671       13.6%      17.7%     15.9%      936    3.97     25.1%    94.5*     0    0.000       0
     1.94        1427     901      6735       13.4%      37.5%     36.8%     1052    1.86     53.0%    82.6*     0    0.000       0
     1.74        1608    1000      7585       13.2%      96.0%     99.1%     1216    0.65    135.8%    41.8*     0    0.000       0
     1.59        1701    1071      8383       12.8%     276.6%    265.0%     1260    0.25    391.2%    14.8*     0    0.000       0
     1.47        1899    1193      9089       13.1%     624.2%    628.0%     1412    0.09    882.7%     6.7      0    0.000       0
     1.37        1938    1207      9755       12.4%    2637.8%   2646.0%     1462    0.03   3730.4%     6.1      0    0.000       0
     1.30        2081    1309     10344       12.7%     -99.9%    -99.9%     1544    0.00    -99.9%     2.7      0    0.000       0
    total       13398    8378     64528       13.0%      19.2%     19.0%    10020    2.13     27.2%    99.2*   -50    0.767       8


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     13561
 NUMBER OF REJECTED MISFITS                               152
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                        13409
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                   8386



 ******************************************************************************
    WILSON STATISTICS OF DATA SET  "XDS_ASCII.HKL"
 ******************************************************************************

 Data is divided into resolution shells and a straight line 
 A - 2*B*SS is fitted to log<I>, where
   RES    = mean resolution (Angstrom) in shell
   SS     = mean of (sin(THETA)/LAMBDA)**2 in shell
   <I>    = mean reflection intensity in shell
   BO     = (A - log<I>)/(2*SS)
    #     = number of reflections in resolution shell

   WILSON LINE (using all data) : A=   5.181 B=  27.953 CORRELATION=  0.98
      #      RES      SS        <I>       log(<I>)       BO
     108     7.042   0.005  1.2568E+02       4.834      34.4
     186     4.520   0.012  1.6124E+02       5.083       4.0
     268     3.537   0.020  1.0175E+02       4.622      14.0
     288     3.001   0.028  3.8144E+01       3.641      27.7
     334     2.654   0.035  1.9577E+01       2.974      31.1
     382     2.402   0.043  1.3294E+01       2.587      29.9
     404     2.211   0.051  8.9623E+00       2.193      29.2
     423     2.058   0.059  6.7405E+00       1.908      27.7
     455     1.933   0.067  3.7786E+00       1.329      28.8
     473     1.829   0.075  2.1526E+00       0.767      29.5
     500     1.740   0.083  1.3316E+00       0.286      29.6
     490     1.663   0.090  8.1187E-01      -0.208      29.8
     530     1.596   0.098  5.8613E-01      -0.534      29.1
     561     1.536   0.106  1.6928E-01      -1.776      32.8
     582     1.482   0.114  5.1791E-01      -0.658      25.6
     573     1.434   0.122  3.1798E-01      -1.146      26.0
     587     1.389   0.130 -5.1213E-02    -999.999   -1000.0
     617     1.349   0.137  5.2704E-02      -2.943      29.6
     615     1.313   0.145  1.5709E-01      -1.851      24.2


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

      14     7.042       1.7006       1.5755       1.1799    
      14     4.520      0.96656      0.55290      0.24655    
      17     3.537      0.46162      0.24603      0.10775    
      12     3.001      0.60677      0.35569      0.15861    
      17     2.654       2.5977       3.4637       3.6238    
      14     2.402      0.60914      0.31292      0.13004    
      15     2.211       5.3846       9.1600       12.285    
      15     2.058       1.5791       1.7183       1.5989    
      16     1.933       4.4979       7.6244       10.969    
      15     1.829       2.7241       2.9006       2.7420    
      16     1.740       6.1510       5.9961       11.606    
      11     1.663       19.787       4.3672       71.900    
      15     1.482       25.699       50.185       139.67    
      10     1.434       66.479       204.95       870.46    
      15     1.389       5.0448       1.0817       6.7243    
      17     1.349       2844.4      -17747.      0.18037E+07
      15     1.313       340.56       3685.1       36948.    
     279   overall       218.84      -834.75      0.11361E+06


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      94     7.042      0.97405       1.3546       2.1794    
     172     4.520      0.92944      0.78883      0.59435    
     251     3.537       1.3205       1.9950       3.0792    
     276     3.001       1.4180       2.5369       5.0029    
     317     2.654       1.1747       1.4585       1.7587    
     368     2.402       1.2214       1.5661       2.0452    
     389     2.211       1.3366       2.1312       3.6768    
     408     2.058       1.6753       3.1441       6.2672    
     439     1.933       1.6286       1.8901       2.4347    
     458     1.829       2.8550       4.5280       9.1914    
     484     1.740       4.5119       6.6245       16.857    
     479     1.663       11.521       17.174       84.644    
     515     1.596       19.885       39.611       263.82    
     545     1.536       208.01       262.19       29744.    
     567     1.482       24.947       36.844       397.45    
     563     1.434       65.950       32.208       2513.2    
     600     1.349       2783.8       4918.8      0.40428E+07
     600     1.313       315.55       80.333       56462.    
    8097   overall       252.84       397.77      0.30599E+06

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      14     7.042  0.58  0.83  0.86  0.91  0.82  0.76  0.84  0.80  0.76  0.84
      14     4.520  0.86  1.03  0.86  1.06  0.96  0.89  0.84  0.80  0.76  0.73
      17     3.537  1.42  1.19  1.27  1.12  1.24  1.15  1.08  1.12  1.07  1.03
      12     3.001  1.34  1.45  1.60  1.41  1.44  1.34  1.26  1.19  1.14  1.10
      17     2.654  1.18  1.19  0.99  0.87  0.79  0.73  0.89  0.94  0.99  0.95
      14     2.402  1.15  0.83  0.69  0.75  0.69  0.89  0.84  0.91  1.09  1.05
      15     2.211  0.54  0.77  0.64  0.56  0.51  0.83  0.89  0.85  0.81  0.78
      15     2.058  1.07  0.77  0.64  0.70  0.64  0.71  0.67  0.74  0.81  0.88
      16     1.933  0.50  0.36  0.30  0.26  0.24  0.22  0.31  0.30  0.29  0.27
      15     1.829  1.07  0.77  0.64  0.70  0.64  0.59  0.56  0.64  0.61  0.59
      16     1.740  1.01  0.72  0.75  0.66  0.84  0.78  0.73  0.70  0.67  0.73
      11     1.663  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
      15     1.482  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
      10     1.434  1.61  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
      15     1.389  3.22  2.31  1.92  1.69  1.54  1.43  1.34  1.27  1.22  1.17
      17     1.349  2.13  1.53  1.27  1.12  1.02  0.94  0.89  0.84  0.81  0.78
      15     1.313  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
     279   overall  1.18  0.97  0.86  0.80  0.76  0.75  0.73  0.72  0.71  0.70


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      94     7.042  0.78  1.12  1.19  1.06  1.03  1.06  1.01  1.00  0.97  0.98
     172     4.520  0.98  0.99  1.12  1.08  1.05  1.06  1.03  1.02  1.00  0.98
     251     3.537  1.47  1.23  1.17  1.10  1.02  1.06  1.03  1.06  1.04  1.04
     276     3.001  1.33  1.26  1.26  1.16  1.14  1.11  1.09  1.09  1.05  1.08
     317     2.654  1.62  1.27  1.28  1.22  1.17  1.12  1.10  1.08  1.09  1.06
     368     2.402  1.71  1.26  1.17  1.13  1.04  1.03  1.05  1.02  1.03  1.02
     389     2.211  2.40  1.47  1.30  1.17  1.12  1.08  1.08  1.10  1.07  1.07
     408     2.058  2.32  1.57  1.28  1.15  1.07  1.11  1.09  1.07  1.09  1.08
     439     1.933  2.97  1.73  1.34  1.15  1.04  1.00  0.95  0.95  0.94  0.92
     458     1.829  3.53  2.04  1.51  1.27  1.13  1.03  0.96  0.92  0.88  0.85
     484     1.740  4.10  2.27  1.61  1.29  1.12  1.03  0.96  0.92  0.87  0.84
     479     1.663  4.61  2.46  1.78  1.41  1.20  1.10  1.01  0.94  0.89  0.84
     515     1.596  4.67  2.47  1.76  1.41  1.20  1.06  0.97  0.89  0.83  0.80
     545     1.536  5.32  2.79  1.96  1.55  1.31  1.14  1.02  0.94  0.87  0.82
     567     1.482  4.71  2.50  1.78  1.41  1.19  1.05  0.96  0.89  0.84  0.80
     563     1.434  4.83  2.55  1.80  1.43  1.21  1.07  0.96  0.89  0.84  0.79
     600     1.349  5.22  2.75  1.92  1.52  1.27  1.11  1.00  0.91  0.85  0.80
     600     1.313  5.01  2.64  1.85  1.46  1.24  1.08  0.97  0.90  0.83  0.78
    8097   overall  3.45  1.96  1.48  1.23  1.08  1.00  0.93  0.89  0.86  0.83


 List of    500 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES        Z       Intensity    Sigma

    1   31   20    1.34   231.4      0.2439E+02  0.8817E+01 "alien"
    4   30   21    1.34   221.2      0.2332E+02  0.9115E+01 "alien"
   65   11   13    1.34   217.1      0.2288E+02  0.1070E+02 "alien"
   10   10   29    1.39   213.9     -0.2191E+02  0.8637E+01 "alien"
   72    2    8    1.40   180.2     -0.1845E+02  0.6999E+01 "alien"
   15   11   29    1.35   179.2      0.1889E+02  0.6745E+01 "alien"
    0   28   21    1.40   178.2     -0.1826E+02  0.7871E+01 "alien"
   63    9   13    1.39   177.7     -0.1820E+02  0.1030E+02 "alien"
   67    3   14    1.33   175.2      0.1846E+02  0.7437E+01 "alien"
    0   24   25    1.37   174.7      0.1842E+02  0.7920E+01 "alien"
   39   35   -4    1.35   174.3      0.1837E+02  0.6597E+01 "alien"
    2   28   23    1.33   173.2      0.1826E+02  0.6083E+01 "alien"
   10   32  -19    1.36   173.0      0.1824E+02  0.7545E+01 "alien"
   51   27    6    1.39   172.5     -0.1766E+02  0.8805E+01 "alien"
   71    1   10    1.37   171.9     -0.1760E+02  0.7259E+01 "alien"
   67   13    9    1.38   171.3     -0.1755E+02  0.7061E+01 "alien"
    6   30  -22    1.33   169.9      0.1791E+02  0.8118E+01 "alien"
   69   17    7    1.33   164.2      0.1731E+02  0.7481E+01 "alien"
   10   26   23    1.35   163.8      0.1727E+02  0.9260E+01 "alien"
   70   16    6    1.35   163.1      0.1719E+02  0.9908E+01 "alien"
   19    3   29    1.37   161.3     -0.1652E+02  0.6884E+01 "alien"
   50   30    1    1.38   158.8     -0.1626E+02  0.6080E+01 "alien"
   13   19   25    1.41   158.4     -0.1623E+02  0.6058E+01 "alien"
   57   25    5    1.38   158.4     -0.1622E+02  0.9409E+01 "alien"
   21    1   29    1.36   158.2      0.1668E+02  0.7254E+01 "alien"
    2   30  -19    1.41   156.6     -0.1604E+02  0.7543E+01 "alien"
   24    4   29    1.33   155.9      0.1643E+02  0.6991E+01 "alien"
   19    7   28    1.40   154.2     -0.1579E+02  0.6713E+01 "alien"
   69   17    6    1.35   152.2      0.1604E+02  0.9792E+01 "alien"
   62   20    6    1.40   151.9     -0.1555E+02  0.6774E+01 "alien"
   59   21    9    1.37   151.5     -0.1552E+02  0.6158E+01 "alien"
   20   10   27    1.41   150.6     -0.1543E+02  0.9530E+01 "alien"
   17    3   29    1.39   149.8     -0.1535E+02  0.6358E+01 "alien"
    3   31   19    1.36   148.0      0.1561E+02  0.8593E+01 "alien"
   16   18   25    1.41   146.7     -0.1502E+02  0.6199E+01 "alien"
   69    5   11    1.37   146.4     -0.1500E+02  0.7032E+01 "alien"
   44   32    2    1.37   144.6     -0.1481E+02  0.7700E+01 "alien"
   14   18   26    1.38   143.7     -0.1472E+02  0.6352E+01 "alien"
   49   31    0    1.37   141.5      0.1491E+02  0.6298E+01 "alien"
   64   22    4    1.36   141.4      0.1490E+02  0.9999E+01 "alien"
   25    1   28    1.37   141.1     -0.1445E+02  0.7073E+01 "alien"
    2   24   25    1.36   140.6      0.1482E+02  0.8366E+01 "alien"
   14   24   24    1.34   140.4      0.1480E+02  0.9247E+01 "alien"
   74   10    7    1.34   139.1      0.1466E+02  0.9647E+01 "alien"
    1   23   25    1.39   135.9     -0.1391E+02  0.7265E+01 "alien"
   56   26    7    1.34   135.5      0.1428E+02  0.9459E+01 "alien"
   13   15   27    1.40   135.4     -0.1387E+02  0.5733E+01 "alien"
   59   19   10    1.39   135.3     -0.1386E+02  0.6272E+01 "alien"
   14    4   29    1.41   135.3     -0.1385E+02  0.6291E+01 "alien"
   71   13    8    1.34   134.3      0.1415E+02  0.1005E+02 "alien"
   59   25    4    1.37   133.4      0.1406E+02  0.9989E+01 "alien"
   64   20    5    1.39   133.3     -0.1365E+02  0.9405E+01 "alien"
   18    8   28    1.40   132.4     -0.1356E+02  0.6449E+01 "alien"
   48   30    4    1.37   132.4     -0.1356E+02  0.8294E+01 "alien"
   59   27    3    1.33   131.1      0.1382E+02  0.7014E+01 "alien"
    7   27   21    1.40   131.0     -0.1342E+02  0.5833E+01 "alien"
   19   15   27    1.36   129.8      0.1368E+02  0.6689E+01 "alien"
   65    7   12    1.40   128.8     -0.1319E+02  0.1052E+02 "alien"
    2   26  -24    1.36   127.7      0.1346E+02  0.7624E+01 "alien"
   37   35   -4    1.36   127.4      0.1343E+02  0.6200E+01 "alien"
   15   21   24    1.40   127.1     -0.1302E+02  0.6028E+01 "alien"
   52   28    2    1.40   126.6     -0.1296E+02  0.8754E+01 "alien"
    7   23   25    1.36   126.4      0.1333E+02  0.6021E+01 "alien"
   67   11    9    1.40   125.4     -0.1285E+02  0.6987E+01 "alien"
   10   24   24    1.36   125.4      0.1321E+02  0.9312E+01 "alien"
   18   10   28    1.38   122.9     -0.1259E+02  0.6726E+01 "alien"
    4   28   21    1.39   122.1     -0.1250E+02  0.6013E+01 "alien"
   65    5   12    1.41   120.8     -0.1238E+02  0.9526E+01 "alien"
    1   27   24    1.33   120.5      0.1270E+02  0.6072E+01 "alien"
   51   29    1    1.39   119.4     -0.1223E+02  0.7959E+01 "alien"
   51   29    0    1.40   118.6     -0.1215E+02  0.8394E+01 "alien"
   71   15    6    1.35   118.1      0.1245E+02  0.9314E+01 "alien"
   56   22    8    1.41   118.1     -0.1209E+02  0.6091E+01 "alien"
   66   18    6    1.38   117.7     -0.1205E+02  0.6905E+01 "alien"
   55   27    4    1.37   117.6      0.1240E+02  0.6507E+01 "alien"
   13   19   26    1.37   116.9     -0.1197E+02  0.6293E+01 "alien"
   37   35   -6    1.36   116.5      0.1228E+02  0.6651E+01 "alien"
   22   14   27    1.35   116.5      0.1228E+02  0.1004E+02 "alien"
    9   27   21    1.39   116.0     -0.1188E+02  0.8611E+01 "alien"
   75    5    8    1.34   115.2      0.1215E+02  0.7670E+01 "alien"
   41   33    0    1.38   115.0     -0.1178E+02  0.6427E+01 "alien"
    1   31  -18    1.40   113.9     -0.1166E+02  0.8076E+01 "alien"
   10   34  -15    1.38   113.6     -0.1163E+02  0.6670E+01 "alien"
   63   23    5    1.34   113.0      0.1191E+02  0.7115E+01 "alien"
   10   16   28    1.36   112.9      0.1190E+02  0.8137E+01 "alien"
   25    3   28    1.37   112.3      0.1184E+02  0.7271E+01 "alien"
   69    9    9    1.39   112.1     -0.1148E+02  0.1041E+02 "alien"
   51   31    1    1.34   110.7      0.1167E+02  0.6227E+01 "alien"
   11   35  -16    1.33   110.7      0.1166E+02  0.6251E+01 "alien"
   21   15   27    1.34   110.4      0.1164E+02  0.9093E+01 "alien"
   68    4   13    1.34   109.8      0.1158E+02  0.7259E+01 "alien"
   17   15   27    1.37   109.7     -0.1123E+02  0.9295E+01 "alien"
   52   28    4    1.38   109.0     -0.1116E+02  0.6055E+01 "alien"
   52   26    6    1.40   108.8     -0.1115E+02  0.5808E+01 "alien"
    5   29   20    1.39   108.4     -0.1110E+02  0.8283E+01 "alien"
    9   23   25    1.35   108.4      0.1143E+02  0.6369E+01 "alien"
   66    4   13    1.37   108.3     -0.1109E+02  0.6753E+01 "alien"
   46   32    1    1.36   108.2      0.1141E+02  0.5701E+01 "alien"
    1   25   25    1.34   107.4      0.1132E+02  0.5532E+01 "alien"
   21    1   28    1.41   107.3     -0.1099E+02  0.6688E+01 "alien"
   22   10   28    1.35   107.1      0.1129E+02  0.6951E+01 "alien"
   21    3   28    1.40   106.9     -0.1095E+02  0.6706E+01 "alien"
   71    5    9    1.38   106.6     -0.1092E+02  0.7139E+01 "alien"
    8   18   28    1.34   106.1      0.1118E+02  0.8557E+01 "alien"
   60   22    5    1.40   105.7     -0.1083E+02  0.9567E+01 "alien"
   14   20   26    1.35   104.8      0.1105E+02  0.6426E+01 "alien"
    1   29   22    1.34   104.1      0.1097E+02  0.5969E+01 "alien"
   42   32   -1    1.41   103.8     -0.1063E+02  0.6738E+01 "alien"
   72    8    8    1.37   103.7      0.1094E+02  0.7584E+01 "alien"
   65    9   12    1.38   103.7     -0.1062E+02  0.8661E+01 "alien"
   13    9   30    1.34   103.6      0.1092E+02  0.6600E+01 "alien"
   57   27    2    1.37   102.6      0.1082E+02  0.9673E+01 "alien"
    4   20   26    1.40   101.8     -0.1043E+02  0.7399E+01 "alien"
   64   16   10    1.36   101.5      0.1070E+02  0.6944E+01 "alien"
    7   13   30    1.34   100.6      0.1060E+02  0.8235E+01 "alien"
   46   32    2    1.35   99.63      0.1050E+02  0.8081E+01 "alien"
   65   11   12    1.36   99.08      0.1044E+02  0.9997E+01 "alien"
   24   12   27    1.36   98.86      0.1042E+02  0.9121E+01 "alien"
   61   21    6    1.40   98.70     -0.1011E+02  0.9464E+01 "alien"
   16   18   27    1.33   98.65      0.1040E+02  0.6864E+01 "alien"
   34   34   -6    1.41   98.63     -0.1010E+02  0.5402E+01 "alien"
   21   11   27    1.39   97.82     -0.1002E+02  0.9170E+01 "alien"
   47   29    5    1.40   97.77     -0.1001E+02  0.8317E+01 "alien"
   12   34  -17    1.34   96.74      0.1020E+02  0.4928E+01 "alien"
    1   25   24    1.38   95.89     -0.9822E+01  0.8064E+01 "alien"
   21    9   28    1.37   95.76      0.1009E+02  0.6821E+01 "alien"
   41   35   -4    1.33   95.50      0.1007E+02  0.6881E+01 "alien"
   64   22    5    1.35   95.41      0.1006E+02  0.9910E+01 "alien"
   58   22    8    1.38   94.81     -0.9710E+01  0.8435E+01 "alien"
   60   16   11    1.40   94.67     -0.9696E+01  0.6247E+01 "alien"
   68   12   10    1.35   94.55      0.9967E+01  0.1028E+02 "alien"
    4   34  -17    1.34   94.24      0.9934E+01  0.7785E+01 "alien"
   20   14   27    1.36   94.07      0.9916E+01  0.9583E+01 "alien"
    2   28  -22    1.37   94.01     -0.9629E+01  0.5334E+01 "alien"
   64   20    8    1.34   93.97      0.9906E+01  0.6703E+01 "alien"
    3   29   22    1.34   93.48      0.9853E+01  0.6336E+01 "alien"
   17   13   27    1.40   93.41     -0.9568E+01  0.9186E+01 "alien"
   68    4   12    1.37   93.23      0.9828E+01  0.7405E+01 "alien"
   20   10   28    1.37   93.14      0.9818E+01  0.7027E+01 "alien"
   21   13   27    1.37   92.82      0.9784E+01  0.1012E+02 "alien"
   48   30    0    1.40   92.37     -0.9461E+01  0.8452E+01 "alien"
   19    1   30    1.34   92.00      0.9698E+01  0.7096E+01 "alien"
   66   20    5    1.36   91.93      0.9691E+01  0.9760E+01 "alien"
   52   30    3    1.34   91.88      0.9684E+01  0.8742E+01 "alien"
   63   19    8    1.37   91.77     -0.9400E+01  0.6602E+01 "alien"
   57   27    5    1.33   91.58      0.9653E+01  0.6578E+01 "alien"
   42   34   -4    1.35   91.36      0.9631E+01  0.7610E+01 "alien"
   14   18   27    1.35   91.29      0.9623E+01  0.9673E+01 "alien"
   56   28    3    1.34   90.82      0.9573E+01  0.6457E+01 "alien"
   15    7   30    1.34   90.78      0.9569E+01  0.6477E+01 "alien"
   13   11   28    1.41   90.45     -0.9264E+01  0.8508E+01 "alien"
   58   20   10    1.38   90.41     -0.9260E+01  0.9463E+01 "alien"
    3   29  -21    1.38   90.23     -0.9242E+01  0.8078E+01 "alien"
   23    9   27    1.39   89.87     -0.9205E+01  0.6724E+01 "alien"
   70    6   11    1.35   89.69      0.9454E+01  0.1028E+02 "alien"
   13   19   27    1.34   89.35      0.9418E+01  0.9324E+01 "alien"
   71    7   10    1.35   88.74      0.9354E+01  0.7529E+01 "alien"
   12   14   29    1.34   88.65      0.9345E+01  0.6505E+01 "alien"
   12    8   30    1.36   88.18      0.9295E+01  0.6545E+01 "alien"
   17   21   24    1.38   87.82     -0.8995E+01  0.8714E+01 "alien"
   69   15    6    1.38   86.91     -0.8902E+01  0.9364E+01 "alien"
   18    8   30    1.32   86.91      0.2730E+02  0.1093E+02 "alien"
   47   33   -2    1.34   86.24      0.9090E+01  0.8702E+01 "alien"
   57   21   10    1.38   85.78     -0.8786E+01  0.8996E+01 "alien"
   68    6   13    1.33   85.17      0.8978E+01  0.1003E+02 "alien"
    9   35  -16    1.33   85.09      0.8969E+01  0.6547E+01 "alien"
    4   30   20    1.36   84.68      0.8926E+01  0.8439E+01 "alien"
   11   35  -15    1.35   84.65      0.8923E+01  0.6554E+01 "alien"
   61   17   10    1.39   84.44     -0.8649E+01  0.6473E+01 "alien"
   62   18    9    1.38   83.78     -0.8581E+01  0.6547E+01 "alien"
   53   27    5    1.38   83.77     -0.8580E+01  0.5914E+01 "alien"
    5   19   27    1.38   83.76     -0.8579E+01  0.7865E+01 "alien"
    8   34  -16    1.36   83.22      0.8772E+01  0.7402E+01 "alien"
   58   24    8    1.34   82.78      0.8726E+01  0.8967E+01 "alien"
   36   34   -2    1.40   82.75     -0.8475E+01  0.5157E+01 "alien"
    8   16   28    1.37   82.63     -0.8464E+01  0.7529E+01 "alien"
   52   28    6    1.35   82.34      0.8679E+01  0.8802E+01 "alien"
   12   26   23    1.34   82.28      0.8673E+01  0.8595E+01 "alien"
   14   14   24    1.54   82.10      0.2780E+02  0.8496E+01 "alien"
   45   31    1    1.40   81.45     -0.8343E+01  0.5597E+01 "alien"
   49   31   -1    1.37   80.90     -0.8286E+01  0.8926E+01 "alien"
   44   34   -1    1.33   80.51      0.8486E+01  0.7453E+01 "alien"
   66   22    4    1.34   80.34      0.8469E+01  0.9549E+01 "alien"
   57   29    1    1.33   80.14      0.2518E+02  0.1001E+02 "alien"
   38   34   -3    1.39   79.81     -0.8175E+01  0.4245E+01 "alien"
   10    8   30    1.37   79.69      0.8400E+01  0.6577E+01 "alien"
   18   14   27    1.38   79.33     -0.8126E+01  0.9313E+01 "alien"
   64   18    6    1.41   78.74     -0.8065E+01  0.6680E+01 "alien"
   24    2   28    1.38   77.58     -0.7946E+01  0.6781E+01 "alien"
   17    5   30    1.34   77.45      0.8164E+01  0.6436E+01 "alien"
   68   16    6    1.38   76.99     -0.7886E+01  0.9132E+01 "alien"
   47   31   -1    1.39   76.93     -0.7880E+01  0.8617E+01 "alien"
   24    4   28    1.37   76.73     -0.7859E+01  0.6767E+01 "alien"
   15    1   30    1.36   76.41      0.8054E+01  0.6654E+01 "alien"
   70    4   10    1.38   76.02     -0.7787E+01  0.7181E+01 "alien"
   59   23    5    1.40   75.94     -0.7778E+01  0.6682E+01 "alien"
   68    2   11    1.40   75.61     -0.7744E+01  0.6926E+01 "alien"
   11   33  -17    1.37   75.25     -0.7708E+01  0.6532E+01 "alien"
    5   21   27    1.34   75.09      0.7915E+01  0.6018E+01 "alien"
   13   21   26    1.34   74.86      0.7891E+01  0.6660E+01 "alien"
   24   10   27    1.38   73.61     -0.7539E+01  0.6590E+01 "alien"
   25    9   27    1.38   73.57     -0.7535E+01  0.9110E+01 "alien"
   18   18   25    1.39   73.30     -0.7508E+01  0.6289E+01 "alien"
   20   16   27    1.34   72.83      0.7677E+01  0.9674E+01 "alien"
   12   20   25    1.40   72.62     -0.7438E+01  0.6033E+01 "alien"
   25    9   28    1.34   72.46      0.7638E+01  0.9559E+01 "alien"
   63   15   10    1.39   72.33     -0.7409E+01  0.6323E+01 "alien"
   68   22    4    1.31   72.13      0.2266E+02  0.1058E+02 "alien"
   65    3   14    1.36   71.92      0.7581E+01  0.1000E+02 "alien"
    9   17   28    1.35   71.43      0.7529E+01  0.8149E+01 "alien"
   59   23    7    1.37   70.98      0.7482E+01  0.8876E+01 "alien"
   54   26    6    1.38   70.51     -0.7222E+01  0.5997E+01 "alien"
    8   28   21    1.37   70.49     -0.7220E+01  0.8535E+01 "alien"
   19   35  -13    1.37   70.09     -0.7179E+01  0.4441E+01 "alien"
   12    8   29    1.40   70.09     -0.7179E+01  0.6144E+01 "alien"
    9   25   24    1.35   69.97      0.7376E+01  0.6393E+01 "alien"
   74    2    9    1.34   69.38      0.7313E+01  0.7515E+01 "alien"
    8   20   26    1.38   69.30     -0.7098E+01  0.5711E+01 "alien"
   62   18   10    1.36   69.30      0.7305E+01  0.6726E+01 "alien"
   61   27    2    1.32   69.10      0.2171E+02  0.1041E+02 "alien"
    9   27   24    1.30   68.87      0.2164E+02  0.9214E+01 "alien"
   60   26    4    1.33   68.19      0.7188E+01  0.7109E+01 "alien"
   65    5   13    1.38   67.88     -0.6952E+01  0.6824E+01 "alien"
    5   33  -16    1.39   67.56     -0.6920E+01  0.7164E+01 "alien"
   69    9   10    1.37   67.50      0.7116E+01  0.7415E+01 "alien"
   59   25    6    1.34   67.49      0.7114E+01  0.9651E+01 "alien"
   67    5   12    1.38   67.38     -0.6901E+01  0.7075E+01 "alien"
   19   11   28    1.36   66.93      0.7055E+01  0.6826E+01 "alien"
   40   34   -2    1.37   66.69      0.7030E+01  0.4556E+01 "alien"
   69    5   14    1.30   66.67      0.2095E+02  0.1197E+02 "alien"
   10   18   23    1.54   66.13      0.2239E+02  0.5842E+01 "alien"
    6   18   27    1.39   66.00     -0.6760E+01  0.5856E+01 "alien"
   54   28    3    1.37   65.87      0.6943E+01  0.6342E+01 "alien"
   63   19    6    1.40   65.87     -0.6746E+01  0.6862E+01 "alien"
    4   22   25    1.40   65.66     -0.6725E+01  0.7822E+01 "alien"
   63   21    6    1.37   65.45      0.6899E+01  0.6949E+01 "alien"
   58   24    4    1.40   65.33     -0.6692E+01  0.6490E+01 "alien"
   15   15   23    1.56   65.00      0.2200E+02  0.8766E+01 "alien"
   73    3    8    1.38   64.88     -0.6645E+01  0.7061E+01 "alien"
    3   29  -22    1.35   64.83      0.6834E+01  0.5731E+01 "alien"
   49   31    3    1.34   64.36      0.6785E+01  0.8703E+01 "alien"
    2   22   27    1.33   64.33      0.6781E+01  0.8648E+01 "alien"
   13    5   30    1.37   64.32      0.6780E+01  0.6444E+01 "alien"
   66   14   10    1.36   63.94      0.6739E+01  0.6967E+01 "alien"
   18    4   30    1.34   63.73      0.6717E+01  0.6429E+01 "alien"
   69    9    8    1.41   63.29     -0.6482E+01  0.7161E+01 "alien"
    6   28   22    1.35   63.15      0.6656E+01  0.6031E+01 "alien"
   15    3   30    1.36   62.85      0.6625E+01  0.6263E+01 "alien"
   64    8   13    1.38   62.34     -0.6385E+01  0.9635E+01 "alien"
   75    9    9    1.30   62.00      0.1948E+02  0.1252E+02 "alien"
   58   22    9    1.36   62.00      0.6535E+01  0.8710E+01 "alien"
   47   29    4    1.41   61.97     -0.6347E+01  0.8336E+01 "alien"
   18   10   29    1.34   61.61      0.6494E+01  0.6672E+01 "alien"
    2   30   21    1.34   61.08      0.6439E+01  0.8601E+01 "alien"
   18   20   24    1.40   61.08     -0.6256E+01  0.8766E+01 "alien"
    5   23   27    1.30   60.58      0.1903E+02  0.8837E+01 "alien"
   54   30    1    1.33   60.34      0.6361E+01  0.9639E+01 "alien"
   41   33   -4    1.39   60.22     -0.6168E+01  0.7561E+01 "alien"
   65   17    7    1.39   60.10     -0.6156E+01  0.6656E+01 "alien"
   69   15    8    1.34   60.02      0.6326E+01  0.1033E+02 "alien"
   70    4   12    1.33   59.92      0.6316E+01  0.7480E+01 "alien"
    8   32  -20    1.33   59.89      0.6313E+01  0.5386E+01 "alien"
   67   23    4    1.30   59.66      0.1874E+02  0.1052E+02 "alien"
   15    3   29    1.40   59.36     -0.6080E+01  0.6347E+01 "alien"
   73    9    8    1.34   59.30      0.6251E+01  0.7357E+01 "alien"
   15   17   27    1.36   58.94      0.6213E+01  0.9725E+01 "alien"
   69    3   11    1.38   58.63     -0.6005E+01  0.7021E+01 "alien"
   72    8    9    1.35   58.58      0.6175E+01  0.7624E+01 "alien"
   21    3   29    1.36   58.40      0.6156E+01  0.6879E+01 "alien"
    1   27  -24    1.34   58.13      0.6127E+01  0.7966E+01 "alien"
   26    4   27    1.40   57.91     -0.5932E+01  0.6571E+01 "alien"
   17   35  -15    1.35   57.78      0.6091E+01  0.5797E+01 "alien"
   70   14    7    1.36   57.70      0.6082E+01  0.7339E+01 "alien"
   59   27    4    1.32   57.52      0.1807E+02  0.6871E+01 "alien"
   65   15   11    1.34   57.46      0.6057E+01  0.1011E+02 "alien"
   17    7   30    1.33   57.35      0.6045E+01  0.7004E+01 "alien"
   11   13   28    1.40   57.33     -0.5872E+01  0.8153E+01 "alien"
    9   33  -16    1.40   57.00     -0.5838E+01  0.4831E+01 "alien"
   54   24    8    1.39   56.84     -0.5822E+01  0.8758E+01 "alien"
   14    8   29    1.38   56.49     -0.5786E+01  0.6250E+01 "alien"
   17   13   24    1.53   56.47      0.1912E+02  0.6233E+01 "alien"
   38   34   -5    1.38   55.96     -0.5732E+01  0.6570E+01 "alien"
    7   21   27    1.33   55.85      0.5887E+01  0.6066E+01 "alien"
   25    9   23    1.55   55.83      0.1890E+02  0.9039E+01 "alien"
   10   18   26    1.41   55.67     -0.5702E+01  0.5588E+01 "alien"
   66    2   13    1.38   55.47     -0.5682E+01  0.7016E+01 "alien"
    9   23   24    1.39   55.21     -0.5655E+01  0.5904E+01 "alien"
    4   28  -24    1.32   55.00      0.1728E+02  0.8590E+01 "alien"
   19   35  -14    1.36   54.75      0.5772E+01  0.4920E+01 "alien"
   46   24    7    1.52   54.58      0.1848E+02  0.9390E+01 "alien"
    8   32  -18    1.38   54.22     -0.5553E+01  0.5190E+01 "alien"
   23   11   29    1.30   54.12      0.1700E+02  0.7350E+01 "alien"
   54   26    5    1.39   53.53     -0.5483E+01  0.5883E+01 "alien"
   11   33  -18    1.35   53.28      0.5617E+01  0.7253E+01 "alien"
   67   15    9    1.35   53.12      0.5599E+01  0.6987E+01 "alien"
    0   30   20    1.37   53.10     -0.5439E+01  0.5495E+01 "alien"
   15   17   28    1.32   52.98      0.1665E+02  0.9773E+01 "alien"
   11   21   25    1.38   52.73     -0.5401E+01  0.5879E+01 "alien"
   67   19    5    1.36   52.73      0.5558E+01  0.9297E+01 "alien"
   65   15    8    1.40   52.52     -0.5379E+01  0.7074E+01 "alien"
   38   34   -1    1.38   52.51     -0.5378E+01  0.5553E+01 "alien"
   64   16   11    1.34   52.31      0.5514E+01  0.6938E+01 "alien"
   62   16   10    1.39   52.14     -0.5340E+01  0.6470E+01 "alien"
   46   32   -2    1.38   52.08     -0.5334E+01  0.8079E+01 "alien"
   64   24    6    1.30   52.05      0.1635E+02  0.1103E+02 "alien"
   13    1   30    1.38   52.00     -0.5326E+01  0.6350E+01 "alien"
   57    9   11    1.55   51.70      0.1750E+02  0.6599E+01 "alien"
   14    4   26    1.56   51.50      0.1744E+02  0.6391E+01 "alien"
   60   18    9    1.41   51.19     -0.5243E+01  0.6418E+01 "alien"
   66   16    7    1.39   50.83     -0.5206E+01  0.6998E+01 "alien"
    0   26   24    1.36   50.79      0.5353E+01  0.8417E+01 "alien"
    4   24   25    1.35   50.76      0.5350E+01  0.8501E+01 "alien"
   23    1   25    1.53   50.44      0.1708E+02  0.6800E+01 "alien"
   62   18    8    1.40   50.37     -0.5159E+01  0.6477E+01 "alien"
    8   28   20    1.40   50.04     -0.5125E+01  0.8213E+01 "alien"
   70   10   10    1.34   50.03      0.5274E+01  0.7337E+01 "alien"
   72    4   10    1.35   49.71      0.5239E+01  0.7248E+01 "alien"
    5   23   25    1.37   49.66     -0.5087E+01  0.8032E+01 "alien"
   10   20   26    1.37   49.51     -0.5072E+01  0.5770E+01 "alien"
   26    2   29    1.32   49.32      0.1550E+02  0.6796E+01 "alien"
   76    4    8    1.33   49.20      0.1546E+02  0.7752E+01 "alien"
   69    5   10    1.39   48.98     -0.5017E+01  0.1024E+02 "alien"
   63    7    8    1.54   48.49      0.1642E+02  0.9621E+01 "alien"
   14   10   25    1.55   48.36      0.1637E+02  0.6320E+01 "alien"
    6   26   24    1.34   48.34      0.5096E+01  0.5985E+01 "alien"
   65    9   11    1.41   48.15     -0.4931E+01  0.6876E+01 "alien"
   16   16   27    1.36   48.07      0.5067E+01  0.9083E+01 "alien"
   41   35   -3    1.33   47.80      0.5038E+01  0.6490E+01 "alien"
   12    4   30    1.38   47.65     -0.4880E+01  0.6132E+01 "alien"
   68    2   13    1.35   47.37      0.4993E+01  0.7095E+01 "alien"
    5   25   25    1.33   47.26      0.1485E+02  0.8659E+01 "alien"
   70    8    9    1.38   47.16     -0.4830E+01  0.1031E+02 "alien"
   23    1   28    1.39   46.85     -0.4799E+01  0.6984E+01 "alien"
   11   19   26    1.38   46.84     -0.4797E+01  0.5783E+01 "alien"
   16   18   26    1.37   46.62     -0.4775E+01  0.6566E+01 "alien"
   15   35  -14    1.37   46.39      0.4890E+01  0.5807E+01 "alien"
   66   12    9    1.41   46.29     -0.4742E+01  0.7072E+01 "alien"
   69   21    4    1.31   46.11      0.1449E+02  0.1002E+02 "alien"
    3   27   23    1.35   46.06      0.4855E+01  0.5999E+01 "alien"
   67    1   12    1.39   45.76     -0.4687E+01  0.6948E+01 "alien"
    6   34  -19    1.30   45.66      0.1434E+02  0.8748E+01 "alien"
   11   11   30    1.34   45.53      0.4799E+01  0.6640E+01 "alien"
   61   15   12    1.38   45.20     -0.4630E+01  0.9397E+01 "alien"
   69    7   12    1.34   45.14      0.4758E+01  0.7671E+01 "alien"
   24    8   28    1.35   44.89      0.4732E+01  0.6620E+01 "alien"
   12    4   31    1.34   44.89      0.4731E+01  0.6336E+01 "alien"
   10   22   27    1.30   44.85      0.1409E+02  0.6355E+01 "alien"
   57   21    9    1.40   44.73     -0.4582E+01  0.8465E+01 "alien"
    6   20   28    1.32   44.58      0.1401E+02  0.8719E+01 "alien"
   73    5    8    1.37   44.55     -0.4563E+01  0.7179E+01 "alien"
   19    5   28    1.41   44.42     -0.4550E+01  0.6814E+01 "alien"
   67    7   12    1.37   44.22      0.4661E+01  0.7441E+01 "alien"
   72   10    9    1.33   44.11      0.4650E+01  0.7606E+01 "alien"
   68   14    8    1.37   43.96     -0.4503E+01  0.7187E+01 "alien"
   69    7    9    1.40   43.61     -0.4467E+01  0.7081E+01 "alien"
   67    5   13    1.35   43.59      0.4595E+01  0.7371E+01 "alien"
   53   29    4    1.34   43.40      0.4575E+01  0.6220E+01 "alien"
   14    6   29    1.40   43.36     -0.4441E+01  0.6392E+01 "alien"
   12    6   30    1.37   43.32      0.4566E+01  0.6373E+01 "alien"
    6   26   23    1.37   43.25      0.4559E+01  0.9117E+01 "alien"
    6   34  -16    1.36   43.19      0.4552E+01  0.7279E+01 "alien"
   74    2   10    1.32   43.09      0.1354E+02  0.7671E+01 "alien"
   47   33    2    1.32   43.09      0.1354E+02  0.8154E+01 "alien"
   17   17   26    1.38   43.08     -0.4413E+01  0.6370E+01 "alien"
    5   19   29    1.30   43.01      0.1351E+02  0.5802E+01 "alien"
   48   32    1    1.34   43.01      0.4533E+01  0.5964E+01 "alien"
   33   35   -7    1.38   42.88     -0.4392E+01  0.5087E+01 "alien"
   74    6    9    1.33   42.49      0.1335E+02  0.1053E+02 "alien"
   69   17    8    1.32   42.39      0.1332E+02  0.1107E+02 "alien"
    0   24   27    1.30   42.29      0.1329E+02  0.8647E+01 "alien"
    6   20   27    1.35   42.24      0.4453E+01  0.5846E+01 "alien"
   70   18    7    1.31   42.23      0.1327E+02  0.1073E+02 "alien"
    8   22   25    1.38   41.92     -0.4294E+01  0.5737E+01 "alien"
   12   12   29    1.36   41.88      0.4414E+01  0.6505E+01 "alien"
   74   14    6    1.32   41.47      0.1303E+02  0.1042E+02 "alien"
   72    4    9    1.37   41.40     -0.4241E+01  0.7333E+01 "alien"
   68   10   13    1.31   41.34      0.1299E+02  0.7352E+01 "alien"
    4   20   27    1.36   41.34      0.4357E+01  0.8288E+01 "alien"
   19   17   22    1.53   41.30      0.1398E+02  0.8534E+01 "alien"
   13   21   25    1.37   41.16     -0.4215E+01  0.6411E+01 "alien"
    7   29   22    1.32   40.84      0.1283E+02  0.9234E+01 "alien"
   62   24    6    1.32   40.79      0.1282E+02  0.6837E+01 "alien"
   16    6   29    1.38   40.56     -0.4155E+01  0.6436E+01 "alien"
    6   26   18    1.54   40.56      0.1373E+02  0.5366E+01 "alien"
   24    4   24    1.56   40.48      0.1371E+02  0.6825E+01 "alien"
   37   29    2    1.55   40.47      0.1370E+02  0.7047E+01 "alien"
   21    5   28    1.39   40.27     -0.4125E+01  0.6715E+01 "alien"
   15   11   25    1.53   40.09      0.1357E+02  0.6280E+01 "alien"
   17   11   29    1.34   40.09      0.4225E+01  0.6839E+01 "alien"
   20   12   27    1.39   39.99     -0.4096E+01  0.6678E+01 "alien"
   65   11   11    1.39   39.98     -0.4095E+01  0.7156E+01 "alien"
   55   23    9    1.38   39.89     -0.4086E+01  0.9038E+01 "alien"
   57   27    3    1.36   39.59      0.4173E+01  0.6841E+01 "alien"
   68    6   12    1.36   39.41      0.4154E+01  0.1043E+02 "alien"
   74   12    7    1.32   39.30      0.1235E+02  0.7502E+01 "alien"
    9   11   31    1.31   39.08      0.1228E+02  0.6565E+01 "alien"
   62    6   10    1.51   39.02      0.1321E+02  0.7112E+01 "alien"
    3   29   23    1.31   39.00      0.1225E+02  0.6552E+01 "alien"
   53   27    3    1.40   38.89     -0.3983E+01  0.6352E+01 "alien"
   51   29    5    1.35   38.61      0.4069E+01  0.8997E+01 "alien"
   16   16   26    1.40   38.49     -0.3942E+01  0.6347E+01 "alien"
   66   22    6    1.31   38.24      0.1201E+02  0.7492E+01 "alien"
   46   30    4    1.39   38.22     -0.3914E+01  0.7870E+01 "alien"
   18   10   30    1.30   38.20      0.1200E+02  0.7245E+01 "alien"
   67   13   11    1.34   38.19      0.4025E+01  0.1016E+02 "alien"
    8   34  -17    1.34   38.16      0.4023E+01  0.7245E+01 "alien"
   50   26    7    1.41   38.15     -0.3907E+01  0.8177E+01 "alien"
    1   29  -21    1.38   38.01     -0.3893E+01  0.7985E+01 "alien"
   70    6   12    1.33   38.00      0.1194E+02  0.1055E+02 "alien"
   50   28    4    1.40   37.91     -0.3883E+01  0.5683E+01 "alien"
   23    3   28    1.38   37.66     -0.3857E+01  0.6754E+01 "alien"
   17    3   31    1.31   37.65      0.1183E+02  0.6984E+01 "alien"
   26    6   27    1.39   37.64     -0.3855E+01  0.8770E+01 "alien"
   66    8   12    1.37   37.64     -0.3855E+01  0.9543E+01 "alien"
   18   14   29    1.30   37.32      0.1173E+02  0.7012E+01 "alien"
   18   36  -16    1.30   37.29      0.1172E+02  0.6485E+01 "alien"
   73   15    6    1.32   37.17      0.1168E+02  0.1054E+02 "alien"
   11    5   30    1.38   37.05     -0.3795E+01  0.6296E+01 "alien"
   49   21    7    1.56   36.75      0.1244E+02  0.5930E+01 "alien"
   34   30   -2    1.56   36.70      0.1242E+02  0.6678E+01 "alien"
    3   33  -20    1.30   36.50      0.1147E+02  0.8862E+01 "alien"
    3   27  -23    1.37   36.49     -0.3737E+01  0.7590E+01 "alien"
   66   10   11    1.38   35.53     -0.3639E+01  0.1021E+02 "alien"
   17    1   31    1.31   35.43      0.1113E+02  0.7019E+01 "alien"
    5   29   23    1.30   35.37      0.1111E+02  0.9362E+01 "alien"
   11   23   26    1.31   35.35      0.1111E+02  0.9463E+01 "alien"
    1   29   20    1.40   35.30     -0.3616E+01  0.5663E+01 "alien"
   54   14   11    1.53   35.24      0.1193E+02  0.9575E+01 "alien"
   63    1   14    1.40   35.13     -0.3598E+01  0.7197E+01 "alien"
   61    5   10    1.54   35.03      0.1186E+02  0.6870E+01 "alien"
   25    1   29    1.33   35.03      0.3692E+01  0.7034E+01 "alien"
   63   21    8    1.33   35.00      0.3689E+01  0.6686E+01 "alien"
    1   29  -16    1.52   34.90      0.1182E+02  0.7173E+01 "alien"
   24    8   27    1.39   34.82     -0.3567E+01  0.6761E+01 "alien"
   59   13    7    1.56   34.67      0.1174E+02  0.9281E+01 "alien"
   64   14   11    1.37   34.65      0.3653E+01  0.6945E+01 "alien"
    4   30   22    1.31   34.64      0.1088E+02  0.9007E+01 "alien"
   16    4   31    1.31   34.55      0.1085E+02  0.6576E+01 "alien"
   41   33   -1    1.39   34.54     -0.3537E+01  0.4605E+01 "alien"
   13   25   23    1.36   34.43      0.3629E+01  0.8868E+01 "alien"
    6   24   26    1.31   34.30      0.1078E+02  0.6108E+01 "alien"
   11   11   28    1.42   34.29      0.2181E+02  0.6999E+01 "alien"
   10    6   30    1.38   34.21     -0.3504E+01  0.6119E+01 "alien"
   43   35   -5    1.32   34.12      0.1072E+02  0.7964E+01 "alien"
   11    5   32    1.30   34.00      0.1068E+02  0.6932E+01 "alien"
   13   15   29    1.32   33.99      0.1068E+02  0.6128E+01 "alien"
   13   25   24    1.33   33.95      0.1067E+02  0.9220E+01 "alien"
   61   15   11    1.40   33.81     -0.3463E+01  0.6486E+01 "alien"
    5   29   18    1.44   33.74      0.2146E+02  0.8467E+01 "alien"
   61   19    8    1.40   33.68     -0.3450E+01  0.6563E+01 "alien"
   56   30    3    1.30   33.42      0.1050E+02  0.1089E+02 "alien"
   12   20   22    1.52   33.36      0.1129E+02  0.5876E+01 "alien"
   25    3   24    1.55   33.27      0.1126E+02  0.6444E+01 "alien"
    8   30  -15    1.52   33.26      0.1126E+02  0.6352E+01 "alien"
   18    4   31    1.30   33.26      0.1045E+02  0.7409E+01 "alien"
    9   19   27    1.36   33.13      0.3492E+01  0.5890E+01 "alien"
   55   15   10    1.52   33.04      0.1119E+02  0.9514E+01 "alien"
   59   23    8    1.35   33.01      0.3480E+01  0.6619E+01 "alien"
   49   25    2    1.53   32.99      0.1117E+02  0.8694E+01 "alien"
   40   34   -1    1.37   32.92      0.3470E+01  0.6061E+01 "alien"
   58   12    9    1.55   32.88      0.1113E+02  0.8700E+01 "alien"
   10    4   32    1.31   32.70      0.1027E+02  0.6637E+01 "alien"
   72    2   10    1.35   32.67      0.3444E+01  0.7412E+01 "alien"
   21   15   22    1.55   32.59      0.1103E+02  0.8703E+01 "alien"
    5   25   24    1.36   32.58      0.3434E+01  0.5864E+01 "alien"
   70    8    8    1.40   32.53     -0.3332E+01  0.7316E+01 "alien"
   46   34   -2    1.32   32.31      0.1015E+02  0.5520E+01 "alien"
   69   17    9    1.30   32.24      0.1013E+02  0.9411E+01 "alien"
   62   16   11    1.37   32.21     -0.3299E+01  0.6358E+01 "alien"
    9   19   26    1.40   32.19     -0.3297E+01  0.5589E+01 "alien"
    3   29  -23    1.32   32.01      0.1006E+02  0.5885E+01 "alien"
    1   29  -22    1.35   31.97      0.3370E+01  0.5846E+01 "alien"
   48   30    3    1.38   31.95     -0.3273E+01  0.5801E+01 "alien"
    7   35  -17    1.31   31.84      0.1000E+02  0.7633E+01 "alien"
   61    9    8    1.56   31.83      0.1078E+02  0.6932E+01 "alien"
   40   28    4    1.52   31.62      0.1071E+02  0.7626E+01 "alien"
   51   17   10    1.55   31.62      0.1070E+02  0.6221E+01 "alien"
   70    4   13    1.31   31.59      0.9925E+01  0.7342E+01 "alien"
   74    2   11    1.30   31.51      0.9898E+01  0.7239E+01 "alien"
   63   15   11    1.37   31.49     -0.3226E+01  0.6531E+01 "alien"
   46   34   -3    1.32   31.46      0.9885E+01  0.8264E+01 "alien"
   55    9   13    1.52   31.46      0.1065E+02  0.9319E+01 "alien"
   72    2   11    1.33   31.44      0.3314E+01  0.6977E+01 "alien"
   66   14    9    1.38   31.25     -0.3201E+01  0.7013E+01 "alien"
   14   12   25    1.52   31.23      0.1057E+02  0.6068E+01 "alien"
   15   19   21    1.56   31.19      0.1056E+02  0.6048E+01 "alien"
   18    2   29    1.38   31.03     -0.3178E+01  0.6329E+01 "alien"
   12   18   28    1.32   31.00      0.9738E+01  0.6580E+01 "alien"
    9   15   29    1.34   30.92      0.3259E+01  0.8870E+01 "alien"
    0   32   18    1.37   30.73      0.3239E+01  0.7531E+01 "alien"
    6   28   21    1.38   30.72     -0.3147E+01  0.6049E+01 "alien"
    0   28   16    1.56   30.70      0.1039E+02  0.5000E+01 "alien"
   45   25    6    1.52   30.67      0.1038E+02  0.5816E+01 "alien"
   58    2   12    1.55   30.67      0.1038E+02  0.6854E+01 "alien"
   11    5   28    1.47   30.60      0.3169E+02  0.1006E+02 "alien"
   74   12    8    1.30   30.39      0.9547E+01  0.7700E+01 "alien"
   40   32   -3    1.43   30.16      0.1918E+02  0.7560E+01 "alien"
   23    7   26    1.45   30.09      0.1913E+02  0.6973E+01 "alien"
   16   36  -16    1.30   30.07      0.9449E+01  0.6286E+01 "alien"
   12    8   26    1.54   30.05      0.1017E+02  0.6101E+01 "alien"


 List of    500 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES        Z       Intensity    Sigma

   56   30    3    1.30   33.42      0.1050E+02  0.1089E+02
   75    9    9    1.30   62.00      0.1948E+02  0.1252E+02
   64   24    6    1.30   52.05      0.1635E+02  0.1103E+02
   23   11   29    1.30   54.12      0.1700E+02  0.7350E+01
   18   36  -16    1.30   37.29      0.1172E+02  0.6485E+01
   11    5   32    1.30   34.00      0.1068E+02  0.6932E+01
   18    4   31    1.30   33.26      0.1045E+02  0.7409E+01
   69   17    9    1.30   32.24      0.1013E+02  0.9411E+01
   69    5   14    1.30   66.67      0.2095E+02  0.1197E+02
    6   34  -19    1.30   45.66      0.1434E+02  0.8748E+01
    3   33  -20    1.30   36.50      0.1147E+02  0.8862E+01
   16   36  -16    1.30   30.07      0.9449E+01  0.6286E+01
   10   22   27    1.30   44.85      0.1409E+02  0.6355E+01
    5   29   23    1.30   35.37      0.1111E+02  0.9362E+01
   18   14   29    1.30   37.32      0.1173E+02  0.7012E+01
    9   27   24    1.30   68.87      0.2164E+02  0.9214E+01
   74    2   11    1.30   31.51      0.9898E+01  0.7239E+01
    0   24   27    1.30   42.29      0.1329E+02  0.8647E+01
   18   10   30    1.30   38.20      0.1200E+02  0.7245E+01
    5   19   29    1.30   43.01      0.1351E+02  0.5802E+01
   74   12    8    1.30   30.39      0.9547E+01  0.7700E+01
   67   23    4    1.30   59.66      0.1874E+02  0.1052E+02
    5   23   27    1.30   60.58      0.1903E+02  0.8837E+01
   70   18    7    1.31   42.23      0.1327E+02  0.1073E+02
   68   10   13    1.31   41.34      0.1299E+02  0.7352E+01
   10    4   32    1.31   32.70      0.1027E+02  0.6637E+01
    3   29   23    1.31   39.00      0.1225E+02  0.6552E+01
   17    3   31    1.31   37.65      0.1183E+02  0.6984E+01
    4   30   22    1.31   34.64      0.1088E+02  0.9007E+01
    9   11   31    1.31   39.08      0.1228E+02  0.6565E+01
   66   22    6    1.31   38.24      0.1201E+02  0.7492E+01
   68   22    4    1.31   72.13      0.2266E+02  0.1058E+02
   17    1   31    1.31   35.43      0.1113E+02  0.7019E+01
   11   23   26    1.31   35.35      0.1111E+02  0.9463E+01
   16    4   31    1.31   34.55      0.1085E+02  0.6576E+01
   70    4   13    1.31   31.59      0.9925E+01  0.7342E+01
   69   21    4    1.31   46.11      0.1449E+02  0.1002E+02
    7   35  -17    1.31   31.84      0.1000E+02  0.7633E+01
    6   24   26    1.31   34.30      0.1078E+02  0.6108E+01
   43   35   -5    1.32   34.12      0.1072E+02  0.7964E+01
   74   14    6    1.32   41.47      0.1303E+02  0.1042E+02
   69   17    8    1.32   42.39      0.1332E+02  0.1107E+02
   47   33    2    1.32   43.09      0.1354E+02  0.8154E+01
   18    8   30    1.32   86.91      0.2730E+02  0.1093E+02
   61   27    2    1.32   69.10      0.2171E+02  0.1041E+02
   46   34   -2    1.32   32.31      0.1015E+02  0.5520E+01
   73   15    6    1.32   37.17      0.1168E+02  0.1054E+02
    6   20   28    1.32   44.58      0.1401E+02  0.8719E+01
   15   17   28    1.32   52.98      0.1665E+02  0.9773E+01
    4   28  -24    1.32   55.00      0.1728E+02  0.8590E+01
   46   34   -3    1.32   31.46      0.9885E+01  0.8264E+01
   13   15   29    1.32   33.99      0.1068E+02  0.6128E+01
   59   27    4    1.32   57.52      0.1807E+02  0.6871E+01
   74   12    7    1.32   39.30      0.1235E+02  0.7502E+01
   12   18   28    1.32   31.00      0.9738E+01  0.6580E+01
   26    2   29    1.32   49.32      0.1550E+02  0.6796E+01
    7   29   22    1.32   40.84      0.1283E+02  0.9234E+01
   62   24    6    1.32   40.79      0.1282E+02  0.6837E+01
    3   29  -23    1.32   32.01      0.1006E+02  0.5885E+01
   74    2   10    1.32   43.09      0.1354E+02  0.7671E+01
   57   29    1    1.33   80.14      0.2518E+02  0.1001E+02
   13   25   24    1.33   33.95      0.1067E+02  0.9220E+01
   76    4    8    1.33   49.20      0.1546E+02  0.7752E+01
   70    6   12    1.33   38.00      0.1194E+02  0.1055E+02
   74    6    9    1.33   42.49      0.1335E+02  0.1053E+02
    5   25   25    1.33   47.26      0.1485E+02  0.8659E+01
    6   30  -22    1.33   169.9      0.1791E+02  0.8118E+01
    1   27   24    1.33   120.5      0.1270E+02  0.6072E+01
   17    7   30    1.33   57.35      0.6045E+01  0.7004E+01
   60   26    4    1.33   68.19      0.7188E+01  0.7109E+01
   25    1   29    1.33   35.03      0.3692E+01  0.7034E+01
   41   35   -3    1.33   47.80      0.5038E+01  0.6490E+01
    8   32  -20    1.33   59.89      0.6313E+01  0.5386E+01
   41   35   -4    1.33   95.50      0.1007E+02  0.6881E+01
   59   27    3    1.33   131.1      0.1382E+02  0.7014E+01
   44   34   -1    1.33   80.51      0.8486E+01  0.7453E+01
   57   27    5    1.33   91.58      0.9653E+01  0.6578E+01
   72    2   11    1.33   31.44      0.3314E+01  0.6977E+01
    7   21   27    1.33   55.85      0.5887E+01  0.6066E+01
   68    6   13    1.33   85.17      0.8978E+01  0.1003E+02
   24    4   29    1.33   155.9      0.1643E+02  0.6991E+01
   72   10    9    1.33   44.11      0.4650E+01  0.7606E+01
   69   17    7    1.33   164.2      0.1731E+02  0.7481E+01
   54   30    1    1.33   60.34      0.6361E+01  0.9639E+01
   11   35  -16    1.33   110.7      0.1166E+02  0.6251E+01
    2   22   27    1.33   64.33      0.6781E+01  0.8648E+01
    9   35  -16    1.33   85.09      0.8969E+01  0.6547E+01
   67    3   14    1.33   175.2      0.1846E+02  0.7437E+01
    2   28   23    1.33   173.2      0.1826E+02  0.6083E+01
   63   21    8    1.33   35.00      0.3689E+01  0.6686E+01
   16   18   27    1.33   98.65      0.1040E+02  0.6864E+01
   70    4   12    1.33   59.92      0.6316E+01  0.7480E+01
    4   30   21    1.34   221.2      0.2332E+02  0.9115E+01
    3   29   22    1.34   93.48      0.9853E+01  0.6336E+01
   66   22    4    1.34   80.34      0.8469E+01  0.9549E+01
   69    7   12    1.34   45.14      0.4758E+01  0.7671E+01
   12    4   31    1.34   44.89      0.4731E+01  0.6336E+01
    1   27  -24    1.34   58.13      0.6127E+01  0.7966E+01
    6   26   24    1.34   48.34      0.5096E+01  0.5985E+01
   25    9   28    1.34   72.46      0.7638E+01  0.9559E+01
   13   19   27    1.34   89.35      0.9418E+01  0.9324E+01
   67   13   11    1.34   38.19      0.4025E+01  0.1016E+02
   20   16   27    1.34   72.83      0.7677E+01  0.9674E+01
   11   11   30    1.34   45.53      0.4799E+01  0.6640E+01
   19    1   30    1.34   92.00      0.9698E+01  0.7096E+01
    7   13   30    1.34   100.6      0.1060E+02  0.8235E+01
   47   33   -2    1.34   86.24      0.9090E+01  0.8702E+01
   13   21   26    1.34   74.86      0.7891E+01  0.6660E+01
   71   13    8    1.34   134.3      0.1415E+02  0.1005E+02
   52   30    3    1.34   91.88      0.9684E+01  0.8742E+01
   18    4   30    1.34   63.73      0.6717E+01  0.6429E+01
   75    5    8    1.34   115.2      0.1215E+02  0.7670E+01
   12   14   29    1.34   88.65      0.9345E+01  0.6505E+01
   56   26    7    1.34   135.5      0.1428E+02  0.9459E+01
   64   20    8    1.34   93.97      0.9906E+01  0.6703E+01
   59   25    6    1.34   67.49      0.7114E+01  0.9651E+01
   51   31    1    1.34   110.7      0.1167E+02  0.6227E+01
   17   11   29    1.34   40.09      0.4225E+01  0.6839E+01
   74   10    7    1.34   139.1      0.1466E+02  0.9647E+01
   65   11   13    1.34   217.1      0.2288E+02  0.1070E+02
    4   34  -17    1.34   94.24      0.9934E+01  0.7785E+01
   12   26   23    1.34   82.28      0.8673E+01  0.8595E+01
   58   24    8    1.34   82.78      0.8726E+01  0.8967E+01
   17    5   30    1.34   77.45      0.8164E+01  0.6436E+01
    5   21   27    1.34   75.09      0.7915E+01  0.6018E+01
   65   15   11    1.34   57.46      0.6057E+01  0.1011E+02
   68    4   13    1.34   109.8      0.1158E+02  0.7259E+01
   14   24   24    1.34   140.4      0.1480E+02  0.9247E+01
    2   30   21    1.34   61.08      0.6439E+01  0.8601E+01
   48   32    1    1.34   43.01      0.4533E+01  0.5964E+01
    1   29   22    1.34   104.1      0.1097E+02  0.5969E+01
   53   29    4    1.34   43.40      0.4575E+01  0.6220E+01
   18   10   29    1.34   61.61      0.6494E+01  0.6672E+01
   64   16   11    1.34   52.31      0.5514E+01  0.6938E+01
    8   18   28    1.34   106.1      0.1118E+02  0.8557E+01
    9   15   29    1.34   30.92      0.3259E+01  0.8870E+01
   49   31    3    1.34   64.36      0.6785E+01  0.8703E+01
   13    9   30    1.34   103.6      0.1092E+02  0.6600E+01
   70   10   10    1.34   50.03      0.5274E+01  0.7337E+01
   12   34  -17    1.34   96.74      0.1020E+02  0.4928E+01
   21   15   27    1.34   110.4      0.1164E+02  0.9093E+01
   69   15    8    1.34   60.02      0.6326E+01  0.1033E+02
    8   34  -17    1.34   38.16      0.4023E+01  0.7245E+01
   63   23    5    1.34   113.0      0.1191E+02  0.7115E+01
   74    2    9    1.34   69.38      0.7313E+01  0.7515E+01
   56   28    3    1.34   90.82      0.9573E+01  0.6457E+01
    1   25   25    1.34   107.4      0.1132E+02  0.5532E+01
   15    7   30    1.34   90.78      0.9569E+01  0.6477E+01
    1   31   20    1.34   231.4      0.2439E+02  0.8817E+01
   73    9    8    1.34   59.30      0.6251E+01  0.7357E+01
    9   25   24    1.35   69.97      0.7376E+01  0.6393E+01
   17   35  -15    1.35   57.78      0.6091E+01  0.5797E+01
   68    2   13    1.35   47.37      0.4993E+01  0.7095E+01
   14   18   27    1.35   91.29      0.9623E+01  0.9673E+01
   39   35   -4    1.35   174.3      0.1837E+02  0.6597E+01
    1   29  -22    1.35   31.97      0.3370E+01  0.5846E+01
   71   15    6    1.35   118.1      0.1245E+02  0.9314E+01
   22   14   27    1.35   116.5      0.1228E+02  0.1004E+02
   72    8    9    1.35   58.58      0.6175E+01  0.7624E+01
   59   23    8    1.35   33.01      0.3480E+01  0.6619E+01
   72    4   10    1.35   49.71      0.5239E+01  0.7248E+01
   70    6   11    1.35   89.69      0.9454E+01  0.1028E+02
   70   16    6    1.35   163.1      0.1719E+02  0.9908E+01
   64   22    5    1.35   95.41      0.1006E+02  0.9910E+01
   14   20   26    1.35   104.8      0.1105E+02  0.6426E+01
   69   17    6    1.35   152.2      0.1604E+02  0.9792E+01
   71    7   10    1.35   88.74      0.9354E+01  0.7529E+01
   10   26   23    1.35   163.8      0.1727E+02  0.9260E+01
    6   28   22    1.35   63.15      0.6656E+01  0.6031E+01
    3   29  -22    1.35   64.83      0.6834E+01  0.5731E+01
   11   33  -18    1.35   53.28      0.5617E+01  0.7253E+01
   24    8   28    1.35   44.89      0.4732E+01  0.6620E+01
   51   29    5    1.35   38.61      0.4069E+01  0.8997E+01
   22   10   28    1.35   107.1      0.1129E+02  0.6951E+01
    9   17   28    1.35   71.43      0.7529E+01  0.8149E+01
   67    5   13    1.35   43.59      0.4595E+01  0.7371E+01
   11   35  -15    1.35   84.65      0.8923E+01  0.6554E+01
   52   28    6    1.35   82.34      0.8679E+01  0.8802E+01
   15   11   29    1.35   179.2      0.1889E+02  0.6745E+01
   46   32    2    1.35   99.63      0.1050E+02  0.8081E+01
   67   15    9    1.35   53.12      0.5599E+01  0.6987E+01
   42   34   -4    1.35   91.36      0.9631E+01  0.7610E+01
   68   12   10    1.35   94.55      0.9967E+01  0.1028E+02
   72    2   10    1.35   32.67      0.3444E+01  0.7412E+01
    9   23   25    1.35   108.4      0.1143E+02  0.6369E+01
    6   20   27    1.35   42.24      0.4453E+01  0.5846E+01
    4   24   25    1.35   50.76      0.5350E+01  0.8501E+01
    3   27   23    1.35   46.06      0.4855E+01  0.5999E+01
   24   12   27    1.36   98.86      0.1042E+02  0.9121E+01
   57   27    3    1.36   39.59      0.4173E+01  0.6841E+01
   15   17   27    1.36   58.94      0.6213E+01  0.9725E+01
   12    8   30    1.36   88.18      0.9295E+01  0.6545E+01
    0   26   24    1.36   50.79      0.5353E+01  0.8417E+01
   37   35   -6    1.36   116.5      0.1228E+02  0.6651E+01
   10   32  -19    1.36   173.0      0.1824E+02  0.7545E+01
   68    6   12    1.36   39.41      0.4154E+01  0.1043E+02
   13   25   23    1.36   34.43      0.3629E+01  0.8868E+01
   19   15   27    1.36   129.8      0.1368E+02  0.6689E+01
    9   19   27    1.36   33.13      0.3492E+01  0.5890E+01
   66   20    5    1.36   91.93      0.9691E+01  0.9760E+01
   70   14    7    1.36   57.70      0.6082E+01  0.7339E+01
   19   35  -14    1.36   54.75      0.5772E+01  0.4920E+01
   21    3   29    1.36   58.40      0.6156E+01  0.6879E+01
   66   14   10    1.36   63.94      0.6739E+01  0.6967E+01
   46   32    1    1.36   108.2      0.1141E+02  0.5701E+01
   12   12   29    1.36   41.88      0.4414E+01  0.6505E+01
   62   18   10    1.36   69.30      0.7305E+01  0.6726E+01
   37   35   -4    1.36   127.4      0.1343E+02  0.6200E+01
   15    3   30    1.36   62.85      0.6625E+01  0.6263E+01
   10   16   28    1.36   112.9      0.1190E+02  0.8137E+01
    2   26  -24    1.36   127.7      0.1346E+02  0.7624E+01
    2   24   25    1.36   140.6      0.1482E+02  0.8366E+01
   67   19    5    1.36   52.73      0.5558E+01  0.9297E+01
   10   24   24    1.36   125.4      0.1321E+02  0.9312E+01
    5   25   24    1.36   32.58      0.3434E+01  0.5864E+01
    4   30   20    1.36   84.68      0.8926E+01  0.8439E+01
   64   22    4    1.36   141.4      0.1490E+02  0.9999E+01
    6   34  -16    1.36   43.19      0.4552E+01  0.7279E+01
   64   16   10    1.36   101.5      0.1070E+02  0.6944E+01
    4   20   27    1.36   41.34      0.4357E+01  0.8288E+01
   21    1   29    1.36   158.2      0.1668E+02  0.7254E+01
   20   14   27    1.36   94.07      0.9916E+01  0.9583E+01
   58   22    9    1.36   62.00      0.6535E+01  0.8710E+01
    7   23   25    1.36   126.4      0.1333E+02  0.6021E+01
   65    3   14    1.36   71.92      0.7581E+01  0.1000E+02
    8   34  -16    1.36   83.22      0.8772E+01  0.7402E+01
   65   11   12    1.36   99.08      0.1044E+02  0.9997E+01
   16   16   27    1.36   48.07      0.5067E+01  0.9083E+01
    3   31   19    1.36   148.0      0.1561E+02  0.8593E+01
   19   11   28    1.36   66.93      0.7055E+01  0.6826E+01
   15    1   30    1.36   76.41      0.8054E+01  0.6654E+01
   40   34   -1    1.37   32.92      0.3470E+01  0.6061E+01
   59   23    7    1.37   70.98      0.7482E+01  0.8876E+01
   68    4   12    1.37   93.23      0.9828E+01  0.7405E+01
   49   31    0    1.37   141.5      0.1491E+02  0.6298E+01
   57   27    2    1.37   102.6      0.1082E+02  0.9673E+01
   59   25    4    1.37   133.4      0.1406E+02  0.9989E+01
   67    7   12    1.37   44.22      0.4661E+01  0.7441E+01
   15   35  -14    1.37   46.39      0.4890E+01  0.5807E+01
   54   28    3    1.37   65.87      0.6943E+01  0.6342E+01
   69    9   10    1.37   67.50      0.7116E+01  0.7415E+01
   20   10   28    1.37   93.14      0.9818E+01  0.7027E+01
   13    5   30    1.37   64.32      0.6780E+01  0.6444E+01
   63   21    6    1.37   65.45      0.6899E+01  0.6949E+01
    6   26   23    1.37   43.25      0.4559E+01  0.9117E+01
   10    8   30    1.37   79.69      0.8400E+01  0.6577E+01
   25    3   28    1.37   112.3      0.1184E+02  0.7271E+01
   21   13   27    1.37   92.82      0.9784E+01  0.1012E+02
   21    9   28    1.37   95.76      0.1009E+02  0.6821E+01
   40   34   -2    1.37   66.69      0.7030E+01  0.4556E+01
   55   27    4    1.37   117.6      0.1240E+02  0.6507E+01
   64   14   11    1.37   34.65      0.3653E+01  0.6945E+01
   72    8    8    1.37   103.7      0.1094E+02  0.7584E+01
    0   24   25    1.37   174.7      0.1842E+02  0.7920E+01
   12    6   30    1.37   43.32      0.4566E+01  0.6373E+01
    0   32   18    1.37   30.73      0.3239E+01  0.7531E+01
   48   30    4    1.37   132.4     -0.1356E+02  0.8294E+01
   59   21    9    1.37   151.5     -0.1552E+02  0.6158E+01
   63   19    8    1.37   91.77     -0.9400E+01  0.6602E+01
   69    5   11    1.37   146.4     -0.1500E+02  0.7032E+01
   63   15   11    1.37   31.49     -0.3226E+01  0.6531E+01
   72    4    9    1.37   41.40     -0.4241E+01  0.7333E+01
   73    5    8    1.37   44.55     -0.4563E+01  0.7179E+01
   62   16   11    1.37   32.21     -0.3299E+01  0.6358E+01
   68   14    8    1.37   43.96     -0.4503E+01  0.7187E+01
   49   31   -1    1.37   80.90     -0.8286E+01  0.8926E+01
   16   18   26    1.37   46.62     -0.4775E+01  0.6566E+01
   71    1   10    1.37   171.9     -0.1760E+02  0.7259E+01
    3   27  -23    1.37   36.49     -0.3737E+01  0.7590E+01
   17   15   27    1.37   109.7     -0.1123E+02  0.9295E+01
   25    1   28    1.37   141.1     -0.1445E+02  0.7073E+01
   66    4   13    1.37   108.3     -0.1109E+02  0.6753E+01
    8   28   21    1.37   70.49     -0.7220E+01  0.8535E+01
   10   20   26    1.37   49.51     -0.5072E+01  0.5770E+01
    5   23   25    1.37   49.66     -0.5087E+01  0.8032E+01
    8   16   28    1.37   82.63     -0.8464E+01  0.7529E+01
   44   32    2    1.37   144.6     -0.1481E+02  0.7700E+01
   13   21   25    1.37   41.16     -0.4215E+01  0.6411E+01
   13   19   26    1.37   116.9     -0.1197E+02  0.6293E+01
   11   33  -17    1.37   75.25     -0.7708E+01  0.6532E+01
   24    4   28    1.37   76.73     -0.7859E+01  0.6767E+01
    0   30   20    1.37   53.10     -0.5439E+01  0.5495E+01
    2   28  -22    1.37   94.01     -0.9629E+01  0.5334E+01
   66    8   12    1.37   37.64     -0.3855E+01  0.9543E+01
   19   35  -13    1.37   70.09     -0.7179E+01  0.4441E+01
   19    3   29    1.37   161.3     -0.1652E+02  0.6884E+01
    1   29  -21    1.38   38.01     -0.3893E+01  0.7985E+01
   24   10   27    1.38   73.61     -0.7539E+01  0.6590E+01
   61   15   12    1.38   45.20     -0.4630E+01  0.9397E+01
   46   32   -2    1.38   52.08     -0.5334E+01  0.8079E+01
   25    9   27    1.38   73.57     -0.7535E+01  0.9110E+01
   57   21   10    1.38   85.78     -0.8786E+01  0.8996E+01
   50   30    1    1.38   158.8     -0.1626E+02  0.6080E+01
    1   25   24    1.38   95.89     -0.9822E+01  0.8064E+01
    5   19   27    1.38   83.76     -0.8579E+01  0.7865E+01
   69    3   11    1.38   58.63     -0.6005E+01  0.7021E+01
   66    2   13    1.38   55.47     -0.5682E+01  0.7016E+01
   67    5   12    1.38   67.38     -0.6901E+01  0.7075E+01
   73    3    8    1.38   64.88     -0.6645E+01  0.7061E+01
   54   26    6    1.38   70.51     -0.7222E+01  0.5997E+01
   12    4   30    1.38   47.65     -0.4880E+01  0.6132E+01
   18   14   27    1.38   79.33     -0.8126E+01  0.9313E+01
   66   18    6    1.38   117.7     -0.1205E+02  0.6905E+01
   57   25    5    1.38   158.4     -0.1622E+02  0.9409E+01
   33   35   -7    1.38   42.88     -0.4392E+01  0.5087E+01
   13    1   30    1.38   52.00     -0.5326E+01  0.6350E+01
   69   15    6    1.38   86.91     -0.8902E+01  0.9364E+01
   52   28    4    1.38   109.0     -0.1116E+02  0.6055E+01
   24    2   28    1.38   77.58     -0.7946E+01  0.6781E+01
   67   13    9    1.38   171.3     -0.1755E+02  0.7061E+01
   68   16    6    1.38   76.99     -0.7886E+01  0.9132E+01
    8   22   25    1.38   41.92     -0.4294E+01  0.5737E+01
   53   27    5    1.38   83.77     -0.8580E+01  0.5914E+01
    3   29  -21    1.38   90.23     -0.9242E+01  0.8078E+01
   11    5   30    1.38   37.05     -0.3795E+01  0.6296E+01
    6   28   21    1.38   30.72     -0.3147E+01  0.6049E+01
    8   32  -18    1.38   54.22     -0.5553E+01  0.5190E+01
   48   30    3    1.38   31.95     -0.3273E+01  0.5801E+01
   70    8    9    1.38   47.16     -0.4830E+01  0.1031E+02
   64    8   13    1.38   62.34     -0.6385E+01  0.9635E+01
   17   17   26    1.38   43.08     -0.4413E+01  0.6370E+01
   10    6   30    1.38   34.21     -0.3504E+01  0.6119E+01
   70    4   10    1.38   76.02     -0.7787E+01  0.7181E+01
   62   18    9    1.38   83.78     -0.8581E+01  0.6547E+01
   38   34   -1    1.38   52.51     -0.5378E+01  0.5553E+01
   66   14    9    1.38   31.25     -0.3201E+01  0.7013E+01
   58   20   10    1.38   90.41     -0.9260E+01  0.9463E+01
   55   23    9    1.38   39.89     -0.4086E+01  0.9038E+01
   18   10   28    1.38   122.9     -0.1259E+02  0.6726E+01
    8   20   26    1.38   69.30     -0.7098E+01  0.5711E+01
   58   22    8    1.38   94.81     -0.9710E+01  0.8435E+01
   71    5    9    1.38   106.6     -0.1092E+02  0.7139E+01
   66   10   11    1.38   35.53     -0.3639E+01  0.1021E+02
   65    9   12    1.38   103.7     -0.1062E+02  0.8661E+01
   65    5   13    1.38   67.88     -0.6952E+01  0.6824E+01
   41   33    0    1.38   115.0     -0.1178E+02  0.6427E+01
   17   21   24    1.38   87.82     -0.8995E+01  0.8714E+01
   38   34   -5    1.38   55.96     -0.5732E+01  0.6570E+01
   14   18   26    1.38   143.7     -0.1472E+02  0.6352E+01
   11   21   25    1.38   52.73     -0.5401E+01  0.5879E+01
   10   34  -15    1.38   113.6     -0.1163E+02  0.6670E+01
   11   19   26    1.38   46.84     -0.4797E+01  0.5783E+01
   16    6   29    1.38   40.56     -0.4155E+01  0.6436E+01
   18    2   29    1.38   31.03     -0.3178E+01  0.6329E+01
   23    3   28    1.38   37.66     -0.3857E+01  0.6754E+01
   14    8   29    1.38   56.49     -0.5786E+01  0.6250E+01
   38   34   -3    1.39   79.81     -0.8175E+01  0.4245E+01
   59   19   10    1.39   135.3     -0.1386E+02  0.6272E+01
    5   29   20    1.39   108.4     -0.1110E+02  0.8283E+01
    1   23   25    1.39   135.9     -0.1391E+02  0.7265E+01
   20   12   27    1.39   39.99     -0.4096E+01  0.6678E+01
    4   28   21    1.39   122.1     -0.1250E+02  0.6013E+01
   64   20    5    1.39   133.3     -0.1365E+02  0.9405E+01
   26    6   27    1.39   37.64     -0.3855E+01  0.8770E+01
   41   33   -1    1.39   34.54     -0.3537E+01  0.4605E+01
   63    9   13    1.39   177.7     -0.1820E+02  0.1030E+02
   65   11   11    1.39   39.98     -0.4095E+01  0.7156E+01
   67    1   12    1.39   45.76     -0.4687E+01  0.6948E+01
   69    9    9    1.39   112.1     -0.1148E+02  0.1041E+02
    9   23   24    1.39   55.21     -0.5655E+01  0.5904E+01
   23    1   28    1.39   46.85     -0.4799E+01  0.6984E+01
   51   27    6    1.39   172.5     -0.1766E+02  0.8805E+01
    6   18   27    1.39   66.00     -0.6760E+01  0.5856E+01
   17    3   29    1.39   149.8     -0.1535E+02  0.6358E+01
   21   11   27    1.39   97.82     -0.1002E+02  0.9170E+01
   46   30    4    1.39   38.22     -0.3914E+01  0.7870E+01
   47   31   -1    1.39   76.93     -0.7880E+01  0.8617E+01
   54   24    8    1.39   56.84     -0.5822E+01  0.8758E+01
   61   17   10    1.39   84.44     -0.8649E+01  0.6473E+01
   65   17    7    1.39   60.10     -0.6156E+01  0.6656E+01
   54   26    5    1.39   53.53     -0.5483E+01  0.5883E+01
   24    8   27    1.39   34.82     -0.3567E+01  0.6761E+01
   63   15   10    1.39   72.33     -0.7409E+01  0.6323E+01
   62   16   10    1.39   52.14     -0.5340E+01  0.6470E+01
   66   16    7    1.39   50.83     -0.5206E+01  0.6998E+01
   41   33   -4    1.39   60.22     -0.6168E+01  0.7561E+01
   23    9   27    1.39   89.87     -0.9205E+01  0.6724E+01
   69    5   10    1.39   48.98     -0.5017E+01  0.1024E+02
    9   27   21    1.39   116.0     -0.1188E+02  0.8611E+01
   51   29    1    1.39   119.4     -0.1223E+02  0.7959E+01
   10   10   29    1.39   213.9     -0.2191E+02  0.8637E+01
    5   33  -16    1.39   67.56     -0.6920E+01  0.7164E+01
   21    5   28    1.39   40.27     -0.4125E+01  0.6715E+01
   18   18   25    1.39   73.30     -0.7508E+01  0.6289E+01
    9   33  -16    1.40   57.00     -0.5838E+01  0.4831E+01
    9   19   26    1.40   32.19     -0.3297E+01  0.5589E+01
   68    2   11    1.40   75.61     -0.7744E+01  0.6926E+01
   61   21    6    1.40   98.70     -0.1011E+02  0.9464E+01
   47   29    5    1.40   97.77     -0.1001E+02  0.8317E+01
   18   20   24    1.40   61.08     -0.6256E+01  0.8766E+01
   11   13   28    1.40   57.33     -0.5872E+01  0.8153E+01
   72    2    8    1.40   180.2     -0.1845E+02  0.6999E+01
   15   21   24    1.40   127.1     -0.1302E+02  0.6028E+01
   57   21    9    1.40   44.73     -0.4582E+01  0.8465E+01
   59   23    5    1.40   75.94     -0.7778E+01  0.6682E+01
   26    4   27    1.40   57.91     -0.5932E+01  0.6571E+01
    4   22   25    1.40   65.66     -0.6725E+01  0.7822E+01
   65    7   12    1.40   128.8     -0.1319E+02  0.1052E+02
   17   13   27    1.40   93.41     -0.9568E+01  0.9186E+01
   12   20   25    1.40   72.62     -0.7438E+01  0.6033E+01
   13   15   27    1.40   135.4     -0.1387E+02  0.5733E+01
   12    8   29    1.40   70.09     -0.7179E+01  0.6144E+01
   61   19    8    1.40   33.68     -0.3450E+01  0.6563E+01
   14    6   29    1.40   43.36     -0.4441E+01  0.6392E+01
   19    7   28    1.40   154.2     -0.1579E+02  0.6713E+01
   63    1   14    1.40   35.13     -0.3598E+01  0.7197E+01
    1   29   20    1.40   35.30     -0.3616E+01  0.5663E+01
   18    8   28    1.40   132.4     -0.1356E+02  0.6449E+01
   60   16   11    1.40   94.67     -0.9696E+01  0.6247E+01
   61   15   11    1.40   33.81     -0.3463E+01  0.6486E+01
    0   28   21    1.40   178.2     -0.1826E+02  0.7871E+01
   36   34   -2    1.40   82.75     -0.8475E+01  0.5157E+01
   45   31    1    1.40   81.45     -0.8343E+01  0.5597E+01
   52   28    2    1.40   126.6     -0.1296E+02  0.8754E+01
   51   29    0    1.40   118.6     -0.1215E+02  0.8394E+01
    4   20   26    1.40   101.8     -0.1043E+02  0.7399E+01
    1   31  -18    1.40   113.9     -0.1166E+02  0.8076E+01
   62   20    6    1.40   151.9     -0.1555E+02  0.6774E+01
   70    8    8    1.40   32.53     -0.3332E+01  0.7316E+01
   62   18    8    1.40   50.37     -0.5159E+01  0.6477E+01
    8   28   20    1.40   50.04     -0.5125E+01  0.8213E+01
   67   11    9    1.40   125.4     -0.1285E+02  0.6987E+01
   50   28    4    1.40   37.91     -0.3883E+01  0.5683E+01
   21    3   28    1.40   106.9     -0.1095E+02  0.6706E+01
    7   27   21    1.40   131.0     -0.1342E+02  0.5833E+01
   58   24    4    1.40   65.33     -0.6692E+01  0.6490E+01
   52   26    6    1.40   108.8     -0.1115E+02  0.5808E+01
   16   16   26    1.40   38.49     -0.3942E+01  0.6347E+01
   65   15    8    1.40   52.52     -0.5379E+01  0.7074E+01
   69    7    9    1.40   43.61     -0.4467E+01  0.7081E+01
   53   27    3    1.40   38.89     -0.3983E+01  0.6352E+01
   15    3   29    1.40   59.36     -0.6080E+01  0.6347E+01
   60   22    5    1.40   105.7     -0.1083E+02  0.9567E+01
   48   30    0    1.40   92.37     -0.9461E+01  0.8452E+01
   63   19    6    1.40   65.87     -0.6746E+01  0.6862E+01
   21    1   28    1.41   107.3     -0.1099E+02  0.6688E+01
    2   30  -19    1.41   156.6     -0.1604E+02  0.7543E+01
   66   12    9    1.41   46.29     -0.4742E+01  0.7072E+01
   13   11   28    1.41   90.45     -0.9264E+01  0.8508E+01
   65    9   11    1.41   48.15     -0.4931E+01  0.6876E+01
   10   18   26    1.41   55.67     -0.5702E+01  0.5588E+01
   64   18    6    1.41   78.74     -0.8065E+01  0.6680E+01
   14    4   29    1.41   135.3     -0.1385E+02  0.6291E+01
   47   29    4    1.41   61.97     -0.6347E+01  0.8336E+01
   20   10   27    1.41   150.6     -0.1543E+02  0.9530E+01
   16   18   25    1.41   146.7     -0.1502E+02  0.6199E+01
   65    5   12    1.41   120.8     -0.1238E+02  0.9526E+01
   56   22    8    1.41   118.1     -0.1209E+02  0.6091E+01
   69    9    8    1.41   63.29     -0.6482E+01  0.7161E+01
   34   34   -6    1.41   98.63     -0.1010E+02  0.5402E+01
   42   32   -1    1.41   103.8     -0.1063E+02  0.6738E+01
   19    5   28    1.41   44.42     -0.4550E+01  0.6814E+01
   13   19   25    1.41   158.4     -0.1623E+02  0.6058E+01
   60   18    9    1.41   51.19     -0.5243E+01  0.6418E+01
   50   26    7    1.41   38.15     -0.3907E+01  0.8177E+01
   11   11   28    1.42   34.29      0.2181E+02  0.6999E+01
   40   32   -3    1.43   30.16      0.1918E+02  0.7560E+01
    5   29   18    1.44   33.74      0.2146E+02  0.8467E+01
   23    7   26    1.45   30.09      0.1913E+02  0.6973E+01
   11    5   28    1.47   30.60      0.3169E+02  0.1006E+02
   62    6   10    1.51   39.02      0.1321E+02  0.7112E+01
   12   20   22    1.52   33.36      0.1129E+02  0.5876E+01
   46   24    7    1.52   54.58      0.1848E+02  0.9390E+01
   55    9   13    1.52   31.46      0.1065E+02  0.9319E+01
   40   28    4    1.52   31.62      0.1071E+02  0.7626E+01
    1   29  -16    1.52   34.90      0.1182E+02  0.7173E+01
   55   15   10    1.52   33.04      0.1119E+02  0.9514E+01
    8   30  -15    1.52   33.26      0.1126E+02  0.6352E+01
   45   25    6    1.52   30.67      0.1038E+02  0.5816E+01
   14   12   25    1.52   31.23      0.1057E+02  0.6068E+01
   49   25    2    1.53   32.99      0.1117E+02  0.8694E+01
   23    1   25    1.53   50.44      0.1708E+02  0.6800E+01
   54   14   11    1.53   35.24      0.1193E+02  0.9575E+01
   19   17   22    1.53   41.30      0.1398E+02  0.8534E+01
   15   11   25    1.53   40.09      0.1357E+02  0.6280E+01
   17   13   24    1.53   56.47      0.1912E+02  0.6233E+01
   10   18   23    1.54   66.13      0.2239E+02  0.5842E+01
    6   26   18    1.54   40.56      0.1373E+02  0.5366E+01
   61    5   10    1.54   35.03      0.1186E+02  0.6870E+01
   63    7    8    1.54   48.49      0.1642E+02  0.9621E+01
   14   14   24    1.54   82.10      0.2780E+02  0.8496E+01
   12    8   26    1.54   30.05      0.1017E+02  0.6101E+01
   37   29    2    1.55   40.47      0.1370E+02  0.7047E+01
   58    2   12    1.55   30.67      0.1038E+02  0.6854E+01
   58   12    9    1.55   32.88      0.1113E+02  0.8700E+01
   57    9   11    1.55   51.70      0.1750E+02  0.6599E+01
   51   17   10    1.55   31.62      0.1070E+02  0.6221E+01
   25    3   24    1.55   33.27      0.1126E+02  0.6444E+01
   21   15   22    1.55   32.59      0.1103E+02  0.8703E+01
   14   10   25    1.55   48.36      0.1637E+02  0.6320E+01
   25    9   23    1.55   55.83      0.1890E+02  0.9039E+01
   14    4   26    1.56   51.50      0.1744E+02  0.6391E+01
    0   28   16    1.56   30.70      0.1039E+02  0.5000E+01
   61    9    8    1.56   31.83      0.1078E+02  0.6932E+01
   24    4   24    1.56   40.48      0.1371E+02  0.6825E+01
   49   21    7    1.56   36.75      0.1244E+02  0.5930E+01
   15   15   23    1.56   65.00      0.2200E+02  0.8766E+01
   59   13    7    1.56   34.67      0.1174E+02  0.9281E+01
   15   19   21    1.56   31.19      0.1056E+02  0.6048E+01
   34   30   -2    1.56   36.70      0.1242E+02  0.6678E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT       772
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"       13561
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)   12208
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))     1353


 cpu time used                  0.3 sec
 elapsed wall-clock time        0.2 sec
