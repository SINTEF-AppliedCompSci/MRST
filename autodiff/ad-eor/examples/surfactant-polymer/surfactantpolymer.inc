NOECHO
PLYVISC
0      1.0
0.5    4.0
1.0    8.0
1.5   13.0
2.0   26.0
3.0   52.0 /


PLYROCK
--IPV  RRF  dens  AI  max ads
 0.05  1.3  2600  2  0.000025 /
 0.05  1.3  2600  2  0.000025 /


PLYADS
-- conc  adsorb-conc
 0.000   0.000000
 0.250   0.000014
 0.500   0.000018
 0.750   0.000021
 1.000   0.000022
 1.250   0.000023
 1.500   0.000023
 1.750   0.000027
 3.000   0.000027 /
  0.000   0.000000
 0.250   0.000014
 0.500   0.000018
 0.750   0.000021
 1.000   0.000022
 1.250   0.000023
 1.500   0.000023
 1.750   0.000027
 3.000   0.000027 /

PLMIXPAR
-- 1.0  /
 0.2 /

PLYMAX
 3.00  0  /

PLYSHEAR
-- Vw.    VRF
0           1
0.00002   0.9
0.00008   0.8
0.0004    0.74
0.0008    0.69
0.0017    0.67
0.0043    0.63
0.011     0.60
0.015     0.58
0.019     0.58
0.22      0.48
2.2       0.41   /


SURFST
-- Concentration   Surface Tension
0                  0.05
0.1                0.0005
0.5                1e-5
1                  1e-6
30                 1e-6
100                1e-6 /

SURFCAPD
-- Log10(capillary number)   Capillary number interpolant
-10                          0
-5.5                         0
-4                           0.5
-3                           1
2                            1   /
-10                          0
-5.5                         0
-4                           0.5
-3                           1
2                            1   /

SURFVISC
-- Concentration   Water viscisity
0                  0.61
30                 0.8
100                1 /

SURFADS
-- Concentration   Adsorption coefficient
0                  0.
1                  5e-4
30                 5e-4
100                5e-4 /
0                  0
1                  5e-4
30                 5e-4
100                5e-4 /

SURFROCK
-- AI   Density
 1      2650 /
 1      2650 /

ECHO
