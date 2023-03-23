NOECHO
-- Surfactant properties
-- Values are taken from SPE paper 145036,
-- Title   : Case Study: Numerical Simulation of Surfactant Flooding in Low Permeability Oil Field,
-- Authors : Xu, Guo, Wang, Zhang, Jia, Wang.

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
