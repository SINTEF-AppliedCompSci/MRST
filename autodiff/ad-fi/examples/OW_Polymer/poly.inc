NOECHO
-- Polymer properties (keywords starting with PLY)

-- PLYVISC: Viscosity multiplier. Defined as PLYVISC=u/uw, where uw is water viscosity.
PLYVISC
0    1.0
0.5  3.0
1.0  5.0
1.5  15.0
2.0  30.0
3.0  60.0 
4.0  65.0 /

PLYROCK
--IPV  RRF  dens  AI  max ads
 0.05 1.1  2000   2 0.000050 /


--PLYADS: Adsorption isotherm
PLYADS
-- conc  adsorb-conc
 0.000   0.000000       
 0.250	 0.00001
 0.500	 0.000013
 0.750	 0.000015
 1.000	 0.000017
 1.250	 0.000018
 1.500	 0.000019
 2.000	 0.000023
 2.750	 0.000025
 3.000	 0.000025
 4.000	 0.000025 /
--Polymer: Todd-Longstaff mixing parameter. 1 = full mixing
TLMIXPAR
 1.0  /

--Max. concentration to be used in the mixing calculation
PLYMAX
 4.00  0  /


ECHO
