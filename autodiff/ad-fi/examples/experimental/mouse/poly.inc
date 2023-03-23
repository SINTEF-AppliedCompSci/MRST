NOECHO
-- Polymer properties (keywords starting with PLY)

-- PLYVISC: Viscosity multiplier. Defined as PLYVISC=u/uw, where uw is water viscosity. uw@50 C is 0.5 mPa*s (cP)
--          Tables defined twice as NTSFUN=2
PLYVISC
0    1.0
0.5  3.0
1.0  5.0
1.5  15.0
2.0  30.0
3.0  50.0 /


-- PLYROCK gives fraction of inaccessible pore volume (IPV, must be less than Swi)
-- residual resistivity factor (RRF), density of rock incl. PV (P-sandgrain=2.7 g/cm3 + poro=0.185)
-- desorpsjon av polymer (paa=1,av=2), og RRF. ved maksimalt mengde adsorbert polymer
-- IPV: Eclipse krever at IPV<Swcr. Vi har derfor satt denne til 0.199. Our recommended base-case value is 0.25 
-- RRF: 33% reduction in water. perm (33% was observed for Oerrell)
-- AI: adsorbtion index: 1=desorption may occur, 2= NO desorption
-- Max. adsorbert polymer (kg/kg), 

PLYROCK
--IPV  RRF  dens  AI  max ads
 0.05 1.3  2650   2 0.000020 /


--PLYADS: Adsorption isotherm
PLYADS
-- conc  adsorb-conc
 0.000   0.000000       
 0.250	 0.00001
 0.500	 0.00002
 0.750	 0.000015
 1.000	 0.000017
 1.250	 0.000018
 1.500	 0.000019
 1.750	 0.000020
 2.000	 0.000020
 3.000	 0.000020 /
  
--Polymer: Todd-Longstaff mixing parameter. 1 = full mixing
TLMIXPAR
 1.0  /

--Max. concentration to be used in the mixing calculation (set equal to pol. inj. conc.)
PLYMAX
 3.00  0  /



ECHO
