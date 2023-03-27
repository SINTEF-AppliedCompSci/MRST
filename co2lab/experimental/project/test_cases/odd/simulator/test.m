


Pt = 100 * barsa;

CO2_props = CO2props;

CO2.beta = @CO2_props.beta;
CO2.bder = @CO2_props.betader;
CO2.rho  = @CO2_props.rho;

Wat.beta = @(p, t) 4.2e-10;
Wat.bder = @(p, t) 0;
Wat.rho = @(p, t) 1000 * (1 + Wat.beta(p,t) * p);

h = 40;
H = 100;
T = 310;
theta = 0;

[Pb, Pi] = Pt2PbPi(Pt, CO2, Wat, h, H, T, theta);

[Pt2, Pi2] = Pb2PtPi(Pb, CO2, Wat, h, H, T, theta);