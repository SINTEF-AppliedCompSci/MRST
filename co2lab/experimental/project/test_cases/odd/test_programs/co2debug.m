CO2 = CO2props;

gravity on;
[pc, tc] = CO2CriticalPoint;
pRefMin = 5.0e5;%0.8e7; %1.0e5;
pRefMax = pc - (pc - pRefMin) * 0.05;%7e7;%1.0e8;%1.0e8;
tRefMin = 280;
tRefMax = tc - (tc - tRefMin) * 0.05;;%400;

bounds = [tRefMin tRefMax pRefMin pRefMax];

p_res = 200;
t_res = 200;

tvals = linspace(tRefMin, tRefMax, t_res);
pvals = linspace(pRefMin, pRefMax, p_res);

[T_top, P_top] = meshgrid(tvals, pvals);

top_rho  = reshape(CO2.rho    (P_top(:), T_top(:)), p_res, t_res);
top_der  = reshape(CO2.rhoder (P_top(:), T_top(:)), p_res, t_res);
top_beta = reshape(CO2.beta   (P_top(:), T_top(:)), p_res, t_res);
top_bder = reshape(CO2.betader(P_top(:), T_top(:)), p_res, t_res);

