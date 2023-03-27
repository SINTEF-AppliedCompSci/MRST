% ============ SCRIPT SPECIFYING SYSTEM UNITS AND CHARACTERISTICS ============

% ----------------------------------- Units -----------------------------------
km     = 1000; % base unit is _meter_
m      = 1;    % base unit is _meter_
cm     = 0.01; % base unit is _meter_
Pa     = 1;    % base unit is _Pascal_
MPa    = 1e6;  % base unit is _Pascal_
kg     = 1;    % base unit is _kilogram_
sec    = 1;    % base unit is _second_
centi  = 1e-2;
milli  = 1e-3;
poise  = 0.1 * Pa * sec; 
atm    = 101325 * Pa;
Darcy  = (cm^3)/sec/(cm^2) * (centi*poise) / (atm/cm);
minute = 60 * sec;
hour   = 60 * minute;
day    = 24 * hour;
year   = 365 * day;

% ------------------------ Defining physics and domain ------------------------

L       = 393216 * m    ;            % Physical dimension of (square) simulation domain
l       = 3;                         % grid refinement level (ref. level 9 should take about 2 hours....)
[dx_vec, N, h] = define_grid(l,L);

Bpress  = 31.36 * MPa;               % (fixed) boundary presure
Sres    = 0.29;                      % Residual brine saturation in CO2 phase
H       = 38;                        % Thickness of aquifer (meter)
k       = 23 * milli * Darcy;        % Taken from the Basal sandstone in Janzen Thesis

%k = k * 100;
phi     = 0.12;                      % Taken from the Basal sandstone in Janzen Thesis              
rhoc    = 768 * kg/m^3;              % Taken from the Basal sandstone in Janzen Thesis              
rhob    = 1161 * kg/m^3;             % Taken from the Basal sandstone in Janzen Thesis
krc     = 0.54;                      % Rel. perm. of CO2 at saturation (1 - Sres)
muc     = 0.067*milli*Pa/sec;        % CO2 viscosity
mub     = 0.624*milli*Pa/sec;        % brine viscosity
lambdac = krc/muc;                   % mobility of CO2 at saturation (1 - Sres)
lambdab = 1/mub;                     % brine mobility at full saturation
Q       = 60 * kg / sec;             % Injection rate (well centered at middle of domain)
Qvol    = Q / rhoc;                  % volumetric injection rate (when only
                                     % injecting CO2) 
c_c     = 7.25e-9;                   % CO2 compressibility (@@ unused)
c_b     = 4.6e-10;                   % brine compressibility (@@ unused)
c_r     = 1e-10;                     % aquifer rock compressibility