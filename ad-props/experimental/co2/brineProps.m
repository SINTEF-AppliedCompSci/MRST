function brine = brineProps(salinity)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

if(nargin<1)
    salinity= 0.03;
end
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
h2o=h2oProps();

brine.viscosity = @(p,T) viscosity(T,p,salinity);
brine.density = @(p,T) density(T,p,h2o,salinity);
brine.enthalpy = @(p,T) enthalpy(T,p,h2o,salinity);
brine.name='brine';
%brine=viscosity(T,p,salinity);
%brine=density(T,p,h2o_s,salinity);
%brine=enthalpy(T,p,h2o_s,salinity);
end

function v =density(T, p,h2o, salinity)
TempC = T - 273.15;
pMPa = p/1.0E6;

rhow = h2o.density(p, T);
v=rhow +...
    1000*salinity.*(...
    0.668 +...
    0.44*salinity +...
    1.0E-6.*(...
    300.*pMPa -...
    2400.*pMPa.*salinity +...
    TempC.*(...
    80.0 +...
    3.*TempC -...
    3300*salinity -...
    13.*pMPa +...
    47.*pMPa*salinity)));
end
function v = enthalpy(T,p,h2o, salinity)

%/*Numerical coefficents from PALLISER*/
f = [2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10];


%/*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
a = [
    -9633.6, -4080.0, +286.49 ;
    +166.58, +68.577, -4.6856 ;
    -0.90963, -0.36524, +0.249667E-1 ;
    +0.17965E-2, +0.71924E-3, -0.4900E-4
    ];

%Scalar theta, h_NaCl;
%Scalar m, h_ls, h_ls1, d_h;
%Scalar S_lSAT, delta_h;
%int i, j;
%Scalar hw;

theta = T - 273.15;

S = salinity;
%S_lSAT = f(1) + f(2)*theta + f(3)*pow(theta,2) + f(4)*pow(theta,3);
S_lSAT=0;
for i=1:4
    S_lSAT=S_lSAT+f(i)*power(theta,i-1);
end
%/*Regularization*/
if (S>S_lSAT)
    S = S_lSAT;
end

hw = h2o.enthalpy(p, T)/1E3;% /* kJ/kg */

%/*DAUBERT and DANNER*/
%/*U=*/
h_NaCl = (3.6710E4*T + 0.5*(6.2770E1).*T.*T - ((6.6670E-2)/3).*T.*T.*T...
                +((2.8000E-5)/4)*power(T,4))/(58.44E3)- 2.045698e+02;% /* kJ/kg */

m = (1E3/58.44)*(S./(1-S));
%i = 0;
%j = 0;
d_h = 0;

for i=1:4
    for j=1:3
        d_h = d_h + a(i,j) * power(theta, i-1) * power(m, j-1);
    end
end

delta_h = (4.184./(1E3 + (58.44 * m))).*d_h;

%/* Enthalpy of brine */

h_ls1 =(1-S).*hw + S.*h_NaCl + S.*delta_h; %/* kJ/kg */

h_ls = h_ls1.*1E3; %/*J/kg*/

v = h_ls;
end
function v= viscosity(T, p, salinity)
    
        if(T <= 275.) %// regularisation
           T = 275;
        end
        T_C = T - 273.15;

        A = (0.42*power((power(salinity, 0.8)-0.17), 2) + 0.045).*power(T_C, 0.8);
        mu_brine = 0.1 + 0.333*salinity + (1.65+91.9*salinity*salinity*salinity).*exp(-A);

        v=mu_brine/1000.0;% /* unit: Pa s */
end
