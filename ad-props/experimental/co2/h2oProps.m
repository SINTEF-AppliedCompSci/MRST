function h2o = h2oProps()
% from dumux h2o

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

    h2o.isvalid =@(p,T) T <= 623.15 & p<= 100e6;
    h2o.density =@(p,T) R1_density(p,T);
    h2o.enthalpy =@(p,T) R1_entalpy(p,T);
    h2o.viscosity = @(p,T) R1_viscosity(p,T);
    h2o.name='h2o';
end
function v= R1_viscosity(p,T)
    rho=R1_density(p,T);
    v =common_viscosity(T,rho);
end
function v = common_viscosity(T,rho)
rhoc=322;
rhoBar = rho/rhoc;
cT=647.096;
TBar = T/cT;
Hij = Visc_const();
tmp3 = 1;
muBar = 0;
for i = 1:6
    tmp = 0;
    tmp2 = 1;
    for j = 1:7
        tmp = tmp + Hij(i,j).*tmp2;
        tmp2 =tmp2.* (rhoBar - 1);
    end
    muBar =muBar+ tmp3.* tmp;
    tmp3 =tmp3.*(1.0./TBar - 1);
end
muBar =muBar.*rhoBar;
muBar = exp(muBar);

muBar  =muBar.*100.*power(TBar,0.5);
H = [1.67752, 2.20462, 0.6366564, -0.241605];
tmp = 0; tmp2 = 1;
for i = 1:4
    tmp =tmp + H(i)./tmp2;
    tmp2 =tmp2.* TBar;
end
muBar =muBar./ tmp;
v= 1e-6*muBar;
end
function d=R1_density(p,T)
    d=1./volumeRegion1(T,p);
end
function v = R1_entalpy(p,T)
    molarMass = 18.01518e-3;
    R=8.314472;
    Rs=R/molarMass;
    v= R1_tau(T).*R1_dGamma_dTau(T,p).*Rs.*T;
end

function v = volumeRegion1(T,p)
molarMass = 18.01518e-3;
R=8.314472;
Rs=R/molarMass;
    v= R1_pi(p).*R1_dGamma_dPi(T,p).*Rs.*T./p;
end
function v = R1_pi(p)
     v=p ./ 16.53e6;
end
function v = R1_tau(T)
     v=1386.0 ./ T;
end
function result = R1_dGamma_dTau(T,p)
     tau_ =R1_tau(T);
     pi_ = R1_pi(p);
     [n,I,J] = R1_constants();
     result = 0.0;
     for i = 1:34
         result = result + ...
             n(i) .*...
             power(7.1 - pi_, I(i)) .*...
             power(tau_ - 1.222,  J(i)-1) .*...
                J(i);
      end
     
end

function result = R1_dGamma_dPi(T,p)
     tau_ =R1_tau(T);
     pi_ = R1_pi(p);
     [n,I,J] = R1_constants();
     result=0;
     for  i = 1:34 
            result = result-...
                n(i).*...
                I(i).*...
                power(7.1 - pi_, I(i) - 1).*...
                power(tau_ - 1.222, J(i));
     end    
end

function [n,I,J] = R1_constants();

n=[         0.14632971213167, -0.84548187169114, -0.37563603672040e1,...
    0.33855169168385e1, -0.95791963387872, 0.15772038513228,...
    -0.16616417199501e-1, 0.81214629983568e-3, 0.28319080123804e-3,...
    -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1,...
    -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3,...
    -0.30001780793026e-3, 0.47661393906987e-4, -0.44141845330846e-5,...
    -0.72694996297594e-15,-0.31679644845054e-4, -0.28270797985312e-5,...
    -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6,...
    -0.14341729937924e-12,-0.40516996860117e-6, -0.12734301741641e-8,...
    -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19,...
    0.26335781662795e-22,-0.11947622640071e-22, 0.18228094581404e-23,...
    -0.93537087292458e-25...
    ];


I = [
    0, 0, 0,...
    0, 0, 0,...
    0, 0, 1,...
    1, 1, 1,...
    1, 1, 2,...
    2, 2, 2,...
    2, 3, 3,...
    3, 4, 4,...
    4, 5, 8,...
    8, 21, 23,...
    29, 30, 31,...
    32 ...
    ];

J =[
    -2, -1, 0,...
    1, 2, 3,...
    4, 5, -9,...
    -7, -1, 0,...
    1, 3, -3,...
    0, 1, 3,...
    17, -4, 0,...
    6, -5, -2,...
    10, -8, -11,...
    -6, -29, -31,...
    -38, -39, -40,...
    -41....
    ];
end
function H = Visc_const()
 H = [5.20094e-1, 2.22531e-1,-2.81378e-1, 1.61913e-1,-3.25372e-2, 0, 0;...
             8.50895e-2, 9.99115e-1,-9.06851e-1, 2.57399e-1, 0, 0, 0 ;
             -1.08374 , 1.88797 ,-7.72479e-1, 0, 0, 0, 0;...
            -2.89555e-1, 1.26613 ,-4.89837e-1, 0, 6.98452e-2, 0,-4.35673e-3;...
            0, 0,-2.57040e-1, 0, 0, 8.72102e-3, 0;...
            0, 1.20573e-1, 0, 0, 0, 0,-5.93264e-4];
end

