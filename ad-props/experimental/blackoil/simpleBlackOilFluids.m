function fluid = simpleBlackOilFluids(mycomp,mycomp)
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

fluid=initSimpleADIFluid('mu',[0.5 5 1e-1]*centi*poise,'rho',[1000 1000 1000],'n',[2 1 1]);

% fluid comp
switch mycomp
    case 'comp'
        cR= 0.00045/barsa;pRef=235*barsa;
        fluid.pvMultR = @(p)(1+cR.*(p-pRef));
        cRw= 5.43E-05/barsa;pRefw=234.46*barsa;rfw= 1.0042;
        fluid.bW   =@(p) (1+cRw.*(p-pRefw))/rfw;
        cRo= 6.65E-05/barsa;pRefo=234*barsa;rfo=1.065;
        fluid.bO   =@(p) (1+cRo.*(p-pRefo))/rfo;
    case 'incomp'
        cR= 0.000/barsa;pRef=235*barsa;
        fluid.pvMultR = @(p)(1+cR.*(p-pRef));
        cRw=  0/barsa;pRefw=234.46*barsa;rfw= 1;
        fluid.bW   =@(p) (1+cRw.*(p-pRefw))/rfw;
        cRo=  0/barsa;pRefo=234*barsa;rfo=1;
        fluid.bO   =@(p) (1+cRo.*(p-pRefo))/rfo;
        
    otherwise
        error();
end

%% presure dependent viscosity

%% energy
cW       = 4.1813*10.^6/1e3; % energy density per mass
cR       = 2.17*10.^6;% energy density per volume

%my_fluid='linear'
%my_fluid='simple';
switch my_fluid
    case 'linear'        
        fluid.uW = @(p,T) cW.*T;
        fluid.uG = @(p,T) 0.1*cW.*T;
        fluid.uO = @(p,T) 0.7*cW.*T;
        fluid.uR = @(T) cR.*T;

        fluid.hW = @(p,T) cW.*T;
        fluid.hG = @(p,T) 0.1*cW.*T;
        fluid.hO = @(p,T) 0.7*cW.*T;
        %fluid.hO = @(p,T) 1*cW.*T;
    case 'simple'
        fluid.uW = @(p,T) cW.*T;
        fluid.uG = @(p,T) 0.1*cW.*T;
        fluid.uO = @(p,T) 0.7*cW.*T;
        fluid.uR = @(T) cR.*T;
        fluid.hW = @(p,T) cW.*T + (1-T.*207e-6).*p/(cW);
        fluid.hW = @(p,T) cW.*T + 2.2e-7.*p.*(cW);
        cO=2.2e6/1e3;
        fluid.hO = @(p,T) cO.*T+(1-T.*950e-6).*p/(cO);% gasonline        
    otherwise
        error()
end

end
