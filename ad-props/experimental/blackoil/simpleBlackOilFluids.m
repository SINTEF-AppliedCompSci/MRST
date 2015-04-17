function fluid = simpleBlackOilFluids(mycomp,mycomp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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

