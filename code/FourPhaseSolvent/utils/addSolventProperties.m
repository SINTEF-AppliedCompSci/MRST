function fluid = addSolventProperties(fluid, varargin)
% Add solvent phase and properties to fluid

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

opt = struct('mu'              , 1   , ...
             'rho'             , 1   , ...
             'n'               , 1   , ...
             'b'               , 1   , ...
             'c'               , []  , ...
             'pRef'            , 0   , ...
             'mixPar'          , 1   , ...
             'sOres_m'         , 0   , ...
             'sOres_i'         , 0   , ...
             'sSGres_m'        , 0   , ...
             'sSGres_i'        , 0   , ...
             'sWres'           , 0   , ...
             'Msat'            , []  , ...
             'Mpres'           , []  , ...
             'useTabulatedMsat', true);

opt = merge_options(opt, varargin{:});

%%

% Standard properties (b, kr, rho, mu)

b = opt.b;
if isempty(opt.c)
    % Constant value (incompressible phase)
    bf = @(p, varargin) b*constantReciprocalFVF(p, varargin{:});
else
    % Compressibility on the form
    % b = b_ref exp((p-p_ref)*c)
    c = opt.c;
    if c < 0
        warning('Negative compressibility detected.')
    end
    bf = @(p, varargin) b*exp((p-opt.pRef)*c);
end
    
fluid.rhoSS = opt.rho;
fluid.bS = bf;
fluid.muS = @(p, varargin) constantViscosity(opt.mu, p, varargin{:});

% Extra model properties

fluid.mixPar   = opt.mixPar;
fluid.sOres_m  = opt.sOres_m;
fluid.sOres_i  = opt.sOres_i;
fluid.sSGres_m = opt.sSGres_m;
fluid.sSGres_i = opt.sSGres_i;
fluid.sWres    = opt.sWres;


nW = opt.n(1); nO = opt.n(2); nG = opt.n(3);
fluid.krW_i  = @(sW, sWres, sOres, sSGres) coreyRelperm(sW, nW, sWres, fluid.krW(1-(sOres + sSGres)) , sWres + sOres + sSGres);
fluid.krO_i  = @(sO, sWres, sOres, sSGres) coreyRelperm(sO, nO, sOres , fluid.krO(1-(sWres + sSGres)), sWres + sOres + sSGres);
fluid.krGT_i = @(sGT, sWres, sOres, sSGres) coreyRelperm(sGT, nG, sSGres, fluid.krG(1-(sWres+sOres)), sWres + sOres + sSGres);
fluid.krN    = @(sN, sWres, sOres, sSGres) coreyRelperm(sN, nO, sOres + sSGres, fluid.krO(1-sWres) , sWres + sOres + sSGres);


fluid.krW = @(sW, varargin) fluid.krW_i (sW, fluid.sWres, fluid.sOres_i, fluid.sSGres_i);
fluid.krO = @(sO, varargin) fluid.krO_i (sO, fluid.sWres, fluid.sOres_i, fluid.sSGres_i);
fluid.krG = @(sG, varargin) fluid.krGT_i(sG, fluid.sWres, fluid.sOres_i, fluid.sSGres_i);
fluid.krOW = fluid.krO;
fluid.krOG = fluid.krO;

fluid.Msat = opt.Msat;
if isempty(opt.Msat)
    if opt.useTabulatedMsat
        T = tabulatedSaturationFraction();
        
        if 0
            a = 1e-3;
            b = 1e-10;
            % b = 0;
            smthr = @(s) (s > a+b) + (s >= b & s <= a+b).*((s-b)/a).^2./(((s-b)/a).^2 + (1-(s-b)/a).^2);
        elseif 1
            smin = 1e-12;
            smax = 1e-3;
            smthr = @(s)  min(max((s - smin)/(smax - smin), 0), 1);
        else
            tol = 1e-12;
            smthr = @(sS) sS>tol;
        end
        tol = 1e-12;
        fluid.Msat = @(sG, sS) interp2DTable(T, sG, sS).*smthr(sS);
%         fluid.Msat = @(sG, sS) interp2DTable(T, sG, sS).*(sS>tol);
        fluid.satFrac = @(sA, sB) interp2DTable(T, sB, sA).*(sA>tol);
    else
        fluid.Msat = @(sG, sS) linearSaturationMiscibility(sG, sS);
    end
end

fluid.Mpres = opt.Mpres;
if isempty(opt.Mpres)
    fluid.Mpres = @(p) constantPressureMiscibility(p);
end



end

function [nW, nO, nG] = getRelpermExponent(fluid)

    tol = 0.2;
    sWres = fluid.sWres;
    sOres = fluid.sOres_i;
    sSGres = fluid.sSGres_i;

    n = -100:0.1:0;
    s = 10.^n;
    
    sW = s(s>sWres & s < 1 - (sOres + sSGres));
    nW = polyfit(log(sW), log(fluid.krW(sW)), 1); nW = nW(1);    

    sO = s(s>sOres & s < 1 - (sWres + sSGres));
    nO = polyfit(log(sO), log(fluid.krO(sO)), 1); nO = nO(1);    

    sG = s(s>sSGres & s < 1 - (sWres + sOres));
    nG = polyfit(log(sG), log(fluid.krG(sG)), 1); nG = nG(1);    
    
end

function B = constantReciprocalFVF(p, varargin)
    B = p*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
    mu = p*0 + mu;
end

function kr = coreyRelperm(s, n, sr, kwm, sr_tot)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end
    
    kr = kwm.*sat.^n;
end

function T = tabulatedSaturationFraction()
    
    n  = 400;
    ds = 1e-4;
    s  = linspace(-ds,1+ds,n)';
    [ss, sg] = meshgrid(s,s);
    M = ss./(ss + sg);
    
    assert(~any(any(isnan(M))));
    M(M<0) = 0; M(M>1) = 1;
    
    T = struct('x', s, 'y', s, 'data', M);
    
end

function Msat = linearSaturationMiscibility(sG, sS)

    tol = 0;
    sS(sS < tol) = 0;
    sG(sG < tol) = 0;
    
%     tol = 1e-2;
%     Msat = sS./(sS + sG + tol).*(1+tol);
    
    Msat = sS./(sS + sG);
    Msat(isnan(double(Msat))) = 0.5;

end

function Mpres = constantPressureMiscibility(p)

    Mpres = p.*0 + 1;
    
end
    