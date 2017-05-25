function fluid = addSolventProperties(fluid, varargin)

opt = struct('mu'    , 1                      , ...
             'rho'   , 1                      , ...
             'n'     , 1                      , ...
             'b'     , 1                      , ...
             'c'     , []                     , ...
             'pRef'  , 0                      , ...
             'mixPar', 1                      , ...
             'sOres_m', 0                     , ...
             'sOres_i', 0                     , ...
             'sSGres_m', 0                    , ...
             'sSGres_i', 0                    , ...
             'sWres',    0                    , ...
             'Msat'  , [], ...
             'Mpres' , []);

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
kr = @(s) s.^opt.n;
    
fluid.rhoSS = opt.rho;
fluid.bS = bf;
fluid.muS = @(p, varargin) constantViscosity(opt.mu, p, varargin{:});
fluid.krS = kr;

% Extre model properties

fluid.mixPar   = opt.mixPar;
fluid.sOres_m  = opt.sOres_m;
fluid.sOres_i  = opt.sOres_i;
fluid.sSGres_m = opt.sSGres_m;
fluid.sSGres_i = opt.sSGres_i;
fluid.sWres    = opt.sWres;

fluid.Msat = opt.Msat;
if isempty(opt.Msat)
    fluid.Msat = @(sG, sS) linearSaturationMiscibility(sG, sS);
end

fluid.Mpres = opt.Mpres;
if isempty(opt.Mpres)
    fluid.Mpres = @(p) constantPressureMiscibility(p);
end



end

function B = constantReciprocalFVF(p, varargin)
    B = p*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
    mu = p*0 + mu;
end

function Msat = linearSaturationMiscibility(sG, sS)
% Miscibility is dependent on fraction of solvent in the total gas
% saturation. We add small value tol to avoid 0/0, and specifically assign
% Msat = 0 if both sG and sS is very small.

%     tol = 100*eps;
%     zz = sG < tol & sS < tol;
%     Msat = (sS./(sG + sS + (sS < tol).*tol)).*(~zz) + 0.*zz;

%     tol = eps;
%     sS  = sS - sS.*(sS < tol);
%     sG  = sG - sG.*(sG < tol);
    
    
%     sS = trimSaturations(sS);
%     sG = trimSaturations(sG);

%     sS = max(sS,0);

%     tol = 1e-2;
%     sG(sG < tol) = tol;
%     sS(sS < tol) = tol;
%     sS(sS < 10*eps) = 0;

%     Msat = saturationFraction(sS, sS + sG);

%      Msat = saturationFraction(sS, sG);
     Msat = 1 - saturationFraction(sG, sS);

%     tol = 10*eps;
%     zz = sG < tol & sS < tol;
%     sS  = sS - sS.*(sS < tol);
%     sG  = sG - sG.*(sG < tol);
%     
%     Msat = (sS./(sG + sS + (sS < tol).*tol)).*(~zz) + 0.*zz;
    
end

function Mpres = constantPressureMiscibility(p)

    Mpres = p.*0 + 1;
    
end