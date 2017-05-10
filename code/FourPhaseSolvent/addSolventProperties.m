function fluid = addSolventProperties(fluid, varargin)

opt = struct('mu'    , 1                      , ...
             'rho'   , 1                      , ...
             'n'     , 1                      , ...
             'b'     , 1                      , ...
             'c'     , []                      , ...
             'pRef'  , 0                      , ...
             'mixPar', 1                      , ...
             'sOres_m', 0                     , ...
             'sOres_i', 0                     , ...
             'sSGres_m', 0                     , ...
             'sSGres_i', 0                     , ...
             'sWres',    0                    , ...
             'Msat'  , @(sG, sS) sS./(sS + sG), ...
             'Mpres' , @(p) ones(numel(p),1)   );
             

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
fluid.Msat     = opt.Msat;
fluid.Mpres    = opt.Mpres;

end

function B = constantReciprocalFVF(p, varargin)
    B = p*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
    mu = p*0 + mu;
end

