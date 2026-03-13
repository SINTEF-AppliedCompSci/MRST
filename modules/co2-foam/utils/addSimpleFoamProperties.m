function fluid = addSimpleFoamProperties(fluid, varargin)
%Add foam and extra properties for the CO2 foam model to fluid
%
% SYNOPSIS:
%       fluid = addSimpleFoamProperties(fluid)
%       fliud = addSimpleFoamProperties(fluid, 'pn1', 'pv1', ...)
%
% PARAMETERS:
%   fluid  - Model fluid where the parameters for foam effect shall be added.
%
% OPTIONAL ARGUMENTS:
%   adsMax          - Maximum adsorbed concentration, kg surfactant per kg rock.
%                     Defaults to 2.0e-4
%   adsSat          - Approximate surfactant concentration (mass fraction) where
%                     adsorption will start to saturate. Defaults to 1e-3.
%   rhoRSft         - Rock density, kg rock per bulk m3. Defaults to 2000.
%   foam            - Struct with foam parameters:
%                     sF    - Critical water saturation
%                     s1    - Slope in capillary regime
%                     s2    - Slope in pressure gradient regime
%                     F0    - Scale factor for large saturations
%                     v0    - Reference gas velocity
%   cCrit           - Minimum surfactant concentration for full foam strength,
%                     kg surfactant per kg fluid (solution). Defaults to 1.0e-3
%   cDecl           - Foam strength essentially zero for concentration below
%                     cCrit/cDecl. Defaults to 100.
%   Cpart           - Partitioning constant for surfactant concentration in
%                     gas and water. No default. Undefined if not specified.
%   surfingas       - True if surfactant is carried only in the gas phase. Default false.
%   useSmoother     - Smoothen mobility reduction function at critical water
%                     saturation. Default false.
%   mobMultIsLinear - True if mobility reduction factor does not depend on gas velocity.
%                     Default true.
%   noShear         - True if no shear dependence is included. Default false.
%   usePermDep      - True if permeability effect depends on absolute permeabiliyt. Default false.
%   mobMult         - Handle for mobility reduction factor function
%
% RETURNS:
%   fluid   - Fluid model with additional parameters added.

%{ 
Copyright 2009-2023 SINTEF

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

    opt = struct('adsMax'         , 2.0e-4 , ...
                 'adsSat'         , 1e-3   , ...
                 'rhoRSft'        , 2000   , ...
                 'foam'     , []     , ...
                 'fmax'  , []     , ...
                 'csmax'          , []     , ...    
                 'cCrit'          , 1.0e-3 , ...
                 'cDecl'          , 100    , ...
                 'surfingas'      , false  , ...
                 'useSmoother'    , false  , ...
                 'mobMultIsLinear', true   , ...
                 'usePermDep'     , false  , ...
                 'noShear'        , false  , ...
                 'mobMult'        , []     ); 

    [opt, extra] = merge_options(opt, varargin{:});

    %
    fluid.adsMax = opt.adsMax;
    fluid.adsSat = opt.adsSat;
    fluid.surfads = @(c) opt.adsMax*(2./(1 + exp(-c/opt.adsSat)) - 1.0); % Will saturate at about 7*opt.cads
    fluid.rhoRSft = opt.rhoRSft;
    fluid.csmax = opt.csmax;
    fluid.cCrit = opt.cCrit;
    fluid.cDecl = opt.cDecl;
    fluid.usePermDep = opt.usePermDep;

    % Check if partitioning is specified in input parameters. Leave undefined 
    % if not specified.
    if any(strcmpi('cpart',extra))
        idx = find(strcmpi('cpart',extra(1:2:end)));
        tmp = extra{2*idx(1)};
        if ~isempty(tmp)
            fluid.Cpart = tmp;
        end
    end
    
    fluid.surfingas = opt.surfingas;

    % Get foam model parameters
    if isempty(opt.foam)
        opt.foam = getFoam(1);
    end
    fluid.foam = opt.foam;
    
    % mobMult can be user-defined. If linear, we assume it returns two
    % values [A,B] so that F = A + B*|vG|.
    fluid.mobMultIsLinear = opt.mobMultIsLinear;
    if isempty(opt.mobMult)
        if opt.mobMultIsLinear
            % If mobility multiplicator F is on the form F = Fcap +
            % Fgrad*|vG|, we can solve for |vG| (absolute gas velocity)
            % explicitly to obtain F
            fluid.mobMult = @(sW, c) mobMultFunc(sW, c, [], fluid, opt);
        else
            % If not, we must solve a nonlinear problem. If no nonlinear
            % mobMult function is given, we simply return the linear
            % function
            fluid.mobMult = @(sW, c, vG) mobMultFunc(sW, c, vG, fluid, opt);
        end
    end
    fluid.noShear = opt.noShear;
end

function varargout = mobMultFunc(sW, c, vG, fluid, opt)
% Calculate gas mobility reduction factor based on water saturation,
% surfactant concentration and foam property parameters.

    foam = fluid.foam;

    F0 = foam.F0;
    s1 = foam.s1;
    s2 = foam.s2;
    sF = foam.sF;
    cfrac = foam.cfrac;
    v0 = foam.v0;
    
    cf = fluid.cCrit; % Concentration threshold
    cd = fluid.cDecl; % Foam strength decline interval
    
    if sF > 0 && opt.useSmoother
        % Hack to elimiate jump in F at sW = sF
        alpha = 5;
        C0 = (s1/s2)*1e-2;
        beta = C0/(sF*s1);
        smoother = sigmoid(sW, sF, alpha, beta);
%         smoother = ((sW - sF)./(0.2*sF)).^3.*(sW < 1.2*sF) + (sW >= 1.2*sF);
    else 
        smoother = 1;
    end
 
    % Mobility reduction factor is given by F = Fcap + Fgrad*|vG|
    % Capillary pressure part
    Fcap  = exp(-(sW-sF).*s1);
    % Pressure gradient part (to be multiplied by gas velocity)
    Fgrad = exp(-(sW-sF).*s2).*(F0./v0).*smoother;

    % Calculate low concentration correction factor
    if cf > 0
        % Increase from zero to one in the log-interval cf/cd to cf.
        G = 1./(1+exp(-5*(1+2*(log(c+eps)-log(cf))/log(cd))));
        %a = 5/log(10);
        %b = 5-a*log(cf);
        %G = 1./(1+exp(-a*log(c+eps)-b)); % Add eps to avoid log(0)
        % Alternative sigmoid variant. Does not go to zero for c->0
        %a=10;
        %b=2000;
        %G = exp(-a*exp(-b*c));
    else
        G = 1;
    end
    % This needs to be done last, otherwise F==2 where foam strength is
    % zero.
    %Fcap  = 1 - G.*(1-Fcap );
    %Fgrad = 1 - G.*(1-Fgrad);

    % Cap value at 1 for sW < sF. Mainly a problem for Fcap if F0<<1.
    Fcap = min(Fcap,1).*(sW >= sF) + (sW < sF);
    
    
    if opt.mobMultIsLinear
        varargout = {Fcap, Fgrad, G};
    else
        varargout = {Fcap + Fgrad.*vG};
    end
    
end
