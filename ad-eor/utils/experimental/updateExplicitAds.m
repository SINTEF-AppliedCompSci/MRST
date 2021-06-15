function state = updateExplicitAds(state, model, varargin)
%
%
% SYNOPSIS:
%   function state = updateExplicitAds(state, model, varargin)
%
% DESCRIPTION: Update of the adsorption which ensures total mass
% conservation. Used with an explicit transport solver.
%
% PARAMETERS:
%   state    - State at current time-step.
%   model    - Model instance
%   varargin - Optional arguments
%
% RETURNS:
%   state - State at current time-step with updated adsorption

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

    opt = struct('desorptionThreshold', -inf);
    opt = merge_options(opt, varargin{:});

    s = model.operators; % shortcut
    fluid = model.fluid; % shortcut

    [p, sW, cs, csmax, ads, adsmax] = model.getProps(state, 'pressure', 'water', 'surfactant', ...
                                                           'surfactantmax', 'ads', 'adsmax');

    new_ads = computeEffAds(c, opt.desorptionThreshold, fluid);

    % Multipliers for properties
    [pvMult, ~, ~] = getMultipliers(model.fluid, p, p);
    bW   = fluid.bW(p);
    poro = model.rock.poro;

    % Solve the equation
    %
    %   (s.pv/dt).*(pvMult.*bW.*sW.*(new_c - c) + fluid.rhoRSft.*((1-poro)./poro).*(new_ads - ads))= 0
    %
    % We remove values where sW is zero. At those points, the concentration is not defined and we
    % let it unchanged.
    %
    ind = sW > 1e-12;
    coef1 = fluid.rhoRSft.*((1-poro)./poro);
    coef2 = pvMult.*bW.*sW;
    coef = coef1./coef2;
    coef(~ind) = 0; % if sW=0, concentration is in fact not well-defined and we let it unchanged.
    new_c  = c - coef.*(new_ads - ads);

    negc = new_c < 0;
    % For the terms where the concentration becomes negative, we set the adsorption term so that
    % surfactant mass conservation is ensured, that is
    %
    %   (s.pv/dt).*(pvMult.*bW.*sW.*(0 - c) + fluid.rhoRSft.*((1-poro)./poro).*(new_ads - ads) )= 0
    %
    if any(negc)
        new_ads(negc) = ads(negc) + coef2(negc)./coef1(negc).*c(negc);
        new_c = max(0, new_c);
    end

    state = model.setProp(state, 'ads', new_ads);
    state = model.setProp(state, 'adsmax', max(new_ads, adsmax));
    state = model.setProp(state, 'surfactant', max(new_c, 0));
    state = model.setProp(state, 'surfactantmax', max(new_cs, csmax));

end
