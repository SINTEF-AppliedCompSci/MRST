function fluid_incomp = convertToIncompressibleFluid(model, varargin)
% Simple fluid-coverter setting pressure-dependent values constant  
%
% SYNOPSIS:
%   fluid_incomp = convertToIncompressibleFluid(model, pn1, pv1, ...)
%
% DESCRIPTION:
%  
% REQUIRED PARAMETERS:
%   model - model including AD-type fluid
%
% OPTIONAL PARAMETERS:
%   state    -  state structure which is used to compute mean pressure 
%              (and mean rs,rv if relevant)
%   pressure - pressure value or vector (used only if state is empty)
%   rs       - rs value or vector (used only if state is empty)
%   rv       - rv value or vector (used only if state is empty)

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('state',           [], ..., 
             'pressure', 200*barsa, ...
             'rs',             nan, ...
             'rv',             nan, ...
             'unitBFactors',  nan);

opt = merge_options(opt, varargin{:});

if ~isempty(opt.state)
    p = mean(opt.pressure);    
else
    p = mean(opt.pressure);
end

if model.disgas
    if ~isempty(opt.state)
        rs = mean(state.rs);
    else
        rs = mean(opt.rs);
    end
    assert(~isempty(rs), 'RS can''t be defaulted for model with disolved gas')
end

if model.vapoil
    if ~isempty(opt.state)
        rv = mean(state.rv);
    else
        rv = mean(opt.rv);
    end
    assert(~isempty(rs), 'RV can''t be defaulted for model with vaporized oil')
end

fluid = model.fluid;
fluid_incomp = fluid;

if model.disgas
    fluid_incomp.bO  = @(x1, x2)fluid.bO(p, rs) + 0*x1;
    fluid_incomp.muO = @(x1, x2)fluid.muO(p, rs) + 0*x1;
elseif model.oil
    fluid_incomp.bO  = @(x1)fluid.bO(p) + 0*x1;
    fluid_incomp.muO = @(x1)fluid.muO(p) + 0*x1;
end

if model.vapoil
    fluid_incomp.bG  = @(x1, x2)fluid.bG(p, rv) + 0*x1;
    fluid_incomp.muG = @(x1, x2)fluid.muG(p, rv) + 0*x1;
elseif model.gas
    fluid_incomp.bG  = @(x1)fluid.bG(p) + 0*x1;
    fluid_incomp.muG = @(x1)fluid.muG(p) + 0*x1;
end

if model.water
    fluid_incomp.bW  = @(x1)fluid.bW(p) + 0*x1;
    fluid_incomp.muW = @(x1)fluid.muW(p) + 0*x1;
end

if isfield(fluid, 'pvMultR')
    fluid_incomp.pvMultR = @(x1)fluid.pvMultR(p) + 0*x1;
end
    
end
