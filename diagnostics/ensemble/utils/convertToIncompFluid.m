function fluid_incomp = convertToIncompFluid(model, varargin)
% Simple fluid-coverter setting pressure-dependent values constant  
%
% SYNOPSIS:
%   fluid_incomp = convertToIncompFluid(model, pn1, pv1, ...)
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

opt = struct('state',           [], ..., 
             'pressure', 200*barsa, ...
             'rs',             nan, ...
             'rv',             nan);

opt = merge_options(opt, varargin{:});

if ~isempty(opt.state)
    p = mean(opt.state.pressure);
else
    p = mean(opt.pressure);
end

if model.disgas
    if ~isempty(opt.state)
        rs = mean(opt.state.rs);
    else
        rs = mean(opt.rs);
    end
    assert(any(isnan(rs)), 'RS can''t be defaulted for model with disolved gas')
end

if model.vapoil
    if ~isempty(opt.state)
        rv = mean(opt.state.rv);
    else
        rv = mean(opt.rv);
    end
    assert(any(isnan(rv)), 'RV can''t be defaulted for model with vaporized oil')
end

fluid = model.fluid;
fluid_incomp = fluid;

if model.water
    fluid_incomp.bW  = setConstant(fluid.bW, p);
    fluid_incomp.muW = setConstant(fluid.muW, p);
end

if model.oil
    arg = {p};
    if model.disgas
        arg = [arg, {rs}];
    end
    fluid_incomp.bO  = setConstant(fluid.bO, arg{:});
    fluid_incomp.muO = setConstant(fluid.muO, arg{:});
end

if model.gas
    arg = {p};
    if model.vapoil
        arg = [arg, {rv}];
    end
    fluid_incomp.bG  = setConstant(fluid.bG, arg{:});
    fluid_incomp.muG = setConstant(fluid.muG, arg{:});
end

if isfield(fluid, 'pvMultR')
    fluid_incomp.pvMultR = setConstant(fluid.pvMultR, p); 
end
    
end


function fn_incomp = setConstant(fn, varargin)
if ~iscell(fn)
    if nargin == 1
        fn_incomp = @(x)fn(varargin{:})+0*x;
    else
        fn_incomp = @(x1, x2)fn(varargin{:})+0*x1;
    end
else
    fn_incomp = cell(1, numel(fn));
    for k = 1:numel(fn)
        if nargin == 1
            fn_incomp{k} = @(x)fn{k}(varargin{:})+0*x;
        else
            fn_incomp{k} = @(x1, x2)fn{k}(varargin{:})+0*x1;
        end
    end
end
end
