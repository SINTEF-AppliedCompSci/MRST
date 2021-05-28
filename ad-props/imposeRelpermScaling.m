function model = imposeRelpermScaling(model, varargin)
% Utility to impose rel-perm end-point scaling for a model not associated 
% with a deck. 
%
% SYNOPSIS:
%   model = imposeRelpermScaling(model, ...)
%
% REQUIRED PARAMETERS:
%  model - Standard model where model.fluid has field 'krPts' present.
%          Oil rel-perm functions should have names 'krOW'/'krOG'.  
%
% OPTIONAL PARAMETERS:
%   'nPoints'   - 2/3-point scaling (default: 2)
%   '[kw]'      - where kw is any of the fields SWL, SWCR, SWU, SGL, SGCR, 
%                 SGU, SOWCR, SOGCR, KRW, KRO, KRG. Value is either scalar 
%                 or nCells. 
% RETURNS:
% model         - model with updated rock.krscale and the subsequent 
%                 required fields of a 'fake' deck in model.inputdata.
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
assert(isfield(model.fluid, 'krPts'), ...
    'To impose rel-perm scaling, the fluid must contain field ''krPts''');

opt   = struct('nPoints',  2);
[opt, argScale] = merge_options(opt, varargin{:});
assert(any(opt.nPoints == [2,3]), 'Only 2- or 3-point scaling is supported')

scale = struct(argScale{:});
valid = {'SWL', 'SWCR', 'SWU', 'SGL', 'SGCR', 'SGU', 'SOWCR', 'SOGCR', ...
         'KRW', 'KRO', 'KRG'};
     
fn      = fieldnames(scale);
fn      = cellfun(@upper, fn, 'UniformOutput', false);
isValid = cellfun(@(x)ismember(x, valid), fn);
if ~all(isValid)
    warnProblem(prob)
    fn = fn(isValid);
end

nc = model.G.cells.num;
for k = 1:numel(fn)
    nv = numel(scale.(fn{k}));
    assert(any(nv == [1,nc]), ...
        'Scaling values does not match number of grid cells');
    if numel(scale.(fn{k})) == 1  % expand single value
        scale.(fn{k}) = scale.(fn{k})*ones(nc,1);
    end
end

% setup model
model.rock.krscale = initRelpermScaling(struct('PROPS', scale), nc);
% add fake deck to model to get correct setup of FlowPropertyFunctions
str = 'NO';
if opt.nPoints == 3
    str = 'YES';
end
model.inputdata = struct(...
    'RUNSPEC', struct('ENDSCALE', {{'NODIR', 'REVERS', 1, 20, 0}}), ...
    'PROPS',   struct('SCALECRS', {{str}}), ...
    'GRID', [], 'SOLUTION', []);

% check rel-perm names
if model.oil &&  isfield(model.fluid, 'krO')
    fluid = model.fluid;
    if ~model.gas && model.water && ~isfield(fluid, 'krOW')
        fluid.krOW = fluid.krO;
        if isfield(fluid.krPts, 'o')
            fluid.krPts.ow = fluid.krPts.o;
        end
    end
    if ~model.water && model.gas && ~isfield(fluid, 'krOG')
        fluid.krOG = fluid.krG;
        if isfield(fluid.krPts, 'o')
            fluid.krPts.og = fluid.krPts.o;
        end
    end
    if (~model.water || isfield(fluid, 'krOW')) && ...
       (~model.gas   || isfield(fluid, 'krOG')) 
        fluid = rmfield(fluid, 'krO');  
    else
        warning('Scaling of ''krO''-function ignored');
    end
    model.fluid = fluid;
end
end

%--------------------------------------------------------------------------

function warnProblem(prob)
for k = 1:numel(prob)
    warning('Ignoring unrecognized/unsupported scaling keyword: %s', prob{k});
end
end



