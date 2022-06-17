function W = addThermalWellProps(W, G, rock, fluid, varargin)
%Add thermal properties to an existing well structure.
% 
% SYNOPSIS: 
%   W = addThermalWellProps(W, G, rock, fluid, 'pn1', pv1)
% 
% PARAMETERS: 
%   W - Well structure created with e.g. addWell.
% 
%   G, rock, fluid - Grid, rock, and fluid structure of the model
%
% OPTIONAL PARAMETERS:
%   T    - Well temperature (scalar or one per well). Will only be used
%            for injeciton wells. Default: 40 C.
%
%   WIth - Thermal well index used for computing conductive heat flux from
%          the wellbore to the perforated cell. Computed by Peaceman
%          approximation if left empty (see `computeWellIndex.m`).
% 
% RETURNS:
%   W - valid well structure with temperature and thermal well index 
% 
% SEE ALSO:
% 'addWell'.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('T'   , convertFromCelcius(50), ...
                 'WIth', []                    );
    opt = merge_options(opt, varargin{:});
    
    fNames = setdiff(fieldnames(opt), 'WIth');
    nW = numel(W);
    for fNo = 1:numel(fNames)
        fn = fNames{fNo};
        v = opt.(fn);
        if numel(v) < nW
            v = repmat(v, nW, 1);
        end
        for wNo = 1:nW
            W(wNo).(fn) = v(wNo);
        end
    end
    
    W = computeThermalWellIndex(W, G, rock, fluid, opt);
    
end

%-------------------------------------------------------------------------%
function W = computeThermalWellIndex(W, G, rock, fluid, opt)

    lambdaR = rock.lambdaR.*(1-rock.poro);
    lambdaF = fluid.lambdaF.*rock.poro;
    lambda  = lambdaR + lambdaF;

    givenWIth = ~isempty(opt.WIth);
    [W.WIth] = deal(W.WI);
    for i = 1:numel(W)
        cells  = W(i).cells;
        if givenWIth
            W(i).WIth = repmat(WIth, numel(cells),1);
            continue
        end
        % Assume thermal conductivity tensors are diagonal for now
        rock   = struct('perm', ones(numel(value(lambda)), 1));
        radius = W(i).r;
        if isfield(G.faces, 'nodePos')
            % Compute well index
            % Ensure equivalent radius greater than well radius
            [dx, dy, dz] = cellDims(G, cells);
            re = 2*0.14*sqrt(dx.^2 + dy.^2)/2;
            C  = max(1.1*radius(1)./re,1);
            dx = [dx.*C, dy.*C, dz];
        else
            dx = [];
        end
        wi = computeWellIndex(G, rock, radius, cells, 'cellDims', dx);
        wi = wi.*lambda(cells);
        W(i).WIth = wi;
    end
    
end