function uplift = computeUpliftForState(model, state, topnode, varargin)
%
%
% SYNOPSIS:
%   function uplift = computeUpliftForState(model, state, topnode, varargin)
%
% DESCRIPTION: Compute the vertical displacement at the node given by
% topnode. This function is used in the example runAdjointExample
%
% PARAMETERS:
%   model    - poroelastic model (fully coupled type e.g. MechWaterModel or MechOilWaterModel)
%   state    - state structure, which should have the same structure as the one
%              used in poroelastic model (e.g. MechWaterModel)
%   topnode  - index of the node where the uplift will be computed
%   varargin - 
%
% RETURNS:
%   uplift - value of the vertical displacement at the node topnode
%
% EXAMPLE:
%
% SEE ALSO:
%
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


    opt = struct('ComputePartials', false);
    opt = merge_options(opt, varargin{:});

    
    % get index in the mechanical displacement field xd where the vertical
    % displacement of the top node is stored.
    G = model.G;
    nx = G.cartDims(1);
    ny = G.cartDims(2);
    isdirdofs = model.mechModel.operators.isdirdofs;
    u = (1 : (G.griddim*G.nodes.num))';
    indlift = G.griddim*(topnode - 1) + 2;
    u = u(~isdirdofs);
    indlift = find(u == indlift);
    

    wellSol = state.wellSol;
        
    % This is specific to oil water fluid model
    % TODO: extend to other fluid models
    [p, sW, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                 'water', ...
                                                 'wellSol', ...
                                                 'xd');
    [wellVars, wellVarNames, wellMap] = ...
        model.fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);
    
    if opt.ComputePartials
        [p, sW, wellVars{:}, xd] = initVariablesADI(p, sW, wellVars{:}, xd);
    end
    uplift = xd(indlift);

    
end