classdef DualContinuumReservoirModel < DualPorosityReservoirModel
% Base dual-continuum flow models
%
% SYNOPSIS:
%   model = DualContinuumReservoirModel(G, rock, fluid)
%
% DESCRIPTION:
%   Base flow class for all reservoirs that show dual-permeability behaviour.
%
% REQUIRED PARAMETERS:
%   G     - Simulation grid.
%   rock  - Rock cell for fracture / matrix
%   fluid - Fluid cell for fracture / matrix
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   DualPorosityReservoirModel, DualContMechFluidModel,
%   DualContMechFluidSplitModel

    properties    
        operators_matrix % required to consider flow in matrix
    end

    methods
        % --------------------------------------------------------------------%
        function model = DualContinuumReservoirModel(G, rock, fluid, varargin)
            % Parent class handles everything for us
            model = model@DualPorosityReservoirModel(G, rock, fluid);
            model.operators_matrix = setupOperatorsTPFA(G, model.rock_matrix, 'deck', model.inputdata);
        end
    end
end

%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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

