classdef DualContinuumWaterModel < DualContinuumReservoirModel
%
% SYNOPSIS:
%   model = DualContinuumWaterModel(G, rock, fluid, varargin)
%
% DESCRIPTION: 
%   Single phase water model that will be required by the fully coupled 
%   and fixed stress dual-continuum water models.
%
%
% PARAMETERS:
%   G     - Simulation grid.
%   rock  - Rock cell for fracture / matrix
%   fluid - Fluid cell for fracture / matrix
%
% RETURNS:
%   class instance
%
% EXAMPLE: 
%
% SEE ALSO: 
%
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
    properties

    end
    
    methods
        function model = DualContinuumWaterModel(G, rock, fluid, varargin)
            model = model@DualContinuumReservoirModel(G, rock, fluid);
            % Only enable water
            model.oil = false;
            model.gas = false;
            model.water = true;
            model.useCNVConvergence = false;   
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            % we can just keep the equations for the single-phase instance, as
            % they will not be used anyway. 
            [problem, state] = equationsWater(state0, state, model,...
                                              dt, drivingForces,...
                                              varargin{:});           
        end
        
        function [frac_index, mat_index] = findDCWells(model, wellSol)
        % By definition the dual-continuum model assumes flow from both
        % porosity systems. This function helps to find the wells in
        % wellSol that correspond to matrix and fracture continua. It
        % enforces that when we add in a well, we add a well for BOTH
        % continua. Rows of the wellSol struct that correspond to fracture
        % and matrix wells are denoted frac_index and mat_index resp. 
            if ~ isempty(wellSol) 
                frac_index = find(contains({wellSol.name}, 'frac'));
                mat_index = find(contains({wellSol.name}, 'mat'));
                assert(~ isempty(frac_index) && ~ isempty(mat_index), ...
                    ['well names are required to include fracture and matrix labels e.g. '...
                     'Prod_frac & Prod_mat. If a single well connection is required '...
                     'for example to a fracture, create a dummy "matrix" well with a '...
                     'dummy structure containing a zero permeability field entry'])
            else
                % dummy indices
                frac_index = 1; mat_index = 2;
            end
        end
        
    end
end

