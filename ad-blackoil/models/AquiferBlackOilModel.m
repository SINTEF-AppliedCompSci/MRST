classdef AquiferBlackOilModel < ThreePhaseBlackOilModel
    
    properties
        aquifer % aquifer data table
        aquind % structure which explains the fields in the aquifer data table
    end
    
    methods
        function model = AquiferBlackOilModel(G, rock, fluid, aquifer, aquind, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});
            model.aquifer = aquifer;
            model.aquind  = aquind;
        end
        
        % --------------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            switch(lower(name))
              case {'aquiferpressures', 'aquifervolumes'}
                fn = lower(name);
                index = 1;
              otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ThreePhaseBlackOilModel(model, name);
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = updateAfterConvergence@ThreePhaseBlackOilModel(model, ...
                                                              state0, state, ...
                                                              dt, ...
                                                              drivingForces);
            sW = model.getProp(state, 'sW');
            p = model.getProp(state, 'pressure');
            q = computeAquiferFluxes(model, p, sW, state, dt);
        end
        
        function eqs = addAquiferContribution(model, eqs, names, state, p, sW, dt)
            
            q = computeAquiferFluxes(model, p, sW, state, dt);
            wind = strcmpi('water', names);
            conn = model.aquifer.conn;
            eqs{wind}(conn) = eqs{wind}(conn) - q;

        end
        
    end
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
