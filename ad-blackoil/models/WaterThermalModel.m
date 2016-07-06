classdef WaterThermalModel < ReservoirModel
    % Two phase oil/water system without dissolution
    properties

    end
    
    methods
        function model = WaterThermalModel(G, rock, fluid, varargin)
            
            model = model@ReservoirModel(G, rock, fluid);
            
            % This is the model parameters for oil/water
            model.oil = false;
            model.gas = false;
            model.water = true;
            %model.addflux=true;
            
            % Blackoil -> use CNV style convergence 
            model.useCNVConvergence = false;
            
            model.saturationVarNames = {'sw'};
            model.wellVarNames = {'qWs', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
            rock_heat=struct('perm',rock.lambdaR);
            T_r=computeTrans(G,rock_heat);
            cf = G.cells.faces(:,1);
            nf = G.faces.num;
            T_r  = 1 ./ accumarray(cf, 1./T_r, [nf, 1]);
            model.operators.T_r_all=T_r;
            intInx = all(G.faces.neighbors ~= 0, 2);
            model.operators.T_r = T_r(intInx);
            % Setup operators
            %model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function forces = getValidDrivingForces(model)
        forces = getValidDrivingForces@ReservoirModel(model);
        %
        forces.bcT = [];
    end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsWaterThermal(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Parent class handles almost everything for us
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            
            % Update wells based on black oil specific properties
            saturations = model.saturationVarNames;
            wi = strcmpi(saturations, 'sw');
            oi = strcmpi(saturations, 'so');
            gi = strcmpi(saturations, 'sg');

            W = drivingForces.W;
            state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);

        end
    end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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