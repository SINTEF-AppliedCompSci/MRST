classdef AquiferBlackOilModel < ThreePhaseBlackOilModel
    
    properties
        aquifers % aquifers data table
        aquind % structure which explains the fields in the aquifers data table
        aquiferprops
        modeltype
    end
    
    methods
        function model = AquiferBlackOilModel(G, rock, fluid, aquifers, aquind, ...
                                              aquiferprops, varargin)
            opt = struct('modeltype', 'blackoil');
            [opt, rest] = merge_options(opt, varargin{:});
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, rest{:});

            model.modeltype = opt.modeltype;
            switch model.modeltype
              case 'blackoil'
                model.oil = true;
                model.gas = true;
                model.water = true;
              case 'oilwater'
                model.oil = true;
                model.gas = false;
                model.water = true;
              otherwise
                error('modeltype not recognized');
            end
            
            model.aquifers     = aquifers;
            model.aquind       = aquind;
            model.aquiferprops = aquiferprops;
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ThreePhaseBlackOilModel(model);
            forces.aquifer = true;
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            switch model.modeltype
              case 'blackoil'
                [problem, state] = equationsBlackOil(state0, state, model, dt, ...
                                                     drivingForces, ...
                                                     varargin{:});
              case 'oilwater'
                [problem, state] = equationsOilWater(state0, state, model, dt, ...
                                                     drivingForces, ...
                                                     varargin{:});
              otherwise
                error('modeltype not recognized');
            end

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
            sW   = model.getProp(state, 'sW');
            p    = model.getProp(state, 'pressure');
            p_aq = model.getProp(state, 'aquiferpressures');
            V_aq = model.getProp(state, 'aquifervolumes');
            
            q = computeAquiferFluxes(model, p, sW, state, dt);
            
            aquifers     = model.aquifers;
            aquind       = model.aquind;
            aquiferprops = model.aquiferprops;
            conn = aquifers(:, aquind.conn);
            aquid = aquifers(:, aquind.aquid);
            nconn = size(conn, 1);
            naq = max(aquid);
            aquid2conn = sparse(aquid, (1 : nconn)', 1, naq, nconn)';
           
            Q = aquid2conn'*q;
            Q = dt.*Q;
            C = aquiferprops.C;
            p_aq = p_aq - Q./(C.*V_aq);
            V_aq = V_aq + Q;
        
            state = model.setProp(state, 'aquiferpressures', p_aq);
            state = model.setProp(state, 'aquifervolumes'  , V_aq);
        end
        
        function eqs = addAquifersContribution(model, eqs, names, state, p, sW, dt)
            
            q = computeAquiferFluxes(model, p, sW, state, dt);
            wind = strcmpi('water', names);
            aquifers = model.aquifers;
            aquind   = model.aquind;
            conn = aquifers(:, aquind.conn);
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
