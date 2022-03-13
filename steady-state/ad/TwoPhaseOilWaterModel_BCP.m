classdef TwoPhaseOilWaterModel_BCP < TwoPhaseOilWaterModel
    % Two phase oil/water system without dissolution
properties
    
end

methods
    
    function model = TwoPhaseOilWaterModel_BCP(G, rock, fluid, varargin)
        %model = model@TwoPhaseOilWaterModel(G, rock, fluid, varargin{:});
        warning('off', 'mrst:ReservoirModel');
        model = model@TwoPhaseOilWaterModel(G, [], [], varargin{:});
        warning('on', 'mrst:ReservoirModel');
        model.rock  = rock;
        model.fluid = fluid;
        
        isPeriodic = isfield(G, 'parent');
        if isPeriodic
            % Periodic grids must be set up in a special way
            N  = double(G.faces.neighbors);
            intInx = (prod(N,2)~=0);
            N  = N(intInx, :);
            
            % Transmissibility
            T  = computeTrans(G.parent, rock);
            hfmap = G.cells.faces(:,end);
            T  = T(hfmap);
            cf = G.cells.faces(:,1); 
            nf = G.faces.num;
            T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
            T_all = T;
            T  = T(intInx);
            
            s = setupOperatorsTPFA(G.parent, model.rock, ...
                'deck', model.inputdata, 'neighbors', N, 'trans', T);
            s.T = T;
            s.T_all = T_all;
            s.internalConn = intInx;
            
            model.operators = s;
        else
            model.operators = setupOperatorsTPFA(G, model.rock, ...
                'deck', model.inputdata);
        end
        
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsOilWater_BCP(state0, state, model,...
                        dt, ...
                        drivingForces,...
                        varargin{:});

    end
    
    % --------------------------------------------------------------------%
    function forces = getValidDrivingForces(model)
        forces = getValidDrivingForces@TwoPhaseOilWaterModel(model);
        % Support for periodic boundary conditions
        forces.bcp = [];
    end
end
end

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
