classdef ScheduleConsistencyTest < matlab.unittest.TestCase
    properties
    end
    
    methods
    end
    
    methods (Test)
       function [ws, s, rep] = ChangingWellNumbers(test)
           mrstModule reset
            mrstModule add deckformat ad-fi ad-core ad-props ad-unittest ad-blackoil

            G = cartGrid([11 1 3]);
            G = computeGeometry(G);

            fluid = initSimpleADIFluid();
            ons = ones(G.cells.num, 1);
            rock = struct('perm', ons*darcy, 'poro', ons);

            model = ThreePhaseBlackOilModel(G, rock, fluid);

            state = initResSol(G, 0, [1 0 0]);
            state.rs = 0;
            state.rv = 0;
            
            [ii, jj, kk] = gridLogicalIndices(G);

            W = [];
            W = addWell(W, G, rock, find(ii == 1), 'Val', 1*barsa, 'sign', -1);
            W = addWell(W, G, rock, find(ii == 6), 'Val', 2*barsa, 'sign', 1);
            W = addWell(W, G, rock, find(ii == 11), 'Val', 3*barsa, 'sign', 1);

            schedule = struct();

            schedule.step.val = ones(3, 1)*day;
            schedule.step.control = [1; 2; 3];
            
            % Make three controls that are the subsets of the same superset
            W_1 = W([1, 3]);
            W_2 = W([1, 2, 3]);
            W_3 = W([2, 3]);

            schedule.control(1).W = W_1;
            schedule.control(2).W = W_2;
            schedule.control(3).W = W_3;
            
            % Strip a perforation from one of the control's first well
            subs = [1; 3];
            w = schedule.control(1).W(1);
            w.dZ = w.dZ(subs);
            w.cells = w.cells(subs);
            w.WI = w.WI(subs);
            w.cstatus = w.cstatus(subs);
            
            schedule.control(1).W(1) = w;
            schedule = makeScheduleConsistent(schedule);

            % This  will probably error out if the code isn't working
            [ws, s, rep] = simulateScheduleAD(state, model, schedule);
       end
    end
end

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
