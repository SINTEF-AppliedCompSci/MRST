classdef ScheduleConsistencyTest < matlab.unittest.TestCase
    methods (Test)
       function testReorderStrategies(test)
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

           [ii, ~, ~] = gridLogicalIndices(G);

           W = [];
           W = addWell(W, G, rock, find(ii == 1),  'Val', 1*barsa, 'sign', -1);
           W = addWell(W, G, rock, find(ii == 6),  'Val', 2*barsa, 'sign',  1);
           W = addWell(W, G, rock, find(ii == 11), 'Val', 3*barsa, 'sign',  1);

           schedule = struct();
           schedule.step.val = ones(3, 1)*day;
           schedule.step.control = [1; 2; 3];

           W_1 = W([1, 3]);
           W_2 = W([1, 2, 3]);
           W_3 = W([2, 3]);

           schedule.control(1).W = W_1;
           schedule.control(2).W = W_2;
           schedule.control(3).W = W_3;

           % Strip a perforation from the first control's first well
           subs = [1; 3];
           w = schedule.control(1).W(1);
           w.dZ = w.dZ(subs);
           w.cells = w.cells(subs);
           w.WI = w.WI(subs);
           w.cstatus = w.cstatus(subs);
           schedule.control(1).W(1) = w;

           strategies = {'origin', 'depth', 'depth-origin', 'none', 'direction-xyz'};

           for k = 1:numel(strategies)
               strat = strategies{k};
               s = schedule; 
               s = makeScheduleConsistent(s, 'G', G, 'ReorderStrategy', strat);
               % Should run without error and produce outputs
               [ws, simState, report] = simulateScheduleAD(state, model, s);

               test.verifyNotEmpty(ws);
               test.verifyNotEmpty(simState);
               test.verifyNotEmpty(report);
           end
       end
    end
end