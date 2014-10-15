classdef TestSimplePolymer < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSimplePolymer(varargin)
            test = test@ScheduleTest();
            test = merge_options(test, varargin{:});

            mrstModule add deckformat ad-fi ad-blackoil ad-core

            fn = fullfile('SINTEF', 'simplePolymer', 'POLYMER.DATA');
            
            [deck, schedule, model, rock] = setupADcase(fn);

            gravity on
            
            G = model.G;
            ijk = gridLogicalIndices(G);
            state0 = initResSol(G, deck.PROPS.PVCDO(1), [ .9, .1]);
            state0.s(ijk{3} == 1, 2) = .9;
            state0.s(ijk{3} == 2, 2) = .8;
            state0.s(:,1) = 1 - state0.s(:,2);
            
            % Add zero polymer concentration to the state.
            state0.c    = zeros(G.cells.num, 1);
            state0.cmax = zeros(G.cells.num, 1);
            
            test.state0 = state0;
            
            test.schedule = schedule;
            test.model = model;
            test.rock = rock;
            
        end
        function s = getIdentifier(test, name) %#ok
            s = [mfilename('class'), '_', name];
        end
    end
    
    methods (Test)
        % Add your own, test specific tests here
    end
    
end

