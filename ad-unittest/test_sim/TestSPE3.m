classdef TestSPE3 < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSPE3(varargin)
            test = test@ScheduleTest();
            
            
            test = merge_options(test, varargin{:});

            mrstModule add deckformat ad-fi ad-refactor
            
            fn = fullfile('SPE', 'SPE3', 'BENCH_SPE3.DATA');
            
            [G, rock, fluid, deck, schedule] = setupADcase(fn);
            
            % The case includes gravity
            gravity on
            
            p0  = deck.SOLUTION.PRESSURE;
            sw0 = deck.SOLUTION.SWAT;
            sg0 = deck.SOLUTION.SGAS;
            s0  = [sw0, 1-sw0-sg0, sg0];
            rv0 = deck.SOLUTION.RV;
            rs0 = 0;
            test.state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
            
            test.schedule = schedule;
            test.model = ThreePhaseBlackOilModel(G, rock, fluid, 'inputdata', deck);
            test.model.disgas = false;
        end
        function s = getIdentifier(test, name) %#ok
            s = [mfilename('class'), '_', name];
        end
    end
    
    methods (Test)
        % Add your own, test specific tests here
    end
    
end

