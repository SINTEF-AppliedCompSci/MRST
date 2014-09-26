classdef TestSPE9 < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSPE9(varargin)
            test = test@ScheduleTest();
            
            
            test = merge_options(test, varargin{:});

            mrstModule add ad-fi deckformat
            
            fn = fullfile('SPE', 'SPE9', 'BENCH_SPE9.DATA');
            
            [deck, schedule, model, rock] = setupADcase(fn);
            
            gravity on
            
            p0  = deck.SOLUTION.PRESSURE;
            sw0 = deck.SOLUTION.SWAT;
            sg0 = deck.SOLUTION.SGAS;
            s0  = [sw0, 1-sw0-sg0, sg0];
            rs0 = deck.SOLUTION.RS;
            rv0 = 0;
            
            
            test.state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
            
            test.schedule = schedule;
            test.model = model;
            test.rock = rock;
            
            test.model.drsMaxRel = .2;
            test.model.dpMaxRel  = .2;
            test.model.dsMaxAbs  = .05;

        end
        function s = getIdentifier(test, name) %#ok
            s = [mfilename('class'), '_', name];
        end
    end
    
    methods (Test)
        % Add your own, test specific tests here
    end
    
end

