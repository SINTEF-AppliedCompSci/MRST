classdef TestSimpleOW < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSimpleOW(varargin)
            test = test@ScheduleTest();
            test = merge_options(test, varargin{:});

            mrstModule add deckformat ad-fi ad-blackoil ad-core

            fn = fullfile('SINTEF', 'simpleOW', 'simple10x1x10.data');
            
            [deck, schedule, model, rock] = setupADcase(fn);

            gravity on
            
            test.state0 = initResSol(model.G, deck.PROPS.PVCDO(1), [.15, .85]);
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

