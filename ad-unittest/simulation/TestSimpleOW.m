classdef TestSimpleOW < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSimpleOW(varargin)
            test = test@ScheduleTest();
            test = merge_options(test, varargin{:});

            mrstModule add deckformat ad-fi ad-refactor

            fn = fullfile('SINTEF', 'simpleOW', 'simple10x1x10.data');
            
            [G, rock, fluid, deck, schedule] = test.setupADcase(fn);

            gravity on
            
            test.state0 = initResSol(G, deck.PROPS.PVCDO(1), [.15, .85]);
            test.schedule = schedule;
            test.model = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);
            
        end
        function s = getIdentifier(test, name) %#ok
            s = [mfilename('class'), '_', name];
        end
    end
    
    methods (Test)
        function baseline(test)
            name = test.getIdentifier('baseline');
            test.runSchedule(name);
        end
        
        function CPR_mldivide(test)
            name = test.getIdentifier('cpr_mldivide');
            test.runSchedule(name, 'useCPR', true, 'useAGMG', false, 'cprTol', 1-3);
        end
                
        function CPR_AGMG(test)
            name = test.getIdentifier('cpr_agmg');
            mrstModule add agmg
            try
                agmg(speye(3), ones(3, 1));
            catch
                test.verifyFail( ...
                    'AGMG is not installed properly, test cannot proceed')
                return
            end
            test.runSchedule(name, 'useCPR', true, 'useAGMG', true);
        end
    end
    
end

