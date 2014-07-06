classdef TestSimpleOW < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSimpleOW(varargin)
            test = test@ScheduleTest();
            test = merge_options(test, varargin{:});

            mrstModule add deckformat ad-fi ad-refactor
            
            moddir = mrstPath('query', 'ad-unittest');
            fn = fullfile(moddir, 'data', 'simpleOW', 'simple10x1x10.data');
            
            deck = readEclipseDeck(fn);
            
            % The deck is given in field units, MRST uses metric.
            deck = convertDeckUnits(deck);
            
            G = initEclipseGrid(deck);
            G = computeGeometry(G);
            
            rock  = initEclipseRock(deck);
            rock  = compressRock(rock, G.cells.indexMap);
            
            % Create a special ADI fluid which can produce differentiated fluid
            % properties.
            fluid = initDeckADIFluid(deck);
            
            % The case includes gravity
            gravity on
            
            test.state0 = initResSol(G, deck.PROPS.PVCDO(1), [.15, .85]);
            
            test.schedule = convertDeckScheduleToMRST(G, rock, deck);
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

