classdef TestSPE1 < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSPE1(varargin)
            test = test@ScheduleTest();
            
            
            test = merge_options(test, varargin{:});

            mrstModule add deckformat ad-fi ad-refactor
            
            moddir = mrstPath('query', 'ad-unittest');
            fn = fullfile(moddir, 'data', 'b1-SPE1', 'BENCH_SPE1.DATA');
            
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
            
            p0  = deck.SOLUTION.PRESSURE;
            sw0 = deck.SOLUTION.SWAT;
            sg0 = deck.SOLUTION.SGAS;
            s0  = [sw0, 1-sw0-sg0, sg0];
            rs0 = deck.SOLUTION.RS;
            rv0 = 0;
            test.state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
            
            test.schedule = convertDeckScheduleToMRST(G, rock, deck);
            test.model = ThreePhaseBlackOilModel(G, rock, fluid, 'inputdata', deck);
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
            test.runSchedule(name, 'useCPR', true, 'useAGMG', false);
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

