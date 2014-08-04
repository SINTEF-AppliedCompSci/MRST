classdef TestEGG < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestEGG(varargin)
            test = test@ScheduleTest();
            test = merge_options(test, varargin{:});

            mrstModule add deckformat ad-fi ad-refactor

            fn = fullfile('external', 'TUDelft-EGG', 'BENCH_EGG.DATA');
            
            [deck, schedule, model] = setupADcase(fn);

            gravity on

            % Approximate initial conds:
            pr   = 400*barsa;
            rz   = G.cells.centroids(1,3);
            dz   = G.cells.centroids(:,3) - rz;
            rhoO    = fluid.bO(400*barsa)*fluid.rhoOS;
            rhoW    = fluid.bW(400*barsa)*fluid.rhoWS;
            rhoMix  = .1*rhoW + .9*rhoO;
            p0   = pr + norm(gravity)*rhoMix*dz;      
            
            test.state0 = initResSol(G, p0, [0.1, .90]);
            test.schedule = schedule;
            test.model = model;
            
        end
        function s = getIdentifier(test, name) %#ok
            s = [mfilename('class'), '_', name];
        end
    end
    
    methods (Test)
        % Add your own, test specific tests here
    end
    
end

