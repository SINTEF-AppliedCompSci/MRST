classdef TestGenericBO < matlab.unittest.TestCase
    properties (TestParameter)
        phases = {'wog', 'wo', 'og', 'w', 'o', 'g'};
        backend = {'diagonal', 'diagonal-mex', 'sparse'};
    end
    
    methods
        function test = TestGenericBO()
            mrstModule reset
            mrstModule add ad-unittest ad-core ad-blackoil ad-props
        end

        function testMultiPhase(test, phases, backend, varargin)
            gravity reset on
            opt = struct('dims', [2, 2, 2], ...
                         'physDims', [1000, 1000, 10], ...
                         'useAMGCL', false);
            opt = merge_options(opt, varargin{:});

            G = cartGrid(opt.dims, opt.physDims);
            G = computeGeometry(G);
            rock = makeRock(G, 1*darcy, 0.5);

            nph = numel(phases);

            inj = 1:nph;
            inj = inj./sum(inj);
            
            s = nph:-1:1;
            s = s./sum(s);
            
            p = 100*barsa;
            p_w = p/2;
            
            state0 = initResSol(G, p, s);
            
            nkr = [2, 3, 4];
            mu = [1, 5, 0.1]*centi*poise;
            rhoS = [1000, 700, 10];
            c = [0, 1e-7, 1e-4]/barsa;
            
            ph = 'wog';
            act = ismember(ph, phases);
            
            fluid = initSimpleADIFluid('phases', ph(act), ...
                                       'mu', mu(act), ...
                                       'rho', rhoS(act), ...
                                       'c', c(act), ...
                                       'n', nkr(act));
            
            switch lower(backend)
                case 'diagonal'
                    auto = DiagonalAutoDiffBackend('useMex', false);
                case 'diagonal-mex'
                    auto = DiagonalAutoDiffBackend('useMex', true);
                case 'sparse'
                    auto = AutoDiffBackend();
                otherwise
                    error('Unknown backend %s', opt.backend);
            end
            
            model = GenericBlackOil(G, rock, fluid, ...
                'disgas', false, 'vapoil', false, ...
                'water', any(phases == 'w'), ...
                'oil', any(phases == 'o'), ...
                'gas', any(phases == 'g'), ...
                'AutoDiffBackend', auto);
            
            time = 1*year;
            irate = sum(model.operators.pv)/time;
            
            
            W = [];
            W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', inj, ...
                            'val', irate, 'type', 'rate');
            W = verticalWell(W, G, rock, opt.dims(1), opt.dims(1), [], ...
                            'comp_i', inj, 'val', p_w, 'type', 'bhp');
            
            schedule = simpleSchedule(time, 'W', W);
            
            [ws, states, rep] = simulateScheduleAD(state0, model, schedule);
        end
    end
    
    methods (Test)
        function immiscibleMultiPhaseTest(test, phases, backend)
            testMultiPhase(test, phases, backend);
        end
    end
end
