classdef TestGenericCompositional < matlab.unittest.TestCase
    properties (TestParameter)
        includeWater = {true, false};
        backend = {'diagonal', 'diagonal-mex', 'sparse'};
        fluidSystem = {'simple'};
        modelType = {'natural'};
    end
    
    methods
        function test = TestGenericCompositional()
            mrstModule reset
            mrstModule add ad-unittest ad-core ad-blackoil ad-props compositional
        end

        function testMultiPhase(test, modelType, includeWater, fluidSystem, backend, varargin)
            gravity reset on
            opt = struct('dims', [2, 2, 2], ...
                         'physDims', [1000, 1000, 10], ...
                         'useAMGCL', false);
            opt = merge_options(opt, varargin{:});

            G = cartGrid(opt.dims, opt.physDims);
            G = computeGeometry(G);
            rock = makeRock(G, 1*darcy, 0.5);

            nph = 2 + includeWater;

            inj = 1:nph;
            inj = inj./sum(inj);
            
            s = nph:-1:1;
            s = s./sum(s);
            
            p = 100*barsa;
            p_w = p/2;
            
            
            [compfluid, info] = getCompositionalFluidCase(fluidSystem);
            eos = EquationOfStateModel(G, compfluid);
            state0 = initCompositionalState(G, info.pressure, info.temp, s, info.initial, eos);
            
            nkr = [2, 3, 4];
            mu = [1, 5, 0.1]*centi*poise;
            rhoS = [1000, 700, 10];
            c = [0, 1e-7, 1e-4]/barsa;
            
            ph = 'wog';
            act = [includeWater, true, true];
            
            fluid = initSimpleADIFluid('phases', ph(act), ...
                                       'mu', mu(act), ...
                                       'rho', rhoS(act), ...
                                       'cR', 1e-12/barsa, ...
                                       'c', c(act), ...
                                       'n', nkr(act));
            fluid = rmfield(fluid, 'bO');
            fluid = rmfield(fluid, 'bG');
            fluid = rmfield(fluid, 'muO');
            fluid = rmfield(fluid, 'muG');
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
            
            model = GenericNaturalVariables(G, rock, fluid, eos, ...
                'water', includeWater, ...
                'AutoDiffBackend', auto);
            
            time = 1*year;
            irate = sum(model.operators.pv)/time;
            
            
            W = [];
            W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', inj, ...
                            'val', irate, 'type', 'rate');
            W = verticalWell(W, G, rock, opt.dims(1), opt.dims(1), [], ...
                            'comp_i', inj, 'val', p_w, 'type', 'bhp');
            for i = 1:numel(W)
                W(i).components = info.injection;
            end
            schedule = simpleSchedule(time, 'W', W);
            
            [ws, states, rep] = simulateScheduleAD(state0, model, schedule);
        end
    end
    
    methods (Test)
        function compositionalTest(test, modelType, includeWater, fluidSystem, backend)
            testMultiPhase(test, modelType, includeWater, fluidSystem, backend);
        end
    end
end
