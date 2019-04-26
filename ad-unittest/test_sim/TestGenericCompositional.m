classdef TestGenericCompositional < matlab.unittest.TestCase
    properties (TestParameter)
        includeWater = {true, false};
        backend = {'diagonal', 'diagonal-mex', 'sparse'};
        fluidSystem = {'simple'};
        modelType = {'natural','natural-legacy'};
    end
    
    methods
        function test = TestGenericCompositional()
            mrstModule reset
            mrstModule add ad-unittest ad-core ad-blackoil ad-props compositional
        end

        function [state0, model, schedule] = buildTestCase(test, modelType, includeWater, fluidSystem, backend, varargin)
            opt = struct('dims', [2, 2, 2], ...
                         'physDims', [1000, 1000, 10], ...
                         'gravity', true, ...
                         'useAMGCL', false);
            opt = merge_options(opt, varargin{:});
            if opt.gravity
                gravity reset on
            else
                gravity reset off
            end

            G = cartGrid(opt.dims, opt.physDims);
            G = computeGeometry(G);
            rock = makeRock(G, 1*darcy, 0.5);

            nph = 2 + includeWater;

            inj = 1:nph;
            inj = inj./sum(inj);
            
            s = nph:-1:1;
            s = s./sum(s);
            
            [compfluid, info] = getCompositionalFluidCase(fluidSystem);
            p = info.pressure;
            p_w = p/2;

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
            
            arg = {G, rock, fluid, eos, ...
                'water', includeWater, ...
                'AutoDiffBackend', auto};
            switch lower(modelType)
                case 'natural'
                    model = GenericNaturalVariables(arg{:});
                case 'natural-legacy'
                    model = NaturalVariablesCompositionalModel(arg{:});
                otherwise
                    error('%s is unknown', modelType);
            end

            time = 1*year;
            irate = 0.1*sum(model.operators.pv)/time;

            W = [];
            W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', inj, ...
                            'val', irate, 'type', 'rate');
            W = verticalWell(W, G, rock, opt.dims(1), opt.dims(2), [], ...
                            'comp_i', inj, 'val', p_w, 'type', 'bhp');
            for i = 1:numel(W)
                W(i).components = info.injection;
            end
            schedule = simpleSchedule(time, 'W', W);
        end
        
        function [ws, states, reports] = testMultiPhase(test, modelType, includeWater, fluidSystem, backend, varargin)
            fprintf('Testing fluid "%s" with %s model and %s backend\n', fluidSystem, modelType, backend);
            [state0, model, schedule] = test.buildTestCase(modelType, includeWater, fluidSystem, backend, varargin{:});
            [ws, states, reports] = simulateScheduleAD(state0, model, schedule);
        end
    end
    
    methods (Test)
        function varargout = compositionalTest(test, modelType, includeWater, fluidSystem, backend)
            varargout = cell(nargout, 1);
            [varargout{:}] = testMultiPhase(test, modelType, includeWater, fluidSystem, backend);
        end
    end
end
