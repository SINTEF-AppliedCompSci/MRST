classdef TestGenericCompositional < matlab.unittest.TestCase
    properties (TestParameter)
        includeWater = {true, false};
        backend = {'diagonal', 'diagonal-mex', 'sparse'};
        fluidSystem = {'simple'};
        modelType = {'natural',...
                     'natural-legacy',...
                     'overall',...
                     'natural-si', ...
                     'overall-si', ...
                     'overall-legacy',...
                     'natural-legacy-si', ...
                     'overall-legacy-si'};
    end
    
    methods
        function test = TestGenericCompositional()
            mrstModule reset
            mrstModule add ad-unittest ad-core ad-blackoil ad-props compositional blackoil-sequential
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
            rhoS = [1000, 10, 10];
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
                    model = GenericNaturalVariablesModel(arg{:});
                case 'natural-legacy'
                    model = NaturalVariablesCompositionalModel(arg{:});
                case 'overall'
                    model = GenericOverallCompositionModel(arg{:});
                case 'overall-legacy'
                    model = OverallCompositionCompositionalModel(arg{:});
                case 'natural-legacy-si'
                    model = NaturalVariablesCompositionalModel(arg{:});
                    model = getSequentialModelFromFI(model);
                case 'overall-legacy-si'
                    model = OverallCompositionCompositionalModel(arg{:});
                    model = getSequentialModelFromFI(model);
                case 'natural-si'
                    parent = GenericNaturalVariablesModel(arg{:});
                    pmodel = PressureModel(parent);
                    tmodel = TransportModel(parent);
                    model = SequentialPressureTransportModel(pmodel, tmodel);
                case 'overall-si'
                    parent = GenericOverallCompositionModel(arg{:});
                    pmodel = PressureModel(parent);
                    tmodel = TransportModel(parent);
                    model = SequentialPressureTransportModel(pmodel, tmodel);
                otherwise
                    error('%s is unknown', modelType);
            end

            time = 1*year;
            pv = poreVolume(G, rock);
            irate = 0.1*sum(pv)/time;

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

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
