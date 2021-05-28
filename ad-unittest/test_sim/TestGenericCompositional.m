classdef TestGenericCompositional < matlab.unittest.TestCase
    properties (TestParameter)
        includeWater = {true, false};
        backend = {'diagonal', 'diagonal-mex', 'sparse'};
        fluidSystem = {'simple'};
        modelType = {'natural',...
                     'natural-legacy',...
                     'overall',...
                     'overall-si', ...
                     'overall-legacy',...
                     'natural-legacy-si', ...
                     'overall-legacy-si'};
        phases = {'OG', 'GO', 'WO', 'WG'};                 
    end
    
    methods
        function test = TestGenericCompositional()
            mrstModule reset
            mrstModule add ad-unittest ad-core ad-blackoil ad-props compositional sequential
        end

        function [state0, model, schedule] = buildTestCase(test, modelType, includeWater, fluidSystem, backend, varargin)
            opt = struct('dims', [2, 2, 2], ...
                         'physDims', [1000, 1000, 10], ...
                         'gravity', true, ...
                         'phases', [], ...
                         'useAMGCL', false);
            opt = merge_options(opt, varargin{:});
            ok_sim = true;
            if opt.gravity
                gravity reset on
            else
                gravity reset off
            end
            defaultedPhases = strcmp(opt.phases, 'OG');

            G = cartGrid(opt.dims, opt.physDims);
            G = computeGeometry(G);
            rock = makeRock(G, 1*darcy, 0.5);

            nph = 2 + includeWater;
            uniform = true;
            if uniform
                inj = ones(1, nph);
            else
                inj = 1:nph;
            end
            inj = inj./sum(inj);
            
%             s = nph:-1:1;
            s = ones(1, nph);
            s = s./sum(s);
            
            [compfluid, info] = getBenchmarkMixture(fluidSystem);
            p = info.pressure;
            p_w = p/2;

            eos = EquationOfStateModel(G, compfluid);
%             state0 = initCompositionalState(G, info.pressure, info.temp, s, info.initial, eos);
            nc = G.cells.num;
            e = @(x) repmat(x, nc, 1);
            state0 = struct('pressure', e(info.pressure), 'T', e(info.T), 's', e(s), 'components', e(info.initial));
            if uniform
                nkr = [2, 2, 2];
                mu = [1, 1, 1]*centi*poise;
                rhoS = [10, 10, 10];
                c = [1e-7, 1e-7, 1e-7]/barsa;
            else
                nkr = [2, 3, 4];
                mu = [1, 5, 0.1]*centi*poise;
                rhoS = [1000, 10, 10];
                c = [0, 1e-7, 1e-4]/barsa;
            end
            ph = 'wog';
            act = [includeWater, true, true];
            
            fluid = initSimpleADIFluid('phases', ph(act), ...
                                       'mu', mu(act), ...
                                       'rho', rhoS(act), ...
                                       'cR', 1e-12/barsa, ...
                                       'c', c(act), ...
                                       'n', nkr(act));
            if any(opt.phases == 'O')
                fluid = rmfield(fluid, 'bO');
                fluid = rmfield(fluid, 'muO');
            end
            if any(opt.phases == 'G')
                fluid = rmfield(fluid, 'bG');
                fluid = rmfield(fluid, 'muG');
            end
            if any(opt.phases == 'W')
                if includeWater
                    fluid = rmfield(fluid, 'bW');
                    fluid = rmfield(fluid, 'muW');
                else
                    ok_sim = false;
                end
            end
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
            isLegacy = false;
            isSeq = false;
            switch lower(modelType)
                case 'natural'
                    model = GenericNaturalVariablesModel(arg{:});
                case 'natural-legacy'
                    model = NaturalVariablesCompositionalModel(arg{:});
                    isLegacy = true;
                case 'overall'
                    model = GenericOverallCompositionModel(arg{:});
                case 'overall-legacy'
                    model = OverallCompositionCompositionalModel(arg{:});
                    isLegacy = true;
                case 'natural-legacy-si'
                    model = NaturalVariablesCompositionalModel(arg{:});
                    model = getSequentialModelFromFI(model);
                    isLegacy = true;
                    isSeq = true;
                case 'overall-legacy-si'
                    model = OverallCompositionCompositionalModel(arg{:});
                    model = getSequentialModelFromFI(model);
                    isLegacy = true;
                    isSeq = true;
                case 'natural-si'
                    parent = GenericNaturalVariablesModel(arg{:});
                    pmodel = PressureModel(parent);
                    tmodel = TransportModel(parent);
                    model = SequentialPressureTransportModel(pmodel, tmodel);
                    isSeq = true;
                case 'overall-si'
                    parent = GenericOverallCompositionModel(arg{:});
                    pmodel = PressureModel(parent);
                    tmodel = TransportModel(parent);
                    model = SequentialPressureTransportModel(pmodel, tmodel);
                    isSeq = true;
                otherwise
                    error('%s is unknown', modelType);
            end

            time = 1*year;
            pv = poreVolume(G, rock);
            irate = 0.1*sum(pv)/time;

            z = min(G.cells.centroids(:, 3));
            W = [];
            W = verticalWell(W, G, rock, 1, 1, 1, 'comp_i', inj, ...
                            'val', irate, 'type', 'rate', 'refDepth', z);
            W = verticalWell(W, G, rock, opt.dims(1), opt.dims(2), opt.dims(3), ...
                            'comp_i', inj, 'val', p_w, 'type', 'bhp', 'refDepth', z);
            for i = 1:numel(W)
                W(i).components = info.injection;
            end
            L = opt.phases(1);
            V = opt.phases(2);
            if isSeq
                if isLegacy
                    model.pressureModel.liquidPhase = L;
                    model.pressureModel.vaporPhase = V;
                    
                    model.pressureModel.liquidPhase = L;
                    model.transportModel.vaporPhase = V;
                else
                    model.pressureModel.parentModel.liquidPhase = L;
                    model.transportModel.parentModel.vaporPhase = V;
                    
                    model.pressureModel.parentModel.liquidPhase = L;
                    model.transportModel.parentModel.vaporPhase = V;
                end
            else
                model.liquidPhase = L;
                model.vaporPhase = V;
            end
            ok_sim = ok_sim && ~(isLegacy && ~defaultedPhases);
            if ok_sim
                schedule = simpleSchedule(time, 'W', W);
            else
                schedule = [];
            end
        end
        
        function [ws, states, reports] = testMultiPhase(test, modelType, includeWater, fluidSystem, backend, phases, varargin)
            fprintf('Testing fluid "%s" with %s model and %s backend\n', fluidSystem, modelType, backend);
            fprintf('Water: %d\n', double(includeWater));
            fprintf('Phases: %s\n', phases);
            [state0, model, schedule] = test.buildTestCase(modelType, includeWater, fluidSystem, backend, 'phases', phases, varargin{:});
            if ~isempty(schedule)
                [ws, states, reports] = simulateScheduleAD(state0, model, schedule);
            end
        end
    end
    
    methods (Test)
        function varargout = compositionalTest(test, modelType, includeWater, fluidSystem, backend, phases)
            varargout = cell(nargout, 1);
            [varargout{:}] = testMultiPhase(test, modelType, includeWater, fluidSystem, backend, phases);
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
