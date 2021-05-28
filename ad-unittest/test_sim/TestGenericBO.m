classdef TestGenericBO < matlab.unittest.TestCase
    properties (TestParameter)
        strategy = {'fi', 'fi-legacy', 'si', 'si-legacy'};
        phases = {'wog', 'wo', 'og', 'w', 'o', 'g'};
        backend = {'diagonal', 'diagonal-mex', 'sparse'};
    end
    
    methods
        function test = TestGenericBO()
            mrstModule reset
            mrstModule add ad-unittest ad-core ad-blackoil ad-props sequential
        end

        function [ws, states, rep] = testMultiPhase(test, strategy, phases, backend, varargin)
            gravity reset on
            opt = struct('dims', [2, 2, 2], ...
                         'physDims', [1000, 1000, 10], ...
                         'pvi', 1, ...
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
                                       'cR', 1e-12/barsa, ...
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
            
            water = any(phases == 'w');
            oil = any(phases == 'o');
            gas = any(phases == 'g');
            
            arg = {G, rock, fluid, 'disgas', false, 'vapoil', false, 'AutoDiffBackend', auto};
            switch strategy
                case 'fi'
                    model = GenericBlackOilModel(arg{:}, ...
                        'water', water, ...
                        'oil', oil, ...
                        'gas', gas);
                    pv = model.operators.pv;
                case 'fi-legacy'
                    model = getLegacy(water, oil, gas, arg);
                    if ~isa(model, 'PhysicalModel')
                        return
                    end
                    pv = model.operators.pv;
                case 'si'
                    model = GenericBlackOilModel(arg{:}, ...
                        'water', water, ...
                        'oil', oil, ...
                        'gas', gas);
                    pv = model.operators.pv;

                    if water + oil + gas  < 2
                        return
                    end
                    pmodel = PressureModel(model);
                    tmodel = TransportModel(model);
                    model = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
                    
                case 'si-legacy'
                    if water + oil + gas  < 2
                        return
                    end
                    model_fi = getLegacy(water, oil, gas, arg);
                    if ~isa(model_fi, 'PhysicalModel')
                        return
                    else
                        model = getSequentialModelFromFI(model_fi);
                    end
                    model.pressureModel.AutoDiffBackend = auto;
                    model.transportModel.AutoDiffBackend = auto;
                    pv = model.pressureModel.operators.pv;
            end
            time = 1*year;
            irate = opt.pvi*sum(pv)/time;
            
            
            W = [];
            W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', inj, ...
                            'val', irate, 'type', 'rate');
            W = verticalWell(W, G, rock, opt.dims(1), opt.dims(2), [], ...
                            'comp_i', inj, 'val', p_w, 'type', 'bhp');
            
            schedule = simpleSchedule(time, 'W', W);
            
            [ws, states, rep] = simulateScheduleAD(state0, model, schedule);
        end
    end
    
    methods (Test)
        function immiscibleMultiPhaseTest(test, strategy, phases, backend)
            testMultiPhase(test, strategy, phases, backend);
        end
    end
end

function model = getLegacy(water, oil, gas, arg)
    if water && oil && gas
        model = ThreePhaseBlackOilModel(arg{:});
    elseif water && oil
        model = TwoPhaseOilWaterModel(arg{:});
    elseif water && gas
        mrstModule add co2lab
        model = TwoPhaseWaterGasModel(arg{:});
    elseif water
        model = WaterModel(arg{:});
    else
        disp('Pseudocomponent combination not implemented in legacy solvers. Test automatically successful.')
        model = nan;
        return
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
