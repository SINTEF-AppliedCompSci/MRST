classdef ReservoirStateQoI < BaseQoI
    
    properties
        name = 'sW'
        cells = ':'
        time  = inf
        
        dt = nan
        dx
    end
    
    methods
        %-----------------------------------------------------------------%
        function qoi = ReservoirStateQoI(varargin)
            qoi = qoi@BaseQoI();
            qoi = merge_options(qoi, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            qoi = validateQoI@BaseQoI(qoi, problem);
            qoi.dx = problem.SimulatorSetup.model.G.cells.volumes;
            if strcmpi(qoi.time, ':')
                qoi.time = cumsum(problem.SimulatorSetup.schedule.step.val);
            end
            if numel(qoi.time) == 1
                if isinf(qoi.time)
                    qoi.time = sum(problem.SimulatorSetup.schedule.step.val);
                end
            else
                qoi.dt = diff([0,qoi.time]);
            end
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            nt = numel(qoi.time);
            u  = cell(nt,1);
            states   = problem.OutputHandlers.states;
            model    = problem.SimulatorSetup.model.setupStateFunctionGroupings();
            schedule = problem.SimulatorSetup.schedule;
            for i = 1:nt
                n     = qoi.findTimestep(schedule, qoi.time(i));
                state = states{n};
                u{i}  = qoi.getStateValue(model, state);
            end
        end
        
        %-----------------------------------------------------------------%
        function n = findTimestep(qoi, schedule, time)
            [~, n] = min(abs(cumsum(schedule.step.val) - time));
        end
        
        %-----------------------------------------------------------------%
        function u = getStateValue(qoi, model, state)
            u = model.getProp(state, qoi.name);
        end
        
        %-----------------------------------------------------------------%
        function n = norm(qoi, u)
            n = cellfun(@(u) sum(abs(u).*qoi.dx), u);
            if ~isnan(qoi.dt)
                n = n.*qoi.dt;
            end
            n = sum(n);
        end
        
        %-----------------------------------------------------------------%
        function plotQoI(qoi, ensemble, u, varargin)
            opt = struct('isMean'   , false      , ...
                         'threshold', [-inf, inf], ...
                         'contour'  , false      , ...
                         'numLines' , 10         );
            [opt, extra] = merge_options(opt, varargin{:});
            if iscell(u)
                cellfun(@(u) qoi.plotQoI(u), u);
                return
            end
            if isa(ensemble, 'MRSTEnsemble')
                setup = ensemble.setup;
            elseif isa(ensemble, 'MRSTExample')
                setup = ensemble;
            else
                error(['Input ensemble must either be an instance of ', ...
                      'class ''MRSTEnsemble'' or ''MRSTExample'''     ]);
            end
            if ~opt.contour
                subs = u >= opt.threshold(1) & u <= opt.threshold(2);
                keep = qoi.cells & subs;
                plotCellData(setup.model.G, u(keep), keep, extra{:});
            else
                assert(setup.model.G.griddim == 2 && ...
                       any(strcmpi(setup.model.G.type, 'cartGrid')), ...
                       'Contour plotting requires Cartesian 2D grid')
                uc = nan(setup.model.G.cells.num, 1);
                uc(qoi.cells) = u;
                uc = reshape(uc, setup.model.G.cartDims);
                xMin = min(setup.model.G.cells.centroids);
                xMax = max(setup.model.G.cells.centroids);
                x = linspace(xMin(1), xMax(1), setup.model.G.cartDims(1));
                y = linspace(xMin(2), xMax(2), setup.model.G.cartDims(2));
                color = [1,1,1]*0.8*(1-opt.isMean); % Plot mean in distinct color
                contour(x, y, uc, opt.numLines, extra{:}, 'color', color);
            end
            setup.setAxisProperties(gca);
        end
        
    end
    
end