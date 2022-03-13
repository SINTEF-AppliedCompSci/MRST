classdef ReservoirStateQoI < BaseQoI
     % Class for extracting state related quantity of interests from (an
     % ensemble of) reservoir simulation problems.
    %
    % DESCRIPTION:
    %   This class is used within a MRSTEnsemble to extract, store, and 
    %   work with quantities of interest related to the reservoir state
    %
    % SYNOPSIS
    %   qoi = ReservoirStateQoI(...);
    %
    % OPTIONAL PARAMETERS
    %   'names' - Name of the reservoir state variable we are interested in.
    %            Default: 'sW'.
    %   'cells' - range of cells that we are interested in. Default: ':'
    %   'time' - At what simulation time we want to store the QoI.
    % 
    % SEE ALSO:
    %   `BaseQoI`, `WellQoI`, `MRSTExample`, `BaseSamples`
    properties
        cells = inf
        time  = inf
        
        dt = nan
        dx
    end
    
    methods
        %-----------------------------------------------------------------%
        function qoi = ReservoirStateQoI(varargin)
            qoi = qoi@BaseQoI('names', 'sW');
            qoi = merge_options(qoi, varargin{:});
            qoi.plotAllSamples = false;
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            % Check that the configs that are inserted to the constructor
            % makes sense for the base problem for the ensemble, and
            % updates remaining fields.
            %
            % SYNOPSIS:
            %   qoi = qoi.validateQoI(problem)
            %
            % PARAMETERS:
            %   problem - MRST problem that is representative for the
            %             (ensemble member) problems that will be used with
            %             this QoI instance.
            %
            % NOTE:
            %   When used in an MRSTEnsemble, this function is called by
            %   the ensemble constructor.
            qoi = validateQoI@BaseQoI(qoi, problem);
            qoi.dx = problem.SimulatorSetup.model.G.cells.volumes;
            if isinf(qoi.time)
                qoi.time = cumsum(problem.SimulatorSetup.schedule.step.val);
            end
            qoi.time = reshape(qoi.time, 1, []);
            if numel(qoi.time) == 1
                if qoi.time == -1
                    qoi.time = sum(problem.SimulatorSetup.schedule.step.val);
                end
            else
                qoi.dt = diff([0, qoi.time]);
            end
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            % Reads the simulation state solutions from the given problem
            % and extract the relevant data.
            %
            % SYNOPSIS:
            %   u = qoi.computeQoI(problem)
            %
            % PARAMETERS:
            %   problem - The specific problem for which to compute the QoI
            
            states   = problem.OutputHandlers.states;
            model    = problem.SimulatorSetup.model.setupStateFunctionGroupings();
            schedule = problem.SimulatorSetup.schedule;
            nt = numel(qoi.time);
            subs = qoi.cells;
            if isinf(qoi.cells), subs = 1:model.G.cells.num; end
            nc = numel(subs);
            u = cell2struct(repmat({zeros(nc, nt)}, numel(qoi.names), 1), qoi.names);
            for i = 1:nt
                n       = qoi.findTimestep(schedule, qoi.time(i));
                state   = states{n};
                for fn = qoi.names
                    prop = model.getProps(state, fn{1});
                    u.(fn{1})(:,i) = prop(subs,:);
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function n = norm(qoi, u)
            n = sum(abs(u).*qoi.dx,1);
            if ~isnan(qoi.dt)
                n = bsxfun(@times, n, qoi.dt);
            end
            n = sum(n);
        end
        
        %-----------------------------------------------------------------%
        function plotQoI(qoi, ensemble, u, varargin)
            % Plot a single well QoI u in current figure.
            opt = struct('isMean'    , false      , ...
                         'threshold' , [-inf, inf], ...
                         'contour'   , false      , ...
                         'labels'    , true       , ...
                         'title'     , true       , ...
                         'numLines'  , 10         , ...
                         'names'     , {qoi.names});
            [opt, extra] = merge_options(opt, varargin{:});
            if iscell(u)
                cellfun(@(u) qoi.plotQoI(u), u);
                return
            end
            if isa(ensemble, 'MRSTEnsemble') || isa(ensemble, 'MRSTEnsembleV2')
                setup = ensemble.setup;
            elseif isa(ensemble, 'TestCase') 
                setup = ensemble;
            else
                error(['Input ensemble must either be an instance of ', ...
                      'class ''MRSTEnsemble'' or ''TestCase'''     ]);
            end
            for i = 1:numel(opt.names)
                ui = u.(opt.names{i});
                for j = 1:size(ui,2)
                    subplot(1,size(ui,2), j);
                    if ~opt.contour
                        subs = ui(:,j) >= opt.threshold(1) & ui(:,j) <= opt.threshold(2);
                        keep = qoi.cells & subs;
                        plotCellData(setup.model.G, ui(keep,j), keep, extra{:});
                    else
                        assert(setup.model.G.griddim == 2 && ...
                               any(strcmpi(setup.model.G.type, 'cartGrid')), ...
                               'Contour plotting requires Cartesian 2D grid')
                        uc = nan(setup.model.G.cells.num, 1);
                        uc(qoi.cells) = ui(:,j);
                        uc = reshape(uc, setup.model.G.cartDims);
                        xMin = min(setup.model.G.cells.centroids);
                        xMax = max(setup.model.G.cells.centroids);
                        x = linspace(xMin(1), xMax(1), setup.model.G.cartDims(1));
                        y = linspace(xMin(2), xMax(2), setup.model.G.cartDims(2));
                        color = [1,1,1]*0.8*(1-opt.isMean); % Plot mean in distinct color
                        contour(x, y, uc, opt.numLines, extra{:}, 'color', color);
                    end
                    if opt.title
                        title(sprintf('%s after %d simulation days', qoi.names{i}, qoi.time(j)/day));
                    end
                    setup.setAxisProperties(gca);
                end
            end
        end
        
    end
    
    methods (Access = protected)
        %-----------------------------------------------------------------%
        function n = findTimestep(qoi, schedule, time)
            % Utility function to get the solution index of the requested
            % simulation time
            [~, n] = min(abs(cumsum(schedule.step.val) - time));
        end
        
    end
    
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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