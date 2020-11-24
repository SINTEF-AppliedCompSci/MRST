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
    %   'name' - Name of the reservoir state variable we are interested in.
    %            Default: 'sW'.
    %   'cells' - range of cells that we are interested in. Default: ':'
    %   'time' - At what simulation time we want to store the QoI.
    % 
    % SEE ALSO:
    %   `BaseQoI`, `WellQoI`, `MRSTExample`, `BaseSamples`
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
            % Reads the simulation state solutions from the given problem
            % and extract the relevant data.
            %
            % SYNOPSIS:
            %   u = qoi.computeQoI(problem)
            %
            % PARAMETERS:
            %   problem - The specific problem for which to compute the QoI
            
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
        function n = norm(qoi, u)
            n = cellfun(@(u) sum(abs(u).*qoi.dx), u);
            if ~isnan(qoi.dt)
                n = n.*qoi.dt;
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
                         'cellNo'    , 1          , ...
                         'subCellNo' , 1);
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
            if opt.title
                title(sprintf('%s after %d simulation days', qoi.name, qoi.time(opt.cellNo)/day));
            end
            setup.setAxisProperties(gca);
        end
        
    end
    
    methods (Access = protected)
        %-----------------------------------------------------------------%
        function n = findTimestep(qoi, schedule, time)
            % Utility function to get the solution index of the requested
            % simulation time
            [~, n] = min(abs(cumsum(schedule.step.val) - time));
        end
        
        %-----------------------------------------------------------------%
        function u = getStateValue(qoi, model, state)
            % Utility to function to extract the correct field from a given
            % state object
            u = model.getProp(state, qoi.name);
        end
        

    end
    
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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