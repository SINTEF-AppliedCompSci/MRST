classdef SubdomainModel < WrapperModel
    % Simulation-ready model constructed for a subset of a full model.
    % Problem equations are constructed assuming Dirichlet BCs for all
    % boundary faces that are not also a boundary face of the full model
    properties
        mappings              % Mappings from 
        restrictionOperators  % Operators for restricting equations
        prolongationOperators % Operators for prolongation equations
        noflowBC = false;     % Optionally use no-flow BCs
        implicitType
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = SubdomainModel(parent, cells, varargin)
            % Construct subdomain model
            model = model@WrapperModel(parent);
            [model, submodelOpt] = merge_options(model, varargin{:});
            % Set the submodel
            model = model.setSubModel(cells, submodelOpt{:});
            % Construct restriction operators
            model.restrictionOperators  = constructRestrictionOperators(model);
            model.prolongationOperators = constructProlongationOperators(model);
            if isa(parent, 'TransportModel')
                model.implicitType = parent.implicitType;
            end
        end
        
        %-----------------------------------------------------------------%
        function model = setSubModel(model, cells, varargin)
            % Set submodel
            rmodel = getReservoirModel(model);
            % Get submodel of reservoir model from subset
            [submodel, map] = getSubModel(rmodel, cells, varargin{:});
            % Update reservoir model
            model = setReservoirModel(model, submodel);
            model.G = submodel.G; % Replace grid
            model.mappings = map; % Set mappings
        end
        
        %-----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, forces, varargin)
            % Let parent model assemble equations
            [problem, state] = model.parentModel.getEquations(state0, state, dt, forces, varargin{:});
            % Replace rhs with 0 and Jacobian with identity matrix for all
            % external cells
            keep = model.restrictionOperators.keep;
            if numel(keep) == numel(value(state.pressure)) && ~model.noflowBC
                ML      = model.restrictionOperators.ML;
                MR      = model.restrictionOperators.MR;
                I       = model.restrictionOperators.I;
                problem = transformProblem(problem, 'ML', ML, 'MR', MR, 'B', I);
            end
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            [eqs, names, types, state] = model.parentModel.getModelEquations(state0, state, dt, drivingForces);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p0 = state.pressure;
            [state, report] = updateState@WrapperModel(model, state, problem, dx, drivingForces);
            % Recompute relative pressure change using global pressure
            % range in denominator if given
            if isfield(state, 'pressureRange')
                if isfield(state, 'dpRel')
                    state.dpRel = (state.pressure - p0)./state.pressureRange;
                end
                pmodel = model.parentModel;
                isPressure = isa(model, 'PressureModel');
                if isfield(state, 'pressureChange') && isPressure ...
                        && strcmpi(pmodel.pressureIncTolType, 'relative')
                    state.pressureChange = (state.pressure - p0)./state.pressureRange;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            if isfield(problem.state, 'globalIteration') && 0
                % We may want to use the global iteration when checking for
                % subdomain convergence (currently disabled)
                problem.iterationNo = max(problem.iterationNo, problem.state.globalIteration);
            end
            [convergence, values, names] = checkConvergence@WrapperModel(model, problem, varargin{:});
        end
            
        function [state, names, origin] = getStateAD(model, state, init)
            [state, names, origin] = model.parentModel.getStateAD(state, init);
        end
%         %-----------------------------------------------------------------%
%         function pmodel = setReservoirModel(model, pmodel, rmodel)
%             % Set reservoir model to the model
%             if ~isa(pmodel.parentModel, 'ReservoirModel')
%                 pmodel.parentModel = model.setReservoirModel(pmodel.parentModel, rmodel);
%                 pmodel.parentModel.G = rmodel.G;
%                 return
%             end
%             pmodel.parentModel = rmodel;
%         end

    end
end

%-------------------------------------------------------------------------%
function operators = constructRestrictionOperators(model)
    % Get logical mask for values we keep unchanged
    map  = model.mappings.cells;
    v    = map.internal | map.overlap;
    keep = v(map.keep);
    % Construct restriction operator
    rmodel = getReservoirModel(model);
    nc = rmodel.G.cells.num;
    ML = sparse(1:nc, 1:nc, keep, nc, nc); 
    if isa(model.parentModel.AutoDiffBackend, 'DiagonalAutoDiffBackend')
        model = model.validateModel();
        if isa(model.parentModel, 'PressureModel')
            ncomp = 1;
        else
            rmodel = getReservoirModel(model);
            ncomp = rmodel.getNumberOfComponents();
        end
        MR = sparse(1:nc*ncomp, 1:nc*ncomp, repmat(keep, ncomp, 1), nc*ncomp, nc*ncomp);
        I  = @(i,j) sparse(1:nc, (1:nc) + nc*(i-1), ~keep, nc, nc*ncomp);
    else
        MR = ML;
        I  = @(i,j) sparse(1:nc, 1:nc, ~keep, nc, nc).*(i == j);
    end
    % Make operator struct
    operators = struct('ML', ML, 'MR', MR, 'I', I, 'keep', keep);
end

%-------------------------------------------------------------------------%
function operators = constructProlongationOperators(model)
    % Get logical mask for values we keep unchanged
    map  = model.mappings.cells;
    v    = map.internal | map.overlap;
    keep = v(map.keep);
    % Construct restriction operator
    rmodel = getReservoirModel(model);
    Nc = numel(map.keep);
    nc = rmodel.G.cells.num;
%     nc = model.parentModel.G.cells.num; 
    ML0 = sparse(find(map.keep), 1:nc, keep, Nc, nc);
    ML  = sparse(find(map.keep), 1:nc, 1, Nc, nc);
    if isa(model.parentModel.AutoDiffBackend, 'DiagonalAutoDiffBackend')
        model = model.validateModel();
        if isa(model.parentModel, 'PressureModel')
            ncomp = 1;
        else
            rmodel = getReservoirModel(model);
            ncomp = rmodel.getNumberOfComponents();
        end
        jj  = find(repmat(map.keep, ncomp, 1));
        v   = repmat(keep, ncomp, 1);
        MR0 = sparse(1:nc*ncomp, jj, v, nc*ncomp, Nc*ncomp);
        MR  = sparse(1:nc*ncomp, jj, 1, nc*ncomp, Nc*ncomp);
    else
        MR0 = ML0';
        MR  = ML';
    end
    % Make operator struct
    operators = struct('ML', ML, 'ML0', ML0, 'MR', MR, 'MR0', MR0);
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
