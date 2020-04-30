classdef SubdomainModel < WrapperModel
   
    properties
        mappings
        restrictionOperators
        noflowBC = false;
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = SubdomainModel(parent, cells, varargin)
            model = model@WrapperModel(parent);
            [model, submodelOpt] = merge_options(model, varargin{:});
            model = model.setSubModel(cells, submodelOpt{:});
            model.restrictionOperators = constructRestrictionOperators(model);
        end
        
        %-----------------------------------------------------------------%
        function model = setSubModel(model, cells, varargin)
            rmodel = model.getReservoirModel();
            [submodel, map] = getSubModel(rmodel, cells, varargin{:});
            model = model.setReservoirModel(model, submodel);
            model.G = model.parentModel.G;
            model.mappings = map;
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
                problem = transformProblem(problem, 'ML', ML, 'MR', MR, 'B', I, 'zero', ~keep);
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p0 = state.pressure;
            [state, report] = updateState@WrapperModel(model, state, problem, dx, drivingForces);
            % 
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
                problem.iterationNo = max(problem.iterationNo, problem.state.globalIteration);
            end
            [convergence, values, names] = checkConvergence@WrapperModel(model, problem, varargin{:});
        end
            
        %-----------------------------------------------------------------%
        function pmodel = setReservoirModel(model, pmodel, rmodel)
            if ~isa(pmodel.parentModel, 'ReservoirModel')
                pmodel.parentModel = model.setReservoirModel(pmodel.parentModel, rmodel);
                pmodel.parentModel.G = rmodel.G;
                return
            end
            pmodel.parentModel = rmodel;
        end

    end
end

%-------------------------------------------------------------------------%
function operators = constructRestrictionOperators(model)
    % Get logical mask for values we keep unchanged
    map  = model.mappings.cells;
    v    = map.internal | map.overlap;
    keep = v(map.keep);
    % Construct restriction operator
    nc = model.parentModel.G.cells.num;
    ML = sparse(1:nc, 1:nc, keep, nc, nc); 
    if isa(model.parentModel.AutoDiffBackend, 'DiagonalAutoDiffBackend')
        model = model.validateModel();
        rmodel = model.getReservoirModel();
        ncomp = rmodel.getNumberOfComponents();
        MR = sparse(1:nc*ncomp, 1:nc*ncomp, repmat(keep, ncomp, 1), nc*ncomp, nc*ncomp);
        I  = @(i,j) sparse(1:nc, (1:nc) + nc*(i-1), ~keep, nc, nc*ncomp);
    else
        MR = ML;
        I  = @(i,j) sparse(1:nc, 1:nc, ~keep, nc, nc).*(i == j);
    end
    % Make operator struct
    operators = struct('ML', ML, 'MR', MR, 'I', I, 'keep', keep);
end