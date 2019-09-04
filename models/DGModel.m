classdef DGModel < WrapperModel
   
    properties
        disc         = [];
        tryMaxDegree = true;
    end
    
    methods

        % ----------------------------------------------------------------%
        function model = DGModel(parentModel, varargin)
            model = model@WrapperModel(parentModel);

            % If we use reordering, this tells us which cells are actually
            % part of the discretization, and which cells that are included
            % to get fluxes correct
            G = parentModel.G;
            G.cells.ghost = false(parentModel.G.cells.num,1);
            model.G = G;
            [model, discArgs] = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, discArgs{:});
            end
            
        end

        % ----------------------------------------------------------------%
        function state = validateState(model, state)
            state = model.parentModel.validateState(state);
            state = assignDofFromState(model.disc, state);    
        end
        
        % ----------------------------------------------------------------%
        function gdxyz = getGravityGradient(model)
            gdxyz = model.parentModel.getGravityGradient();
        end
        
        % ----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            [problem, state] ...
                = transportEquationOilWaterDG(state0, state, model, dt, drivingForces, ...
                                  'solveForOil'  , model.conserveOil  , ...
                                  'solveForWater', model.conserveWater, ...
                                  varargin{:}                         );
        end
        
        % ----------------------------------------------------------------%
        function [p, state] = getProp(model, state, name)
            
            [fn, index] = model.getVariableField(name, false);
            if isempty(fn)
                % Not known - check property functions
                ok = false;
                containers = model.getStateFunctionGroupings();
                for i = 1:numel(containers)
                    c = containers{i};
                    nms = c.getNamesOfStateFunctions();
                    sub = strcmpi(nms, name);
                    if any(sub)
                        pdof = c.get(model, state, nms{sub});
                        ok = true;
                    end
                end
                if ~ok
                    error('PhysicalModel:UnknownVariable', ...
                        'Unknown variable field %s', name);
                end
            else
                if iscell(state.(fn))
                    if ischar(index)
                        pdof = state.(fn);
                    else
                        pdof = state.(fn){index};
                    end
                else
                    pdof = state.(fn)(:, index);
                end
            end
            
            cells = (1:model.G.cells.num)';
            
            [W, x, c] = model.disc.getCubature(cells, 'volume');
            [x, ~, scaling] = model.disc.transformCoords(x, c);
            p = model.disc.evaluateDGVariable(x, c, state, pdof);
            
        end
        
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch lower(name)
                case {'pressure'}
                    index = 1;
                    fn = 'pressuredof';
                case {'water', 'sw'}
                    index = 1;
                    fn = 'sdof';
                case {'oil', 'so'}
                    index = 2;
                    fn = 'sdof';
                case {'gas', 'sg'}
                    index = 3;
                    fn = 'sdof';
                case{'sdof'}
                    index = ':';
                    fn = 'sdof';
                otherwise
                    [fn, index] = model.parentModel.getVariableField(name, varargin{:});
            end
        end
        
    end
    
end