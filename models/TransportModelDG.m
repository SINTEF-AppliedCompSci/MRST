classdef TransportModelDG < TransportModel
    
    properties
        disc = [];
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = TransportModelDG(parent, varargin)
           
            model = model@TransportModel(parent);
           
            [model, discArgs] = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, discArgs{:});
            end
        end

        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            
            [fn, index] = getVariableField@TransportModel(model, name, varargin{:});
            fn = [fn, 'dof'];
            
        end
        
        % ----------------------------------------------------------------%
        function state = validateState(model, state)
            state    = assignDofFromState(model.disc, state);    
            state    = validateState@TransportModel(model, state);
            state.sT = sum(state.s,2);
        end
        
        %-----------------------------------------------------------------%
        function [p, state] = getProp(model, state, name)
            % Get dofs
            [pdof, state] = getProp@TransportModel(model, state, name);
            % Get cubature
            cells = (1:model.G.cells.num)';
            [W , x, cellNo] = model.disc.getCubature(cells, 'volume');
            [xc, ~, ~     ] = model.disc.transformCoords(x, cellNo);
            % Evaluate at cubature points
            p = model.disc.evaluateDGVariable(xc, cellNo, state, pdof);
        end
        
    end
    
end