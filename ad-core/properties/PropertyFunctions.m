classdef PropertyFunctions
    properties
        
    end
    
    properties (Access = protected)
        structName
        structFields
    end
    
    methods
        function props = PropertyFunctions()
            props.structFields = properties(props);
        end
        
        function [container, name] = getPropertyContainer(props)
            % Set up dynamic container (handle class) for storing
            % properties as we go
            name = props.structName;
            fld = [props.structFields(:)'; cell(1, numel(props.structFields))];
            s = struct(fld{:});
            container = HandleStruct(s);
        end
        
        function state = evaluateProperty(props, model, state, name)
            % Force evaluation of a property, assuming all dependencies are
            % met
            struct_name = props.structName;
            if isstruct(state) && ~isfield(state, struct_name)
                props_struct = props.getPropertyContainer();
            else
                props_struct = state.(struct_name);
            end
            props_struct.(name) = props.(name).evaluateOnGrid(model, state);
            if nargout > 0
                state.(struct_name) = props_struct;
            end
        end
        
        function v = getProperty(props, model, state, name)
            % Get a property with lazy evaluation
            state = props.evaluateDependencies(model, state, {name});
            v = state.(props.structName).(name);

            if isnumeric(v) && size(v, 2) > 1
                n = size(v, 2);
                out = cell(1, n);
                for i = 1:n
                    out{i} = v(:, i);
                end
                v = out;
            end
        end
        
        function ok = isPropertyEvaluated(props, model, state, name)
            % Check if property is present in cache
            if isstruct(state) && ~isfield(state, props.structName)
                % Cache object is missing, we have no properties
                ok = false;
            else
                % Cache is present, but this specific property is not
                % necessarily present
                ok = ~isempty(state.(props.structName).(name));
            end
        end
        
        function state = evaluateDependencies(props, model, state, names)
            % Evaluate dependencies (order dependent)
            for i = 1:numel(names)
                name = names{i};
                if ~isPropertyEvaluated(props, model, state, name)
                    state = props.evaluateProperty(model, state, name);
                end
            end
        end
    end
end