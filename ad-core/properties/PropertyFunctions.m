classdef PropertyFunctions
    properties
        
    end
    
    properties (Access = protected)
        structName
        structFields
        excludedFields
    end
    
    methods
        function props = PropertyFunctions()
            props.structFields = setdiff(properties(props), props.excludedFields);
        end
        
        function names = getPropertyNames(props)
            names = props.structFields;
        end
        
        function name = getPropertyContainerName(props)
            name = props.structName;
        end
        
        function [container, name] = getPropertyContainer(props, state)
            % Set up dynamic container (handle class) for storing
            % properties as we go
            name = props.getPropertyContainerName();
            if nargin > 1 && isfield(state, name)
                container = state.(name);
            else
                fld = [props.structFields(:)'; cell(1, numel(props.structFields))];
                s = struct(fld{:});
                container = HandleStruct(s);
            end
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
            props_struct.(name) = props.(name).evaluateOnDomain(model, state);
            if nargout > 0
                state.(struct_name) = props_struct;
            end
        end
        
        function v = get(props, model, state, name)
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
            if isfield(state, props.structName)
                % Cache is present, but this specific property is not
                % necessarily present
                ok = structPropEvaluated(state.(props.structName), name);
            else
                % Cache object is missing, we have no properties
                ok = false;
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
        
        function props = subset(props, cell_subset)
            names = props.structFields;
            for i = 1:numel(names)
                pn = names{i};
                if ~isempty(props.(pn))
                    props.(pn) = props.(pn).subset(cell_subset);
                end
            end
        end
    end
end