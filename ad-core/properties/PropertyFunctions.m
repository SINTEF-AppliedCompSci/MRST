classdef PropertyFunctions
    properties
        
    end
    
    properties (Access = protected)
        structName
        structFields
    end
    
    methods
        function [container, name] = getPropertyContainer(props)
            % Set up dynamic container (handle class) for storing
            % properties as we go
            name = props.structName;
            fld = [props.structFields(:)'; cell(1, numel(props.structFields))];
            s = struct(fld{:});
            container = DynamicStruct(s);
        end
        
        function evaluateProperty(props, model, state, name)
            % Force evaluation of a property, assuming all dependencies are
            % met
            props_struct = state.(props.structName);
            props_struct.(name) = props.(name).evaluateOnGrid(model, state);
        end
        
        function v = getProperty(props, model, state, name)
            % Get a property with lazy evaluation
            props.evaluateDependencies(model, state, {name})
            v = state.(props.structName).(name);
        end
        
        function ok = isPropertyEvaluated(props, model, state, name)
            % Check if property is present in cache
            ok = ~isempty(state.(props.structName).(name));
        end
        
        function evaluateDependencies(props, model, state, names)
            % Evaluate dependencies (order dependent)
            for i = 1:numel(names)
                name = names{i};
                if ~isPropertyEvaluated(props, model, state, name)
                    props.evaluateProperty(model, state, name);
                end
            end
        end
    end
end