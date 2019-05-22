classdef AutoDiffFunctionGrouping
    % A AutoDiffFunctions class is a grouping of interdependent properties.
    % Inside each grouping, dependencies are handled and all grouped
    % properties are stored in the same normalized field on the state.
    %
    % Each grouping has intrinsic propreties which must be defined and
    % typically represent external interfaces. These are given as
    % properties on the inherited subclass, with name corresponding to the
    % name of the property.
    %
    % In addition, property functions can be added to the class instance.
    % These functions represent additional properties which can be
    % dependencies of the intrinsics properties.
    properties
        
    end
    
    properties (Access = protected)
        structName % Name of the struct where the properties should be stored on state
        propertyNames % Base name of all properties (i.e. what they implement)
        propertyTypes % Indicator of property (0 for class member "intrinsic", 1 for stored)
        extraProperties = {}; % Additional properties, not present as class properties
        excludedFields % Class properties which are not intended as functions
    end
    
    methods
        function props = AutoDiffFunctionGrouping()
            props.propertyNames = setdiff(properties(props), props.excludedFields);
            props.propertyTypes = zeros(size(props.propertyNames));
        end
        % ----------------------- Getters --------------------------------%
        function [names, types] = getPropertyNames(props)
            % Get the names of all properties in collection. If a second
            % output argument is requested, it will give the internal
            % indicators if a property is intrinsic to the class (and will
            % always be present) or if it has been added (and cannot be
            % depended upon).
            names = props.propertyNames;
            types = props.propertyTypes;
        end

        function prop = getPropertyFunction(props, name)
            % Get named property function (case-sensitive)
            sub = strcmp(props.propertyNames, name);
            if props.propertyTypes(sub) == 0
                prop = props.(name);
            else
                extrasub = sub(props.propertyTypes == 1);
                prop = props.extraProperties{extrasub};
            end
        end
        
        function name = getPropertyContainerName(props)
            % Get the name of the proprety container used to store
            % evaluated properties on state.
            name = props.structName;
        end

        function [container, name] = getPropertyContainer(props, state)
            % Set up dynamic container (handle class) for storing
            % properties as we go
            name = props.getPropertyContainerName();
            if nargin > 1 && isfield(state, name)
                container = state.(name);
            else
                fld = [props.propertyNames(:)'; cell(1, numel(props.propertyNames))];
                s = struct(fld{:});
                container = HandleStruct(s);
            end
        end
        % ----------------------- Setters --------------------------------%
        function props = setPropertyFunction(props, prop, name)
            % Set or replace a property function
            sub = strcmp(props.propertyNames, name);
            if any(sub)
                % We are replacing an existing property
                ptypes = props.propertyTypes;
                type = ptypes(sub);
                if type == 0
                    % Class property ("intrinsic"), just replace directly
                    props.(name) = prop;
                else
                    % This property was not found, we are adding extending
                    % the list of extra properties
                    assert(type == 1);
                    props.extraProperties{sub(ptypes == 1)} = prop;
                end
            else
                % We are adding a new property and must extend the
                % corresponding internal datastructures.
                props.propertyNames = [props.propertyNames; name];
                props.propertyTypes = [props.propertyTypes; 1];
                props.extraProperties{end+1, 1} = prop;
            end
        end

        function props = removePropertyFunction(props, name)
            % Remove a property function from the list. Only allowed for
            % non-intrinsic functions.
            sub = strcmp(props.propertyNames, name);
            if any(sub)
                ptypes = props.propertyTypes;
                type = ptypes(sub);
                if type == 0
                    error('Cannot remove intrinsic property %s', name);
                else
                    assert(type == 1);
                    extrasub = sub(ptypes == 1);
                    props.extraProperties = props.extraProperties(~extrasub);
                    props.propertyTypes = props.propertyTypes(~sub);
                    props.propertyNames = props.propertyNames(~sub);
                end
            else
                warning('Property %s was not found, and cannot be removed.', name);
            end
        end
        % --------------- Evaluation functions ---------------------------%
        function v = get(props, model, state, name)
            % Get value of a property (possibily triggering several function
            % evaluations if required.
            if ~props.isPropertyEvaluated(model, state, name)
                state = props.evaluatePropertyWithDependencies(model, state, name);
            end
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

        function state = evaluateProperty(props, model, state, name)
            % Force evaluation of a property, assuming all dependencies are
            % met. If all dependencies are not met, use
            % evaluatePropertyWithDependencies or simply get.
            struct_name = props.structName;
            if isstruct(state) && ~isfield(state, struct_name)
                props_struct = props.getPropertyContainer();
            else
                props_struct = state.(struct_name);
            end
            prop = props.getPropertyFunction(name);
            if isempty(prop.structName)
                prop.structName = props.structName;
            end
            props_struct.(name) = prop.evaluateOnDomain(model, state);
            if nargout > 0
                state.(struct_name) = props_struct;
            end
        end
        
        function state = evaluatePropertyWithDependencies(props, model, state, name)
            % Evaluate property, and all required dependencies in state.
            prop = props.getPropertyFunction(name);
            state = props.evaluateDependencies(model, state, prop.dependencies);
            state = props.evaluateProperty(model, state, name);
        end
        
        function ok = isPropertyEvaluated(props, model, state, name)
            % Check if property is present in cache.
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
                    state = props.evaluatePropertyWithDependencies(model, state, name);
                end
            end
        end
        
        function props = subset(props, cell_subset)
            % Take the subset of all the properties (reducing regions etc
            % to the new local domain defined by cell_subset).
            names = props.propertyNames;
            for i = 1:numel(names)
                name = names{i};
                prop = props.getPropertyFunction(name);
                if ~isempty(prop)
                    prop = prop.subset(cell_subset);
                end
                props.setPropertyFunction(prop, name);
            end
        end
        
        function disp(props)
            % Custom display function
            fprintf('  ');
            name = class(props);
            isDesktop = usejava('desktop');
            if isDesktop
                fprintf('<a href="matlab:helpPopup %s">%s</a> (<a href="matlab:edit %s.m">edit</a>)', name, name, name);
            else
                fprintf('%s', class(props));
            end
            fprintf(' property function grouping instance.\n');
            names = props.propertyNames;
            types = props.propertyTypes;
            
            len = max(cellfun(@numel, names))+1;
            
            for type = 0:1
                if type == 0
                    fprintf('\n  Intrinsic properties (Class properties for %s):\n\n', class(props));
                else
                    fprintf('\n  Extra properties (Added to class instance):\n\n');
                end
                typeIndices = find(types == type);
                for i = 1:numel(typeIndices)
                    index = typeIndices(i);
                    fn = props.getPropertyFunction(names{index});

                    cl = class(fn);
                    fprintf('%*s: ', len, names{index});
                    if isDesktop
                        fprintf('<a href="matlab:helpPopup %s">%s</a>', cl, cl);
                        fprintf(' (<a href="matlab:edit %s.m">edit</a>)', cl);
                    else
                        fprintf('%s', class(fn));
                    end
                    fprintf('\n');
                end
            end
        end
    end
end
