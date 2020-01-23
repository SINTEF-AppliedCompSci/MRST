classdef StateFunctionGrouping
    % A StateFunctionGrouping class is a grouping of interdependent properties.
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
        functionNames % Base name of all properties (i.e. what they implement)
        functionTypes % Indicator of property (0 for class member "intrinsic", 1 for stored)
        extraFunctions = {}; % Additional properties, not present as class properties
        excludedFields % Class properties which are not intended as functions
    end
    
    methods
        function props = StateFunctionGrouping(structname)
            if nargin > 0
                props.structName = structname;
            end
            props.functionNames = setdiff(properties(props), props.excludedFields);
            props.functionTypes = zeros(size(props.functionNames));
        end
        % ----------------------- Getters --------------------------------%
        function [names, types, implementation] = getNamesOfStateFunctions(props)
            % Get the names of all properties in collection. If a second
            % output argument is requested, it will give the internal
            % indicators if a property is intrinsic to the class (and will
            % always be present) or if it has been added (and cannot be
            % depended upon). If a third argument is requested, the class
            % name of the currently set implementation will be output as
            % well.
            names = props.functionNames;
            types = props.functionTypes;
            if nargout > 2
                n = numel(names);
                implementation = cell(n, 1);
                for i = 1:n
                    p = props.getStateFunction(names{i});
                    implementation{i} = class(p);
                end
            end
        end

        function prop = getStateFunction(props, name)
            % Get named property function (case-sensitive)
            sub = strcmp(props.functionNames, name);
            if props.functionTypes(sub) == 0
                prop = props.(name);
            else
                extrasub = sub(props.functionTypes == 1);
                prop = props.extraFunctions{extrasub};
            end
        end
        
        function name = getStateFunctionContainerName(props)
            % Get the name of the proprety container used to store
            % evaluated properties on state.
            name = props.structName;
        end

        function state = initStateFunctionContainer(group, state)
            % Set up state function container on a state
            [f, name] = group.getStateFunctionContainer();
            state.(name) = f;
        end

        function [container, name] = getStateFunctionContainer(props, state)
            % Set up dynamic container (handle class) for storing
            % properties as we go
            name = props.getStateFunctionContainerName();
            if nargin > 1 && isfield(state, name)
                container = state.(name);
            else
                fld = [props.functionNames(:)'; cell(1, numel(props.functionNames))];
                s = struct(fld{:});
                container = HandleStruct(s);
            end
        end
        % ----------------------- Setters --------------------------------%
        function props = setStateFunction(props, name, prop)
            % Set or replace a property function.
            % INPUT:
            %   - props (class instance, automatic)
            %   - name. Canonical name of property.
            %   - prop. Property implementation (StateFunction instance)
            assert(isa(prop, 'StateFunction'));
            assert(ischar(name));
            sub = strcmp(props.functionNames, name);
            if any(sub)
                % We are replacing an existing property
                ptypes = props.functionTypes;
                type = ptypes(sub);
                if type == 0
                    % Class property ("intrinsic"), just replace directly
                    props.(name) = prop;
                else
                    % This property was not found, we are adding extending
                    % the list of extra properties
                    assert(type == 1);
                    props.extraFunctions{sub(ptypes == 1)} = prop;
                end
            else
                % We are adding a new property and must extend the
                % corresponding internal datastructures.
                props.functionNames = [props.functionNames; name];
                props.functionTypes = [props.functionTypes; 1];
                props.extraFunctions{end+1, 1} = prop;
            end
        end

        function props = removeStateFunction(props, name)
            % Remove a property function from the list. Only allowed for
            % non-intrinsic functions.
            sub = strcmp(props.functionNames, name);
            if any(sub)
                ptypes = props.functionTypes;
                type = ptypes(sub);
                if type == 0
                    error('Cannot remove intrinsic property %s', name);
                else
                    assert(type == 1);
                    extrasub = sub(ptypes == 1);
                    props.extraFunctions = props.extraFunctions(~extrasub);
                    props.functionTypes = props.functionTypes(~sub);
                    props.functionNames = props.functionNames(~sub);
                end
            else
                warning('Property %s was not found, and cannot be removed.', name);
            end
        end
        % --------------- Evaluation functions ---------------------------%
        function [v, state] = get(props, model, state, name)
            % Get value of a property (possibily triggering several function
            % evaluations if required. Repeated calls to get within the
            % same AD-initialization of state for the same property will
            % not result in evaluations (caching)
            if ~props.isStateFunctionEvaluated(model, state, name)
                state = props.evaluateStateFunctionWithDependencies(model, state, name);
            end
            v = state.(props.structName).(name);
            v = expandIfUniform(v);
        end

        function state = evaluateStateFunction(props, model, state, name)
            % Force evaluation of a property, assuming all dependencies are
            % met. If all dependencies are not met, use
            % evaluateStateFunctionWithDependencies or simply get.
            struct_name = props.structName;
            if isstruct(state) && ~isfield(state, struct_name)
                props_struct = props.getStateFunctionContainer();
            else
                props_struct = state.(struct_name);
            end
            prop = props.getStateFunction(name);
            if isempty(prop.structName)
                prop.structName = props.structName;
            end
            props_struct.(name) = prop.evaluateOnDomain(model, state);
            if nargout > 0
                state.(struct_name) = props_struct;
            end
        end
        
        function state = evaluateStateFunctionWithDependencies(props, model, state, name)
            % Evaluate property, and all required dependencies in state.
            prop = props.getStateFunction(name);
            state = props.evaluateDependencies(model, state, prop.dependencies);
            state = props.evaluateStateFunction(model, state, name);
        end
        
        function ok = isStateFunctionEvaluated(props, model, state, name)
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
            % Internal function for evaluating a list of dependencies, in
            % an ordered fashion.
            % PARAMETERS:
            %   props - class instance
            %   model - model instance used to initialize the state
            %   state - state used to evaluate dependencies
            %   names - ordered cell array of all dependencies
            for i = 1:numel(names)
                name = names{i};
                if ~isStateFunctionEvaluated(props, model, state, name)
                    state = props.evaluateStateFunctionWithDependencies(model, state, name);
                end
            end
        end
        
        function props = subset(props, subset)
            % Take the subset of all the properties (reducing regions etc
            % to the new local domain defined by cell_subset).
            names = props.functionNames;
            for i = 1:numel(names)
                name = names{i};
                prop = props.getStateFunction(name);
                if ~isempty(prop)
                    prop = prop.subset(subset);
                end
                props = props.setStateFunction(name, prop);
            end
        end
        % ----- Display, plotting etc ------------------------------------%
        function disp(props)
            % Custom display function
            iname = inputname(1);
            fprintf('  ');
            name = class(props);
            isDesktop = usejava('desktop');
            if isDesktop
                fprintf(['<a href="matlab:helpPopup %s">%s</a> (<a href="matlab:edit %s.m">edit</a>', ...
                    '|<a href="matlab:figure;plotStateFunctionGroupings(%s)">plot</a>)'], name, name, name, iname);
            else
                fprintf('%s', class(props));
            end
            fprintf(' state function grouping instance.\n');
            names = props.functionNames;
            types = props.functionTypes;
            
            len = max(cellfun(@numel, names))+1;
            
            for type = 0:1
                typeIndices = find(types == type);
                if type == 0
                    if isempty(typeIndices)
                        fprintf('  Class has no intrinsic state functions.\n');
                    else
                        fprintf('\n  Intrinsic state functions (Class properties for %s, always present):\n', class(props));
                    end
                else
                    if isempty(typeIndices)
                        fprintf('  No extra state functions have been added to this class instance.\n');
                    else
                        fprintf('\n  Extra state functions (Added to class instance):\n');
                    end
                end
                for i = 1:numel(typeIndices)
                    index = typeIndices(i);
                    fn = props.getStateFunction(names{index});

                    cl = class(fn);
                    fprintf('    %*s: ', len, names{index});
                    if isDesktop
                        fprintf('<a href="matlab:helpPopup %s">%s</a>', cl, cl);
                        fprintf([' (<a href="matlab:edit %s.m">edit</a>|', ...
                            '<a href="matlab:figure;plotStateFunctionGroupings(%s,''center'',''%s'')">plot</a>)'],...
                            cl, iname, [class(props), '.', names{index}]);
                    else
                        fprintf('%s', class(fn));
                    end
                    fprintf('\n');
                end
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
