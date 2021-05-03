classdef StateFunctionGrouping < StateFunctionDependent
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
        validationLevel = 0; % Level of validation checking performed. Performance hit > 0!
    end
    
    methods
        function group = StateFunctionGrouping(structname)
            if nargin > 0
                group.structName = structname;
            end
            % Check if running Octave. Octave does not support querying
            % the properties of a class under construction. For this reason
            % we perform a check and hope for the best.
            isOctave = mrstPlatform('octave');
            if ~isOctave
                group = group.setInternalNames();
            end
        end

        % ----------------------- Getters --------------------------------%
        function [names, types, implementation, labels] = getNamesOfStateFunctions(props)
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
                labels = cell(n, 1);
                for i = 1:n
                    p = props.getStateFunction(names{i});
                    implementation{i} = class(p);
                    s = p.label;
                    if isempty(s)
                        s = names{i};
                    end
                    labels{i} = s;
                end
            end
        end

        function prop = getStateFunction(props, name)
            % Get named property function (case-sensitive)
            [present, sub] = props.hasStateFunction(name);
            if ~present
                error('%s is not known to this instance of the %s StateFunctionGrouping', name, class(props));
            end
            if props.functionTypes(sub) == 0
                prop = props.(name);
            else
                extrasub = sub(props.functionTypes == 1);
                prop = props.extraFunctions{extrasub};
            end
        end
        
        function [present, sub] = hasStateFunction(group, name)
            sub = strcmp(group.functionNames, name);
            present = any(sub);
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

        function present = isStateInitialized(group, state)
            % Check if the state contains the requisite container
            present = isfield(state, group.structName);
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
        function group = setStateFunction(group, name, prop)
            % Set or replace a property function.
            % INPUT:
            %   - props (class instance, automatic)
            %   - name. Canonical name of property.
            %   - prop. Property implementation (StateFunction instance)
            assert(isa(prop, 'StateFunction'));
            assert(ischar(name));
            group = group.setInternalNames();
            sub = strcmp(group.functionNames, name);
            prop.structName = group.structName;
            if any(sub)
                % We are replacing an existing property
                ptypes = group.functionTypes;
                type = ptypes(sub);
                if type == 0
                    % Class property ("intrinsic"), just replace directly
                    % group.(name) = prop;
                    sr = struct('subs', name, 'type', '.');
                    group = builtin('subsasgn', group, sr, prop);
                else
                    % This property was not found, we are adding extending
                    % the list of extra properties
                    assert(type == 1);
                    group.extraFunctions{sub(ptypes == 1)} = prop;
                end
            else
                % We are adding a new property and must extend the
                % corresponding internal datastructures.
                group.functionNames = [group.functionNames; name];
                group.functionTypes = [group.functionTypes; 1];
                group.extraFunctions{end+1, 1} = prop;
            end
        end
        
        function group = subsasgn(group, sub, val)
            if numel(sub) == 1 && strcmp(sub.type, '.')
                group = group.setStateFunction(sub.subs, val);
            else
                group = builtin('subsasgn', group, sub, val);
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
        function [v, state] = get(props, model, state, name, expand_to_cell)
            % Get value of a property (possibily triggering several function
            % evaluations if required. Repeated calls to get within the
            % same AD-initialization of state for the same property will
            % not result in evaluations (caching)
            if nargin == 4
                expand_to_cell = false;
            end
            if ~props.isStateFunctionEvaluated(model, state, name)
                state = props.evaluateStateFunctionWithDependencies(model, state, name);
            end
            v = state.(props.structName).(name);
            if expand_to_cell
                % Always return a cell array.
                v = expandMatrixToCell(v);
            else
                % Return a cell array if output was a matrix only.
                v = expandIfUniform(v);
            end
        end

        function state = evaluateStateFunction(sfg, model, state, name)
            % Force evaluation of a property, assuming all dependencies are
            % met. If all dependencies are not met, use
            % evaluateStateFunctionWithDependencies or simply get.
            struct_name = sfg.structName;
            if isstruct(state) && ~isfield(state, struct_name)
                props_struct = sfg.getStateFunctionContainer();
            else
                props_struct = state.(struct_name);
            end
            prop = sfg.getStateFunction(name);
            props_struct.(name) = prop.evaluateOnDomain(model, state);
            if sfg.validationLevel > 0
                prop.validateOutput(props_struct.(name), sfg.validationLevel);
            end
            if nargout > 0
                state.(struct_name) = props_struct;
            end
        end
        
        function state = evaluateStateFunctionUnsafe(sfg, model, state, name)
            % Evaluate if missing, assuming that:
            % - Struct is here.
            % - No validation should be performed.
            % - May re-evaluate the state function
            prop = sfg.getStateFunction(name);
            props_struct = state.(sfg.structName);
            props_struct.(name) = prop.evaluateOnDomain(model, state);
        end

        
        function state = evaluateStateFunctionWithDependencies(props, model, state, name)
            % Evaluate property, and all required dependencies in state.
            prop = props.getStateFunction(name);
            state = props.evaluateDependencies(model, state, prop.dependencies);
            state = props.evaluateExternalDependencies(model, state, prop.externals);
            state = props.evaluateStateFunction(model, state, name);
        end
        
        function ok = isStateFunctionEvaluated(props, model, state, dep)
            % Check if property is present in cache.
            if isfield(state, 'evaluated')
                ok = true;
            elseif ischar(dep)
                % Internal dependency - same group
                nm = props.structName;
                if isfield(state, nm)
                    % Cache is present, but this specific property is not
                    % necessarily present
                    if props.validationLevel && ~isfield(state.(nm), dep)
                        error(['Did not find %s in %s field. %s does not appear', ...
                            ' to belong to %s. Check your dependencies.'], ...
                            dep, nm, dep, class(props));
                    end
                    ok = structPropEvaluated(state.(nm), dep);
                else
                    % Cache object is missing, we have no properties
                    ok = false;
                end
            else
                % External dependency - either state or some other function
                % group belonging to the model
                if strcmp(dep.grouping, 'state')
                    ok = true;
                else
                    try
                        ok = model.(dep.grouping).isStateFunctionEvaluated(model, state, dep.name);
                    catch
                        ok = true;
                    end
                end
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
        
        function state = evaluateExternalDependencies(props, model, state, externals)
            % Internal function for evaluating a list external of dependencies, in
            % an ordered fashion.
            % PARAMETERS:
            %   props - class instance
            %   model - model instance used to initialize the state
            %   state - state used to evaluate dependencies
            %   names - ordered struct array of all external dependencies
            for i = 1:numel(externals)
                dep = externals(i);
                if strcmpi(dep.grouping, 'state')
                    % Do nothing
                    continue
                end
                if ~isStateFunctionEvaluated(props, model, state, dep)
                    state = model.(dep.grouping).evaluateStateFunctionWithDependencies(model, state, dep.name);
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
        
        function ok = checkDependencies(group, groups)
            % Checks if dependencies are met for a given group. Internal
            % dependencies are always checked, and externals are checked if
            % a cell array of other groups are provided.
            %
            % If output is requested, this function returns a boolean. If
            % not, it produces an error at first unmet dependency.
            % Additional diagnostic can be enabled by mrstVerbose.
            checkExternal = nargin > 1;
            throwError = nargout == 0;
            if throwError
                dispfn = @(varargin) error(varargin{:});
            else
                dispfn = @(varargin) dispif(mrstVerbose(), varargin{:});
            end
            
            nf = numel(group.functionNames);
            ok = true;
            for i = 1:nf
                % Check all dependencies
                name = group.functionNames{i};
                fn = group.getStateFunction(name);
                % We first check internals (i.e. in same group)
                for j = 1:numel(fn.dependencies)
                    depname = fn.dependencies{j};
                    present = group.hasStateFunction(depname);
                    if ~present
                        dispfn('Unmet internal dependency for %s in group %s: Did not find %s in own group.\n', name, class(group), depname);
                    end
                    ok = ok && present;
                end
                % If groups were provided, we also check the externals for
                % our group
                if checkExternal
                    for j = 1:numel(fn.externals)
                        e = fn.externals(j);
                        extgroupname = e.grouping;
                        if strcmpi(extgroupname, 'state')
                            % State is always fulfilled
                            continue
                        end
                        depname = e.name;
                        candidates = find(cellfun(@(x) isa(x, extgroupname), groups));
                        if isempty(candidates)
                            ok = false;
                            dispfn('Unmet external dependency for %s in group %s: Did not find any groups of class %s that may contain %s.\n',...
                                    name, class(group), extgroupname, depname);
                        else
                            for k = 1:numel(candidates)
                                found = groups{candidates(k)}.hasStateFunction(depname);
                                if found
                                    break;
                                end
                            end
                            ok = ok && found;
                            if ~found
                                dispfn('Unmet external dependency for function %s in group %s: Did not find %s in group of class %s.\n',...
                                    name, class(group), depname, extgroupname);
                            end
                        end
                    end
                end
            end
            if ok
                if checkExternal
                    dparg = 'internal and external';
                else
                    dparg = 'internal';
                end
                dispif(mrstVerbose(), 'All %s dependencies are met for group of class %s.\n', dparg, class(group));
            end
        end
        
        % ----- Display, plotting etc ------------------------------------%
        function disp(props)
            % Custom display function
            iname = inputname(1);
            canPlot = ~isempty(iname); % Make sure that the variable is defined in workspace
            fprintf('  ');
            name = class(props);
            allowRichText = mrstPlatform('richtext');
            if allowRichText
                fprintf('<a href="matlab:helpPopup %s">%s</a> (<a href="matlab:edit %s.m">edit</a>', name, name, name);
                if canPlot
                    fprintf('|<a href="matlab:figure;plotStateFunctionGroupings(%s)">plot</a>', iname);
                end
                fprintf(')');
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
                    if allowRichText
                        fprintf('<a href="matlab:helpPopup %s">%s</a>', cl, cl);
                        fprintf(' (<a href="matlab:edit %s.m">edit</a>', cl);
                        if canPlot
                            fprintf(['|<a href="matlab:figure;plotStateFunctionGroupings(%s,', ...
                                     '''center'',''%s'')">plot</a>'],...
                                iname, [class(props), '.', names{index}]);
                        end
                        fprintf(')');
                    else
                        fprintf('%s', class(fn));
                    end
                    fprintf('\n');
                end
            end
        end
        
        function props = setCompactEvaluation(props, val)
            % Use more compact evaluation for properties where the same
            % value can be evaluated multiple times (e.g. multiple
            % densities in the single-phase region)
            names = props.functionNames;
            for i = 1:numel(names)
                name = names{i};
                fn = props.getStateFunction(name);
                if isprop(fn, 'useCompactEvaluation')
                    fn.useCompactEvaluation = val;
                    props = props.setStateFunction(name, fn);
                end
            end
        end
        
        function group = setValidationLevel(group, level)
            group.validationLevel = level;
        end
    end
    methods (Access = protected)
        function group = setInternalNames(group)
            if ~iscell(group.functionNames)
                group.functionNames = setdiff(propertynames(group), group.excludedFields);
                group.functionTypes = zeros(size(group.functionNames));
            end
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
