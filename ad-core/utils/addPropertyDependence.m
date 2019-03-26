function prop = addPropertyDependence(prop, name, grouping)
    % Document dependencies and external dependencies
    if iscell(name)
        name = reshape(name, [], 1);
    end
    if nargin < 3 || isempty(grouping)
        if isstruct(name)
            prop.externals = [prop.externals; name];
        else
            prop.dependencies = [prop.dependencies; name];
        end
    else
        s = struct('name', name, 'grouping', grouping);
        prop.externals = [prop.externals; s];
    end
end