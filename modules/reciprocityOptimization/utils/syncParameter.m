function setup = syncParameter(setup, params, param_name, new_value)
    % Update a single parameter in `setup` (mimicking setParameter logic)
    idx = find(strcmp({params.name}, param_name));
    if ~isempty(idx)
        params{idx}.value = new_value;  % Update the parameter struct
        setup = params{idx}.setParameter(setup, new_value);
    end
end