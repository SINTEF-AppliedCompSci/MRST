function var = expand_var(var, target_size)
% Ensure that 'var' is of the right dimension and size, as specified by 'target_size'

    if isstr(var)
        var = {var}; % ensure strings are treated with cell arrays below
    end
    
    if isscalar(var) && (prod(target_size) > 1)
        if iscell(var)
            tmp = var;
            var = cell(target_size,1);
            for i = 1:numel(var)
                var{i} = tmp;
            end
        else 
            var = var * ones([target_size, 1]); % ',1' avoids the creation of a
                                                % square matrix in the scalar
                                                % case
        end
    end

    % Sanity check - will catch if the size of res doesn't match the
    % prescribed size after all.
    if iscell(var)
        assert(numel(var) == target_size);
    else
        varsize = size(var);
        % ignoring trivial dimensions
        varsize = varsize(varsize~=1);
        target_size = target_size(target_size~=1);
        % Verifying that the numbers of nontrivial dimensions are the same
        assert(numel(varsize) == numel(target_size));
        % Verify that the sizes of nontrivial dimensions are the same
        assert(isequal(varsize, target_size));
    end
end


