function [ variables ] = addLogToNames( variables )

    for i = 1 : numel(variables)
        variables{i} = ['log', variables{i}];
    end

    
    ind = strcmpi('logCVC',variables);
    variables(ind) = {'CVC'};
    
    
end

