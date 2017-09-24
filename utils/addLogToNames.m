function [ variables ] = addLogToNames( variables )

    for i = 1 : numel(variables)
        variables{i} = ['log', variables{i}];
    end

    ind = strcmpi('logfluidVolumeFraction',variables);
    variables(ind) = {'logFluidVolumeFraction'};
    
    ind = strcmpi('logCVC',variables);
    variables(ind) = {'CVC'};
    
    ind = strcmpi('logCVC',variables);
    variables(ind) = {'CVC'};
    
end

