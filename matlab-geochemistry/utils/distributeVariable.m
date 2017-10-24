function [ variableValues ] = distributeVariable( variableNames, knownNames, unknownNames, knownVal, unknownVal )

    variableValues = cell(1, numel(variableNames));        
    for i = 1 : numel(variableNames)
        ind = strcmpi(unknownNames, variableNames{i});
        if any(ind)
            variableValues{i} = unknownVal{ind};
        end
        ind = strcmpi(knownNames, variableNames{i});
        if any(ind)
            variableValues{i} = knownVal{ind};
        end
    end

end

