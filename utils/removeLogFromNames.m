function [ variables ] = removeLogFromNames( variables )

variables = regexprep(variables, 'log', '');


ind = strcmpi('fluidVolumeFraction',variables);
variables{ind} = 'fluidVolumeFraction';
end

