function [ variables ] = removeLogFromNames( variables )

variables = regexprep(variables, 'log', '');


end

