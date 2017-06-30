function [ states ] = changeUnits( states, unit_conv )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

unit_conv = unit_conv^-1;

if numel(states) == 1
        states.components = states.components*unit_conv;
        states.masterComponents = states.masterComponents*unit_conv;
        if isfield(states, 'activities')
           states.activities = states.activities*unit_conv;
        end
        
        states.logcomponents = log(states.components);
        states.logmasterComponents = log(states.masterComponents);
else

    for i = 1 : numel(states)
    
        states{i}.components = states{i}.components*unit_conv;
        states{i}.masterComponents = states{i}.masterComponents*unit_conv;
        if isfield(states, 'activities')
           states.activities{i} = states.activities{i}*unit_conv;
        end
        states{i}.logcomponents = log(states{i}.components);
        states{i}.logmasterComponents = log(states{i}.masterComponents);
        
    end
end

end

