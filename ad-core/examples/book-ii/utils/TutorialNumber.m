classdef TutorialNumber < StateFunction
    properties
        stateField;
    end
    
    methods
        function tn = TutorialNumber(n)
            tn.stateField = n;
            tn = tn.dependsOn(n, 'state');
        end
        function v = evaluateOnDomain(n, model, state)
            fprintf('Retrieving %s from state.\n', n.stateField);
            v = state.(n.stateField);
        end
    end
end