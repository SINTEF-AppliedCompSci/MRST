classdef TutorialProduct < StateFunction
    properties
        leftNumber
        rightNumber
    end
    methods
        function tp = TutorialProduct(left, right)
            [tp.leftNumber, tp.rightNumber] = deal(left, right);
            tp = tp.dependsOn({left, right});
        end
        function v = evaluateOnDomain(tp, model, state)
            [l, r] = tp.getEvaluatedDependencies(state, tp.leftNumber, tp.rightNumber);
            fprintf('Multiplying %s and %s.\n', tp.leftNumber, tp.rightNumber);
            v = l*r;
        end
    end
end