classdef ThresholdedTransmissibility < Transmissibility
    properties
        pressureThreshold
    end
    
    methods
        function pp = ThresholdedTransmissibility(backend, model)
            pp@Transmissibility(backend);
            equil = model.inputdata.REGIONS.EQLNUM(model.G.cells.indexMap);
            thpres = model.inputdata.SOLUTION.THPRES;
            N = model.operators.N;
            nf = size(N, 1);
            threshold = zeros(nf, 2);
            
            left_equil = equil(N(:, 1));
            right_equil = equil(N(:, 2));
            for i = 1:size(thpres, 1)
                left = left_equil == thpres(i, 1);
                right = right_equil == thpres(i, 2);
                
                threshold(left & right, 1) = thpres(i, 3);
                % Other direction
                left = left_equil == thpres(i, 2);
                right = right_equil == thpres(i, 1);
                threshold(left & right, 2) = thpres(i, 3);
            end
            pp.pressureThreshold = threshold;
        end
        
        function v = evaluateOnDomain(prop, model, state)
            v = model.operators.T;
            p = model.getProps(state, 'pressure');
            if isfield(model.fluid, 'transMult')
                v = model.fluid.transMult(p).*T;
            end
            pv = value(p);
            left = model.operators.N(:, 1);
            right = model.operators.N(:, 2);
            dp = pv(left) - pv(right);
            thres = prop.pressureThreshold;
            act = dp > thres(:, 1) | dp < thres(:, 2);
            v = v.*act;
        end
    end
end