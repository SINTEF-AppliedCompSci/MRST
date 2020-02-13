classdef NumericalPressureReductionFactors < StateFunction
   
    properties
    end
    
    methods
        function prf = NumericalPressureReductionFactors(model, varargin)
            
            prf@StateFunction(model, varargin{:});
            prf.dependsOn({'Density'});
            
        end
        
        function weights = evaluateOnDomain(prop, model, state)
           
            state = model.storeDensities(state, rhoW, rhoO, rhoG);
            
            rho = prop.getEvaluatedDependencies(state, 'Densities');
            
            ncell = model.G.cells.num;
            ndof = ncell*(ncomp+wat);
            wellVarIndices = (ndof+1):(ndof+nwellvar);

            [w, dwdp] = getPartialVolumes(model, state, acc, ...
                'iteration',                  opt.iteration, ...
                'wellVarIndices',             wellVarIndices, ...
                'singlePhaseStrategy',        model.singlePhaseStrategy, ...
                'twoPhaseStrategy',           model.twoPhaseStrategy, ...
                'singlePhaseDifferentiation', model.singlePhaseDifferentiation, ...
                'twoPhaseDifferentiation',    model.twoPhaseDifferentiation);

            weights = cell(ncomp+wat, 1);
            for i = 1:(ncomp+wat)
                wi = w(:, i);
                if any(dwdp)
                    Wp = double2ADI(wi, p);
                    Wp.jac{1} = sparse(1:ncell, 1:ncell, dwdp(:, i), ncell, ncell);
                else
                    Wp = wi;
                end
                weights{i} = Wp;
            end
            
        end

    end
    
end