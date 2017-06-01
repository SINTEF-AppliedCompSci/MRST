classdef SinglephaseFixedStressFluidModel < WaterBiotModel
    
    properties
        pressCoef;
    end

    methods
        function model = SinglephaseFixedStressFluidModel(G, rock, fluid, varargin)
            model = model@WaterBiotModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            
            opt = merge_options(opt, varargin{:});

            [p, wellSol] = model.getProps(state, 'pressure', 'wellsol');
            pBH    = vertcat(wellSol.bhp);
            qWs    = vertcat(wellSol.qWs);
            p0 = model.getProp(state0, 'pressure');

            if ~opt.resOnly,
                [p, qWs, pBH] = initVariablesADI(p, qWs, pBH);
            end
            
            fnew = drivingForces.fixedStressTerms.new;
            mechTerm.new = fnew.pTerm.*p - fnew.sTerm;
            fold = drivingForces.fixedStressTerms.old;
            mechTerm.old = fold.pTerm.*p - fold.sTerm;
            
            otherDrivingForces = rmfield(drivingForces, 'fixedStressTerms');
            
            [eqs, state] = equationsWaterBiot(p0, p, qWs, pBH, state, model, dt, ...
                                              mechTerm, otherDrivingForces, ...
                                              varargin{:});
            primaryVars = {'pressure','qWs', 'bhp'};
            names = {'water'};
            types = {'cell'};
            W = drivingForces.W;
            if ~isempty(W)
                names(2:3) = {'waterWells', 'closureWells'};
                types(2:3) = {'perf', 'well'};
            end
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@WaterModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.fixedStressTerms = [];
        end
        
        
    end
    
end
    