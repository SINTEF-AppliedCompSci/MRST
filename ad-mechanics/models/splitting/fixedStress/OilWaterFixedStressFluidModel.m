classdef OilWaterFixedStressFluidModel < OilWaterBiotModel
    
    properties
        pressCoef;
    end

    methods
        function model = OilWaterFixedStressFluidModel(G, rock, fluid, varargin)
            model = model@OilWaterBiotModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            
            opt = merge_options(opt, varargin{:});


            [p, sW, wellSol] = model.getProps(state, 'pressure', 'sw', 'wellsol');
            pBH    = vertcat(wellSol.bhp);
            qOs    = vertcat(wellSol.qOs);
            qWs    = vertcat(wellSol.qWs);

            [p0, sW0] = model.getProps(state0, 'pressure', 'sw');

            if ~opt.resOnly,
                [p, sW, qWs, qOs, pBH] = initVariablesADI(p, sW, qWs, qOs, pBH);
            end
            
            fnew = drivingForces.fixedStressTerms.new;
            mechTerm.new = fnew.pTerm.*p - fnew.sTerm;
            fold = drivingForces.fixedStressTerms.old;
            mechTerm.old = fold.pTerm.*p - fold.sTerm;
            
            otherDrivingForces = rmfield(drivingForces, 'fixedStressTerms');
            
            [eqs, state] = equationsOilWaterBiot(p0, sW0, p, sW, pBH, qWs, ...
                                                 qOs, state, model, dt, mechTerm, ...
                                                 otherDrivingForces, ...
                                                 varargin{:});

            primaryVars = {'pressure', 'sw', 'qWs', 'qOs', 'bhp'};
            names = {'water', 'oil' 'waterWells', 'oilWells', 'closureWells'};
            types = {'cell', 'cell', 'perf', 'perf', 'well'};
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@TwoPhaseOilWaterModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.fixedStressTerms = [];
        end
        
        
    end
    
end
    