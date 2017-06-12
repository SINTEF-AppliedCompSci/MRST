classdef OilWaterFixedStressFluidModel < TwoPhaseOilWaterModel
    
    methods
        function model = OilWaterFixedStressFluidModel(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
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

            [wellVars, wellVarNames, wellMap] = ...
                model.FacilityModel.getAllPrimaryVariables(wellSol);

            if ~opt.resOnly,
                [p, sW, wellVars{:}] = initVariablesADI(p, sW, wellVars{:});
            end
            
            fnew = drivingForces.fixedStressTerms.new;
            mechTerm.new = fnew.pTerm.*p - fnew.sTerm;
            fold = drivingForces.fixedStressTerms.old;
            mechTerm.old = fold.pTerm.*p - fold.sTerm;
            
            otherDrivingForces = rmfield(drivingForces, 'fixedStressTerms');
            
            [eqs, names, types, state] = equationsOilWaterMech(state0, p, sW, ...
                                                              wellVars, state, ...
                                                              model, dt, ...
                                                              mechTerm, ...
                                                              otherDrivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            primaryVars = {'pressure', 'sw', wellVarNames{:}};
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@TwoPhaseOilWaterModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.fixedStressTerms = [];
        end
        
        function fds = getAllVarsNames(model)
            fds = {'wellSol', 'pressure', 's'};
        end
        
    end
    
end
    