classdef WaterFixedStressFluidModel < WaterModel
    

    methods
        function model = WaterFixedStressFluidModel(G, rock, fluid, varargin)
            model = model@WaterModel(G, rock, fluid);
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

            [wellVars, wellVarNames, wellMap] = ...
                model.FacilityModel.getAllPrimaryVariables(wellSol);

            if ~opt.resOnly,
                [p, wellVars{:}] = initVariablesADI(p, wellVars{:});
            end
            
            fnew = drivingForces.fixedStressTerms.new;
            mechTerm.new = fnew.pTerm.*p - fnew.sTerm;

            fold = drivingForces.fixedStressTerms.old;
            mechTerm.old = fold.pTerm.*p - fold.sTerm;
            
            otherDrivingForces = rmfield(drivingForces, 'fixedStressTerms');
            
            [eqs, names, types, state] = equationsWaterMech(state0, p, ...
                                                            wellVars, state, ...
                                                            model, dt, ...
                                                            mechTerm, ...
                                                            otherDrivingForces, ...
                                                            'iteration', ...
                                                            opt.iteration);
            
            primaryVars = {'pressure', wellVarNames{:}};
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@WaterModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.fixedStressTerms = [];
        end
        
        function fds = getAllVarsNames(model)
            fds = {'wellSol', 'pressure'};
        end
    
    end
    
end
    