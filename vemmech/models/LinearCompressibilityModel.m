classdef LinearCompressibilityModel < ReservoirModel
    properties
        thermal = false;
    end
    
    methods
        function model = LinearCompressibilityModel(G, rock, fluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            % Only enable water
            model.oil = false;
            model.gas = false;
            model.water = true;
            model.useCNVConvergence = false;
            model.saturationVarNames = {'sw'};
            
            model = merge_options(model, varargin{:});
            
            model.stepFunctionIsLinear = true;
        end

        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, iteration, varargin)

        T = model.operators.T_all;
        
        solfn = @(A, b) linsolve.solveLinearSystem(A, b);
        state = lincompTPFA(dt, state0, model.G, T, model.operators.pv, model.fluid, model.rock,...
            'src', drivingForces.src,...
            'bc', drivingForces.bc, ...
            'Wells', drivingForces.W, ...
            'LinSolve', solfn,...
            'use_trans', true);
        
        if model.thermal
            assert(all(isfinite(model.rock.lambdaR)), ...
                ['Thermal conductivity is not finite and thermal is enabled.'...
                ' Aborting.']);
            rock_heat = struct('perm', model.rock.lambdaR);
            T_r = computeTrans(model.G, rock_heat);
            state = linearTransport(dt, state, model.G, T_r, model.operators.pv, model.fluid, model.rock,...
                            'src', drivingForces.src,...
                            'bc', drivingForces.bc, ...
                            'Wells', drivingForces.W ...
                            );
        end
            
        report = model.makeStepReport(...
                        'Failure',      false, ...
                        'Converged',    true, ...
                        'Residuals',    0, ...
                        'ResidualsConverged', true);
        end
    end
end