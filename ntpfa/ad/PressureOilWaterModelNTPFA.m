classdef PressureOilWaterModelNTPFA < PressureOilWaterModel
    properties
        fluxType
        fixPointIteration
    end
    
    methods
        function model = PressureOilWaterModelNTPFA(G, rock, fluid, varargin)
            if isempty(fluid)
                fluid = initSimpleADIFluid();
            end
            
            model = model@PressureOilWaterModel(G, rock, fluid);
            model.incTolPressure = 1e-8;
            model.fixPointIteration = false;
            model = merge_options(model, varargin{:});
            model.fluxType = 'ntpfa';
            
            model.operators.collocationSet = getCollactionSet(G, model.rock);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = pressureEquationOilWaterNTPFA(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@PressureOilWaterModel(model, state, problem, dx, drivingForces);
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = checkConvergence@PressureOilWaterModel(model, problem);
        end
        
        function [kgrad, model] = getPermGradient(model, p, p0, forces, transMult)
        % Get gradient operator for K grad(p)
        switch lower(model.fluxType)
            case 'tpfa'
                T = model.operators.T.*transMult;
                C = model.operators.C;
                Ct = bsxfun(@times, C, T);                
                kgrad = @(x) -Ct*x;
            case 'mpfa'
                error('NotImplemented');
            case 'ntpfa'
                assert(isfield(model.operators, 'collocationSet'));
                if model.fixPointIteration
                    p = double(p);
                end
                nc = numel(double(p));
                N = model.operators.N;
                intx = model.operators.internalConn;

                nf = size(N, 1);
                T = computeNonLinearTrans(model.G, model.operators.collocationSet, p);
%                 T = bsxfun(@times, T, transMult);
                T1 = T{1}(intx);
                T2 = T{2}(intx);
                kgrad = @(x) T2.*x(N(:, 2)) - x(N(:,1)).*T1;
%                 C = sparse( [(1:nf)'; (1:nf)'], N, [-T(intx, 1), T(intx, 2)], nf, nc);
%                 kgrad = @(x) C*x;
            otherwise
                error(['Unknown discretization ''', model.fluxType, '''']);
        end
    end
    end

end