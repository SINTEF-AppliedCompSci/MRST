classdef WaterThermalModel < ReservoirModel
    % Two phase oil/water system without dissolution
    properties

    end
    
    methods
        function model = WaterThermalModel(G, rock, fluid, varargin)
            
            model = model@ReservoirModel(G, rock, fluid);
            
            % This is the model parameters for oil/water
            model.oil = false;
            model.gas = false;
            model.water = true;
            %model.addflux=true;
            
            % Blackoil -> use CNV style convergence 
            model.useCNVConvergence = false;
            
            model.saturationVarNames = {'sw'};
            model.wellVarNames = {'qWs', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
            rock_heat=struct('perm',rock.lambdaR);
            T_r=computeTrans(G,rock_heat);
            cf = G.cells.faces(:,1);
            nf = G.faces.num;
            T_r  = 1 ./ accumarray(cf, 1./T_r, [nf, 1]);
            model.operators.T_r_all=T_r;
            intInx = all(G.faces.neighbors ~= 0, 2);
            model.operators.T_r = T_r(intInx);
            % Setup operators
            %model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function [vararg, driving] = getDrivingForces(model, control) %#ok
            % Setup and pass on driving forces
            vararg = {};
            driving = struct('Wells', [], 'bc', [], 'src', []);
            
            if isfield(control, 'W') && ~isempty(control.W)
                vararg = [vararg, 'Wells', control.W];
                driving.Wells = control.W;
            end
            
            if isfield(control, 'bc') && ~isempty(control.bc)
                vararg = [vararg, 'bc', control.bc];
                driving.bc = control.bc;
            end
            
            if isfield(control, 'src') && ~isempty(control.src)
                vararg = [vararg, 'src', control.src];
                driving.src = control.src;
            end
            
            if isfield(control, 'bcT') && ~isempty(control.bcT)
                vararg = [vararg, 'bcT', control.bcT];
                driving.bcT = control.bcT;
            end
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsWaterThermal(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Parent class handles almost everything for us
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            
            % Update wells based on black oil specific properties
            saturations = model.saturationVarNames;
            wi = strcmpi(saturations, 'sw');
            oi = strcmpi(saturations, 'so');
            gi = strcmpi(saturations, 'sg');

            W = drivingForces.Wells;
            state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);

        end
    end
end
