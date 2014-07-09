classdef TwoPhaseOilWaterModel < PhysicalModel
    % two phase oil / water system
    properties

    end
    
    methods
        function model = TwoPhaseOilWaterModel(G, rock, fluid, varargin)
            
            model = model@PhysicalModel(G, rock, fluid);
            
            % This is the model parameters for oil/water
            model.oil = true;
            model.gas = false;
            model.water = true;
            
            model.phaseNames = {'water', 'oil'};
            
            model = merge_options(model, varargin{:});
            
            % Name and operators
            model.name = 'OilWater_2ph';
            model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            dsMax = model.dsMax;
            dpMax = model.dpMax;
            
            pInd = problem.indexOfPrimaryVariable('pressure');
            wInd = problem.indexOfPrimaryVariable('sW');
            
            dp = dx{pInd};
            ds = dx{wInd};

            ds = sign(ds).*min(abs(ds), dsMax);
            dp = sign(dp).*min(abs(dp), abs(dpMax.*state.pressure));

            state.pressure = state.pressure + dp;
            sw = state.s(:,1) + ds;
            % Cap values
            sw = min(sw, 1); sw = max(sw, 0);

            state.s = [sw, 1-sw];
            
            dqWs    = dx{problem.indexOfPrimaryVariable('qWs')};
            dqOs    = dx{problem.indexOfPrimaryVariable('qOs')};
            dpBHP   = dx{problem.indexOfPrimaryVariable('bhp')};
            
            if ~isempty(dpBHP)
                dpBHP = sign(dpBHP).*min(abs(dpBHP), abs(dpMax.*vertcat(state.wellSol.bhp)));
                for w = 1:numel(state.wellSol)
                    state.wellSol(w).bhp      = state.wellSol(w).bhp + dpBHP(w);
                    state.wellSol(w).qWs      = state.wellSol(w).qWs + dqWs(w);
                    state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
                end
            end
            
            report = struct();
        end
        
    end
end
