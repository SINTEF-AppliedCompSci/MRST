classdef ThreePhaseBlackOilSolventModel < ThreePhaseBlackOilModel
    
    properties
        solvent
        
    end
    
    methods
    
        function model = ThreePhaseBlackOilSolventModel(G, rock, fluid, ...
                                                                  varargin)
            
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);
            
            model.solvent = true;
            
            model.wellVarNames = {'qWs', 'qOs', 'qGs', 'qSs', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                                               dt, drivingForces, varargin)
                                           
            [problem, state] = equationsThreePhaseBlackOilSolvent(state0, state, ...
                                    model, dt, drivingForces, varargin{:});
            
        end
        
        function state = validateState(model, state)
            state = validateState@ThreePhaseBlackOilModel(mode, state, problem, dx, drivingForces);
            
            model.checkProperty(
        
    end
    
end
    