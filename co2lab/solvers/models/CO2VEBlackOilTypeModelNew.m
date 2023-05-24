classdef CO2VEBlackOilTypeModelNew < ReservoirModel & GenericReservoirModel
    
    properties
    end

    methods
        function model = CO2VEBlackOilTypeModelNew(Gt, rock2D, fluid, varargin)
            
        end
        
        function [eqs, names ,types, state] = ...
                getModelEquations(model, state0, state, dt, drivingForces)
        end
        
        function model = setupOperators(model, G, rock, varargin)
        end
        
        function model = validateModel(model, varargin)
        end
        
        function state = validateState(model, varargin)
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
        end
        
        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
        end
        
        function [state, report] = ...
                updateAfterConvergence(model, state0, state, dt, drivingForces)
        end

        function [fn, index] = getVariableField(model, name, varargin)
        end
        
        function names = getComponentNames(model)
        end
    
        function [vars, names, origin] = getPrimaryVariables(model, state)
        end
        
        function [eqs, state, src] = ...
                addBoundaryConditionsAndSources(model, eqs, names, types, ...
                                                state, forces)
        end
        
        function ctrl = validateDrivingForces(model, ctrl, index)
        end
        
        function forces = getValidDrivingForces(model)
        end
        
        function state = initStateAD(model, state, vars, names, origin)
        end
            
        function [v_eqs, tolerances, names] = ...
                getConvergenceValues(model, problem, varargin)
        end
        
        function [model, state] = updateForChangedControls(model, state, forces)
        end

        function [eq, src] = ...
                addComponentContributions(model, cname, eq, component, src, force)        
        end
        
    end
end



---------------------------------------

function model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, varargin)
function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
function [fn, index] = getVariableField(model, name, varargin)
function [state, report] = updateState(model, state, problem, dx, drivingForces)
function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
function model = setupOperators(model, Gt, rock, varargin)
    
    
function gdz = getGravityGradient(model)
function rhoS = getSurfaceDensities(model)
function [problem, state] = getAdjointEquations(model, state0, state, dt, drivingForces, varargin)
