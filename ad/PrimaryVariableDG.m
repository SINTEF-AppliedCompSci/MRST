classdef PrimaryVariableDG < PrimaryVariable
   
    properties
        parentProp
    end
    
    methods
        function prop = PrimaryVariableDG(parentProp)

            model = PhysicalModel([], 'AutoDiffBackend', parentProp.AutoDiffBackend);
            prop = prop@PrimaryVariable(model);
            prop.dofname = [parentProp.dofname, 'dof'];
            prop.parentProp = parentProp;
        end
        
        function value = evaluateOnDomain(prop, model, state)
            % Given state, evaluate the canonical representation for the
            % current model.
            dof = state.(prop.dofname);
            if ~isfield(state, 'type')
                state.type = 'cell';
            end
            value = model.disc.evaluateProp(state, dof, state.type);
        end
        
    end
    
end