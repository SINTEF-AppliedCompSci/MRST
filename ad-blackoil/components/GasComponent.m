classdef GasComponent < ImmiscibleComponent
    properties
    end
    
    methods
        function c = GasComponent(name, gasIndex)
            c@ImmiscibleComponent(name, gasIndex);
        end
        
        function c = getComponentDensity(component, model, state)
            c = getComponentDensity@ComponentImplementation(component, model, state);
            if model.disgas
                phasenames = model.getPhaseNames();
                gix = phasenames == 'G';
                oix = phasenames == 'O';
                [b, rs] = model.getProps(state, 'ShrinkageFactors', 'rs');
                rhoS = model.getSurfaceDensities();
                rhoGS = rhoS(gix);
                bO = b{oix};
                bG = b{gix};
                c{oix} = rs.*rhoGS.*bO;
                c{gix} = rhoGS.*bG;
            end
        end
    end
end