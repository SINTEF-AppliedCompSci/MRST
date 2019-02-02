classdef GasComponent < ImmiscibleComponent
    properties
    end
    
    methods
        function c = GasComponent(name, gasIndex)
            c@ImmiscibleComponent(name, gasIndex);
        end
        
        function c = getComponentDensity(component, model, state)
            c = getComponentDensity@ComponentImplementation(component, model, state);
            if model.disgas || model.vapoil
                phasenames = model.getPhaseNames();
                gix = phasenames == 'G';
                oix = phasenames == 'O';
                [b, rs] = model.getProps(state, 'ShrinkageFactors', 'rs');
                rhoS = model.getSurfaceDensities();

                if model.vapoil
                    rhoGS = rhoS(gix);
                    bG = b{gix};
                    c{gix} = rhoGS.*bG;
                end

                if model.disgas
                    bO = b{oix};
                    rs = model.getProp(state, 'rs');
                    c{oix} = rs.*rhoGS.*bO;
                end
            end
        end
    end
end