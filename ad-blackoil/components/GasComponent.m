classdef GasComponent < ImmiscibleComponent
    properties
        disgas
        vapoil
    end
    
    methods
        function c = GasComponent(name, gasIndex, disgas, vapoil)
            c@ImmiscibleComponent(name, gasIndex);
            c.disgas = disgas;
            c.vapoil = vapoil;
            c = c.dependsOn('ShrinkageFactors');
            if disgas
                c = c.dependsOn('rs', 'state');
            end
        end
        
        function c = getComponentDensity(component, model, state)
            c = getComponentDensity@ImmiscibleComponent(component, model, state);
            if component.disgas || component.vapoil
                phasenames = model.getPhaseNames();
                gix = phasenames == 'G';
                oix = phasenames == 'O';
                b = model.getProps(state, 'ShrinkageFactors');
                rhoS = model.getSurfaceDensities();
                rhoGS = rhoS(gix);
                if component.vapoil
                    bG = b{gix};
                    c{gix} = rhoGS.*bG;
                end
                if component.disgas
                    bO = b{oix};
                    rs = model.getProp(state, 'rs');
                    c{oix} = rs.*rhoGS.*bO;
                end
            end
        end
    end
end