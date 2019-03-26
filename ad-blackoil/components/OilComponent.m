classdef OilComponent < ImmiscibleComponent
    properties
        disgas
        vapoil
    end
    
    methods
        function c = OilComponent(name, gasIndex, disgas, vapoil)
            c@ImmiscibleComponent(name, gasIndex);
            c.disgas = disgas;
            c.vapoil = vapoil;
            c = c.dependsOn('ShrinkageFactors');
            if vapoil
                c = c.dependsOn('rv', 'state');
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
                rhoOS = rhoS(oix);
                if component.disgas
                    bO = b{oix};
                    c{oix} = rhoOS.*bO;
                end
                if component.vapoil
                    bG = b{gix};
                    rv = model.getProp(state, 'rv');
                    c{gix} = rv.*rhoOS.*bG;
                end
            end
        end
    end
end