classdef RsGasComponent < ImmiscibleComponent
    properties
    end
    
    methods
        function c = RsGasComponent(name)
            c@ImmiscibleComponent(name, []);
        end
        
        function c = getComponentDensity(component, model, state)
            c = getComponentDensity@ComponentImplementation(component, model, state);
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