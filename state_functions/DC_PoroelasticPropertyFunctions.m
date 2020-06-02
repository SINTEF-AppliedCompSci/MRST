classdef DC_PoroelasticPropertyFunctions < StateFunctionGrouping
    % Default grouping for describing rock properties affected by
    % poroelasticity. Namely porosity and permeability. 
    properties
        PoroelasticMatrixPoro % scalar
        PoroelasticFracturePoro %
        Strain % tensor
        EffectiveStress % tensor
        TotalStress % tensor
        PoroelasticFracturePerm % scalar, for now
        PoroelasticFractureTransMult % scalar, for now
    end

    methods
        function props = DC_PoroelasticPropertyFunctions(model)
            props@StateFunctionGrouping();
            props.PoroelasticMatrixPoro = PoroelasticMatrixPoro(model);
            props.PoroelasticFracturePoro = PoroelasticFracturePoro(model);
            props.Strain = Strain(model);
            props.EffectiveStress = EffectiveStress(model);
            props.TotalStress = TotalStress(model);
            props.PoroelasticFracturePerm = PoroelasticFracturePerm(model);
            props.PoroelasticFractureTransMult = PoroelasticFractureTransMult(model);
            
            % Define storage field in state
            props.structName = 'DC_PoroelasticProps';
        end
    end
end