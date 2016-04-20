classdef CompositionalFluid
    properties
        % Critical temperature in Kelvin
        Tcrit
        % Critical pressure in Pascal
        Pcrit
        % Critical volume in m^3 / mol
        Vcrit
        % Acentric factors (dimensionless)
        acentricFactors
        % Component mass (kg / mol)
        molarMass
        T_ref
        P_ref
        % Names of each component. Must be unique.
        names
    end
    
    properties ( Access = private)
        % Binary interaction coefficients
        bic
    end
    
    methods
        function fluid = CompositionalFluid(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass, T_ref, P_ref)
            if nargin == 0
                disp('Empty fluid - no validation');
                return
            end
            
            fluid.names = names;
            fluid.Tcrit = Tcrit;
            fluid.Pcrit = Pcrit;
            fluid.Vcrit = Vcrit;
            fluid.acentricFactors = acentricFactors;
            fluid.molarMass = molarMass;
            fluid.T_ref = T_ref;
            fluid.P_ref = P_ref;
            ncomp = numel(fluid.names);
            
            fluid.bic = zeros(ncomp, ncomp);
        end
        
        function n = getNumberOfComponents(fluid)
            n = numel(fluid.names);
        end
        
        function bic = getBinaryInteraction(fluid)
            bic = fluid.bic;
        end
        
        function fluid = setBinaryInteraction(fluid, input)
            % Set BIC via a matrix. Must be symmetric and ncomp by ncomp
            ncomp = fluid.getNumberOfComponents();
            assert(size(input, 1) == ncomp)
            assert(size(input, 1) == size(input, 2));
            assert(all(all(input == input')));
            fluid.bic = input;
        end
    end
    
end