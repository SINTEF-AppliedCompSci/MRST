classdef FracturedDomainFlowDiscretization < FlowDiscretization
    % Flow discretization for fractured domains
    
    methods
        function props = FracturedDomainFlowDiscretization(model)
            % Call parent constructor
            props@FlowDiscretization(model);
            % Replace state functions for gravity potential difference,
            % face mobility and phase flux with fractured domain variants
            props = props.setStateFunction('GravityPotentialDifference', FracturedDomainGravityPotentialDifference(model));
            props = props.setStateFunction('FaceMobility', FracturedDomainFaceMobility(model));
            props = props.setStateFunction('PhaseFlux', FracturedDomainPhaseFlux(model));
        end
    end
    
end