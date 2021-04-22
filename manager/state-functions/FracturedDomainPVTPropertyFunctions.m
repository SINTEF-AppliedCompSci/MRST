classdef FracturedDomainPVTPropertyFunctions < PVTPropertyFunctions
    % PVT property functions for fractured domains
    
    methods
        function props = FracturedDomainPVTPropertyFunctions(model)
            % Add in pv for fracture cells to please PoreVolume constructor
            pv0 = model.operators.pv;
            model.operators.pv = [model.operators.pv; 
                   ones(model.G.cells.num - numel(model.operators.pv), 1)];
            % Call parent constructor
            props@PVTPropertyFunctions(model);
            % Set pv back to what it was
            model.operators.pv = pv0;
            % Replace pore volume
            pvt = props.getRegionPVT(model);
            props = props.setStateFunction('PoreVolume', FracturedDomainPoreVolume(model, pvt));
        end
    end
    
end