classdef BiotTpfaCompositionalModel < BiotCompositionalModel

    methods
        
        function model = setupStateFunctionGroupings(model, varargin) 
            
            model.PVTPropertyFunctions = []; % make sure this ir reset
            model = setupStateFunctionGroupings@GenericOverallCompositionModel(model, varargin{:});
           
            biotprops = model.BiotPropertyFunctions; 
            pvtprops  = model.PVTPropertyFunctions; 
            mprops  = model.MechPropertyFunctions;
            
            pv = pvtprops.getStateFunction('PoreVolume');
            
            biotprops = biotprops.setStateFunction('BasePoreVolume'               , pv);
            biotprops = biotprops.setStateFunction('Dilatation'                   , BiotBlackOilDilatation(model));
            pvtprops  =  pvtprops.setStateFunction('PoreVolume'                   , BiotPoreVolume(model));
            mprops    =    mprops.setStateFunction('FaceNodeDisplacement'         , BiotFaceNodeDisplacement(model));
            
            model.BiotPropertyFunctions = biotprops;
            model.PVTPropertyFunctions  = pvtprops;
            model.MechPropertyFunctions = mprops;
            
        end
        
    end
end

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}

