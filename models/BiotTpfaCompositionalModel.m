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
