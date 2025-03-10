classdef BactConvertionRate <  StateFunction
    % The bacterial growth rate, given per cell
    properties
    end
    
    methods
        function gp = BactConvertionRate(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PsiGrowthRate', 'PsiDecayRate'}, 'FlowDiscretization');
           %gp = gp.dependsOn({'nbact', 's'}, 'state');
           %gp = gp.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');
            gp.label = 'Q_biot';
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            qbiot = 0;
         if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase
             
             ncomp = model.ReservoirModel.EOSModel.getNumberOfComponents;

             Psigrowth = model.getProps(state, 'PsiGrowthRate'); 

             Y_H2 = model.ReservoirModel.Y_H2;
             gammak =model.ReservoirModel.gammak;
             
             qbiot_temp =  (Psigrowth)./Y_H2;
             qbiot = cell(ncomp,1);
             mc = model.ReservoirModel.EOSModel.CompositionalMixture.molarMass;
             for c = 1:ncomp            
                qbiot{c} = gammak(c)/2.*qbiot_temp.*mc(c) +0;
             end
         end         
        end
    end
end