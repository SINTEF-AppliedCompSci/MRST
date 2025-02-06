classdef BactConvertionRate <  StateFunction
    % The bacterial growth rate, given per cell
    properties
    end
    
    methods
        function gp = BactConvertionRate(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PsiGrowthRate'}, 'FlowDiscretization');
            gp = gp.dependsOn({'nbact'}, 'state');
            gp.label = 'Q_biot';
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            qbiot = 0;
         if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase
             
             ncomp = model.ReservoirModel.EOSModel.getNumberOfComponents;

             Psigrowth = model.getProps(state, 'PsiGrowthRate'); 

             nbact = model.ReservoirModel.getProps(state, 'nbact');
             pv = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'PoreVolume');

             Y_H2 = model.ReservoirModel.Y_H2;
             gammak =model.ReservoirModel.gammak;
             m_rate =model.ReservoirModel.m_rate;
             nbactMax = model.ReservoirModel.nbactMax;
             bact_limit = 1 - (nbact./nbactMax).^0.5;
             qbiot_temp =  (Psigrowth)./Y_H2;
             qbiot =cell(ncomp,1);
             
             for c = 1:ncomp            
                qbiot{c} = pv.*gammak(c).*qbiot_temp +0;
            end
         end         
        end
    end
end