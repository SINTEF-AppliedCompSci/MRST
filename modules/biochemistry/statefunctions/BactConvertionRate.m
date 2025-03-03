classdef BactConvertionRate <  StateFunction
    % The bacterial growth rate, given per cell
    properties
    end
    
    methods
        function gp = BactConvertionRate(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PsiGrowthRate', 'PsiDecayRate'}, 'FlowDiscretization');
            gp = gp.dependsOn({'nbact', 's'}, 'state');
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            gp.label = 'Q_biot';
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            qbiot = 0;
         if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase
             
             ncomp = model.ReservoirModel.EOSModel.getNumberOfComponents;

             Psigrowth = model.getProps(state, 'PsiGrowthRate'); 
             Psidecay = model.getProps(state, 'PsiDecayRate'); 


             nbact = model.ReservoirModel.getProps(state, 'nbact');
             pv = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'PoreVolume');
             rho = model.ReservoirModel.PVTPropertyFunctions.get(model, state, 'Density');

             rhoS = model.ReservoirModel.getSurfaceDensities();

             sat = model.ReservoirModel.getProps(state, 's');

             Y_H2 = model.ReservoirModel.Y_H2;
             gammak =model.ReservoirModel.gammak;
             m_rate =model.ReservoirModel.m_rate;
             nbactMax = model.ReservoirModel.nbactMax;
             bact_limit = 1 - (nbact./nbactMax).^0.5;
             qbiot_temp =  (Psigrowth)./Y_H2;
             qbiot =cell(ncomp,1);
             L_ix = model.ReservoirModel.getLiquidIndex();
             N=sum(nbact.*sat{L_ix}.*pv);
             Yp =1;
             Yc = 0.25;
             for c = 1:ncomp            
                qbiot{c} = rho{L_ix}.*gammak(c)/2.*qbiot_temp +0;
             end
           %   qbiot{3} = qbiot{3} -  rho{L_ix}.*Psidecay./1.0e9;
         end         
        end
    end
end