classdef ComponentMolecularDiffPhaseFlux < StateFunction
    % Flux of each component, in each phase
    properties
    end

    methods
        function gp = ComponentMolecularDiffPhaseFlux(model, varargin)
            gp@StateFunction(model);                
            gp = gp.dependsOn({'Density', 'PoreVolume'}, 'PVTPropertyFunctions');
             gp = gp.dependsOn('s', 'state');
             gp = gp.dependsOn('x', 'state');

            gp.label = 'J_{i,\alpha}';
        end

        function J = evaluateOnDomain(prop, model, state)            
            ncomp = model.getNumberOfComponents;            
            nph = model.getNumberOfPhases;
            J = cell(ncomp, nph);
            J = cellfun(@(x) 0, J, 'UniformOutput', false);      
            if isfield(state,'x')
               nm = model.getPhaseNames();
               tau = [1,1];
               rho = prop.getEvaluatedExternals(model, state, 'Density');
               pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');

               avg = model.operators.faceAvg;
               poro = pv./model.G.cells.volumes;
               L_ix = model.getLiquidIndex();
               V_ix = model.getVaporIndex();
               % Define diffusion coefficients in m²/s for liquid and gas phases
               % These are example values, please replace them with actual data as needed
               % Format: [liquid_diff gas_diff] for each component
               for c = 1:ncomp
                   for ph = 1:nph
                       s = model.getProp(state, ['s', nm(ph)]);
                       D_diff = avg(s.*rho{ph}.*model.mol_diff(c,ph).*tau(ph).*poro);%Millington and Quirk model
                       if (ph==L_ix)                    
                           J{c, ph} = - D_diff.*model.operators.Grad(state.x{c});
                       elseif (ph==V_ix)                                              
                           J{c, ph} = - D_diff.*model.operators.Grad(state.y{c});
                       end
                   end
               end
            end
        end
    end
end
