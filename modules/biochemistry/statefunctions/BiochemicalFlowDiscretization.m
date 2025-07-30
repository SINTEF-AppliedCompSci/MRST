classdef BiochemicalFlowDiscretization < FlowDiscretization
%Discretization and state function grouping for bio-chemistry flow
    
    properties
              
        MolecularDiffPhaseFlux
        PsiGrowthRate
        BactConvRate
        PsiDecayRate
        BactFlux
    end
    
    methods
        %-----------------------------------------------------------------%
        function props = BiochemicalFlowDiscretization(model)
            % Inherit most of the state functions from FLuxDiscretization
            props = props@FlowDiscretization(model);
            % Fluid flow transmissibility
            if model.dynamicFlowTrans
                % Dynamic transmissibility. Conductivities are caluclated
                % at each nonlinear iteration, and used to compute
                % corresponding transmissibilities
                props = props.setStateFunction('Permeability', BactPermeability(model));
                props = props.setStateFunction('Transmissibility', ...
                                                DynamicFlowTransmissibility(model, 'Permeability'));
            else
                % Static transmissibilities already been set up by parent
                % or set up BactTransmissibility
            end

            if model.dynamicFlowPv
                % Dynamic transmissibility. Conductivities are caluclated
                % at each nonlinear iteration, and used to compute
                % corresponding transmissibilities
                props = props.setStateFunction('Porosity', BactPorosity(model));
                props = props.setStateFunction('PoreVolume', ...
                    DynamicFlowPoreVolume(model, 'Porosity'));
            else
                % Static transmissibilities already been set up by parent
            end
            props = props.setStateFunction('MolecularDiffPhaseFlux', ComponentMolecularDiffPhaseFlux(model));
            props = props.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
            props = props.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
            props = props.setStateFunction('BactConvRate', BactConvertionRate(model));
            props = props.setStateFunction('BactFlux', DiffusiveBactFlux(model));
            if model.bDiffusionEffect
                props = props.setStateFunction('BactFlux'          , DiffusiveBactFlux(model)          );
            end
            if model.moleculardiffusion
                % Set molecular diffusion flux
                props = props.setStateFunction('MolecularDiffusivity'       , ComponentMolecularDiffPhaseFlux(model));
            end
        end
        
        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
            flowState = fd.buildFlowState(model, state, state0, dt);
            act = model.getActivePhases();
            ncomp = model.getNumberOfComponents();
            nph = sum(act); 
            if model.moleculardiffusion           
                J = model.getProps(flowState, 'MolecularDiffPhaseFlux');
                for c = 1:ncomp
                    for ph = 1:nph                   
                        flux{c} = flux{c} + J{c,ph};                     
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function [acc, bflux, name, type] = bacteriaConservationEquation(fd, model, state, state0, dt)
            bactmass = model.getProps(state, 'BacterialMass');
            bactmass0 = model.getProps(state0, 'BacterialMass');
            flowState = fd.buildFlowState(model, state, state0, dt);
            bflux = [];
            if model.bDiffusionEffect
                bflux      = model.getProp(flowState, 'BactFlux');
            end
            % Add up accumulation
            acc = (bactmass  - bactmass0 )./dt;
            % Convert to cell arrays to comply with AD framework
            acc  = {acc};
            bflux = {bflux};
            % Name and type
            name = 'bacteria';
            type = 'cell';
        end

        function c = getPopulationMass(model, state, extra)
            c = component.getPhaseComposition(model, state);
            if nargin < 4
                rho = model.getProps(state, 'Density');
            else
                rho = extra.rho;
            end
            for ph = 1:numel(c)
                if ~isempty(c{ph})
                    c{ph} = rho{ph}.*c{ph};
                end
            end
        end
    end
    
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}