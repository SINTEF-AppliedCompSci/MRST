classdef FluxDiscretization < PropertyFunctions
    properties
        PermeabilityPotentialGradient % K * (grad(p) + rho g dz)
        PressureGradient % Gradient of phase pressures
        GravityPotentialDifference % rho * g * dz term
        PhasePotentialDifference % (grad p_alpha + dpdz)
        PhaseFlux % Phase volumetric fluxes
        PhaseUpwindFlag
        ComponentFlux
        Transmissibility
    end
    
    methods
        function props = FluxDiscretization(model)
            backend = model.AutoDiffBackend;
            upstr = UpstreamFunctionWrapper(model.operators.faceUpstr);
            tpfa = TwoPointFluxApproximation(model);

            props@PropertyFunctions();
            % Darcy flux
            props.Transmissibility = Transmissibility(backend);
            props.PermeabilityPotentialGradient = PermeabilityPotentialGradient(backend, tpfa);
            props.PressureGradient = PressureGradient(backend);
            props.GravityPotentialDifference = GravityPotentialDifference(backend);
            
            % Phase flux
            props.PhasePotentialDifference = PhasePotentialDifference(backend);
            props.PhaseUpwindFlag = PhaseUpwindFlag(backend);
            
            % Fluxes - these are upwinded properties
            props.ComponentFlux = ComponentFlux(backend, upstr);
            props.PhaseFlux = PhaseFlux(backend, upstr);

            % Define storage
            props.structName = 'FluxProps';
        end

        function state = evaluateProperty(props, model, state, name)
            switch name

            end
            state = evaluateProperty@PropertyFunctions(props, model, state, name);
        end
        
        function [acc, div, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;

            acc = cell(ncomp, 1);
            div = cell(ncomp, 1);
            names = model.getComponentNames();
            types = cell(ncomp, 1);
            [types{:}] = deal('cell');
            
            [mass, v] = model.getProps(state, 'ComponentMass', 'ComponentFlux');
            mass0 = model.getProps(state, 'ComponentMass');
            
            Div = model.operators.Div;
            for c = 1:ncomp
                % Loop over phases where the component may be present
                for ph = 1:nph
                    % Check if present
                    if ~isempty(mass{c, ph})
                        accumulated = (mass{c, ph} - mass0{c, ph});
                        dflux = Div(v{c, ph});
                        
                        if isempty(acc{c})
                            acc{c} = accumulated;
                            div{c} = dflux;
                        else
                            acc{c} = acc{c} + accumulated;
                            div{c} = div{c} + dflux;
                        end
                    end
                end
                acc{c} = acc{c}./dt;
            end
        end
    end
end