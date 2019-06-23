classdef EquationOfStateComponent < ComponentImplementation
    properties
        componentIndex % Global component numbering
        surfacePhaseMassFractions % Mass fraction for each phase
        surfacePhaseDensityPure % Density
        pressure % pressure
        T % Temp
    end
    
    methods
        function c = EquationOfStateComponent(name, p, T, cindex, surfaceMassFractions, density)
            c@ComponentImplementation(name);
            c.componentIndex = cindex;
            c.pressure = p;
            c.T = T;
            c.surfacePhaseMassFractions = surfaceMassFractions;
            c.surfacePhaseDensityPure = density;
            c = c.dependsOn({'Density', 'ComponentPhaseMassFractions'});
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            c = getComponentDensity@ComponentImplementation(component, model, state, varargin{:});
            [rho, massFractions] = model.getProps(state, 'Density', 'ComponentPhaseMassFractions');
            ix = component.componentIndex;
            for ph = 1:numel(c)
                mf = massFractions{ix, ph};
                if ~isempty(mf)
                    c{ph} = rho{ph}.*mf;
                end
            end
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            for i = 1:nph
                c{i} = component.surfacePhaseMassFractions(i);
            end
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            c = component.getPhaseComposition(model, state);
        end

        function c = getPurePhaseDensitySurface(component, model, state, pressure, temperature)
            % Surface density, for a pure component
            rho = component.surfacePhaseDensityPure;
            % Ideal gas scaling
            scale = (pressure./component.pressure).*(component.T./temperature);
            rho = bsxfun(@times, rho, scale);
            c = arrafyn(@(x) x, rho, 'UniformOutput', false);
        end

        function c = getPhaseComponentFractionWell(component, model, state, W)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(W.components);
            comp_i = model.EOSModel.getMassFraction(comp_i);
            index = component.componentIndex - model.water;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                for i = (1+model.water):nph
                    c{i} = ci;
                end
            end
        end
    end
end