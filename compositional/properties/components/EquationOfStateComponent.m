classdef EquationOfStateComponent < ComponentImplementation
    properties
        componentIndex % Global component numbering
    end
    
    methods
        function c = EquationOfStateComponent(name, cindex)
            c@ComponentImplementation(name);
            c.componentIndex = cindex;
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
            mw = model.EOSModel.fluid.molarMass;
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            ix = component.componentIndex - model.water;
            % Just put "lighter" components in gas phase by default
            c{component.componentIndex} = double(mw(ix) > mean(mw));
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            c = component.getPhaseComposition(model, state);
        end
        
        function c = getPhaseComponentFractionWell(component, model, state, W)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(W.components);
            comp_i = model.EOSModel.getMassFraction(comp_i);
            index = component.componentIndex - model.water;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                c{index} = ci;
            end
        end
    end
end