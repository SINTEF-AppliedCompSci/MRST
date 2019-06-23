classdef BasicSeparator
    properties
        pressure
        T
    end
    
    methods
        function sep = BasicSeparator(varargin)
            sep = merge_options(sep, varargin{:});
        end
        
        function [phaseMassStreams, massFractions, densities] = separateComponentMassStream(sep, model, massStream)
            assert(iscell(massStream));
            ncomp = numel(massStream);
            nph = model.getNumberOfPhases();
            % Total mass stream of each phase
            phaseMassStreams = cell(1, nph);
            % Mass density of each phase
            densities = cell(1, nph);
            [phaseMassStreams{:}, densities{:}] = deal(0);
            % Mass fractions
            massFractions = cell(1, nph);
            [massFractions{:}] = deal(cell(1, ncomp));
            % Compute phase split based on simplified assumptions
            % (components act independently of each other)
            for c = 1:ncomp
                comp = model.Components{c};
                frac = comp.getPhaseCompositionSurface(model, [], sep.pressure, sep.T);
                for ph = 1:nph
                    result = massStream{c}*frac{ph};
                    massFractions{ph}{c} = result;
                    phaseMassStreams{ph} = phaseMassStreams{ph} + result;
                end
            end
            % We can now use a simple linear mixing rule to get densities
            for c = 1:ncomp
                rho = comp.getPurePhaseDensitySurface(model, [], sep.pressure, sep.T);
                for ph = 1:nph
                    densities{ph} = densities{ph} + massFractions{ph}{c}.*rho{ph};
                end
            end
        end
        
        function [moleStreams, moleFractions, molardensities] = separateComponentMoleStream(sep, model, molestream)
            assert(false, 'Not implemented in class');
        end
    end
end