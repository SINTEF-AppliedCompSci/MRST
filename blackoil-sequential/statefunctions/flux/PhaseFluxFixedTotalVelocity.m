classdef PhaseFluxFixedTotalVelocity < StateFunction
    properties
    end
    
    methods
        function pf = PhaseFluxFixedTotalVelocity(varargin)
            pf@StateFunction(varargin{:});
            pf = pf.dependsOn({'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility'});
        end
        function f = evaluateOnDomain(prop, model, state)
            [vT, G, T, fmob, mobT] = prop.getEvaluatedDependencies(state,...
                 'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility');
            nph = numel(fmob);
            f = cell(1, nph);
            for i = 1:nph
                mobG = 0;
                for j = 1:nph
                    if i ~= j
                        mobG = mobG + fmob{j}.*(G{i} - G{j});
                    end
                end
                f{i} = (fmob{i}./mobT).*(vT + T.*mobG);
            end
        end
    end
end