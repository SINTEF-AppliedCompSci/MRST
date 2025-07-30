classdef BactConvertionRate < StateFunction
    % Calculates the bacterial growth rate per cell.
    
    properties
    end

    methods
        function gp = BactConvertionRate(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PsiGrowthRate', 'PsiDecayRate'}, 'FlowDiscretization');
            gp = gp.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');
            gp.label = 'Q_biot';
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            qbiot = 0; % Default initialization

            % Retrieve component names and indices for H2 and CO2
            compn = model.ReservoirModel.EOSModel.getComponentNames();
            idxH2 = find(strcmp(compn, 'H2'));  % Index for H2
            idxCO2 = find(strcmp(compn, 'CO2')); % Index for CO2

            % Check if bacterial modeling is enabled and components are valid
            if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase && ...
               ~isempty(idxH2) && ~isempty(idxCO2)

                ncomp = model.ReservoirModel.EOSModel.getNumberOfComponents();
                qbiot = cell(ncomp, 1);

                % Retrieve required properties from the model
                pv = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'PoreVolume');
                s = model.ReservoirModel.getProps(state, 's');
                nbact = model.ReservoirModel.getProps(state, 'nbact');
                L_ix = model.ReservoirModel.getLiquidIndex();
                x = model.ReservoirModel.getProps(state, 'x');

                % Handle liquid phase properties
                if iscell(x)
                    xH2 = x{idxH2};
                    xCO2 = x{idxCO2};
                    sL = s{L_ix};
                else
                    xH2 = x(:, idxH2);
                    xCO2 = x(:, idxCO2);
                    sL = s(:, L_ix);
                end

                if iscell(sL)
                    Voln = sL{1};
                else
                    Voln = sL;
                end

                Voln = max(Voln, 1.0e-8); % Ensure non-zero volume

                % Growth rate parameters
                alphaH2 = model.ReservoirModel.alphaH2;
                alphaCO2 = model.ReservoirModel.alphaCO2;
                Psigrowthmax = model.ReservoirModel.Psigrowthmax;

                % Calculate Psigrowth using mole fractions
                axH2 = (xH2 ./ (alphaH2 + xH2));
                axCO2 = (xCO2 ./ (alphaCO2 + xCO2));

                Psigrowth = pv .* Psigrowthmax .* axH2 ...
                            .* axCO2 .* nbact .* Voln;

                % Conversion factors
                Y_H2 = model.ReservoirModel.Y_H2;
                gammak = model.ReservoirModel.gammak;
                nbactMax = model.ReservoirModel.nbactMax;
                mc = model.ReservoirModel.EOSModel.CompositionalMixture.molarMass;

                % Calculate the biotic reaction rate per component
                qbiot_temp = Psigrowth ./ Y_H2 ./ abs(gammak(idxH2));
                for c = 1:ncomp
                    qbiot{c} = gammak(c) .* qbiot_temp .* mc(c) .* nbactMax;
                end
            end
        end
    end
end