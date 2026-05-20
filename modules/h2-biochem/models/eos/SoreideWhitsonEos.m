classdef SoreideWhitsonEos < EquationOfStateModel
    % Soreide-Whitson equation of state model for H2O-containing mixtures.
    %
    % This EOS extends the Peng-Robinson model with modified mixing rules
    % for the liquid phase to account for water-salt interactions. It uses
    % different binary interaction coefficients for liquid and vapor phases,
    % and a modified alpha function for the H2O component in the liquid
    % phase that depends on salt molality.
    %
    % SYNOPSIS:
    %   eos = SoreideWhitsonEos(G, fluid)
    %   eos = SoreideWhitsonEos(G, fluid, 'msalt', 5)
    %
    % REQUIRED PARAMETERS:
    %   G     - Simulation grid (can be empty [])
    %   fluid - CompositionalMixture instance
    %
    % OPTIONAL PARAMETERS (key/value):
    %   msalt - Salt molality (default: 0.0)
    %
    % SEE ALSO:
    %   EquationOfStateModel

    properties
        msalt = 0.0; % Salt molality for Soreide-Whitson EoS
    end

    methods
        function model = SoreideWhitsonEos(G, fluid, varargin)
            % Constructor: set up as Peng-Robinson, then adjust omegaB
            model = model@EquationOfStateModel(G, fluid, 'pr');
            model.omegaB = 0.0777961;
            model = merge_options(model, varargin{:});
        end

        function t = shortname(model) %#ok<MANU>
            t = 'sw';
        end

        function [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = getMixtureFugacityCoefficients(model, P, T, x, y, acf)
            % For Soreide-Whitson, use H2O-corrected mixing for liquid
            % and standard PR mixing for vapor
            [A_ij, Bi] = model.getMixingParametersH2O(P, T, acf, iscell(x));
            [Si_L, A_L, B_L] = model.getPhaseMixCoefficients(x, A_ij, Bi);
            [A_ij, Bi] = model.getMixingParameters(P, T, acf, iscell(y));
            [Si_V, A_V, B_V] = model.getPhaseMixCoefficients(y, A_ij, Bi);
        end

        function [A_ij, Bi] = getLiquidMixingParameters(model, P, T, acf, useCell)
            % Override to use H2O-corrected mixing for liquid phase
            namecp = model.CompositionalMixture.names;
            indH2O = find(strcmp(namecp, 'Water') | strcmp(namecp, 'H2O'), 1);
            if ~isempty(indH2O)
                [A_ij, Bi] = model.getMixingParametersH2O(P, T, acf, useCell);
            else
                [A_ij, Bi] = model.getMixingParameters(P, T, acf, useCell);
            end
        end

        function L = singlePhaseLabel(eos, p, T, z, Z_L, Z_V)
            % Fugacity-based phase labeling for Soreide-Whitson
            Vc = eos.CompositionalMixture.Vcrit;
            Tc = eos.CompositionalMixture.Tcrit;
            Pc = eos.CompositionalMixture.Pcrit;

            Vz = bsxfun(@times, Vc, z);

            Tc_est = sum(bsxfun(@times, Vz, Tc), 2)./sum(Vz, 2);
            Pc_est = sum(bsxfun(@times, Vz, Pc), 2)./sum(Vz, 2);

            Tr = T ./ Tc_est;
            Pr = p ./ Pc_est;

            if ~isempty(p)
                if Tr > 1
                    if Pr > 1
                        L = 0; % Supercritical
                    else
                        if abs(Z_V - Z_L) < 0.1
                            if Z_V < 0.5
                                L = 1;
                            else
                                L = 0;
                            end
                        elseif Z_V > Z_L
                            L = 0;
                        else
                            L = 1;
                        end
                    end
                else
                    if abs(Z_V - Z_L) < 0.1
                        if Z_V < 0.5
                            L = 1;
                        else
                            L = 0;
                        end
                    elseif Z_V > Z_L
                        L = 0;
                    else
                        L = 1;
                    end
                end
            else
                L = [];
            end
        end

        function [A_ij, Bi] = getMixingParametersH2O(model, P, T, acf, useCell)
            % H2O-corrected mixing parameters for the liquid phase.
            % Uses the Soreide-Whitson alpha function for H2O and
            % liquid-phase BICs that depend on salt molality.
            if nargin < 5
                useCell = true;
            end

            ncomp = model.getNumberOfComponents();
            [Pr, Tr] = model.getReducedPT(P, T, useCell);

            mSalt = model.msalt;
            coef_msalt = mSalt^1.1;
            namecomp = model.CompositionalMixture.names;
            indH2O = find(strcmp(namecomp, 'H2O'));

            if useCell
                [sAi, Bi] = deal(cell(1, ncomp));
                [oA, oB] = deal(cell(1, ncomp));
                [oB{:}] = deal(model.omegaB);
            else
                oB = model.omegaB;
            end

            % In the liquid phase, use PR alpha except for H2O where
            % the Soreide-Whitson salt-corrected alpha is used
            if useCell
                for i = 1:ncomp
                    ai = acf(i);
                    if ai > 0.49
                        oA{i} = model.omegaA.*(1 + (0.379642 + 1.48503.*ai - 0.164423.*ai.^2 + 0.016666.*ai.^3).*(1-Tr{i}.^(1/2))).^2;
                    else
                        oA{i} = model.omegaA.*(1 + (0.37464 + 1.54226.*ai - 0.26992.*ai.^2).*(1-Tr{i}.^(1/2))).^2;
                    end
                end
                % Override alpha for H2O component
                TrH2O = Tr{indH2O};
                oA{indH2O} = model.omegaA.*(1 + 0.4530.*(1-(1-0.0103*coef_msalt).*TrH2O)+0.0034.*(TrH2O.^(-3)-1)).^2;
            else
                tmp = bsxfun(@times, 0.37464 + 1.54226.*acf - 0.26992.*acf.^2, 1-Tr.^(1/2));
                oA = model.omegaA.*(1 + tmp).^2;
                % Override alpha for H2O component
                TrH2O = Tr(:, indH2O);
                tmp1 = bsxfun(@power, 1 + 0.4530.*(1-(1-0.0103*coef_msalt).*TrH2O)+0.0034.*(TrH2O.^(-3)-1), 2);
                oA(:, indH2O) = model.omegaA.*tmp1;
            end

            bic = model.getBinaryInteractionLiquidWater(T, mSalt);

            if useCell
                A_ij = cell(ncomp, ncomp);
                for i = 1:ncomp
                    sAi{i} = ((oA{i}.*Pr{i}).^(1/2))./Tr{i};
                    Bi{i} = oB{i}.*Pr{i}./Tr{i};
                end
                if iscell(bic)
                    for i = 1:ncomp
                        for j = i:ncomp
                            A_ij{i, j} = (sAi{i}.*sAi{j}).*(1 - bic{i, j});
                            A_ij{j, i} = A_ij{i, j};
                        end
                    end
                else
                    for i = 1:ncomp
                        for j = i:ncomp
                            A_ij{i, j} = (sAi{i}.*sAi{j}).*(1 - bic(i, j));
                            A_ij{j, i} = A_ij{i, j};
                        end
                    end
                end
            else
                Ai = oA.*Pr./Tr.^2;
                Bi = oB.*Pr./Tr;
                A_ij = cell(ncomp, 1);
                if iscell(bic)
                    for j = 1:ncomp
                        bicc = [];
                        for i = 1:ncomp
                            bicc(:, i) = bic{j, i};
                        end
                        A_ij{j} = bsxfun(@times, bsxfun(@times, Ai, Ai(:, j)).^(1/2), 1 - bicc);
                    end
                else
                    for j = 1:ncomp
                        A_ij{j} = bsxfun(@times, bsxfun(@times, Ai, Ai(:, j)).^(1/2), 1 - bic(j, :));
                    end
                end
            end
        end

        function bic = getBinaryInteractionGasWater(model, T)
            % Gas-phase BIC calculation for components with H2O.
            % References:
            %   Ref1: A. Afanasyev & al., https://doi.org/10.1016/j.jngse.2021.103988
            %   Ref2: S. Chabab & al., Energies 2021,14, 5239. https://doi.org/10.3390/en14175239
            %   Ref3: S. Chabab & al., https://doi.org/10.1016/j.ijhydene.2023.10.290
            fluid = model.CompositionalMixture;
            namecp = fluid.names;
            indices = struct('H2', find(strcmp(namecp, 'H2')), ...
                'C1', find(strcmp(namecp, 'C1')), ...
                'CO2', find(strcmp(namecp, 'CO2')), ...
                'H2S', find(strcmp(namecp, 'H2S')), ...
                'N2', find(strcmp(namecp, 'N2')), ...
                'C2', find(strcmp(namecp, 'C2')), ...
                'C3', find(strcmp(namecp, 'C3')), ...
                'NC4', find(strcmp(namecp, 'NC4')), ...
                'H2O', find(strcmp(namecp, 'H2O')));

            nComponents = numel(namecp);
            bic = cell(nComponents);

            coeffs = struct(...
                'H2',  [0.01993, 0.042834], ...
                'C1', [0.494435, 0.0], ...
                'H2S', [0.19031, -0.05965], ...
                'CO2', [-2.066623464504e-2, 0.207440935], ...
                'N2',  [0.385438, 0.0], ...
                'C2',  [0.4920, 0.0], ...
                'C3',  [0.5525, 0.0], ...
                'NC4', [0.5091, 0.0]);

            fields = fieldnames(indices);
            for i = 1:numel(fields)
                comp = fields{i};
                indComp = indices.(comp);
                indH2O = indices.H2O;
                if ~isempty(indComp) && ~isempty(indH2O) && isfield(coeffs, comp)
                    Tr = T ./ fluid.Tcrit(indComp);
                    bic{indComp, indH2O} = coeffs.(comp)(1) + coeffs.(comp)(2) .* Tr;
                    bic{indH2O, indComp} = bic{indComp, indH2O};
                end
            end
            for i = 1:nComponents
                Tr = T ./ fluid.Tcrit(i);
                bic{i, i} = fluid.getBinaryInteraction();
                bic{i, i} = bic{i, i}(i, i) + 0.*Tr;
            end
            for i = 1:nComponents
                for j = 1:nComponents
                    if isempty(bic{i, j})
                        Tr = T ./ fluid.Tcrit(i);
                        bic{i, j} = 0 .* Tr;
                    end
                end
            end
        end

        function bic = getBinaryInteractionLiquidWater(model, T, msalt)
            % Liquid-phase BIC calculation for components with H2O.
            % References:
            %   Ref1: A. Afanasyev & al., https://doi.org/10.1016/j.jngse.2021.103988
            %   Ref2: S. Chabab & al., Energies 2021,14, 5239. https://doi.org/10.3390/en14175239
            %   Ref3: S. Chabab & al., https://doi.org/10.1016/j.ijhydene.2023.10.290
            fluid = model.CompositionalMixture;
            namecp = fluid.names;
            indices = struct('H2', find(strcmp(namecp, 'H2')), ...
                'C1', find(strcmp(namecp, 'C1')), ...
                'CO2', find(strcmp(namecp, 'CO2')), ...
                'H2S', find(strcmp(namecp, 'H2S')), ...
                'N2', find(strcmp(namecp, 'N2')), ...
                'C2', find(strcmp(namecp, 'C2')), ...
                'C3', find(strcmp(namecp, 'C3')), ...
                'NC4', find(strcmp(namecp, 'NC4')), ...
                'H2O', find(strcmp(namecp, 'H2O')));

            nComponents = numel(namecp);
            bic = cell(nComponents);

            coeffs = struct('H2', [-2.11917, 0.14888, -13.01835, -0.43946, -2.26322e-2, -4.4736e-3, 0, 0, 1], ...
                'C1', [-1.625685, 1.114873, 0,0, 8.590105e-21, 1.812763e-3, -0.169968, -4.198569e-2,1], ...
                'CO2', [-1.709096, 0.450487, 0, 0,  1.792130e-2,  0.066426, 0,0,1], ...
                'H2S', [-4.2619e-1, 6.73586E-1, 0, 0,  -0.0575,  -0.078823343,-2.16250E-1,-0.160085087,1], ...
                'N2',  [-1.702359096, 0.450487, 0, 0, 1.792130E-2,0.066426, 0,0,0.8], ...
                'C2', [ 1.1120-1.7369.*0.0990.^(-0.1),  1.1001 + 0.8360.*0.0990,0,0,0.017407,0.033516,-0.15742-1.0988.*0.0990,0.011478,1], ...
                'C3', [ 1.1120-1.7369.*0.1520.^(-0.1),  1.1001 + 0.8360.*0.1520,0,0,0.017407,0.033516,-0.15742-1.0988.*0.1520,0.011478,1], ...
                'NC4', [ 1.1120-1.7369.*0.200810094644.^(-0.1),1.1001+0.8360.*0.200810094644,0,0,0.017407,0.033516,-0.15742-1.0988.*0.200810094644,0.011478,1]);

            fields = fieldnames(indices);
            for i = 1:numel(fields)
                comp = fields{i};
                indComp = indices.(comp);
                indH2O = indices.H2O;
                if ~isempty(indComp) && ~isempty(indH2O) && isfield(coeffs, comp)
                    Tr = T ./ fluid.Tcrit(indComp);
                    csalt = coeffs.(comp)(9);
                    bic{indComp, indH2O} = coeffs.(comp)(1) .* (1 + coeffs.(comp)(5) .* msalt.^csalt) + ...
                        coeffs.(comp)(2) .* Tr .* (1 + coeffs.(comp)(6) .* msalt.^csalt) + ...
                        coeffs.(comp)(3) .* exp(coeffs.(comp)(4) .* Tr)+...
                        +coeffs.(comp)(7) .* Tr.^2.*(1+coeffs.(comp)(8) .* msalt.^csalt);
                    bic{indH2O, indComp} = bic{indComp, indH2O};
                end
            end
            indCO2 = indices.CO2;
            indH2O = indices.H2O;
            if ~isempty(indCO2) && ~isempty(indH2O)
                a = 0.43575155;
                b = -5.766906744e-2;
                c = 8.26464849e-3;
                d = 1.29539193e-3;
                e = -1.6698848e-3;
                f = -0.47866096;
                Tr_CO2 = T ./ fluid.Tcrit(indices.CO2);
                bic{indH2O, indCO2} = Tr_CO2 .* ...
                    (a + b .* Tr_CO2 + c .* Tr_CO2 .* msalt) + ...
                    msalt^2 .* (d + e .* Tr_CO2) + ...
                    f;
                bic{indCO2, indH2O} = bic{indH2O, indCO2};
            end
            baseBic = fluid.getBinaryInteraction();
            for i = 1:nComponents
                Tr = T ./ fluid.Tcrit(i);
                bic{i, i} = baseBic(i, i) + 0.*Tr;
            end
            for i = 1:nComponents
                for j = 1:nComponents
                    if isempty(bic{i, j})
                        bic{i, j} = 0 .* Tr;
                    end
                end
            end
        end
    end
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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
