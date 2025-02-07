classdef CompositionalMixture
    properties
        Tcrit % Critical temperature in Kelvin
        Pcrit % Critical pressure in Pascal
        Vcrit % Critical volume in m^3 / mol
        acentricFactors % Acentric factors (dimensionless)
        molarMass % Component mass (kg / mol)
        names % Names of each component. Each name must be unique.
        name % Name of the fluid mixture itself
    end
    
    properties ( Access = protected)
        bic % Binary interaction coefficients
    end
    
    methods
        function fluid = CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass, varargin)
            if nargin == 0
                disp('Empty fluid - no validation');
                return
            end
            
            for i = 1:numel(names)
                cts = sum(strcmp(names{i}, names));
                if cts > 1
                    warning(['Component ', num2str(i), ': ', names{i}, ' occurs multiple times.']);
                end
            end
            fluid.names = names;
            fluid.Tcrit = Tcrit;
            fluid.Pcrit = Pcrit;
            fluid.Vcrit = Vcrit;
            fluid.acentricFactors = acentricFactors;
            fluid.molarMass = molarMass;
            ncomp = numel(fluid.names);
            
            fluid.bic = zeros(ncomp, ncomp);
        end
        
        function n = getNumberOfComponents(fluid)
            n = numel(fluid.names);
        end
        
        function bic = getBinaryInteraction(fluid)
            bic = fluid.bic;
        end


function bic = getBinaryInteractionGasWater(fluid, T)
    % Generalized gas-phase BIC calculation for components with H2O.
    %
    % Inputs:
    %   - fluid: structure containing component names, critical temperatures, and BICs
    %   - T: temperature (scalar or vector)
    %
    % Outputs:
    %   - bic: binary interaction coefficients matrix
    % References: 
    %   -Ref1: A. Afanasyev & al., https://doi.org/10.1016/j.jngse.2021.103988
    %   -Ref2: S. Chabab & al., Energies 2021,14, 5239. https://doi.org/10.3390/en14175239
    %   -Ref3: S. Chabab & al., https://doi.org/10.1016/j.ijhydene.2023.10.290
    

    % Extract component indices dynamically
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

    % Initialize BIC matrix
    nComponents = numel(namecp);
    bic = cell(nComponents);

    % Define gas-phase BIC coefficients
    coeffs = struct(...
        'H2',  [0.01993, 0.042834], ...%Ref3
        'C1', [0.494435, 0.0], ... %Ref2
        'H2S', [0.19031, -0.05965], ...%Ref1
        'CO2', [-2.066623464504e-2, 0.207440935], ... %Ref2 (coeffs2=coeffs2(Ref2)*Tc(CO2)
        'N2',  [0.385438, 0.0], ... %Ref2 % 'N2',  [0.4778, 0.0], ... %Ref1
        'C2',  [0.4920, 0.0], ... %Ref1
        'C3',  [0.5525, 0.0], ... %Ref1
        'NC4', [0.5091, 0.0]); %Ref1

    % Assign BIC values for each component interacting with H2O
    fields = fieldnames(indices);
    for i = 1:numel(fields)
        comp = fields{i};
        indComp = indices.(comp);
        indH2O = indices.H2O;

        if ~isempty(indComp) && ~isempty(indH2O) && isfield(coeffs, comp)
            Tr = T ./ fluid.Tcrit(indComp);
            bic{indComp, indH2O} = coeffs.(comp)(1) + coeffs.(comp)(2) .* Tr;
            bic{indH2O, indComp} = bic{indComp, indH2O}; % Symmetry
        end
    end
    % Handle self-interactions for H2, CO2, CH4, N2, and H2O
    for i=1:nComponents
        Tr = T ./ fluid.Tcrit(i);
        bic{i,i} = fluid.bic(i,i)+0.*Tr;
    end
    for i =1:nComponents
        for j=1:nComponents
            if isempty(bic{i,j})
                Tr = T ./ fluid.Tcrit(i);
                bic{i,j} = 0 .* Tr;
            end
        end
    end
end
function bic = getBinaryInteractionLiquidWater(fluid, T, msalt)
    % Generalized liquid-phase BIC calculation for components with H2O.
    %
    % Inputs:
    %   - fluid: structure containing component names, critical temperatures, and BICs
    %   - T: temperature (scalar or vector)
    %   - msalt: salinity (scalar or vector)
    %
    % Outputs:
    %   - bic: binary interaction coefficients matrix
    % References: 
    %   -Ref1: A. Afanasyev & al., https://doi.org/10.1016/j.jngse.2021.103988
    %   -Ref2: S. Chabab & al., Energies 2021,14, 5239. https://doi.org/10.3390/en14175239
    %   -Ref3: S. Chabab & al., https://doi.org/10.1016/j.ijhydene.2023.10.290
    

    % Extract component indices dynamically
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

    % Initialize BIC matrix
    nComponents = numel(namecp);
    bic = cell(nComponents);

    % Define liquid-phase BIC coefficients with salinity influence
    coeffs = struct('H2', [-2.11917, 0.14888, -13.01835, -0.43946, -2.26322e-2, -4.4736e-3, 0, 0, 1], ...%Ref3
                    'C1', [-1.625685, 1.114873, 0,0, 8.590105e-21, 1.812763e-3, -0.169968, -4.198569e-2,1], ...%Ref2
                    'CO2', [-1.709096, 0.450487, 0, 0,  1.792130e-2,  0.066426, 0,0,1], ...%?
                    'H2S', [-4.2619e-1, 6.73586E-1, 0, 0,  -0.0575,  -0.078823343,-2.16250E-1,-0.160085087,1], ...%Ref1
                    'N2',  [-1.702359096, 0.450487, 0, 0, 1.792130E-2,0.066426, 0,0,0.8], ...%Ref2
                    'C2', [ 1.1120-1.7369.*0.0990.^(-0.1),  1.1001 + 0.8360.*0.0990,0,0,0.017407,0.033516,-0.15742-1.0988.*0.0990,0.011478,1], ...%Ref1
                    'C3', [ 1.1120-1.7369.*0.1520.^(-0.1),  1.1001 + 0.8360.*0.1520,0,0,0.017407,0.033516,-0.15742-1.0988.*0.1520,0.011478,1], ...%Ref1
                    'NC4', [ 1.1120-1.7369.*0.200810094644.^(-0.1),1.1001+0.8360.*0.200810094644,0,0,0.017407,0.033516,-0.15742-1.0988.*0.200810094644,0.011478,1]);%Ref1


    % Assign BIC values for each component interacting with H2O
    fields = fieldnames(indices);
    for i = 1:numel(fields)
        comp = fields{i};
        indComp = indices.(comp);
        indH2O = indices.H2O;

        if ~isempty(indComp) && ~isempty(indH2O) && isfield(coeffs, comp)
            Tr = T ./ fluid.Tcrit(indComp);
            csalt=coeffs.(comp)(9);
             bic{indComp, indH2O} = coeffs.(comp)(1) .* (1 + coeffs.(comp)(5) .* msalt.^csalt) + ...
                                   coeffs.(comp)(2) .* Tr .* (1 + coeffs.(comp)(6) .* msalt.^csalt) + ...
                                   coeffs.(comp)(3) .* exp(coeffs.(comp)(4) .* Tr)+...
                                   +coeffs.(comp)(7) .* Tr.^2.*(1+coeffs.(comp)(8) .* msalt.^csalt);
            bic{indH2O, indComp} = bic{indComp, indH2O}; % Symmetry
        end
    end
    indCO2 = indices.CO2;
    indH2O = indices.H2O;
    if ~isempty(indCO2)&&~isempty(indH2O)%Ref2
        a = 0.43575155;
        b = -5.766906744e-2;
        c = 8.26464849e-3;
        d = 1.29539193e-3;
        e = -1.6698848e-3;
        f = -0.47866096;
        Tr_CO2 = T ./fluid.Tcrit(indices.CO2);
        bic{indH2O, indCO2} = Tr_CO2 .* ...
            (a + b .* Tr_CO2 + c .* Tr_CO2 .* msalt) + ...
            msalt^2 .* (d + e .* Tr_CO2) + ...
            f;
        bic{indCO2, indH2O} = bic{indH2O, indCO2};
    end
    % Handle self-interactions for H2, CO2, CH4, N2, and H2O
    for i=1:nComponents
        Tr = T ./ fluid.Tcrit(i);
        bic{i,i} = fluid.bic(i,i)+0.*Tr;
    end

    for i =1:nComponents
        for j=1:nComponents
            if isempty(bic{i,j})
                bic{i,j} = 0 .* Tr;
            end
        end
    end
end


        function fluid = setBinaryInteraction(fluid, input)
            % Set BIC via a matrix. Must be symmetric and ncomp by ncomp
            ncomp = fluid.getNumberOfComponents();
            
            if isvector(input)
                n_el = ncomp*(ncomp-1)/2;
                assert(numel(input) == n_el);
                
                tmp = zeros(ncomp, ncomp);
                offset = 0;
                for i = 2:ncomp
                    cts = i-1;
                    tmp(i, 1:i-1) = input(offset + (1:cts));
                    offset = offset + cts;
                end
                input = tmp + tmp';
            end
            assert(size(input, 1) == ncomp)
            assert(size(input, 1) == size(input, 2));
            assert(all(all(input == input')));
            fluid.bic = input;
        end
        
        function disp(mixture)
            % builtin('disp', mixture);
            cn = class(mixture);
            isDesktop = usejava('desktop');
            if isDesktop
                fprintf('  <a href="matlab:helpPopup %s">%s</a>:\n', cn, cn);
            else
                fprintf('  %s:\n', cn);
            end
            fprintf('\n');
            nc = mixture.getNumberOfComponents();
            cnames = mixture.names;
            fprintf('  %d component mixture', nc);
            if ~isempty(mixture.name)
                fprintf(' (%s)', mixture.name)
            end
            fprintf(':\n');

            ml = max(max(cellfun(@numel, cnames)), 4);
            sep = repmat('-', 55 + ml, 1);
            fprintf('  %*s | p_c [Pa] | T_c [K] | V_c [m^3] |  acf  | mw [kg/mol] \n', ml, 'Name');
            fprintf('  %s\n', sep)
            for i = 1:nc
                pc = mixture.Pcrit(i);
                tc = mixture.Tcrit(i);
                vc = mixture.Vcrit(i);
                acf = mixture.acentricFactors(i);
                mw = mixture.molarMass(i);
                fprintf('  %*s | %1.2e | %3.1f K | %2.3e | %1.3f | %1.7f \n', ...
                        ml, cnames{i}, pc, tc, vc, acf, mw);
            end
            fprintf('  %s\n', sep);
            bi = mixture.getBinaryInteraction();
            fprintf('  ');
            if any(bi(:))
                fprintf('Binary interaction coefficients:\n');
                disp(bi);
            else
                fprintf('No non-zero binary interaction coefficients.\n')
            end
            fprintf('\n');
            known = {'Pcrit', 'Tcrit', 'Vcrit', 'acentricFactors', 'molarMass', 'names', 'name'};
            props = propertynames(mixture);
            extra = setdiff(props, known);
            if numel(extra)
                fprintf('\n  Additional properties:\n');
                for i = 1:numel(extra)
                    e = extra{i};
                    fprintf('  %s (%s)\n', e, class(mixture.(e)));
                end
                fprintf('\n');
            end
        end
    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
