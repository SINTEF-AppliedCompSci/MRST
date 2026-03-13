classdef ECPACompositionalMixture
    properties
        Tcrit % Critical temperature in Kelvin
        Pcrit % Critical pressure in Pascal
        Vcrit % Critical volume in m^3 / mol
        acentricFactors % Acentric factors (dimensionless)
        molarMass % Component mass (kg / mol)
        names % Names of each component. Each name must be unique.
        name % Name of the fluid mixture itself
        b % Co-volume parameter of the EoS (m^3 / mol)
        a0 % Energy parameter of the EoS (Pa * m^6 / mol)
        c1 % pure-compound parameter of the EoS
        epsilonAB % parameter for associating compounds
        betaAB % parameter for associating compounds
        alpha0  % molecular polarizability (C^2 * m^2 / J)
        mu0 % vacuum dipole moment (C * m)
        d % parameter of ions
        rborn % Born radius
        hyd % hydration
        charge % charge
        parachor % parachor
        c % Peneloux volume translation
    end
     
    properties ( Access = protected)
        bic % Binary interaction coefficients
        ap
    end
    
    methods
        function fluid = ECPACompositionalMixture(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass,...
                b, a0, c1, epsilonAB, betaAB, alpha0, mu0, d, rborn, hyd, charge, parachor, c, varargin)
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
            fluid.b = b;
            fluid.a0 = a0;
            fluid.c1 = c1;
            fluid.epsilonAB = epsilonAB;
            fluid.betaAB = betaAB;
            fluid.alpha0 = alpha0;
            fluid.mu0 = mu0;
            fluid.d = d;
            fluid.rborn = rborn;
            fluid.hyd = hyd;
            fluid.charge = charge;
            fluid.parachor = parachor;
            fluid.c = c;
            ncomp = numel(fluid.names);
            
            fluid.bic = zeros(ncomp, ncomp);
            fluid.ap = zeros(ncomp, ncomp);
        end
        
        function n = getNumberOfComponents(fluid)
            n = numel(fluid.names);
        end
        
        function n = getNumberOfIons(fluid)
            n = sum(~isnan(fluid.d));
        end
        
        function n = getNumberOfMolecules(fluid)
            n = sum(~isnan(fluid.Tcrit));
        end

        function bic = getBinaryInteraction(fluid)
            bic = fluid.bic;
        end

        function ap = getAssociationParameter(fluid)
            ap = fluid.ap;
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
        
        function fluid = setAssociationParameter(fluid, input)
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
            fluid.ap = input;
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
            sep = repmat('-', 198 + ml, 1);
            fprintf('  %*s | p_c [Pa] | T_c [K] | V_c [m^3] |  acf  | mw [kg/mol] | b [m^3/mol] | a0 [Pa*m^6/mol] |   c1   | epsilonAB | betaAB |   alpha0   |   mu0   |  d  [m]  | rborn [m] | hyd | charge | parachor | c \n', ml, 'Name');
            fprintf('  %s\n', sep)
            for i = 1:nc
                pc = mixture.Pcrit(i);
                tc = mixture.Tcrit(i);
                vc = mixture.Vcrit(i);
                acf = mixture.acentricFactors(i);
                mw = mixture.molarMass(i);
                bi = mixture.b(i);
                a0i = mixture.a0(i);
                c1i = mixture.c1(i);
                epsilonABi = mixture.epsilonAB(i);
                betaABi = mixture.betaAB(i);
                alpha0i = mixture.alpha0(i);
                mu0i = mixture.mu0(i);
                di = mixture.d(i);
                rborni = mixture.rborn(i);
                hydi = mixture.hyd(i);
                chargei = mixture.charge(i);
                parachori = mixture.parachor(i);
                ci = mixture.c(i);
                fprintf('  %*s | %8.2e | %5.1f K | %9.3e | %5.3f | %1.9f | %11.3e | %15.6f | %1.4f | %9.1f | %6.4f | %10.3e | %7.1e | %8.2e | %9.3e | %3.0f | %6.0f | %8.1f | %3.2f \n', ...
                        ml, cnames{i}, pc, tc, vc, acf, mw, bi, a0i, c1i, epsilonABi, ...
                        betaABi, alpha0i, mu0i, di, rborni, hydi, chargei, parachori, ci);
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
            acp = mixture.getAssociationParameter();
            fprintf('  ');
            if any(bi(:))
                fprintf('Association parameters:\n');
                disp(acp);
            else
                fprintf('No non-zero association parameters.\n')
            end
            fprintf('\n');
            known = {'Pcrit', 'Tcrit', 'Vcrit', 'acentricFactors', 'molarMass', 'names', 'name', 'b','a0','c1','epsilonAB','betaAB','alpha0','mu0','d','rborn','hyd','charge','parachor','c'};
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
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
