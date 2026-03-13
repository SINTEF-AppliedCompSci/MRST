classdef ECPATableCompositionalMixture < ECPACompositionalMixture
    % Create MRST compositional mixture from stored tables
    properties
        
    end
    
    methods
        function mix = ECPATableCompositionalMixture(names, outnames)
            % Create fluid from names
            %
            % SYNOPSIS:
            %   f = ECPATableCompositionalMixture({'Water','CarbonDioxide','Na+','Cl-'});
            %
            % PARAMETERS:
            %   names    - Cell array of valid names. See `getFluidList` for
            %              valid names.
            %   outnames - The names the components will have in output
            %              (optional).
            % RETURNS:
            %   mix - Initialized compositional mixture.
            if ischar(names)
                names = {names};
            end
            if nargin == 1
                outnames = names;
            else
                if ischar(outnames)
                    outnames = {outnames};
                end
                assert(numel(names) == numel(outnames))
            end
            
            ncomp = numel(names);
            [Tcrit, Pcrit, rhocrit, acc, molarMass, b, a0, c1, epsilonAB,...
                betaAB, alpha0, mu0, d, rborn, hyd, charge, parachor, c] ...
                = deal(zeros(1, ncomp));
                        
            fluids = eCPAPropFluidsStructs();
            validChoices = ECPATableCompositionalMixture.getFluidList();
            
            ok = ismember(lower(names), lower(validChoices));
            
            if ~all(ok)
                s ='Unable to create fluid. The following names were not known to CoolProps table: ';
                msg = [s, sprintf('%s ', names{~ok})];
                error(msg);
            end
            
            for i = 1:numel(names)
                isF = strcmpi(names{i}, validChoices);
                str = fluids(isF);
                Tcrit(i) = str.Tcrit;
                Pcrit(i) = str.Pcrit;
                rhocrit(i) = str.rhocrit;
                acc(i) = str.acentric;
                molarMass(i) = str.molarmass;
                b(i) = str.b;
                a0(i) = str.a0;
                c1(i) = str.c1;
                epsilonAB(i) = str.epsilonAB;
                betaAB(i) = str.betaAB;
                alpha0(i) = str.alpha0;
                mu0(i) = str.mu0;
                d(i) = str.d;
                rborn(i) = str.rborn;
                hyd(i) = str.hyd;
                charge(i) = str.charge;
                parachor(i) = str.parachor;
                c(i) = str.c;
            end
            Vcrit = molarMass./rhocrit;
            mix = mix@ECPACompositionalMixture(outnames, Tcrit, Pcrit, ...
                Vcrit, acc, molarMass, b, a0, c1, epsilonAB, betaAB, ...
                alpha0, mu0, d, rborn, hyd, charge, parachor, c);
            mix.name = 'eCPA - Tabulated';
        end
    end
    methods (Static)
        function varargout = getFluidList()
            fluids = eCPAPropFluidsStructs();
            names = {fluids.name};
            if nargout > 0
                varargout{1} = names;
            else
                disp('Possible fluid choices are:');
                fprintf('%s\n', names{:});
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
