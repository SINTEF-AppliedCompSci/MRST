classdef WellboreFrictionLoss < StateFunction
% State function for wellbore friction loss in WellboreMode

    properties
        assumeTurbulent = false; % Assume turbulent flow only
    end
    
    methods
        function dp = evaluateOnDomain(prop, model, state)

            % Get dependencies
            v = model.getProp(state, 'massFlux');
            [s, rho, mu, flag] = model.parentModel.getProps(state, ...
                's'              , ...
                'Density'        , ...
                'Viscosity'      , ...
                'PhaseUpwindFlag'  ...
            );
                                    
            % Compute wellbore mixture density/viscosity
            [rhoMix, muMix] = deal(0);
            if ~iscell(s), s = expandMatrixToCell(s); end
            nph = model.getNumberOfPhases();
            for ph = 1:nph
                rhoMix = rhoMix + rho{ph}.*s{ph};
                muMix  = muMix  + mu{ph}.*s{ph};
            end
            % Face average density and and upstream-weight viscosity
            rhoMix = model.parentModel.operators.faceAvg(rhoMix);
            muMix  = model.parentModel.operators.faceUpstr(flag{1,1}, muMix);
            % Get wellbore segment inner and outer diameters
            [Di, Do] = deal(0, model.G.faces.radius.*2);
            if size(Do, 2) == 2
                Di = Do(:,1);
                Do = Do(:,2);
            end
            % Get wellbore segment lengths
            L = model.G.faces.length;
            % Get wellbore segment roughness
            [~, roughness] = model.getSegmentRoughness();
            
            % Convert mass flux to velocity
            v = v./(pi*rhoMix.*((Do/2).^2 - (Di/2).^2));
            % Compute Reynold's number
            Re = abs(rhoMix.*v.*(Do-Di)./muMix);
            % Fanning friction factor
            f = (-3.6*log(6.9./Re+(roughness./(3.7*(Do))).^(10/9))/log(10)).^(-2);

            [lam, int] = deal(false);
            if ~prop.assumeTurbulent
                % Divide into laminar, intermediate and tubulent regions
                [Re1, Re2] = deal(2000, 4000);
                lam = Re <= Re1;
                tur = Re >= Re2;
                int = ~or(lam, tur);
                f1 = 16/Re1;
                f2 = (-3.6*log10(6.9./Re2+(roughness./(3.7*(Do))).^(10/9))).^(-2);
                
            end

            if any(int)
                % Intermediate flow regime - interpolate between f1 and f2
                if numel(f2) > nnz(int), f2 = f2(int); end
                f(int) = f1 + ((f2-f1)/(Re2-Re1)).*(Re(int)-Re1);
            end
            
            % Compute friction loss
            dp = -(2*sign(value(v)).*L./(Do-Di)).*(f.*rhoMix.*v.^2);
            
            % Explicitly calculate laminar flow friction loss to ensure
            % correct derivatives at zero rate
            dpLam   = -32*muMix.*(L./(Do - Di).^2).*v;
            dp(lam) = dpLam(lam);
           
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