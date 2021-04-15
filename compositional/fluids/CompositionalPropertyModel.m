classdef CompositionalPropertyModel < PropertyModel
    % Property model for compositional models
    properties
        volumeShift
    end
    
    methods
        function model = CompositionalPropertyModel(varargin)
            model = model@PropertyModel(varargin{:});
        end
        function rho = computeDensity(model, EOS, p, x, Z, T, isLiquid)
            % Predict mass density from EOS Z-factor
            ncomp = numel(model.fluid.names);
            M = 0;
            if iscell(x)
                for i = 1:ncomp
                    M = M + x{i}.*model.fluid.molarMass(i);
                end
            else
                M = sum(bsxfun(@times, x, model.fluid.molarMass), 2);
            end
            rho = model.computeMolarDensity(EOS, p, x, Z, T, isLiquid);
            rho = rho.*M;
        end
        
        function rho = computeMolarDensity(model, EOS, p, x, Z, T, isLiquid)
            % Predict molar density from EOS Z-factor
            R = 8.3144598;
            V = (R.*T.*Z)./p;
            if ~isempty(model.volumeShift)
                shift = reshape(model.volumeShift, 1, []);
                Tc = model.fluid.Tcrit;
                Pc = model.fluid.Pcrit;
                
                beta = EOS.omegaB.*R.*Tc./Pc;
                factors = shift.*beta;
                if iscell(x)
                    corr = 0;
                    for i = 1:numel(x)
                        corr = corr + factors(i).*x{i};
                    end
                else
                    corr = sum(bsxfun(@times, factors, x), 2);
                end
                V = V - corr;
            end
            rho = 1./V;
        end
        
        function mu = computeViscosity(model, eos, P, x, Z, T, isLiquid)
            % Compute viscosity using the Lohrenz, Bray and Clark
            % correlation for hydrocarbon mixtures (LBC viscosity)
            if ~iscell(x)
                x = expandMatrixToCell(x);
            end
            ncomp = numel(x);
            molfactor = 1/gram;
            
            rho = model.computeMolarDensity(eos, P, x, Z, T, isLiquid);
            
            % We first compute an estimated low pressure (surface)
            % viscosity for the mixture
            MW = (molfactor*model.fluid.molarMass).^(1/2);
            [a, b] = deal(0);
            for i = 1:ncomp
                tr = T./model.fluid.Tcrit(i);
                Tc = model.fluid.Tcrit(i)./Rankine();
                Pc = model.fluid.Pcrit(i)./psia();

                mwi = MW(i);
                e_i = (5.4402*Tc.^(1/6))./(mwi.*Pc.^(2/3).*(centi*poise));
                
                large = double(tr) > 1.5;
                % Different estimates based on how far above we are from
                % the critical temp
                mu_i = (~large.*34e-5.*tr.^(0.94) + large.*17.78e-5.*(4.58*tr - 1.67).^0.625)./e_i;
                a = a + x{i}.*mu_i.*mwi;
                b = b + x{i}.*mwi;
            end
            % Final atmospheric / low pressure viscosity is the ratio of a and b
            mu_atm = a./b;
            % Compute critical properties and coefficient for final
            % expression
            [P_pc, T_pc, Vc, mwc] = model.computePseudoCriticalPhaseProperties(x);            
            e_mix = 5.4402*(T_pc./Rankine()).^(1/6)./((molfactor*mwc).^(1/2).*(P_pc./psia()).^(2/3).*(centi*poise));
            % Reduced density via definition of critical density 
            rhor = Vc.*rho;
            % Final adjusted viscosity at current conditions
            
            % From Jossi et al
            % coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093724, 1e-4];
            % From LBC paper
            coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093324, -1e-4];
            mu = mu_atm + ((coeffs(1) + coeffs(2).*rhor + coeffs(3).*rhor.^2 ...
                + coeffs(4).*rhor.^3 + coeffs(5).*rhor.^4).^4 + coeffs(6))./e_mix;
        end

        function [P_pc, T_pc, Vc, mw] = computePseudoCriticalPhaseProperties(model, x)
            % Molar fraction weighted aggregated properties to get
            % pseudocritical properties of mixture.
            f = model.fluid;
            ncomp = numel(x);
            [T_pc, P_pc, Vc, mw] = deal(0);
            for i = 1:ncomp
                T_pc = T_pc + f.Tcrit(i)*x{i};
                P_pc = P_pc + f.Pcrit(i)*x{i};
                mw = mw + f.molarMass(i)*x{i};
                Vc = Vc + f.Vcrit(i)*x{i};
            end
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
