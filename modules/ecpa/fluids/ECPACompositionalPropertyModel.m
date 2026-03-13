classdef ECPACompositionalPropertyModel < PropertyModel
    % Property model for compositional models
    properties
        
    end
    
    methods
        function model = ECPACompositionalPropertyModel(varargin)
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
            R = 8.314462618;
            V = (R.*T.*Z)./p;
            c = model.fluid.c;
            if iscell(x)
                corr = 0;
                for i = 1:numel(x)
                    corr = corr + c(i).*x{i};
                end
            else
                corr = sum(bsxfun(@times, c, x), 2);
            end
            V = V - 1e-6 .* corr;
            rho = 1./V;
        end
        
        function mu = computeViscosity(model, eos, P, x, Z, T, isLiquid)
            if 0
                % Compute viscosity using free volume theory (FV viscosity)
                R = 8.3144598;
                Tc = model.fluid.Tcrit;
                Tc = Tc(~isnan(Tc));
                nmole = numel(Tc);
                Ace = model.fluid.acentricFactors;
                names = model.fluid.names;
                MW = model.fluid.molarMass .* 1e3;
                Vc = model.fluid.Vcrit .* 1e6;
                Vc = Vc(1:nmole);
                Tr = bsxfun(@rdivide, T, Tc);
                ncomp = numel(MW);
                
                T1 = 1.2593 .* Tr;
                q = 1.16145 ./ T1 .^ 0.14874 + 0.52487 ./ exp(0.7732 .* T1) ...
                    + 2.16178 ./ exp(2.43787 .* T1) - 6.435e-4 .* T1 .^ 0.14874 ...
                    .* sin(18.0323 .* T1 .^ (-0.7683) - 7.27371);
                Fc = 1 - 0.2756 .* Ace(1:nmole);
                if strcmp(names(1), {'Water'})
                    Fc(1) = 1.139494;
                end
                mu0 = 40.785 .* bsxfun(@times, sqrt(bsxfun(@times, MW(1:nmole), T)), Fc) ...
                    ./ bsxfun(@times, Vc.^(2/3), q);
                if iscell(x)
                    mu0_m = 0;
                    for i = 1:nmole
                        mu0_m = mu0_m + exp(x{i} .* log(mu0(:,i)));
                    end
                    M = 0;
                    for i = 1:ncomp
                        M = M + x{i}.*MW(i);
                    end
                    [L, alp, B] = model.getFVparam(x, Tr, names, nmole);
                    rho = model.computeMolarDensity(eos, P, x, Z, T, isLiquid);
                else
                    mu0_m = exp(sum(x(:,1:nmole) .* log(mu0), 2));
                    M = sum(bsxfun(@times, x, MW), 2);
                    [L, alp, B] = model.getFVparam(x, Tr, names, nmole);
                    rho = model.computeMolarDensity(eos, P, x, Z, T, isLiquid);
                end
                E = 1e-3 .* M .* alp .* rho  + P ./ rho;
                d_mu = 1e-6.*rho.*L.*E.* (1000.*M./(3.*R.*T)).^0.5 .*exp(B.*(E./(R.*T)).^1.5);
                mu = 1e-7 .* (mu0_m + d_mu); % Pa¡¤s
            else
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
        end
        
        function [Lm, alpm, Bm] = getFVparam(model, x, Tr, names, nmole)
            nc = numel(Tr(:,1));
            ncomp = numel(names);
            L = zeros(nc, ncomp);alp = zeros(nc, ncomp);B = zeros(nc, ncomp);
            for i = 1 : nmole
                switch(lower(names{i}))
                    case {'water', 'h2o'}
                        L(:,i) = 1.56866857e3.*Tr(:,i).^4 - 4.0307305e3.*Tr(:,i).^3 ...
                            + 3.86452737e3.*Tr(:,i).^2 - 1.63961155e3.*Tr(:,i) + 260.50228;
                        alp(:,i) = 2016.2.*Tr(:,i).^2 - 2746.*Tr(:,i) + 985.06;
                        B(:,i) = 0.0043.*Tr(:,i) - 0.0042;
                    case {'carbondioxide', 'co2'}
                        L(:,i) = 0.55;
                        alp(:,i) = 27.4;
                        B(:,i) = 0.0026;
                    case {'methane', 'ch4'}
                        L(:,i) = 0.76;
                        alp(:,i) = 29.7;
                        B(:,i) = 0.0114;
                    case {'hydrogensulfide', 'h2s'}
                        L(:,i) = 0.59;
                        alp(:,i) = 39.6;
                        B(:,i) = 0.0049;
                    case {'nitrogen', 'n2'}
                        L(:,i) = 0.65;
                        alp(:,i) = 11.8;
                        B(:,i) = 0.013;
                    case {'oxygen', 'o2'}
                        L(:,i) = 0.53;
                        alp(:,i) = 14.1;
                        B(:,i) = 0.0048;
                    case {'argon', 'ar'}
                        L(:,i) = 0.52;
                        alp(:,i) = 12.6;
                        B(:,i) = 0.0043;
                    case {'sulfurdioxide', 'so2'}
                        L(:,i) = 0.55;
                        alp(:,i) = 28.8;
                        B(:,i) = 0.0088;
                end
            end
            if iscell(x)
                Lm = 0; alpm = 0; Bm = 0;
                for i = 1 : ncomp
                    Lm = Lm + x{i} .* L(:,i);
                    alpm = alpm + x{i} .* alp(:,i);
                    Bm = Bm + x{i} .* B(:,i);
                end
            else
                Lm = sum(x .* L, 2);
                alpm = sum(x .* alp, 2);
                Bm = sum(x .* B, 2);
            end
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
