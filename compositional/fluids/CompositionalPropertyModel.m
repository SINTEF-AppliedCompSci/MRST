classdef CompositionalPropertyModel < PropertyModel
    % Property model for compositional models
    properties
    end
    
    methods
        function model = CompositionalPropertyModel(varargin)
            model = model@PropertyModel(varargin{:});
        end
        function rho = computeDensity(model, p, x, Z, T, isLiquid)
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
            R = 8.3144598;
            rho = p.*M./(R.*T.*Z);
        end
        
        function rho = computeMolarDensity(model, p, x, Z, T, isLiquid)
            % Predict molar density from EOS Z-factor
            R = 8.3144598;
            rho = p./(R.*T.*Z);
        end
        
        function mu = computeViscosity(model, P, x, Z, T, isLiquid)
            % Compute viscosity using the Lohrenz, Bray and Clark
            % correlation for hydrocarbon mixtures (LBC viscosity)
            ncomp = numel(x);
            molfactor = 1/gram;
            
            rho = model.computeMolarDensity(P, x, Z, T, isLiquid);
            
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
            mu = mu_atm + ((0.1023 + 0.023364.*rhor + 0.058533.*rhor.^2 ...
                - 0.040758.*rhor.^3 + 0.0093324.*rhor.^4).^4 - 1e-4)./e_mix;
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