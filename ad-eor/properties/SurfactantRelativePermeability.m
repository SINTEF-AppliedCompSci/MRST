classdef SurfactantRelativePermeability < BaseRelativePermeability

    properties
        zeroSurf
        fullSurf
    end
    
        
    methods
        function prop = SurfactantRelativePermeability(model, satreg, surfreg, varargin)
            prop@BaseRelativePermeability(model, varargin{:});
            prop = prop.dependsOn('surfactant', 'state');
            prop = prop.dependsOn({'CapillaryNumber'});
            prop.zeroSurf = StateFunction(model, satreg);
            prop.fullSurf = StateFunction(model, surfreg);
        end
        
        function kr = evaluateOnDomain(prop, model, state)

            fluid = model.fluid;
            [cs, Nc] = model.getProps(state, 'surfactant', 'CapillaryNumber');

            isSft = value(cs) > 0;
            m = zeros(model.G.cells.num, 1);
            m = model.AutoDiffBackend.convertToAD(m, cs);
            if nnz(isSft) > 0
                logNc = log(Nc(isSft))/log(10);  % log10 is not yet implemented in ADI
                % We cap logNc (as done in Eclipse)
                logNc = min(max(-20, logNc), 20);
                m(isSft) = fluid.miscfact(logNc, 'cellInx', find(isSft));
            end
            
            if numel(state.s) > 2
                [sW, sO, sG] = model.getProps(state, 'sw', 'so', 'sg');
            else
                [sW, sO] = model.getProps(state, 'sw', 'so');
            end

            satreg  = model.rock.regions.saturation; 
            surfreg = model.rock.regions.surfactant;

            sWcon_noSft = fluid.krPts.w(satreg   , 2); % Residual water saturation   without surfactant
            sOWres_noSft = fluid.krPts.ow(satreg , 2); % Residual oil saturation     without surfactant
            sWcon_Sft   = fluid.krPts.w(surfreg  , 2); % Residual water saturation   with    surfactant
            sOWres_Sft   = fluid.krPts.ow(surfreg, 2); % Residual oil saturation     with    surfactant
            
            % Interpolated water/oil residual saturations
            sNcWcon = m.*sWcon_Sft + (1 - m).*sWcon_noSft;
            sNcOWres = m.*sOWres_Sft + (1 - m).*sOWres_noSft;

            sNcWEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOWres);
            sNcOWEff = (sO - sNcOWres)./(1 - sNcWcon - sNcOWres);

            % Rescaling of the saturation - without surfactant
            sW_noSft  = (1 - sWcon_noSft - sOWres_noSft).*sNcWEff + sWcon_noSft;
            sOW_noSft  = (1 - sWcon_noSft - sOWres_noSft).*sNcOWEff + sOWres_noSft;
            % Compute rel perm - without surfactant
            krW_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krW, sW_noSft);
            if isfield(fluid, 'krO')
                krO_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krO, ...
                                                                 sOW_noSft);
            else
                krOW_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krOW, ...
                                                                 sOW_noSft);
                krO_noSft = krOW_noSft;
            end
            
            % Rescaling of the saturation - with surfactant
            sW_Sft  = (1 - sWcon_Sft - sOWres_Sft).*sNcWEff + sWcon_Sft;
            sOW_Sft  = (1 - sWcon_Sft - sOWres_Sft).*sNcOWEff + sOWres_Sft;
            % Compute rel perm - with surfactant
            krW_Sft = prop.fullSurf.evaluateFunctionOnDomainWithArguments(fluid.krW, sW_Sft);
            if isfield(fluid, 'krO')
                krO_Sft = prop.fullSurf.evaluateFunctionOnDomainWithArguments(fluid.krO, ...
                                                               sOW_Sft);
            else
                krOW_Sft = prop.fullSurf.evaluateFunctionOnDomainWithArguments(fluid.krOW, ...
                                                               sOW_Sft);
                krO_Sft = krOW_Sft;
            end
            
            if model.gas
                sOGres_noSft   = fluid.krPts.og(surfreg, 2); % Residual oil saturation without surfactant
                sGres_noSft = fluid.krPts.g(satreg, 2);  % Residual gas saturation without surfactant
                sOGres_Sft   = fluid.krPts.og(surfreg, 2); % Residual oil saturation with surfactant
                sGres_Sft   = fluid.krPts.g(surfreg, 2); % Residual gas saturation with surfactant
                
                sNcOGres = m.*sOGres_Sft + (1 - m).*sOGres_noSft;
                sNcGres = m.*sGres_Sft + (1 - m).*sGres_noSft;
                
                sNcGEff = (sG - sNcGres)./(1 - sNcGres - sNcOGres);
                sNcOGEff = (sO - sNcOGres)./(1 - sNcGres - sNcOGres);
                
                sG_noSft  = (1 - sGres_noSft - sOGres_noSft).*sNcGEff + sGres_noSft;
                sOG_noSft  = (1 - sGres_noSft - sOGres_noSft).*sNcOGEff + sOGres_noSft;
                
                krOG_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krOG, ...
                                                              sOG_noSft);
                krG_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krG, ...
                                                              sG_noSft);
                                                          
                sWcon_noSft = min(sWcon_noSft, value(sW_noSft)-1e-5);
                d  = (sG_noSft - sGres_noSft + sW_noSft - sWcon_noSft);
                ww = (sW_noSft - sWcon_noSft)./d;
                wg = 1 - ww;
                krO_noSft = wg.*krOG_noSft + ww.*krOW_noSft;
                
                sOG_Sft  = (1 - sGres_Sft - sOGres_Sft).*sNcOGEff + sOGres_Sft;
                sG_Sft  = (1 - sGres_Sft - sOGres_Sft).*sNcGEff + sGres_Sft;
                
                krOG_Sft = prop.fullSurf.evaluateFunctionOnDomainWithArguments(fluid.krOG, ...
                                                              sOG_Sft);
                krG_Sft = prop.fullSurf.evaluateFunctionOnDomainWithArguments(fluid.krG, ...
                                                              sG_Sft);
                
                sWcon_Sft = min(sWcon_Sft, value(sW_Sft)-1e-5);
                d  = (sG_Sft - sGres_Sft + sW_Sft - sWcon_Sft);
                ww = (sW_Sft - sWcon_Sft)./d;
                wg = 1 - ww;
                krO_Sft = wg.*krOG_Sft + ww.*krOW_Sft;
                
                krG  = m.*krG_Sft + (1 - m).*krG_noSft;
            end

            % Interpolate relperm, with and without surfactant
            
            krW  = m.*krW_Sft + (1 - m).*krW_noSft;
            krO  = m.*krO_Sft + (1 - m).*krO_noSft;
            
            if model.gas
                kr = {krW, krO, krG};
            else
                kr = {krW, krO};
            end
            
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
