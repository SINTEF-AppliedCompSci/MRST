classdef SurfactantRelativePermeability < BaseRelativePermeability
% The surfactant relative permeability model corresponds to a scaling of the
% permeability in both the saturation and relative permeability axis. The
% scaling with respect to saturation enables in particular to change the
% value of residual saturations.
    
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
            assert(all(isfield(model.fluid,{'miscfact','krPts'})));
            assert(any(isfield(model.fluid,{'krO','krOW'})));
            assert(isfield(model.fluid,'krW'));
        end
        
        function kr = evaluateOnDomain(prop, model, state)
            fluid      = model.fluid;
            immiscEval = @(a,b) prop.zeroSurf.evaluateFunctionOnDomainWithArguments(a, b);
            miscEval   = @(a,b) prop.fullSurf.evaluateFunctionOnDomainWithArguments(a, b);
            satreg     = model.rock.regions.saturation; 
            surfreg    = model.rock.regions.surfactant;
            
            [cs, Nc] = model.getProps(state, 'surfactant', 'CapillaryNumber');

            % Compute interpolation paramter
            m = zeros(model.G.cells.num, 1);
            m = model.AutoDiffBackend.convertToAD(m, cs);
            if nnz(value(cs) > 0) > 0
                logNc = log(Nc)/log(10);          % ADI does not implement log10
                logNc = min(max(-20, logNc), 20); % We cap logNc as in ECLIPSE
                m     = fluid.miscfact(logNc);
            end
            
            if model.gas
                [sW, sO, sG] = model.getProps(state, 'sw', 'so', 'sg');
            else
                [sW, sO] = model.getProps(state, 'sw', 'so');
            end

            % Residual saturations
            sWc_ns  = fluid.krPts.w(satreg  , 2); % Water without surfactant
            sWc_s   = fluid.krPts.w(surfreg , 2); % Water with    surfactant
            sOWr_ns = fluid.krPts.ow(satreg , 2); % Oil   without surfactant
            sOWr_s  = fluid.krPts.ow(surfreg, 2); % Oil   with    surfactant
            
            % Interpolated water/oil residual saturations
            sNcWc    = m.*sWc_s + (1 - m).*sWc_ns;
            sNcOWr   = m.*sOWr_s + (1 - m).*sOWr_ns;
            sNcWeff  = (sW - sNcWc)./(1 - sNcWc - sNcOWr);
            sNcOWeff = (sO - sNcOWr)./(1 - sNcWc - sNcOWr);

            % Rescaling of the saturation - without surfactant
            sW_ns   = (1 - sWc_ns - sOWr_ns).*sNcWeff + sWc_ns;
            sOW_ns  = (1 - sWc_ns - sOWr_ns).*sNcOWeff + sOWr_ns;
            % Compute rel perm - without surfactant
            krW_ns = immiscEval(fluid.krW, sW_ns);
            if isfield(fluid, 'krO')
                krO_ns = immiscEval(fluid.krO, sOW_ns);
            else
                krOW_ns = immiscEval(fluid.krOW, sOW_ns);
                krO_ns  = krOW_ns;
            end
            
            % Rescaling of the saturation - with surfactant
            sW_s  = (1 - sWc_s - sOWr_s).*sNcWeff + sWc_s;
            sOW_s = (1 - sWc_s - sOWr_s).*sNcOWeff + sOWr_s;
            % Compute rel perm - with surfactant
            krW_s = miscEval(fluid.krW, sW_s);
            if isfield(fluid, 'krO')
                krO_s = miscEval(fluid.krO, sOW_s);
            else
                krOW_s = miscEval(fluid.krOW, sOW_s);
                krO_s  = krOW_s;
            end
            
            if model.gas
                % Residual saturations
                sOGr_ns = fluid.krPts.og(surfreg, 2); % Oil without surfactant
                sOGr_s  = fluid.krPts.og(surfreg, 2); % Oil with    surfactant
                sGr_ns  = fluid.krPts.g(satreg,   2); % Gas without surfactant
                sGr_s   = fluid.krPts.g(surfreg,  2); % Gas with    surfactant
                
                sNcOGr   = m.*sOGr_s + (1 - m).*sOGr_ns;
                sNcGr    = m.*sGr_s + (1 - m).*sGr_ns;
                sNcGeff  = (sG - sNcGr)./(1 - sNcGr - sNcOGr);
                sNcOGeff = (sO - sNcOGr)./(1 - sNcGr - sNcOGr);
                
                sG_ns  = (1 - sGr_ns - sOGr_ns).*sNcGeff + sGr_ns;
                sOG_ns = (1 - sGr_ns - sOGr_ns).*sNcOGeff + sOGr_ns;
                
                krOG_ns = immiscEval(fluid.krOG, sOG_ns);
                krG_ns  = immiscEval(fluid.krG,   sG_ns);
                                                          
                sWc_ns = min(sWc_ns, value(sW_ns)-1e-5);
                d      = (sG_ns - sGr_ns + sW_ns - sWc_ns);
                ww     = (sW_ns - sWc_ns)./d;
                wg     = 1 - ww;
                krO_ns = wg.*krOG_ns + ww.*krOW_ns;
                sOG_s  = (1 - sGr_s - sOGr_s).*sNcOGeff + sOGr_s;
                sG_s   = (1 - sGr_s - sOGr_s).*sNcGeff + sGr_s;
                
                krOG_s = miscEval(fluid.krOG, sOG_s);
                krG_s  = miscEval(fluid.krG,   sG_s);
                
                sWc_s = min(sWc_s, value(sW_s)-1e-5);
                d     = (sG_s - sGr_s + sW_s - sWc_s);
                ww    = (sW_s - sWc_s)./d;
                wg    = 1 - ww;
                krO_s = wg.*krOG_s + ww.*krOW_s;
                krG   = m.*krG_s + (1 - m).*krG_ns;
            end

            % Interpolate relperm, with and without surfactant and return
            % the result
            krW  = m.*krW_s + (1 - m).*krW_ns;
            krO  = m.*krO_s + (1 - m).*krO_ns;
            if model.gas
                kr = {krW, krO, krG};
            else
                kr = {krW, krO};
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
