classdef DiagnosticsNPV < DiagnosticsObjective
    properties
        timeHorizon    = 10;
        timeUnit       = 'year';
        smoothing      = true;
        smoothFac      = 0.05;
        ro             = 50/stb;%45/stb;
        rwi            = -5/stb;%-2/stb;
        rwp            = -3/stb;%-6/stb;
        rgi            = 0;
        rgp            = 0;
        discount       = 0.1;
        phaseWeights    = [];
        lengthCost     = 0*10e3/meter;  
        scaleProduction   = 1;
        scaleBreakthrough = 1;
        pRef           = 200*barsa;
        %scaleTOFBreakThrough = [];
    end
    
    methods
        function obj = DiagnosticsNPV(model, varargin)
            assert(isa(model, 'DiagnosticsModel'), ...
                'Expected model of class DiagnosticsModel, got %s', class(model));
            obj = obj@DiagnosticsObjective(model);
            obj = merge_options(obj, varargin{:});
            obj.model = makeCompatibleForObjective(obj.model);
        end
        
        function varargout = evaluate(obj, state, D, W)
            computeDerivatives = nargout > 1;
            [fstate, varNames] = obj.model.getStateAD(state, 'initPressure', computeDerivatives, 'initTOFBackward', computeDerivatives);
            gboModel = obj.model.parentModel.parentModel;
            % function only of backward tof
            tof_b = obj.model.getProp(fstate, 'tof_backward');
            qs    = fstate.FacilityState.surfacePhaseRates;
            cqs   = gboModel.FacilityModel.getComponentSources(fstate);
           
            pv = gboModel.operators.pv;
            
            w = obj.phaseWeights;
            if isempty(w)
                w = obj.model.phaseWeights;
            elseif isa(w, 'function_handle')
                w = w(state.s);
            end
            if size(w,1) == 1
                w = repmat(w, [gboModel.G.cells.num, 1]);
            end
            
            isInj = vertcat(state.wellSol.sign) > 0;
            [injValue, prodValue, injProdValues] = deal(0);
            
            fluid = gboModel.fluid;
            
            wix = gboModel.getPhaseIndex('W');
            if ~isempty(wix)
                injWater  = sign( qs{wix} ) == 1;
                injValue  = injValue  + obj.rwi*sum( qs{wix}.*injWater);
                prodValue = prodValue + fluid.bW(obj.pRef)*obj.rwp*w(:,wix);
                injWaterCon = sign(cqs.value{wix}) == 1;
                rhow = gboModel.fluid.rhoWS;
                injProdValues = injProdValues + obj.rwp*cqs.value{wix}.*injWaterCon/rhow;
            end
            
            oix = gboModel.getPhaseIndex('O');
            if ~isempty(oix)
                injOil     = sign( qs{oix} ) == 1;
                if any(injOil)
                    warning('Injecting oil!');
                    injValue  = injValue  - obj.ro*sum( qs{oix}.*injOil);
                end
                
                prodValue = prodValue + fluid.bO(obj.pRef)*obj.ro*w(:,oix);
            end
            
            gix = gboModel.getPhaseIndex('G');
            if ~isempty(gix)
                injValue  = injValue  + obj.rgi*sum( qs{gix}.*isInj);
                prodValue = prodValue + obj.rgp*w(:,gix);
                rhog = gboModel.fluid.rhoGS;
                injGasCon = sign(cqs.value{gix}) == 1;
                injProdValues = injProdValues + obj.rgp*cqs.value{gix}.*injGasCon/rhog;
            end
            
            dfac    = (1+obj.discount).^(-tof_b/year);
            
            % handle cutoff
            tunit = obj.timeUnit;
            if ischar(tunit)
                if strcmpi(tunit, 'pvi')
                    tunit = sum(pv)./sum( qs{1}.*isInj);
                else
                    tunit = eval(tunit);
                end
            end    
            
            for nt = 1:numel(obj.timeHorizon)
                T = convertFrom(obj.timeHorizon(nt), tunit);
                if ~obj.smoothing
                    cutoff = tof_b <= T;
                else
                    dT  = T -tof_b;
                    sf  = obj.smoothFac*T;
                    cutoff = .5*dT.*(dT.^2+sf^2).^(-.5)+.5;
                end
                
                % injection value must be multiplied with integral (1+d)^(-tof) from 0 to T:
                if obj.discount > 0
                    rr = log(1+obj.discount)/tunit;
                    injFac = (1-exp(-rr*T))/rr;
                else
                    injFac = T;
                end
                % account for production of injeceted water/gas
                tof_i = tof_b(cqs.cells)/obj.scaleBreakthrough;
                %cut_i = cutoff(cqs.cells);
                cut_i = tof_i <= T;
                if obj.discount > 0
                    rr = log(1+obj.discount)/tunit;
                    injProdVal = sum( (exp(-rr*tof_i) - exp(-rr*T)).*cut_i.*injProdValues )/rr;
                else
                    injProdVal = sum((T-tof_i).*cut_i.*injProdValues);
                end
                
                v = injFac*injValue + sum( (prodValue.*pv.*cutoff).*dfac )/obj.scaleProduction + injProdVal;
                dvdpos = cell(1, numel(W));
                if obj.lengthCost ~= 0 && isfield(W(1), 'posControl')
                    % loop over position controls (only these are added to v)
                    pcontrols = {W.posControl};
                    wp = find(~cellfun(@isempty, pcontrols));
                    
                    for k = 1:numel(wp)
                        pc  = pcontrols{wp(k)};
                        tr  = pc.getTrajectory();
                        dtr = diff(tr);
                        v = v - obj.lengthCost*sum(sqrt( dot(dtr, dtr, 2) ));
                        % partials wrt perturbation directions
                        dvdpos{wp(k)} = -obj.lengthCost*pc.parameters.lengthDerivative;
                    end
                end
                if numel(obj.timeHorizon) > 1
                    vals(nt) = v;
                end
            end
            if numel(obj.timeHorizon) > 1
                varargout{1} = vals;
                return;
            end
            
            if ~computeDerivatives
                varargout{1} = v;
            else
                varargout{1} = v.val;
                % dFdx
                tmp = struct('tof', [], 'tracerIx', [], 'tracer', []);
                partials = struct('pressure', [], 'forward', [], 'backward', tmp);
                partials.pressure = v.jac(1:(end-1));
                partials.backward.tof = v.jac{end};
                varargout{2} = partials;
                varargout{3} = struct('well', [], 'position', {dvdpos});
            end
        end
    end
end


function dm = makeCompatibleForObjective(dm)
dm.computeForward  = false;
dm.computeBackward = true;
dm.tracerWells     = 'none';
dm.computeWellTOFs = false;
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
