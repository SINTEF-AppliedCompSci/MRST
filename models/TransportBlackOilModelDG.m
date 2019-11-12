classdef TransportBlackOilModelDG < TransportModelDG
    % Two phase oil/water system without dissolution with discontinuous
    % Galerking discretization
    
    properties
    end

    methods
        % ----------------------------------------------------------------%
        function model = TransportBlackOilModelDG(G, rock, fluid, varargin)
            model = model@TransportModelDG(G, rock, fluid);
        end

        function state = initStateAD(model, state, vars, names, origin)
            removed = false(size(vars));
            parent = model.parentModel;
            if parent.disgas || parent.vapoil
                % Black-oil specific variable switching
                if parent.water
                    isw   = strcmpi(names, 'swdof');
                    sWdof = vars{isw};
                    sW    = model.disc.getCellMean(state, sWdof);
                    removed = removed | isw;
                else
                    sW = 0;
                end

                isx = strcmpi(names, 'x');
                x = vars{isx};
                sG = model.getProps(state, 'sg');
                st  = model.getCellStatusVO(state, 1-sW-sG, sW, sG);
                sG = st{2}.*(1-sW) + st{3}.*x;
                sO = st{1}.*(1-sW) + ~st{1}.*(1 - sW - sG);
                if model.water
                    sat = {sW, sO, sG};
                else
                    sat = {sO, sG};
                end
                removed(isx) = true;
            else
                % Without variable switching
                phases = parent.getPhaseNames();
                nph = numel(phases);
                sat = cell(1, nph);
                fill = model.disc.getFillSat(state);
                removed_sat = false(1, nph);
                for i = 1:numel(phases)
                    sub = strcmpi(names, ['s', phases(i), 'dof']);
                    if any(sub)
                        fill = fill - vars{sub};
                        removed = removed | sub;
                        removed_sat(i) = true;
                        sat{i} = vars{sub};
                    end
                end
                if any(~removed_sat)
                    sat{~removed_sat} = fill;
                end
            end
            state = model.setProp(state, 'sdof', sat);
            state = model.initStateFunctionContainers(state);
            
            if not(isempty(parent.FacilityModel))
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(parent.FacilityModel);
                isF = strcmp(origin, fm);
                state = parent.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end

            % Set up state with remaining variables
%             state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
            % Account for dissolution changing variables
            if parent.disgas
                rsSat = model.getProp(state, 'RsMax');
                rs = ~st{1}.*rsSat + st{1}.*x;
                % rs = rs.*(value(sO) > 0);
                state = model.setProp(state, 'rs', rs);
            end

            if parent.vapoil
                rvSat = model.getProp(state, 'RvMax');
                rv = ~st{2}.*rvSat + st{2}.*x;
                % rv = rv.*(value(sG) > 0);
                state = model.setProp(state, 'rv', rv);
                % No rv, no so -> zero on diagonal in matrix
                bad_oil = value(sO) == 0 & value(rv) == 0;
                if any(bad_oil)
                    sO(bad_oil) = 1 - sW(bad_oil) - value(sG(bad_oil));
                    state = model.setProp(state, 'sO', sO);
                end
            end
        end
        
    end
    
end
        
%{
Copyright 2009-2019 SINTEF ICT, Applied Mathematics.

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