classdef BaseRelativePermeability < StateFunction & SaturationProperty
    % Phase relative permeability (per cell)
    properties
        relpermPoints = 2;
        immobileChop = false;
        regions_imbibition = [];
    end

    methods
        function gp = BaseRelativePermeability(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'s'}, 'state');
            gp.label = 'k_\alpha';
            if isfield(model.fluid, 'ehystr')
                warning('Hysteresis is not supported by %s', class(gp));
            end
        end

        function kr = evaluateOnDomain(prop, model, state)
            nph = model.water + model.oil + model.gas;
            if nph < 2
                kr = {ones(model.G.cells.num, 1)};
            else
                if model.oil
                    if model.gas && model.water
                        kr = prop.relPermWOG(model, state);
                    elseif model.gas
                        kr = prop.relPermOG(model, state);
                    elseif model.water
                        kr = prop.relPermWO(model, state);
                    end
                else
                    kr = prop.relPermUnified(model, state);
                end
            end
        end

        % --------------------------------------------------------------------%
        function kr = relPermWOG(prop, model, state)
            % Three-phase system
            [sw, so, sg] = model.getProps(state, 'sw', 'so', 'sg');
            krW = prop.evaluatePhaseRelativePermeability(model, 'w', sw);
            krG = prop.evaluatePhaseRelativePermeability(model, 'g', sg);
            % Oil rel perm is special
            if isfield(model.fluid, 'krO')
                krO = prop.evaluateFluid(model, 'krO', so);
            else
                swcon = prop.getConnateWater(model, state);
                swcon = min(swcon, value(sw)-1e-5);
                d  = (sg+sw-swcon);
                ww = (sw-swcon)./d;
                krow = prop.evaluatePhaseRelativePermeability(model, 'ow', so);
                krog = prop.evaluatePhaseRelativePermeability(model, 'og', so);
                krO  = (1-ww).*krog + ww.*krow;
            end
            kr = {krW, krO, krG};
        end

        function kr = relPermWO(prop, model, state)
            % Water-oil system
            [sw, so] = model.getProps(state, 'sw', 'so');
            f = model.fluid;
            krW = evaluatePhaseRelativePermeability(prop, model, 'w', sw);
            if isfield(f, 'krO')
                krO = prop.evaluateFluid(model, 'krO', so);
            else
                krO = prop.evaluatePhaseRelativePermeability(model, 'ow', so);
            end
            kr = {krW, krO};
        end

        function kr = relPermOG(prop, model, state)
            % Oil-gas system
            [sg, so] = model.getProps(state, 'sg', 'so');
            f = model.fluid;
            krG = evaluatePhaseRelativePermeability(prop, model, 'g', sg);
            if isfield(f, 'krO')
                krO = prop.evaluateFluid(model, 'krO', so);
            else
                krO = prop.evaluatePhaseRelativePermeability(model, 'og', so);
            end
            kr = {krO, krG};
        end

        function kr = relPermUnified(prop, model, state)
            % Rel.perm without special oil treatment
            phases = model.getPhaseNames();
            snames = arrayfun(@(x) ['s', x], phases, 'UniformOutput', false);
            nph = numel(phases);
            s = cell(1, nph);
            kr = cell(1, nph);
            [s{:}] = model.getProps(state, snames{:});
            for i = 1:nph
                kr{i} = prop.evaluatePhaseRelativePermeability(model, 'g', s{i});
            end
        end

        function kr = evaluatePhaseRelativePermeability(prop, model, phase, s, cells)
            if nargin < 5
                cells = ':';
            end
            fn = model.fluid.(['kr', upper(phase)]);
            f = model.fluid;
            if prop.scalingActive
                reg = prop.regions;
                if ~isempty(reg)
                    reg = reg(cells);
                end
                [ss, kr_max_m] = prop.scaleSaturation(model.rock.krscale.drainage, phase, reg, f, s, cells);
                kr = kr_max_m.*prop.evaluateFunctionCellSubset(fn, cells, ss);
            else
                kr = prop.evaluateFunctionCellSubset(fn, cells, s);
            end
        end

        function [s_scale, k_max_m] = scaleSaturation(prop, pts, phase, reg, f, s, cells)
            if nargin  < 7
                cells = ':';
            end
            if ischar(cells)
                cells = prop.cell_subset;
            else
                assert(isnumeric(cells));
                cells = cells(prop.cell_subset);
            end
            if prop.relpermPoints == 2 % 2-point
                [m, c, p, k] = prop.getTwoPointScalers(pts, phase, reg, f, cells);
                ix1 = s < p{1};
                ix2 = s >= p{2};
                ix  = ~or(ix1,ix2);
                s_scale = (ix.*m).*s + (ix.*c + ix2);
            elseif prop.relpermPoints == 3
                [m, c, p, k] = prop.getThreePointScalers(pts, phase, reg, f, cells);
                ix1 = and(s >= p{1}, s < p{2});
                ix2 = and(s >= p{2}, s < p{3});
                ix3 = s >= p{3};
                a = ix1.*m{1} + ix2.*m{2};
                s_scale = a.*s + ix1.*c{1} + ix2.*c{2} + ix3;
            else
                error('Unknown number of scaling points');
            end
            k_max_m = k{2}./k{1};
        end

        function [state, chopped] = applyImmobileChop(prop, model, state, state0)
            chopped = false;
            if prop.immobileChop
                s_min = prop.getCriticalPhaseSaturations(model, state);
                s = state.s;
                s0 = state0.s;
                chopped = false;
                for ph = 1:size(s, 2)
                    sm = s_min(:, ph);
                    tol = 1e-8;

                    toMobile = s0(:, ph) < sm + tol & s(:, ph) > sm - tol;
                    toImmobile = s0(:, ph) > sm - tol & s(:, ph) < sm + tol;
                    s(toMobile, ph) = tol;
                    if size(sm, 1) > 1
                        s(toImmobile, ph) = sm(toIoImmobile) - 2*tol;
                        s(toMobile, ph) = sm(toMobile) + 2*tol;
                    else
                        s(toImmobile, ph) = sm - 2*tol;
                        s(toMobile, ph) = sm + 2*tol;
                    end
                    chopped = chopped | (toImmobile | toMobile);
                end
                state.s = s;
            end
        end

        function property = subset(property, subs)
            property = subset@StateFunction(property, subs);
            property.cell_subset = subs;
        end
    end
end



%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
