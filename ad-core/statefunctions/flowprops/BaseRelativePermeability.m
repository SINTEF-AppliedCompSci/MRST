classdef BaseRelativePermeability < StateFunction
    properties
        relpermScaling = false;
        relpermPoints = 2;
        immobileChop = false;
        regions_imbibition = [];
    end
    
    properties (Access = private)
        cell_subset = ':';
    end
    
    methods
        function gp = BaseRelativePermeability(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'s', 'sMax'}, 'state');
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
        f = model.fluid;
        swcon = 0;
        if prop.relpermScaling && ...
                isfield(model.rock, 'krscale') && ...
                isfield(model.rock.krscale.drainage, 'w')
            % Connate water in rock, endpoint scaling
            cix = prop.cell_subset;
            swcon = model.rock.krscale.drainage.w(cix, 1);
            % check for defaulted (nan) swcon -> use table values
            nix = isnan(swcon);
            if any(nix)
                swcon(nix) = reshape(f.krPts.w(prop.regions(nix), 1), [], 1);
            end
        elseif isfield(f, 'krPts')
            % Connate water from rel perm table
            swcon = reshape(f.krPts.w(prop.regions, 1), [], 1);
        end
        swcon = min(swcon, value(sw)-1e-5);

        d  = (sg+sw-swcon);
        ww = (sw-swcon)./d;
        
        krW = prop.evaluatePhaseRelativePermeability(model, 'w', sw);
        krG = prop.evaluatePhaseRelativePermeability(model, 'g', sg);

        wg = 1-ww;
        if isfield(f, 'krO')
            krO = prop.evaluateFunctionOnDomainWithArguments(f.krO, so);
        else
            krow = prop.evaluatePhaseRelativePermeability(model, 'ow', so);
            krog = prop.evaluatePhaseRelativePermeability(model, 'og', so);
            krO  = wg.*krog + ww.*krow;
        end
        kr = {krW, krO, krG};
    end
    
    function kr = relPermWO(prop, model, state)
        % Water-oil system
        [sw, so] = model.getProps(state, 'sw', 'so');
        f = model.fluid;
        krW = evaluatePhaseRelativePermeability(prop, model, 'w', sw);
        if isfield(f, 'krO')
            krO = prop.evaluateFunctionOnDomainWithArguments(f.krO, so);
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
            krO = prop.evaluateFunctionOnDomainWithArguments(f.krO, so);
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
        if prop.relpermScaling
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
            [m, c, p, k] = getTwoPointScalers(pts, phase, reg, f, cells);
            ix1 = s < p{1};
            ix2 = s >= p{2};
            ix  = ~or(ix1,ix2);
            s_scale = (ix.*m).*s + (ix.*c + ix2);
        elseif prop.relpermPoints == 3
            [m, c, p, k] = getThreePointScalers(pts, phase, reg, f, cells);
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
    
    function s_min = getCriticalPhaseSaturations(prop, model, state, drainage)
        if nargin < 4
            drainage = true;
        end
        
        if drainage
            satnum = prop.regions_imbibition;
        else
            satnum = prop.regions;
        end
        
        if isempty(satnum)
            satnum = 1;
        end
        
        phases = model.getPhaseNames();
        nph = numel(phases);
        s_min = zeros(size(satnum));
        if isfield(model.fluid, 'krPts')
            pts = model.fluid.krPts;
            for i = 1:nph
                ph = lower(phases(i));
                if strcmp(ph, 'o') && ~isfield(pts, 'o')
                    if model.water && model.gas
                        s_min(i) = max(pts.ow(satnum, 2), pts.og(satnum, 2));
                    elseif model.water
                        s_min(i) = pts.ow(satnum, 2);
                    elseif model.gas
                        s_min(i) = pts.og(satnum, 2);
                    end
                else
                    s_min(i) = pts.(ph)(satnum, 2);
                end
            end
        end
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

function [m, c, p, k] = getTwoPointScalers(pts, ph, reg, f, cells)
    [get, CR, U, L, KM] = getSatPointPicker(f, pts, reg, cells);
    switch ph
        case {'w', 'g'}
            [su, SU] = get(ph, U);
        case {'ow', 'og'}
            if isfield(f.krPts, 'w')
                [swl, SWL] = get('w', L);
            else
                swl = 0; SWL = 0;
            end
            if isfield(f.krPts, 'g')
                [sgl, SGL] = get('g', L);
            else
                sgl = 0; SGL = 0;
            end
            SU = 1 - SWL - SGL;
            su = 1 - swl - sgl;
    end
    [scr, SCR] = get(ph, CR);
    p = cell(1, 2);
    p{1} = SCR;
    m    = (su-scr)./(SU-SCR);
    c    = scr - SCR.*m;
    p{2} = SU;
    [k{1}, k{2}] = get(ph, KM);
end

function [m, c, p, k] = getThreePointScalers(pts, ph, reg, f, cells)
    [get, CR, U, L, KM] = getSatPointPicker(f, pts, reg, cells);
    switch ph
        case 'w'
            [sowcr, SOWCR] = get('ow', CR);
            [sgl, SGL] = get('g', L);
            
            SR = 1 - SOWCR - SGL;
            sr = 1 - sowcr - sgl;
            
            [su, SU] = get('w', U);  
        case 'g'
            [sogcr, SOGCR] = get('og', CR);
            [swl, SWL] = get('w', L);
            SR = 1 - SOGCR - SWL;
            sr = 1 - sogcr - swl;
            
            [su, SU] = get('g', U);
        case 'ow'
            [swcr, SWCR] = get('w', CR);
            [sgl, SGL] = get('g', L);
            
            SR = 1 - SWCR - SGL;
            sr = 1 - swcr - sgl;

            [swl, SWL] = get('w', L);
            [sgl, SGL] = get('g', L);
            
            SU = 1 - SWL - SGL;
            su = 1 - swl - sgl;
        case 'og'
            [sgcr, SGCR] = get('g', CR);
            [swl, SWL] = get('w', L);
            [sgl, SGL] = get('g', L);
            SR   = 1 - SGCR - SWL;
            sr   = 1 - sgcr - swl;

            SU   = 1 - SGL - SWL;
            su = 1 - sgl - swl;
    end
    [scr, SCR] = get(ph, CR);

    m = cell(1, 2);
    p = cell(1, 3);
    c = cell(1, 2);

    p{1} = SCR;
    m{1} = (sr-scr)./(SR-SCR);
    c{1} = scr - SCR.*m{1};
    p{2} = SR;

    p{3} = SU;

    m{2} = (su-sr)./(SU-SR);
    c{2} = sr - SR.*m{2};
    ix = SU <= SR;
    if nnz(ix) > 0
        m{2}(ix) = 0;
        c{2}(ix) = 0;
    end
    [k{1}, k{2}] = get(ph, KM);
end

function [get, CR, U, L, KM] = getSatPointPicker(f, pts, reg, cells)
    L  = 1;
    CR = 2;
    U  = 3;
    KM = 4;
    
    tbl = @(phase, index) f.krPts.(phase)(reg, index);
    scal = @(phase, index) pts.(phase)(cells, index);
    get = @(phase, index) getPair(phase, index, tbl, scal);
    
end

function [v1, v2] = getPair(phase, index, fn1, fn2)
    v1 = fn1(phase, index);
    v2 = fn2(phase, index);
    % When scaling values are not given, they fall back to tables values
    ind = isnan(v2);
    v2(ind) = v1(ind);
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
