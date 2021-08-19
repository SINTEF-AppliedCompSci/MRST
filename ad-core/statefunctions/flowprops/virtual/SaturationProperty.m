classdef SaturationProperty
    % Virtual class for utilities related to saturation function endpoint
    % scaling
    properties
        scalingActive = false;
        swcon = [];
        scalers
    end
    
    properties (Access = protected)
        cell_subset = ':';
    end

    methods
        function prop = storeConnateWater(prop, model, varargin)
            prop.swcon = prop.getConnateWater(model, varargin{:});
        end

        function prop = storeScalers(prop, model, varargin)
            % Store scalers in state function for faster evaluation
            r = model.rock;
            if isfield(r, 'krscale')
                tmp = struct();
                for type = reshape(fieldnames(r.krscale), 1, [])
                    t = type{1};
                    fn = fieldnames(r.krscale.(t));
                    fn_f = fieldnames(model.fluid.krPts);
                    fn = intersect(fn, fn_f);
                    s = struct();
                    for i = 1:numel(fn)
                        f = fn{i};
                        s.(f) = prop.getScalers(model, f, ':', t);
                    end
                    tmp.(t) = s;
                end
                prop.scalers = tmp;
            end
        end
        
        function swcon = getConnateWater(prop, model, state)
            % Get the connate water in each cell of the domain. This may
            % come from either the fluid or the rock, depending if the
            % scaling is active.
            if ~isempty(prop.swcon)
                swcon = prop.swcon;
            elseif prop.scalingActive && ...
                    isfield(model.rock, 'krscale') && ...
                    isfield(model.rock.krscale.drainage, 'w')
                % Connate water in rock, endpoint scaling
                cix = prop.cell_subset;
                swcon = model.rock.krscale.drainage.w(cix, 1);
                % check for defaulted (nan) swcon -> use table values
                nix = isnan(value(swcon));
                if any(nix)
                    swcon(nix) = reshape(model.fluid.krPts.w(prop.regions(nix), 1), [], 1);
                end
            elseif isfield(model.fluid, 'krPts')
                % Connate water from rel perm table
                swcon = reshape(model.fluid.krPts.w(prop.regions, 1), [], 1);
            else
                % No found - zero swcon in every cell
                swcon = zeros(model.G.cells.num, 1);
            end
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

        function scaler = getScalers(prop, model, phase, cells, type)
            % Get the scalers for a given phase(pair), cells and type
            if nargin < 5
                type = 'drainage';
            end
            if nargin < 4
                cells = ':';
            end
            if isempty(prop.scalers)
                pts = model.rock.krscale.(type);
                reg = prop.regions(cells);
                npts = prop.relpermPoints;
                f = model.fluid;
                if npts == 2
                    [m, c, p, k] = prop.getTwoPointScalers(pts, phase, reg, f, cells);
                elseif npts == 3
                    [m, c, p, k] = prop.getThreePointScalers(pts, phase, reg, f, cells);
                else
                    error('Unknown number of scaling points');
                end
                scaler = struct('m', {m},...
                                'c', {c}, ...
                                'p', {p}, ...
                                'k', {k}, ...
                                'k_max_m', k{2}./k{1}, ...
                                'points', npts);
            else
                % We already had it stored, we can just retrieve and pick
                % the subset.
                scaler = prop.scalers.(type).(phase);
                if ~ischar(cells)
                    for f = {'m', 'c', 'p', 'k'}
                        fi = f{1};
                        scaler.(fi) = applyFunction(@(x) x(cells, :), scaler.(fi));
                    end
                    scaler.k_max_m = scaler.k_max_m(cells);
                end
            end
        end
    end

    methods (Static)
        function [m, c, p, k] = getTwoPointScalers(pts, ph, reg, f, cells)
            % Get scaling factors for two-point rel.perm. scaling
            [get, CR, U, L, KM] = SaturationProperty.getSatPointPicker(f, pts, reg, cells);
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
            % Get scaling factors for three-point rel.perm. scaling
            [get, CR, U, L, KM] = SaturationProperty.getSatPointPicker(f, pts, reg, cells);
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
                otherwise
                    error('No valid scalers for phase %s', ph);
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

        function [getter, CR, U, L, KM] = getSatPointPicker(f, pts, reg, cells)
            % Get function handle for getting saturation-based scaling
            % points
            L  = 1; % Connate (always present)
            CR = 2; % Critical (first point where phase becomes mobile)
            U  = 3; % Saturation at which maximum rel. perm occurs
            KM = 4; % Maximum rel. perm.

            tbl = @(phase, index) f.krPts.(phase)(reg, index);
            scal = @(phase, index) pts.(phase)(cells, index);
            getter = @(phase, index) SaturationProperty.getPair(phase, index, tbl, scal);
        end

        function [v1, v2] = getPair(phase, index, fn1, fn2)
            v1 = fn1(phase, index);
            v2 = fn2(phase, index);
            % When scaling values are not given, they fall back to tables values
            ind = isnan(value(v2));
            v2(ind) = v1(ind);
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
