function kr = eCPARelativePermeabilityScaling(model, state)
% Phase relative permeability (per cell)

kr = relPermUnified(model, state);

end

% --------------------------------------------------------------------%

function kr = relPermUnified(model, state)
% Rel.perm without special oil treatment
phases = model.getPhaseNames();
snames = arrayfun(@(x) ['s', x], phases, 'UniformOutput', false);
nph = numel(phases);
s = cell(1, nph);
kr = cell(1, nph);
[s{:}] = model.getProps(state, snames{:});
if strcmpi(phases, 'OG')
    phases = {'og', 'g'};
else
    phases = {'ow', 'w'};
end
for i = 1:nph
    kr{i} = evaluatePhaseRelativePermeability(model, phases{i}, s{i});
end
end

function kr = evaluatePhaseRelativePermeability(model, phase, s, cells)
if nargin < 4
    cells = ':';
end

fn = model.fluid.(['kr', upper(phase)]);

[ss, kr_max_m] = scaleSaturation(model, s, phase, cells, 'drainage');
kr = kr_max_m.*evaluateFunctionCellSubset(fn, cells, ss);
end

function [s_scale, k_max_m] = scaleSaturation(model, s, phase, cells, type)
if nargin  < 4
    cells = ':';
end
if nargin  < 5
    type = 'drainage';
end

sv = value(s);
scaler = getScalers(model, phase, cells, type);
[p, c, m] = deal(scaler.p, scaler.c, scaler.m);
n_pts = numel(p);
if n_pts == 2 % 2-point
    ix1 = sv < p{1};
    ix2 = sv >= p{2};
    ix  = ~(ix1 | ix2);
    s_scale = (ix.*m).*s + (ix.*c + ix2);
elseif n_pts == 3
    ix1 = sv >= p{1} & sv < p{2};
    ix2 = sv >= p{2} & sv < p{3};
    ix3 = sv >= p{3};
    a = ix1.*m{1} + ix2.*m{2};
    s_scale = a.*s + ix1.*c{1} + ix2.*c{2} + ix3;
else
    error('Unknown number of scaling points');
end
k_max_m = scaler.k_max_m;
end

function scaler = getScalers(model, phase, cells, type)
% Get the scalers for a given phase(pair), cells and type
if nargin < 4
    type = 'drainage';
end
if nargin < 3
    cells = ':';
end
if 1
    pts = model.rock.krscale.(type);
    regions = ones(model.G.cells.num, 1);
    reg = regions(cells);
    npts = 2;
    f = model.fluid;
    if npts == 2
        [m, c, p, k] = getTwoPointScalers(pts, phase, reg, f, cells);
    elseif npts == 3
        [m, c, p, k] = getThreePointScalers(pts, phase, reg, f, cells);
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
    scaler = prop.scalers.(type).(lower(phase));
    if ~ischar(cells)
        for f = {'m', 'c', 'p', 'k'}
            fi = f{1};
            scaler.(fi) = applyFunction(@(x) x(cells, :), scaler.(fi));
        end
        scaler.k_max_m = scaler.k_max_m(cells);
    end
end
end

function v = evaluateFunctionCellSubset(fn, subset, varargin)
% Evaluate specific function on a given subset
if iscell(fn)
    local_region = prop.regions;
    if isempty(local_region)
        error(['The function provided is a cell array, which', ...
            ' indicates that multiple regions are present.', ...
            ' This instance of %s has empty .regions.', ...
            ' An region must be provided when', ...
            ' the input function is a cell array'], class(prop));
    end
    if prop.isSingleRegion
        lr = local_region(1);
        v = fn{lr}(varargin{:});
    else
        if isnumeric(subset) || islogical(subset)
            local_region = local_region(subset);
        end
        % We have multiple regions and have to evaluate for each
        % subregion
        nc = size(prop.regions, 1);
        isCell = cellfun(@(x) numelValue(x) == nc, varargin);
        [sample, isAD] = getSampleAD(varargin{:});
        v = zeros(numel(local_region), 1);
        if isAD
            v = prop.AutoDiffBackend.convertToAD(v, sample);
        end
        for reg = 1:numel(fn)
            act = local_region == reg;
            arg = varargin;
            carg = cellfun(@(x) x(act), arg(isCell), 'UniformOutput', false);
            [arg{isCell}] = carg{:};
            if any(act)
                v(act) = fn{reg}(arg{:});
            end
        end
    end
else
    v = fn(varargin{:});
end
end
function [m, c, p, k] = getTwoPointScalers(pts, ph, reg, f, cells)
% Get scaling factors for two-point rel.perm. scaling
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
% Get scaling factors for three-point rel.perm. scaling
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
getter = @(phase, index) getPair(phase, index, tbl, scal);
end

function [v1, v2] = getPair(phase, index, fn1, fn2)
v1 = fn1(phase, index);
v2 = fn2(phase, index);
% When scaling values are not given, they fall back to tables values
ind = isnan(value(v2));
v2(ind) = v1(ind);
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