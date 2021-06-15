function fluid = initDeckADIFluid(deck, varargin)
%Initialize AD-solver fluid from ECLIPSE-style input deck
%
% SYNOPSIS:
%   f = initDeckADIFluid(deck)
%
% REQUIRED PARAMETERS:
%   deck - Output structure from readEclipseDeck
%
% OPTIONAL PARAMETERS:
%   G    - Grid to be used for the fluid model. Only required for the
%          situation when a deck has more than one region, and the intended
%          simulation grid is different from the one defined by the deck
%          (more specifically, the deck.GRID.ACTNUM field).
%
%   region_method - Method for defining regions. Either 'deck' (use exactly as
%           prescribed in the input deck [DEFAULT] or 'override' which
%           allows individual fields to be overwritten. Since the region
%           support in MRST is abstracted away, this option should only be
%           used if you are comfortable with the internal structure of
%           getRegMap and interpReg.
%
%   regionOverride - Struct containing fields to be overwritten. Only used
%           if 'method' equals 'override'. The supported fields include
%           PVTNUM, SATNUM, ROCKNUM, IMBNUM and SURFNUM. For each of these
%           fields, if present, the value can be formatted in three
%           different ways:
%               * Single numerical value. Will be repeated for all grid
%                 cells.
%               * The char ':'. Interpreted as using the first region for
%                 all cells. This is normally used when no regions are
%                 present.
%               * A Nx1 array containing one region indicator per cell in
%                 the grid. The region indicator is a single number
%                 indicating which table is to be used for that cell.
%
%   singleRegion - Internal debug option. This region will be used for ALL
%                 fields and any other value than 1 or ':' will usually
%                 result in exceptions being thrown.
%
% RETURNS:
%   f - Fluid model suitable for several different MRST models.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, initEclipseDeck

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

    % this only work for full deck or first region.
    opt = struct('G',               [], ...
                 'region_method',   'deck', ...
                 'singleRegion',    [], ...
                 'resamplePoly',    0, ...
                 'pvtMethodOil',    'parallel', ...
                 'pvtMethodGas',    'linshift', ...
                 'prange',          [], ...
                 'ignoreErrors',    true, ...
                 'useMex',          false, ...
                 'optimize',        false, ...
                 'interp1d',        [], ...
                 'regionOverride',  struct());

    opt = merge_options(opt,varargin{:});
    reg = getRegions(deck, opt);
    % Pass interpolation along with regions to avoid breaking interface
    if opt.optimize && opt.resamplePoly == 0
        opt.resamplePoly = 20;
    end
    if isempty(opt.interp1d)
        if opt.useMex
            reg.interp1d = @interpTableMEX;
            reg.interp1d_uniform = @(X, F, x) interpTableMEX(X, F, x, 3);
        else
            reg.interp1d = @interpTable;
            reg.interp1d_uniform = reg.interp1d;
        end
    else
        reg.interp1d = opt.interp1d;
        reg.interp1d_uniform = reg.interp1d;
    end
    reg.pvtMethodOil = opt.pvtMethodOil;
    reg.pvtMethodGas = opt.pvtMethodGas;
    reg.prange = opt.prange;
    reg.useMex = opt.useMex;
    reg.optimize = opt.optimize; % Perform additional optimizations for large cases
    if isempty(reg.prange)
       reg.prange = get_pressures(deck, opt);
    end

    fluid = struct();

    % Properties
    props = deck.PROPS;

    for fld = assignable_fields(prioritise(fieldnames(props)))
        try
            fluid = feval(['assign', fld{1}], fluid, props.(fld{1}), reg);
        catch ME
            if opt.ignoreErrors
                warning(msgid('Assign:Failed'), ...
                    ['Could not assign property ''%s''. ', ...
                    'Encountered error: ''%s'''], fld{1}, ME.message);
            else
                rethrow(ME);
            end
        end
    end

    fn = fieldnames(fluid);
    for i = 1:numel(fn)
        f = fn{i};
        if iscell(fluid.(f)) && numel(fluid.(f)) == 1
            fluid.(f) = fluid.(f){1};
        end
    end
    if isfield(deck.RUNSPEC, 'GAS') && any(fluid.rhoGS == 0)
        fluid.rhoGS(fluid.rhoGS == 0) = 1;
    end
    if isfield(deck.RUNSPEC, 'OIL') && any(fluid.rhoOS == 0)
        fluid.rhoOS(fluid.rhoOS == 0) = 1;
    end
end

%--------------------------------------------------------------------------

function reg = getRegions(deck, opt)
    % Legacy regions
    reg = handleRegions(deck, opt.G);
    % Modern region treatment
    reg.sat = 1;
    reg.pvt = 1;
    if isfield(deck.RUNSPEC, 'TABDIMS')
        tab = deck.RUNSPEC.TABDIMS;
        reg.sat = tab(1);
        reg.pvt = tab(2);
    else
        if ~isempty(reg.PVTNUM)
            reg.pvt = numel(reg.PVTINX);
            warning('RUNSPEC does not include TABDIMS, number of PVT-tables set to %d', reg.pvt)
        end
        if ~isempty(reg.SATNUM)
            reg.sat = numel(reg.SATINX);
            warning('RUNSPEC does not include TABDIMS, number of saturation-tables set to %d', reg.sat)
        end
    end
end
%--------------------------------------------------------------------------
function prange = get_pressures(deck, opt)
    pmax = -inf;
    pmin = inf;
    props = deck.PROPS;
    f = {'PVDG', 'PVDO'};
    for i = 1:numel(f)
        fn = f{i};
        if isfield(props, fn)
            d = props.(fn);
            for j = 1:numel(d)
                pi = d{j}(:, 1);
                pmax = max(pmax, max(pi));
                pmin = min(pmin, min(pi));
            end
        end
    end
    if ~isfinite(pmax)
        pmax = 1000*barsa;
    end
    if ~isfinite(pmin)
        pmin = 0;
    end
    prange = linspace(pmin, pmax, opt.resamplePoly)';
end
%--------------------------------------------------------------------------

function fnames = prioritise(fnames)
% Pull saturation functions for water to front of list to ensure that SWCON
% exists when processing saturation functions for gas.

   sfunc_wat = { 'SWFN', 'SWOF' };
   [in, pos] = ismember(sfunc_wat, fnames);

   if any(in)
      shape = size(fnames);

      fnames(pos(in)) = [];
      fnames = [reshape(sfunc_wat(in), [], 1) ; ...
                reshape(fnames       , [], 1)];

      fnames = reshape(fnames, shape);
   end
end

%--------------------------------------------------------------------------

function flds = assignable_fields(fnames)
   excpt = reshape(exclude_properties(), [], 1);
   n     = numel(excpt);

   [i, j] = blockDiagIndex(numel(fnames), n);

   mtch = strcmp(reshape(fnames(i), [], n), ...
                 reshape(excpt (j), [], n));

   flds = reshape(fnames(~ any(mtch, 2)), 1, []);
end

%--------------------------------------------------------------------------

function excpt = exclude_properties()
   % Properties not resulting in individual functions
   excpt = {'KRW'   , 'KRO'   , 'KRG',   ...
            'SWL'   , 'SWCR'  , 'SWU',   ...
            'SGL'   , 'SGCR'  , 'SGU',   ...
            'SOWCR' , 'SOGCR' ,          ...
            'ISWL'  , 'ISWCR' , 'ISWU',  ...
            'ISGL'  , 'ISGCR' , 'ISGU',  ...
            'ISOWCR', 'ISOGCR',          ...
            'CNAMES', 'BIC'   , 'ACF',   ...
            'PCRIT' , 'TCRIT' , 'VCRIT', ...
            'MW'    , 'ZCRIT' , 'EOS',    ...
            'ZMFVD' , 'TEMPVD', 'SSHIFT', ...
            'SCALECRS', 'SWATINIT', ...
            'PRCORR', 'STCOND' ...
            };
end
