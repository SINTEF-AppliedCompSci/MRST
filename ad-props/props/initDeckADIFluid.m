function fluid = initDeckADIFluid(deck, varargin)
%Initialize AD-solver fluid from ECLIPSE-style input deck
%
% SYNOPSIS:
%   f = initDeckADIFluid(deck)
%
% REQUIRED PARAMETERS:
%   deck - Output structure from readEclipseDeck
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs, ('pn'/pv)):
%   G    - Grid to be used for the fluid model. Only required for the
%          situation when a deck has more than one region, and the intended
%          simulation grid is different from the one defined by the deck
%          (more specifically, the deck.GRID.ACTNUM field).
%
%  region_method - Method for defining regions. Either 'deck' (use exactly as
%           prescribed in the input deck [DEFAULT] or 'override' which
%           allows individual fields to be overwritten. Since the region
%           support in MRST is abstracted away, this option should only be
%           used if you are comfortable with the internal structure of
%           getRegMap and interpReg.
%
%  regionOverride - Struct containing fields to be overwritten. Only used
%           if 'method' equals 'override'. The supported fields include
%           PVTNUM, SATNUM, ROCKNUM, IMBNUM and SURFNUM. For each of these
%           fields, if present, the value can be formatted in three
%           different ways:
%               * Single numerical value. Will be repeated for all grid
%               cells.
%               * The char ':'. Interpreted as using the first region for
%               all cells. This is normally used when no regions are
%               present.
%               * A Nx1 array containing one region indicator per cell in
%               the grid. The region indicator is a single number
%               indicating which table is to be used for that cell.
%
%  singleRegion - Internal debug option. This region will be used for ALL
%                 fields and any other value than 1 or ':' will usually
%                 result in exceptions being thrown.
%
% RETURNS:
%   f - Fluid model suitable for several different MRST models.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, initEclipseDeck

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
opt = struct('G',              [], ...
             'region_method', 'deck', ...
             'singleRegion',   [], ...
             'regionOverride', struct());
 
opt = merge_options(opt,varargin{:});
reg = handleRegions(deck, opt.G);
switch opt.region_method
     case 'deck'
         % Default behavior - do nothing
     case 'override'
         % In override mode, it is possible to overwrite 
         nc = getNumberOfCells(opt.G, deck);
         flds = {'PVTNUM', 'SATNUM', 'SURFNUM', 'IMBNUM', 'ROCKNUM'};
         for i = 1:numel(flds)
             f = flds{i};
             if isfield(opt.regionOverride, f)
                 % First check if this specific field has been overridden
                 d = expandValue(opt.regionOverride.(f), nc);
             elseif ~isempty(opt.singleRegion)
                 % Otherwise just use the expanded single region
                 d = expandValue(opt.singleRegion, nc);
             else
                 % No override found, just use whatever was in the deck
                 continue
             end
             % Get the maximum allowable value
             switch f
                 case 'PVTNUM'
                     maxval = deck.RUNSPEC.TABDIMS(2);
                 case 'SATNUM'
                     maxval = deck.RUNSPEC.TABDIMS(1);
                 case {'IMBNUM', 'ROCKNUM', 'SURFNUM'}
                     if isfield(deck.REGIONS, f)
                         maxval = deck.REGIONS.(f);
                     end
                 otherwise
                     if ischar(d)
                         maxval = inf;
                     else
                         maxval = max(d);
                     end
             end
             % Sanity check
             if isnumeric(d)
                assert(all(d <= maxval), ...
                    ['Input override for ''', f, ''' exceeds maximum allowable value ''', num2str(maxval),'''']);
             end
             % Set INX values, used internally for the region interpolator
             reg.(f) = d;
             inx = [f(1:end-3), 'INX'];
             inx_v = arrayfun(@(x) find(x==d), 1:maxval, 'UniformOutput', false);

             reg.(inx) = inx_v;
         end
     otherwise
         error('No such region method')
end
fluid = struct();
% Properties
props = deck.PROPS;
fns = fieldnames(props);
for k = 1:numel(fns)
    fn = fns{k};
    if doAssign(fn)
        asgn = str2func(['assign',fn]);
        try
            fluid = asgn(fluid, props.(fn), reg);
        catch ME
            warning(msgid('Assign:Failed'), ...
            'Could not assign property ''%s''. Encountered error: ''%s''',...
                fn, ME.message);
        end
    end
end
fluid = assignRelPerm(fluid);
end

function flag = doAssign(propNm)
% Properties not resulting in individual functions
excpt = {'SWL'   ,'SWCR'   ,'SWU' , ...
         'SGL'   ,'SGCR'   ,'SGU' , ...
         'SOWCR' ,'SOGCR'  , ...
         'CNAMES','BIC'   , 'ACF', ...
         'PCRIT' ,'TCRIT' , 'VCRIT',...
         'MW',    'ZCRIT', ...
         'ISWL'  ,'ISWCR'  ,'ISWU', ...
         'ISGL'  ,'ISGCR'  ,'ISGU', ...
         'ISOWCR','ISOGCR'};
flag = ~any( strcmp(propNm , excpt) );
end

function nc = getNumberOfCells(G, deck)
    if ~isempty(G)
        nc = G.cells.num;
    else
        nc = sum(deck.GRID.ACTNUM);
    end
end

function v = expandValue(v, nc)
    if ischar(v)
        assert(strcmp(v, ':'), 'Only valid char region is '':''.');
    else
        assert(isnumeric(v));
        if numel(v) == 1
            v = repmat(v, nc, 1);
        end
        assert(numel(v) == nc, 'Values must be provided per cell or one for entire grid.');
    end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
