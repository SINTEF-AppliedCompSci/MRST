function W = addWellShale(W, G, rock, cellInx, varargin)
%Insert a well into the simulation model.
%
% SYNOPSIS:
%   W = addWellShale(W, G, rock, cellInx)
%   W = addWellShale(W, G, rock, cellInx, 'pn', pv, ...)
%
% REQUIRED PARAMETERS:
%   W       - Well structure or empty if no other wells exist.
%             Updated upon return.
%
%   G       - Grid data structure.
%
%   rock    - Rock data structure.  Must contain valid field `perm`.
%
%   cellInx - Perforated well cells (vector of cell indices or a logical
%             mask with length equal to G.cells.num).
%
% OPTIONAL PARAMETERS:
%   Type -   Well control type. String. Supported values depend on the
%            solver in question. Most solvers support at least two options: 
%              - 'bhp': Controlled by bottom hole pressure (DEFAULT)
%              - 'rate': Controlled by total rate.
%
%   Val    - Well control target value.  Interpretation of this value is
%            dependent upon `Type`.  Default value is 0.  If the well
%            `Type` is 'bhp', then `Val` is given in unit Pascal and if
%            the `Type` is 'rate', then `Val` is given in unit m^3/second.
%
%   Radius - Well bore radius (in unit of meters).  Either a single,
%            scalar value which applies to all perforations or a vector of
%            radii (one radius value for each perforation).
%            Default value: Radius = 0.1 (i.e., 10 cm).
%
%   Dir    - Well direction.  Single CHAR or CHAR array.  A single CHAR
%            applies to all perforations while a CHAR array defines the
%            direction of the corresponding perforation.  In other words,
%            Dir(i) is the direction in perforation cellInx(i) in the CHAR
%            array case.
%
%            Supported values for a single perforation are 'x', 'y', or
%            'z' (case insensitive) meaning the corresponding cell is
%            perforated in the X, Y, or Z direction, respectively.
%            Default value: Dir = 'z' (vertical well).
%
%   InnerProduct - The inner product with which to define the mass matrix.
%            String.  Default value = 'ip_tpf'.
%            Supported values are 'ip_simple', 'ip_tpf', 'ip_quasitpf',
%            and 'ip_rt'.
%
%   WI     - Well productivity index. Vector of length `nc=numel(cellInx)`.
%            Default value: `WI = repmat(-1, [nc, 1])`, whence the
%            productivity index will be computed from available grid block
%            data in grid blocks containing well completions.
%
%   Kh     - Permeability thickness.  Vector of length `nc=numel(cellInx)`.
%            Default value: `Kh = repmat(-1, [nc, 1])`, whence the thickness
%            will be computed from available grid block data in grid
%            blocks containing well completions.
%
%   Skin   - Skin factor for computing effective well bore radius.  Scalar
%            value or vector of length `nc=numel(cellInx)`.
%            Default value: 0.0 (no skin effect).
%
%   Comp_i - Fluid composition for injection wells.  Vector of saturations.
%            Default value:  `Comp_i = [1, 0, 0]` (water injection)
%
%   Sign   - Well type: Production (Sign = -1) or Injection (Sign = 1).
%            Default value: 0 (Undetermined sign. Will be derived from
%            rates if possible).
%
%   Name   - Well name (string).
%            Default value: `sprintf('W%d', numel(W) + 1)`
%
%   refDepth - Reference depth for the well, i.e. the value for which
%            bottom hole pressures are defined.
%
%   calcReprRad - Whether or not to compute the representative radius of
%            each perforation.  The representative radius is needed to
%            derive shear thinning calculations in the context of polymer
%            injection.  LOGICAL.  Default value: `calcReprRad = true` (do
%            calculate the representative radius).  If set to false, the
%            resulting well structure cannot be used to simulate polymer
%            injection.
%
%   cellDims - optional cellDims of grid cells
%
% RETURNS:
%   W - Updated (or freshly created) well structure, each element of which
%       has the following fields:
%         - cells:   Grid cells perforated by this well (== cellInx).
%         - type:    Well control type (== Type).
%         - val:     Target control value (== Val).
%         - r:       Well bore radius (== Radius).
%         - dir:     Well direction (== Dir).
%         - WI:      Well productivity index.
%         - dZ: Displacement of each well perforation measured from
%           'highest' horizontal contact (i.e., the 'TOP' contact with the
%           minimum 'Z' value counted amongst all cells perforated by this
%           well).
%         - name:    Well name (== Name).
%         - compi:   Fluid composition--only used for injectors (== Comp_i).
%         - sign:    Injection (+) or production (-) flag.
%         - status:  Boolean indicating if the well is open or shut.
%         - cstatus: One entry per perforation, indicating if the completion is open.
%         - lims:    Limits for the well. Contains subfields for the types
%           of limits applicable to the well (bhp, rate, orat, ...)
%           Injectors generally have upper limits, while producers have
%           lower limits.
%         - rR:      The representative radius for the wells, which is used
%           in the shear thinning calculation when polymer is involved in
%           the simulation.  Empty if 'calcReprRad' is false.
%
% EXAMPLE:
%   incompTutorialWells
%
% NOTE:
%   Wells in two-dimensional grids are not well defined in terms of
%   computing well indices. However, such wells are often useful for
%   simulation scenarios where the exact value of well indices is not of
%   great importance. For this reason, we make the following approximations
%   when addWell is used to compute e.g. horizontal wells in 2D:
%
%       - K_z is assumed to be the harmonic average of K_x and K_y:
%         K_z = 1/(1/K_x + 1/K_y).
%       - The depth of a grid block is assumed to be unit length (1 meter)
%
%   This generally produces reasonable ranges for the WI field, but it is
%   the user's responsibility to keep these assumptions in mind.
%
% SEE ALSO:
%   `verticalWell`, `addSource`, `addBC`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.
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


if ~isempty(W) && ~isfield(W, 'WI'),
   error(msgid('CallingSequence:Changed'), ...
        ['The calling sequence for function ''addWell'' has changed\n', ...
         'Please use\n\tW = addWell(W, G, rock, cellInx, ...)\n', ...
         'from now on.']);
end

mrstNargInCheck(4, [], nargin);
if islogical(cellInx)
    assert(numel(cellInx) == G.cells.num, ...
    'Logical mask does not match number of grid cells.');
    cellInx = find(cellInx);
end
if(isfield(G,'FracGrid'))%EDFM On
    numC = size(cellInx,1);
else
    numC = numel(cellInx);
end


opt = struct('InnerProduct', 'ip_tpf',                     ...
             'Dir'         , 'z',                          ...
             'Name'        , sprintf('W%d', numel(W) + 1), ...
             'Radius'      , 0.1,                          ...
             'Type'        , 'bhp',                        ...
             'Val'         , 0,                            ...
             'Comp_i'      , [1, 0, 0],                    ...
             'WI'          , -1*ones(numC,1),              ...
             'Kh'          , -1*ones(numC,1),              ...
             'Skin'        , zeros(numC, 1),               ...
             'refDepth'    , 0,                            ...
             'lims'        , [],                           ...
             'vfp_index'   , 0,                            ...
             'c'           , [],                           ...
             'Sign'        , 0,                            ...
             'calcReprRad' , true,                         ...
             'cellDims',   []);

opt = merge_options(opt, varargin{:});

WI = reshape(opt.WI, [], 1);

assert (numel(WI)       == numC, ...
    'Provided WI should be one entry per perforated cell.')
assert (numel(opt.Kh)   == numC, ...
    'Provided Kh should be one entry per perforated cell.')
assert (numel(opt.Skin) == numC || numel(opt.Skin) == 1, ...
    ['Provided Skin should be one entry per perforated cell or a single', ...
    ' entry for all perforated cells']);

if numel(opt.Skin) == 1, opt.Skin = opt.Skin(ones([numC, 1]));  end

% Set reference depth default value.
if isempty(opt.refDepth)
   g_vec = gravity();
   dims  = G.griddim;
   if norm(g_vec(1:dims)) > 0
      g_vec = g_vec ./ norm(g_vec);
      opt.refDepth = min(G.nodes.coords * g_vec(1:dims)');
   else
      opt.refDepth = 0;
   end
end
ip = opt.InnerProduct;

if(isfield(G,'FracGrid'))%EDFM On
   rR = EDFMradiusRep(G, cellInx,opt.Radius); %OMO modification
elseif opt.calcReprRad
   % Compute the representative radius for the grid block in which the well
   % is completed.  It is needed for computing the shear rate of the wells.
   rR = radiusRep(G, opt.Radius, opt.Dir, reshape(cellInx, [], 1));
else
   rR = [];
end

% Initialize Well index - WI. ---------------------------------------------
% Check if we need to calculate WI or if it is supplied.

compWI = WI < 0;

if any(compWI) % calculate WI for the cells in compWI
   if(isfield(G,'FracGrid'))%EDFM On
       WI(compWI) = computeWIShale(G,opt.Radius,cellInx, ...
                                   'Skin', opt.Skin,...
                                   're', rR(compWI));
   else
       WI(compWI) = computeWellIndex(G, rock, opt.Radius, reshape(cellInx, [], 1), ...
                        'Dir', opt.Dir, ...
                        'Skin', opt.Skin, ...
                        'cellDims', opt.cellDims, ...
                        'Kh', opt.Kh, ...
                        'InnerProduct', ip, ...
                        'Subset', compWI);
   end
end
if(isfield(G,'FracGrid')), cellInx=cellInx(:,1); end

% Set well sign (injection = 1 or production = -1)
% for bhp wells or rate controlled wells with rate = 0.
if opt.Sign ~= 0
   if sum(opt.Sign == [-1, 1]) ~= 1
      error(msgid('Sign:NonUnit'), 'Sign must be -1 or 1');
   end
   if strcmp(opt.Type, 'rate') && (sign(opt.Val) ~= 0) ...
         && (opt.Sign ~= sign(opt.Val))
      warning(msgid('Sign'), ...
             ['Given sign does not match sign of given value. ', ...
              'Setting w.sign = sign( w.val )']);
      opt.Sign = sign(opt.Val);
   end
else
   if strcmp(opt.Type, 'rate')
      if opt.Val == 0
         warning(msgid('Sign'), 'Given value is zero, prod or inj ???');
      else
         opt.Sign = sign(opt.Val);
      end
   end
end

% Add well to well structure. ---------------------------------------------
%
W  = [W; struct('cells'    , cellInx(:),           ...
                'type'     , opt.Type,             ...
                'val'      , opt.Val,              ...
                'r'        , opt.Radius,           ...
                'dir'      , opt.Dir,              ...
                'rR'       , rR,                   ...
                'WI'       , WI,                   ...
                'dZ'       , getDeltaZ(G, cellInx(:), opt.refDepth), ...
                'name'     , opt.Name,             ...
                'compi'    , opt.Comp_i,           ...
                'refDepth' , opt.refDepth,         ...
                'lims'     , opt.lims,             ...
                'sign'     , opt.Sign,             ...
                'c'        , opt.c,                ...
                'status'   , true,                 ...
                'vfp_index', opt.vfp_index,        ...
                'cstatus'  , true(numC,1))];

if numel(W(end).dir) == 1
   W(end).dir = repmat(W(end).dir, [numel(W(end).cells), 1]);
end
assert (numel(W(end).dir) == numel(W(end).cells));

%--------------------------------------------------------------------------
% Private helper functions follow
%--------------------------------------------------------------------------

function dZ = getDeltaZ(G, cells, refDepth)
direction = gravity();
dims      = G.griddim;
if norm(direction(1:dims)) > 0
   direction = direction ./ norm(direction(1:dims));
else
   direction      = zeros(1, dims);
   direction(end) = 1;
end
Z = G.cells.centroids(cells, :) * direction(1:dims).';
dZ = Z - refDepth;
if any(dZ < 0)
    warning('Negative distance between perforation and reference depth.')
end
if isfield(G, 'nodes')
    % Grid with nodes
    delta = max(G.nodes.coords)-min(G.nodes.coords);
elseif isfield(G.faces, 'centroids')
    % Possibly a coarse grid
    delta = max(G.faces.centroids)-min(G.nodes.centroids);
else
    % Skip check
    delta = inf;
end

if max(dZ) > delta*direction(1:dims).'
    msg = 'Pressure drop from refDepth to perforation exceeds depth of model.';
    if refDepth == 0
        msg = [msg, ' refDepth is defaulted to zero. Consider setting a different value.'];
    else
        msg = [msg, ' Please check refDepth value.'];
    end
    warning(msg)
end

%--------------------------------------------------------------------------
% A function to compute the representative radius of the grid block in
% which the well is completed.
% rR = sqrt(re * rw).
% Here, rw is the wellbore radius, re is the area equivalent radius of the
% grid block where the well is completed, which means the circle with radius
% re has the same area as the cross section of the grid block in the
% completion's orthogonal direction.
% The current formulation theoretically only works for Cartisian grids,
% while it has been working well for the cases we have,
% including some corner-point grids.
% TODO: REMAIN TO BE VERIFIED for really twisted grids.
function rr = radiusRep(G, radius, welldir, cells)

if(isfield(G,'nodes'))
   [dx, dy, dz] = cellDims(G, cells);
else
   [dx, dy, dz] = cellDimsCG(G, cells);
end

welldir = lower(welldir);
if numel(welldir) == 1
    welldir = repmat(welldir, numel(cells), 1);
end
re = zeros(size(welldir, 1), 1);

% The following formulation only works for Cartisian mesh
ci = welldir == 'x';
re(ci) = sqrt(dy(ci) .* dz(ci) / pi);

ci = welldir == 'y';
re(ci) = sqrt(dx(ci) .* dz(ci) / pi);

ci = welldir == 'z';
re(ci) = sqrt(dy(ci) .* dx(ci) / pi);

rr = sqrt( re .* radius);


%--------------------------------------------------------------------------
% A function to compute the representative radius of the frac grid block 
% in EDFM
function rr = EDFMradiusRep(G, cells,radius)
    NumFracs=length(fieldnames(G.FracGrid));
    
    rr=zeros(size(cells,1),1);
    for i=1:NumFracs
        mask=find(cells(:,2)==i);
        startID=G.FracGrid.(['Frac' num2str(i)]).cells.start;
        localFracCells=cells(mask)-startID+1;
        
        FracAperture=G.FracGrid.(['Frac' num2str(i)]).aperture;
        WellCellVolume=G.FracGrid.(['Frac' num2str(i)]).cells.volumes(localFracCells);
        Areas=WellCellVolume/FracAperture;
        
        re(mask)=sqrt(Areas/pi);
        rr = sqrt( re .* radius); %OMO modification
        
    end
