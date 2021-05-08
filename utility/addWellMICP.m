function W = addWellMICP(W, G, rock, cellInx, varargin)
% Function to insert a well into the simulation model.
% 
% This function is extended from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see 
%   mrst/params/wells_and_bc/addWell.m
%
% We refer to that function for a complete commented version of the file. 
% In this file we comment on the new added lines.

%{
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

mrstNargInCheck(4, [], nargin);

if isempty(cellInx) || (islogical(cellInx) && ~any(cellInx))
   warn_empty_cell_index_range(numel(W), varargin{:});
   return
end

if islogical(cellInx)
    assert(numel(cellInx) == G.cells.num, ...
    'Logical mask does not match number of grid cells.');
    cellInx = find(cellInx);
end
numC = numel(cellInx);


opt = struct('InnerProduct', 'ip_tpf',                     ...
             'Dir'         , 'z',                          ...
             'Name'        , sprintf('W%d', numel(W) + 1), ...
             'Radius'      , 0.1,                          ...
             'Type'        , 'bhp',                        ...
             'Val'         , 0,                            ...
             'Comp_i'      , [],                           ...
             'compi'       , [],                           ...
             'WI'          , repmat(-1, [numC, 1]),        ...
             'Kh'          , repmat(-1, [numC, 1]),        ...
             'Skin'        , zeros(numC, 1),               ...
             'refDepth'    , [],                           ...
             'lims'        , [],                           ...
             'vfp_index'   , [],                           ...
             'o'           , [],                           ...% oxygen
             'u'           , [],                           ...% urea
             'm'           , [],                           ...% microbes
             'b'           , [],                           ...% biofilm
             'c'           , [],                           ...% calcite
             'cp'          , [],                           ...
             'components'  , [],                           ...
             'heel'        , [],                           ...
             'Sign'        , 0,                            ...
             'calcReprRad' , true,                         ...
             'cellDims'    , [],                           ...
             'lineSegments', []);

opt = merge_options(opt, varargin{:});

if ~isempty(opt.Comp_i)
    if numel(opt.compi) ~= numel(opt.Comp_i)|| any(opt.compi ~= opt.Comp_i)
        assert(isempty(opt.compi), ...
         ['Ambigious injection phase composition. Both ''compi'' and ', ...
                       '''Comp_i'' was specified with different values.']);
        opt.compi = opt.Comp_i;
    end
    opt.Comp_i = [];
end

if isempty(opt.compi)
    opt.compi = [1, 0, 0];
end

if numel(opt.WI) == 1
    opt.WI = repmat(opt.WI, numC, 1);
end
if numel(opt.Skin) == 1
    opt.Skin = repmat(opt.Skin, numC, 1);
end
if numel(opt.Kh) == 1
    opt.Kh = repmat(opt.Kh, numC, 1);
end

WI = reshape(opt.WI, [], 1);

assert (numel(WI) == numC, ...
        'Provided WI should be one entry per perforated cell.')
assert (numel(opt.Kh) == numC, ...
        'Provided Kh should be one entry per perforated cell.')
assert (numel(opt.Skin) == numC || numel(opt.Skin) == 1, ...
       ['Provided Skin should be one entry per perforated cell or ', ...
        'a single entry for all perforated cells']);

if numel(opt.Skin) == 1
   opt.Skin = opt.Skin(ones([numC, 1]));
end

if isempty(opt.refDepth)
   dispif(mrstVerbose(), ...
 ['Defaulting reference depth to top of formation for well %s. Please', ...
' specify ''refDepth'' optional argument to ''addWell'' if you',...
'require bhp  at a specific depth.\n'], opt.Name);
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

if opt.calcReprRad
   rR = radiusRep(G, opt.Radius, opt.Dir, reshape(cellInx, [], 1));
else
   rR = [];
end

compWI = WI < 0;

if any(compWI) 
    if isempty(opt.lineSegments)
        WI(compWI) = computeWellIndex(G, rock, opt.Radius, ...
            reshape(cellInx, [], 1), ...
            'Dir', opt.Dir, ...
            'Skin', opt.Skin, ...
            'cellDims', opt.cellDims, ...
            'Kh', opt.Kh, ...
            'InnerProduct', ip, ...
            'Subset', compWI);
    else
        len = opt.lineSegments;
        assert(all(size(opt.lineSegments)== [numC, 3]), ...
      'Expected format of input ''segments'' is [%d, 3], got [%d, %d]', ...
            numC, size(len,1), size(len,2));
        d = {'x', 'y', 'z'}; 
        tmp = 0; 
        if isempty(opt.cellDims)
            [dx, dy, dz] = cellDims(G, cellInx);
            opt.cellDims = [dx, dy, dz];
        end
        for k = 1:3
            WI_d = computeWellIndex(G, rock, opt.Radius, ...
                   reshape(cellInx, [], 1), ...
                   'Dir', d{k}, ...
                   'Skin', opt.Skin, ...
                   'cellDims', opt.cellDims, ...
                   'Kh', opt.Kh, ...
                   'InnerProduct', ip, ...
                   'Subset', compWI);         
        tmp = tmp + ( (len(compWI, k)./opt.cellDims(compWI, k)).*WI_d ).^2;
        end
        WI(compWI) = sqrt(tmp);
    end
end

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

if numel(opt.Radius) == 1
    opt.Radius = repmat(opt.Radius, numC, 1);
end
W  = [W; struct('cells'     , reshape(cellInx, [], 1), ...
                'type'      , opt.Type,             ...
                'val'       , opt.Val,              ...
                'r'         , opt.Radius,           ...
                'dir'       , opt.Dir,              ...
                'rR'        , rR,                   ...
                'WI'        , WI,                   ...
                'dZ'        , getDeltaZ(G, reshape(cellInx, [], 1), ...
                                           opt.refDepth, opt.Name), ...
                'heel'      , opt.heel,             ...
                'name'      , opt.Name,             ...
                'compi'     , opt.compi,            ...
                'refDepth'  , opt.refDepth,         ...
                'defaulted' , [],                   ...
                'lims'      , opt.lims,             ...
                'o'         , opt.o,                ...%oxygen
                'u'         , opt.u,                ...%urea
                'm'         , opt.m,                ...%microbes
                'b'         , opt.b,                ...%biofilm
                'c'         , opt.c,                ...%calcite
                'cp'        , opt.cp,               ...
                'components', opt.components,      ...
                'sign'      , opt.Sign,             ...
                'status'    , true,                 ...
                'vfp_index' , opt.vfp_index,        ...
                'cstatus'   , true(numC,1))];

if numel(W(end).dir) == 1
   W(end).dir = repmat(W(end).dir, [numel(W(end).cells), 1]);
end
assert (numel(W(end).dir) == numel(W(end).cells));

function dZ = getDeltaZ(G, cells, refDepth, wellName)
direction = gravity();
dims      = G.griddim;
if norm(direction(1:dims)) > 0
   direction = direction ./ norm(direction(1:dims));
else
   direction      = zeros(1, dims);
   if dims > 2
      direction(end) = 1;
   end
end
xyz = G.cells.centroids(cells, :);
bad = any(~isfinite(xyz), 1);
if norm(direction(bad))
    error('Non-finite centroids in dimensions: %s', num2str(find(bad)));
end
xyz(isnan(xyz)) = 0;
Z = xyz * direction(1:dims).';
dZ = Z - refDepth;
if any(dZ < -max(eps(Z), eps(refDepth)))
    warning('RefDepth:BelowTopConnection', ...
           ['Reference depth for well BHP in well ''%s'' is set ', ...
            'below well''s top-most reservoir connection'], wellName);
end
if isfield(G, 'nodes')
    delta = max(G.nodes.coords)-min(G.nodes.coords);
elseif isfield(G.faces, 'centroids')
    delta = max(G.faces.centroids)-min(G.faces.centroids);
else
    delta = inf(1, G.griddim);
end

if max(dZ) > delta*direction(1:dims).' && norm(gravity()) > 0
    msg = ['Pressure drop distance from BHP reference depth to ', ...
           'selected perforation exceeds total thickness of model ', ...
           'in well ''%s''.'];
    if refDepth == 0
        msg = [msg, ' BHP reference depth is defaulted to zero. ', ...
               'Consider setting a different value.'];
    else
        msg = [msg, ' Please check refDepth value.'];
    end
    warning('RefDepth:PDropExceedsThickness', msg, wellName)
end

function rr = radiusRep(G, radius, welldir, cells)

if(isfield(G,'nodes'))
   [dx, dy, dz] = cellDims(G, cells);
elseif isfield(G.faces, 'centroids')
   [dx, dy, dz] = cellDimsCG(G, cells);
else
   warning('RepRad:FaceGeomMissing', ...
           'Face geometry missing, will not compute representative radius')
   [dx, dy, dz] = deal(nan(G.cells.num, 1));
end

welldir = lower(welldir);
if numel(welldir) == 1
    welldir = repmat(welldir, numel(cells), 1);
end
re = zeros(size(welldir, 1), 1);

ci = welldir == 'x';
re(ci) = sqrt(dy(ci) .* dz(ci) / pi);

ci = welldir == 'y';
re(ci) = sqrt(dx(ci) .* dz(ci) / pi);

ci = welldir == 'z';
re(ci) = sqrt(dy(ci) .* dx(ci) / pi);

rr = sqrt( re .* radius);

function warn_empty_cell_index_range(nW, varargin)
[opt, ign] = merge_options(struct('Name', ''), varargin{:});   

name = opt.Name;
if isempty(name), name = sprintf('%d', nW + 1); end

warning('Empty:CellIndexSelection', ...
        'Empty Cell Selection in Well ''%s''. Well Ignored.', name);
