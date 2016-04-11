function bc = VEM_addBC(bc, G, f, t, v, varargin)

if isempty(f),
   warning('MRST:addBC', 'Empty list of boundary faces.');
   return;
end

opt = struct('sat', []);
opt = merge_options(opt, varargin{:});
s   = opt.sat;

if isempty(bc),
   bc = struct('face', [], 'type', {{}}, 'value', [], 'sat', []);
end

% Validate boundary condition type.
% Convert to (one-element) cell array if valid.
%
if ischar(t),
   t = repmat({ t }, [1, numel(f)]);
else
   assert(iscellstr(t) && numel(t) == numel(f), ...
      ['Boundary condition type should be either a string or a cell', ...
       'array of strings']);
    t = reshape(t, 1, []);
end

t = lower(t);

assert (all(strcmpi(t, 'pressure') | strcmpi(t, 'flux')), ...
        'Boundary condition type should be either ''pressure'' of ''flux''');

% Validate saturation input.
%  - If ISEMPTY(bc.sat), then 's' may be any numeric array (including []).
%  - Otherwise, 's' must be numeric and have the same number of columns as
%    the existing 'bc.sat' array.
%
assert (isnumeric(s));
assert (isempty(bc.sat) || (size(bc.sat,2) == size(s,2)));

% Verify that boundary condition is not already set
bc_given = false(max([f(:); bc.face]), 1);
bc_given(bc.face) = true;
assert(~any(bc_given(f)), ...
   'New boundary condition overlaps with conditions already in bc-struct.');

nf = numel(f);

% Expand single-element saturations and BC values to cover all faces.
%
if size(s,1) == 1, s = s(ones([nf, 1]), :); end

if ~isa(v,'function_handle')
    if size(v,1) == 1, v = v'; end
    if numel(v) == 1, v = repmat(v(ones([nf, 1])),1,3);
    else v = repmat(v,1,3);
    end
else
    nodeNum = mcolon(G.faces.nodePos(f),G.faces.nodePos(f+1)-1);
    nodes = G.faces.nodes(nodeNum);
    nodes = reshape(nodes,2,[])';
    v = [v(G.nodes.coords(nodes(:,1),:)), ...
         v(G.nodes.coords(nodes(:,2),:)), ...
         v(G.faces.centroids(f,:))      ];
end

% Verify that v and s are same length as faces (or s empty).
assert (size(v,1)     ==     nf  );
assert (any(size(s,1) == [0, nf]));

% Boundary conditions structurally verified.  Append to existing structure.
bc.face  = [bc.face ; f(:)];
bc.type  = [bc.type , t];
bc.value = [bc.value; v];
bc.sat   = [bc.sat  ; s];
