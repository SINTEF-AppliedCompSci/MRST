function bc = VEM3D_addBC(bc, faces, type, g)
%   Adds boundary condition to (new or existing) BC object. Usage differs
%   from that of the MRST function addBC in the following way:
%
%       - Input values v is a function handle, and understood to be the
%         normal derivative of the unknown function. Can also be scalars as
%         in addBC, and is the understood to be these function values at
%         the faces.
%
%   See MRST function addBC for details and copyright info.
%-----------------------------------------------------------------ØSK-2016-

assert(strcmp(type, 'pressure') | strcmp(type, 'flux'));
assert( isa(g, 'function_handle') | numel(g) == 1);

if isempty(bc),
   bc = struct('faces', [], 'type', {{}}, 'func', [], 'face2Func', []);
end

if ~isa(g, 'function_handle')
    g = @(X) g*ones(size(X,1),1);
end

nF = numel(faces);

type  = repmat({type}, 1, nF);
g     = repmat({g   }, 1, nF);

bc.faces     = [bc.faces    ; faces(:) ];
bc.type      = [bc.type     , type     ];
bc.func      = [bc.func     , g        ];

end