function bc = VEM3D_addBC(bc, faces, type, g)

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

if isempty(bc.face2Func)
    face2Func = ones(nF,1);
else
    face2Func = repmat(bc.face2Func(end) + 1, nF, 1);
end

bc.faces     = [bc.faces    ; faces(:) ];
bc.type      = [bc.type     , type     ];
bc.func      = [bc.func     , g        ];

end