function bc = VEM3D_addBC(bc, faces, type, g)

assert(strcmp(type, 'pressure') | strcmp(type, 'flux'));
assert( isa(g, 'function_handle') | numel(g) == 1);

if isempty(bc),
   bc = struct('faces', [], 'type', {{}}, 'func', [], funcPos, 0);
end

if ~isa(g, 'function_handle')
    g = @(X) g*ones(size(X,1),1);
end

nF = numel(faces);

type  = repmat({type}, 1, nF);
g     = repmat({g   }, 1, nF);

face2Func = repmat(funcPos(end) + 1, nF, 1);

bc.faces     = [bc.faces    ; faces(:) ];
bc.type      = [bc.type     ; type     ];
bc.func      = [bc.func     ; g        ];
bc.face2Func = [bc.face2Func; face2Func];

end