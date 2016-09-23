function bc = addBCFunc(bc, f, t, g)
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

assert(strcmp(t, 'pressure') | strcmp(t, 'flux'));
assert( isa(g, 'function_handle') | numel(g) == 1);

if isempty(bc),
   bc = struct('face', [], 'type', {{}}, 'func', []);
end

if ~isa(g, 'function_handle')
    g = @(X) g*ones(size(X,1),1);
end

nf = numel(f);

t  = repmat({t}, 1, nf);
g  = repmat({g}, 1, nf);

bc.face = [bc.face; f(:)];
bc.type = [bc.type, t   ];
bc.func = [bc.func, g   ];

end