function bc_nfvm = convertBC(G, bc)

% Convert a std mrst bc structure to bc_nfvm.
% Fill the bc_nfvm by setting default homogeneous Neumann conditions,
% then copy from bc_std. We should set BC explicitly on all boundary
% faces, which means that bc_nfvm members should be of size
% numel(boundaryFaces(G)).
bf = boundaryFaces(G);
nf = numel(bf);
bc_nfvm.face = zeros(nf, 1);
bc_nfvm.type = repmat({'flux'}, [nf, 1]);
bc_nfvm.value = repmat({@(x) 0},[nf, 1]);

for i = 1:numel(bc.face)
    faceno = bc.face(i);
    bc_nfvm.face(faceno) = faceno;
    bc_nfvm.type(faceno) = bc.type(i);

    % Can't do this with @(x) since the value is not evaluated
    %bc_nfvm.value(faceno) = {@(x) bc_std.value(i)};
    
    if strcmpi(bc.type(i), 'flux')
        bc_nfvm.value(faceno) = {@(x) 1}; % flux value is 1
    else
        bc_nfvm.value(faceno) = {@(x) 1}; % pressure value is 0
        [i,faceno,bc_nfvm.type(faceno),bc_nfvm.value(faceno)]
    end
end

end
