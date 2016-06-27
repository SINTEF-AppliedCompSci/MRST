function [div] = VEM_div(G)
% for no define all the primary operators with out units i.e without
% volum, area, length

%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 
if(G.griddim == 3)
    div = VEM3D_div(G);
    %div = calVemDiv_vec(G);
else
    assert(G.griddim == 2)
    div = VEM2D_div(G);
end
end
%{
old code for 3D
function vdiv = calVemDiv(G)
%ncl = diff(G.cells.nodePos
%svdiv = struct('iind', [], 'jind', [], 'matel', []);
nn = sum(diff(G.faces.nodePos))*G.griddim;
svdiv = struct('iind', zeros(nn, 1), 'jind', zeros(nn, 1), 'matel', zeros(nn, 1))
pos = 1;
for cell = 1:G.cells.num
    ifaces = G.cells.facePos(cell):(G.cells.facePos(cell+1)-1);
    faces = G.cells.faces(ifaces, 1);
    for j = 1:numel(faces)
        iedges = G.faces.edgePos(faces(j)):(G.faces.edgePos(faces(j)+1)-1);
        inodes = G.faces.nodePos(faces(j)):(G.faces.nodePos(faces(j)+1)-1);
        edges = G.faces.edges(iedges);
        nodes = G.faces.nodes(inodes);
        eSign = G.faces.edgeSign(iedges);
        % cheking
        assert(numel(nodes) == numel(edges)); %number
        ienodes = mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1);
        enodes = G.edges.nodes(ienodes);
        % test if nodes of faces and nodes of edges of faces is consist
        %assert(all(unique(enodes, 'sorted') == sort(nodes)));
        % test sign of  edges compeared to orientation of faces
        assert(all(enodes([1:2:numel(enodes)]'+(1-eSign)/2) == nodes))
        nn = numel(nodes);
        assert(nn == 4);
        nvec = G.faces.normals(faces(j), :)*(2*(G.faces.neighbors(faces(j), 1) == cell)-1);        
        matel = repmat(nvec', nn, 1)/(nn);
        jind = (rldecode(nodes, G.griddim)-1)*G.griddim+repmat([1:G.griddim]', numel(nodes), 1);
        iind = repmat(cell, numel(nodes)*G.griddim, 1);
        nn = numel(iind);
        %svdiv.iind = [svdiv.iind;iind];svdiv.jind = [svdiv.jind;jind];svdiv.matel = [svdiv.matel;matel];
        ind = pos:pos+nn-1;
        svdiv.iind(ind) = iind;
        svdiv.jind(ind) = jind;
        svdiv.matel(ind) = matel;
        pos = pos+nn;
    end   
end
vdiv = sparse(svdiv.iind, svdiv.jind, svdiv.matel, G.cells.num, G.nodes.num.*G.griddim);
end
function vdiv = calVemDiv_vec(G)
%ncl = diff(G.cells.nodePos
%svdiv = struct('iind', [], 'jind', [], 'matel', []);
nn = sum(diff(G.faces.nodePos))*G.griddim;

cells = 1:G.cells.num;
ifaces = mcolon(G.cells.facePos(cells), (G.cells.facePos(cells+1)-1));
faces = G.cells.faces(ifaces, 1);
%iedges = mcolon(G.faces.edgePos(faces), (G.faces.edgePos(faces+1)-1));
inodes = mcolon(G.faces.nodePos(faces), (G.faces.nodePos(faces+1)-1));
%edges = G.faces.edges(iedges);
nodes = G.faces.nodes(inodes);
cellfaces = rldecode(cells', diff(G.cells.facePos));
nvec = bsxfun(@times, G.faces.normals(faces, :), (2*(G.faces.neighbors(faces, 1) == cellfaces)-1));
nfl = G.faces.nodePos(faces+1)-G.faces.nodePos(faces);
matel = reshape(rldecode(nvec, nfl)', [], 1);
%facenodes = rldecode(faces, diff(G.faces.nodePos));
jind = mcolon(G.griddim*(nodes-1)+1, G.griddim*(nodes-1)+G.griddim);
iind = rldecode(cellfaces, nfl*3);
%assert(nn == numel(matel(:)));
vdiv = sparse( iind', jind', matel(:), G.cells.num, G.nodes.num.*G.griddim);

end
%}