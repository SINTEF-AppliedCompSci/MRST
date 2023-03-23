function G=mrstGridWithFullMappings(G)
% Add all mappings to mrst grid. This is needed for VEM of finte element
% type methods. 
%
% SYNOPSIS:
%   G=mrstGridWithFullMappings(G)
%   
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
% OPTIONAL PARAMETERS:
%
%   'hingenodes'  - A struct with fields 'faces' and 'nodes'.  A hinge node
%                   is an extra center node for a face, that is used to
%                   triangulate the face geometry.  For each face number F
%                   in 'faces' there is a row in 'nodes' which holds the
%                   node coordinate for the hinge node belonging to face F.
%
%
% COMMENTS:
%   NB at the moment may mappings which are not needed is made. The mapping
%   from local nodes to global is not calcualted and need to be made localy
%   if needed. Forexample by a sparse matrix.
%   The methods make the assumtion of sorted edges and faces   
%
% SEE ALSO:
%   grid_structure. calulateGeometryCalc, primaryMimeticOp, CC_linElast,
%   VEM_linElast
%

%G=sortEdges(G);

%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 
if(G.griddim==3)
    ind1=mcolon(G.faces.nodePos(1:end-1),G.faces.nodePos(2:end)-1,1);
    ind2=mcolon(G.faces.nodePos(1:end-1)+1,G.faces.nodePos(2:end),1);
    ind2(G.faces.nodePos(2:end)-1)=G.faces.nodePos(1:end-1);
    eh=[G.faces.nodes(ind1),G.faces.nodes(ind2)];
    ehs=sort(eh,2);
    [e,~,j]=unique(ehs,'rows','sorted');
    assert(all(all(e(j,:)==ehs)))
    
    
    G.faces.edges=j;
    G.faces.edgePos=[1;cumsum(diff(G.faces.nodePos))+1];
    G.faces.edgeSign=2*(e(j,1)==eh(:,1))-1;
    
    % make ting si pos format
    %eo=e;
    e=e';ne=size(e',1);
    G.edges.nodes=e(:);
    nPos=cumsum(repmat(2,ne,1));
    G.edges.nodePos=[1;nPos+1];
    G.edges.num=ne;
    %%
    we=[rldecode([1:ne]',diff(G.edges.nodePos)),G.edges.nodes(:)];%#ok
    n2e=unique(we,'rows','sorted');
    [pos,val] = map2Pos(n2e);
    G.nodes.edgePos=pos;
    G.nodes.edges=val;
else
    %G=sortEdges(G);
    %warning('Hope for sorted edges')
    %ifaces=mcolon(G.cells.facePos(1:end-1),G.cells.facePos(2:end));
    %faces=G.cells.faces(ifaces);
    faces=G.cells.faces(:,1);
    inodes=mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1);
    assert(2*numel(faces)==numel(inodes))
    cell=rldecode([1:G.cells.num]',diff(G.cells.facePos));%#ok
    sign=2*(cell==G.faces.neighbors(faces,1))-1;
    nf=[1:2:numel(inodes)]'-(sign-1)/2;%1:2:numel(inodes)-(sign-1)/2;
    nodes=G.faces.nodes(inodes(nf));
    G.cells.nodes=nodes;
    G.cells.nodePos=G.cells.facePos;
end
w =createGridMappings(G);

%e2n=sparse(e(1:2:end),e(2:2:end),1:ne,G.nodes.num,G.nodes)...
%    +sparse(e(2:2:end),e(1:2:end),1:ne,G.nodes.num,G.nodes.num)


% find face2edge map


%find node to face mapping
n2f=unique([w(:,2),w(:,4)],'rows','sorted');
[pos,val] = map2Pos(n2f);
G.nodes.facePos=pos;
G.nodes.faces=val;


%find cell to node mapping
%c2n=sortrows([w(:,1),w(:,2)]);
if(G.griddim==3)
    c2n = unique([w(:,1),w(:,2)],'rows','sorted');
    [pos,val] = map2Pos(c2n);
    G.cells.nodePos=pos;
    G.cells.nodes=val;
end

% find node to cell mapping
%n2c=sortrows([w(:,2),w(:,1)]);
n2c=unique([w(:,2),w(:,1)],'rows','sorted');
[pos,val] = map2Pos(n2c);
G.nodes.cellPos=pos;
G.nodes.cells=val;



end

function [pos,val]=map2Pos(nn)
ind=find(diff(nn(:,1))==1);
pos=[1;ind+1;size(nn,1)+1];
val=nn(:,2);
end

