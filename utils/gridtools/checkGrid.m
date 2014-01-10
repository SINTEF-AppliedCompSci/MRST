function ok=checkGrid(G)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ok=true;
internal=sum(G.faces.neighbors~=0,2)>1;
c=G.cells.centroids(G.faces.neighbors(internal,2),:)-G.cells.centroids(G.faces.neighbors(internal,1),:);
n=G.faces.normals(internal,:);
if ~all(sum(n.*c,2)>0)
    ok=false;
   disp('Something wrong between neigbours and sgn of face normal')
end
cellno=rldecode([1:G.cells.num]',diff(G.cells.facePos));
cc=G.faces.centroids(G.cells.faces(:,1),:)-G.cells.centroids(cellno,:);
ncc=bsxfun(@times,G.faces.normals(G.cells.faces(:,1),:),(1-2*(G.faces.neighbors(G.cells.faces(:,1),2)==cellno)));
if ~all(sum(ncc.*cc,2)>0)
    ok=false;
   disp('Something wrong between cellface normals')
end
%{
vv=volumeByGaussGreens(G);
if ~all(vv>0)
    ok=false;
   disp('Greans gauss volumes wrong')
end
%}
if (~all(G.cells.volumes>0))
%end
    ok = false;
end
if (~all(G.faces.areas>0))
%end
ok =false;
end