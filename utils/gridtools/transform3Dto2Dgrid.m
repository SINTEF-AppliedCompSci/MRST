function G2 = transform3Dto2Dgrid(G)
% Transforms a 3D grid into a 2D grid.
% 3D grid must be constant in either x-, y- or z-direction.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


  if (numel(unique(G.cells.centroids(:,1)))==1)
    dir=1;
  elseif (numel(unique(G.cells.centroids(:,2)))==1)
    dir=2;
  elseif (numel(unique(G.cells.centroids(:,3)))==1)
    dir=3;
  end

  % G2.cells
  G2.cells.num         = G.cells.num;
  G2.cells.numFaces    = G.cells.numFaces;
  G2.cells.numFaces(:) = 4;
  if(isfield(G.cells,'indexMap'))
  	G2.cells.indexMap = G.cells.indexMap;
  end

  if dir==1
    dx = max(G.nodes.coords(:,1))-min(G.nodes.coords(:,1));
    G2.cells.volumes   = G.cells.volumes./dx;
    G2.cells.centroids = G.cells.centroids(:,2:3);
    G2.cartDims        = ([numel(unique(G.cells.centroids(:,2))),...
                           numel(unique(G.cells.centroids(:,3)))]);

    % G2.nodes
    nodes           = find(G.nodes.coords(:,1)== min(G.nodes.coords(:,1)));
    G2.nodes.coords = G.nodes.coords(nodes,[2,3]);
    G2.nodes.num    = size(G2.nodes.coords,1);

    % G2.faces
    faces = true(G.faces.num,1);
    faces(G.faces.normals(:,1)~=0) = false;
    I = find(faces);
    G2.faces.num         = numel(I);
    G2.faces.numNodes    = G.faces.numNodes(I);
    G2.faces.numNodes(:) = 2;
    G2.faces.neighbors   = G.faces.neighbors(I,:);
    G2.faces.areas       = G.faces.areas(I)./dx;
    G2.faces.normals     = G.faces.normals(I,2:3);
    G2.faces.centroids   = G.faces.centroids(I,2:3);

    % G2.faces.nodes:
    facenodes    = [rldecode((1 : G.faces.num) .', ...
                             double(G.faces.numNodes)), G.faces.nodes];
    actfaces     = ismember(facenodes(:,1),find(faces));
    facenodes    = facenodes(actfaces,2);
    G2.faces.nodes = facenodes(ismember(facenodes,nodes));
    G2.faces.nodes = compressFaces(G2.faces.nodes);
    facenodes    = [rldecode((1 : G2.faces.num) .', ...
                             double(G2.faces.numNodes)), G2.faces.nodes];
    y_faces      = find(G2.faces.normals(:,2));
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,2));
    y_facenodes  = facenodes(ind_fn,:);
    y_face       = 0;
    for i = 1:size(y_faces)
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;
      node       = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,1), 'descend');
     G2.faces.nodes(ind_fn(y_face)) = G2.faces.nodes(ind_fn(y_face(ind)));
    end
    y_faces      = find(G2.faces.normals(:,1));
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,1));
    y_facenodes  = facenodes(ind_fn,:);
    y_face       = 0;
    for i = 1:size(y_faces)
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;
      node       = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,2), 'ascend');
      G2.faces.nodes(ind_fn(y_face)) = G2.faces.nodes(ind_fn(y_face(ind)));
    end

    % G2.cells.faces
    faces        = true(size(G.cells.faces,1),1);
    faces(G.cells.faces(:,2)==1) = false;
    faces(G.cells.faces(:,2)==2) = false;
    G2.cells.faces = G.cells.faces(faces==true,:);
    G2.cells.faces = compressFaces(G2.cells.faces);
    G2.cells.faces(G2.cells.faces(:,2)==3,2)=1;
    G2.cells.faces(G2.cells.faces(:,2)==4,2)=2;
    G2.cells.faces(G2.cells.faces(:,2)==5,2)=3;
    G2.cells.faces(G2.cells.faces(:,2)==6,2)=4;

  elseif dir==2
    dy = max(G.nodes.coords(:,2))-min(G.nodes.coords(:,2));
    G2.cells.volumes   = G.cells.volumes./dy;
    G2.cells.centroids = G.cells.centroids(:,[1,3]);
    G2.cartDims        = ([numel(unique(G.cells.centroids(:,1))),...
                           numel(unique(G.cells.centroids(:,3)))]);

    % G2.nodes
    nodes           = find(G.nodes.coords(:,2)== min(G.nodes.coords(:,2)));
    G2.nodes.coords = G.nodes.coords(nodes,[1,3]);
    G2.nodes.num    = size(G2.nodes.coords,1);

    % G2.faces
    faces = true(G.faces.num,1);
    faces(G.faces.normals(:,2)~=0) = false;
    I = find(faces);
    G2.faces.num         = numel(I);
    G2.faces.numNodes    = G.faces.numNodes(I);
    G2.faces.numNodes(:) = 2;
    G2.faces.neighbors   = G.faces.neighbors(I,:);
    G2.faces.areas       = G.faces.areas(I)./dy;
    G2.faces.normals     = G.faces.normals(I,[1,3]);
    G2.faces.centroids   = G.faces.centroids(I,[1,3]);

    % G2.faces.nodes:
    facenodes    = [rldecode((1 : G.faces.num) .', ...
                             double(G.faces.numNodes)), G.faces.nodes];
    actfaces     = ismember(facenodes(:,1),find(faces));
    facenodes    = facenodes(actfaces,2);
    G2.faces.nodes = facenodes(ismember(facenodes,nodes));
    G2.faces.nodes = compressFaces(G2.faces.nodes);
    facenodes    = [rldecode((1 : G2.faces.num) .', ...
                             double(G2.faces.numNodes)), G2.faces.nodes];
    y_faces      = find(G2.faces.normals(:,2));
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,2));
    y_facenodes  = facenodes(ind_fn,:);
    y_face       = 0;
    for i = 1:size(y_faces)
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;
      node       = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,1), 'descend');
      G2.faces.nodes(ind_fn(y_face)) = G2.faces.nodes(ind_fn(y_face(ind)));
    end

    % G2.cells.faces
    faces        = true(size(G.cells.faces,1),1);
    faces(G.cells.faces(:,2)==3) = false;
    faces(G.cells.faces(:,2)==4) = false;
    G2.cells.faces = G.cells.faces(faces==true,:);
    G2.cells.faces = compressFaces(G2.cells.faces);
    G2.cells.faces(G2.cells.faces(:,2)==5,2) = 3;
    G2.cells.faces(G2.cells.faces(:,2)==6,2) = 4;

elseif dir==3
    dz = max(G.nodes.coords(:,3))-min(G.nodes.coords(:,3));

    G2.cells.volumes   = G.cells.volumes./dz;
    G2.cells.centroids = G.cells.centroids(:,1:2);
    G2.cartDims        = ([numel(unique(G.cells.centroids(:,1))),...
                           numel(unique(G.cells.centroids(:,2)))]);

    % G2.nodes
    nodes           = find(G.nodes.coords(:,3)== min(G.nodes.coords(:,3)));
    G2.nodes.coords = G.nodes.coords(nodes,1:2);
    G2.nodes.num    = size(G2.nodes.coords,1);

    % G2.faces
    faces = true(G.faces.num,1);
    faces(G.faces.normals(:,3)~=0) = false;
    I = find(faces);
    G2.faces.num         = numel(I);
    G2.faces.numNodes    = G.faces.numNodes(I);
    G2.faces.numNodes(:) = 2;
    G2.faces.neighbors   = G.faces.neighbors(I,:);
    G2.faces.areas       = G.faces.areas(I)./dz;
    G2.faces.normals     = G.faces.normals(I,1:2);
    G2.faces.centroids   = G.faces.centroids(I,1:2);

    % G2.faces.nodes:
    facenodes    = [rldecode((1 : G.faces.num) .', ...
                             double(G.faces.numNodes)), G.faces.nodes];
    actfaces     = ismember(facenodes(:,1),find(faces));
    facenodes    = facenodes(actfaces,2);
    G2.faces.nodes = facenodes(ismember(facenodes,nodes));
    facenodes    = [rldecode((1 : G2.faces.num) .', ...
                                     double(G2.faces.numNodes)), G2.faces.nodes];
    y_faces      = find(G2.faces.normals(:,2));
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,2));
    y_facenodes  = facenodes(ind_fn,:);
    y_face       = 0;
    for i=1:size(y_faces)
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;
      node = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,1), 'descend');
      G2.faces.nodes(ind_fn(y_face)) = G2.faces.nodes(ind_fn(y_face(ind)));
    end

    % G2.cells.faces
    faces        = true(size(G.cells.faces,1),1);
    faces(G.cells.faces(:,2)==5) = false;
    faces(G.cells.faces(:,2)==6) = false;
    G2.cells.faces = G.cells.faces(faces==true,:);
    G2.cells.faces = compressFaces(G2.cells.faces);
  end
  G = [G.type, mfilename];
end

function cf = compressFaces(cf)
  f=cf(:,1);
  active = find(accumarray(f, 1) > 0);
  compr  = zeros([max(f), 1]);
  compr(active) = 1 : numel(active);
  f = compr(f);
  cf(:,1)=f;
end
