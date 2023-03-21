function geometry=mrstGrid2NG(G)
%Undocumented Utility Function

%{ 
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
%}

geometry = struct(); 
geometry.Nc = G.cells.num;
geometry.Nf = G.faces.num;
geometry.Ncorn = G.nodes.num;
geometry.center = rowMat2cell(G.cells.centroids); 
geometry.corner = rowMat2cell(G.nodes.coords);

[cno, nno, hfno, fno,w]  = createVEMappings(G);
%w=[cells(i), nodes(i), hfaces(i), faces(i)];

%find node to face mapping
n2f=sortrows([w(:,2),w(:,4)]);
[pos,val] = map2Pos(n2f);
geometry.corner_faces = clean_list(pos2list(pos,val));

%find cell to node mapping
c2n=sortrows([w(:,1),w(:,2)]);
[pos,val] = map2Pos(c2n);
geometry.cell_corners = clean_list(pos2list(pos,val));

% find node to cell mapping
n2c=sortrows([w(:,2),w(:,1)]);
[pos,val] = map2Pos(n2c);
geometry.corner_cells = clean_list(pos2list(pos,val));
assert(G.nodes.num==numel(geometry.corner_cells));
%geometry.cell_corners = [];%clean_list(cell_corners);
%geometry.corner_cells = [];%clean_list(corner_cells);


norm_o=bsxfun(@rdivide,G.faces.normals,G.faces.areas);
[geometry.face_norm,geometry.face_cells] =   normAndNeigh(norm_o,G.faces.neighbors);
geometry.cell_faces = pos2list(G.cells.facePos,G.cells.faces(:,1));%clean_list(cell_faces); 
geometry.face_corners = pos2list(G.faces.nodePos,G.faces.nodes);%clean_list(corner_faces); 

%geometry.face_corners = [];%clean_list(face_corners); 
geometry.face_center = rowMat2cell(G.faces.centroids); 
%geometry.face_norm = rowMat2cell(bsxfun(@rdivide,G.faces.normals,G.faces.areas)); 
geometry.face_tangent = [];%x_t; 
geometry.area = rowMat2cell(G.faces.areas); 
geometry.volume = rowMat2cell(G.cells.volumes); 
geometry.dx = [];%dx; 

geometry.Nd = G.griddim;  

geometry.isbnd = [];%face_bcn;
face_bcn = cell(G.faces.num,1); 
[face_bcn{:}] = deal('no');
BC='Dirichlet';
% Orient sign on normal vectors to point from "low" to "high" cell, out otherwise 
for iter1 = 1:G.faces.num
    lfc = geometry.face_cells{iter1};         
    if length(lfc) == 1; 
        face_bcn{iter1} = BC; 
    end    
end
geometry.isbnd =face_bcn;

end
function list=pos2list(pos,val)
 nl=numel(pos)-1;
 list=cell(nl,1);
 for i=1:nl;
     list{i}=sort(val(pos(i):pos(i+1)-1));
     assert(size(list{i},2)==1)
 end
end
function [norm,N] = normAndNeigh(norm_o,N_o)
    [N,j] =  sort(N_o,2);
    sign_cang=2*(N(:,1)==N_o(:,1))-1;
    sign_cang(N_o(:,2)==0)=1;
    norm=bsxfun(@times,norm_o,sign_cang);
    % boundary is always out 
    norm(N_o(:,1)==0,:)=-1.*norm(N_o(:,1)==0,:);
    norm=rowMat2cell(norm);
    list = cell(size(N_o,1),1);
    for i=1:numel(list)
        if(any(N(i,:)==0))
            list{i}=sum(N(i,:));
        else
            list{i}=N(i,:)';
        end        
    end
    N=list;
    
end

function [pos,val]=map2Pos(nn)
 ind=find(diff(nn(:,1))==1);
 pos=[1;ind+1;size(nn,1)+1];
 val=nn(:,2);
end


function cc=rowMat2cell(mat)
nc=size(mat,1);
cc=cell(nc,1);
for i=1:nc;
    cc{i}=mat(i,:)';
end
end

function ilist =  invert_list(list, maxval)

ilist = cell(maxval,1);
for iter1 = 1:length(list)
    for iter2 = list{iter1}'
        ilist{iter2} = [ilist{iter2}; iter1];
    end
end
for iter1 = 1:maxval
    ilist{iter1} = unique(ilist{iter1});
end
end

function list =  clean_list(list)

for iter1 = 1:length(list)
    list{iter1} = unique(list{iter1});
end
end
