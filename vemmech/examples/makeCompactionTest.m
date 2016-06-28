function [el_bc,load] =makeCompactionTest(G,opt)
sides={'XMin','XMax','YMin','YMax','ZMin','ZMax'}
for j=1:G.griddim
    Lmax=max(G.faces.centroids(:,j));
    Lmin=min(G.faces.centroids(:,j));
    x = [Lmin, Lmax];
    for i = 1:2      
        if(false)
            faces = find(abs(G.faces.centroids(:, j)-x(i))<tol);
        else
            mside=sides{2*(j-1)+i};
            tmp=pside([],G,mside,0);
            faces=tmp.face;
        end
        bc{i+(j-1)*2} = addBC([], faces, 'pressure', 0);
    end
end

% finde node of the differens sides and prepare elastisity boundary
% conditions

for i=1:2*G.griddim
    inodes=mcolon(G.faces.nodePos(bc{i}.face),G.faces.nodePos(bc{i}.face+1)-1);
    nodes=unique(G.faces.nodes(inodes));
    disp_bc=struct('nodes',nodes,...
                    'uu',0,...
                    'faces',bc{i}.face,...
                    'uu_face',0,...
                    'mask',true(numel(nodes),G.griddim));
    bc{i}.el_bc=struct('disp_bc', disp_bc,'force_bc',[]);
end
%% define load as gravity
density=3000;% kg/m^3
grav=10;% gravity 
gdir=zeros(1,G.griddim);gdir(end)=1;
if(opt.islinear)
    origo=mean(G.cells.centroids,1);
    fac=(max(G.cells.centroids(:,G.griddim))-min(G.cells.centroids(:,G.griddim)))*1000;
    bcdisp=@(x) bsxfun(@minus,bsxfun(@times,x,[0 0 1]),origo)./fac;
else
  bcdisp=@(x) x*0.0;  
end
if(opt.gravity_load)
    load=@(x) -(grav*density)*repmat(gdir,size(x,1),1);
    
else    
    load=@(x) -(0*density)*repmat(gdir,size(x,1),1);
end

clear bc_el_sides
% set direclet boundary conditions at selected sides
    % on left side nod displace ment in x direction only, this is done by
    % mask
bc_el_sides=bc;
if(~opt.hanging || opt.islinear)
for i=1:2*(G.griddim-1);
    bc_el_sides{i}.el_bc.disp_bc.mask(:,G.griddim)=false;
end

for i=2*G.griddim;
    bc_el_sides{i}.el_bc.disp_bc.mask(:,1:(G.griddim-1))=false;
end
for i=2*G.griddim-1;
    if(opt.free_top)
        bc_el_sides{i}=[]; 
    else    
       bc_el_sides{i}.el_bc.disp_bc.mask(:,1:G.griddim-1)=false; 
    end   
end
else
  for i=2*G.griddim-2:2*G.griddim;
    bc_el_sides{i}=[];
  end  
end
% collect bounary conditions
nodes=[];
faces=[];
mask=[];

for i=1:numel(bc)
    if(~isempty(bc_el_sides{i}))
        nodes=[nodes;bc_el_sides{i}.el_bc.disp_bc.nodes];%#ok
        faces=[faces;bc_el_sides{i}.el_bc.disp_bc.faces];%#ok
        mask=[mask;bc_el_sides{i}.el_bc.disp_bc.mask];%#ok
    end
end

disp_node=bcdisp(G.nodes.coords(nodes,:));
disp_faces=bcdisp(G.faces.centroids(faces,:));
disp_bc=struct('nodes',nodes,'uu',disp_node,'faces',faces,'uu_face',disp_faces,'mask',mask);
if(opt.top_load)
    fvec=zeros(1,G.griddim);
    fvec(G.griddim)=1;
    H=max(G.nodes.coords(:,G.griddim))-min(G.nodes.coords(:,G.griddim));
    face_force =@(x) H*10*3000*repmat(-fvec,size(x,1),1);
    faces=bc{2*G.griddim-1}.face;
    % make force boundary structure NB force is in units Pa/m^3
    force_bc=struct('faces',faces,'force',face_force(G.faces.centroids(faces,:)));
else
   force_bc=[]; 
end
% final structure fo boundary conditions
%}
%el_bc=struct('disp_bc',disp_bc,'force_bc',force_bc);

el_bc=struct('disp_bc',disp_bc,'force_bc',force_bc);
end