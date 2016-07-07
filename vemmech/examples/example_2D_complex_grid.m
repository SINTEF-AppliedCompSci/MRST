%% A short example solving linear elasticity on a complex grids
% 
%% Define parameters
opt=struct('L',[1 1],...
    'cartDims',[4 4],...
    'grid_type','square',...
    'disturb',0.02,... %parameter for disturbing grid
    'E',4e7,...  %youngs modolo
    'nu',0.44);% poiso ratio

%% make a mixed type of gird 
G= squareGrid(opt.cartDims,opt.L,'grid_type','mixed4','disturb',opt.disturb);
G=removeCells(G,[140:151,233:235,249:250,93,117:118]')
G=createAugmentedGrid(G);
G=computeGeometry(G);
clf,plotGrid(G);




%% Find sides of domain
% find sides normal sides is not normaly defined on complex grids
[Lx,Ly] = deal(opt.L(1),opt.L(2));
assert(G.griddim==2);
x=[0,Lx];
for i=1:2
    clear find
    faces=find(abs(G.faces.centroids(:,1)-x(i))<eps);
    bc{i}=addBC([], faces, 'pressure', 0);
    bc{i}= rmfield(bc{i},'type');
    bc{i}= rmfield(bc{i},'sat');
end
y=[0,Ly];
for i=1:2
    faces=find(abs(G.faces.centroids(:,2)-y(i))<eps);
    bc{i+2}=addBC([], faces, 'pressure', 0);
    bc{i+2}= rmfield(bc{i+2},'type');
    bc{i+2}= rmfield(bc{i+2},'sat');
end

for i=1:4
    inodes=mcolon(G.faces.nodePos(bc{i}.face),G.faces.nodePos(bc{i}.face+1)-1);
    nodes=unique(G.faces.nodes(inodes));
    disp_bc=struct('nodes',nodes,...
        'uu',0,...
        'faces',bc{i}.face,...
        'uu_face',0,...
        'mask',true(numel(nodes),G.griddim));
    bc{i}.el_bc=struct('disp_bc', disp_bc,'force_bc',[]);
end

%% Set up loading terms and boundary conditions
% define load as gravity
density=3000*0;% kg/m^3
grav=0;% gravity
load=@(x) -(grav*density)*repmat([0,1],size(x,1),1);
% set boundary dispace ment function to zeros
bcdisp=@(x) x*0.0;


% set up boundary conditions for each side
clear bc_el_sides
% set direclet boundary conditions at selected sides
bc_el_sides{1}=bc{1};% x side is fixed
bc_el_sides{2}=bc{2};% y side is fixed
bc_el_sides{3}=[];% bottom  is free
bc_el_sides{4}=[];% top is free


% collect boundary conditions
% logical cartesian boundary conditions to proper boundary structure
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
% define forces at boundary
% find midpoint face set all force corresponding to the "weight fo the ting
%a t limited area
sigma=opt.L(2)/10;force=50*barsa;
%face_force =@(x) force*exp(-(((x(:,1)-opt.L(1)/2))./sigma).^2)+100*barsa;
face_force =@(x) force*sign(x(:,1)-opt.L(1)/2)+100*barsa;
faces=bc{4}.face;
% make force boundary structure NB force is in units Pa/m^3
force_bc=struct('faces',faces,'force',bsxfun(@times,face_force(G.faces.centroids(faces,:)),[0 -1]));
% final structure fo boundary conditions
el_bc=struct('disp_bc',disp_bc,'force_bc',force_bc);


% define rock parameters
Ev=repmat(opt.E,G.cells.num,1);nuv=repmat(opt.nu*0+0.4,G.cells.num,1);
C=Enu2C(Ev,nuv,G);

%% solve system with profiler on
% do profiling
profile off;profile on
lsolve=@mldivide;
[uu,extra]=VEM_linElast(G,C,el_bc,load,'linsolve',lsolve);
profile off;profile on

%% Plot displacement in y direction
plotopts={'EdgeAlpha',0.0,'EdgeColor','none'};
plotopts1={'EdgeAlpha',0.01};
clf,plotNodeData(G,uu(:,2),plotopts{:}),colorbar(),

%% Plot displacement in x direction
clf,plotNodeData(G,uu(:,1),plotopts{:}),colorbar(),

%% plot the deformed grid
fac=1;
clf,plotGridDeformed(G,uu*fac)
%%
%% calulate divergence and plot
vdiv=VEM_div(G);
mdiv=vdiv*reshape(uu',[],1)./G.cells.volumes;
clf,plotCellDataDeformed(G,mdiv,uu*fac,plotopts1{:}),colorbar()


%% calulate stress and strain
op=extra;
stress=reshape(op.D*op.WC'*op.assemb'*reshape(uu',[],1),3,[])';
% Note the ordering is not traditional voits and the of diagnol terms have
% have a different




