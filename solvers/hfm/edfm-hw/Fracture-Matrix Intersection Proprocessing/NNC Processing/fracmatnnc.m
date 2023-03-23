function [cells,CI,area,type]=fracmatnnc(F,Gm,tol)
% INDIVIDUAL FRACTURE-MATRIX NNC GENERATOR - This function takes in the
% fracture grid of an individual fracture and the matrix grid. It then
% generates outputs cells (which contains global indices in each row for
% every frac-mat nnc), CI (which contains the CI for each frac-mat nnc,
% corresponding to the cells array) and area (which contains the area for
% each intersection, corresponding to the cells array).

% METHODOLOGY - We retrieve the matrix cells which the fracture F
% intersects. Then for every matrix cell, we calculate d_avg and then
% prescribe a CI accordingly

% Prelims
mcells=F.matrix_connection.cells;
area=F.matrix_connection.area;
typenum=F.matrix_connection.type;
Nm=length(mcells);
fcellstart=F.cells.start;
Nf=F.cells.num;
assert((Nf==Nm)||(Nm==2*Nf),'Number of intersected matrices don"t make sense. Check code.')

% First output
cells=[mcells,repmat(((1:Nf)+fcellstart-1)',Nm/Nf,1)];

% node data.
% [fn,fpos]=gridCellNodes(F,1:Nf);
[cn,cpos]=gridCellNodes(Gm,1:Gm.cells.num);
% fcellnodes_par=cell(Nm,1);
mcellnodes_par=cell(Nm,1);
for j=1:Nm
%     i=mod(j,Nf);
%     fcellnodeind=fn(fpos(i):(fpos(i+1)-1));
%     Nfcellnodes=size(fcellnodeind,1); % this is always even    
%     fcellnodes=F.nodes2D.coords(fcellnodeind(1:Nfcellnodes/2),:); % Nx3 matrix, each row containing xyz coords of nodes
% 
%     fcellnodes_par(j)={fcellnodes};
    
    mcellnodeind=cn(cpos(mcells(j)):(cpos(mcells(j)+1)-1));
    mcellnodes=Gm.nodes.coords(mcellnodeind,:); % 8x3 matrix, each row containing xyz coords of nodes
    mcellnodes_par(j)={mcellnodes};
end



%% Calculate CI

CI=ones(Nm,1); % preallocate
planepoint=F.points(1,:);
planenormal=getnormal(F.points,tol);
parfor i=1:Nm 
    mcellnodes=mcellnodes_par{i}; % we can do this because first Nm elements in mcellnodes_par follows the original arrangement
    davg=calcdavg(mcellnodes,planenormal,planepoint,tol);
    xarea=area(i);
    
    [xCI,~]=calcfracmatCI(mcellnodes,xarea,planenormal,planepoint,davg,tol);
    CI(i)=xCI;   
end

type=cell(size(typenum));
for i=1:size(typenum,1)
    switch typenum(i)
        case 1
            type{i}='fracmat boundary';
        case 2
            type{i}='fracmat interior';
    end
end

end


