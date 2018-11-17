function davg=calcdavg(nodes,planenormal,planepoint,tol)
% calculates davg for fracture-matrix CI calculation. Refer to Moinfar et
% al (2013) for details.

% METHODOLOGY - We generate a cartesian grid for the matrix cell in
% question. Then for every subgridcell, we calculate the point-plane
% distance, multiply it by the subgridcell volume, and add it to a sum.
% This will be the numerator that we need to calculate davg.
%
%
% We divide the cell into 1000 subgridcells and compute davg in a discrete
% manner.
%
% Alternatively: We use a while loop with error calculation (SLOW)
% We will need a starting point for generating davg values. Then while
% successive iterations differ by more than a tolerance, we refine the
% subgrid.

origin=min(nodes); % location of origin of cell
planepoint=planepoint-origin; % move planepoint
lx=max(nodes(:,1))-min(nodes(:,1));
ly=max(nodes(:,2))-min(nodes(:,2));
lz=max(nodes(:,3))-min(nodes(:,3));

% Alternative method
% starting subgrid - make subgridcells approximately cubic
% minside=min([lx,ly,lz]);
% gridscale=minside/10;
% nx=ceil(lx/gridscale);
% ny=ceil(ly/gridscale);
% nz=ceil(lz/gridscale);

nx=10;ny=10;nz=10;
celldim=[nx,ny,nz];
physdim=[lx,ly,lz];
G=cartGrid(celldim,physdim);
G=computeGeometry(G);
N = G.cells.num;

cellvol=sum(G.cells.volumes); % only need to calculate this once

normdist=pointplanedistance(G.cells.centroids,repmat(planenormal,N,1),repmat(planepoint,N,1));
davg = sum(normdist.*G.cells.volumes)/cellvol;

% davg=numericaldavg(G,cellvol,planenormal,planepoint); % first approximation

% Alternative method
% error=2*tol;
% while error>tol
%     davg_old=davg; % save previous davg
%     
%     celldim=celldim*2;
%     G=cartGrid(celldim,physdim); G=computeGeometry(G);
%     davg=numericaldavg(G,cellvol,planenormal,planepoint); % next approximation
%     
%     error=abs(davg-davg_old); % compare davg
% end


end

% function d=numericaldavg(G,cellvol,planenormal,planepoint)
% numcells=G.cells.num;
% sum=0;
% for i=1:numcells
%     normdist=pointplanedistance(G.cells.centroids(i,:),planenormal,planepoint);
%     sum=sum+normdist*G.cells.volumes(i);
% end
% 
% 
% d=sum/cellvol; 
% end


