function [G,new_grdecl]=makeMatchingGridFromGrdel(grdecl)
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


%grdecl=simpleGrdecl( [5 6 3])
% shange shape of coord and zcorn for easier indexing
[xyz,zcorn]= grdeclXYZ(grdecl);
d=grdecl.cartDims;
new_cartDim=[d(1:2)+2,d(3)];
% extend z corn by moving and flipping cells form on side
% to an other
%new_zcorn=nan(2*(d+[2,2,0]));
new_zcorn=nan(2*(new_cartDim));
new_zcorn(3:end-2,3:end-2,:)=zcorn;
% flip in xdirection
new_zcorn(1:2,3:end-2,:)=zcorn(end-1:end,:,:);
new_zcorn(end-1:end,3:end-2,:)=zcorn(1:2,:,:);
% flip in ydirection
new_zcorn(3:end-2,1:2,:)=zcorn(:,end-1:end,:);
new_zcorn(3:end-2,end-1:end,:)=zcorn(:,1:2,:);

% extend xyz
new_xyz=nan(6,new_cartDim(1)+1,new_cartDim(2)+1);
new_xyz(:,2:end-1,2:end-1)=xyz;
% extend x direction
new_xyz(:,1,2:end-1)=xyz(:,1,:);
new_xyz(:,end,2:end-1)=xyz(:,end,:);
new_xyz([1,4],1,2:end-1)=xyz([1,4],1,:)-1;
new_xyz([1,4],end,2:end-1)=xyz([1,4],end,:)+1;

% extend y direction
new_xyz(:,2:end-1,1)=xyz(:,:,1);
new_xyz(:,2:end-1,end)=xyz(:,:,end);
new_xyz([2,5],2:end-1,1)=xyz([2,5],:,1)-1;
new_xyz([2,5],2:end-1,end)=xyz([2,5],:,end)+1;
% actnum corners with nan
ACTNUM=ones(new_cartDim);
ACTNUM(1,1,:)=0;
ACTNUM(1,end,:)=0;
ACTNUM(end,end,:)=0;
ACTNUM(end,1,:)=0;
% define final domain
mydomain=zeros(new_cartDim);
mydomain(2:end-1,2:end-1,:)=1;
new_grdecl=struct('cartDims', new_cartDim,...
                  'COORD'  , new_xyz(:),...
                  'ZCORN'   , new_zcorn(:),...
                  'ACTNUM',int32(ACTNUM(:)));
% make grid with faults on outer domain of original grid
% introduced by oposite sides
G=processGRDECL(new_grdecl);
% remove cells introduced to make correct faces on outer boundary
ind=find(mydomain(:)==1);
cells=zeros(prod(G.cartDims),1);
cells(ind)=1;
cells=cells(G.cells.indexMap);
rm_cells=find(cells==0);
G=removeCells(G,rm_cells);
end

