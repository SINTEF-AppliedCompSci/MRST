function h=plotSlice(g,data_inn,slice_ind,dim)
%Plot Cartesian slices of cell data on faces
%
% SYNOPSIS:
%    h=plotSlice(g,data_inn,slice_ind,dim)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   data   - cell data to be plotted
%
%   slice_ind - index of slice
%
%   dim - dimension of slice
%
%  RETURNS:
%    h - Handle to resulting PATCH object.  The patch object is added to the
%       current AXES object.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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


[IX,IY,IZ]=ndgrid(1:g.cartDims(1),1:g.cartDims(2),1:g.cartDims(3));
if(dim==1)
   cellind=sub2ind(g.cartDims,IX(slice_ind,:,:),IY(slice_ind,:,:),IZ(slice_ind,:,:));
   face_tag=1;
elseif(dim==2)
   cellind=sub2ind(g.cartDims,IX(:,slice_ind,:),IY(:,slice_ind,:),IZ(:,slice_ind,:));
   face_tag=3;
elseif(dim==3)
   cellind=sub2ind(g.cartDims,IX(:,:,slice_ind),IY(:,:,slice_ind),IZ(:,:,slice_ind));
   face_tag=5;
else
   error('not valid dim');
end
cellind=cart2active(g,cellind(:));
cellno = rldecode(1:g.cells.num, diff(g.cells.facePos), 2)';
mask_c = false([g.cells.num, 1]);
mask_f = false([6, 1]);
mask_c(cellind) = true;
mask_f(face_tag) = true;
mask = mask_c(cellno) & mask_f(g.cells.faces(:,2));
faces = g.cells.faces(mask, 1);
data  = data_inn(cellind);
fn=accumarray(cellno, double(mask_f(g.cells.faces(:,2))));
data=rldecode(data,fn(cellind));
h=plotFaces(g, faces, data);axis tight;shading flat;
vv=zeros(1,3);
vv(dim)=1;
view(vv);
end

