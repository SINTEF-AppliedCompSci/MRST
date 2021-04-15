function h=plotSlice(G, data_in, slice_ind, dim)
% Plot Cartesian slices of cell data on faces
%
% SYNOPSIS:
%    h = plotSlice(G, data_inn, slice_ind, dim)
%
% PARAMETERS:
%   G         - Grid data structure.
%
%   data_in   - Cell data to be plotted
%
%   slice_ind - Index of slice
%
%   dim       - Dimension of slice
%
% RETURNS:
%    h - Handle to resulting `patch` object.  The patch object is added to
%        the current `axes` object. 

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


[IX,IY,IZ]=ndgrid(1:G.cartDims(1),1:G.cartDims(2),1:G.cartDims(3));
if(dim==1)
   cellind=sub2ind(G.cartDims,IX(slice_ind,:,:),IY(slice_ind,:,:),IZ(slice_ind,:,:));
   face_tag=1;
elseif(dim==2)
   cellind=sub2ind(G.cartDims,IX(:,slice_ind,:),IY(:,slice_ind,:),IZ(:,slice_ind,:));
   face_tag=3;
elseif(dim==3)
   cellind=sub2ind(G.cartDims,IX(:,:,slice_ind),IY(:,:,slice_ind),IZ(:,:,slice_ind));
   face_tag=5;
else
   error('not valid dim');
end
cellind=cart2active(G,cellind(:));
cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
mask_c = false([G.cells.num, 1]);
mask_f = false([6, 1]);
mask_c(cellind) = true;
mask_f(face_tag) = true;
mask = mask_c(cellno) & mask_f(G.cells.faces(:,2));
faces = G.cells.faces(mask, 1);
data  = data_in(cellind);
fn=accumarray(cellno, double(mask_f(G.cells.faces(:,2))));
data=rldecode(data,fn(cellind));
h=plotFaces(G, faces, data);axis tight;shading flat;
vv=zeros(1,3);
vv(dim)=1;
view(vv);
end

