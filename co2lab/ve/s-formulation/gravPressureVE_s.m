function ff = gravPressureVE_s(G, omega)
%Computes innerproduct of (face_centroid - cell_centroid) * g for each face
%
% SYNOPSIS:
%   ff = gravPressureVE_s(g, omega)
%
% DESCRIPTION:
%   This function is an alternate gravity contribution for
%   computePressureRHS which is used with s-formulation VE.
%
% PARAMETERS:
%   G     - Top surface grid as defined by topSurfaceGrid
%
%   omega - Accumulated phase densities \rho_i weighted by fractional flow
%           functions f_i -- i.e., omega = \sum_i \rho_i f_i.  One scalar
%           value for each cell in the discretised reservoir model, G.
%
% RETURNS:
%   ff = omega*(face_centroid - cell_centroid)*g for each face for use in
%        construction of right hand systems for VE models.
% 

%{
#COPYRIGHT#
%}

   g_vec = gravity();   % Must be a 1-by-3 row vector for subsequent code.

   if norm(g_vec) > 0,       
      dim = size(G.nodes.coords,2);

      assert (1 < dim && dim < 4);
      assert (all(size(g_vec) == [1,3]));
      cellno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';     
      
      if (any(strcmp(G.type, 'topSurfaceGrid')))
         % VE - 2D (i.e. g is created using topSurfaceGrid)
         cvec   = (G.faces.z(G.cells.faces(:,1), :) - ...
                   G.cells.z(cellno            , :));
         ff     = omega(cellno) .* cvec*norm(gravity);
      else
         % VE - 3D
         cvec   = G.faces.centroids(G.cells.faces(:,1), :) - ...
                  G.cells.centroids(cellno            , :);
         cellF2D = (1:size(G.cells.faces,1))';
         
         
         hftb = any(bsxfun(@eq,G.cells.faces(cellF2D,2),[5,6]),2);
         % remove contribution from top and bottom, but only for 2D (VE)
         % cells i.e. not for real 3D cells (in the coupled version)         
         cvec(cellF2D(hftb),:) = 0;
         ff     = omega(cellno) .* (cvec * g_vec(1:dim).');
         
      end
   else
      ff     = zeros([size(G.cells.faces,1), 1]);
   end
end
