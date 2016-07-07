function D = C2D(C, G)
% Compute the matrix of normalized strain energies (D) from the elasticity
% tensor (C) associated with grid (G).
%
% convert form voit notation to matrix in VEM format of vem formulation

%{ 
   Copyright 2009 - 2014 SINTEF ICT, Applied Mathematics
%} 

   if(G.griddim == 3)
      nlin = 6; 
   else
      assert(G.griddim == 2)
      nlin = 3; 
   end

   factor = [ones(G.griddim, G.griddim),            ones(G.griddim, nlin - G.griddim) * 2; ...
             ones(nlin - G.griddim, G.griddim) * 2, ones(nlin - G.griddim, nlin - G.griddim) * ...
             4]; 
   factor = reshape(factor, 1, []); 
   D = bsxfun(@times, C, factor);
end

