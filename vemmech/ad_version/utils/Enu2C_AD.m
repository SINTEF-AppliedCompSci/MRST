function C = Enu2C_AD(E, nu, G, ixnames)

   assert(G.griddim == 3) % @@ only 3D supported for now
   if ~exist('ixnames', 'var')
      ixnames = {'i', 'j', 'c'}; % default index names
   end
   
   
   % scalar factor
   T1 = SparseMultiArray(E ./ (1 + nu) ./ (1 - 2 * nu));
   
   % upper-left block
   T2 = SparseMultiArray(nu * ones(9, 1), ...
                     vertcat(repmat([1,2,3], 1, 3), ...
                             [ones(1,3), 2*ones(1,3), 3*ones(1,3)])', ...
                     ixnames(1:2));
   
   T3 = SparseMultiArray((1 - 2 * nu) * ones(3, 1), ...
                     [1, 1; 2, 2; 3, 3], ixnames(1:2));
   
   % lower-right block
   T4 = SparseMultiArray((1 - 2 * nu)/2 * ones(3, 1), ...
                     [4, 4; 5, 5; 6, 6], ixnames(1:2));
   
   % tensor representing cell-space
   Tcells = SparseMultiArray([], (1:G.cells.num)', ixnames(3));
   
   % resulting C tensor
   
   C = (T1 * (T2 + T3 + T4)) * Tcells;
   
end
