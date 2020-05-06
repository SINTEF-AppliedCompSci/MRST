function bi = blockInverter(opt)
   if ~ischar(opt.invertBlocks)
      dispif(opt.verbose, ...
            ['Unsupported option value of type ''%s'' in ', ...
             'option ''invertBlocks''. Reset to default ' , ...
             '(''matlab'')\n'], class(opt.invertBlocks));
      opt.invertBlocks = 'matlab';
   end

   switch lower(opt.invertBlocks)
      case {'matlab', 'm', 'builtin'}
         bi = @invertDiagonalBlocks;
      case {'mex', 'c', 'accelerated'}
         bi = @invertDiagonalBlocksMex;
      otherwise
         dispif(opt.verbose, ...
               ['Unsupported value ''%s'' in option ', ...
                '''invertBlocks''.\nMust be one of ''matlab'' or ', ...
                '''mex''.\nReset to default (''matlab'').'], ...
                opt.invertBlocks);

         bi = @invertDiagonalBlocks;
   end
end

% Matlab code calling mex function invertSmallMatrices.c
function iA = invertDiagonalBlocksMex(A, sz)
   sz     = int32(sz);
   blocks = matrixBlocksFromSparse(A, sz);
   iA     = blockDiagMatrix(invv(blocks, sz), sz);
end

%--------------------------------------------------------------------------

% Pure Matlab code using inv
function iA = invertDiagonalBlocks(A, sz)
   V = zeros([sum(sz .^ 2), 1]);
   [p1, p2] = deal(0);

   for b = 1 : numel(sz)
      n  = sz(b);
      n2 = n * n;
      i  = p1 + (1 : n);

      V(p2 + (1 : n2)) = inv(full(A(i, i)));

      p1 = p1 + n;
      p2 = p2 + n2;
   end

   iA = blockDiagMatrix(V, sz);
end

%--------------------------------------------------------------------------

function A = blockDiagMatrix(V, sz)
   [I, J] = blockDiagIndex(sz);
   A      = sparse(I, J, V);
end
