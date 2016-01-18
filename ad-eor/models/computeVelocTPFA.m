function veloc = computeVelocTPFA(G, rock, varargin)

   opt = struct('verbose', false, 'operators', [] );
   opt = merge_options(opt, varargin{:});

   if isempty(operator)
      operator = setupOperatorsTPFA(G, rock);
   end
   T = operators.T;
   N = operators.N;
   intInx = operators.internalConn;
   K = rock.perm;

   % We compute value of product K by outer normal n on each half-faces. Obtain array of
   % dimension intInx-by-dim
   cellNo = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2).';
   sgn    = 2*double(G.faces.neighbors(G.cells.faces(:,1), 1) == cellNo) - 1;
   % cfn  = cell-face normal = *outward* face normals on each cell.
   cfn    = bsxfun(@times, G.faces.normals(G.cells.faces(:,1), :), sgn);
   dim    = size(G.nodes.coords, 2);
   [K, r, c] = permTensor(rock, dim);
   Kn = zeros(size(cellNo, 1), dim);
   Kn = K(cellNo,:).*cfn(:, c);
   Kn = reshape(Kn',3, []);
   Kn = sum(Kn);
   Kn = reshape(Kn, 3, [])';
   Kn = Kn(intInx);
   Af = G.faces.areas(G.cells.faces(:, 1));
   vol = G.cells.volume(CellNo);
   wKn = 1./vol(intInx).*bsxfun(@times, Kn, Af./T);
   sumFaceOper = sparse(cellNo(intInx), 1:numel(intInx), 1, G.cells.num, numel(intInx));
   veloc = @(v)(sumFaceOper*bsxfun(wKn, v)); 
   
end