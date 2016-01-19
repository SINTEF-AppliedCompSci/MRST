function veloc = computeVelocTPFA(G, rock, varargin)

   opt = struct('verbose', false, 'operators', [] );
   opt = merge_options(opt, varargin{:});

   if isempty(opt.operators)
      operators = setupOperatorsTPFA(G, rock);
   else
      operators = opt.operators;
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
   Kn = K(cellNo,:).*cfn(:, c);
   Kn = reshape(Kn', dim, []);
   Kn = sum(Kn);
   Kn = reshape(Kn, dim, [])';
   Kn = Kn(intInx, :);
   Af = G.faces.areas(G.cells.faces(:, 1));
   Af = Af(intInx);
   vol = G.cells.volumes(cellNo);
   wKn = bsxfun(@times, Kn, 1./vol(intInx).*Af./T);
   sumFaceOper = sparse(cellNo(intInx), 1 : nnz(intInx), 1, G.cells.num, nnz(intInx));
   dim = G.griddim;
   veloc = cell(dim, 1);
   for i = 1 : dim
      veloc{i} = @(v)(sumFaceOper*(wKn(:, i).*v));
   end

end