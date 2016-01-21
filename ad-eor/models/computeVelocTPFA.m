function veloc = computeVelocTPFA(G, rock, varargin)

% Compute approximation of velocity for TPFA
% Note: The approximation method has flaws at well cells.

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

   % We compute value of product K by outer normal n on each half-faces.
   cellNo = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2).';
   cfn    = G.faces.normals(G.cells.faces(:,1), :);
   dim    = size(G.nodes.coords, 2);
   [K, r, c] = permTensor(rock, dim);
   Kn = K(cellNo, :).*cfn(:, c);
   Kn = reshape(Kn', dim, []);
   Kn = sum(Kn);
   Kn = reshape(Kn, dim, [])';
   vol = G.cells.volumes(cellNo);

   % Define mapping from internal faces to half faces.
   nhf = size(G.cells.faces, 1); % number of half faces
   nf  = G.faces.num;            % number of faces
   nif = nnz(intInx);            % number of internal faces
   nc  = G.cells.num;            % number of cells
   fromIntfacesToFaces = sparse(find(intInx), 1 : nif, 1, nf, nif);
   fromFacesToHalffaces = sparse(1 : nhf, G.cells.faces(:, 1), 1, nhf, nf);
   fromIntfacesToHalffaces = fromFacesToHalffaces*fromIntfacesToFaces;

   hfinvT = fromIntfacesToHalffaces*(1./T);
   wKn = bsxfun(@times, Kn, hfinvT./vol);
   sumHalffaces = sparse(cellNo, 1 : nhf, 1, nc, nhf);
   veloc = cell(dim, 1);
   for i = 1 : dim
      veloc{i} = @(v)(sumHalffaces*(wKn(:, i).*(fromIntfacesToHalffaces*v)));
   end

end