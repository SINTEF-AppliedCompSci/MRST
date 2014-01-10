g         = computeGeometry(cartGrid([4, 4]));
rock.perm = ones([g.cells.num, 1]);

[BI, pconn, conn] = mex_ip_simple(g, rock);

nconn  = double(diff(pconn));
[i, j] = blockDiagIndex(nconn, nconn);

S = struct('BI', sparse(i, j, BI), 'type', 'hybrid', 'ip', 'ip_simple');
x = initResSol(g, 0);


fluid = initSingleFluid('mu', 1, 'rho', 1);
mu = fluid.properties(x);
s  = fluid.saturation(x);
kr = fluid.relperm(s, x);

mob = sum(bsxfun(@rdivide, kr, mu), 2);

p = mex_partition_ui(double(g.cells.indexMap), g.cartDims, [2, 2]);

[b2c_pos, b2c] = mex_partition_invert(p);

cg = mex_generate_coarsegrid(g, p, 5);
cs = generateCoarseSystem(g, rock, S, cg, mob);

BPsi = basisMatrixHybrid(g, cg, cs);
Psi  = S.BI * BPsi;


cellno = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
a = sortrows([double(p(cellno)), (1:size(g.cells.faces,1))']);
blkcf_pos = cumsum([1; accumarray(a(:,1), 1)]);


blkno = rldecode(1:cg.cells.num, diff(cg.cells.facePos), 2) .';
t = [blkno(cs.activeCellFaces), (1 : numel(cs.activeCellFaces)) .'];
tpos = cumsum([1; accumarray(t(:,1), 1)]);

Psi_vals = [];
for b = 1 : numel(tpos) - 1,
   for cf = reshape(t(tpos(b) : tpos(b+1) - 1, 2), 1, []),
      Psi_vals = [Psi_vals; ...
                  full(Psi(a(blkcf_pos(b) : blkcf_pos(b+1) - 1, 2), cf))];
   end
end


dof                 = zeros([cg.faces.num, 1]);
dof(cs.activeFaces) = 1 : numel(cs.activeFaces);


[cell_ip, Binv] = mex_compute_coarse_contrib(BI, Psi_vals, p, pconn, ...
                                             tpos, diff(blkcf_pos));

[ss, rr, ff, ll] = ...
   mex_schur_comp_symm(Binv, tpos, dof(cg.cells.faces(cs.activeCellFaces)));
