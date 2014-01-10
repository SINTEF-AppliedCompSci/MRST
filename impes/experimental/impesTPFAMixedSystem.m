function [A, b, dofPos, fmob, rho] = ...
      impesTPFAMixedSystem(state, G, T, fluid, dt, pv, varargin)
   opt = struct('bc', [], 'wells', []);
   opt = merge_options(opt, varargin{:});

   OP = tpfaConnectionOperators(G, opt.wells, size(state.z,2));

   wdp = impesTPFADefaultWellModel(state, opt.wells, fluid);
   [p, z, fp] = impesAssembleStateVars(state, opt.bc, opt.wells, wdp);
   [rho, rho] = fluid.pvt(p, z);                                       %#ok

   [cmob, dcmob] = impesComputeMobility(state, fluid, ...
                                        opt.bc, opt.wells, wdp);

   [fmob, fz] = tpfaUpwindStateVars(G, p, z, rho, cmob, dcmob, ...
                                    opt.bc, opt.wells, OP);

   press = struct('cell', p, 'face', fp);
   mass  = struct('cell', z, 'face', fz);

   trans = 1 ./ accumarray(OP.connno, 1 ./ T, [G.faces.num, 1]);
   vt    = state.flux;

   [WI, wcellno, wconnno, wother, wsgn, pickp, dFw, dCw, nw] = ...
      well_connections(opt.wells, G.cells.num, OP);

   if ~isempty(opt.wells),
      vt    = [ vt    ; vertcat(state.wellSol.flux) ];
      trans = [ trans ; WI                          ];
   end

   trans = trans .* sum(fmob, 2);

   NF = ~OP.active;

   [dF, dC] = deal([]);
   if ~isempty(opt.bc),
      if ~all(strcmpi('pressure', opt.bc.type)),
         error('Boundary conditions other than Dirichlet not supported.');
      end

      dF = opt.bc.face;
      c  = sum(G.faces.neighbors(dF,:), 2);

      hf = OP.hf(dF, c);

      NF(hf) = false;
      dC     = - OP.sgn(hf) .* fp(dF);
   end

   numNF = sum(NF);

   IC = [ OP.connno ; wconnno ; wconnno(pickp) ];
   JC = [ OP.cellno ; wcellno ; wother( pickp) ];
   VC = [ OP.sgn    ; wsgn    ; -wsgn(  pickp) ];

   B = sparse(1 : numel(trans), 1 : numel(trans), 1 ./ trans);
   C = sparse(IC, JC, VC, numel(trans), G.cells.num + nw);
   D = sparse(OP.connno(NF), 1 : numNF, ...
              OP.sgn   (NF), numel(trans), numNF);

   [F, Fp, Fv] = impesMixedContinuityResidual(pv, dt, fluid, fmob, vt, ...
                                              press, mass, opt.wells, OP);

   NullCt = sparse(size(F,1), numNF);
   NullDt = sparse(numNF, size(F,1) + numNF);

   A = [B  ,  C , D     ; ...
        Fv , -Fp, NullCt; ...
        D.', NullDt     ];

   b1                    = zeros([numel(trans), 1]);
   b1(dF)                = dC;
   b1(G.faces.num + dFw) = dCw;

   sp = state.pressure;
   if ~isempty(opt.wells),
      i  = ~strcmpi('bhp', { opt.wells.type });
      sp = [sp ; vertcat(state.wellSol(i).pressure) ];
   end
   b1 = b1 - B*vt + C*sp - D*state.facePressure(OP.connno(NF));
   b3 = - D' * vt;

   b  = [ b1 ; -F ; b3 ];

   numFlux = numel(trans);   % Number of flux equations
   numRes  = numel(F);   % Number of residual equations (nc + sum(~is_bhp))

   dofPos = cumsum([1; numFlux; numRes; numNF]);
end

%--------------------------------------------------------------------------

function [WI, cellno, connno, other, sgn, pickp, dF, dC, nw] = ...
      well_connections(W, nc, OP)
   if isempty(W),

      [WI, cellno, connno, other, sgn, pickp, dF, dC] = deal([]);
      nw = 0;

   else

      is_bhp = reshape(strcmpi('bhp', { W.type }), [], 1);

      %if all(is_bhp),

       %  [WI, cellno, connno, other, sgn, pickp, dF, dC] = deal([]);

      %else

         nperf  = reshape(cellfun('prodofsize', { W.cells }), [], 1);
         pickw  = ~is_bhp;
         pickp  = rldecode(pickw, nperf);

         WI     = vertcat(W.WI);
         cellno = OP.wcellno;
         connno = OP.wconnno;
         other  = rldecode(nc + (1 : sum(pickw)), nperf(pickw), 2) .';
         sgn    = OP.wsgn;

         dF     = find(~pickp);
         dC     = rldecode(vertcat(W(is_bhp).val), nperf(is_bhp));

      %end

      nw = numel(W) - sum(is_bhp);
   end
end
