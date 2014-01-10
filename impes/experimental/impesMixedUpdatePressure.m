function state = impesMixedUpdatePressure(state, soln, dofPos, G, T, ...
                                          mob, rho, OP, wells, bc)
   flx = state.flux;
   p   = state.pressure;
   if ~isempty(wells),
      is_bhp = strcmpi('bhp', { wells.type });

      flx = [ flx ; vertcat(state.wellSol.flux)              ];
      p   = [ p   ; vertcat(state.wellSol(~is_bhp).pressure) ];
   end

   cf  = G.cells.faces(:,1);

   flx = flx + soln(dofPos(1) : dofPos(1 + 1) - 1);
   p   = p   - soln(dofPos(2) : dofPos(2 + 1) - 1);
   lam =       soln(dofPos(3) : dofPos(3 + 1) - 1);

   dG  = bsxfun(@times, OP.gpot, rho(OP.cellno, :));

   num = bsxfun(@plus, p(OP.cellno), dG);
   Tm  = bsxfun(@times, T, mob(cf, :));
   fp  = accumarray(cf, sum(Tm .* num, 2)) ./ accumarray(cf, sum(Tm, 2));

   NF = ~OP.active;
   if ~isempty(bc),
      di      = strcmpi(bc.type, 'pressure');
      bcf     = bc.face(di);
      fp(bcf) = bc.value(di);

      c       = sum(G.faces.neighbors(bcf, :), 2);
      hf      = OP.hf(bcf, c);
      NF(hf)  = false;
   end

   fp(OP.connno(NF)) = state.facePressure(OP.connno(NF)) + lam;

   i = find(~isfinite(fp));
   if ~isempty(i),
      % Assertion:
      %
      %   This situation occurs only in multi-phase flow and only if *all*
      %   (upwind) phase mobilites are identically zero on the faces 'i' in
      %   which case the above defintion produces fp=0/0 (==NaN).
      %   Moreover, these faces *must* be internal (not boundary).
      %
      %   A typical case in which this may happen is if the interfaces 'i'
      %   represent *sharp* separators between mobile phases of different
      %   mass densities with the heaviest at "the bottom".  Fairly
      %   academic, but nevertheless possible as a result of equilibrium
      %   based intialisation or as a stable state at T_\infty.
      %
      % Use a simple transmissibility weighted average of the cell
      % pressures in this case as there is effectively *no* connection
      % across the interface.  In other words, the exact interface pressure
      % value does not matter, we need only compute something vaguely
      % reasonable in order to subsequently compute A(fp).
      %
      % On the other hand, as there is no connection across the interface,
      % the Jacobian matrix will effectively decouple into sub-systems and
      % we may experience difficulty solving the system of linear eqns.
      %
      N = G.faces.neighbors(i, :);

      % Verify above assertion.
      %
      assert (size(mob, 2) > 1, ...
             ['Non-finite connection pressure in single-phase ', ...
              'flow is highly unexpected.']);

      assert (all(isnan(fp(i))), ...
              'Inifinite connection pressures highly unexpected.');

      assert (all(all(mob(i, :) == 0)), ...
             ['Non-finite connection pressure even if upwind ', ...
              'phase mobilities are non-zero.']);

      assert (all(all(N ~= 0, 2)), ...
              'Unexpected zero (upwind) total mobility on boundary.');

      % Compute connection pressure by means of transmissibility weighted
      % cell pressure values.
      %
      c  = reshape(N, [], 1);
      hf = OP.hf(repmat(i, [2, 1]), c);
      ix = repmat((1 : numel(i)).', [2, 1]);

      fp(i) = accumarray(ix, T(hf) .* p(c)) ./ accumarray(ix, T(hf));
   end

   if ~isempty(wells),
      pw         = zeros([numel(is_bhp), 1]);
      pw(is_bhp) = vertcat(wells(is_bhp).val);
      wflx       = flx((G.faces.num + 1) : end);

      if any(~is_bhp),
         pw(~is_bhp) = p((G.cells.num + 1) : end);
         p           = p(1 : G.cells.num);
      end

      if ~all(pw > 0),
         error('Non-physical pressure value in wells');
      end

      offset = 0;
      for w = 1 : numel(wells),
         nperf = numel(wells(w).cells);

         state.wellSol(w).pressure = pw(w);
         state.wellSol(w).flux     = wflx(offset + (1 : nperf));

         offset = offset + nperf;
      end
   end

   if ~all(p > 0),
      error('Non-physical pressure value in cells');
   end

   state.pressure     = p;
   state.facePressure = fp;

   state.flux = flx(1 : G.faces.num);
end
