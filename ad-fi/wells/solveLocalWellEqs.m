function sol = solveLocalWellEqs(W, pBH, qs, p, rho, b, rs, m, sol, maxit)

   [sol, sol] = checkLims(sol, pBH, qs{:});                            %#ok

   consts = {double(p), cellfun(@double, rho, 'UniformOutput', false), ...
             cellfun(@double, b, 'UniformOutput', false), double(rs), ...
             cellfun(@double, m, 'UniformOutput', false)};

   X = cell(4,1);
   [X{1}, X{2}, X{3}, X{4}] = ...
      initVariablesADI(double(pBH)  , double(qs{1}), ...
                       double(qs{2}), double(qs{3}));

   % We have to compute residual before starting  the loop in the case we
   % have a solution.  The solution may not be unique (for example if null
   % flux, any alpha is solution).
   %
   [eqs, sol, sol] = ...
      getWellContributionsBO(W, X{1}, X(2:4), consts{:}, sol);         %#ok

   res = norm(double(vertcat(eqs{:})), inf);
   res_diff = inf; its=0;

   while (res_diff > 1e-5) && (its < maxit)
      alpha0 = vertcat(sol.alpha);

      [eqs, sol, sol] = ...
         getWellContributionsBO(W, X{1}, X(2:4), consts{:}, sol);      %#ok

      dx = SolveEqsADI(eqs,[]);
      for k = 1:numel(X)
         if k==1
            dx{k} = sign(dx{k}).*min(abs(dx{k}), 10*barsa);
         end

         X{k}.val = X{k}.val + dx{k};
      end

      res_diff = norm(vertcat(sol.alpha)-alpha0, inf);
      res = norm(double(vertcat(eqs{:})), inf);

      [withinLims, sol] = checkLims(sol, X{:});

      if any(~withinLims)
         res_diff = inf;
      end

      its = its +1;
   end

   fprintf('Well its: %d\n', its);
   for k = 1:numel(sol)
      sol(k).pressure = X{1}.val(k);
      sol(k).qWs      = X{2}.val(k);
      sol(k).qOs      = X{3}.val(k);
      sol(k).qGs      = X{4}.val(k);
   end
end
