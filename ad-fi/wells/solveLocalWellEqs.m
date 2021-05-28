function sol = solveLocalWellEqs(W, pBH, qs, p, rho, b, rs, m, sol, maxit)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
