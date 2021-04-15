function [resSol, report] = newtonRaphson2ph(resSol, tf, F, Jac, update, opt)
%Solve non-linear equation F(s)=0 using Newton-Raphson method.
%
% SYNOPSIS:
%   resSol          = newtonRaphson2ph(resSol, tf, F, Jac, update, opt)
%   [resSol,report] = newtonRaphson2ph(resSol, tf, F, Jac, update, opt)
%
% DESCRIPTION:
%   Solve general nonlinear equation F(s)=0 using a Newton-Raphson
%   iteration with user-specified UPDATE scheme to adjust iterate.
%   Specifically, the k-th iteration is defined schematically as
%
%      ds^k = J \ F
%      s^k  = update(s^k-1, ds^k, ...)
%
% PARAMETERS:
%   resSol - the reservoir state
%   tf     - final solution time
%   F      - discretization on residual form
%   Jac    - Jacobian of F
%   update - user-specified scheme to adjust iterate
%   opt    - structure with optional parameters, see function
%            'implicitTransport' for documentation
%
% RETURNS:
%   resSol - updated solution
%   report - structure reporting success (true/false), number of
%       iterations, vasted iterations, substeps, failed steps, residual and
%       convergence rate

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


   report = struct('success',           true,...
                   'iterations',        0, ...
                   'vasted_iterations', 0, ...
                   'sub_steps',         0, ...
                   'failed_steps',      0, ...
                   'residual',          inf, ...
                   'convergence_rate',  0);

   show_history = opt.verbose && opt.show_convergence;

   mints   = pow2(tf, -opt.tsref);
   dispif(show_history, 'Time interval (s)     iter  relax   residual          rate\n');
   dispif(show_history, [repmat('-', [1, 70]), '\n']);

   [t, dt, dtprev, count] = deal(0.0, tf, -tf, 0);

   %-----------------------------------------------------------------------
   % Main time loop (ideally, just a single step) -------------------------
   %
   while t < tf && dt >= mints,
      dt = min(dt, tf - t);

      %--------------------------------------------------------------------
      % Outer controlling loop (sub step size &c) -------------------------
      %
      redo_newton = true;
      while redo_newton,
         if isempty(opt.init_state),
            sn_0 = resSol;
            sn = resSol; sn.s(:) = min(1,sn.s+0.05);
         else
            sn_0 = resSol;
            sn = resSol; sn.s(:) = opt.init_state.s(:);
         end

         %-----------------------------------------------------------------
         % Initialise NR algorithm ----------------------------------------
         %
         res     = F(sn, sn_0, dt);
         err     = norm(res(:), inf);
         E       = err;
         nwtfail = err > opt.nltol;
         linfail = false;
         it      = 0;
         report.residual = err;

         %-----------------------------------------------------------------
         % Execute inner NR algorithm -------------------------------------
         %
         while nwtfail && ~linfail && it < opt.maxnewt,

            J  =  Jac(sn, sn_0, dt);
            ns =  size(res, 2);
            ds = -reshape(opt.LinSolve(J, reshape(res', [], 1)), ns, [])';
            [sn, res, alph, linfail] = update(sn, sn_0, ds, dt, err);

            it      = it + 1;
            err     = norm(res(:), inf);

            % Append current error
            E    = [E; err];                                           %#ok

            nwtfail = err > opt.nltol;
            dispif(show_history, ...
                   ' [%5.1e, %5.1e]:%4d  %6.2f    %11.5e\t %9.2f\n', ...
                   t, t + dt,it, alph, err, computeConvergenceRate(E));
            report.residual = norm(res, inf);
         end

         %-----------------------------------------------------------------
         % Determine if NR succeeded or not -------------------------------
         %
         count = count - 1;
         if nwtfail,
            dispif(show_history, ...
               '-------------------- Reducing step -------------------\n');

            if count > 0 && dtprev > 0,
               dt = dtprev;
            else
               dt = dt / (1.5 + 0.5);
            end
            count = 5;
            report.failed_steps      = report.failed_steps      + 1;
            report.vasted_iterations = report.vasted_iterations + it;
         else
            redo_newton = false;
            t = t + dt;  dtprev = -dt; %intentionally invalid dtprev.
            if it == 0,
               % ?
               dispif(show_history, '\t%4d\t %s\t %5.5e', it, '-', err);
               dispif(show_history, '\t NB: err <= ntol.');
            end
            if t<tf && it <= 5  && count < 0, % arbitrary threshold
               dispif(show_history, ...
               '-------------------- Increasing step -----------------\n');

               dtprev = dt;
               dt     = min((1 + 0.5) * dt, tf - t);
               count  = 5;
            end
            if (t<tf),
               dispif(show_history, ...
                  '-------------------- Next substep -----------------------\n');
            end
         end
      end

      % Time step [t, t+dt] was successful
      resSol = sn;

      if isfield(resSol, 'extSat'),
         % Save extremal values of saturation for hysteresis purposes
         % First column contains minima, second maxima
           resSol.extSat(:,1) = min(resSol.s(:,1), resSol.extSat(:,1));
           resSol.extSat(:,2) = max(resSol.s(:,1), resSol.extSat(:,2));
      end

      report.iterations = report.iterations + it;
      report.sub_steps  = report.sub_steps  + 1;
   end

   report.success = ~(t<tf) || any(isnan(resSol.s(:,1)));
   report.convergence_rate = computeAverageConvergenceRate(E);
   dispif(show_history, '\n');

   % Update oil saturation if applicable.
   if size(resSol.s,2) > 1,
      resSol.s(:,2) = 1 - resSol.s(:,1);
   end
   report.str = dispReport(report, dt, mints);
end

%--------------------------------------------------------------------------

function R = computeConvergenceRate(E)
% Convergence rate of iteration based on
%
%    |E(n)|            |E(n-1)|
%   --------   = A and --------   = A.
%   |E(n-1)|^R         |E(n-2)|^R
%
%
% Let  a = log|E(n)|, b = log|E(n-1)| and c = log|E(n-2)|, then
%
%    [a]   [1 b][A]
%    [ ] = [   ][ ]
%    [b]   [1 c][R]

   if numel(E) < 3,
      R = nan;
      return;
   end

   a = log(E(end));
   b = log(E(end-1));
   c = log(E(end-2));

   y = [1, b; 1, c]\[a;b];
   R = y(2);
end

%--------------------------------------------------------------------------

function R = computeAverageConvergenceRate(E)
% Least squares approximation to convergence rate of iteration based on
%
%   |E(n+1)|
%   -------- = A
%    |E(n)|^R

    x = log(E);

    N    =  numel(x)-1;
    xbm  = sum(x(1:end-1))/N;
    xbp  = sum(x(2:end))/N;
    x2bm = sum(x(1:end-1).*x(1:end-1))/N;
    x2b  = sum(x(1:end-1).*x(2:end))/N;

    den = x2bm - xbm*xbm;
    if abs(den) > 0,
       R = -(xbm*xbp - x2b) / den;
    else
       R = NaN;
    end
end

%--------------------------------------------------------------------------

function str = dispReport(r, dt, mints)
   str = sprintf([...
      ' Iterations        :%4d            Wasted iterations : %d\n',...
      ' Sub steps         :%4d            Failed steps      : %d\n',...
      ' Final residual    :%11.2e     Convergence rate  : %3.1f\n',...
      ], ...
      r.iterations, r.vasted_iterations, r.sub_steps, r.failed_steps, ...
      r.residual, r.convergence_rate);

   if ~r.success, %?
      str = sprintf('implicitTransport: FAILED due to timestep %g < %g.\n', ...
                dt, mints);
   end
end
