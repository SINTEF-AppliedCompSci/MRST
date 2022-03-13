function [v, u, history, status] = optimizeSR1(u0, f, varargin)
% 
% Unconstrained optimization based on a quasi-Newton approach refered to as
% the SR1-method, combined with a trust-region approach to determine step
% direction and distance.   See Chapter 6 of "Numerical Optimization" by
% Nocedal and Wright.
% 
% SYNOPSIS:
%   function [v, u, history, status] = optimizeSR1(u0, f, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   u0       - Initial guess of the vector of variables
%   f        - Function to minimize
%   varargin - Various optional arguments with default values
%              - 'B_init' - starting point for quasi-Hessian
%              - 'B_scale - scale of quasi-hessian (only relevant if 'B_init'
%                           is not explicitly provided
%              - 'delta' - initial trust region radius 
%              - 'eta' - parameter for when to update 
%              - 'epsilon' - convergence tolerance
%              - 'funval_tol' - function change  tolerance
%              - 'r' - threshold for when to update quasi-Hessian
%              - 'maxIt' - maximum number of search directions
%              - 'backup_file' - file to backup (save) temporary result along
%                                the way
%              - 'plotEvolution' - if true, a progress window will appear and
%                                  stay updated during simulation
%
% RETURNS:
%   v       - optimal value identified for f
%   u       - optimal vector of variables
%   history - history information from the optimizatoin procedure
%   status  - 1: convergence; -1: non-convergence; 2: unable to make further
%             progress 
%
% SEE ALSO:
% unitBoxBFGS

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

% NB: this function _minimizes_, not _mazimizes_

   opt.B_init = []; % starting point for quasi-Hessian
   opt.B_scale = 1; % scale of quasi-hessian (if B_init not set)
   opt.delta = 1; % trust region radius
   opt.eta = 1e-3; % parameter for when to update
   opt.epsilon = 1e-8; % convergence tolerance
   opt.funval_tol = 1e-8; % function change tolerance
   opt.r = 1e-3; % threshold for when to update quasi-Hessian
   opt.maxIt = 100; % max number of search directions
   opt.backup_file = []; % file to back up results along the way
   opt.plotEvolution = true;
   opt.rat_lim = 0.75;  % parameter controlling when to adjust trust region
   opt.delta_fac = 0.8; % parameter controlling when to adjust trust region
   
   opt = merge_options(opt, varargin{:});
   
   % set initial guess for quasi-Hessian matrix, and radius for trust region
   delta = opt.delta;
   B = opt.B_scale * eye(numel(u0));
   if ~isempty(opt.B_init)
      B = opt.B_init;
   end
   
   % check if there is an existing backup to use
   if ~isempty(opt.backup_file) && ...
          ( exist(opt.backup_file) == 2  || exist([opt.backup_file, '.mat']) )
   
      str = input('Continue from existing backup? Y/N [Y]: ', 's');
      if isempty(str)
         str = 'Y';
      end
      if strcmpi(str, 'Y')
         fprintf('Loading from existing backup.\n');
         
         backup = load(opt.backup_file);
         u = backup.history.u{end};
      
      else
         fprintf('Ignoring existing backup.  Starting from scratch.\n');
          u = u0;
      end
   else
      u = u0; % starting point for search
   end

   % compute starting point
   [v0, g0] = f(u); g0 = g0(:);
   
   history = gatherInfo([], v0, u0, norm(g0), g0, nan, nan, B);
   
   % main loop
   status = 1;
   iter_count = 0;
   mini_iter = 1;
   skipped_last_update = false;
   if norm(g0) <= opt.epsilon
      v = v0;
      g = g0;
      return;
   end
   while norm(g0) > opt.epsilon 
      if iter_count > opt.maxIt
         status = -1; % signal non-convergence
         return;
      end
      
      [s, Bstar] = compute_step(g0, B, delta);
       
      [v, g] = f(u + s);
      
      dg = g(:) - g0;                         % change in gradient 
      act_red = v0 - v;                       % actual fct. reduction
      pre_red = -(g0' * s + 0.5 * s' * Bstar * s); % predicted fct. reduction
      
      % check if we shall update the current point u
      ratio = act_red / pre_red;
      if (ratio > opt.eta)
         u = u + s;

         [v0, g0] = deal(v, g(:)); % previously calculated values
         iter_count = iter_count + 1;
         % register point in history and (optionally) plot
         history = gatherInfo(history, v, u, norm(g, inf), g(:), mini_iter, nan, B);
         history.B = B;
         if opt.plotEvolution
            plotInfo(10, history);
         end
         if ~isempty(opt.backup_file)
            save(opt.backup_file, 'history');
         end
         if abs(act_red) < opt.funval_tol
            break; % stop criteria on function reduction reached
         end
         mini_iter = 1;
      else
         if abs(act_red) < opt.funval_tol && mini_iter > 4
            % unable to make progress
            status = 2;
            break;
         end
         mini_iter = mini_iter + 1;
      end
         
      % adjust trust region if necessary
      if ratio > opt.rat_lim
         if norm(s) > delta * opt.delta_fac
            delta = delta * 2;
         end
      elseif ratio < 0.1
         
         if (norm(s) < delta) && skipped_last_update
            % fast-forward shrinking of trust region
            delta = norm(s)/2;
         else
            delta = delta / 2;
         end
      end
      
      % update quasi-hessian if applicable  
      tmp = dg(:) - B * s;
      if abs(s' * tmp) >= opt.r * norm(s) * norm(tmp)
         B = B + (tmp * tmp') / (s' * tmp);
         skipped_last_update = false;
      else
         skipped_last_update = true;
      end
      
      % fprintf(['iteration is: %i.  Norm of B is %f.  min(eig) is: %f; max(eig) ' ...
      %          'is %f; cond is: %f\n\n'], iter_count, norm(B), min(eig(B)), ...
      %         max(eig(B)), cond(B));
      
   end
end

% ----------------------------------------------------------------------------
function [e1, e2, Bstar] = determine_subspace(g, B)

   % check whether B is positive  - if not, work on a modified version
   Bstar = B;
   lmin = min(eig(Bstar));
   if lmin < 0;
      % @@ factor 1.5 could be refined (should theoretically be between 1 and
      % 2, see page 77 of Nocedal: "Numerical Optimization"
      factor = 1.5;
      Bstar = Bstar + abs(lmin) * factor * eye(size(Bstar)); 
   end
   
   g = full(g); % avoid sparse vectors (not handled by 'orth' below)
   v1 = g(:);
   v2 = Bstar\g(:); 
   if any(isinf(v2))
      v2 = v1; % Bstar was singular
   end
   
   q = orth([v1, v2]);
   e1 = q(:,1);
   if size(q, 2) == 2
      e2 = q(:,2);
   else 
      e2 = nan; % degenerate subspace
   end
   
end

% ----------------------------------------------------------------------------
function [s, Bstar] = compute_step(g, B, delta)
% * see page 76-77 in Nocedal: "Numerical optimization"

   % determine subspace, and possibly modified B*, for which we search the
   % minimum of g'p + 1/2 p'B*p
   [e1, e2, Bstar] = determine_subspace(g, B);

   p_int = [];
   s = -Bstar\g(:);
   if norm(s) < delta
      % the extremum was found inside the trust region.  
      if min(eig(B)) > 0
         % we know this is a minimum and not a maximum or saddle point
         return;
      else
         % this could be a maximum or saddle point.  Search boundary of trust
         % region too and return the best point
         p_int = s;
      end
   end
   
   % Minimum is located outside trust region.  We must constrain search to
   % boundary of this region.
   if isnan(e2)
      % degenerate subspace, use gradient descent
      s = - g(:) / norm(g(:)) * delta;
      return;
   end
   
   % we seek p on the form: p(theta) = delta (e1 cos(theta) + e2 sin(theta))
   
   p = @(theta) delta * (e1 * cos(theta) + e2 * sin(theta));
   dp = @(theta) delta * (-e1 * sin(theta) + e2 * cos(theta));
   fun = @(p) g(:)'*p + 1/2 * p' * Bstar * p;

   %fun = @(theta) g(:)'*p(theta) + 1/2 * p(theta)' * B * p(theta);
   dfun = @(theta) (g(:)' + p(theta)' * Bstar) * dp(theta);
   
   % the derivative of fun should have two zeros, corresponding to maximum
   % and minimum.  Choose the minimum.  This can be tricky, as the curve may
   % have multiple minima, each very close

   % @@ The following code is a workaround and ought to be reimplemented in a
   % more proper way.
   samples = 200;
   dt = 2*pi/(samples-2);
   theta = linspace(-dt, 2*pi+dt, samples);
   val = zeros(samples,1);
   for i = 1:samples
      val(i) = fun(p(theta(i)));
   end
   der = val(2:end) - val(1:end-1);
   dleft = der(1:end-1);
   dright = der(2:end);
   min_ixs = dleft <=0 & dright > 0; % indices close to minima
   
   theta_min_candidates = theta(find(min_ixs)+1);
   
   t_exacts = arrayfun(@(t) fzero(dfun, [t-2*dt, t+2*dt]), theta_min_candidates);
   f_exacts = zeros(size(t_exacts));
   for i = 1:numel(t_exacts)
      f_exacts(i) = fun(p(t_exacts(i))); % convert angle to function value
   end
   [~, ix] = min(f_exacts);
   theta_opt = t_exacts(ix);
   
   if ~isempty(p_int) && fun(p_int) < fun(p(theta_opt))
      s = p_int;
   else
      s = p(theta_opt);
   end
end

%--------------------------------------------------------------------------
function hst = gatherInfo(hst, val, u, pg, grad, lsit, lsfl, hess)   
% obj.val | contr | norm proj.grad | ls-step | ls-its | ls-flag | hessian 
if isempty(hst)
    hst = struct('val', val, 'u', {{u}}, 'pg', pg, ...
                 'grad', grad, 'lsit', lsit, 'lsfl', lsfl, ...
                 'hess', {hess});
else
    hst.val   = [hst.val  , val  ];
    hst.u     = [hst.u    , {u}  ];
    hst.pg    = [hst.pg   , pg   ];
    hst.grad  = [hst.grad , grad ];
    hst.lsit  = [hst.lsit , lsit ];
    hst.lsfl  = [hst.lsfl , lsfl ];
    hst.hess  = [hst.hess ,{hess}];
end
end

%--------------------------------------------------------------------------      
function [] = plotInfo(fig, hst)
if ~ishandle(fig)
    figure(fig)
else
    % Avoid stealing focus if figure already exists
    set(0, 'CurrentFigure', fig);
end
xt = 0:(numel(hst.val)-1);
ch = [0, hst.val(2:end)-hst.val(1:end-1)];
subplot(5,1,1), plot(xt, hst.val, '.-','LineWidth', 2, 'MarkerSize', 20), title('Objective');
subplot(5,1,2), semilogy(xt,hst.pg,'.-','LineWidth', 2, 'MarkerSize', 20), title('Gradient norm');
subplot(5,1,3), semilogy(xt,abs(ch),'.-','LineWidth', 2, 'MarkerSize', 20), title('Objective change');
subplot(5,1,4), bar(xt,hst.lsit), title('Trust region size iterations');
subplot(5,1,5), bar(hst.u{end}), title('Current scaled controls');
drawnow
end
