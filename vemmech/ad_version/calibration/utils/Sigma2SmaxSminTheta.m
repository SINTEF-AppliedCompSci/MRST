function [Smax, Smin, theta] = Sigma2SmaxSminTheta(o_xx, o_yy, o_xy)
%Undocumented Utility Function

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

   tic; fprintf('Starting Sigma2SmaxSminTheta\n');
   
   % NB: The definition of Smax and Smin is chosen to be consistent with a
   % compressive regime, i.e. Smax and Smin are both negative.  As such, Smax
   % is chosen to be the smallest value in the arithmetric sense (the "most
   % negative" value)
   
   %% check if we are considering a compressive or expansive stress regime
   
   % sgn = sign(o_xx + o_xy);  % negative sign indicates compressive stress
   %                           % regime, in which case we consider Smax to be
   %                           % the maximum _negative_ stress

   sgn = -1;
   
   %% define intermediary variables
   s = o_xx + o_yy;
   d = o_xx - o_yy;
   g = 2 * o_xy;
   a = (d.^2 + g.^2).^(1/2);
   
   if isa(a, 'ADI')
      a(a.val==0) = 0; % remove broken derivatives (NaNs)
   end

   %% compute Smax and Smin
   Smax = 0.5 * (s + sgn .* a);  % @@ we assume negative (compressive) stress here
   Smin = 0.5 * (s - sgn .* a);  
   
   %% determine angles for which acos is the natural choice (otherwise, use asin)
   
   % d and g will only be used for angle calculations from here on.  We switch
   % sign according to stress regime
   d = sgn.*d;
   g = sgn.*g;
   
   acos_ix = abs(d) < abs(g);

   theta2 = 0 * Smax; % In case of ADI variables, let theta2 be one
      
   tmp = acos(d(acos_ix)./a(acos_ix)); 
   tmp_g = g(acos_ix);
                                          
   if numel(value(tmp)) > 0
      theta2(acos_ix & g >= 0) =  tmp(tmp_g >= 0);
      theta2(acos_ix & g  < 0) = -tmp(tmp_g < 0);
   end
   
   tmp = asin(g(~acos_ix)./a(~acos_ix));
   tmp_d = d(~acos_ix);
   
   % removing nans that might arise from a completely stressless system
   tmp(isnan(value(tmp))) = 0;
   
   if numel(value(tmp)) > 0
      theta2(~acos_ix & d >= 0) = tmp(tmp_d >= 0);
      theta2(~acos_ix & d <  0) = pi * ones(sum(tmp_d<0), 1) - tmp(tmp_d < 0);
   end
      
   if sum(theta2 < 0) > 0
      theta2(theta2 < 0) = theta2(theta2 < 0) + 2*pi;
   end
   
   theta = theta2./2;
   
   %% Correct for cases with isotropic stress (equal principal stresses)
   % same result, but ADI-friendly
   % margin = abs(100 * eps * s);
   % isotropic_ixs = a < margin; % derivatives of a will be corrupt when a=0
   %                                           % due to square root in its definition above.
   % % isotropic_ixs = a < abs((eps * s)); % derivatives of a will be corrupt when a=0
   % %                                % due to square root in its definition above.
   
   % % Smax(isotropic_ixs) = 0.5 * s(isotropic_ixs);  % @@ can we really ignore a here for derivatives?
   % % Smin(isotropic_ixs) = 0.5 * s(isotropic_ixs);
   
   % %theta(isotropic_ixs) = 0;  % indeterminate, so we can set it to any value
   
   % theta(isotropic_ixs) = value(theta(isotropic_ixs));
   
   % % scale theta to 0 (along with the possibly arbitrarily large derivatives
   % %damper = a(isotropic_ixs) ./ margin(isotropic_ixs);
   % %theta(isotropic_ixs) = theta(isotropic_ixs) .* damper.^2;
   
   

   %% In case of ambiguous sign of stress tensor, set angle to zero
   % this shoudl usually not happen in real situations, but if there is a
   % mixture of compressive stress and extensive stress, the definition of
   % Smax and Smin becomes ambiguous, and we set the angle to 0
   if any(sgn<0) && any(sgn>0)
      theta = 0 * theta;
   end
   fprintf('Finished Sigma2SmaxSminTheta\n');
   toc;
end
