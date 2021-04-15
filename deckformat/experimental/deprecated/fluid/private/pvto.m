function [B, dB, R, dR, mu, P] = pvto(tab, p, volratio, varargin)
%Interpolate live oil PVT properties.
%
% SYNOPSIS:
%   [B, dB, R, dR, mu, P] = pvto(tab, p, volratio)
%
% PARAMETERS:
%
% RETURNS:
%   B, dB - Volume formation factor (B) and derivative of B with respect to
%           pressure (dB).
%
%   R, dR - Miscibility ratio (R) and derivative of R with respect to
%           pressure (dR).
%
%   mu    - Oil viscosity.
%
%   P     - Predicate for identifying which cells are in a saturated state.
%
% SEE ALSO:
%   `pvtg`, `pvdx`.

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


   % Declare size
   sz = size(p);
   B  = nan(sz);
   dB = nan(sz);
   dR = nan(sz);
   mu = nan(sz);
   P  = true(sz);
   % Find row numbers such that each table section i
   % spans rows row(i):row(i+1)-1
   k       = find(~isnan(tab(:,1)));
   row     = [find(~isnan(tab(:,1))); length(tab(:,1))+1];

   % --- saturated case
   R       = interpTable (tab(k,2), tab(k,1), p);
   maxR    = volratio;
   %  R       = min (R, maxR);

   % If R < maxR, we have saturated case
   is      = find(R <  maxR);

   dR(is) = dinterpTable(tab(k,2), tab(k,1), p(is));
   [mu(is), B(is), dB(is)] = interpViscVolOil(tab(k,2), tab(k,[4,3]), p(is));

   % undersaturated case
   iu      = find(~(R < maxR));
   if any(iu)

      R(iu)  = maxR(iu);
      dR(iu) = 0;
      %
      r      = maxR(iu);
      p      = p(iu);

      % group points together that refer to same table section
      delim   = [0;tab(k,1);inf];
      nbin    = numel(delim)-1;
      [n,bin] = histc  (r, delim);
      groups  = unique (bin);

      % Compute Weighted average of Bo, dBo and mu for each rl
      for i=1:length(groups)
         g      = groups(i);
         I      = bin==g;
         p1     = p(I);
         iuI    = iu(I);


         if g <= 1
            % Extrapolate from first table section
%{
            dispif(false, 'Extrapolate from first table section\n');
            tab_r = [tab(row(1),1), tab(row(2),1)];
            w     = [tab_r(2)-r(I), r(I)-tab_r(1)] ./ diff(tab_r);


            i2 = row(1) : row(1+1) - 1;
            [mu(iuI), B(iuI), dB(iuI)] = ...
               interpViscVol(tab(i2,2), tab(i2, [4, 3]), p1);
%}
            % Extrapolate from first table section
            dispif(false, 'Extrapolate from last table section\n');
            tab_r = [tab(row(1),1), tab(row(2),1)];
            w     = [tab_r(2)-r(I), r(I)-tab_r(1)] ./ diff(tab_r);

            %i1 = row(g-2) : row(g-1) - 1;
            i1 = row(1) : row(2) - 1;
            [mu1, B1, dB1] = interpViscVolOil(tab(i1,2), tab(i1, [4, 3]), p1);

            %i2 = row(g-1) : row(g) - 1;
            i2 = row(2) : row(3) - 1;

            [mu2, B2, dB2] = interpViscVolOil(tab(i2,2), tab(i2, [4, 3]), p1);

            mu(iuI) = w(:,1) .* mu1 + w(:,2) .* mu2;
            B (iuI) = 1 ./ (w(:,1) ./ B1  + w(:,2) ./ B2);
            dB(iuI) = w(:,1) .* dB1 + w(:,2) .* dB2;
         elseif g == nbin
            % Extrapolate from last table section
            dispif(false, 'Extrapolate from last table section\n');
            tab_r = [tab(row(g-2),1), tab(row(g-1),1)];
            w     = [tab_r(2)-r(I), r(I)-tab_r(1)] ./ diff(tab_r);

            i1 = row(g-2) : row(g-1) - 1;
            [mu1, B1, dB1] = interpViscVolOil(tab(i1,2), tab(i1, [4, 3]), p1);

            i2 = row(g-1) : row(g) - 1;

            [mu2, B2, dB2] = interpViscVolOil(tab(i2,2), tab(i2, [4, 3]), p1);

            mu(iuI) = w(:,1) .* mu1 + w(:,2) .* mu2;
            B (iuI) = 1 ./ (w(:,1) ./ B1  + w(:,2) ./ B2);
            dB(iuI) = w(:,1) .* dB1 + w(:,2) .* dB2;

            %[mu(iuI), B(iuI), dB(iuI)] = ...
            %   interpViscVol(tab(i1,2), tab(i1, [4, 3]), p1);

         else

            % Interpolate between table sections
            tab_r = [tab(row(g-1),1), tab(row(g),1)];
            w     = [tab_r(2)-r(I), r(I)-tab_r(1)] ./ diff(tab_r);

            i1    = row(g-1) : row(g) - 1;
            [mu1, B1, dB1] = interpViscVolOil(tab(i1,2), tab(i1, [4, 3]), p1);

            i2      = row(g):row(g+1)-1;
            [mu2, B2, dB2] = interpViscVolOil(tab(i2,2), tab(i2, [4, 3]), p1);

            mu(iuI) = w(:,1) .* mu1 + w(:,2) .* mu2;
            B (iuI) = 1 ./ (w(:,1) ./ B1  + w(:,2) ./ B2);
            dB(iuI) = w(:,1) .* dB1 + w(:,2) .* dB2;
         end
      end
   end
   P(iu) = false;
end
%--------------------------------------------------------------------------

function [mu, B, dB] = interpViscVolOil(p, t, pi)
   if isempty(pi),

      [mu, B, dB] = deal([]);

   elseif numel(p) == 1,
      % No undersaturated data available.
      % Assume incompressible.

      mu = repmat(t(:,1), size(pi));
      B  = repmat(t(:,2), size(pi));
      dB = zeros(size(pi));

   else

      mu =            interpTable(p,      t(:,1), pi);
      B  =   1    ./  interpTable(p, 1 ./ t(:,2), pi);
      dB = - B.^2 .* dinterpTable(p, 1 ./ t(:,2), pi);

   end
end
