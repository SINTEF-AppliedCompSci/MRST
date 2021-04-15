function [B, dB, R, dR, mu, P] = pvtg(tab, p, volratio, varargin)
%Interpolate live gas PVT properties.
%
% SYNOPSIS:
%   [B, dB, R, dR, mu, P] = pvtg(tab, p, volratio)
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
%   mu    - Gas viscosity.
%
%   P     - Predicate for identifying which cells are in a saturated state.
%
% SEE ALSO:
%   `pvto`, `pvdx`.

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
   sz  = size(p);
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
   R       = interpTable (tab(k,1), tab(k,2), p);
   maxR    = volratio;
   %R      = min (R, maxR);

   % If R < maxR, we have saturated case
   is      = find(R <  maxR);

   dR(is) = dinterpTable(tab(k,1), tab(k,2), p(is));
   mu(is) = interpTable (tab(k,1), tab(k,4), p(is));
   B (is) = interpTable (tab(k,1), tab(k,3), p(is));
   dB(is) = dinterpTable(tab(k,1), tab(k,3), p(is));

   % undersaturated case
   iu      = find(R >= maxR);
   if any(iu)

      R(iu)  = maxR(iu);
      dR(iu) = 0;
      %
      r      = R(iu); %maxR(iu);
      p      = p(iu);

      % group points together that refer to same table section
      delim   = [0;tab(k,1);inf];
      nbin    = numel(delim)-1;
      [n,bin] = histc  (p, delim);
      groups  = unique (bin);

      % Compute Weighted average of Bo, dBo and mu for each rl
      for i=1:length(groups)
         g      = groups(i);
         I      = bin==g;
         r1     = r(I);
         p1     = p(I);
         iuI    = iu(I);


         if g==1
            % Extrapolate from first table section
            dispif(false, 'Extrapolate from first table section\n');
            i2      = row(g):row(g+1)-1;
            mu(iuI) = interpTable (tab(i2,2), tab(i2,4), r1);
            B (iuI) = interpTable (tab(i2,2), tab(i2,3), r1);
            %dB(iuI) = dinterpTable(tab(i2,2), tab(i2,3), r1); % GALT

            i3      = row(g+1):row(g+2)-1;
            B3      = interpTable (tab(i3,2), tab(i3,3), r1);
            dB(iuI) = (B3-B(iuI))./(tab(row(g+1),1)- tab(row(g),1));

         elseif g==nbin
            % Extrapolate from last table section
            dispif(false, 'Extrapolate from last table section\n');
            i1      = row(g-1):row(g)-1;
            mu(iuI) = interpTable (tab(i1,2), tab(i1,4), r1);
            B (iuI) = interpTable (tab(i1,2), tab(i1,3), r1);
            %dB(iuI) = dinterpTable(tab(i1,2), tab(i1,3), r1); % GALT

            i0      = row(g-2):row(g-1)-1;
            B0      = interpTable (tab(i0,2), tab(i0,3), r1);
            dB(iuI) = (B(iuI)-B0)./(tab(row(g-1),1)- tab(row(g-2),1));
         else

            % Interpolate between table sections
            tab_p   = [tab(row(g-1),1), tab(row(g),1)];
            w       = [tab_p(2)-p1, p1-tab_p(1)]./diff(tab_p);

            i1      = row(g-1):row(g)-1;
            mu1     = interpTable (tab(i1,2), tab(i1,4), r1);
            B1      = interpTable (tab(i1,2), tab(i1,3), r1);
            %dB1     = dinterpTable(tab(i1,2), tab(i1,3), r1); % GALT

            i2      = row(g):row(g+1)-1;
            mu2     = interpTable (tab(i2,2), tab(i2,4), r1);
            B2      = interpTable (tab(i2,2), tab(i2,3), r1);
            %dB2     = dinterpTable(tab(i2,2), tab(i2,3), r1); % GALT

            mu(iuI) = w(:,1) .* mu1 + w(:,2) .* mu2;
            B (iuI) = w(:,1) .* B1  + w(:,2) .* B2;
            %dB(iuI) = w(:,1) .* dB1 + w(:,2) .* dB2;

            % Use simple FD approximation
            dB(iuI) = (B2-B1)./(tab_p(2)-tab_p(1));
         end
      end
   end
   P(iu) = false;
end
%--------------------------------------------------------------------------
%{
function [mu, B, dB] = interpViscVol(p, t, pi)
   if isempty(pi),
      [mu, B, dB] = deal([]);
   else
      mu =            interpTable(p,      t(:,1), pi);
      B  =   1    ./  interpTable(p, 1 ./ t(:,2), pi);
      dB = - B.^2 .* dinterpTable(p, 1 ./ t(:,2), pi);
   end
end
%}
