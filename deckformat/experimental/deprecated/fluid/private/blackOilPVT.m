function [c, rho, mu, u, R, B, A, dA, dB, dR] = ...
      blackOilPVT(pvtfun, rhos, p, z)
%Evaluate blackoil pvt functions and return phase properties.
%
% SYNOPSIS:
%   [c, rho, mu, u, R, B] = blackOilPVT(pvtfun, surfacedensity, p, z)
%
% PARAMETERS:
%
%   pvtfun- Cell array containing one pvt function for each phase.
%
%   rhos  - Surface denisty for each primary component, i.e., oil, water
%           or gas.
%
%   p     - pressure.
%
%   z     - surface volume (per porevolume) of each component.
%
% RETURNS:
%   c     - phase compressibilities.
%
%   rho   - phase densities.
%
%   mu    - phase viscosities.
%
%   u     - reservoir volume (per porevolume) for each phase.
%
%   R     - block-diagonal matrix of dissolution or vaporization ratios.
%           Each diagonal block has dimensions np x np where np is the
%           number of phases present.
%
%   u     - block-diagonal matrix of formation volume factors, i.e., ratio
%           of phase reservoir volume to surface volume of the primary
%           component. Each diagonal block is diagonal and has dimensions
%           np x np where np is the number of phases present.
%
%
% SEE ALSO:
%   `pvtg`, `pvto`.

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


% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
% To compute the effective phase properties, we use auxillary matrices R
% and B defined by the black-oil model to be
%
%         [ 1   Rl  0]                 [Bl  0   0 ]
%     R = [ Rv  1   0],     and    B = [0   Bv  0 ].
%         [ 0   0   1]                 [0   0   Ba]
%
% Then the reservoir phase volumes z = R·inv(B) u = A·u, and we define the
% phase compressibilities as
%
%               1   dc_i
%      c_i = - ---  ----,
%               c_i  dp
%
% for each phase which translates to
%                      dA
%      c  =  - inv(A)· --
%                      dp
%
% in vector notation, where dA/dp = (dR/dp·B - R·dB/dp)/( B^2 ).  The
% reservoir phase density is computed as rho = A'*rhos.

   %% Hard-wired sequence, aqua (water), liquid (oil), vapor (gas)
   phases = reshape(~cellfun(@isempty, pvtfun), [], 1);

   [aqua, liquid, vapor] = deal(1, 2, 3);

   if ~any(phases), error('Huh!?, No phases?'); end

   if numel(p) ~= size(z, 1),
      error('There must be one pressure for each row in mass table');
   end
   if size(z, 2) ~= sum(phases),
      error('There must be one column in mass table for each phase');
   end

   %% Initialize return values to zero, evaluate pvt functions
   [B, dB, R, dR, mu, saturated] = deal(zeros([numel(p), numel(pvtfun)]));

   for i = reshape(find(phases), 1, []),
      [B(:,i), dB(:,i), ...
       R(:,i), dR(:,i), mu(:,i), saturated(:,i)] = pvtfun{i}(p, z);
   end

   % Convert z from an (n x numphases) to an (n x 3) array (numphases <=3).
   z = [zeros([size(z, 1), 1]), z];
   z = z(:, 1 + cumsum(phases).*phases);

   %% Enforce undersaturated constraints
   % If the aquaic phase is missing, no gas can dissolve in the aquaic
   % phase.
   i = (z(:, aqua)>0);
   R(~i,  aqua ) = 0;
   R( i,  aqua ) = min(R( i, aqua ), z(i, vapor)./z(i, aqua));

   % If the liquid phase is undersaturated, the vapor phase is missing.
   % Consequently, no oil is allowed to evaporate.
   i = (z(:, vapor)>0);
   R(~saturated(:, liquid) | ~i,  vapor ) = 0;
   R(i, vapor) = min(R(i, vapor), z(i, liquid)./z(i, vapor));

   % Likewise, if the vapor phase is undersaturated, the liquid phase is
   % missing. Consequently, no gas can be dissolved in liquid.
   i = (z(:, liquid)>0);
   R(~saturated(:, vapor)  | ~i, liquid) = 0;
   R(i, liquid) = min(R(i, liquid), z(i, vapor)./z(i, liquid));

   %% Make sparse matrices of B, dB, R, dR
   Borig = B;
   B    = spdiags(reshape(B',[], 1), 0, numel(B), numel(B));
   dB   = spdiags(reshape(dB',[], 1), 0, numel(dB), numel(dB));

   null = zeros(size(R,1), 1);
   RgL  = reshape([null, null,        R(:,liquid)]',[],1);
   RoV  = reshape([null, R(:,vapor),  null  ]',[],1);
   RgA  = reshape([null, null,       R(:,aqua)]',[],1);

   dRgL = reshape([null, null,        dR(:,liquid)]',[],1);
   dRoV = reshape([null, dR(:,vapor), null   ]',[],1);
   dRgA = reshape([null, null,        dR(:,aqua)]',[],1);

   R    = speye(numel(R)) + ...
          spdiags([RgA,  RgL,  RoV], [2, 1, -1], numel(R), numel(R))';
   dR   = spdiags([dRgA, dRgL, dRoV], [2, 1, -1], numel(dR), numel(dR))';

   i    = repmat(phases, [numel(p), 1]);
   B    = B( i,i);
   dB   = dB(i,i);
   R    = R( i,i);
   dR   = dR(i,i);

   %% Phase compressibilities
   A    = R/B;
   dA   = (dR - A*dB) / B;
   e    = ones(size(A,1),1);
   c    = reshape(e'*(A\dA), sum(phases), [])';

   %% Reservoir phase densities
   rho  = reshape(A'*repmat(rhos', [size(p, 1), 1]), sum(phases), []) .';
   mu = mu(:,phases);

   %% Reservoir phase volumes
   u    = reshape(A\reshape(z(:,phases)', [], 1), sum(phases), []) .';

   % Fundamental requirement for PVT model for mixtures of fluids (exact
   % range of fluids uncertain... isothermal?)
   if any(sum(u.*c, 2)<0),
      warning(msgid('TotalCompressibility:Negative'), ...
              'Huh?! Negative total compressibility in ''%s''.', ...
              mfilename);
   end

   if any(any(u>0 & Borig(:,phases)<0)),
      warning(msgid('FVF:Negative'), ...
              'Huh?! Negative FVF in ''%s''...', mfilename);
   end

   if any(any(R<0)),
      warning(msgid('xxR:Negative'), ...
              'Huh?! Negative GOR, OGR or GAR in ''%s''...', mfilename);
   end

   % Sanity check from Coats paper.
   %if ~all(B(:, liquid) - B(:, vapor).*R(:, liquid) > 0),
   %   warning('Basic PVT sanity check failed!');
   %end
end
