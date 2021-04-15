function fluid = artificialLiveOilModel(varargin)
%
% Three-phase fluid with linear relperm and constant compressibility
% live-oil PVT (from Lee, Wolfsteiner and Tchelepi).
%
   opt = struct('c_water',     convertFrom( 2.37e-5, 1/psia), ...'
                'c_gas',       convertFrom( 6.00e-4, 1/psia), ...
                'c_satoil',    convertFrom(-4.80e-5, 1/psia), ...
                'c_usatoil',   convertFrom( 1.37e-5, 1/psia), ...
                'rho_surface', [1000, 800, 100],              ...
                'mu',          [1, 2, 0.1]*centi*Poise);
   opt = merge_options(opt, varargin{:});

   fluid.info    = 'PVT model from Lee, Wolfsteinier & Tchelepi 2008.';
   fluid.pvt     = pvtLeeWolfsteinerTchelepi(opt);
   fluid.relperm = @relperm;
   fluid.names = {'aqua','liquid','vapor'};
   fluid.pc = [];
   fluid.surfaceDensity = opt.rho_surface;
end

function [kr, dkr] = relperm(s, varargin)
   assert(size(s, 2) == 3);

   kr = s.^2;
   dkr = zeros(size(s, 1), 9);
   dkr(:,[1,5,9]) = 2*s;%ones(size(dkr,1), 3);
end

function pvt = pvtLeeWolfsteinerTchelepi(opt)

   pvtfun = {...
   @(p,z) pvtwater(p, z, opt.c_water,  opt.mu(1)), ...
   @(p,z)   pvtoil(p, z, opt.c_satoil, opt.c_usatoil, opt.c_gas, opt.mu(2)),   ...
   @(p,z)   pvtgas(p, z, opt.c_gas,    opt.mu(3))};

   pvt    = @(p, z) blackOilPVT(pvtfun, opt.rho_surface, p, z);
end


function [B, dB, R, dR, mu, P] = pvtoil(p, z, c_soil, c_uoil, c_gas, mu)
%
% Live oil, i.e., oil in which gas can be dissolved.  The capacity to
% dissole gas is pressure dependant, given by the function R(p). At the
% bubble point pressure pb, all gas is dissoved in the oil. Thus, R(pb) =
% zg/zo. The bubble point pressure is a function of zg/zo.  If pressure
% falls below pb, there may is free gas and the gas saturation is non-zero.
%
% This model has constant compressibilities for water, saturated oil,
% under saturated oil and gas like in Lee, Wolfsteiner and Tchelepi
% (CompGeosci 2008).

   pr  = convertFrom(3300, psia);
   b   = stb/(ft^3);
   a   = 0.1*b/psia;
   Rs  = @(p)       a*(p-pr) + b*90;
   pb  = @(z) pr + 1/a*z(:,3)./z(:,2) - b/a*90;
   dRs = @(p) ones(size(p))*a;

   % c_soil = -1/Bo·dBo/dp + dR/dp·Bg/Bo = constant, B1=1 at 1 athmosphere.
   p0 = convertFrom(14.7, psia);
   B1  = @(p)     (1 +  a/(c_soil-c_gas)*(exp((c_soil-c_gas)*(p-p0)) - 1)) .* exp(-c_soil*(p-p0));
   dB1 = @(p)     -c_soil*B1(p) + a*exp(-c_gas*(p-p0));

   % c_uoil = -1/Bo·dBo/dp = constant, B1 = B2 at pb(z).
   B2  = @(p, z)         B1(pb(z)).*exp(-c_uoil*(p-pb(z)));
   dB2 = @(p, z) -c_uoil*B1(pb(z)).*exp(-c_uoil*(p-pb(z)));


%{
   Bo  = @(p, z) (p<pb(z)).* B1(p) + ~(p<pb(z)).* B2(p,z);
   dBo = @(p, z) (p<pb(z)).*dB1(p) + ~(p<pb(z)).*dB2(p,z);
   figure(1)
   plot(p, Bo(p, z));
   figure(2)
   drs = @(p, z) (p<pb(z)).*(dRs(p));
   Bg = @(p)exp(-c_gas *(p-p0));
   plot(p, -1./Bo(p,z).*dBo(p, z) + drs(p, z).*Bg(p)./Bo(p, z))
%}

   i = p<pb(z);
   if any(i),
      B(i,:)  = B1(p(i));
      dB(i,:) = dB1(p(i));
      R(i)    = Rs(p(i,:));
      dR(i)   = dRs(p(i,:));
   end
   if any(~i),
      B(~i,:)  = B2(p(~i), z(~i,:));
      dB(~i,:) =dB2(p(~i), z(~i,:));
      R(~i)    = Rs(pb(z(~i,:)));
      dR(~i)   = zeros(sum(~i),1) ;
   end
   mu = repmat(mu, size(p));


   P  = (p<pb(z));
   B(~isfinite(B)) = 0;
end

function [B, dB, R, dR, mu, P] = pvtwater(p, z, c_water, mu)  %#ok
   %cw =  convertFrom(2.37e-5, 1/psia);
   B  =     exp(-c_water *(p-1*atm));
   dB = -c_water*exp(-c_water *(p-1*atm));
   R  = zeros(size(p));
   dR = zeros(size(p));
   mu = repmat(mu, size(p));
   P = ones(size(p));
      %figure(3)
   %plot(p, B)

end

function [B, dB, R, dR, mu, P] = pvtgas(p, z, c_gas, mu) %#ok
   %cg =  convertFrom(60.0e-5, 1/psia);
   B  =     exp(-c_gas *(p-1*atm));
   dB = -c_gas*exp(-c_gas *(p-1*atm));
   R  = zeros(size(p));
   dR = zeros(size(p));
   mu = repmat(mu, size(p));
   P = ones(size(p));
      %figure(4)
   %plot(p, B)
end

function [c, rho, mu, u, R, B, A, dA] = blackOilPVT(pvtfun, rhos, p, z)
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
   phases= reshape(~cellfun(@isempty, pvtfun), [], 1);

   aqua   = 1;
   liquid = 2;
   vapor  = 3;

   if sum(phases) < 1,
      error('Huh!?, No phases?');
   end


   if numel(p) ~= size(z, 1),
      error('There must be one pressure for each row in mass table');
   end
   if size(z, 2) ~= sum(phases),
      error ('There must be one column in mass table for each phase');
   end


   %% Initialize return values to zero, evaluate pvt functions
   [B, dB, R, dR, mu, saturated] = deal(zeros([numel(p), numel(pvtfun)]));

   for i = 1 : numel(pvtfun),
      if isempty(pvtfun{i}), phases(i)=false; continue; end
      [B(:,i), dB(:,i), R(:,i), dR(:,i), mu(:,i), saturated(:,i)]= ...
         pvtfun{i}(p, z);
   end

   % Convert z from an (n x numphases) to an (n x 3) array (numphases <=3).
   z = [zeros(size(z, 1), 1), z];z = z(:,1+cumsum(phases).*phases);

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
   dA   = (dR*B - R*dB)/B/B;
   e    = ones(size(A,1),1);
   c    = reshape(e'*(A\dA), sum(phases), [])'  ;


   %% Reservoir phase densities
   rho  = reshape(A'*repmat(rhos', [size(p, 1), 1]), sum(phases), []) .';
   mu = mu(:,phases);

   %% Reservoir phase volumes
   u    = reshape(A\reshape(z(:,phases)', [], 1), sum(phases), []) .';

   % Fundamental requirement for PVT model for mixtures of fluids (exact
   % range of fluids uncertain... isothermal?)
   if any(sum(u.*c, 2)<0),
      warning(['Huh?! Negative total compressibility in ', mfilename,'.']); %#ok
   end

   if any(any(u>0 & Borig(:,phases)<0)),
      warning('Huh?! Negative FVF in blackOilPVT...'); %#ok
   end

   if any(any(R<0)),
  %    warning('Huh?! Negative GOR, OGR or GAR in blackOilPVT...');
   end

   % Sanity check from Coats paper.
   %if ~all(B(:, liquid) - B(:, vapor).*R(:, liquid) > 0),
   %   warning('Basic PVT sanity check failed!');
   %end
end
