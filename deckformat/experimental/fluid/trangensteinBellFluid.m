function fluid = initTrangensteinBellFluid(rho0, varargin)
%Construct PVT and relperm functions, constant compressibility three-phase
%mode a la Juanes&Patzek
% SYNOPSIS:
%   fluid = initCompressibleFluid(rho0, p0, c0, mu0)
%
% PARAMETERS:
%   rho0 - Fluid reference density, in units of kg/m^3.
%
%   p0   - Pressure, in units of Pascal, at which fluid reference
%          parameters (density, compressibility, and viscosity) are
%          defined.
%
%   c0   - Fluid reference compressibility at pressure 'p0'.
%
%   mu0  - Fluid reference viscosity, in units of Pa s, at pressure 'p0'.
%
% RETURNS:
%   fluid - A 1-by-1 structure array of PVT and relative permeability
%           evaluation functions.  The individual evaluation functions
%           (structure fields) are,
%             - pvt -- Evaluates pvt data for the single phase compressible
%                      fluid.  Specifically, given an n-by-1 vector of
%                      pressures p and an n-by-1 vector of surface volumes
%                      z, the call
%
%                         [c, rho, mu, u] = fluid.pvt(p, z)
%
%                      computes, respectively, n-by-1 phase
%                      compressibilities (c), densities (rho), viscosities
%                      (mu), and reservoir volumes (u).
%
%              - relperm --
%                      Evaluates one-phase relative permeability curves.
%                      Specifically, given an n-by-1 vector of fluid
%                      saturations, s (implicitly assumed s = ones([n,1])),
%                      the call
%
%                         kr = fluid.relperm(s)
%
%                      computes an n-by-1 array of relative permeability
%                      values.  As 'fluid' is a single phase fluid model,
%                      all relative permeability values are one, i.e.,
%
%                         kr = ones(size(s))
%
% EXAMPLE:
%
% SEE ALSO:
%   `initBlackoilFluid`.

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


if nargin == 0, rho0 = [1,1,1]; end

   function [c, rho, mu, u, R] = pvtfun(varargin)
      % Poor man's class replacement
      if nargin == 1,
         str = varargin{1};
         assert(ischar(str));

         if strcmpi(str, 'info'),
            assert (nargout == 0);
            disp('single compressible fluid');
            return;
         elseif strcmpi(str, 'surfacedensity'),
            assert (nargout <= 1);
            c = rho0;
            return
         elseif strcmpi(str, 'miscible'),
            assert (nargout <= 1);
            c = [true, true, true];
            return;
         elseif strcmpi(str, 'names'),
            c = {'aqua', 'liquid', 'vapour'};
            return
         else
            error('unsupported mode');
         end
      elseif nargin == 2 && isnumeric([varargin{1:2}]),
         p = varargin{1};
         z = varargin{2};
      else
         error('Unsupported input.');
      end

      if numel(p) ~= size(z, 1),
         error('There must be one pressure for each row in mass table');
      end

      % convert p from Pascal to Psia
      p = convertTo(p, psia());

      Rmax = [z(:,2)./z(:,1), z(:,1)./z(:,2), z(:,2)./z(:,3)];
      Rmax(z(:,2)==0 | z(:,1)==0, 1)=0;
      Rmax(z(:,1)==0 | z(:,2)==0, 2)  =0;
      Rmax(z(:,2)==0 | z(:,3)==0, 3)  =0;
      pb   = bubble(Rmax);


      mu          = viscosity(p, pb);
      [B, dB]     = formationVolumeFactors(p, z, Rmax);
      [R, dR]     = mix(p);
      dR(R>Rmax)= 0;
      R(R>Rmax)=Rmax(R>Rmax);
      [c, rho, u] = convertBlackoil(z, B, dB, R, dR, rho0, {3,1,2});

      % Convert from psia and cenitpoise to Pascal and Pascal·s
      c   = convertFrom(c,   1/psia());
      %rho = convertFrom(rho, 1); %??
      mu = convertFrom(mu, centi*poise());

   end

   function kr = relpermfun(s)
      bg = 0.0;

      % Juanes&Patzek "Three-Phase Displacement Theory: An Improved
      % Description of Relative Permeabilities", eqs. 32--34.
      kr = [s(:,1).^2, ...
            s(:,2).*(1-s(:,1)).*(1-s(:,3)), ...
            bg*s(:,3) + (1-bg)*s(:,3).^2];
   end



   fluid.pvt     = @pvtfun;
   fluid.relperm = @relpermfun;
end

function [R, dR] = mix(p)
   %regn om til Pascal.
   R = [0.05*p, 9.0e-5-6.0e-8*p+1.6e-11*p.^2, 0.005*p];
   dR= [repmat(0.05, size(p)), -6.0e-8+2*1.6e-11*p, repmat(0.005, size(p))];
end

function pb = bubble(Rmax)
   pb = zeros(size(Rmax));
   pb(:,1) = Rmax(:,1)./0.05; % bubble point for gas in liquid
   pb(:,2) = Rmax(:,2);        % vapour point oil in vapour is not uniquie.
   pb(:,3) = Rmax(:,3)./0.005; % bubble point for gas in aqua
end

function mu = viscosity(p, pb)
   mu = zeros(size(pb));

   s = p<pb(:,1);
   mu( s,1) = 0.8-1e-4*p(s,1);
   mu(~s,1) = (0.8-1e-4*pb(~s,1)).*(1+6.78e-5.*(p(~s,1)-pb(~s,1)));

   mu(:,2) = 0.012 + 3e-5*p;

   s = p<pb(:,3);
   mu(s,3) = 0.35;
   mu(~s,3) = 0.35*(1+678e-5.*(p(~s,1)-pb(~s,3)));
end

function [B, dB] = formationVolumeFactors(p, z, Rmax)

   B = zeros(size(z));
   %Rmax = [z(:,3)./z(:,1), z(:,3)./z(:,2), z(:,2)./z(:,3)];
   pb   = bubble(Rmax);
   R = mix(p);

   s = R(:,1)<Rmax(:,1);
   r = abs(R(:,1)) < sqrt(eps);


   B(s &  r, 1) = 1 - 2.31e-5*p(s & r);
   B(s & ~r, 1) = 1 + 1.50e-4*p(s &~r);
   B(~s,     1) =(1 + 1.50e-4*pb(~s))./(1+2.31e-5*(p(~s)-pb(~s, 1)));
   dB(s &  r, 1) =  -2.31e-5+0*p(s & r);
   dB(s & ~r, 1) =   1.50e-4+0*p(s &~r);
   dB(~s,     1) =  -B(~s, 1).*2.31e-5./(1+2.31e-5*(p(~s)-pb(~s, 1)));

   s = R(:,2)<Rmax(:,2);
   B( s, 2) = 1./(6+0.06*p(s));
   B(~s, 2) = 1./(7+0.06*p(~s,1))+Rmax(~s,2)./R(~s,2).*(1./(6+0.06*p(~s,1)) -1./(7+0.06*p(~s,1)));
   dB( s, 2) = -B(s, 1).*0.06./(6+0.06*p(s,1));
   dB(~s, 2) = -0.06./(7+0.06*p(~s,:)).^2 +Rmax(~s,2)./R(~s,2).*(-0.06./(6+0.06*p(~s,1)).^2 -0.06./(7+0.06*p(~s,1)).^2);

   s = R(:,2)<Rmax(:,2);
   r = abs(R(:,3))<sqrt(eps);
   B(s &  r, 3) = 1-1.8e-5*p(s & r);
   B(s & ~r, 3) = 1-3e-6*p(s & ~r);
   B(~s, 3)     = (1-3e-6*pb(~s, 3))./(1+1.8e-5*(p(~s,1)-pb(~s, 3)));
   dB(s &  r, 3) = -1.8e-5+0*p(s & r);
   dB(s & ~r, 3) = -3e-6+0*p(s & ~r);
   dB(~s, 3)     = -B(~s, 3).*(1.8e-5+0*p(~s,1))./(1+1.8e-5*(p(~s,1)-pb(~s, 3)));

end


