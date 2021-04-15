function state = initEclipseState(G, deck, fluid)
%Initialise reservoir state using hydrostatic equilibration
%
% SYNOPSIS:
%   state = initEclipseState(G, deck, fluid)
%
% PARAMETERS:
%   G     - A grid structure discretising a specific reservoir geometry.
%
%   deck  - An ECLIPSE/FrontSim input deck defining a particular run.  It
%           is assumed that the grid 'G' is created from the information
%           contained in 'deck'.
%
%   fluid - A fluid structure as defined by function 'initEclipseFluid'.
%
% RETURNS:
%   state - An initialised reservoir state object in hydrostatic
%           equilibrium.  The pressure values correspond to the initial oil
%           pressure distributions while saturations/masses are defined to
%           balance capillary pressure forces.
%
% RESTRICTION:
%   This function is presently only applicable to three-phase black-oil
%   problems for which the equilibrium datum depth is in the oil zone.
%
% SEE ALSO:
%   `readEclipseDeck`, `convertDeckUnits`, `initEclipseGrid`, `initEclipseFluid`.

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


   assert (isstruct(deck) && ...
           all(isfield(deck, {'RUNSPEC', 'SOLUTION'})),   ...
          ['Parameter ''deck'' does not appear to be a ', ...
           'valid ECLIPSE deck']);

   assert (isfield(deck.RUNSPEC, 'OIL') && logical(deck.RUNSPEC.OIL), ...
          ['Function ''%s'' is only supported in runs that declare ', ...
           'OIL as an active phase.'], mfilename);

   pix = phase_index(fluid, {'WATER', 'OIL', 'GAS'});

   if isfield(deck.SOLUTION, 'EQUIL'),

      [ppress, s, rs, rv] = equilibration(G, deck, fluid, pix);

   elseif any(isfield(deck.SOLUTION, {'PRESSURE', 'PRVD'})) && ...
          any(isfield(deck.SOLUTION, {'SGAS', 'SOIL', 'SWAT'})),

      [press, s, rs, rv] = direct_assignment(G, deck, fluid, pix);

      if isfield(fluid, 'pc'),
         ppress = bsxfun(@plus, press, fluid.pc(s));
      else
         ppress = repmat(press, [1, size(s, 2)]);
      end

   else

      error(msgid('Scheme:Unknown'), ...
           ['Initialisation scheme specified in ', ...
            '''deck'' is not supported.']);

   end

   assert (~ any(any((s < 0) | (s > 1))), ...
           'Meaningless initial saturations outside [0,1]');

   pressure = ppress(:, pix.OIL);
   z        = derive_masses(pix, ppress, s, rs, rv, fluid);

   cellno   = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
   fp       = sparse(double(G.cells.faces(:,1)), 1 : numel(cellno), 1) ...
              * [ pressure(cellno), ones([numel(cellno), 1]) ];

   % Check initial volume discrepancy.
   [u, u, u, u] = fluid.pvt(pressure, z);                              %#ok
   e = sqrt(eps(1));
   if any(abs(sum(u, 2) - 1) > e),
      if mrstVerbose, s0 = s; end

      s = bsxfun(@rdivide, u, sum(u, 2));

      if mrstVerbose,
         ds    = s - s0;
         ds    = [ min(ds) ; max(ds) ];
         i     = 1 + double(abs(ds(1, :)) < abs(ds(2, :)));
         maxds = ds(sub2ind(size(ds), i, 1 : size(ds, 2)));

         maxds_str = [ '[', sprintf(' %.2g', maxds), ' ]' ];

         warning(msgid('VolDiscr:Large'), ...
                ['Initial volume discrepancy too large.\n -> Using '  , ...
                 'saturations derived from initial mass distribution.', ...
                 '\n -> Largest saturation change: ds = %s.'], maxds_str);
      end
   end

   assert (~ any(any(abs(bsxfun(@rdivide, u, sum(u,2)) - s) > e)), ...
          ['Initial mass distribution does not match ', ...
           'initial saturation distribution.']);

   state = struct('pressure'    , ppress(:, pix.OIL)     , ...
                  'facePressure', fp(:,1) ./ fp(:,2)     , ...
                  'flux'        , zeros([G.faces.num, 1]), ...
                  's'           , s                      , ...
                  'z'           , z                      , ...
                  'rs'          , rs                     , ...
                  'rv'          , rv                     );
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function pix = phase_index(fluid, pname)
   a      = [reshape(lower(fluid.names), 1, []); ...
             num2cell(1 : numel(fluid.names))];
   s      = struct(a{:});

   pix    = num2cell(zeros([1, numel(pname)]));
   i      = isfield(s, lower(pname));
   pix(i) = cellfun(@(f) s.(lower(f)), pname(i), 'UniformOutput', false);

   pix    = cell2struct(pix, reshape(pname, 1, []), 2);
end

%--------------------------------------------------------------------------

function z = derive_masses(pix, p, s, rs, rv, fluid)
% Compute z = R*inv(B)*s

   np = numel(fluid.names);

   assert (np == sum(structfun(@(x) x > 0, pix)));
   assert (all(structfun(@(x) x, pix) <= np));

   if all([pix.OIL, pix.GAS] > 0),

      rs(~ (s(:, pix.OIL) > 0)) = 0;
      rv(~ (s(:, pix.GAS) > 0)) = 0;

   else

      assert (all(all([rs, rv] == 0)));

   end

   %-----------------------------------------------------------------------
   % Compute B(:,p) only where saturation of phase p is nonzero (and if the
   % phase is declared active in the run).

   % Preallocation.
   B   = ones([size(p, 1), np]);
   fvf = @(p, z, o) compute_fvf(fluid, p, z, [o, np]);

   %-----------------------------------------------------------------------

   if pix.WATER > 0,
      % FVF for the water phase.
      kw = pix.WATER;
      pw = p(:, kw);

      i  = isfinite(pw);
      if any(i),
         z       = zeros([1, np]);
         z(kw)   = 1;
         z       = repmat(z, [sum(i), 1]);

         B(i,kw) = fvf(pw(i), z, kw);                          clear pw z i
      end
   end

   %-----------------------------------------------------------------------

   if pix.OIL > 0,
      % FVF for the oil phase.
      ko = pix.OIL;
      po = p(:, ko);

      i = find(isfinite(po));
      if ~isempty(i),
         z       = zeros([1, np]);
         z(:,ko) = 1;
         z       = repmat(z, [numel(i), 1]);

         if pix.GAS > 0,
            % OIL-GAS or WATER-OIL-GAS system.
            %
            % Account for miscibility/dissolution of gas into oil.
            kg = pix.GAS;

            % Saturated oil.
            j1        = s(i, kg) > 0;
            z(j1, kg) = inf;                                   clear     j1

            % Undersaturated oil.
            j2        = ~ (s(i, kg) > 0);   % == ~j1
            z(j2, kg) = rs(i(j2));                             clear     j2
         end

         B(i,ko) = fvf(po(i), z, ko);                          clear po z i
      end
   end

   %-----------------------------------------------------------------------

   if pix.GAS > 0,
      % FVF for the gas phase.
      kg = pix.GAS;
      pg = p(:, kg);

      i = find(isfinite(pg));
      if ~isempty(i),
         z       = zeros([1, np]);
         z(:,kg) = 1;
         z       = repmat(z, [numel(i), 1]);

         if pix.OIL > 0,
            % OIL-GAS or WATER-OIL-GAS system.
            %
            % Account for miscibility/evaporation of oil into gas.
            ko = pix.OIL;

            % Saturated gas.
            j1        = s(i, ko) > 0;
            z(j1, ko) = inf;                                   clear     j1

            % Undersaturated gas.
            j2        = ~ (s(i, ko) > 0);   % == ~j1
            z(j2, ko) = rv(i(j2));                             clear     j2
         end

         B(i,kg) = fvf(pg(i), z, kg);                          clear pg z i
      end
   end

   %-----------------------------------------------------------------------
   % Final result:
   %
   %    z = R*inv(B)*s
   %

   assert (all(all(isfinite(B) & (B > 0))), ...
           'Non-finite or non-positive initial FVF in run.');

   z = s ./ B;

   if all([pix.OIL, pix.GAS]),
      k      = [pix.OIL, pix.GAS];
      z(:,k) = z(:,k) + z(:, fliplr(k)).*[rv, rs];
   end
end

%--------------------------------------------------------------------------

function B = compute_fvf(fluid, p, z, os)
   % [c, rho, mu, u, R, B]
   [  B,  B , B , B, B, B] = fluid.pvt(p, z);                          %#ok

   offset = os(1);  assert (offset > 0);
   stride = os(2);  assert (stride > 0);

   B = spdiags(B, 0);
   B = B(offset : stride : end);
end
