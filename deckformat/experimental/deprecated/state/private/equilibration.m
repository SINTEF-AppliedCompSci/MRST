function [press, s, rs, rv] = equilibration(G, deck, fluid, pix)
%Equilibration facility to handle EQUIL keyword.
%
% Internal helper function.  Interface subject to change.
%
% SEE ALSO:
%   `initEclipseState`.

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

   assert (isfield(G.cells, 'centroids'), ...
           'Input grid must be equipped with geometric information.');

   pressure = zeros([G.cells.num, numel(fluid.names)]);
   rs       = zeros([G.cells.num, 1]);
   rv       = zeros([G.cells.num, 1]);

   [eqlnum, pvtnum] = regions(deck, G);

   act_region = @(I) find(accumarray(I(:), 1) > 0);

   for eqreg = reshape(act_region(eqlnum), 1, [])
      cells  = eqlnum == eqreg;
      pvtreg = act_region(pvtnum(cells));

      if numel(pvtreg) > 1
         error(id('EquilReg:InconsistentPVTReg')            , ...
              ['Cells in equilibration region %d reference ', ...
               '%d PVT tables.'], eqreg, numel(pvtreg));
      end

      [pressure(cells, :), rs(cells, :), rv(cells, :)] = ...
         equilibrate_region(deck, G, cells, fluid, pix, eqreg, pvtreg);
   end

   s = initsat(deck, pressure, pix, ...
               G.cells.centroids(:,3), G.cells.indexMap);

   press = pressure;

%{
   kr = fluid.relperm(s);

   % Don't pretend there's a meaningful phase pressure unless the
   % corresponding phase is mobile.
   pressure(~ (kr > 0.0)) = nan;

   liq      = strcmpi(fluid.names, 'oil');
   i        = isfinite(pressure(:, liq));
   press    = zeros([G.cells.num, 1]);
   press(i) = pressure(i, liq);
%}
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [eqlnum, pvtnum] = regions(deck, G)
   if isfield(deck, 'REGIONS') && isfield(deck.REGIONS, 'EQLNUM')
      eqlnum = reshape(deck.REGIONS.EQLNUM(G.cells.indexMap), [], 1);
   else
      eqlnum = ones([G.cells.num, 1]);
   end

   if isfield(deck, 'REGIONS') && isfield(deck.REGIONS, 'PVTNUM')
      pvtnum = reshape(deck.REGIONS.PVTNUM(G.cells.indexMap), [], 1);
   else
      pvtnum = ones([G.cells.num, 1]);
   end
end

%--------------------------------------------------------------------------

function [press, rs, rv] = ...
      equilibrate_region(deck, G, cells, fluid, pix, eqreg, pvtreg)
   equil = deck.SOLUTION.EQUIL(eqreg, :);

   Z0   = equil(1);
   P0   = equil(2);
   Zwoc = equil(3);
   Zgoc = equil(5);

   Zc   = G.cells.centroids(cells, 3);
   Zmin = min(Zc);  Zmax = max(Zc);

   Rs = define_dissolution(deck, Z0, P0, Zgoc, eqreg, pvtreg);
   Rv = define_evaporation(deck, Z0, P0, Zgoc, eqreg, pvtreg);

   if ~((Zgoc > Z0) || (Z0 > Zwoc))
      % Datum depth in oil zone  (Zgoc <= Z0 <= Zwoc)
      press = equilibrate_OWG(Zc, [Zmin, Zmax], ...
                              fluid, pix, Rs, Rv, equil(1:6));
   else
      error('Datum not in oil zone')
   end

   [rs, rv] = deal(zeros(size(Zc)));

   if all([pix.OIL, pix.GAS] > 0)
      rs = Rs(Zc, press(:, pix.OIL)); % min(Rs(z), Rs(p_oil))
      rv = Rv(Zc, press(:, pix.GAS)); % min(Rv(z), Rv(p_gas))
   end
end

%--------------------------------------------------------------------------

function s = initsat(deck, press, pix, Zc, indexMap)
   % Invert capillary pressure function to derive initial saturations.

   if isfield(deck.REGIONS, 'SATNUM')
      satnum = deck.REGIONS.SATNUM(indexMap);
   else
      satnum = ones([numel(indexMap), 1]);
   end

   if isfield(deck.REGIONS, 'EQLNUM')
      eqlnum = deck.REGIONS.EQLNUM(indexMap);
   else
      eqlnum = ones([numel(indexMap), 1]);
   end

   pairs = sortrows([eqlnum, satnum, (1 : numel(indexMap)) .']);

   [p, n] = rlencode(pairs(:, [1, 2]));
   pos    = cumsum([1 ; n]);

   [Pcow, Pcog] = deal(zeros(size(Zc)));
   [ow, go]     = deal(false);
   if all([pix.WATER, pix.OIL] > 0)
      % Oil and water active
      Pcow = press(:, pix.OIL) - press(:, pix.WATER);   % P_o - P_w
      ow   = true;
   end

   if all([pix.OIL, pix.GAS] > 0)
      % Oil and gas active
      Pcog = press(:, pix.GAS) - press(:, pix.OIL);      % P_g - P_o
      go   = true;
   end

   s = zeros(size(press));

   [inv_Pcow, inv_Pcog] = pc_inv(deck);

   for i = 1 : size(p, 1)
      c = pairs(pos(i) : pos(i + 1) - 1, end);

      [er, sr] = deal(p(i,1), p(i,2));

      if ow, s(c, pix.WATER) = inv_Pcow{er, sr}(Pcow(c), Zc(c)); end
      if go, s(c, pix.GAS)   = inv_Pcog{er, sr}(Pcog(c), Zc(c)); end
   end

   if (ow && go) && any(sum(s, 2) > 1)
      % Three-phase problem (O/W/G) for which the water transition zone
      % extends into the gas transition zone, causing negative oil
      % saturations from simplified view of initial distribution.
      %
      % Recalculate saturation distribution from gas-water capillary
      % pressure defined as Pcgw(Sw) = Pcog(Sg=1-Sw) + Pcow(Sw)) and set
      % initial oil saturation to zero in these cells.  Assume GWC=GOC
      % (gas-water contact coincides with gas-oil contact) in relevant
      % equilibration regions.
      %
      c      = sum(s, 2) > 1;
      Pcgw   = press(c, pix.GAS) - press(c, pix.WATER);
      pairs  = sortrows([eqlnum(c), satnum(c), find(c), (1 : sum(c)) .']);
      [p, n] = rlencode(pairs(:, [1, 2]));
      pos    = cumsum([1 ; n]);

      % Define derived gas-water capillary pressure curves.  This is a
      % convenience definition only that helps simplify the creation of
      % reverse interpolators through the use of function 'create_invpc'.
      %
      if all(isfield(deck.PROPS, { 'SWOF', 'SGOF' }))
         % Family I
         [wtbl, gtbl] = deal(deck.PROPS.SWOF, deck.PROPS.SGOF);
      else
         assert (all(isfield(deck.PROPS, { 'SWFN', 'SGFN' })));
         % Family II
         [wtbl, gtbl] = deal(deck.PROPS.SWFN, deck.PROPS.SGFN);
      end

      deck.PROPS.SLGOF = cell([1, p(end,2)]);  % p(end,2) == MAX(p(:,2)).
      for r = reshape(find(accumarray(p(:, 2), 1) > 0), 1, [])  % Act. reg.
         tw = wtbl{r}(:, [1, end]);
         tg = gtbl{r}(:, [1, end]);

         % Gas-oil capillary pressure as a function of Sw (Sg = 1 - Sw).
         pcgo = interp1(tg(:, 1), tg(:, 2), ...
                        1 - tw(:, 1), 'linear', 'extrap');

         % Tabulate P_{c,gw} as a function of Sw.
         deck.PROPS.SLGOF{r} = [ tw(:,1) , pcgo + tw(:,2) ];
      end

      kw = { 'SLGOF', 'foo' };  % The first table exists by construction.
      phase         = 'gas/water';
      cntct_ix      = 5;        % GWC=GOC
      is_increasing = false;    % d(P_{c,gw})/dSw <= 0

      iPcgw = cell([max(eqlnum), max(satnum)]);
      iPcgw = create_invpc(iPcgw, deck, phase, kw, ...
                           cntct_ix, is_increasing);

      for r = 1 : size(p, 1)
         ix       = pairs(pos(r) : pos(r + 1) - 1, [end - 1, end]);
         [er, sr] = deal(p(r,1), p(r,2));

         sw = iPcgw{er, sr}(Pcgw(ix(:,2)), Zc(ix(:,1)));

         s(ix(:,1), [pix.WATER, pix.GAS]) = [sw, 1 - sw];
      end
   end

   % Define oil saturation to fill pore space.
   s(:, pix.OIL) = 1 - sum(s, 2);

   assert (~any(any(s < 0)), 'Negative saturations during equilibrium');
end

%--------------------------------------------------------------------------

function Rs = define_dissolution(deck, Z0, P0, Zgoc, eqreg, pvtreg)
   if isfield(deck.PROPS, 'PVTO')
      PVTO = deck.PROPS.PVTO{pvtreg};

      interp_Rs = @(p) interp1(PVTO.data(PVTO.pos(1:end-1), 1), ...
                               PVTO.key                       , ...
                               p, 'linear', 'extrap');

      if     isfield(deck.SOLUTION, 'RSVD')

         RSVD = deck.SOLUTION.RSVD{eqreg};
         Rs   = @(z, p) interp1(RSVD(:,1), RSVD(:,2), ...
                                z, 'linear', 'extrap');

      elseif isfield(deck.SOLUTION, 'PBVD')

         PBVD = deck.SOLUTION.PBVD{eqreg};
         Rs   = @(z, p) interp_Rs(min(p, interp1(PBVD(:,1), PBVD(:,2), ...
                                                 z, 'linear', 'extrap')));

      else
         if Z0 ~= Zgoc
            error(id('DefaultRS:DepthMismatch')                  , ...
                 ['Table %d: Datum depth must coincide with GOC ', ...
                  'in absence of explicit solubility data for '  , ...
                  'equilibration.'], eqreg);
         end

         Rsmax = interp_Rs(P0);
         Rs    = @(z, p) min(interp_Rs(p), Rsmax);
      end

   else
      % Immiscible
      Rs = @(varargin) 0;
   end
end

%--------------------------------------------------------------------------

function Rv = define_evaporation(deck, Z0, P0, Zgoc, eqreg, pvtreg)
   if isfield(deck.PROPS, 'PVTG')
      PVTG = deck.PROPS.PVTG{pvtreg};

      interp_Rv = @(p) interp1(PVTG.key                       , ...
                               PVTG.data(PVTG.pos(1:end-1), 1), ...
                               p, 'linear', 'extrap');

      if     isfield(deck.SOLUTION, 'RVVD')

         RVVD = deck.SOLUTION.RVVD{eqreg};
         Rv   = @(z, p) interp1(RVVD(:,1), RVVD(:,2), ...
                                z, 'linear', 'extrap');

      elseif isfield(deck.SOLUTION, 'PDVD')

         PDVD = deck.SOLUTION.PDVD{eqreg};
         Rv   = @(z, p) interp_Rv(min(p, interp1(PDVD(:,1), PDVD(:,2), ...
                                                 z, 'linear', 'extrap')));

      else
         if Z0 ~= Zgoc
            error(id('DefaultRV:DepthMismatch')                  , ...
                 ['Table %d: Datum depth must coincide with GOC ', ...
                  'in absence of explicit solubility data for '  , ...
                  'equilibration.'], eqreg);
         end

         Rv = @(z, p) interp_Rv(min(p, P0));
      end

   else
      % Immiscible
      Rv = @(varargin) zeros([numel(varargin{1}), 1]);
   end
end

%--------------------------------------------------------------------------

function press = equilibrate_OWG(Zc, bdry, fluid, pix, Rs, Rv, equil)
   odeopts = odeset('AbsTol', 1.0e-10, 'RelTol', 5.0e-8);

   np      = sum([pix.WATER, pix.OIL, pix.GAS] > 0);
   press   = zeros([numel(Zc), np]);

   Z0   = equil(1);   P0       = equil(2);
   Zwoc = equil(3);   Pcow_woc = equil(4);  % Water-oil contact
   Zgoc = equil(5);   Pcog_goc = equil(6);  % Gas-oil contact

   assert (~(Zwoc < Zgoc), 'WOC above GOC?!?');

   Zmin = bdry(1);
   Zmax = bdry(2);

   O_up = [Z0  , Zmin];  O_down = [Z0  , Zmax];
   W_up = [Zwoc, Zmin];  W_down = [Zwoc, Zmax];
   G_up = [Zgoc, Zmin];  G_down = [Zgoc, Zmax];

   norm_g = norm(gravity());
   if ~ (norm_g > 0)
       warning('Gravity:Zero', ...
              ['Gravity is zero. Initialization results may ', ...
               'differ from expected behaviour.']);
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 1a) Oil pressure.
   %
   assert (pix.OIL > 0, 'Oil must be active phase');
   if norm_g > 0
      dp = @(z, p) dp_o(z, p, fluid.pvt, np, pix, norm_g, Rs);
   else
      dp = @(z, p) 0;
   end

   sol_u = [];
   if ~(Z0 < Zmin) && abs(diff(O_up)) > 0
      sol_u = ode45(dp, O_up, P0, odeopts);
   end

   sol_d = [];
   if ~(Z0 > Zmax) && abs(diff(O_down)) > 0
      sol_d = ode45(dp, O_down, P0, odeopts);
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 1b) Derive phase pressures at WOC and GOC
   %
   if     diff(W_up  ) > 0   % WOC above reservoir   (unexpected)
      P0w =  inf;
   elseif diff(W_down) < 0   % WOC below reservoir
      P0w = -inf;
   else                      % WOC in    reservoir
      if Zwoc < Z0
         assert (~isempty(sol_u), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set in W-zone                       (unexpected)
         P0w = deval(sol_u, Zwoc) - Pcow_woc;
      else
         assert (~isempty(sol_d), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set above W-zone
         P0w = deval(sol_d, Zwoc) - Pcow_woc;
      end
   end

   if     diff(G_up  ) > 0   % GOC above reservoir
      P0g = -inf;
   elseif diff(G_down) < 0   % GOC below reservoir   (unexpected)
      P0g =  inf;
   else                      % GOC in    reservoir
      if Zgoc < Z0
         assert (~isempty(sol_u), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set below G-zone
         P0g = deval(sol_u, Zgoc) + Pcog_goc;
      else
         assert (~isempty(sol_d), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set in G-zone                       (unexpected)
         P0g = deval(sol_d, Zgoc) + Pcog_goc;
      end
   end

   cu = Zc < Z0;
   if any(cu)
      assert (~isempty(sol_u), 'Internal error.');
      press(cu, pix.OIL) = deval(sol_u, Zc( cu));
   end

   if ~all(cu)
      cdown = Zc > Z0;

      if any(cdown)
         assert (~isempty(sol_d), 'Internal error.');
         press(cdown, pix.OIL) = deval(sol_d, Zc(cdown));
      end

      press(~(cu | cdown), pix.OIL) = P0;
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 2) Water pressure.
   %
   if pix.WATER > 0
      if norm_g > 0
         dp = @(z, p) dp_w(p, fluid.pvt, np, pix, norm_g);
      else
         dp = @(z, p) 0;
      end

      if isfinite(P0w)
         % Water present during initialisation/equilibration.

         sol_u = ode45(dp, W_up  , P0w, odeopts);
         sol_d = ode45(dp, W_down, P0w, odeopts);

         cu = Zc < Zwoc;
         if  any(cu), press( cu, pix.WATER) = deval(sol_u, Zc( cu)); end
         if ~all(cu), press(~cu, pix.WATER) = deval(sol_d, Zc(~cu)); end
      else
         press(:, pix.WATER) = P0w;
      end
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 3) Gas pressure.
   %
   if pix.GAS > 0
      if norm_g > 0
         dp = @(z, p) dp_g(z, p, fluid.pvt, np, pix, norm_g, Rv);
      else
         dp = @(z, p) 0;
      end

      if isfinite(P0g)
         % Gas present during initialisation/equilibration.

         sol_u = ode45(dp, G_up  , P0g, odeopts);
         sol_d = ode45(dp, G_down, P0g, odeopts);

         cu = Zc < Zgoc;
         if  any(cu), press( cu, pix.GAS) = deval(sol_u, Zc( cu)); end
         if ~all(cu), press(~cu, pix.GAS) = deval(sol_d, Zc(~cu)); end
      else
         press(:, pix.GAS) = P0g;
      end
   end
end

%--------------------------------------------------------------------------

function [iPcow, iPcog] = pc_inv(deck)
   if isfield(deck.REGIONS, 'EQLNUM')
      nequilreg = max(deck.REGIONS.EQLNUM);
   else
      nequilreg = 1;
   end

   if isfield(deck.REGIONS, 'SATNUM')
      nsatreg = max(deck.REGIONS.SATNUM);
   else
      nsatreg = 1;
   end

   water = isfield(deck.RUNSPEC, 'WATER') && logical(deck.RUNSPEC.WATER);
   gas   = isfield(deck.RUNSPEC, 'GAS'  ) && logical(deck.RUNSPEC.GAS  );

   [iPcow, iPcog] = deal(cell([nequilreg, nsatreg]));
   if water
      kw    = { 'SWOF', 'SWFN' };
      cntct = 3;  % WOC in EQUIL(:,3)
      iPcow = create_invpc(iPcow, deck, 'Water', kw, cntct, true);
   end

   if gas
      kw    = { 'SGOF', 'SGFN' };
      cntct = 5;  % GOC in EQUIL(:,5)
      iPcog = create_invpc(iPcog, deck, 'Gas'  , kw, cntct, false);
   end
end

%--------------------------------------------------------------------------

function ipc = create_invpc(ipc, deck, phase, kw, cntct_ix, is_increasing)
   i = isfield(deck.PROPS, kw);

   if sum(i) ~= 1
      error(id('Pc:Unsupp'), ...
           ['%s capillary pressure curves must be specified ', ...
            'through either ''%s'' or ''%s'' keywords.'], ...
            phase, kw{:});
   end

   nsatreg = size(ipc, 2);

   sfunc = deck.PROPS.(kw{i});
   assert (nsatreg <= numel(sfunc), ...
          ['Water capillary pressure must be specified in each ', ...
           'saturation region.']);

   cntct = deck.SOLUTION.EQUIL(:, cntct_ix);
   for sr = 1 : nsatreg
      t = sfunc{sr}(:, [1, end]);

      ipc(:, sr) = ...
         arrayfun(@(z) @(pc, depth) inv_interp(pc, depth, z, t, ...
                                               is_increasing), ...
                  cntct, 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['EQUIL:', s];
end

%--------------------------------------------------------------------------

function dp = dp_o(z, p, pvt, np, pix, norm_g, Rs)
   svol = zeros([1, np]);
   svol(pix.OIL) = 1;

   if pix.GAS > 0
      svol(pix.GAS) = Rs(z, p);
   end

   [rho, rho] = pvt(p, svol);  %#ok

   dp = rho(pix.OIL) * norm_g;
end

%--------------------------------------------------------------------------

function dp = dp_w(p, pvt, np, pix, norm_g)
   svol = zeros([1, np]);
   svol(pix.WATER) = 1;

   [rho, rho] = pvt(p, svol);  %#ok

   dp = rho(pix.WATER) * norm_g;
end

%--------------------------------------------------------------------------

function dp = dp_g(z, p, pvt, np, pix, norm_g, Rv)
   svol = zeros([1, np]);
   svol(pix.GAS) = 1;

   if pix.OIL > 0
      svol(pix.OIL) = Rv(z, p);
   end

   [rho, rho] = pvt(p, svol);  %#ok

   dp = rho(pix.GAS) * norm_g;
end

%--------------------------------------------------------------------------

function s = inv_interp(p, d, dc, tbl, is_increasing)
   assert (all(diff(tbl(:, 1)) > 0), ...
           'Saturation values are not monotonically increasing.');

   s = nan(size(p));

   if norm(diff(tbl(:,2))) > 0
      [min_p, min_i] = min(tbl(:,2));
      [max_p, max_i] = max(tbl(:,2));

      i = ~((p < min_p) | (p > max_p));

      [t2, ii] = unique(tbl(:,2), 'first');
      t1       = tbl(ii, 1);

      if max_p ~= min_p
         s(i)      = interp1(t2, t1, p(i));
      end
      s(p < min_p) = tbl(min_i, 1);
      s(p > max_p) = tbl(max_i, 1);

   else
      % Constant capillary pressure function.  Use step definition.

      m = min(tbl(:,1));   % Minimum (connate) saturation
      M = max(tbl(:,1));   % Maximum saturation

      i = d < dc;

      if is_increasing
         s( i) = m; s(~i) = M;
      else
         s(~i) = m; s( i) = M;
      end
   end

   assert (~any(isnan(s)), ...
           'Internal error during capillary pressure inversion.');
end
