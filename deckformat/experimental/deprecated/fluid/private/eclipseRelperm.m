function [relperm, pcap, info] = eclipseRelperm(deck, varargin)
%Construct relperm and capillary pressure evaluators from tabulated data.
%
% SYNOPSIS:
%   [relperm, pcap] = eclipseRelperm(deck)
%   [relperm, pcap] = eclipseRelperm(deck, 'key', value, ...)
%
% PARAMETERS:
%   deck - An ECLIPSE/FrontSim input deck as defined by functions
%          'readEclipseDeck' and 'convertDeckUnits'.
%
% OPTIONAL PARAMETERS:
%   verbose - Display (short) summary of fluid relperm model.
%
% RETURNS:
%   relperm - Function handle supporting the calling sequence
%
%                  kr       = relperm(s)
%                 [kr, dkr] = relperm(s)
%
%             The first call computes the relative permeability curves at
%             the current saturations, 's'.  The second call additionally
%             computes the first partial derivatives with respect to phase
%             saturations.
%
%   pcap    - Function handle supporting the calling sequence
%
%                  pc       = pcap(s)
%                 [pc, dpc] = pcap(s)
%
%             The first call computes the capillary pressure curves at
%             the current saturations, 's'.  The second call additionally
%             computes the first partial derivatives with respect to phase
%             saturations.
%
% SEE ALSO:
%   `readEclipseDeck`, `convertDeckUnits`.

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


   opt = struct('verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   % Hard-wired sequence, aqua (water), liquid (oil), vapor (gas)
   A = 1;   L = 2;    V = 3;

   % Declaration of which phases are present.
   %                               A=1      L=2    V=3
   phase = isfield(deck.RUNSPEC, {'WATER', 'OIL', 'GAS'});
   pos   = cumsum(phase);

   % User readable phase names
   names = {'water', 'oil', 'gas'};
   names = names(phase);

   family_I  = any(isfield(deck.PROPS, {'SWOF', 'SGOF', 'SLGOF'}));
   family_II = any(isfield(deck.PROPS, {'SWFN', 'SGFN', 'SOF2', 'SOF3'}));

   assert (~ (family_I && family_II), ...
          ['Cannot mix keywords from families I and II. ', ...
           'See ECLIPSE documentation.']);

   assert (family_I || family_II, ...
           'No (known) relperm keywords specified.');

   if family_I,
      [kr, pc, connate, smax, scrit] = ...
         process_family_one(deck, phase, A, L, V);
   else
      [kr, pc, connate, smax, scrit] = ...
         process_family_two(deck, phase, A, L, V);
   end

   if all(phase),
      % Three-phase system.

      [relperm, pcap, relpermname] = ...
         build_3p_satfunc(deck, kr, pc, connate, smax, scrit);

   elseif sum(phase) == 2,
      % Two-phase system.

      [relperm, pcap] = build_2p_satfunc(kr, pc, phase, pos, A, L, V);

   elseif sum(phase) == 1,
      % Single-phase system.

      relperm = @singlephase;
      pcap    = @singlephase;

   else

      error(msgid('PhaseComb:NotImplemented'), ...
            'MRST only supports 1-3 phases.');

   end


   % Display summary
   info = '';
   switch sum(phase),
      case 1, info = [info, sprintf('Single phase %s model.', names{1})];
      case 2, info = [info, sprintf('Two-phase %s-%s model.', names{:})];
      case 3,
         info = [info, ...
                 'Three-phase model using ', relpermname, ' oil relperm.'];
   end

   if family_I,
      info = [info, sprintf('\n'), ...
            sprintf('Relperm curves specified with family I keywords.')];
   else
      info = [info, sprintf('\n'), ...
            sprintf('Relperm curves specified with family II keywords.')];
   end

   info = [info, sprintf('\n')];
   if opt.verbose,
      fprintf('%s\n', info);
   end
end

%--------------------------------------------------------------------------

function [kr, pc, connate, smax, scrit] = ...
      process_family_one(deck, phase, A, L, V)
   warn_multi = 0;

   Swco = 0;
   if isfield(deck.PROPS, 'SWOF'),
      if numel(deck.PROPS.SWOF) > 1,
         warn_multi = warn_multi + 1;
      end

      assert (phase(A) && phase(L), ...
             ['The SWOF-type saturation function is not supported ', ...
              'unless both water and oil are declared as active phases.']);

      [krw, krow, krocw, Swco, Swcr, Sowcr, Swmax, pcow] = ...
         swof(deck.PROPS.SWOF{1});

      kr.w  = krw {1};
      kr.ow = krow{1};
      pc.ow = pcow{1};

      connate.water = Swco(1);
      connate.krocw = krocw(1);

      smax.w   = Swmax(1);
      scrit.w  = Swcr(1);
      scrit.ow = Sowcr(1);
   end

   if isfield(deck.PROPS, 'SLGOF'), % convert to SGOF

      assert (~isfield(deck.PROPS, 'SGOF'), ...
             ['Gas relative permeability must be specified using ', ...
              'either ''SGOF'' or ''SLGOF'', not both.']);

      convert = @(t) flipud([max(1 - t(:,1), 0), t(:, 2:end)]);

      deck.PROPS.SGOF = cellfun(convert, deck.PROPS.SLGOF, ...
                                'UniformOutput', false);

      deck.PROPS = rmfield(deck.PROPS, 'SLGOF');

      clear convert
   end

   if isfield(deck.PROPS, 'SGOF'),
      if numel(deck.PROPS.SGOF) > 1,
         warn_multi = warn_multi + 1;
      end

      assert (phase(L) && phase(V), ...
              'Oil and gas must be declared active phases to use SGOF');

      [krg, krog, krocg, Sgco, Sgcr, Sogcr, Sgmax, pcog] ...
         = sgof(deck.PROPS.SGOF{1}, Swco);

      kr.g  = krg {1};
      kr.og = krog{1};
      pc.og = pcog{1};

      connate.gas   = Sgco(1);
      connate.krocg = krocg(1);

      smax.g   = Sgmax(1);
      scrit.g  = Sgcr(1);
      scrit.og = Sogcr(1);
   end

   if warn_multi > 0,
      warning('eclipseRelPerm:MultiRegion:NotSupported', ...
             ['Multiple saturation function regions are not ', ...
              'supported in Family I.\nUsing first region only.']);
   end
end

%--------------------------------------------------------------------------

function [kr, pc, connate, smax, scrit] = ...
      process_family_two(deck, phase, A, L, V)

   warn_multi = 0;

   connate.water = 0;

   if isfield(deck.PROPS, 'SWFN'),
      if numel(deck.PROPS.SWFN) > 1,
         warn_multi = warn_multi + 1;
      end

      assert (phase(A), ...
             ['The SWFN-type saturation function is not supported ', ...
              'unless water is declared as an active phase.']);

      [krw, Swco, pcow, Swcr, Swmax] = swfn(deck.PROPS.SWFN{1});

      kr.w  = krw {1};
      pc.ow = pcow{1};

      connate.water = Swco(1);
      smax.w        = Swmax(1);
      scrit.w       = Swcr(1);
   end

   if isfield(deck.PROPS, 'SGFN'),
      if numel(deck.PROPS.SGFN) > 1,
         warn_multi = warn_multi + 1;
      end

      assert (phase(V), ...
             ['The SGFN-type saturation function is not supported ', ...
              'unless gas is declared as an active phase.']);

      [krg, Sgco, Sgcr, Sgmax, pcog] = ...
         sgfn(deck.PROPS.SGFN{1}, connate.water);

      kr.g  = krg {1};
      pc.og = pcog{1};

      connate.gas = Sgco(1);
      smax.g      = Sgmax(1);
      scrit.g     = Sgcr(1);
   end

   if isfield(deck.PROPS, 'SOF2'),
      if numel(deck.PROPS.SOF2) > 1,
         warn_multi = warn_multi + 1;
      end

      assert (phase(L), ...
             ['The SOF2-type saturation function is not supported ', ...
              'unless oil is declared as an active phase.']);

      assert (xor(phase(A), phase(V)), ...
             ['The SOF2-type saturation function is not supported ', ...
              'unless *either* gas *or* water is declared as an '  , ...
              'active phase.']);

      [kro, Soco, Socr, Somax, kroSomax] = sof2(deck.PROPS.SOF2{1});

      if phase(V),
         kr.og = kro{1};
      else
         kr.ow = kro{1};
      end

      connate.oil = Soco(1);
      scrit.o     = Socr(1);
      smax.o      = Somax(1);

      if smax.o > 1 - connate.water,

         error('SoMax:NoSwco', ...
              ['Maximum oil saturation does not account for presence ', ...
               'of connate water\nMaximum oil (=%.2f) should equal ',   ...
               '1-Swco (=%.2f)'], smax.o, 1 - connate.water);

      elseif smax.o < 1 - connate.water,

         error('SoMax:NoSwco', ...
              ['Maximum oil saturation does not match connate water\n', ...
               'Maximum oil (=%.2f) should equal 1-Swco (=%.2f)'],      ...
               smax.o, 1 - connate.water);
      end

      connate.krocw = kroSomax(1);
   end

   if isfield(deck.PROPS, 'SOF3'),
      if numel(deck.PROPS.SOF3) > 1,
         warn_multi = warn_multi + 1;
      end

      assert (phase(A) && phase(L) && phase(V), ...  % == all(phase)
             ['The SOF3-type saturation function is not supported ', ...
              'unless water, oil and gas are all declared active phases.']);

      assert (deck.PROPS.SOF3{1}(end,2) == deck.PROPS.SOF3{1}(end, 3) , ...
             ['Oil relative permeability at maximum oil saturation '  , ...
              'should be equal in krow (oil-water) and krog (oil-gas)', ...
              ' coulmns (SOF3 table).']);

      [krow, krog, Soco, Socr, Sorw, Sorg, Somax, kroSomax] = ...
         sof3(deck.PROPS.SOF3{1});

      assert (Somax(1) == 1 - connate.water, ...
             ['The maximum oil saturation (=%f) should equal '   , ...
             '1-(connate water saturation) (=%f) in SOF3 table.'], ...
             Somax(1), 1 - connate.water);

      kr.ow = krow{1};
      kr.og = krog{1};

      connate.oil = Soco(1);
      scrit.o     = Socr(1);
      scrit.ow    = Sorw(1);
      scrit.og    = Sorg(1);

      connate.krocw = kroSomax(1);
      connate.krocg = kroSomax(1);
   end

   if warn_multi > 0,
      warning('eclipseRelPerm:MultiRegion:NotSupported', ...
             ['Multiple saturation function regions are not ', ...
              'supported in Family II.\nUsing first region only.']);
   end
end

%--------------------------------------------------------------------------

function [relperm, pcap, relpermname] = ...
      build_3p_satfunc(deck, kr, pc, connate, smax, scrit)

   krocw = connate.krocw;
   krocg = connate.krocg;
   Swco  = connate.water;
   Sgmax = smax.g;
   Sowcr = scrit.ow;
   Sogcr = scrit.og;

   krw  = kr.w;
   krow = kr.ow;
   krg  = kr.g;
   krog = kr.og;

   if all(isfield(deck.PROPS, {'SWOF', 'SGOF'})) && (krocw ~= krocg),
      warning(msgid('Kro:Inconsistent'), ...
             ['Oil relperm at maximum oil saturation in ', ...
              'SWOF and SGOF tables differ.']);
   end

   if Sgmax ~= 1 - Swco,
      warning (msgid('Sgmax:Inconsistent'), ...
              ['The maximum gas saturation (=%.2f) should equal ', ...
               '1-(connate water saturation) (=%.2f).'], Sgmax, 1 - Swco);
   end

   % --- SIMPLE (MRST extension) ---
   if isfield(deck.PROPS, 'SIMPLE'),
      relpermname = 'Simple';
      relperm     = @(s, varargin) Simple(s, Swco, krw, krow, krg);

   % --- STONE I -------------------
   elseif isfield(deck.PROPS, 'STONE1'),
      relpermname = 'Stone I';
      relperm     = @(s, varargin) ...
         StoneI(s, Swco, krw, krow, krg, krog, krocw, Sowcr, Sogcr);

   % --- STONE II ------------------
   elseif any(isfield(deck.PROPS, {'STONE', 'STONE2'})),
      relpermname = 'Stone II';
      relperm     = @(s, varargin) ...
         StoneII(s, Swco, krw, krow, krg, krog, krocw);

   % --- VERTICAL AVERAGE IN CELL --
   else
      relperm = @(s, varargin) threephase(s, Swco, krw, krow, krg, krog);
      relpermname = 'ECLIPSE default';
   end

   pcap = @(s) threephasecap(s, pc.ow, @null, pc.og);
end

%--------------------------------------------------------------------------

function [relperm, pcap] = build_2p_satfunc(kr, pc, phase, pos, A, L, V)
   if phase(A) && phase(L),
      % -------------------------------------------------------------------
      % Water/Oil system
      % -------------------------------------------------------------------

      relperm = @(s, varargin) ...
         twophase(s, @(s) kr.w (s(:, pos(A))), ...
                     @(s) kr.ow(s(:, pos(L))));

      pcap = @(s) twophasecap(s(:, pos([A, L])), pc.ow, @null);

   elseif phase(L) && phase(V),
      % -------------------------------------------------------------------
      % Gas/Oil system
      % -------------------------------------------------------------------

      relperm = @(s, varargin) ...
         twophase(s, @(s) kr.og(s(:, pos(L))), ...
                     @(s) kr.g (s(:, pos(V))));

      pcap = @(s) twophasecap(s(:, pos([L, V])), @null, pc.og);

   elseif phase(A) && phase(V),
      % -------------------------------------------------------------------
      % Water/Gas system
      % -------------------------------------------------------------------

      error('Water/gas systems are not supported at this time.');
   end
end

%--------------------------------------------------------------------------

function varargout = singlephase(s, varargin)
   varargout{1} = s;
   if nargout > 1,
      varargout{2} = ones(size(s));
   end
   if nargout > 2,
      varargout{3} = zeros(size(s));
   end
end

%--------------------------------------------------------------------------

function varargout = threephase(s, Swco, krw, krow, krg, krog)
% This function implements the default Eclipse 3-phase model. The relative
% permeability for oil is approximated from input relative permeabilities
% for two-phase oil in water and oil in gas systems (plus connate water).
% The model assumes complete segregation of water and gas in each cell,
% with water an gas occupying separate fractions of space xw and xg,
% referred to as the water zone and the gas zone respectively,  given by
%
%          Sw - Swco                    Sg
%    xw = -------------        xg = --------------.
%         Sw -Swco + Sg             Sw -Swco + Sg
%
%
% Assuming oil is distributed throughtout the cell, the oil relative
% permeability is approximated by
%
%    ko(Sw, Sg, So) = xw(Sw, Sg)·krow(So) + xg(Sw, Sg)·krog(So).
%
% For volumes of water and gas to be preserved, the water saturation in the
% water zone is Sw* = Sw+Sg, likewise Sg* = Sw+Sg.  Thus, the water and oil
% relative permeabilities are given by
%
%    kw(Sw, Sg) = xw(Sw, Sg)·krw(Sw*)
%
% and
%
%    kg(Sw, Sg) = wg(Sw, Sg)·krg(Sg*)
%
% The derivatives of phase relative permeabilities with respect to each
% saturation can then be computed,
%
%    dkw   dxw          dkrw    xg          dkrw
%    --- = ---·krw + xw·---- =  --·krw + xw·----,  with d = Sw -Swco + Sg,
%    dSw   dSw          dSw*    d           dSw*
%
%    dkw
%    --- = 0,
%    dso
%
%    dkw   dxw          dkrw      xw          dkrw
%    --- = ---·krw + xw·---- = -  --·krw + xw·----,
%    dSg   dSg          dSg*      d           dSw*
%
%
%    dko   xg
%    --- = --·(krow - krog),
%    dSw   d
%
%    dkw      dkrow        dkrog
%    --- = xw·-----  +  xg·-----,
%    dSo       dSo          dSo
%
%    dko   xw
%    --- = --·(krog - krow),
%    dSg   d
%
%
%    dkg     xg          dkrg
%    --- = - --·krg + xg·----,
%    dSw     d           dSg*
%
%    dkg
%    --- = 0,
%    dSo
%
%    dkg   xw           dkrg
%    --- = --·krg + xg· ----,
%    dSg   d            dSg*
%
%
%
%
%
% REMARK:
%    For cells with only oil (d==0), we do not differentiate the water and
%    gas zone fractions since these are nonsensical in this case.
%
%
%
%
% The Jacobian matrix is arranged as
%
%
%           /                            \
%          |   dkrw      dkrw      dkrw   |
%          |   ----      ----      ----   |
%          |   dsw       dso       dsg    |
%          |                              |
%          |                              |
%          |   dkro      dkro      dkro   |
%    J =   |   ----      ----      ----   |   (*).
%          |   dsw       dso       dsg    |
%          |                              |
%          |                              |
%          |   dkrg      dkrg      dkrg   |
%          |   ----      ----      ----   |
%          |   dsw       dso       dsg    |
%           \                            /
%
%
% BEWARE!!
% In the return value varargout{2}, the entries in (*) are arranged
% column-wise in, i.e., [J11, J21, J31, J12, J22, J32, J13, J23, J33].

   sw   = max(Swco, s(:,1));
   so   = s(:, 2);
   sg   = s(:, 3);
   null = zeros(size(s, 1), 1);

   x = [(sw - Swco)./max(sg + sw - Swco, eps), null, null];
   x(:,end) = 1-x(:,1);
   s = [sg + sw, so, sw - Swco + sg];

   nout = nargout;

   [kw{1:nout}]  = krw (s(:,1));
   [kg{1:nout}]  = krg (s(:,3));
   [kwo{1:nout}] = krow(s(:,2));
   [kgo{1:nout}] = krog(s(:,2));

   varargout{1} = [x(:,1).*kw{1}, ...
                   x(:,1).*kwo{1} + x(:,3).*kgo{1}, ...
                   x(:,3).*kg{1}];
   assert(all(varargout{1}(:)>=-1e-5))
   if nargout == 2,

      d = sg + sw - Swco;
      i = abs(d)>0;

      varargout{2} = zeros(size(x, 1), 9);

      % If d <> 0, differentiate water zone and gas zone fractions, ...
      if any(i),
         varargout{2}(i,:)  = [                                   ...
            x(i,3).* kw{1}(i)./d(i)  + x(i,1).*kw{2}(i),         ...
            x(i,3).*(kwo{1}(i)       -         kgo{1}(i))./d(i), ...
           -x(i,3).* kg{1}(i)./d(i)  + x(i,3).*kg{2}(i),         ...
            ...
            null(i),                                             ...
            x(i,1).*kwo{2}(i)  + x(i,3).*kgo{2}(i),              ...
            null(i),                                             ...
            ...
           -x(i,1).* kw{1}(i)./d(i)  + x(i,1).*kw{2}(i),         ...
            x(i,1).*(kgo{1}(i)       -         kwo{1}(i))./d(i), ...
            x(i,1).* kg{1}(i)./d(i)  + x(i,3).*kg{2}(i)          ...
            ];
      end

      % otherwise, don't.
      i = ~i;
      if any(i),
      varargout{2}(i,:)  = [                                   ...
          kw{2}(i),...   x(i,1).*kw{2}(i),                                    ...
          null(i),                                             ...
          null(i), ...    x(i,3).*kg{2}(i),                                    ...
                                                               ...
          null(i),                                             ...
          x(i,1).*kwo{2}(i)  + x(i,3).*kgo{2}(i),              ...
          null(i),                                             ...
                                                               ...
          null(i), ...     x(i,1).*kw{2}(i),                                    ...
          null(i),                                             ...
          kg{2}(i), ...    x(i,3).*kg{2}(i)                                     ...
         ];
      end
   end
end

%--------------------------------------------------------------------------

function varargout = ...
      StoneI(s, Swco, krw, krow, krg, krog, krocw, Sowcr, Sogcr)
% This function implements a version of the first model suggested by Stone
% that is close to the ECLIPSE version. The oil relative permeability is
% approximated from from input relative permeabilities for two-phase oil in
% water and oil in gas systems (plus connate water).
%
% The oil relatve permeability kro is given by
%
%    kro = krocw·SSo·Fw·Fg,
%
% where krocw is the relative permeability of oil in the presence of
% connate water only. At this saturation, the krow and krog relative
% permeabilities must coincide.  Furthermore,
%
%            krow                         krog
%    Fw = ------------ ,         Fg = --------------
%         krocw·(1-SSw)                krocw·(1-SSg)
%
% and the scaled saturations SSw, SSo and SSg are given by
%
%            Sw - Swco
%    SSw = --------------- ,
%           1 - Swco - Som
%
%            So - Som
%    SSo = --------------- ,
%           1 - Swco - Som
%
%                Sg
%    SSg = --------------- .
%           1 - Swco - Som
%
% The minimal residual oil saturation Som is by default the minimum of the
% residual oil saturations in the oil-water and the oil-gas cases, i.e.,
% Som = min(Sowcr, Sogcr).
%
% The functions krw and krg are functions of Sw and Sg, respectively,
% whereas krow and krog are functions of So only.  The Jacobian of this
% relperm model is straightforward to compute.  The entries in the cell
% Jacobians are arranged column-wise as before.

   sw   = max(Swco, s(:,1));
   so   = s(:, 2);
   sg   = s(:, 3);
   null = zeros(size(s, 1), 1);

   nout = nargout;


   [kw{1:nout}]  = krw (sw);
   [kg{1:nout}]  = krg (sg);
   [kwo{1:nout}] = krow(so);
   [kgo{1:nout}] = krog(so);

   Som = min(Sowcr, Sogcr);

   SSo = (so-Som) ./(1-Swco-Som);
   SSw = (sw-Swco)./(1-Swco-Som);
   SSg =  sg      ./(1-Swco-Som);

   Fw  = kwo{1}./krocw./(1-SSw);   Fw(SSw==1) = 1;
   Fg  = kgo{1}./krocw./(1-SSg);   Fg(SSg==1) = 1;

   varargout{1} = [kw{1}, ...
                   krocw.*SSo.*Fw.*Fg,...
                   kg{1}];

   if nargout > 1,

      dFwdSw = kwo{1}./krocw./(1-Swco-Som)./( (1-SSw).^2 );
      dFwdSo = kwo{2}./krocw./(1-SSw);
      dFgdSg = kgo{1}./krocw./(1-Swco-Som)./( (1-SSg).^2 );
      dFgdSo = kgo{2}./krocw./(1-SSg);

      dFwdSw(SSw==1) = 0;
      dFwdSo(SSw==1) = 0;
      dFgdSg(SSg==1) = 0;
      dFgdSo(SSg==1) = 0;

      % BEWARE: entries in 3x3 Jacobian are arranged column-wise
      varargout{2} = [kw{2},                                     ...
                      krocw.*SSo.*dFwdSw.*Fg,                    ...
                      null,                                      ...
                                                                 ...
                      null,                                      ...
                      (krocw.*Fw.*Fg + krocw.*SSo.*dFwdSo.*Fg +  ...
                                       krocw.*SSo.*Fw.*dFgdSo) , ...
                      null,                                      ...
                                                                 ...
                      null,                                      ...
                      krocw.*SSo.*Fw.*dFgdSg,                    ...
                      kg{2}                                      ...
         ];
   end
end

%------------------------------------------------------------------------

function varargout = StoneII(s, Swco, krw, krow, krg, krog, krocw)
% This function implements the second method to compute oil relative
% permeabilities suggested by Stone.
%
% The oil relatve permeability kro is given by
%
%                [   krow          krog                     ]
%    kro = krocw·[ ( ----- + krw)·(----- + krg) - rkw - krg ]
%                [   krocw         krocw                    ]
%
% As in the previous three-phase models, the entries in the cell Jacobians
% are arranged column-wise as before.

   sw   = max(Swco, s(:,1));
   so   = s(:, 2);
   sg   = s(:, 3);
   null = zeros(size(s, 1), 1);

   nout = nargout;

   [kw{1:nout}]  = krw (sw);
   [kg{1:nout}]  = krg (sg);
   [kwo{1:nout}] = krow(so);
   [kgo{1:nout}] = krog(so);

   Fw  = kwo{1}./krocw + kw{1};
   Fg  = kgo{1}./krocw + kg{1};

   varargout{1} = [kw{1}, ...
                   max(krocw.*(Fw.*Fg - kw{1} - kg{1}), 0), ...
                   kg{1}];

   if nargout > 1,

      dFwdSo = kwo{2}./krocw;
      dFgdSo = kgo{2}./krocw;

      % BEWARE: entries in 3x3 Jacobian are arranged column-wise
      varargout{2} = [kw{2},                             ...
                      krocw.*(Fg.*kw{2} - kw{2}),        ...
                      null,                              ...
                      ...
                      null,                              ...
                      krocw.*(dFwdSo.*Fg + Fw.*dFgdSo) , ...
                      null,                              ...
                      ...
                      null,                              ...
                      krocw.*(Fw.*kg{2} - kg{2}),        ...
                      kg{2}                              ...
         ];

      % Set derivatives of oil relperm to zero if oil relperm is zero. The
      % relperm function is not smooth at this point.
      i = varargout{1}(:,2) == 0;
      varargout{2}(i, [2, 5, 8]) = 0;
   end
end

%--------------------------------------------------------------------------

function varargout = Simple(s, Swco, krw, krow, krg)
% This function implements the second method to compute oil relative
% permeabilities suggested by Stone.
%
% The oil relatve permeability kro is given by
%
%                [   krow          krog                     ]
%    kro = krocw·[ ( ----- + krw)·(----- + krg) - rkw - krg ]
%                [   krocw         krocw                    ]
%
% As in the previous three-phase models, the entries in the cell Jacobians
% are arranged column-wise as before.

   sw   = max(Swco, s(:,1));
   so   = s(:, 2);
   sg   = s(:, 3);
   null = zeros(size(s, 1), 1);

   nout = nargout;

   [kw{1:nout}]  = krw (sw);
   [kg{1:nout}]  = krg (sg);
   [kwo{1:nout}] = krow(so);

   varargout{1} = [kw{1}, kwo{1}, kg{1}];

   if nargout > 1,

      % BEWARE: entries in 3x3 Jacobian are arranged column-wise
      varargout{2} = [kw{2},  ...
                      null,   ...
                      null,   ...
                      ...
                      null,   ...
                      kwo{2}, ...
                      null,   ...
                      ...
                      null,   ...
                      null,   ...
                      kg{2}   ...
         ];

      % Set derivatives of oil relperm to zero if oil relperm is zero. The
      % relperm function is not smooth at this point.
      %i = varargout{1}(:,2) == 0;
      %varargout{2}(i, [2, 5, 8]) = 0;
   end

end

%--------------------------------------------------------------------------

function varargout = threephasecap(s, f1, f2, f3)
   [pc1{1 : nargout}] = f1(s(:,1));
   [pc2{1 : nargout}] = f2(s(:,2));
   [pc3{1 : nargout}] = f3(s(:,3));

   varargout{1}    = [ pc1{1}, pc2{1}, pc3{1} ];

   if nargout > 1,
      varargout{2} = [ pc1{2}, pc2{2}, pc3{2} ];
   end
   if nargout > 2,
      error(msgid('Pcap:NoHighDerivatives')                          , ...
           ['You appear to be implementing a highly advanced '       , ...
             'numerical scheme.\nUnfortunately, MRST does presently ', ...
             'not support higher order derivatives of '              , ...
             'capillary pressure'])
   end
end

%--------------------------------------------------------------------------

function varargout = twophase(s, f1, f2)
   varargout = cell([1, nargout]);

   [kr1{1 : nargout}] = f1(s);
   [kr2{1 : nargout}] = f2(s);

   varargout{1}    = [ kr1{1}, kr2{1} ];

   if nargout > 1,
      null         = zeros(size(s, 1), 1);
      varargout{2} = [ kr1{2}, null, null, kr2{2} ];
   end

   if nargout > 2,
      varargout{3} = [ kr1{3}, kr2{3} ];
   end
end

%--------------------------------------------------------------------------

function varargout = twophasecap(s, f1, f2)
   [pc1{1 : nargout}] = f1(s(:,1));
   [pc2{1 : nargout}] = f2(s(:,2));

   varargout{1}    = [ pc1{1}, pc2{1} ];

   if nargout > 1,
      varargout{2} = [ pc1{2}, pc2{2} ];
   end
   if nargout > 2,
      error(msgid('Pcap:NoHighDerivatives')                          , ...
           ['You appear to be implementing a highly advanced '       , ...
             'numerical scheme.\nUnfortunately, MRST does presently ', ...
             'not support higher order derivatives of '              , ...
             'capillary pressure'])
   end
end

%--------------------------------------------------------------------------

function varargout = null(s)
   varargout(1 : nargout) = { zeros(size(s)) };
end
