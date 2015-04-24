function [relperm, pc] = initSatnumRelPerm(T, satnum, imbnum, varargin)
% Construct two-phase relperm evaluation function.
%
% SYNOPSIS:
%   relperm = initSatnumRelPerm(table)
%
% DESCRIPTION:
%   Construct two-phase relperm evaluation function while recognizing
%   Eclipse keywords satnum and imbnum. Satnum specifies the primary
%   drainage curves, while imbnum specifies the primary imbibition curves.
%   If the drainage or imbibition process i reversed at some point, the
%   relperm for the non-wetting phase is computed from a scanning curve
%   using Killoughs hysteresis model. We assume relative permeability
%   hysteresis only in the non-wetting phase and the relperm for the
%   wetting phase is always computed from the drainage curve.
%
%   No hysteresis will occur if satnum = imbnum or if only the drainage
%   curve is provided (i.e. the imbnum is not specified).
%
% PARAMETERS:
%   table   - PVT/Relperm table structure as defined by function
%             'readpvt'.
%
% RETURNS:
%   relperm - Function for evaluating relative permeability curves
%             specified by SWOF or SWFN and SGFN.  Specifically, the call
%
%                [kr, dkr] = relperm(s)
%
%             evaluates the relative permeabilities (kr) and the
%             differentiated relative permeabilities (dkr) of the current
%             fluid composition (s, saturation).
%
% REMARK:
%   When modeling hysteresis it is assumed that rSol.s is the wetting phase
%   saturation. Hysteresis is only modeled in the non-wetting phase.
%   dkr is at the moment computed from the drainage curves.
%
% SEE ALSO:
%   initSatnumFluid, readRelPermTable, swof, sgof, swfn, sgfn.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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



%opt = struct('verbose', mrstVerbose);
%opt = merge_options(opt, varargin{:});

if isfield(T, 'swof')
   %[krw, kro, Swco, krocw, Sorw, pcw] = swof(opt.table);
   %[krw, krn, Swco, tmp, Sncr, pc] = swof(T.swof);
   [krw, krn, krocw, Swco, Swcr, Sncr, Swmax, pc] = swof(T.swof);

   % kr_vec{1}, kr_vec{2}, -(kr_vec{4}-1), .., kr_vec{3}

   Snmax = 1 - Swco;


elseif all(isfield(T, 'sgfn'), isfield(T, 'swfn'));
   [krw, Swco, pc] = swfn(T.swfn);
   %[kr_vec{1}, ..]

   [krn, Sncr, Snmax]  = sgfn(T.sgfn, Swco);
   %[kr_vec{2}, kr_vec{3}, kr_vec{4} ..]
   Snmax = min(1-Swco, Snmax); % OK? end points must meet, see manual

else
   error(msgid('..:MissingKeyword'), ...
      'Table must contain ''swof'' field or ''swfn'' and ''sgfn'' field');
end

relperm = @(s, x) twophaseRelPerm(s, x, satnum, imbnum, ...
                               krw, krn, Sncr, Snmax);

end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------
% Return relperm according to satnum of the cells and
% drainage/imbibition/hysteresis
function varargout = twophaseRelPerm(s, sol, satnum, imbnum, krw, krn, ...
      Sncr, Snmax)

   if(size(s,1)~=numel(satnum))
      % only allow use of non-matching satnum regions for estimating dt in
      % 'explicitTransport'
      %warning('saturation and satnum do not match we use region 1')
      st = dbstack;
      assert(any(cellfun(@(x) strcmp(x, 'estimate_dt'), {st.name})))
      satnum = ones(size(s,1),1);
      imbnum = satnum;
      sol.extSat = sol.s(:,1);
   end


   if nargout > 1,
      dkr = zeros([size(sol.s,1), 2]);  %NB: endre
   end

   scan  = false(numel(sol.s(:,1)),1);
   imb   = false(numel(sol.s(:,1)),1);

   %       if ~isfield(sol, 'satnum') %for use in e.g initTransport
   %          if opt.verbose,
   %             disp('Using default value: first saturation table'),end;
   %          sol.satnum = true(size(sol.s));
   %          sol.satnumActive = 1;
   %   elseif isfield(sol, 'minSat')
   %    sn_max = snmax(sol.satnumActive);

   isMin = sol.s(:,1) == sol.extSat(:,1);
   isEnd = sol.extSat(:,1) < 1 - Snmax(satnum)' + eps; % <= ?

   scan(~isMin & ~isEnd) = true;
   imb(~isMin & isEnd)   = true;

   % no hysteresis if relperm drainage = relperm imibition
   scan(satnum == imbnum) = false;
   imb(satnum == imbnum) = false;

   % DRAINAGE CURVE:
   % Always compute curves for drainage: used for wetting phase and as
   % default value for non-wetting phase. At the moment dkr is also
   % computed from the drainage curves for cells where the relperm is taken from the
   % scanning curve.

   %NB: use modified saturations?!

   [krwd{1 : nargout}] = evalMultipleRegions(krw, satnum, sol.s(:,1));
   [krnd{1 : nargout}] = evalMultipleRegions(krn, satnum, 1-sol.s(:,1));
   kr = [krwd{1}, krnd{1}];

   if nargout > 1,
      dkr(:, 1) = krwd{2};
      dkr(:, 2) = krnd{2};
   end

   % SCANNING CURVE:
   % Compute non-wetting relative permeability for cells where the
   % drainage or imbibition process has been reversed at some point.
   %
   % Use Killough's method for relative permeability hysteresis to
   % compute relative permeabilites from scanning curves for the
   % non-wetting phase. See ECLIPSE techincal documentation ch. 32
   % Hysterisis for details.

   if any(scan)
      ix = scan;

      sncrd = Sncr(satnum(ix))';  % critical saturation drainage curve
      sncri = Sncr(imbnum(ix))';  % critical saturation imbibition curve
      snmax = Snmax(satnum(ix))'; % maximum saturation value for n-w phase (should be same for satnum and imbnum)
      shy = 1-sol.extSat(ix,1);     % max non-wetting saturation in run

      % Calculate trapped critical saturation sncrt:
      C = 1./(sncri-sncrd)-1./(snmax-sncrd);
      a = 0.1; % for robustness
      A = 1 +a.*(snmax-shy);
      sncrt = sncrd+(shy-sncrd)./(A+C.*(shy-sncrd));

      % Normalized saturation
      snorm = sncri+(1-sol.s(ix)-sncrt).* ...
         (snmax-sncri)./(shy-sncrt);

      krni_s     = evalMultipleRegions(krn, imbnum(ix), snorm);
      krnd_hy  = evalMultipleRegions(krn, satnum(ix), shy);
      krnd_max = evalMultipleRegions(krn, satnum(ix), snmax);

      % use formula to calculate relperm for non-wetting phase
      kr(ix, 2) = krni_s.*krnd_hy./krnd_max;

      %% NB: does not compute derivatives for points on scanning curve.
      %% TODO?

   end

   % IMBIBITION CURVE
   % Compute relative permeability in cells where the imibition process
   % has started at Swmin and has not been reversed
   if any(imb)
      ix = imb;
      [krni{1 : nargout}] = evalMultipleRegions(krn, imbnum(ix), 1-sol.s(ix,1));
      kr(ix, 2) = krni{1};

      if nargout > 1,
         dkr(ix, 2) = -krni{2};
      end
   end
   varargout{1}   = kr ;

   if nargout > 1,  varargout{2} = dkr;

   end
end

