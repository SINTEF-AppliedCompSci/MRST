function [wellSols, time] = convertSummaryToWellSols(fn, unit)
% [wellSols, time] = convertSummaryToWellSols(fn, unit)
% Create wellSols with fields qOs, qWs, qGs and bhp from from eclipse
% summary file fn. Supperted units are 'metric', 'field'.

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

    [smry, spec] = get_summary(fn);

    if nargin == 1, unit = []; end

    [qOs, qWs, qGs, bhp, wns, time] = ...
       extract_quantities(smry, extract_units(spec, unit));

    wellSols = assign_wellsols(smry, qOs, qWs, qGs, bhp, wns);
end

%--------------------------------------------------------------------------

function [smry, spec] = get_summary(fn)
    if isstruct(fn)
        [smry, spec] = deal(fn, []);
    else
        [smry, spec] = readEclipseSummaryUnFmt(fn);
    end
end

%--------------------------------------------------------------------------

function u = extract_units(spec, unit)
    % units:
    if ischar(unit)
        u = getUnits(unit);
    elseif ~isempty(spec)
        u = getInteheadUnit(spec);
    elseif ~isempty(unit)
        u = unit;
    else
        warning('Unit:Defaulted', 'No unit given, assuming metric');
        u = getUnits('metric');
    end
end

%--------------------------------------------------------------------------

function [qOs, qWs, qGs, bhp, wns, time] = extract_quantities(smry, u)
   tf   = ':+:+:+:+';   % special field with time-info

   wns  = get_well_names(smry); % well-names
   t = smry.get(tf, 'TIME', ':');
   time = reshape(convertFrom(t, u.t), [], 1);

   nw = numel(wns);
   nt = max(numel(time), size(t, 2));

   [qOs, qWs, qGs, bhp] = deal(zeros([nt, nw]));

   foundOil = true;
   foundGas = true;
   foundWat = true;
   for k = 1 : nw
      wn  = wns{k};
      akw = smry.getKws(wn);

      if ismember('WBHP', akw)
         bhp(:,k) = convertFrom(smry.get(wn, 'WBHP', ':'), u.p);
      end

      if ismember('WOPR', akw)
         qOs(:,k) = - convertFrom(smry.get(wn, 'WOPR', ':'), u.ql);
      else
         foundOil = false;
      end

      if ismember('WGPR', akw)
         qGs(:,k) = - convertFrom(smry.get(wn, 'WGPR', ':'), u.qg);
      elseif ismember('WGOR', akw)
          % We got gas-oil surface ratio, recompute gas-rate from that
         qGs(:,k) = qOs(:,k).*reshape(smry.get(wn, 'WGOR', ':'), [], 1);
      else
         foundGas = false;
      end

      if ismember('WWPR', akw)
         qWs(:,k) = - convertFrom(smry.get(wn, 'WWPR', ':'), u.ql);
      elseif ismember('WWCT', akw)
         wcut = reshape(smry.get(wn, 'WWCT', ':'), [], 1);
         qWs(:,k) = wcut.*qOs(:,k)./(1-wcut);
      else
          foundWat = false;
      end

      if ismember('WWIR', akw)
         qWs(:,k) = qWs(:,k) + ...
            reshape(convertFrom(smry.get(wn, 'WWIR', ':'), u.ql), [], 1);
      end

      if ismember('WGIR', akw)
         qGs(:,k) = qGs(:,k) + ...
            reshape(convertFrom(smry.get(wn, 'WGIR', ':'), u.qg), [], 1);
      end
      if mrstVerbose()
          if ~foundWat
              fprintf('I was not able to reconstruct water surface rate for well %d. Will return zeros.\n', k);
          end
          if ~foundOil
              fprintf('I was not able to reconstruct oil surface rate for well %d. Will return zeros.\n', k);
          end
          if ~foundGas
              fprintf('I was not able to reconstruct gas surface rate for well %d. Will return zeros.\n', k);
          end
      end
   end
end

%--------------------------------------------------------------------------

function wns = get_well_names(smry)
   wkw = smry.KEYWORDS(~cellfun('isempty', regexp(smry.KEYWORDS, '^W')));

   wns = {};
   for kw = reshape(wkw, 1, [])
      wns = [ wns ; reshape(smry.getNms(kw{1}), [], 1) ];       %#ok<AGROW>
   end

   wns = unique(wns);
   wns = setdiff(wns, ':+:+:+:+');
end

%--------------------------------------------------------------------------

function wellSols = assign_wellsols(smry, qOs, qWs, qGs, bhp, wns)
   nw = numel(wns);
   nt = size(smry.data, 2);

   ws = struct('name', '', 'bhp', 0, 'qOs', 0, 'qWs', 0, 'qGs', 0);

   wellSols = repmat({repmat(ws, [1, nw])}, [nt, 1]);

   for kt = 1:nt
      for kw = 1:nw
         wellSols{kt}(kw).name = wns{kw};
         wellSols{kt}(kw).bhp  = bhp(kt, kw);
         wellSols{kt}(kw).qOs  = qOs(kt, kw);
         wellSols{kt}(kw).qWs  = qWs(kt, kw);
         wellSols{kt}(kw).qGs  = qGs(kt, kw);
         wellSols{kt}(kw).sign = sign(qWs(kt, kw)+qOs(kt, kw)+qGs(kt, kw));
      end
   end
end

%--------------------------------------------------------------------------

function u = getUnits(unit)
   switch lower(unit)
      case 'metric'
         u.p  = barsa;
         u.ql = meter^3/day;
         u.qg = meter^3/day;
         u.t  = day;

      case 'field'
         u.p  = psia;
         u.ql = stb/day;
         u.qg = 1000*ft^3/day;
         u.t  = day;

      case 'lab'
         u = struct('p' , atm, ...
                    'ql', (centi*meter)^3/hour, ...
                    'qg', (centi*meter)^3/hour, ...
                    't' , hour);

      case {'pvt-m', 'pvt_m'}
         u = struct('p' , atm, ...
                    'ql', meter^3/day, ...
                    'qg', meter^3/day, ...
                    't' , day);

      otherwise
         error(['Unit ', unit, ' not supported']);
   end
end

%--------------------------------------------------------------------------

function u = getInteheadUnit(spec)
   ustring = { 'METRIC', 'FIELD', 'LAB', 'PVT-M' };
   u = getUnits(ustring{ spec.INTEHEAD.values(1) });
end
