function [wellSols, time] = convertSummaryToWellSols(fn, unit)
% [wellSols, time] = convertSummaryToWellSols(fn, unit)
% Create wellSols with fields qOs, qWs, qGs and bhp from from eclipse 
% summary file fn. Supperted units are 'metric', 'field'.

   if nargin < 2
      warning('No unit given, assuming metric')
      unit = 'metric';
   end

   smry = readEclipseSummaryUnFmt(fn);

   % units:
   if ischar(unit)
      u = getUnits(unit);
   else
      u = unit;
   end

   [qOs, qWs, qGs, bhp, wns, time] = extract_quantities(smry, u);

   wellSols = assign_wellsols(smry, qOs, qWs, qGs, bhp, wns);
end

%--------------------------------------------------------------------------

function [qOs, qWs, qGs, bhp, wns, time] = extract_quantities(smry, u)
   tf   = ':+:+:+:+';   % special field with time-info

   wns  = setdiff(smry.WGNAMES, {tf, 'FIELD', ''}); % well-names
   time = reshape(convertFrom(smry.get(tf, 'TIME', ':'), u.t), [], 1);

   nw = numel(wns);
   nt = numel(time);

   [qOs, qWs, qGs, bhp] = deal(zeros([nt, nw]));

   for k = 1 : nw
      wn  = wns{k};
      akw = smry.getKws(wn);

      if ismember('WBHP', akw)
         bhp(:,k) = convertFrom(smry.get(wn, 'WBHP', ':'), u.p);
      end

      if ismember('WOPR', akw)
         qOs(:,k) = - convertFrom(smry.get(wn, 'WOPR', ':'), u.ql);
      end

      if ismember('WWPR', akw)
         qWs(:,k) = - convertFrom(smry.get(wn, 'WWPR', ':'), u.ql);
      end

      if ismember('WWIR', akw)
         qWs(:,k) = qWs(:,k) + ...
            reshape(convertFrom(smry.get(wn, 'WWIR', ':'), u.ql), [], 1);
      end

      if ismember('WGPR', akw)
         qGs(:,k) = - convertFrom(smry.get(wn, 'WGPR', ':'), u.qg);
      end

      if ismember('WGIR', akw)
         qGs(:,k) = qGs(:,k) + ...
            reshape(convertFrom(smry.get(wn, 'WGIR', ':'), u.qg), [], 1);
      end
   end
end

%--------------------------------------------------------------------------

function wellSols = assign_wellsols(smry, qOs, qWs, qGs, bhp, wns)
   nw = numel(wns);
   nt = size(smry.data, 2);

   ws = struct('name', '', 'bhp', 0, 'qOs', 0, 'qWs', 0, 'qGs', 0);

   wellSols = repmat({repmat(ws, [nw, 1])}, [nt, 1]);

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

      otherwise
         error(['Unit ', unit, ' not supported']);
   end
end
