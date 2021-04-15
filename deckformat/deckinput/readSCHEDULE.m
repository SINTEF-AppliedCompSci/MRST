function [deck, varargout] = readSCHEDULE(fid, dirname, deck, varargin)
% Read schedule

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

   [schd, miss_kw] = get_state(deck);

   [start, ctrl, def_ctrl, cno] = setup(deck, varargin{:});

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section
      switch kw
         case 'BOX'
            warning('BOX:Ignored:Schedule', ...
                    'BOX keyword currently ignored in SCHEDULE section');

         case 'ENDBOX'
            warning('ENDBOX:Ignored:Schedule', ...
                    'ENDBOX keyword currently ignored in SCHEDULE section');

         case {'COMPDAT' , ...
               'COMPSEGS', ...
               'GCONINJE', ...
               'GCONPROD', ...
               'GECON'   , ...
               'GRUPTREE', 'GRUPNET', ...
               'RPTSCHED', 'RPTRST', ...
               'WCONHIST', 'WCONINJ', 'WCONINJE', 'WCONINJH', ...
               'WCONPROD', ...
               'WELSEGS' , ...
               'WELSPECS', ...
               'WELOPEN' , 'WELLOPEN', ...
               'WELTARG' , 'WELLTARG', ...
               'WGRUPCON', ...
               'WPOLYMER', ...
               'WSURFACT', ...
               'WSOLVENT', ...
               'WTEMP'}

            if ~def_ctrl
               def_ctrl = true;
               cno      = cno + 1;
               ctrl     = defaultControl(ctrl);
            end
            ctrl = readWellKW(fid, ctrl, kw);

         case 'VFPINJ'
            if ~def_ctrl
               def_ctrl = true;
               cno      = cno + 1;
               ctrl     = defaultControl(ctrl);
            end

            [tid, vfpinj]    = readVFPINJ(fid);
            ctrl.VFPINJ{tid} = vfpinj;   clear tid vfpinj

         case 'VAPPARS'
            if ~def_ctrl
               def_ctrl = true;
               cno      = cno + 1;
               ctrl     = defaultControl(ctrl);
            end
            data = readDefaultedRecord(fid, {'0', '0'});
            data = cellfun(@str2num, data);
            ctrl.VAPPARS = data;
             
         case 'VFPPROD'
            if ~def_ctrl
               def_ctrl = true;
               cno      = cno + 1;
               ctrl     = defaultControl(ctrl);
            end

            [tid, vfpprod]    = readVFPPROD(fid);
            ctrl.VFPPROD{tid} = vfpprod;   clear tid vfpprod

         case {'DATES', 'TSTEP'}
            if def_ctrl
               def_ctrl = false;
               schd.control = [schd.control; ctrl];
            end
            reader = str2func(['read', kw]);
            data   = reader(fid, start, schd.step.val);
            schd.step.control = [schd.step.control; ...
                                 repmat(cno, [numel(data), 1])];
            schd.step.val     = [schd.step.val; data];

         case {'ECHO', 'NOECHO'}
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         case 'DRSDT'
            data    = readDefaultedRecord(fid, {'NaN', 'ALL'});
            data{1} = sscanf(regexprep(data{1}, '[Dd]', 'e'), '%f');

            if ~def_ctrl
               def_ctrl = true;
               cno      = cno + 1;
               ctrl     = defaultControl(ctrl);
            end
            ctrl.(kw) = data;                                    clear data

            if isfinite(ctrl.DRSDT{1}) && ctrl.DRSDT{1}
               fprintf([...
                  '\n\n', ...
                  '@ The black-oil model currently implemented ' ,  ...
                  'does not support slower than\n@ instantaneous ' , ...
                  'dissolution of gas into oil.\n\n']);

               warning(msgid('DRSDT:Finite'), ...
                      ['Keyword DRSDT is currently ignored in MRST ', ...
                       'blackoil codes.']);
            end

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END'
            % Logical end of input deck.
            %
            in_section    = false;
            deck.SCHEDULE = schd;

         case 'INCLUDE'
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, schd, miss_kw);

            save_cno = cno;
            [deck, def_ctrl, cno] = ...
               readEclipseIncludeFile(@readSCHEDULE, fid, dirname, ...
                                      deck.RUNSPEC,                ...
                                      deck, 'ctrl', ctrl,          ...
                                      'def_ctrl', def_ctrl,        ...
                                      'cno', cno);

            % Prepare for additional reading.
            [schd, miss_kw] = get_state(deck);

            if def_ctrl || (cno > save_cno)
               ctrl = schd.control(end);
            end

            if def_ctrl
               schd.control = schd.control(1:end-1);
            end

         otherwise
            if ischar(kw)
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   if def_ctrl
      % Reading (or, possibly, file) ended whilst defining the current
      % 'control'.  Time stepping information is (presumably) elsewhere.
      %
      % Need to preserve 'control' information (well specification &c) for
      % caller.
      %
      schd.control = [schd.control; ctrl];
   end

   deck = set_state(deck, schd, miss_kw);

   if nargout > 1
      % Return def_ctrl and cno settings back to caller (i.e., ourselves)
      % to prepare for additional SCHEDULE reading.
      %
      varargout{1} = def_ctrl;
      varargout{2} = cno;
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function [start, ctrl, def_ctrl, cno] = setup(deck, varargin)
   opt = struct('def_ctrl', false,          ...
                'cno'     , 0,              ...
                'ctrl'    , defaultControl());
   opt = merge_options(opt, varargin{:});

   start    = deck.RUNSPEC.START;
   ctrl     = opt.ctrl;
   def_ctrl = opt.def_ctrl;
   cno      = opt.cno;
end

%--------------------------------------------------------------------------

function control = defaultControl(varargin)
   control = struct('WELSPECS', [], ...
                    'COMPDAT' , [], ...
                    'COMPSEGS', [], ...
                    'WCONINJ' , [], 'WCONINJE', [], 'WCONINJH', [], ...
                    'WCONPROD', [], 'WCONHIST', [], ...
                    'GCONINJE', [], 'GCONPROD', [], 'GECON'   , [], ...
                    'GRUPTREE', [], 'GRUPNET' , [], ...
                    'RPTSCHED', [], ...
                    'RPTRST'  , [], ...
                    'WELSEGS' , [], ...
                    'WGRUPCON', [], ...
                    'WPOLYMER', [], ...
                    'WSURFACT', [], ...
                    'WSOLVENT', [], ...
                    'VFPINJ'  , [], ...
                    'VFPPROD' , [], ...
                    'WTEMP'   , [], ...
                    'VAPPARS',  [0, 0], ...
                    'DRSDT'   , {{ inf, 'ALL' }});
   if nargin > 0
      control.WELSPECS = [control.WELSPECS; varargin{1}.WELSPECS];
      control.COMPDAT  = [control.COMPDAT;  varargin{1}.COMPDAT ];
      control.DRSDT    = varargin{1}.DRSDT;

      % Copy controls from existing configuration to honour 'WELTARG' &c.
      for ctrl = { 'WCONINJE', 'WCONINJH', 'WCONINJ', ...
                   'WCONPROD', 'WCONHIST', ...
                   'GCONINJE', 'GCONPROD', 'GECON' , ...
                   'GRUPTREE', 'GRUPNET' , ...
                   'COMPSEGS', ...
                   'WELSEGS' , ...
                   'RPTSCHED', 'RPTRST'  , ...
                   'VFPINJ'  , 'VFPPROD' }

         if ~ isempty(varargin{1}.(ctrl{1}))
            control.(ctrl{1}) = varargin{1}.(ctrl{1});
         end
      end
   end
end

%--------------------------------------------------------------------------

function data = readTSTEP(fid, varargin)
   data = readVector(fid, 'TSTEP', inf);
end

%--------------------------------------------------------------------------

function data = readDATES(fid, start, timespec)
   empty_record = @(rec) isempty(rec) || all(isspace(rec));

   getDate = @() strtrim(removeQuotes(readRecordString(fid)));

   dates = {};
   date  = getDate();
   while ~empty_record(date)
      dates = [dates; { date }];     %#ok  % Willfully ignore MLINT advice.
      date  = getDate();
   end

   % Months 'JLY' and 'JUL' are interchangeable in ECLIPSE.
   % Function 'datenum', however, requires 'JUL'.
   %
   dates = strrep(dates, 'JLY', 'JUL');
   data  = datenum(dates);  % dd mmm yyyy format automagically detected.
   data  = diff([start + sum(timespec); data]);  assert (all(data > 0));
end

%--------------------------------------------------------------------------

function [schd, miss_kw] = get_state(deck)
   schd    = deck.SCHEDULE;
   miss_kw = deck.UnhandledKeywords.SCHEDULE;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, schd, miss_kw)
   deck.SCHEDULE                   = schd;
   deck.UnhandledKeywords.SCHEDULE = unique(miss_kw);
end
