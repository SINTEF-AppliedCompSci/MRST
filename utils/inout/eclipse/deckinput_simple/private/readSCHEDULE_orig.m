function sched = readSCHEDULE_orig(fid, dirname, start, varargin)
%Read (simplified version of) SCHEDULE section of ECLIPSE input deck.
%
% SYNOPSIS:
%   sched = readSCHEDULE_orig(fid, dirname, start)
%
% PARAMETERS:
%   fid     - Valid file identifier obtained from FOPEN pointing to the
%             beginning of the SCHEDULE data section of an ECLIPSE input
%             deck.
%
%   dirname - Complete directory name of file from which the input file
%             identifier 'fid' was derived through FOPEN.
%
%   start   - Simulation start date, typically entered using the 'RUNSPEC'
%             keyword 'START'.  Must be represented as a MATLAB DATENUM.
%             This will only be inspected if any report steps are specified
%             using the 'DATES' keyword.
%
% RETURNS:
%   sched -
%
% SEE ALSO:
%   `datenum`, `readGRDECL`.

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


   opt = struct('sched',          initSchedule(), ...
                'define_control', false,          ...
                'control_no',     0,              ...
                'control',        defaultControl());
   opt = merge_options(opt, varargin{:});

   sched          = opt.sched;
   define_control = opt.define_control;
   control_no     = opt.control_no;
   control        = opt.control;

   kw = getEclipseKeyword(fid);
   while ischar(kw),
      switch kw,
         case {'COMPDAT', 'WCONHIST', 'WCONINJE', 'WCONPROD', 'WELSPECS'},
            if ~define_control,
               define_control = true;
               control_no     = control_no + 1;
               control        = defaultControl(control);
            end
            control = readWellKW(fid, control, kw);

         case {'DATES', 'TSTEP'},
            if define_control,
               define_control = false;
               sched.control  = [sched.control; control];
            end
            reader = str2func(['read', kw]);
            data   = reader(fid, start, sched.step.val);
            sched.step.control = [sched.step.control; ...
                                  repmat(control_no, [numel(data), 1])];
            sched.step.val     = [sched.step.val; data];

         case 'INCLUDE',
            sched = readEclipseIncludeFile(@readSCHEDULE_orig, fid, dirname, start, ...
                                           'sched',          sched,            ...
                                           'define_control', define_control,   ...
                                           'control_no',     control_no,       ...
                                           'control',        control);

         case 'END',
            % (Logical) end of SCHEDULE section.  Often coincides with EOF.
            return;
      end

      kw = getEclipseKeyword(fid);
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function sched = initSchedule()
   sched.control = [];
   sched.step    = struct('control', [], 'val', []);
end

%--------------------------------------------------------------------------

function control = defaultControl(varargin)
   control = struct('WELSPECS', [], ...
                    'COMPDAT',  [], ...
                    'WCONHIST', [], ...
                    'WCONINJE', [], ...
                    'WCONPROD', []);
   if nargin > 0,
      control.WELSPECS = [control.WELSPECS; varargin{1}.WELSPECS];
      control.COMPDAT  = [control.COMPDAT;  varargin{1}.COMPDAT ];
   end
end

%--------------------------------------------------------------------------

function data = readTSTEP(fid, varargin) %#ok
   data = readVector(fid, 'TSTEP', inf);
end

%--------------------------------------------------------------------------

function data = readDATES(fid, start, timespec) %#ok
   empty_record = @(rec) isempty(rec) || all(isspace(rec));

   dates = {};
   date  = readRecordString(fid);
   while ~empty_record(date),
      dates = [dates; date];         %#ok  % Willfully ignore MLINT advice.
      date  = readRecordString(fid);
   end

   % Months 'JLY' and 'JUL' are interchangeable in ECLIPSE.
   % Function 'datenum', however, requires 'JUL'.
   %
   dates = strrep(dates, 'JLY', 'JUL');
   data  = datenum(dates);  % dd mmm yyy format automagically detected.
   data  = diff([start + sum(timespec); data]);  assert (all(data > 0));
end
