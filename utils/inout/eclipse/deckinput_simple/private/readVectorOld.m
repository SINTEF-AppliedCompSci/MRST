function v = readVectorOld(fid, field, nel)
%Input vector of floating point numbers from GRDECL file.
%
% SYNOPSIS:
%   v = readVector(fid, field, nel)
%
% PARAMETERS:
%   fid   - File identifier (as defined by FOPEN) of GRDECL file open for
%           reading.  Assumed to point to a seekable (i.e. physical) file
%           on disk, not (e.g.) a POSIX pipe.
%
%   field - Name (string) identifying the GRDECL keyword (field) currently
%           being processed.  Used for error identification/messages only.
%
%   nel   - Number of elements to read from input stream.  As a special
%           case, the caller may pass nel==INF (or nel=='inf') to read as
%           much as possible.  In this case it is imperative that the input
%           vector be terminated by a '/' character.
%
% RETURNS:
%   v     - Vector, length 'nel', of floating point numbers defining the
%           contents of GRDECL keyword 'field'.  If nel==INF, NUMEL(v) is
%           the number of vector elements read before the terminating slash
%           character.

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


% Simple vector read.  Assume the vector is given explictly by listing all
% of its values separated by whitespace, possibly interspersed by comments.
%

if ischar(nel) && strcmpi(nel, 'inf'), nel = inf; end
scanner = getScanner();
opt = getScanOpt(nel);
cs  = opt{end};   % CommentStyle -- relies on knowledge of 'getScanOpt'.
C   = scanner(fid, '%f', opt{:});

pos = ftell(fid);
lin = fgetl(fid);
while ischar(lin) && (isempty(lin) || matches(lin, ['^\*|^\s+$|^\s*', cs]))
   if strncmp('*', lin, 1)
      % Repeat description of the form
      %      N*value
      % with 'N' a positive integer.  There must be no whitespace
      % immediately before or after the literal asterisk ('*').
      %
      % This format allows considerable economy of file space, e.g., when
      % listing uniform porosity values such as
      %
      %   PORO
      %     5000*0.3 /
      %
      % for a vector of five thousand repetitions of the value 0.3.  The
      % MATLAB equivalent is
      %
      %   PORO = REPMAT(0.3, [5000, 1]);
      %
      % This format is most likely encountered when processing GRID section
      % keywords such as 'DXV', 'PERMX', 'PORO', 'ACTNUM' and, possibly,
      % 'ZCORN' and SCHEDULE section keywords such as 'TSTEP'.
      %

      % Validate repeat description format.
      %
      % Algorithm:
      %   1) Read character string starting from (pos - 1) up to first
      %      blank (white-space) or termination ('/') character.
      %   2) Check if first two characters are exactly one digit and an
      %      asterisk respectively.
      %   3) If so, check that the remainder of the character string can be
      %      parsed as a floating point value without any parse failures.
      %
      fseek(fid, pos - 1, 'bof');
      check_string = fscanf(fid, '%[^/ \n\r\v\t\f]', 1);

      [count, count, err, nchar] = ...
         sscanf(check_string(3:end), '%f');                     %#ok<ASGLU>

      if ~ matches(check_string(1:2), '\d\*') || ...
         (count ~= 1) || ~ isempty(err)       || ...
         (nchar ~= (numel(check_string) - 1))

         error('readGRDECL:RepeatDescr:Malformed',                    ...
               'Incorrect repeat description detected in keyword %s', ...
               field);
      end

      % Validate repeat count value.
      nrep = C{1}(end);
      if ~isnumeric(nrep)             || ...
         nrep < 1                     || ...
         nrep > nel - numel(C{1}) + 1 || ...
         nrep ~= fix(nrep)
         error('readGRDECL:RepeatCount:OutOfBounds',             ...
               ['Unexpected repeat count in keyword %s.\n',      ...
                'Expected integral repeat count in [1 .. %d], ', ...
                'but found %d.'], field, nel - numel(C{1}) + 1, nrep);
      end

      % Append repeated value to vector.
      % 1) Position file descriptor directly after literal asterisk.
      %
      fseek(fid, pos + 1, 'bof');

      % 2) Read a single floating point value (whence FTELL(fid) will
      %    be on first whitespace, termination character ('/'), or,
      %    if the input is malformed, EOF).
      %
      [val, count] = fscanf(fid, '%f', 1);   assert (count == 1);

      % 3) Append this value to the vector whilst excluding the repeat
      %    present in C{1}(end).
      %
      C = { [reshape(C{1}(1 : end-1), [], 1); ...
             repmat(val, [nrep, 1])]        };

      % 4) Count new number of elements and update expected read count.
      %
      nv  = numel(C{1});
      opt = getScanOpt(nel, nv);

      % Read next portion of vector data (possibly all remaining data).
      C1 = scanner(fid, '%f', opt{:});

      % Append this data to existing array.
      C  = { [ reshape(C{1}, [], 1);   reshape(C1{1}, [], 1) ] };
   end
   pos = ftell(fid);
   lin = fgetl(fid);
end

iseof           = feof(fid);
[errmsg, errno] = ferror(fid);
readAll = (~isinf(nel) && (numel(C{1}) == nel)) | isinf(nel);

if (iseof || (ischar(lin) && matches(lin, '^\s*/'))) && readAll
   v = C{1};
elseif ischar(lin) && matches(lin, '^\s*/')
   %fclose(fid);
   v  = C{1};
   warning(['readGRDECL:Vector_', field, ':PrematureTermination'],    ...
         ['Detected termination character ''/'' at position %d.\n', ...
          'Scanning of keyword ''%s'' not completed.'], pos, field);
elseif ischar(lin) && matches(lin, '^[A-Z]+')
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':OtherKeyword'],           ...
         ['Encountered non-numeric data ''%s'' at position %d.\n', ...
          'Scanning of keyword ''%s'' not completed.\n',           ...
          'Missing slash (/) in input?'],                          ...
          lin, pos, field);
elseif ischar(lin) && matches(lin, '^[a-z!;=]+')
   error(['readGRDECL:Vector_', field, ':Unexpected'],             ...
         ['Encountered non-numeric data ''%s'' at position %d.\n', ...
          'Scanning of keyword ''%s'' not completed.'],            ...
          lin, pos, field);
elseif iseof
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':EOF'],               ...
         'End of file before complete read of keyword ''%s''.', field);
elseif errno ~= 0
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':SystemError'], ...
         ['Read error while processing keyword %s.\n',  ...
          'Input system reports ''%s''.'], field, errmsg);
else
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':Malformed'],     ...
         ['Unexpected input while reading keyword %s.\n', ...
          'Found ''%s'' at position %d.'], field, lin, pos);
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function scan = getScanner()
%{
if exist('textscan', 'builtin'),
   scan = @textscan;
else
%}
   scan = @myscanner;
%end

%--------------------------------------------------------------------------

function opt = getScanOpt(nel, varargin)
opt = {'CommentStyle', '--'};
if nargin > 1 && isnumeric(varargin{1})
   nel = nel - varargin{1};
end
if ~isinf(nel)
   opt = [{nel}, opt];
end

%--------------------------------------------------------------------------

function b = matches(str, pat)
b = ~isempty(regexp(str, pat, 'once'));

%--------------------------------------------------------------------------

function C = myscanner(fid, format, varargin)
%Scanner function supporting a subset of TEXTSCAN features.

% This function is needed when TEXTSCAN is unavailable (e.g., Octave or old
% (pre-R14) releases of M).
%
% Algorithm:
%   1) Attempt to read as much data as possible using FSCANF.  This handles
%      most common cases (all data listed separately, no comments).
%   2) If we encounter a comment, as identified by 'opt.CommentStyle', then
%      scan individual lines until all comments have been scanned.
%   3) If we encounter 'end-of-record' (i.e., a '/' character), then
%      terminate all scanning and trust caller to detect any errors.

% Setup:
%   Define comment style (usually '--' to exclude a single line of input)
%   Define termination strings (regexp's).
nel = inf;
if mod(numel(varargin), 2) == 1
   if readnum_ok(varargin{1})
      nel = varargin{1};
   end
   varargin = varargin(2 : end);
end
opt = struct('CommentStyle', '--');
opt = merge_options(opt, varargin{:});
if iscell(opt.CommentStyle)
   assert (numel(opt.CommentStyle) == 1);
   opt.CommentStyle = opt.CommentStyle{1};
end
assert (ischar(opt.CommentStyle));
patt_comment   = ['^', opt.CommentStyle];
patt_terminate = '^\*|^/';

% Algorithm step 1.  Consume as much data as possible.
[v, cnt] = fscanf(fid, format, nel);
n = 0;
C = reshape(v, [], 1);
read_complete = n + cnt == nel;
while ~read_complete
   % Short read.  Encountered comment or "premature" termination (usually
   % when ISINF(nel)).  Must either consume a (block of) comment line(s)
   % and all subsequent data or detect end-of-record in which case we need
   % to terminate.

   % Algorithm:
   %   1) Record position in file to allow caller to inspect input stream
   %      if we need to terminate.
   %   2) Inspect the next few character to determine if this is a comment
   %      or 'end-of-record'.
   %   3) Act accordingly.
   pos = ftell(fid);
   lin = fgetl(fid);
   if strncmp(opt.CommentStyle, '--', 2) && matches(lin, '^-')
      % CommentStyle is '--' (i.e., following Eclipse convention).  In
      % MATLABs prior to 7.10 (R2010a), FSCANF consumed the first '-'
      % character thinking that the next token would be (negative) floating
      % point number and failed when the second '-' character was
      % encountered.  In MATLABs 7.10 and later, this behaviour has been
      % corrected and FSCANF no longer consumes characters from the input
      % stream that cannot be converted to outputs according to the scan
      % format.
      %
      % Support both behaviours by backing up to one character before 'lin'
      % (unless pos==0, in which case we back up to 'pos'), before reading
      % a single character.  This one character should either be a space,
      % as defined by ISSPACE, or a single dash that (likely) introduces a
      % comment.
      fseek(fid, pos - (pos > 0), 'bof');
      c   = fscanf(fid, '%c', 1);
      lin = fgetl(fid);

      assert (isspace(c) || c == '-', ...
             ['Unexpected character ''%c'' at position %d in file ', ...
              '''%s''. Expected space or comment character.'], ...
              c, pos, fopen(fid));

      if c == '-'
         lin = [c, lin];  %#ok
      end
   end
   n = n + cnt;
   if matches(lin, patt_comment)
      % This is a comment.  The above FGETL has already consumed all
      % character up to (and including) EOL.  Attempt to read new data with
      % regular FSCANF approach.  This will fail if the next line is
      % another comment, but then the enclosing WHILE will bring us back
      % here on the next iteration...
      [v, cnt] = fscanf(fid, format, nel - n);
      C = [C; reshape(v, [], 1)];  %#ok
      read_complete = n + cnt == nel;
   elseif matches(lin, patt_terminate)
      % This is a termination character (i.e., the '/' character) or a
      % repeat character (i.e., the '*' character).  Back up to beginning
      % of 'lin' to allow caller to rediscover the special character and
      % take appropriate action.  Our work here is done.
      fseek(fid, pos, 'bof');
      read_complete = true;
   else
      fname = fopen(fid);
      fclose(fid);
      error('readVector:Input:Unexpected', ...
           ['Unexpected input ''%s'' encountered at position ''%d''', ...
            ' whilst scanning file ''%s''.'], lin, pos, fname);
   end
end

% Convert to CELL to conform to semantics of TEXTSCAN.
C = { C };

%--------------------------------------------------------------------------

function b = readnum_ok(n)
b = isnumeric(n) && ~isempty(n) && isfinite(n);
