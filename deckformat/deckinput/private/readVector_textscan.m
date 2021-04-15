function vec = readVector_textscan(fid, field, nel)
%Input vector of floating point numbers from ECLIPSE input file.
%
% SYNOPSIS:
%   v = readVector_textscan(fid, field, nel)
%
% PARAMETERS:
%   fid   - File identifier (as defined by FOPEN) of ECLIPSE input file
%           open for reading.  Assumed to point to a seekable (i.e.,
%           physical) file on disk, not (e.g.) a POSIX pipe.
%
%   field - Name (string) identifying the keyword (field) currently being
%           processed.  Used for error identification/messages only.
%
%   nel   - Number of elements to read from input stream.  As a special
%           case, the caller may pass nel==INF (or nel=='inf') to read as
%           much as possible.  In this case it is imperative that the input
%           vector be terminated by a '/' character.
%
% RETURNS:
%   v     - Vector, length 'nel', of floating point numbers defining the
%           contents of ECLIPSE keyword 'field'.  If nel==INF, NUMEL(v) is
%           the number of vector elements read before the terminating slash
%           character.
%
% SEE ALSO:
%   `fopen`, `fseek`, `textscan`.

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

   if isnan(nel)
      error('Argument:Invalid', 'Element Count Cannot be NaN');
   end

   if ischar(nel) && strcmpi(nel, 'inf')
      nel = inf;
   end

   % Optimistic approach.  Try to parse vector that consists entirely of
   % separate floating-point values (i.e., no repeat counts).  This is the
   % fastest approach and often, though not always, useful for the very
   % largest vectors such as ZCORN.
   [vec, pos, complete] = read_vector_elements_simple(fid, nel);

   if ~complete
      % Failed to read a complete vector.  Most likely cause is that we
      % encountered a repeat count character ('*') which signals that we
      % need to perform more elaborate input parsing.  First verify that we
      % didn't encounter a read error.
      [msg, errno] = ferror(fid);

      if errno == 0
         % No read error.  This is likely to be the case of repeat counts.
         % We must, however, verify that we were able to read at least one
         % floating-point value (the repeat count).
         assert (isfloat(vec) && ~isempty(vec), 'Internal Logic Error');

         % Tokenize the remainder of the vector's elements.
         elm = read_vector_elements(vec(end), pos, fid);

         % Parse those elements to numeric values and concatenate to
         % previously read elements.  Recall that the last read element of
         % 'vec' is the first repeat count, so we must NOT include that
         % value in the final output returned to the caller.
         vec = [ vec(1 : (end - 1))         ; ...
                 parse_vector_elements(elm) ];
      else
         error('Read:Error', 'Input Error (TEXTSCAN): %s', msg);
      end
   end

   finalize_reading_and_check_state(vec, fid, field, nel)
end

%--------------------------------------------------------------------------

function [v, pos, complete] = read_vector_elements_simple(fid, nel)
   % Simple vector read.  Try to consume a vector for which each element
   % value is listed separately without any repeat counts.

   % Defer the heavy lifting and comment handling to TEXTSCAN.  In
   % principle this is FSCANF(FID, '%f'), but TEXTSCAN is easier to use.
   [arr, pos] = textscan(fid, '%f', 'CommentStyle', '--', ...
                         'WhiteSpace', ' \r\n\t', ...
                         'MultipleDelimsAsOne', true);

   % Output ('arr') is single-element cell array of array of floating-point
   % values.  Unpack to get the actual values.
   v = arr{1};

   % Check if we completed the entire vector read.
   if isfinite(nel)
      % Finite number of vector elements.  Check if we got all of them.
      complete = numel(v) == nel;
   elseif feof(fid)
      % Non-finite element count (i.e., called as readVector(..., inf)).
      % Treat EOF as '/' for ISINF(nel).
      complete = true;
   else
      % Non-finite element count (i.e., called as readVector(..., inf)).
      % Not at EOF.  Vector read complete if next character is terminator
      % ('/').
      complete = strcmp(fscanf(fid, ' %c', 1), '/');

      % Reverse back one character to enable surrounding logic to identify
      % termination criterion.
      fseek(fid, pos, 'bof');
   end
end

%--------------------------------------------------------------------------

function tok = read_vector_elements(rpt, pos, fid)
   % Elaborate tokenization of vectors that contain repeat counts of the
   % form
   %      N*value
   % with 'N' a positive integer.  There must be no whitespace immediately
   % before or after the literal asterisk ('*').
   %
   % This format allows considerable economy of file space, e.g., when
   % listing uniform porosity values such as
   %
   %   PORO
   %     5000*0.3 /
   %
   % for a vector of five thousand repetitions of the value 0.3.  The MATLAB
   % equivalent is
   %
   %   PORO = REPMAT(0.3, [5000, 1]);
   %
   % This format is most likely encountered when processing GRID section
   % keywords such as 'DXV', 'PERMX', 'PORO', 'ACTNUM' and, possibly,
   % 'ZCORN' and SCHEDULE section keywords such as 'TSTEP'.

   % Consistency checking: Verify that encountering a repeat count ('*')
   % character actually prompted the more involved input routine.  Input
   % 'pos' is directly before the character that could not be converted
   % during the simplified vector read.  Seek to that position, read a
   % single character and verify that that single character actually is the
   % repeat designator.
   fseek(fid, pos, 'bof');
   ast = fscanf(fid, '%c', 1);
   if ~strcmp(ast, '*') || mod(rpt, 1) ~= 0
      error('RepeatDesignator:Erroneous', ...
           ['Incorrect Repeat Descriptor (''%s'': %d).  ', ...
            'Expected N*, but got ''%f%c'''], fopen(fid), pos, rpt, ast);
   end

   % Tokenize relevant portion of the input stream.  Read everything up to
   % the vector terminator '/' (or FEOF or input error), split on blanks
   % (space, newline, tab) and ignore comment lines ('--' designator).
   tok = textscan(fid, '%[^/ \r\n\t]', 'CommentStyle', '--', ...
                  'MultipleDelimsAsOne', true);

   % Output from TEXTSCAN is single-element cell array of cellstring of
   % tokens.  Get the actual tokens.
   tok = tok{1};

   % Re-insert repeat count that was erroneously intepreted as a vector
   % value into the first token to allow the parsing step to treat this
   % complete token correctly.
   tok{1} = sprintf('%d%c%s', rpt, ast, tok{1});
end

%--------------------------------------------------------------------------

function vec = parse_vector_elements(val)
   % Using TEXTSCAN here means we automatically handle Fortran-style double
   % precision exponent characters (1.2D+3, 4.3d-2).  Function TEXTSCAN
   % requires a single string (character vector), so we need to convert the
   % input cellstring to that format (N-by-C CHAR).  Take special care to
   % separate each token (single blank character between tokens).
   cell1     = @(c) c{1};
   mkchar    = @(c) [ char(c).' ; repmat(' ', [1, numel(c)]) ];
   to_double = @(c) cell1(textscan(mkchar(c), '%f'));

   count = ones([numel(val), 1]);

   % Input (val) is cellstring with each element being a token of the form
   % 1.23D+4 or 123*-.4567e-8 .  We need to identify the tokens that have
   % repeat counts and translate those to RLDECODE arguments.

   rpt = regexp(val, '(\d+)\*(\S+)', 'tokens');
   i   = ~cellfun('isempty', rpt);

   has_rpt = any(i);

   if has_rpt
      % Some vector elements have repeat counts.  Decode those and replace
      % the tokens with their values only.

      % 1) Extract tokens with repeat counts.
      rpt = vertcat(rpt{i});

      % 2) Translate those tokens to M-by-2 cellstring.
      rpt = vertcat(rpt{:});

      % 3) Convert string representations of repeat counts to numeric value
      count(i) = to_double(rpt(:,1));

      % 4) Replace repeat count tokens with the value portion (i.e., the
      %    character sequence after the '*' character) only.
      val(i) = rpt(:,2);
   end

   % Convert value string tokens to numeric type.
   vec = to_double(val);

   if has_rpt
      % Some vector elements have repeat counts.  Expand those.
      vec = rldecode(vec, count, 1);
   end
end

%--------------------------------------------------------------------------

function finalize_reading_and_check_state(vec, fid, field, nel)
   if isfinite(nel) && (numel(vec) ~= nel)
      % Finite expected element count, but we did not get that number.
      % Could be something resembling TOPS that has an implied copy
      % operation built into the specification.  Issue a warning here.
      if ~strcmp(field, 'TOPS')
          warning('VectorSize:Mismatch', ...
                  'Failed to Input Keyword Vector ''%s''', field);
      end
   else
      % Element reading complete for this vector.  Finish up by skipping
      % the terminator ('/') and anything following that up to and
      % including the next EOL mark.
      pos   = ftell(fid);
      lin   = fgetl(fid);
      fname = fopen(fid);

      if ischar(lin)
         % Got a set of characters.  If the first character is the vector
         % terminator ('/'), then everything is fine.  Otherwise, issue a
         % warning that we didn't find the expected terminator.

         if isempty(regexp(lin, '^\s*/', 'once'))
            warning('Unexpected:Termination', ...
                   ['Expected Vector Terminator (''/''), ', ...
                    'but Got ''%s'' (''%s'': %d)'], ...
                    lin, fname, pos);
         end

      elseif ~feof(fid)
         % Didn't find any characters after the vector elements and we're
         % not at EOF.  This is a (highly) unexpected input failure (I/O
         % error condition).  Nothing to do here but to report the issue
         % back to our caller and hope they have a way of dealing with it.

         errmsg = ferror(fid);
         fprintf(['I/O Error Reading Keyword Vector ', ...
                  '''%s'' (''%s'': %d): %s\n'], field, fname, pos, errmsg);
      end
   end
end
