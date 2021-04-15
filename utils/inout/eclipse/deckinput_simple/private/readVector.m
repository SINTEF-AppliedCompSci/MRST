function vec = readVector(fid, field, nel)
%Input vector of floating point numbers from ECLIPSE input file.
%
% SYNOPSIS:
%   v = readVector(fid, field, nel)
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

   if exist('OCTAVE_VERSION', 'builtin')
      % We're targeting Octave.  Use original (slow) rV impl.
      vec = readVectorOld(fid, field, nel);
   else
      % We're targeting M.  Use TEXTSCAN-based implementation.
      vec = readVector_textscan(fid, field, nel);
   end
end
