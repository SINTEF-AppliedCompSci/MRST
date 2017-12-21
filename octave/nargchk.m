## Copyright (C) 2008, 2009 Bill Denney
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{msgstr} =} nargchk (@var{minargs}, @var{maxargs}, @var{nargs})
## @deftypefnx {Function File} {@var{msgstr} =} nargchk (@var{minargs}, @var{maxargs}, @var{nargs}, "string")
## @deftypefnx {Function File} {@var{msgstruct} =} nargchk (@var{minargs}, @var{maxargs}, @var{nargs}, "struct")
## Return an appropriate error message string (or structure) if the
## number of inputs requested is invalid.
##
## This is useful for checking to see that the number of input arguments
## supplied to a function is within an acceptable range.
## @seealso{nargoutchk, error, nargin, nargout}
## @end deftypefn

## Author: Bill Denney <bill@denney.ws>
##modification to not give error with struct

function msg = nargchk (mina, maxa, narg, outtype)
  if (nargin < 3 || nargin > 4)
    print_usage ();
  elseif (mina > maxa)
    error ("nargchk: minargs must be <= maxargs");
  elseif (nargin == 3)
    outtype = "string";
  elseif (! any (strcmpi (outtype, {"string", "struct"})))
    warning ("nargchk: output type must be either string or struct");
  elseif (! (isscalar (mina) && isscalar (maxa) && isscalar (narg)))
    error ("nargchk: mina, maxa, and narg must be scalars");
  endif

  msg = struct ("message", "", "identifier", "");
  if (narg < mina)
    msg.message = "not enough input arguments";
    msg.identifier = "Octave:nargchk:not-enough-inputs";
  elseif (narg > maxa)
    msg.message = "too many input arguments";
    msg.identifier = "Octave:nargchk:too-many-inputs";
  endif

  if (strcmpi (outtype, "string"))
    msg = msg.message;
  elseif (isempty (msg.message))
    %msg = struct ([]);
    msg = [];
  endif
endfunction

