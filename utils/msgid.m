function s = msgid(s)
%Construct Error/Warning message ID by prepending function name.
%
% SYNOPSIS:
%   s = msgid(s)
%
% PARAMETERS:
%   s - A string of the form `[<component>]:<mnemonic>` suitable as an
%       ERROR or WARNING-type message identifier.
%
% RETURNS:
%   s - The same string, though with the name of the function calling
%       `msgid` prepended to `s` (or the string 'BASE' if function `msgid`
%       is called from the base workspace).
%
% SEE ALSO:
%   `error`, `warning`.

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


   st = dbstack(1);

   try
      caller = st(1).name;
   catch %#ok
      caller = 'BASE';
   end

   s = [caller, ':', s];
end
