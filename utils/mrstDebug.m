function debug = mrstDebug(varargin)
%Globally control default settings for MRST debugging information.
%
% SYNOPSIS:
%   Either of the modes
%     1)     mrstDebug arg
%     2) d = mrstDebug
%
% DESCRIPTION:
%   This function provides a centralised facility for individual MRST
%   function to query debugging support default setting.  Specifically, if
%   a function provides debugging support (e.g., additional state checks)
%   in the form of a debugging option, then the function is encouraged to
%   initialise the debugging flag with the output of this function.
%   Debugging output may then be subsequently overridden on an individual,
%   per-function basis.
%
% PARAMETERS:
%     arg - Control mode for debugging support. OPTIONAL.  Must be one of
%             - String, `{'off', 'on'}`, for globally disabling or enabling
%               MRST debugging support.  Actual effects depends on specific
%               setting in individual functions and may usually be controlled
%               more targetly.  The default state is 'off'.
%
%             - String, 'reset', for restoring debugging setting to the
%               default state: mrstDebug off
%
%             - Logical, `{false, true}`, for disabling or enabling MRST
%               debugging support.
%
%             - Numeric (Real) scalar value.  Specifically set debugging
%               level.  A debugging level exceeding zero enables debugging.
%               Whether or not higher values produce more output is at the
%               discretion of individual MRST functions.
%
%
% RETURNS:
%     d - A numeric scalar value.  D==0 turns default state of debugging
%         support off, while d>0 signifies different levels of debugging
%         support.  Individual callers of 'mrstDebug' (typically other
%         MRST functions) must support a Boolean on/off state, but may,
%         optionally, support a notion of debugging levels where higher
%         levels signify more invasive (and expensive) checks.
%
% EXAMPLES:
%   % 1) Demonstrate 'String' form (Command and Function syntax).
%   mrstDebug  on    % Enable debugging support.
%   mrstDebug('off') % Disable debugging (default state).
%   mrstDebug reset  % Restore debugging defaults (off).
%
%   % 2) Demonstrate 'Logical' form of mrstDebug function.
%   %    (only Function syntax supported in this case).
%   mrstDebug(true)  % Enable debugging.
%   mrstDebug(false) % Disable debugging (default state).
%
%   % 3) Demonstrate 'Numeric' form of mrstDebug function.
%   %    (only Function syntax supported in this case).
%   mrstDebug(1)   % Enable debugging.
%   mrstDebug(0)   % Disable debugging.
%   mrstDebug(2)   % Set debugging level to 2.
%
%   % 4) Retrieve current debugging setting.
%   %    (Only function syntax supported in this case).
%   d = mrstDebug()

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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


   persistent DEBUG;

   if isempty(DEBUG),
      DEBUG = false;
   end

   if nargin > 0,
      arg = varargin{1};
      if ischar(arg),
         switch lower(arg),
            case {'on' , 'yes', 'true' }         , DEBUG = true;
            case {'off', 'no' , 'false', 'reset'}, DEBUG = false;
            otherwise,
               error(msgid('Mode:Unsupported'), ...
                     'Unsupported debug control mode ''%s''.', arg);
         end
      elseif islogical(arg) || (isnumeric(arg) && isreal(arg)),
         DEBUG = arg;
      end
   else
      debug = DEBUG;
   end
end
