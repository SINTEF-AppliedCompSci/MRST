function verbose = mrstVerbose(varargin)
%Globally control default settings for MRST verbose information.
%
% SYNOPSIS:
%   Either of the modes
%     1)     mrstVerbose arg
%     2) v = mrstVerbose
%
% DESCRIPTION:
%   This function provides a centralised facility for individual MRST
%   functions to query verbosity default settings.  Specifically, if a
%   function provides verbose output (i.e., additional reporting during
%   computational process) in the form of a 'verbose' option, then the
%   function is encouraged to initialise the verbose flag with the output
%   of this function. Verbose output may then be subsequently overridden on
%   an individual, per-function basis.
%
% PARAMETERS:
%   arg - (Mode 1 only) Control mode for verbose output.  Must be one of
%         the following possibilities:
%           - String, `{'off', 'on'}`, for globally disabling or enabling
%             MRST verbose output.  Actual effects depends on specific
%             setting in individual functions and may usually be controlled
%             more targetly.  The default state is 'off'.
%
%           - String, 'reset', for restoring verbose output setting to the
%             default state: mrstVerbose off
%
%           - Logical, `{false, true}`, for disabling or enabling MRST
%             verbose output.
%
%           - Numeric (Real) scalar value.  Specifically set verbosity
%             level.  A verbosity level exceeding zero enables verbose
%             output. Whether or not higher values produce more output is
%             at the discretion of individual MRST functions.
%
%
% RETURNS:
%     v - (Mode 2 only) A numeric scalar value. `V==0` turns default state
%         of verbose output off, while `v>0` signifies different levels of
%         verbose output.  Individual callers of 'mrstVerbose' (typically
%         other MRST functions) must support a Boolean on/off state, but
%         may, optionally, support a notion of verbosity levels where
%         higher levels signify more extensive output.
%
% EXAMPLES:
%   % 1) Demonstrate 'String' form (Command and Function syntax).
%   mrstVerbose  on    % Enable verbose output.
%   mrstVerbose('off') % Disable verbose ouput (default state).
%   mrstVerbose reset  % Restore verbosity defaults (off).
%
%   % 2) Demonstrate 'Logical' form of mrstVerbose function.
%   %    (only Function syntax supported in this case).
%   mrstVerbose(true)  % Enable verbose output.
%   mrstVerbose(false) % Disable verbose output (default state).
%
%   % 3) Demonstrate 'Numeric' form of mrstVerbose function.
%   %    (only Function syntax supported in this case).
%   mrstVerbose(1)   % Enable verbose output.
%   mrstVerbose(0)   % Disable verbose output.
%   mrstVerbose(2)   % Set verbose output level to 2.
%
%   % 4) Retrieve current verbosity setting.
%   %    (Only function syntax supported in this case).
%   v = mrstVerbose()

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


   persistent VERBOSE;

   if isempty(VERBOSE),
      VERBOSE = false;
   end

   if nargin > 0,
      arg = varargin{1};
      if ischar(arg),
         switch lower(arg),
            case {'on' , 'yes', 'true' }         , VERBOSE = true;
            case {'off', 'no' , 'false', 'reset'}, VERBOSE = false;
            otherwise,
               error(msgid('Mode:Unsupported'), ...
                     'Unsupported verbosity control mode ''%s''.', arg);
         end
      elseif islogical(arg) || (isnumeric(arg) && isreal(arg)),
         VERBOSE = arg;
      end
   else
      verbose = VERBOSE;
   end
end
