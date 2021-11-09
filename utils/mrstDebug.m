function debug = mrstDebug(varargin)
%Globally control default settings for MRST verbose information.
%
% SYNOPSIS:
%   Either of the modes
%     1) debuglevel = mrstDebug(debuglevel)
%     2) v = mrstDebug
%
% DESCRIPTION:
%   This function provides a centralised facility for individual MRST
%   functions to query verbosity default settings.  Specifically, if a
%   function provides verbose output (i.e., additional reporting during
%   computational process) in the form of a 'debug' option, then the
%   function is encouraged to initialise the verbose flag with the output
%   of this function. Debug output may then be subsequently overridden on
%   an individual, per-function basis.
%
% PARAMETERS:
%   arg     - Numeric (Real) scalar value.  Specifically set debug
%             level.  A verbosity level exceeding zero enables verbose
%             output. Whether or not higher values produce more debug options is
%             at the discretion of individual MRST functions.
%
%
% RETURNS:
%     v - A numeric scalar value. 
%
% EXAMPLES:
%   mrstDebug(1)   % Enable verbose output.
%   mrstDebug(0)   % Disable verbose output.
%   mrstDebug(2)   % Set verbose output level to 2.
%
%   % Retrieve current verbosity setting.
%   %    (Only function syntax supported in this case).
%   v = mrstDebug()

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


   persistent DEBUG;

   if isempty(DEBUG),
      DEBUG = 0;
   end

   if nargin > 0,
      arg = varargin{1};
      assert(isnumeric(arg) && isreal(arg))
      DEBUG = arg;
   end
   debug = DEBUG;
end
