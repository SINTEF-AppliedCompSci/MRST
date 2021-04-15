function varargout = mrstTranslateExitCode(varargin)
%Attempt to Translate a System Exit Code to Readable Message
%
% SYNOPSIS:
%   msg = mrstTranslateExitCode(ecode)
%
% PARAMETERS:
%   ecode - System exit code.  Typically the first return value from the
%           SYSTEM function, and especially if that value is non-zero.
%
% RETURNS:
%   msg - Human readable string that (hopefully) contains clues as to what
%         caused a failure.  Implemented in terms of the C++
%         'std::system_category()' facility.
%
% SEE ALSO:
%   system.

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

   if ~has_cxx_compiler()
      varargout{1} = 'Untranslatable (no C++ MEX support)';
      return;
   end

   INCLUDE = {};

   OPTS = { '-O' };
   SRC  = { 'mrstTranslateExitCode.cpp' };

   [CXXFLAGS, LINK, LIBS] = mrstDefaultMexFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});

   [varargout{1:nargout}] = mrstTranslateExitCode(varargin{:});
end

%--------------------------------------------------------------------------

function tf = has_cxx_compiler()
   tf = ~isempty(mex.getCompilerConfigurations('c++'));
end
