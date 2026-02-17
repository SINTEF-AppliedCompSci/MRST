function varargout = camproj(varargin)
%Set or get camera projection mode (Octave stub)
%
% SYNOPSIS:
%   camproj('perspective')
%   camproj('orthographic')
%   mode = camproj()
%
% DESCRIPTION:
%   Octave-compatible stub for camproj. This function provides basic
%   compatibility but does not actually change the camera projection in
%   Octave, as this feature is not fully supported. The function silently
%   accepts the commands without error to allow scripts to run.
%
% NOTE:
%   In Octave, perspective vs orthographic projection may not be fully
%   supported or may behave differently than in MATLAB. This stub allows
%   code to run without errors, but visual output may differ.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

   % Silently accept the command
   % In Octave, camera projection mode is not fully implemented
   % This stub allows scripts to run without errors
   
   if nargout > 0
      % If output is requested, return 'orthographic' as default
      varargout{1} = 'orthographic';
   end
end
