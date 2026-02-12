function T = makehgtform(varargin)
%Create 4x4 transformation matrix for graphics
%
% SYNOPSIS:
%   T = makehgtform('xrotate', angle)
%   T = makehgtform('yrotate', angle)
%   T = makehgtform('zrotate', angle)
%
% DESCRIPTION:
%   Octave-compatible implementation of makehgtform for basic rotations.
%   Creates a 4x4 homogeneous transformation matrix for graphics objects.
%
% PARAMETERS:
%   'xrotate', angle - Rotation around X-axis by angle (in radians)
%   'yrotate', angle - Rotation around Y-axis by angle (in radians)
%   'zrotate', angle - Rotation around Z-axis by angle (in radians)
%
% RETURNS:
%   T - 4x4 homogeneous transformation matrix

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

   % Initialize as identity matrix
   T = eye(4);
   
   % Process arguments
   i = 1;
   while i <= length(varargin)
      arg = varargin{i};
      
      if ischar(arg)
         switch lower(arg)
            case 'xrotate'
               if i+1 <= length(varargin)
                  angle = varargin{i+1};
                  T = T * xrotation_matrix(angle);
                  i = i + 2;
               else
                  error('makehgtform:MissingAngle', ...
                        'xrotate requires an angle argument');
               end
               
            case 'yrotate'
               if i+1 <= length(varargin)
                  angle = varargin{i+1};
                  T = T * yrotation_matrix(angle);
                  i = i + 2;
               else
                  error('makehgtform:MissingAngle', ...
                        'yrotate requires an angle argument');
               end
               
            case 'zrotate'
               if i+1 <= length(varargin)
                  angle = varargin{i+1};
                  T = T * zrotation_matrix(angle);
                  i = i + 2;
               else
                  error('makehgtform:MissingAngle', ...
                        'zrotate requires an angle argument');
               end
               
            otherwise
               error('makehgtform:UnsupportedTransform', ...
                     'Transformation ''%s'' not implemented in Octave version', arg);
         end
      else
         error('makehgtform:InvalidArgument', ...
               'Expected string transformation type');
      end
   end
end

function T = xrotation_matrix(angle)
   % Rotation around X-axis
   c = cos(angle);
   s = sin(angle);
   T = [1,  0,  0, 0;
        0,  c, -s, 0;
        0,  s,  c, 0;
        0,  0,  0, 1];
end

function T = yrotation_matrix(angle)
   % Rotation around Y-axis
   c = cos(angle);
   s = sin(angle);
   T = [ c, 0,  s, 0;
         0, 1,  0, 0;
        -s, 0,  c, 0;
         0, 0,  0, 1];
end

function T = zrotation_matrix(angle)
   % Rotation around Z-axis
   c = cos(angle);
   s = sin(angle);
   T = [c, -s,  0, 0;
        s,  c,  0, 0;
        0,  0,  1, 0;
        0,  0,  0, 1];
end
