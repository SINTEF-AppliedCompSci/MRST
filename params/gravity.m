function g = gravity(varargin)
%Manage effects of gravity.
%
% SYNOPSIS:
%   Either of the modes
%     1)     gravity()
%     2)     gravity(arg)
%     3)     gravity arg1 arg2 ...
%     4) g = gravity( )
%
% PARAMETERS:
%     v - Control mode for effect of gravity.  OPTIONAL. Must be one of
%
%           - String, `{'off', 'on'}`, for disabling or enabling effects
%             of gravity.  The default state is 'off'.
%
%           - String `{'x', 'y', 'z'}`.  The string '<p>' sets the gravity
%             direction to point along the positive physical coordinate
%             direction '<p>' of the underlying grid model.
%
%             The default direction in function 'gravity' is 'z'.  In other
%             words, the gravity by default field points along the positive
%             Z axis of the underlying grid model.
%
%           - String, `'reset'`, for restoring effects of gravity to the
%             default state: gravity off, but acceleration strength
%             nevertheless as at the Tellus equator, in the direction
%             'z'.  Note that the gravity vector returned in call mode 2)
%             will be the zero vector unless effects of gravity are
%             specifically enabled.
%
%           - Logical, `{ false, true }`, for disabling or enabling effects
%             of gravity.
%
%           - Numeric (Real) scalar value.  Specifically set acceleration
%             strength.  The value is assumed to be in units of m/s^2.
%             However, unless effects of gravity are specifically enabled,
%             clients (those using call mode 2) will still retrieve a zero
%             gravity vector.
%
%           - Numeric (Real) two- or three-component real vector giving
%             explicit gravity vector.  In this case, the acceleration
%             strength is NORM(v) in units of m/s^2.
%
%
% RETURNS:
%     g - A three-component *ROW* vector defining the current gravity field.
%         The Euclidian norm of `g` is the acceleration strength, in units
%         of m/s^2, of the gravity field.
%
% NOTE:
%   String arguments may be combined together.  Later control modes
%   override former. The default state, unless changed using a call from
%   `gravity on`, is to disable all effects of gravity.  In other words,
%   `gravity off` is the default state of gravity effects.
%
% EXAMPLES:
%   % 1) Demonstrate 'String' form (Command and Function syntax).
%   gravity x on   % Enable gravity (of default strength).
%   gravity('off') % Disable gravity (default state).
%   gravity reset  % Restore gravity defaults (off).
%   gravity x      % Set gravity direction along x-axis of current model.
%
%   % 2) Demonstrate 'Logical' form of gravity function.
%   %    (only Function syntax supported in this case).
%   gravity(true)  % Enable gravity.
%   gravity(false) % Disable gravity (default state).
%
%   % 3) Demonstrate 'Numeric' form of gravity function.
%   %    (only Function syntax supported in this case).
%   gravity(1.0)   % Set acceleration strength to 1 m/s^2.
%   gravity(0.0)   % Set acceleration strength to 0 m/s^2 (i.e., disable).
%   % Gravity field along direction [1,1] of underlying grid model.
%   % Acceleration strength of SQRT(2) m/s^2.
%   gravity([1, 1])
%   % Set gravity field along direction [1,1,-10] of underlying grid model.
%   % Acceleration strength of SQRT(102) m/s^2.
%   gravity([1; 1; -10])
%
%   % 4) Retrieve current gravity vector field.
%   %    (Only function syntax supported in this case).
%   g = gravity()

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


   persistent G g_vec gravityOn

   if isempty(G)
      [G, g_vec, gravityOn] = reset();
   end

   for k = 1 : nargin
      arg = varargin{k};
      if ischar(arg)
         switch lower(arg)
            case {'off', 'no' , 'false'}, gravityOn = false;
            case {'on' , 'yes', 'true' }, gravityOn = true;
            case 'x'                    , g_vec     = [1, 0, 0];
            case 'y'                    , g_vec     = [0, 1, 0];
            case 'z'                    , g_vec     = [0, 0, 1];
            case 'reset'
               [G, g_vec, gravityOn] = reset();
             otherwise
               error(msgid('Mode:Unsupported'), ...
                     'Unsupported gravity control mode ''%s''.', arg);
         end
      elseif islogical(arg)
         gravityOn = arg;
      elseif isnumeric(arg) && isreal(arg)
         if numel(arg) == 1 && ~(norm(g_vec) > 0)
            error(msgid('GravDir:AllZero'), ...
                 ['Gravity does not point anywhere ',       ...
                  '(direction is [0,0,0]).\n',              ...
                  'Use ''gravity reset'' before changing ', ...
                  'the acceleration strength.']);
         end

         switch numel(arg)
            case 1, G = arg;  g_vec = g_vec ./ norm(g_vec)    ;
            case 2, G = 1  ;  g_vec = [reshape(arg, 1, []), 0];
            case 3, G = 1  ;  g_vec =  reshape(arg, 1, [])    ;

            case num2cell(4 : 26)
               error(['%d-D String Theory is yet to be proven as a ', ...
                      'valid description of our world.'], numel(arg));

             otherwise
               error(msgid('foo:bar'), ...
                    ['This does not compute.\n',                ...
                     'Our three dimensional, visible physical', ...
                     ' world cannot possibly\n',                ...
                     'be reconciled with a %d-component',       ...
                     ' gravitational field.\n\n',               ...
                     'Don''t let me catch you doing this again...'], ...
                     numel(arg));
         end
      else
         error(msgid('Mode:Unsupported'), ...
               'Gravity control does not conform to any supported mode.');
      end
   end

   if nargout || nargin == 0
       if gravityOn
          g = G * g_vec;
       else
          g = zeros([1, 3]);
       end
   end
end

%--------------------------------------------------------------------------

function [G, g_vec, gravityOn] = reset()
   G         = 9.80665;
   g_vec     = [0, 0, 1];
   gravityOn = false;
end
