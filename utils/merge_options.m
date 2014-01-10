function prm = merge_options(prm, varargin)
%Override default control options.
%
% SYNOPSIS:
%   prm = merge_options(prm, 'pn1', pv1, ...)
%
% PARAMETERS:
%   prm       - Original/default control option structure.  The contents of
%               this structure is problem specific and defined by caller.
%
%   'pn1'/pv1 - 'key'/value pairs overriding default options in `prm'.  A
%               WARNING is issued if a given 'key' is not already present
%               in FIELDNAMES(prm).  Key names are case insensitive.
%
%               The ``message identifier'' of this warning is
%
%                   [<FUNCTIONNAME>, ':Option:Unsupported']
%
%               with <FUNCTIONNAME> being the name of the function calling
%               MERGE_OPTIONS or the string 'BASE' if MERGE_OPTIONS is
%               called directly from the base workspace.
%
% EXAMPLE:
%   prm = struct('foo', 1, 'bar', pi, 'baz', true)
%   prm = merge_options(prm, 'foo', 0, 'bar', rand(10), 'FimFoo', @exp)
%
% RETURNS:
%   prm - Modified parameter structure.
%
% SEE ALSO:
%   FIELDNAMES, WARNING, STRUCT.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


if nargin > 1,
   if mod(numel(varargin), 2) == 0 && ...
      all(iscellstr(varargin(1 : 2 : end))),
      st = dbstack(1);
      try
         caller = st(1).name;
      catch  %#ok
         caller = 'BASE';
      end
      ofn = fieldnames(prm);
      nfn = varargin(1 : 2 : end);
      nfv = varargin(2 : 2 : end);

      for n = 1 : numel(nfn),
         ix = find(strcmpi(nfn{n}, ofn));

         if ~isempty(ix),
            oclass = class(prm.(ofn{ix}));
            nclass = class(nfv{n});
            empty = isempty(prm.(ofn{ix})) || isempty(nfv{n});
            if empty || isequal(oclass, nclass),
               prm.(ofn{ix}) = nfv{n};
            else
               error([caller, ':Option:ValueWrongClass'], ...
                     ['Option ''', nfn{n}, ''' has value of ', ...
                      'unsupported type.\nExpected ''', oclass, ...
                      ''', but got ''', nclass, ''' in its place.']);
            end
         else
            warning([caller, ':Option:Unsupported'], ...
                    ['Option `', nfn{n}, ''' is not supported']);
         end
      end
   else
      error(msgid('Input:Huh'), ...
            'Huh? Did you remember to unpack VARARGIN?!?');
   end
end
end
