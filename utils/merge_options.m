function prm = merge_options(prm, varargin)
%Override default control options.
%
% SYNOPSIS:
%   prm = merge_options(prm, 'pn1', pv1, ...)
%
% PARAMETERS:
%   prm -
%      Original/default option structure.  The contents of this structure
%      is problem specific and defined by the caller.
%
%   'pn'/pv -
%      List of 'key'/value pairs overriding default options in 'prm'.
%
%      A warning is issued, and no assignment made, if a particular 'key'
%      is not already present in FIELDNAMES(prm).  The message identifier
%      of this warning is
%
%          [<FUNCTIONNAME>, ':Option:Unsupported']
%
%      with <FUNCTIONNAME> being the name of function MERGE_OPTIONS' caller
%      or the string 'BASE' if MERGE_OPTIONS is used directly from the base
%      workspace (i.e., the Command Window).
%
%      Function MERGE_OPTIONS will fail (and call ERROR) if the new value's
%      class is different from the class of the existing value.
%
% RETURNS:
%   prm - Modified parameter structure.
%
% SPECIAL CASE:
%   If the value of a field of the input parameters ('prm') is a CELL
%   array, then the overriding value of that field can be anything.  If the
%   new value is another CELL array (i.e., if ISCELL returns true) it will
%   simply be assigned.  Otherwise, we wrap the overriding value in a cell
%   array so that the field value is always a CELL array.
%
%   This behaviour allows the user of function MERGE_OPTIONS to implement
%   uniform support for both single elements and heterogeneous collections
%   of data in a single option.  That in turn is useful in, for instance, a
%   visualisation application.
%
% EXAMPLE:
%   % 1) Typical use
%   prm = struct('foo', 1, 'bar', pi, 'baz', true)
%   prm = merge_options(prm, 'foo', 0, 'bar', rand(10), 'FimFoo', @exp)
%
%   % 2) Heterogeneous collection in a CELL array
%   prm = struct('f', {{ (@(x) x.^2) }}) % 'f' is cell array of f-handles
%   prm = merge_options(prm, 'f', @exp)  % Pass a simple function handle
%   fplot(prm.f{1}, [0, 3])              % Reference cell array result
%
%   % 3) Heterogeneous collection in a CELL array
%   prm = struct('d', {{ rand(10) }})   % 'd' is cell array of data points
%
%   % Pass multiple data sets
%   prm = merge_options(prm, 'd', { ones([5, 1]), linspace(0, 1, 11) })
%
%   % Plot "last" data set
%   plot(prm.d{end}, '.-')
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
            iscellstr(varargin(1 : 2 : end)),
         st = dbstack(1);
         try
            caller = regexprep(st(1).name, '\W', '_');
         catch  %#ok
            caller = 'BASE';
         end

         prm = process_options(prm, caller, varargin{:});

      else

         error(msgid('Input:NotKeyValuePairs'), ...
              ['Input arguments do not appear to be a list of ', ...
               '''key''/value pairs.\nDid you unpack VARARGIN?']);
      end
   end
end

%--------------------------------------------------------------------------

function prm = process_options(prm, caller, varargin)
   ofn = fieldnames(prm);
   nfn = varargin(1 : 2 : end);
   nfv = varargin(2 : 2 : end);

   for n = 1 : numel(nfn),
      ix = find(strcmpi(nfn{n}, ofn));

      if ~isempty(ix),
         if iscell(prm.(ofn{ix})) && ~iscell(nfv{n}),
            % Original is CELL -> accept anything by turning "new"
            % into CELL too.
            nfv{n} = nfv(n);
         end

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
end
