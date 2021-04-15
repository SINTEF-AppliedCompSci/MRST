function [prm, varargout] = merge_options(prm, varargin)
%Override default control options.
%
% SYNOPSIS:
%    prm         = merge_options(prm, 'pn1', pv1, ...)
%   [prm, extra] = merge_options(prm, 'pn1', pv1, ...)
%
% PARAMETERS:
%   prm -
%      Original/default option structure.  The contents of this structure
%      is problem specific and defined by the caller.
%
%   'pn'/pv -
%      List of 'key'/value pairs overriding default options in 'prm'.
%
%      A warning is issued, and no assignment made, if a particular `key`
%      is not already present in `fieldnames(prm)`.  The message identifier
%      of this warning is
%
%          [<FUNCTIONNAME>, ':Option:Unsupported']
%
%      with <FUNCTIONNAME> being the name of function `merge_options`'
%      caller or the string 'BASE' if `merge_options` is used directly from
%      the base workspace (i.e., the Command Window).
%
%      Function `merge_options` will fail (and call ERROR) if the new
%      value's class is different from the class of the existing value.
%
%      In the interest of convenience for the typical case of using MRST
%      interactively from the Command Window, `merge_options` matches keys
%      (option names) using case insensitive search (i.e., using function
%      `strcmpi`).  If multiple option fields match a given name, such as 
%      in the case of several fields differing only by capitalisation, the
%      `merge_options` function resorts to exact and case sensitive string
%      matching (`strcmp`) to disambiguate options.
%
% RETURNS:
%   prm -
%      Modified parameter structure.
%
%   extra -
%      Cell array of 'key'/value pairs from the 'pn'/pv list that were not
%      matched by any option in the control structure 'prm'.  This allows
%      using function `merge_options` in an intermediate layer to define a
%      set of options and to pass other options unchanged on to lower-level
%      implementation functions--for instance to wrap a pressure solver in
%      a higher-level structure.
%
%      If the caller requests extra output be returned, then no diagnostic
%      message will be emitted for unsupported/undeclared option pairs in
%      the input list.
%
%      If there are no unsupported options in the input list then 'extra'
%      is an empty cell array.
%
% NOTE:
%   If the value of a field of the input parameters ('prm') is a `cell`
%   array, then the overriding value of that field can be anything.  If the
%   new value is another `cell` array (i.e., if `iscell` returns true) it
%   will simply be assigned.  Otherwise, we wrap the overriding value in a
%   cell array so that the field value is always a `cell` array.
%
%   This behaviour allows the user of function `merge_options` to implement
%   uniform support for both single elements and heterogeneous collections
%   of data in a single option.  That in turn is useful in, for instance, a
%   visualisation application.
%
% EXAMPLE:
%   % 1) Typical use
%   prm = struct('foo', 1, 'bar', pi, 'baz', true)
%   prm = merge_options(prm, 'foo', 0, 'bar', rand(10), 'FimFoo', @exp)
%
%   % 2) Heterogeneous collection in a `cell` array
%   prm = struct('f', {{ (@(x) x.^2) }}) % 'f' is cell array of f-handles
%   prm = merge_options(prm, 'f', @exp)  % Pass a simple function handle
%   fplot(prm.f{1}, [0, 3])              % Reference cell array result
%
%   % 3) Heterogeneous collection in a `cell` array
%   prm = struct('d', {{ rand(10) }})   % 'd' is cell array of data points
%
%   % Pass multiple data sets
%   prm = merge_options(prm, 'd', { ones([5, 1]), linspace(0, 1, 11) })
%
%   % Plot "last" data set
%   plot(prm.d{end}, '.-')
%
% SEE ALSO:
%   `fieldnames`, `warning`, `strcmpi`, `strcmp`.

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

   if nargin > 1
      % Caller passed additional input.  Verify that that extra input can
      % be interpreted as a list of 'key'/value pairs.

      if mod(numel(varargin), 2) == 0 && ...
            iscellstr(varargin(1 : 2 : end))
         % Additional input is structurally sound from the point of view of
         % being an apparent list of 'key'/value pairs.  Process options.

         st = dbstack(1);
         try
            caller = regexprep(st(1).name, '\W', '_');
         catch  %#ok
            caller = 'BASE';
         end

         [prm, extra] = process_options(prm, caller, varargin{:});

         if nargout == 2
            % Caller requested that unused options be returned.  Abide by
            % this request and don't emit any diagnostics.
            varargout{1} = extra;

         elseif ~isempty(extra)
            % Caller did not request that unused options be returned
            % unmodified but there were unused options in the input.
            % Report those unused options to the Command Window.

            pl = ''; if numel(extra(1:2:end)) > 1, pl = 's'; end

            warning([caller, ':Option:Unsupported'], ...
                    'Unsupported option%s in %s%s', pl, caller, ...
                    sprintf('\n * %s', extra{1 : 2 : end}));

         elseif nargout > 2
            error(msgid('NargOut:TooMany'), ...
                  'Too many output arguments: %d', nargout);
         end

      else
         % Additional input does not appear to be a list of 'key'/value
         % pairs.  The most common case of this happening is the caller
         % using the syntax
         %
         %   prm = merge_options(prm, varargin)
         %
         % rather than
         %
         %   prm = merge_options(prm, varargin{:})
         %
         % Alert the caller/user to that possibility.

         error(msgid('Input:NotKeyValuePairs'), ...
              ['Input arguments do not appear to be a list of ', ...
               '''key''/value pairs.\nDid you unpack VARARGIN?']);
      end

   elseif nargout == 2
      % Caller requested that unmatched options be returned, but there were
      % no option pairs in the input list (i.e., ISEMPTY(varargin)).
      %
      % Return an empty CELL array to honour requirements of interface lest
      % we fail with a diagnostic of the form
      %
      %   Output argument "varargout" (and maybe others) not assigned ...

      varargout{1} = {};

   end
end

%--------------------------------------------------------------------------

function [prm, extra] = process_options(prm, caller, varargin)
   ofn = fieldnames(prm);
   nfn = varargin(1 : 2 : end);
   nfv = varargin(2 : 2 : end);

   unused = false(size(nfv));

   for n = 1 : numel(nfn)
      ix = find(strcmpi(nfn{n}, ofn));

      if numel(ix) > 1
         % Case insensitive match hits more than one field name.  Defer to
         % case sensitive matching as arbiter/disambiguator.  If we then
         % don't match *any* option, proceed to next field name.  Note: Due
         % to guarantees implied by ofn = FIELDNAMES(prm) we either match
         % zero or one field name here.

         ix = find(strcmp(nfn{n}, ofn));
      end

      if ~isempty(ix)
         if iscell(prm.(ofn{ix})) && ~iscell(nfv{n})
            % Original is CELL -> accept anything by turning "new"
            % into CELL too.
            nfv{n} = nfv(n);
         end

         oclass = class(prm.(ofn{ix}));
         empty = isempty(prm.(ofn{ix})) || isempty(nfv{n});
         if empty || isa(nfv{n}, oclass)
            prm.(ofn{ix}) = nfv{n};
         else
            nclass = class(nfv{n});
            error([caller, ':Option:ValueWrongClass'], ...
                  ['Option ''', nfn{n}, ''' has value of ', ...
                   'unsupported type.\nExpected ''', oclass, ...
                   ''', but got ''', nclass, ''' in its place.']);
         end
      else
         unused(n) = true;
      end
   end

   extra = reshape([ nfn(unused) ; nfv(unused) ], 1, []);
end
