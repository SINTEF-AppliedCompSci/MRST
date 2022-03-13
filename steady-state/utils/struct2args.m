function c = struct2args(s, varargin)
%Converts a structure to a cell array with both fieldnames and values,
%which may be used as arguments to an MRST function.
%
% SYNOPSIS:
%   c = struct2cellWithNames(s)
% 
% PARAMETERS:
%   s - Single structure.
% 
% OPTIONAL PARAMETERS:
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters. The
%              supported options are:
%   include  - Only include the field names given in this cell array.
%   remove   - Remove the field names given in this cell array.
%              Note: remove will be ignored if include is given.
% 
% RETURNS:
%   c - Cell array of the form {'field1', value1, 'field2', value2, .... }
%       where the 'field' strings are the fieldnames of s. The following
%       entry in c is the corresponding value of that field in s.
% 
% COMMENTS:
%   If s contains N fields, then c will contain 2*N elements.
%   
%   The retured cell array c may be used to create a comma seperated list 
%   by calling c{:}.
% 
% See also STRUCT, CELL, STRUCT2CELL.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

% Options
opt = struct('remove',  [], ... % cell array of fields
             'include', []);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.include)
   assert(iscell(opt.include));
   for i = 1:numel(opt.include)
      ns.(opt.include{i}) = s.(opt.include{i});
   end
   s = ns;
elseif ~isempty(opt.remove)
   assert(iscell(opt.remove));
   s = rmfield(s, opt.remove);
end

if isempty(s)
   c = cell(0);
else
   assert(numel(s) == 1, 'Structure cannot be an array.');
   f = fieldnames(s);
   c = cell(1, 2*numel(f));
   c(1:2:end-1) = f;
   c(2:2:end)   = struct2cell(s);
end

end
