function [output, outputLGR] = readEclipseOutputStream(fid, readField, varargin)
%Read single ECLIPSE output/result stream
%
% SYNOPSIS:
%   output = readEclipseOutputStream(fid, readField)
%   output = readEclipseOutputStream(fid, readField, 'pn1', pv1, ...)
%
% PARAMETERS:
%   fid      - Valid file identifier, as obtained from 'fopen', pointing to
%              a stream of summary or restart data.
%
%   readField -
%             Function handle offering the semantics of extracting a single
%             field (datum) from the summary data stream.  Assumed to
%             support the calling syntax
%
%                  [name, data] = readField(fid, ...)
%
%             where the 'fid' is the same open file stream as input to
%             function 'readEclipseSummary'.  The return value 'name' is
%             assumed to be a single string that names or identifies a
%             particular datum, while 'data' is the associated field data
%             represented as a structure with the fields '.type' and
%             '.values'.
%
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters.
%              This list is passed directly on to the field reader function
%              'readField'.
%
% RETURNS:
%   output - Data structure containing all field data represented in the
%            input stream 'fid'.  One structure field, obtained from
%            function 'readField', for each data field in the file.

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

   [output, outputLGR] = deal(struct());

   processLGRData = nargout > 1;
   lgrInfo        = updateLGRInfo;

   while ~feof(fid)
      [name, field] = readField(fid, varargin{:});

      if ~isempty(field)
         % '+' -> p, '-' -> 'n', and other non-word characters to '_'.
         %
         if ~isvarname(name)
            name = regexprep(name, {'+', '-'}, {'p', 'n'});
            name = genvarname(regexprep(name, '\W', '_'));
         end

         % keep track of lgr-fields
         if strcmp(name, 'LGR')
            lgrInfo = updateLGRInfo(lgrInfo, field);
         end

         % need to check if nnchead is for lgs
         if strcmp(name, 'NNCHEAD')
            gridno = field.values(2);
            if gridno > 0
               lgrInfo = updateLGRInfo(lgrInfo, gridno);
            end
         end

         if any(strcmp(name, {'ENDLGR', 'LGRSGONE'}))
            lgrInfo.inLGRSection = false;
         end

         if ~lgrInfo.inLGRSection
            if ~isfield(output, name), output.(name) = empty_data; end
            output.(name) = append_data(output.(name), field);

         elseif processLGRData
            if ~isfield(outputLGR, name)
               outputLGR.(name) = empty_data;
            end

            outputLGR.(name) = ...
               append_data(outputLGR.(name), field, lgrInfo.curLGRNum);
         end
      end
   end
end

%--------------------------------------------------------------------------

function e = empty_data
   e = struct('type', [], 'values', []);
end

%--------------------------------------------------------------------------

function q = append_data(q, d, n)
   if isempty(q.type), q.type = d.type; end

   assert (strcmp(q.type, d.type), ...
           'Data type of new field data does not match existing.');
   if nargin < 3
       q.values = [q.values; reshape(d.values, [], 1)];
   else
       if numel(q.values) < n, q.values{n} = []; end
       q.values{n} = [q.values{n}; reshape(d.values, [], 1)];
   end
end

%--------------------------------------------------------------------------

function lgrInfo = updateLGRInfo(lgrInfo, input)
   if nargin == 0
      lgrInfo = struct('curLGRNum',         0, ...
                       'inLGRSection',  false, ...
                       'names',          {{}});

   elseif isstruct(input) % input is lgr-field

      curname = input.values{1};
      curnum  = find(strcmp(curname, lgrInfo.names), 1);

      if isempty(curnum)
         curnum = numel(lgrInfo.names) + 1;
         lgrInfo.names{curnum} = curname;
      end

      lgrInfo.curLGRNum    = curnum;
      lgrInfo.inLGRSection = true;

   elseif isnumeric(input) % nnc-info for local grid
      lgrInfo.curLGRNum  = input(1);
      lgrInfo.inLGRSection = true;
   end
end
