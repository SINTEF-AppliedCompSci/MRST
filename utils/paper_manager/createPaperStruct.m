function paper = createPaperStruct(id, title, varargin)
%Create a structure containing information about a document that uses MRST
%
% SYNOPSIS:
%   paper = createPaperStruct(id, title)
%   paper = createPaperStruct(id, title, 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   id    - A short string ID for the paper that can be used to
%           programmatically refer to the same id over multiple revisions.
%
%   title - A string containing the title of the paper.
%
% OPTIONAL PARAMETERS:
%   'authors'   - A string containing the names of the authors.
%
%   'published' - Publication avenue (name of conference, journal name with
%                 issue number, thesis, ...)
%
%   'url'       - URL to the official site of a published paper. Typically,
%                 this is a webpage on the publisher's website where a the
%                 paper can be viewed.
%
%   'fileurl'   - URL to a direct download of a preprint of the paper, if
%                 available from the copyright holders.
%
%   'year'      - Double indicating the publication year.
%
%   'modules'   - Cell array of the names modules where the paper is
%                 relevant as documentation or background information.
%
%   'doi'       - Digital object identifier for the publication.
%
% RETURNS:
%   paper - Structure with defaulted values for keywords not specified.
%
% SEE ALSO:
%   `getAvailablePapers`, `mrstReferencesGUI`

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
    paper = struct('authors', '',...
                   'published', '',...
                   'url', '',...
                   'year', -1, ...
                   'doi',   '', ...
                   'modules', {{}}, ...
                   'fileurl', '');
    paper = merge_options(paper, varargin{:});
    paper.id = id;
    paper.name = title;
end
