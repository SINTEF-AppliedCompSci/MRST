function papers = getAvailablePapers()
%Get structures for all papers known to MRST
%
% SYNOPSIS:
%   papers = getAvailablePapers
%
% DESCRIPTION:
%   Get structures (as defined by createPaperStruct) for all known papers.
%   Papers are known to MRST if they are found as functions named 'paper_*'
%   in the 'papers' subdirectory of the directory hosting function
%   getAvailablePapers.
%
% RETURNS:
%   papers - List of all papers known to MRST, as an array of structures.
%
% SEE ALSO:
%   `createPaperStruct`, `mrstExploreModules`, `mrstReferencesGUI`.

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

    this_dir = fileparts(mfilename('fullpath'));
    sets     = what(fullfile(this_dir, 'papers'));

    % Extract functions whose names start with 'paper_'.
    ispaper  = ~ cellfun(@isempty, regexp(sets.m, '^paper_', 'match'));

    % Loop over candidate functions, adding to array as we go (dynamic
    % expansion is okay since the number of papers is relatively small).
    papers = [];
    for func = reshape(sets.m(ispaper), 1, [])
        [paper, paper] = fileparts(func{1});                    %#ok<ASGLU>

        papers = [papers; feval(paper)];                        %#ok<AGROW>
    end
end
