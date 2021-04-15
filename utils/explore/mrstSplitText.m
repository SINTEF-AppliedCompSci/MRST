function h = mrstSplitText(intro)
%Split example introduction (help text) into header and paragraphs
%
% SYNOPSIS:
%   h = mrstSplitText(intro)
%
% PARAMETERS:
%   intro - Introductory text of example, typically obtained through the
%           expression ::
%
%               help('ExampleName')
%
% RETURNS:
%   h - Cell array of strings.  The H1 line is the first element and the
%       remaining elements represent one paragraph each of the input text,
%       in order of appearance.
%
% NOTE:
%   This function is mainly intended to support the MRST module explorer
%   GUI, `mrstExploreModules`.
%
% SEE ALSO:
%   `help`, `mrstExploreModules`, `normaliseTextBox`.

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

   % Identify header (H1) line
   t = regexp(removeCarriageReturn(intro), '\n', 'once', 'split');

   % Split body into paragraphs and squash each paragraph to a single line
   % with words separated by single space.
   if numel(t)>1
       h = normaliseTextBox([ t{1}, regexp(t{2}, '\n\s*\n', 'split') ]);
   else
       h = normaliseTextBox(t{1});
   end
end
