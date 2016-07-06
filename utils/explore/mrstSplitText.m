function h = mrstSplitText(intro)
%Split example introduction (help text) into header and paragraphs
%
% SYNOPSIS:
%   h = mrstSplitText(intro)
%
% PARAMETERS:
%   intro - Introductory text of example, typically obtained through the
%           expression
%
%               help('ExampleName')
%
% NOTE:
%   This function is mainly intended to support the MRST module explorer
%   GUI, mrstExploreModules.
%
% RETURNS:
%   h - Cell array of strings.  The H1 line is the first element and the
%       remaining elements represent one paragraph each of the input text,
%       in order of appearance.
%
% SEE ALSO:
%   help, mrstExploreModules, normaliseTextBox.

%{
#COPYRIGHT#
%}

   % Identify header (H1) line
   t = regexp(removeCarriageReturn(intro), '\n', 'once', 'split');

   % Split body into paragraphs and squash each paragraph to a single line
   % with words separated by single space.
   h = normaliseTextBox([ t{1}, regexp(t{2}, '\n\s*\n', 'split') ]);
end
