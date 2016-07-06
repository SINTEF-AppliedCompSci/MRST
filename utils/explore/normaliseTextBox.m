function s = normaliseTextBox(s)
%Remove explicit line breaks from string for presentation in GUI text box
%
% SYNOPSIS:
%   s = normaliseTextBox(s)
%
% PARAMETERS:
%   s - String or cell array of strings, possibly containing explicit line
%       breaks and/or words separated by multiple spaces.
%
% RETURNS:
%   s - String or cell array of strings (depending on input) without
%       explicit line breaks and with words separated by single space only.

%{
#COPYRIGHT#
%}

   s = removeCarriageReturn(s);

   s = regexprep(s, { '\n', '^\s+', '\s+', ' $' }, ...
                    { ' ' , ''    , ' '  , ''   });
end
