function c = userConsent(question)
%userConsent Ask user for his/her consent to a Y/N question
%  c = userConsent(question) requests a text response from the user in the
%  form of a Y/N question. The function returns logical 1 (true) if the
%  the first letter in the user's response is a 'y' (or 'Y'), and returns
%  logical 0 otherwise. Default answer is 'Y'. 
%
%  Example:
%  userConsent('Do you want to continue') will produce the following prompt
%  "Do you want to continue? Y/N [Y]: "
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
answ = input([question '? Y/N [Y]: '], 's');
if isempty(answ)
    answ = 'Y';
end
c = strncmpi(answ, 'Y', 1);
end

