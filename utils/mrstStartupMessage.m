function mrstStartupMessage()
%Print a welcome message with helpful commands for new MRST users
%
% SYNOPSIS:
%   mrstStartupMessage
%
% DESCRIPTION:
%   Display the welcome message in the command window, indicating that MRST
%   is activated and ready for use. Some helpful links functions are also
%   provided to help new users getting started.
%
% EXAMPLES:
%    mrstStartupMessage()
%
% NOTE:
%    `mrstStartupMessage` is normally automatically run during the startup
%    process. Seeing the output from `mrstStartupMesssage` indicates that
%    MRST was successfully loaded and is ready for use.
%
% SEE ALSO:
%   `startup`, `mrstExamples`, `mrstModule`, `mrstPath`

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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

    fprintf('Welcome to the Matlab Reservoir Simulation Toolbox (MRST)!\nYou are using ')
    githead = fullfile(ROOTDIR, '.git', 'refs', 'heads', 'master');
    if exist(githead, 'file')
        % User is using the git-version.
        f = fopen(githead, 'r');
        fprintf('the development version at commit %s\n', fgetl(f));
        fclose(f);
    else
        % User is using a specific release. Give a bit of extra output.
        fprintf(['the release version 2018a. To download other versions of MRST\n',...
            'and view examples and relevant publications, please visit',...
            ' <a href="http://www.sintef.no/mrst">www.sintef.no/mrst</a>.\n']);
    end
    fprintf('\nUseful commands for getting started:\n');
    fprintf([' - List all introductory examples:   '...
        '<a href="matlab:mrstExamples()">mrstExamples()</a>\n']);
    fprintf([' - List all modules:                 ',...
        '<a href="matlab:mrstPath(''list'')">mrstPath(''list'')</a>\n']);
    fprintf([' - Load modules using GUI:           ',...
        '<a href="matlab:mrstModule(''gui'')">mrstModule(''gui'')</a>\n']);
    fprintf([' - Explore all available data sets   ',...
        '<a href="matlab:mrstDatasetGUI()">mrstDatasetGUI()</a>\n']);
    fprintf([' - List examples of a module:        ',...
        '<a href="matlab:mrstExamples(''ad-blackoil'')">mrstExamples(''ad-blackoil'')</a>\n']);
    fprintf([' - Show all examples in all modules: ',...
        '<a href="matlab:mrstExamples(''all'')">mrstExamples(''all'')</a>\n']);
    fprintf([' - Explore modules and publications: ',...
        '<a href="matlab:mrstExploreModules()">mrstExploreModules()</a>\n']);
    fprintf([' - Display this message:             ',...
        '<a href="matlab:mrstStartupMessage()">mrstStartupMessage()</a>\n']);
    fprintf(['\nFor assistance and discussions about MRST, please visit our mailing list at\n',...
        '\t<a href="http://www.sintef.no/projectweb/mrst/forum/">www.sintef.no/projectweb/mrst/forum/</a>'...
        ' (<a href="mailto:sintef-mrst@googlegroups.com">sintef-mrst@googlegroups.com</a>)\n']);
    fprintf(['For some common queries, see our FAQ: ',...
        '<a href="http://www.sintef.no/projectweb/mrst/faq/">www.sintef.no/projectweb/mrst/faq/</a>\n']);
end
