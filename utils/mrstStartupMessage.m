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
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    isDesktop = usejava('desktop');

    printHeader(isDesktop)

    fprintf('\nUseful commands for getting started:\n');
    fprintf(' - List all introductory examples:   ');
    printMatlabLink(isDesktop, 'mrstExamples()');
    
    fprintf(' - List all modules:                 ');
    printMatlabLink(isDesktop, 'mrstPath(''list'')');

    fprintf(' - Load modules using GUI:           ');
    printMatlabLink(isDesktop, 'mrstModule(''gui'')');

    fprintf(' - Explore all available data sets   ');
    printMatlabLink(isDesktop, 'mrstDatasetGUI()');

    fprintf(' - List examples of a module:        ');
    printMatlabLink(isDesktop, 'mrstExamples(''ad-blackoil'')');

    fprintf(' - Explore modules and publications: ');
    printMatlabLink(isDesktop, 'mrstExploreModules()');

    fprintf(' - Show all examples in all modules: ');
    printMatlabLink(isDesktop, 'mrstExamples(''all'')');

    fprintf(' - Display this message:             ');
    printMatlabLink(isDesktop, 'mrstStartupMessage()');
    
    fprintf(['\nFor assistance and discussions about MRST, ', ...
             'please visit our mailing list at\n']);
    fprintf('\t');
    printLink(isDesktop, 'www.sintef.no/projectweb/mrst/forum/', ...
              'http://www.sintef.no/projectweb/mrst/forum/');

    fprintf(' (');
    printLink(isDesktop, 'sintef-mrst@googlegroups.com', ...
              'mailto:sintef-mrst@googlegroups.com');
    fprintf(')\n');

    fprintf('For some common queries, see our FAQ: ')
    printLink(isDesktop, 'www.sintef.no/projectweb/mrst/faq/', ...
              'http://www.sintef.no/projectweb/mrst/faq/');
    fprintf('\n');
end

%--------------------------------------------------------------------------

function printHeader(isDesktop)
    fprintf(['Welcome to the Matlab Reservoir Simulation ', ...
             'Toolbox (MRST)!\nYou are using '])

    githead = fullfile(ROOTDIR, '.git', 'refs', 'heads', 'master');
    if exist(githead, 'file')
        % User is using the git-version.
        f = fopen(githead, 'r');
        fprintf('the development version at commit %s\n', fgetl(f));
        fclose(f);
    else
        % User is using a specific release. Give a bit of extra output.
        fprintf(['the release version 2020b. To download other ', ...
                 'versions of MRST\n', ...
                 'and view examples and relevant publications, ', ...
                 'please visit ']);

        printLink(isDesktop, 'www.mrst.no', 'www.mrst.no');
        fprintf('\n');
    end
end

%--------------------------------------------------------------------------

function printMatlabLink(isDesktop, text, endl)
    if nargin == 2
        endl = true;
    end
    printLink(isDesktop, text, ['matlab:', text])
    if endl
        fprintf('\n');
    end
end

%--------------------------------------------------------------------------

function printLink(isDesktop, text, link)
    if isDesktop
        fprintf('<a href = "%s">', link);
    end
    fprintf('%s', text);
    if isDesktop
        fprintf('</a>');
    end
end
