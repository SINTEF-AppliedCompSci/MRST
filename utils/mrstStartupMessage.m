function mrstStartupMessage()
%Print a welcome message with helpful commands for new MRST users
%
% SYNOPSIS:
%   mrstStartupMessage
%
% DESCRIPTION:
%   Display the welcome message in the command window, indicating that MRST
%   is activated and ready for use.  Some helpful links functions are also
%   provided to help new users getting started.
%
% NOTE:
%   `mrstStartupMessage` is normally run automatically during the startup
%   process.  Seeing the output from `mrstStartupMesssage` indicates that
%   MRST was successfully loaded and is ready for use.
%
% SEE ALSO:
%   `startup`, `mrstExamples`, `mrstModule`, `mrstPath`.

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

    isDesktop = usejava('desktop');

    printHeader(isDesktop)

    fprintf('\nUseful commands for getting started:\n');

    [menu, nchar] = getUsefulCommands();
    for row = menu .'
       printCommand(isDesktop, row{:}, nchar);
    end

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
    fprintf(['Welcome to the MATLAB Reservoir Simulation ', ...
             'Toolbox (MRST)!\nYou are using '])

    githead = fullfile(ROOTDIR, '.git', 'refs', 'heads', 'master');
    if exist(githead, 'file')
        % User is using the git-version.
        fid = fopen(githead, 'rt');

        if fid < 0
            commit = '<unknown>';
        else
            commit = fgetl(fid);
            fclose(fid);
        end

        fprintf('the development version at commit %s\n', commit);
    else
        % User is using a specific release. Give a bit of extra output.
        fprintf(['the release version 2021b. To download other ', ...
                 'versions of MRST\n', ...
                 'and view examples and relevant publications, ', ...
                 'please visit ']);

        printLink(isDesktop, 'www.mrst.no', 'www.mrst.no');
        fprintf('\n');
    end
end

%--------------------------------------------------------------------------

function [menu, nchar] = getUsefulCommands(padding)
   if nargin < 1, padding = 1; end

   menu = { ...
      'List all introductory examples',   'mrstExamples()'; ...
      'List all modules',                 'mrstPath(''list'')'; ...
      'Load modules using GUI',           'mrstModule(''gui'')'; ...
      'Explore all available data sets',  'mrstDatasetGUI()'; ...
      'List examples of a module',        'mrstExamples(''ad-blackoil'')'; ...
      'Explore modules and publications', 'mrstExploreModules()'; ...
      'Show all examples in all modules', 'mrstExamples(''all'')'; ...
      'Show settings for MRST',           'mrstSettings()'; ...
      'Display this message',             'mrstStartupMessage()'; ...
   };

   nchar = max(cellfun('prodofsize', menu(:, 1))) + padding;
end

%--------------------------------------------------------------------------

function printCommand(isDesktop, summary, command, nchar)
   fprintf(' - %-*s', nchar + 1, [summary, ':']); % +1 for the colon.
   printMatlabLink(isDesktop, command);
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
