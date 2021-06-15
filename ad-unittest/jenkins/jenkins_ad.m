cd('../../../')
ls
run('mrst-core/startup.m');

% Set up module directories
names = {'autodiff', ...
         'internal', ...
         'multiscale', ...
         'visualization', ...
         'model-io', ...
         'solvers', ...
         'thirdparty-modules'};
     
names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);
mrstPath('addroot', names{:});

names = {'test-datasets'};
names = cellfun(@(x) fullfile(ROOTDIR, '..', x), names, ...
                    'UniformOutput', false);
mrstPath('register', 'co2lab', fullfile(ROOTDIR, '..', 'mrst-co2lab', 'co2lab'));
mrstPath('register', 'agmg', fullfile('..', 'hnil-agmg'));
mrstPath('addroot', names{:});

% Ensure Matlab BGL is available
bgl_dir = fullfile(ROOTDIR(), 'utils', '3rdparty','matlab_bgl');
if ~exist(fullfile(bgl_dir, 'matlab_bgl'), 'dir')
    run(fullfile(bgl_dir,'downloadMBGL'))
end


% Remove previous tests
jfolder = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'XML');
if exist(jfolder, 'dir')
    rmdir(jfolder, 's')
end

tap_folder = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'JUnit');
if exist(tap_folder, 'dir')
    rmdir(tap_folder, 's')
end

mrstModule add ad-unittest ad-core ad-blackoil
runTestsAD('runIntegration',    true, ...
            'runUnit',          true, ...
            'runExamples',      true, ...
            'writeToDisk',      true, ...
            'writeXML',         true);
% Disable matlab
exit();

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
