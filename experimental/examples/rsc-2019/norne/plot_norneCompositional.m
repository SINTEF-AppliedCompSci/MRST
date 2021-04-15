[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems, 'readFromDisk', false);

%%

ws = cellfun(@(w) {w{1:w.numelData}}, ws, 'unif', false);

%%

loadStates = cell(size(states));
for mNo = 1:2
    st = cell(states{mNo}.numelData,1);
    for sNo = 1:states{mNo}.numelData
        fprintf('Loading state %d of %d ... \n', sNo, states{mNo}.numelData);
        s = states{mNo}{sNo};
        s.transportModel = [];
        st{sNo} = s;
    end
    loadStates{mNo} = st;
end

%%

sd = cellfun(@(s1, s2) compareStates(s1, s2), loadStates{1}(1:19), loadStates{2}, 'unif', false);

%%

plotWellSols(ws)

%% Copyright Notice
%
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
