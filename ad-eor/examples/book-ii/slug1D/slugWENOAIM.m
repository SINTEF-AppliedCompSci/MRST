%% Slug injection case - discretization schemes
% This script runs the last part of the grid refinement study discussed in
% the 1D slug injection example from the EOR chapter.
%
% If you want to also plot the fully refined solution, it assumes that you
% have conducted simulations with the setup from the SP_1D_hires.DATA input
% file and transformed the results as follows:
% x  = model.G.cells.centroids(:,1);
% t  = reports.ReservoirTime/day;
% s  = cellfun(@(x) x.s(:,1),  states, 'UniformOutput',false); S =horzcat(S{:});
% cp = cellfun(@(x) x.cp(:,1), states, 'UniformOutput',false); Cp=horzcat(Cp{:});
% cs = cellfun(@(x) x.cs(:,1), states, 'UniformOutput',false); Cs=horzcat(Cs{:});
% save slug1D-2slugs.mat x t s cp cs ;
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui

%% Load detailed simulated data
try
    load slug1D-2slugs.mat;
    refsol = true;
catch
    refsol = false;
end
cols = lines(6);
fig=figure('Position',[450 360 800 420]);

%% Set up and run simulations
bookdir = getDatasetPath('eor_book_ii');
nx = [25, 50, 100, 200];
for i=1:numel(nx)
    figure(fig), leg = {};
    subplot(2,2,i)
    if refsol 
        plot(x,cs(:,4000),'Color',cols(5,:),'LineWidth',2);
        leg{end+1} = 'Reference';
    end
    hold on
    
    % Run SPU-FIM
    fn = fullfile(bookdir,'slug1D',sprintf('SP_1D_%d.DATA',nx(i)));
    [state0, model, schedule, nonlinear] = initEclipseProblemAD(fn);
    try
        [~, states] = ...
            simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
        y  = model.G.cells.centroids(:,1);
        figure(fig); plot(model.G.cells.centroids(:,1), states{end}.cs(:,1),...
            'Color',cols(1,:),'LineWidth',1.5);
        drawnow
        leg{end+1} = 'SPU-FIM';
    catch
    end

    % Run SPU-AIM
    model = setTimeDiscretization(model, 'aim', 'saturationCFL', 0.75);
    try
        [~,states] = ...
            simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
        figure(fig); plot(model.G.cells.centroids(:,1), states{end}.cs(:,1),...
            'Color',cols(2,:),'LineWidth',1.5);
        drawnow
        leg{end+1} = 'SPU-AIM';
    catch
    end

    % Run WENO-FIM
    model = setWENODiscretization(model);
    model = setTimeDiscretization(model, 'fim');
    try
        [~, states] = ...
            simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
        figure(fig); plot(model.G.cells.centroids(:,1), states{end}.cs(:,1),...
            'Color',cols(3,:),'LineWidth',1.5);
        drawnow
        leg{end+1} = 'WENO-FIM'; %#ok<*SAGROW>
    catch
    end
    
    % Run WENO-AIM
    model = setTimeDiscretization(model, 'aim', 'saturationCFL', 0.75);
    try
        [~,states] = ...
            simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
        figure(fig); plot(model.G.cells.centroids(:,1), states{end}.cs(:,1),...
            'Color',cols(4,:),'LineWidth',1.5);
        drawnow
        leg{end+1} = 'WENO-AIM';
    catch
    end
    
    legend(leg{:});
    axis([0 40 0 50]);
    text(2,45,['n=',num2str(nx(i))]);
end
               
%% Copyright notice
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
