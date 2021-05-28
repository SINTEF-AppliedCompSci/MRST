%% Slug injection case
% This script runs part of the grid refinement study discussed in the 1D
% slug injection example from the EOR chapter.
%
% To reproduce the figures exactly as they appear in the book chapter, you
% must first perform a somewhat time-consuming fine-scale simulation with
% the setup from the SP_1D_hires.DATA input file and transform the results
% as follows:
% x  = model.G.cells.centroids(:,1);
% t  = reports.ReservoirTime/day;
% s  = cellfun(@(x) x.s(:,1),  states, 'UniformOutput',false); S =horzcat(S{:});
% cp = cellfun(@(x) x.cp(:,1), states, 'UniformOutput',false); Cp=horzcat(Cp{:});
% cs = cellfun(@(x) x.cs(:,1), states, 'UniformOutput',false); Cs=horzcat(Cs{:});
% save slug1D-2slugs.mat x t s cp cs ;
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui

%% Load fine-scale simulations
cols = lines(6);
try
    load slug1D-2slugs.mat
    figure(1), clf
    subplot(3,1,1), hold on, plot(x,s(:,4000), 'Color',cols(1,:),'LineWidth',2);
    subplot(3,1,2), hold on, plot(x,cs(:,4000),'Color',cols(4,:),'LineWidth',2);
    subplot(3,1,3), hold on, plot(x,cp(:,4000),'Color',cols(3,:),'LineWidth',2);
    figure(2), clf
    plot3(s(:,4000),cs(:,4000), cp(:,4000),'LineWidth',2);
catch
end

%% Set up and run model
bookdir = getDatasetPath('eor_book_ii');
nx = [200, 100, 50, 25, 12];
lt = {'-','--','-.', ':','-o'};
qOvals = zeros(3,2);
for i=1:numel(nx)
    fn = fullfile(bookdir,'slug1D',sprintf('SP_1D_%d.DATA',nx(i)));
    [state0, model, schedule, nonlinear] = initEclipseProblemAD(fn);
    [ws1, states1, rep1] = ...
        simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
    qOs1 = cellfun(@(x) x(2).qOs, ws1,'UniformOutput',false);

    figure(1)
    x  = model.G.cells.centroids(:,1);
    subplot(3,1,1), hold on, plot(x, states1{end}.s (:,1),lt{i},'Color',cols(1,:),'LineWidth',1);
    subplot(3,1,2), hold on, plot(x, states1{end}.cs(:,1),lt{i},'Color',cols(4,:),'LineWidth',1);
    subplot(3,1,3), hold on, plot(x, states1{end}.cp(:,1),lt{i},'Color',cols(3,:),'LineWidth',1);
    drawnow
    
    figure(2), hold on
    plot3(states1{end}.s (:,1), states1{end}.cs(:,1), states1{end}.cp(:,1),...
        lt{i},'LineWidth',1.5,'Color',cols(1,:));
    hold off; view(3); box on; drawnow

    %{
    % Simulation for Strategy 1
    schedule.control(2).W(1).cp = 3;
    schedule.control(3).W(1).cs = 50;
    [ws2, states2, rep2] = ...
        simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
    qOs2 = cellfun(@(x) x(2).qOs, ws2,'UniformOutput',false);
    
    qOvals(i,1) = -sum(vertcat(qOs1{:}).*schedule.step.val);
    qOvals(i,2) = -sum(vertcat(qOs2{:}).*schedule.step.val);
    %}
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
