% Function for the co2 assesment of the 3D leakage problem.
%
%{
Copyright 2021, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module. 

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

function [state,lr,tt,Q1]=co2_3Dfls(statesMICP,G,fluid)

    % Setup schedule
    dt = .5;
    nt = 120+59+23+24*(100-1); % 100 days
    clear schedule
    timesteps = repmat(dt, nt, 1);
    
    % Compute current porosity and permeability
    poro = 0.15-statesMICP.c-statesMICP.b;
    KK = fluid.K(poro);
    rock = makeRock(G, KK, poro);

    % Gravity
    gravity on

    % Create Well
    L=max(G.nodes.coords(:,1));
    Q1 = 1600/day; %m^3/day
    [~,iw]= min(abs((G.cells.centroids(:,1)+100).^2+ ...
                                               G.cells.centroids(:,2).^2));
    cellsWell =  1:1:G.cells.num;
    cellsWell1 = cellsWell(abs(G.cells.centroids(:,1)- ...
    G.cells.centroids(iw,1))<.01 & abs(G.cells.centroids(:,2)- ...
                G.cells.centroids(iw,2))<.01 & G.cells.centroids(:,3)>130);
    W = addWellMICP([], G, rock, cellsWell1, 'Type', 'rate', 'Comp_i', ...
         [eps,1-eps], 'Val', Q1, 'Radius',.15*meter,'name', 'I','dir','z');

    % Create model
    model = CO2Model(G, rock, fluid);

    %Boundary Condition
    f = boundaryFaces(G);
    f = f(abs(G.faces.normals(f,1))>eps & (G.faces.centroids(f,1)<-L+2 ...
                                           | G.faces.centroids(f,1)>L-2 ));
    fp = G.faces.centroids(f,3)*fluid.rhoWS*norm(gravity);
    bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);

    % Well different rates and times
    N = 2; % Number of injection changes
    M = zeros(N,2); % Matrix where entries per row are:time and rate.
    M(1,1) = 121*dt; 
    M(1,2) = Q1;
    M(2,1) = (121+59)*dt; 
    M(2,2) = Q1;
    
    % Make schedule
    schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
    for i=1:N
        schedule.control(i+1)=schedule.control(i);
        schedule.control(i+1).W(1).val=M(i,2);
        schedule.control(i+1).W(1).sign=sign(M(i,2));
        schedule.step.control(M(i,1)/dt:end)=i+1;
        schedule.step.val(M(i,1)/dt:end)=minute;
    end
    schedule.step.val(M(2,1)/dt:end)=hour;

    % Initial state
    state0 = initState(G, W, G.cells.centroids(:,3)*fluid.rhoWS* ...
                                      norm(gravity), [1-eps, eps]);

    % Simulate
    [~, states] = simulateScheduleAD(state0, model, schedule);
    state=states{end};
    
    %Compute leakage
    cellsfa =  1:1:G.faces.num;
    cellsfac = cellsfa(G.faces.centroids(:,3)<80+1000*eps & ...
                                       G.faces.centroids(:,3)>80-1000*eps);
    for i=1:1:nt
        lr(i)=abs(states{i}.flux(cellsfac,2));
    end
    for i=1:1:120
        tt(1,i)=i*.5/day;
    end
    for i=121:1:120+59
        tt(1,i)=120*.5/day+(i-120)*minute/day;
    end
    for i=120+60:1:nt
        tt(1,i)=120*.5/day+59*minute/day+(i-120-59)*hour/day;
    end
end