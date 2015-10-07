%% Example: Interactive Diagnostics Tool for the SAIGUP Model
% In this example, we will use the SAIGUP model as set up in the
% 'saigupWithWells' example to show how to lanch an interactive session for
% doing flow diagnostics.
mrstModule add diagnostics;

%% Set up model
saigupWithWells; 
close all;
clearvars -except G rock W state

% Get rid of LaTeX notation in well names
for i=1:numel(W)
    W(i).name = W(i).name([1 5]);
end

%% Launch flow diagnostics
interactiveDiagnostics(G, rock, W, 'state', state);
set(gcf,'position',[440 317 866 480]); axis normal
view(-80,36)