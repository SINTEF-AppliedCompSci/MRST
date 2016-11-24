%% Example: Interactive Diagnostics Tool for the SAIGUP Model
% In this example, we will use the SAIGUP model as set up in the
% 'saigupWithWells' example to show how to launch an interactive session for
% doing flow diagnostics.
mrstModule add book diagnostics;

%% Set up model
% We start by re-running this case study to set up the simulation model and
% compute a flow field. Henceforth, we only need the geological model, the
% description of the wells, and the reservoir state. Hence, we clear all
% other variables and close all plots produced by the script.
saigupWithWells; 
close all;
clearvars -except G rock W state

% Get rid of LaTeX notation in well names
% In the original example, we used LaTeX syntax to specify well names of
% the form $I_1$, $P_1$, etc. Unfortunately, this does not play very well
% along with MATLAB's tools for  graphical interfaces, and we therefore
% post-process the well names to be on the form I1, P1, etc.
for i=1:numel(W)
    W(i).name = W(i).name([1 5]);
end

%% Launch flow diagnostics
interactiveDiagnostics(G, rock, W, 'state', state, 'computeFlux', false);
set(gcf,'position',[440 317 866 480]); axis normal
view(-85,73)

%% Flux allocation plots
% To produce the flux and well allocation figure in the book:
% - in the 'region selection' tab: set maxTOF to zero to remove plotting of
%   volumetric quantities 
% - in the 'advanced' tab: select 'show grid', 'well pairs', and 'plot all
%   wells'. This will produce the overview of the reservoir
% - use your pointer (mouse) to click on any of the wells to bring up the
%   allocation plots

%% Snapshots of swept regions
% In the 'advanced' tab, deselect 'show grid', 'well pairs', and 'plot all
% wells'. Then, in the 'region selection' tab: 
% - set 'selection' to 'intersection'
% - set 'display' to 'tracer selected injector'
% - select wells (mark I6 and P2 to P6, or I4 and P1 to P4)
% - type 'shading faceted' on the command line
% - rotation, zoom, and translate the plot to get an acceptable view