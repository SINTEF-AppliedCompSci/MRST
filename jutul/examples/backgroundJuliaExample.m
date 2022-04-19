%% Brief example demonstrating how to run the Daemon mode
% We can seamlessly run Jutul in the background to accelerate certain MRST
% simulations. This is done using a background process using this package:
% https://github.com/dmolina/DaemonMode.jl
% The Julia process in the background is persistent, and will compile the
% required codes on the first run after the daemon is started. You will
% need to restart the daemon if it is terminated (e.g. by rebooting or
% closing the terminal window).
mrstModule add ad-core ad-blackoil deckformat ad-props test-suite jutul
%% How to set up
% One time setup. Two options:
% 1. Add DaemonMode to your base Julia environment (easiest)
% 2. Create a separate environment with DaemonMode (advanced users).
%
% For option 1:
% Run julia and run 
%    ]add DaemonMode
% It should then be added to your base environment.
%
% Then, proceed with the following 
%
%    mkdir jutul-daemon
%    cd .\jutul-daemon\
%    julia
%    ]activate .
%    add DaemonMode # Required to run the Daemon, not needed if following (1)
%    dev Jutul      # Could be add once we get it into package manager
%    dev JutulDarcy
%    pwd() # will give you the absolute path of the folder
%
% If you are using dev-ed packages, you may have to go into this
% environment as noted above and run the following:
%     using Jutul, JutulDarcy
%
% This needs to be done once per reboot (it seems). You can restart the
% Daemon freely afterwards.
%
% Finally, the daemon can be spun up by calling
%  julia --project="X" --startup-file=no -e 'using Revise; using DaemonMode; serve()'
% where X is the absolute path of the jutul_daemon folder.
% julia --project="C:\\Users\\olavm\\jutul-daemon" --startup-file=no -e 'using Revise; using DaemonMode; serve()'
%
% If you were following path (2) you will need to specify the directory X of
% your project via the 'project' optional input argument to
% runJutulOnDaemon.
%% We run a small example to demonstrate and verify that it is working
setup = qfs_wo();
[state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);
[wells, states] = runJutulOnDaemon(state0, model, schedule, 'name', 'qfs_wo');
%% Plot the results
mrstModule add mrst-gui
figure;
plotToolbar(model.G, states);
