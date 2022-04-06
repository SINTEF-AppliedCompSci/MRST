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
% One time setup:
% mkdir jutul_daemon
% cd .\julia-daemon\
% julia
% ]activate .
% add DaemonMode # Require to run the Daemon
% dev Jutul      # Could be add once we get it into package manager
% dev JutulDarcy
% pwd() # will give you the absolute path of the folder

% Then the daemon can be spun up by calling
%  julia --project="X" --startup-file=no -e 'using Revise; using DaemonMode; serve()'
% where X is the absolute path of the jutul_daemon folder.
% julia --project="C:\\Users\\olavm\\julia-daemon" --startup-file=no -e 'using Revise; using DaemonMode; serve()'
%% We run a small example to demonstrate and verify that it is working
setup = qfs_wo();
[state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);
[wells, states] = runJutulOnDaemon(state0, model, schedule, 'qfs_wo');
%% Plot the results
mrstModule add mrst-gui
figure;
plotToolbar(model.G, states);
