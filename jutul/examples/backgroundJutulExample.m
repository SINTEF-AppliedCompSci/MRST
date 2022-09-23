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
%    ]add JutulDarcy
% DaemonMode and JutulDarcy will then be added to your base environment.
%
% You can then run the JutulDarcy daemon by executing the following
% command:
%  julia --startup-file=no --color=no -e 'using Revise; using DaemonMode; serve()'
%
% For option 2:
% Go to your favorite terminal and navigate to where you want to place your
% environment. Then, execute the following:
%    mkdir jutul-daemon
%    cd jutul-daemon # Linux command - use cd .\jutul-daemon\ if you are on Windows
%    julia --project="."
%    ]                    # To enter Pkg prompt
%    add DaemonMode       # Required to run the Daemon, not needed if following (1)
%    add Revise           # Recommended
%    add JutulDarcy       # Reservoir simulator itself
%    # hit backspace to exit the prompt
%    pwd() # will give you the absolute path of the folder
%
% If you are using dev-ed packages, you may have to go into this
% environment as noted above and run the following:
%     using JutulDarcy
%
% This needs to be done once per reboot (it seems). You can restart the
% Daemon freely afterwards.
%
% Finally, the daemon can be spun up by calling
%  julia --project="X" --startup-file=no --color=no -e 'using Revise; using DaemonMode; serve()'
% where X is the absolute path of the jutul_daemon folder.
%
% For example, the author's path on Windows looks like this:
% julia --project="C:\\Users\\olavm\\jutul-daemon" --startup-file=no --color=no -e 'using Revise; using DaemonMode; serve()'
%
% If you were following path (2) you will need to specify the directory X of
% your project via the 'project' optional input argument to
% runJutulOnDaemon.
%% We run a small example to demonstrate and verify that it is working
setup = qfs_wo();
[state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);

%% Run in background
[ws, states] = runJutulOnDaemon(state0, model, schedule, 'name', 'qfs_wo');
%% Plot the results
mrstModule add mrst-gui
figure;
plotToolbar(model.G, states);
%% We can also use the high level interface:
% Run on daemon (default)
% [ws, states] = simulateScheduleJutul(state0, model, schedule, 'daemon', true);
% Run manually
% [ws, states] = simulateScheduleJutul(state0, model, schedule, 'daemon', false);