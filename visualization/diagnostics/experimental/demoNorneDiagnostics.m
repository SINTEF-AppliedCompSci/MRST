
prefix = '/data/norne/res/BC0407';
mrstModule add deckformat ad-fi mrst-gui diagnostics
cdims = [46   112    22];
tmp = readEclipseOutputFileUnFmt([prefix '.EGRID']);

grdecl.COORD = tmp.COORD.values;
grdecl.ZCORN = tmp.ZCORN.values;
grdecl.ACTNUM = tmp.ACTNUM.values;
grdecl.cartDims = cdims;

G = processGRDECL(grdecl, 'SplitDisconnected', false);
G = computeGeometry(G);
%%
clc
deck = readEclipseDeck([prefix, '_3.DATA']);
deck = convertDeckUnits(deck);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initEclipseFluid(deck);
%%
[smry, smspec] = readSummaryLocal(prefix);
[rstrt, rsspec] = readRestartLocal(prefix);



%%
% Run this twice, first time will crash because of field messyness
% tmp = struct([]);
for blargh = 1:2
    try
    badcells = [8897,11160,13423,15686,26019,28282,30545,32808,34944,37206,44347];
    ind = true(G.cells.num, 1);
    ind(badcells) = false;

    nt = numel(rstrt.PRESSURE);
    f = fields(rstrt);
    data = [];
    for i = 1:nt
        for j = 1:numel(f)
            fld = f{j};
            d = rstrt.(fld);
            if numel(d{i}) == G.cells.num
                tmp.(fld) = d{i}(ind);
            end
        end
        data = [data; tmp];
    %     clear tmp
    end
    catch
        continue
    end
end
%% Recreate grid with splitting...
rock.perm = rock.perm(ind, :);
rock.poro = rock.poro(ind);
rock.ntg = rock.ntg(ind);
G = processGRDECL(grdecl, 'SplitDisconnected', true);
G = computeGeometry(G(1));

%% Read wells. Remove some wells for clarity.
skipwells = {'B-4DH', 'K-1H', 'B-1H', 'B-1BH', 'D-1CH', 'B-2H'};
grdecl_file        = fullfile('~/simmatlab/branches/ioNorne/data_io', 'BC0407_IO.DATA');
deck2 = convertDeckUnits(readEclipseDeck(grdecl_file));
W_all = processWellsLocal(G, rock, deck2.SCHEDULE.control(1));

W = {};
for i = 1:numel(W_all)
    w = W_all(i);
%     if any(strcmpi(w.name, skipwells))
%         continue
%     end

    if ~strcmpi(w.type, 'bhp')
        w.type = 'rate';
    end
    W{end+1} = w;
end
W = [W{:}] .';

%% Save whole visualization history
states = [];
for i = 1:numel(data)
    sw = data(i).SWAT;
    sg = data(i).SGAS;
    states(i).s = min(max([1 - sw - sg, sg, sw], 0), 1);
    states(i).pressure = data(i).PRESSURE;
    states(i).perm = rock.perm;
    states(i).poro = rock.poro;
end

%% Set up state

sw = data(1).SWAT;
sg = data(1).SGAS;
so = 1 - sw - sg;

fluid.properties = @(x) fakeprops;

state = initResSol(G, data(1).PRESSURE, [sw, so, sg]);
%% Do diagnostics
mrstModule add deckformat ad-fi mrst-gui diagnostics
close all
interactiveDiagnostics(G, rock, W, state, 'state', state, 'fluid', fluid)
view(80, 70)
axis tight off

%% Plot states
close all;
plotToolbar(G, states);
axis tight off
