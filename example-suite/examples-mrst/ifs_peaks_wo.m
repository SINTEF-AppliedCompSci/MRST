function [description, options, state0, model, schedule, plotOptions] = ifs_peaks_wo(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    % One-line description
    description ...
        = ['Inverted five-spot pattern with two-phase fluid and ', ...
           'perm/poro made from repeated pattern of peaks'       ];
    % Optional input arguments
    options = struct('ncells', 51, ... % Number of cells in xy for each tile
                     'tiles' , 1 , ... % Number of tiles using peaks (tiles x tiles)
                     'nkr'   , 2 );    % Brooks-Corey relperm exponent
    options = merge_options(options, varargin{:});
    % Adjust so that center well can be placed exactly in model center
    options.ncells = options.ncells - rem(options.ncells,2) + 1;
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil coarsegrid
    % We generate perm/poro by combining multiple tiles of peaks
    ncells = options.ncells;
    partition = partitionCartGrid([1,1]*ncells, [2,2]);
    m = floor(ncells/2)+1;
    partition(m:ncells:end)          = 5;
    partition(ncells*(m-1):ncells*m) = 5;
    k = zeros(ncells*ncells,1);
    for i = 1:4
        ix = find(partition == i);
        [xl, yl] = ind2sub([ncells,ncells], ix);
        XL = [xl, yl];
        XL = (XL - mean(XL))./(ncells/2)*6;
        if i == 2 || i == 4
            XL(:,1) = -XL(:,1);
        end
        if i > 2
            XL(:,2) = -XL(:,2);
        end
        k(ix) = peaks(XL(:,1), XL(:,2));
    end
    k = k./max(k);
    K = zeros(ncells+1);
    K(1:ncells,1:ncells) = reshape(k, ncells, ncells);
    K = repmat(K, options.tiles*[1,1]);
    K = reshape(K(1:end-1,1:end-1),[],1);
    N = ncells.*options.tiles+options.tiles-1;
    % Make model
    G = computeGeometry(cartGrid(N*[1,1], [1000, 1000]*meter));
    logperm = log10(100*milli*darcy) + K*3;
    poro    = 0.4*(1 + K*0.9);
    rock = makeRock(G, 10.^logperm, poro);
    fluid = initSimpleADIFluid('phases', 'WO'                    , ...
                               'n'     , [1,1].*options.nkr      , ...
                               'mu'    , [1,1]*centi*poise       , ...
                               'rho'   , [1,1]*kilogram/(meter^3));
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    % Wells
    time = 2*year;                        % Injection time
    rate = sum(poreVolume(G, rock))/time; % Injection rate
    bhp  = 50*barsa;                      % Producer bhp
    injector = @(W,i,j) verticalWell(W, G, rock, i, j, [], ...
                                         'type'  , 'rate', ...
                                         'val'   , rate  , ...
                                         'comp_i', [1,0] );
    producer = @(W,i,j) verticalWell(W, G, rock, i, j, [], ...
                                          'type'  , 'bhp', ...
                                          'val'   , bhp  , ...
                                          'comp_i', [1,0]);
    % Injection well in the center, production wells in each corner
    M = floor(N/2)+1;
    W = [];
    W = injector(W, M, M);
    W = producer(W, 1, 1);
    W = producer(W, N, 1);
    W = producer(W, N, N);
    W = producer(W, 1, N);
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(time, 30*day), 'W', W);
    % Initial state
    state0      = initResSol(G, bhp, [0,1]);
    % Defaul plotting
    plotOptions = {};
end