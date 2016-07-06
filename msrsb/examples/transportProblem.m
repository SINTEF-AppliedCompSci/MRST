mrstModule add spe10 multiscale-devel coarsegrid libgeometry mrst-gui mrst-experimental agglom

% caseName = 'ness';
caseName = 'tarbert';
% caseName = 'uniform';
n = lower(caseName);

switch n
    case {'ness', 'tarbert'}
        if strcmp(n, 'ness')
            dt = 30*day;
            layer = 71;
        else
            dt = 100*day;
            layer = 1;
        end
        [G, ~, rock] = SPE10_setup(layer);
        
        dims = G.cartDims;
        pdims = max(G.nodes.coords) - min(G.nodes.coords);
        G = cartGrid(dims(1:2), pdims(1:2));
        
        
        rock.perm = rock.perm(:, 1);
    case {'uniform', 'anisotropy'};
        dt = 1*day;
        
        %         dims = [100, 200];
        dims = [50, 100];
        %         dims = [30, 60];
        G = cartGrid(dims, [50, 100]*meter);
        
        if strcmp(n, 'anisotropy')
            rock.perm = repmat([1/10, 1]*milli*darcy, G.cells.num, 1);
            dt = 100*day;
        else
            rock.perm = repmat(500*milli*darcy, G.cells.num, 1);
        end
        rock.poro = ones(G.cells.num, 1);
    otherwise
        error('Unknown case...')
end
G = computeGeometry(G);
nstep = 150;



useBC = false;
useWells = true;

bc = [];
W = [];
if useBC
    % injrate = poreVolume(G, rock)/(5000*day*G.cartDims(1));
    % bc = fluxside(bc, G, 'YMin', repmat(injrate, G.cartDims(1)), 'sat', [1 0]);
    bc = pside(bc, G, 'YMin', 1000*barsa, 'sat', [1 0]);
    bc = pside(bc, G, 'YMax', 0*barsa, 'sat', [1 0]);
end

if useWells
    offset = 0;
    %     offset = 10;
    W = verticalWell(W, G, rock, 1+offset, 1+offset, [], ...
        'Val', 1000*barsa, 'Type', 'bhp', 'Comp_i', [1 0]);
    W = verticalWell(W, G, rock, G.cartDims(1)-offset, G.cartDims(2)-offset, [],...
        'Val', 0*barsa, 'Type', 'bhp', 'Comp_i', [1 0]);
end


% rock.perm = 0*rock.perm + 0.5*darcy;
% rock.poro = 0*rock.poro + 1;

T = computeTrans(G, rock);


mp = 0.01;
rock.poro(rock.poro < mp) = mp;

figure;
plotToolbar(G, rock); axis equal tight off;
view(90, 90)



%%
% cdims = ceil(G.cartDims./[5 10 5]);
cdims = ceil(G.cartDims./[10 20]);
% cdims = ceil(G.cartDims./[10 10]);

uniformPart = true;

useMetis = false;
if ~useMetis
    p = partitionUI(G, cdims);
else
    p = partitionMETIS(G, T, prod(cdims));
    uniformPart = false;
end

if 1
    if useBC
        [ii, jj] = gridLogicalIndices(G);
        % Add a thin top layer of partitions to honor no-flow correctly
        %     noFlowEdgesJ = jj == max(jj) | jj == min(jj);
        noFlowEdgesJ = jj < 3;
        
        if 1
            % Add a few
            p(noFlowEdgesJ) = p(noFlowEdgesJ) + max(p);
        else
            p_refine = partitionUI(G, G.cartDims);
            % Add many
            p(noFlowEdgesJ) = p_refine(noFlowEdgesJ) + max(p);
        end
        
        %     noFlowEdgesI = ii == max(ii) | ii == min(ii);
        %     p(noFlowEdgesI) = p(noFlowEdgesI) + max(p);
        
        p = compressPartition(p);
        uniformPart = false;
    end
    
    if useWells
        wc = vertcat(W.cells);
        
        p_refine = partitionUI(G, G.cartDims);
        if 1
            p0 = p;
            [isRefined, blorg] = deal(false(G.cells.num, 1));
            fl = ones(G.cells.num, 1);
            for i = 1:numel(W)
                wc = W(i).cells(1);
                
                if 1
                    if i == 1
                        pt_well = min(G.nodes.coords);
                    else
                        pt_well = max(G.nodes.coords);
                    end
                else
                    pt_well = G.cells.centroids(wc, :);
                end
                
                [ii, jj] = gridLogicalIndices(G);
                
                maxi = 25;
                maxj = 50;
                local = abs(ii - ii(wc)) < maxi & abs(jj - jj(wc)) < maxj;
                cells = find(ismember(p, p(local)));
                
                
                pts = G.cells.centroids(cells, :);
                angles = [8 8 16 16];
                radiusBins = numel(angles);
                [out, indSector, indBin] = refineNearWell(pts, pt_well, ...
                    'angleBins', angles, 'radiusBins', radiusBins, 'logBins', true);
                
                inner = indBin < max(indBin);

                p(cells(inner)) = max(p) + out(inner);
                isRefined(cells(inner)) = true;
                
                
                %%%%%%%%
                locPart = ismember(p0, p0(isRefined));
                sec = indSector(:, 1);
                blorg(locPart) = true;
                tmp = ones(G.cells.num, 1);
                tmp(cells) = sec;
                p(locPart) = p(locPart) + tmp(locPart) + max(p);
                
                fl(cells(sec == min(sec))) = 10;
                fl(cells(sec ~= min(sec))) = 1;
                %%%%

            end
            
            if 1
                [ic, jc] = ind2sub(cdims, (1:G.cells.num)');
                
                
                bnd = jc == 1 | jc == max(jc) | ic == 1 | ic == max(ic);
                fl(ismember(p0, find(bnd))) = 1;
                % Hand tuned agglom stuff
                vol = 0.001*ones(G.cells.num, 1);
                vol(isRefined) = 1;
                p = processPartition(G, p);
                % lim = 0.06;
                lim = .8;
                p = mergeBlocks(p, G, vol, fl, lim);
                %
            end
            
            if 0
                p0 = p;
                numBins  = 10;
                NL       = 50;
                NU       = 100;
                tof  = computeTimeOfFlight(ref, G, rock, 'wells', W);
                Tr = computeTimeOfFlight(ref, G, rock, 'wells', W, 'reverse', true);
                I = -log10(tof.*Tr); I = I - min(I) + 1;
                
                p = segmentIndicator(G, I, numBins);
                p = mergeBlocks2(p, G, rock.poro, I, NL, NU);
                p = refineBlocks(p, G, I, NU, @refineGreedy2);
                p = mergeBlocks2(p, G, rock.poro, I, NL, NU);
                max(p)
                
                %                 [bfi, bfc] = boundaryFaces(G, find(isRefined));
                %                 [cellNo, cf, cn_isnnc] = getCellNoFaces(G);
                
                %                 T0 = T;
                %                 FT = 1 ./ accumarray(cf, 1 ./ T0, [G.faces.num, 1]);
                %                 FT(bfi) = 0;
                %                 p = partitionMETIS(G, FT, prod(cdims), 'no2hop', true, 'ufactor', 1.05, 'ncuts', 5);
                p(isRefined) = p0(isRefined);
            end
            
            p(vertcat(W.cells)) = max(p) + 1;
            uniformPart = false;
        end
    end
end

p = processPartition(G, p);
p = compressPartition(p);

figure;
plotCellData(G, mod(p, 13), 'edgec', 'w', 'edgea', .1);
outlineCoarseGrid(G, p)
axis equal tight off
view(90, 90);
%%
name = [caseName, '_', num2str(prod(cdims)), 'dof'];
if useMetis
    name = [name, '_metis'];
end


if 0
    fluid = initSimpleFluid('mu',  [1, 1]*centi*poise, ...
        'rho', [1000, 1000]*kilogram/meter^3, ...
        'n',   [1, 1]);
elseif 0
    sharpfront = false;
    fluid = initSimpleFluid('mu',  [1, 10]*centi*poise, ...
        'rho', [1024, 700]*kilogram/meter^3, ...
        'n',   [2, 2]);
    dt = 100*day;
else
    sharpfront = true;
    name = [name, '_visc10'];
    fluid = initSimpleFluid('mu',  [10, 1]*centi*poise, ...
        'rho', [1024, 700]*kilogram/meter^3, ...
        'n',   [2, 2]);
    dt = 400*day;
    nstep = 300;
end
state0 = initResSol(G, 0, [0 1]);

ref = incompTPFA(state0, G, T, fluid, 'Wells', W);

%%
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
if uniformPart
    CG = storeInteractionRegionCart(CG);
else
    CG = storeInteractionRegion(CG, 'adjustCenters', true, 'localTriangulation', true);
end


A = getIncomp1PhMatrix(G, T, state0, fluid);

sg = setupGridsForMex(A, CG);
[I, I_compressed] = mbasisSmoothed(sg, 'tolerance', 1e-3, 'maxiter', 300);
X = restrictOperator(CG);


basis = struct('B', I, 'R', X, 'B_compressed', I_compressed);

figure;
plotCellData(G, mod(p, 13));
axis equal tight off
view(90, 90);
plotGrid(G, sg.globalBoundary + 1)

%%


psolve = @(state) incompTPFA(state, G, T, fluid, 'bc', bc, 'Wells', W);
psolvems = @(state, basis) solveMultiscaleTPFA(state, G, CG, T, fluid, 'bc', bc, 'Wells', W,...
    'basis', basis, 'solver', 'jacobi', 'reconstruct', true, ...
    'iterations', 0, 'smoothsteps', 10, 'iterator', 'jacobi', 'Verbose', false);

tsolve = @(state, bc, dt) implicitTransport(state, G, dt, rock, fluid, 'bc', bc, 'Wells', W, 'nltol', 1e-5);

% solvems   = @(state, dt, basis) tsolve(psolvems(state, basis), dt);
% solvefine = @(state, dt)        tsolve(psolve(  state),        dt);


%%
basis_adaptive = basis;


computeAdaptive = true;

tmp = cell(nstep+1, 1);
tmp{1} = state0;

[ms, msadaptive, fine] = deal(tmp);

% dt = 100*day;

for i = 1:nstep
    disp(['Step ', num2str(i), ' of ', num2str(nstep)]);
    disp(['[', repmat('*', 1, i), repmat('.', 1, nstep - i), ']', ]);
    % Fine
    state = fine{i};
    state = psolve(state);
    % Transport
    state = tsolve(state, bc, dt);
    fine{i+1} = state;
    
    % MS
    state = ms{i};
    state = psolvems(state, basis);
    % Transport
    [bc_transport, state] = incompressibleFluxBC(G, bc, state);
    state = tsolve(state, bc_transport, dt);
    ms{i+1} = state;
    
    if computeAdaptive
        % MS with adaptive updates
        state = msadaptive{i};
        A = getIncomp1PhMatrix(G, T, state, fluid);
        
        sg = updateMatricesForMexGrid(A, CG, sg);
        
        [basis_adaptive.B, basis_adaptive.B_compressed] = mbasisSmoothed(sg, ...
            'tolerance', 5e-3, 'maxiter', 30, 'Interpolator', basis_adaptive.B_compressed);
        
        state = psolvems(state, basis_adaptive);
        
        % Transport
        [bc_transport, state] = incompressibleFluxBC(G, bc, state);
        state = tsolve(state, bc_transport, dt);
        msadaptive{i+1} = state;
    else
        % Fake it.
        msadaptive{i+1} = ms{i+1};
    end
end
%%
ms = ms(cellfun(@(x) ~isempty(x), ms));
fine = fine(cellfun(@(x) ~isempty(x), fine));
msadaptive = msadaptive(cellfun(@(x) ~isempty(x), msadaptive));

% getd = @(x) cellfun(@(y) y.s, x, 'UniformOutput', false);
% getd = @(x) cellfun(@(y) y.pressure, x, 'UniformOutput', false);
close all

figure;
plotToolbar(G, rock); axis equal tight off;
outlineCoarseGrid(G, CG.partition, 'w')
view(90, 90)

for i = 1:2
    if i == 1
        fn = 's';
    else
        fn = 'pressure';
    end
    
    getd = @(x) cellfun(@(y) y.(fn), x, 'UniformOutput', false);
    
    figure;
    plotToolbar(G, getd(fine(2:end))); axis equal tight off;
    view(90, 90)
    title(['Ref: ', fn]);
    
    figure
    plotToolbar(G, getd(ms(2:end))); axis equal tight off;
    view(90, 90)
    title(['MS: ', fn]);
    
    if computeAdaptive
        figure
        plotToolbar(G, getd(msadaptive(2:end))); axis equal tight off;
        view(90, 90)
        title(['Adaptive: ', fn]);
    end
end
%%
% pv = poreVolume(G, rock)

err = @(fine, ms) cellfun(@(x, y) abs(x.s(:, 1) - y.s(:, 1)).*rock.poro, fine, ms, 'unif', false);
figure
e_ms = err(fine, ms);
plotToolbar(G, e_ms, 'EdgeColor', 'none'); axis equal tight off;
% outlineCoarseGrid(G, CG.partition)
view(90, 90)
colorbar
c = caxis();
title('MS error');

if computeAdaptive
    figure
    e_adaptive = err(fine, msadaptive);
    plotToolbar(G, e_adaptive, 'EdgeColor', 'none'); axis equal tight off;
    view(90, 90)
    colorbar
    caxis(c);
    title('Adaptive error');
    
    figure;
    enorm = [cellfun(@norm, e_ms), cellfun(@norm, e_adaptive)];
    enorm = bsxfun(@rdivide, enorm, cellfun(@(x) norm(x.s(:, 1).*rock.poro), fine));
    plot(enorm, 'linewidth', 2)
    legend('MS', 'MS-adaptive')
    grid
    axis tight
end
%%

mrstModule add export_fig

msdir = mrstPath('query', 'multiscale-devel');

fdir = fullfile(msdir, 'papers', 'mrsb', 'paper', 'figures', 'transport2d');
if ~exist(fdir, 'dir')
    mkdir(fdir)
end

saveplot = @(h, ext, varargin) export_fig(fullfile(fdir, [name, '_', get(h, 'Name'), '.', ext]), h, '-transparent', varargin{:});

%% plot saturation front
close all;

bpoints = [.1, .5, 1];
totTime = nstep*dt;
stepT = cumsum(ones(nstep, 1)*dt);

fieldnames = {'pressure', 's'};
                
df = get(0, 'DefaultFigurePosition');
getfig = @() figure('Position', [df(1)*0.8, df(2), 1000, 400]);

for i = 1:numel(bpoints)
    index = find(stepT >= totTime*bpoints(i), 1, 'first');
    
    for j = 1:numel(fieldnames)
        fn = fieldnames{j};
        for k = 1:3
            if k == 1
                n = 'ref';
                state = fine{index};
                stateref = state;
            elseif k == 2
                n = 'ms';
                state = ms{index};
            else
                n = 'adaptive';
                state = msadaptive{index};
            end
            
            sn = [n, '_', num2str(100*bpoints(i)), '_', fn];
            
            h = getfig();
            plotCellData(G, state.(fn)(:, 1), 'edgec', 'none');
            if k == 1
                c = caxis();
            end
            caxis(c);
%             if i == numel(bpoints)
%                 colorbar;
%             end
            view(90, 90);
            axis equal tight off
            
            set(h, 'Name', sn);
            if k > 1 && strcmpi(fn, 's')
                e = state.(fn)(:, 1) - stateref.(fn)(:, 1);
                h = getfig();
                plotCellData(G, abs(e), 'edgec', 'none');
                caxis([0, 0.25]);
                
                if i == numel(bpoints)
                    colorbar;
                end
                view(90, 90);
                axis equal tight on
                set(h, 'Name', [sn, '_err']);
            end
        end
    end
end
%%
ch = get(0, 'Children');
for i = 1:numel(ch)
    saveplot(ch(i), 'png')
end
%%
close all
x = reshape(G.cells.centroids(:, 1), G.cartDims);
y = reshape(G.cells.centroids(:, 2), G.cartDims);

makeContour = @(data, v, varargin) contour(x, y, reshape(data, G.cartDims), v, varargin{:}); 
for i = 1:numel(bpoints)
    figure;
    index = find(stepT >= totTime*bpoints(i), 1, 'first');
    ref = fine{index}.s(:,1); 
    msa = msadaptive{index}.s(:, 1);
    
    % Draw reference as the background
    plotCellData(G, ref, 'EdgeColor', 'none')
    hold on;
    
    v = 0.1:0.1:1;
    % First drawing - make single color outline for contours
    makeContour(msa, v, 'k', 'linewidth', 4)
    % Second drawing - colorized by same colormap as the background
    % saturation front
    makeContour(msa, v, 'linewidth', 2)
    view(90, 90)
    axis equal tight off
    
    set(gcf, 'Name', [num2str(100*bpoints(i)), '_satfront']);
end

ch = get(0, 'Children');
for i = 1:numel(ch)
    saveplot(ch(i), 'png')
end
%%

cmap = gray(128);
cmap = cmap(65:end, :);
cmap = flipud(cmap);


close all
makeContour = @(data, v, varargin) contour(x, y, reshape(data, G.cartDims), v, varargin{:}); 
for i = 1:numel(bpoints)
    getfig();

    index = find(stepT >= totTime*bpoints(i), 1, 'first');
    ref = fine{index}.s(:,1); 
    msa = msadaptive{index}.s(:, 1);
    ms0 = ms{index}.s(:, 1);
    % Draw reference as the background
%     d = log10(rock.perm(:,1));
%     d = d-min(d); d = d./max(d);
%     plotCellData(G, d, 'edgec', 'none')
%     plotGrid(G, 'facec', 'w', 'edgec', 'none')
    plotCellData(G, ref, 'EdgeColor', 'none')
    hold on;
    
    if sharpfront
        v = 0.5;
        cax = [0 1];
    else
        dv = 0.9/5;
        v = 0.1:dv:1;
        cax = [0 1];
    end
    lc = 2;
    c1 = 'b';
    c2 = 'g';
    c3 = 'r';
%     c1 = [255, 0, 0]./255;
%     c2 = [0, 255, 127]./255;
%     c3 = [15, 53, 232]./255;
    makeContour(msa, v, 'linewidth', lc, 'color', c1)
    makeContour(ms0, v, 'linewidth', lc, 'color', c2)
    makeContour(ref, v, '--', 'linewidth', lc, 'color', c3)
    colormap(cmap);

    caxis(cax)
    view(90, 90)
    axis equal tight off
    
    set(gcf, 'Name', [num2str(100*bpoints(i)), '_satfront']);
end
%%
ch = get(0, 'Children');
for i = 1:numel(ch)
    saveplot(ch(i), 'png')
end
%%
figure;
hold on
plot(1, 1, 'linewidth', lc, 'color', c1)
plot(1, 1, 'linewidth', lc, 'color', c2)
plot(1, 1, '--', 'linewidth', lc, 'color', c3)
legend('MS-adaptive', 'MS-static', 'Reference', 'orientation', 'horizontal')
axis tight off
set(gcf, 'Name', 'contourlegend')
% saveplot(gcf, 'pdf')

%%
figure;
h = colorbar();
set(h, 'Position', [.1 .05 .5 .9 ])
colormap(cmap)
caxis(cax);
% axis tight off
set(gca, 'Visible', 'off')
set(gcf, 'Name', 'caxis_sat_gray')
% saveplot(gcf, 'pdf')

%%
h = getfig();
plotCellData(G, log10(rock.perm(:, 1)), 'edgec', 'none');
outlineCoarseGrid(G, p, 'k', 'LineWidth', 2)
axis equal tight off
view(90, 90);
set(h, 'Name', 'perm')
saveplot(h, 'png')
%%
h = getfig();
plotCellData(G, mod(p, 13), 'edgec', 'w', 'edgea', .1);
outlineCoarseGrid(G, p)
axis equal tight off
view(90, 90);
set(h, 'Name', 'coarsegrid')
saveplot(h, 'png')
%%
if computeAdaptive
    fn = @(x) norm(x, 2);
    
%     h = getfig();
    h = figure;
    enorm = [cellfun(fn, e_ms), cellfun(fn, e_adaptive)];
    enorm = bsxfun(@rdivide, enorm, cellfun(@(x) fn(x.s(:, 1).*rock.poro), fine));
    plot(enorm, 'linewidth', 2)
    legend('Static basis', 'Adaptive basis')
    grid
    set(gcf, 'Name', 'saturation_error');
    axis tight
    saveplot(h, 'pdf')
end
