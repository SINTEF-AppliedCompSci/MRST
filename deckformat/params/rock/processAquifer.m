function output = processAquifer(deck, G)
% return array aquifer. Each row corresponds to a cell connected to the aquifer.
% the aquind structure gives the content of each column:

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

    aquind.aquid     = 1; % Aquifer numer (id) 
    aquind.conn      = 2; % Connections (cell number)
    aquind.pvttbl    = 3; % PVT table number
    aquind.J         = 4; % Productivity index
    aquind.C         = 5; % Aquifer compressibility
    aquind.alpha     = 6; % alpha coefficient
    aquind.depthconn = 7; % Connection's depth
    aquind.depthaq   = 8; % Aquifer's depth

    aquancon = deck.SOLUTION.AQUANCON;
    aqufetp  = deck.SOLUTION.AQUFETP;
    
    % Check that aquifer pressures are not defaulted (currently not supported)
    assert(~any(isnan(aqufetp(:,3))), 'Defaulting aquifer pressure is not supported')
    
    aquid = cell2mat(aquancon(:, 1));
    imin = cell2mat(aquancon(:, 2));
    imax = cell2mat(aquancon(:, 3));
    jmin = cell2mat(aquancon(:, 4));
    jmax = cell2mat(aquancon(:, 5));
    kmin = cell2mat(aquancon(:, 6));
    kmax = cell2mat(aquancon(:, 7));

    aquiferdirstr = aquancon(:, 8);
    % note: should add support for X,Y,Z afterwards
    dirstrs = {'I-', 'I+', 'J-', 'J+', 'K-', 'K+'};
    dirflags = {1, 2, 3, 4, 5, 6};
    nentries = numel(aquiferdirstr);
    faceflag = zeros(nentries, 1);
    for i = 1 : numel(dirstrs)
        ind = strcmp(dirstrs{i}, aquiferdirstr);
        faceflag(ind) = dirflags{i};
    end

    influxcoef = cell2mat(aquancon(:, 9));
    influxmultcoef = cell2mat(aquancon(:, 10));
    
    aqutbl = [imin, imax, jmin, jmax, kmin, kmax, ...
             aquid, influxcoef, influxmultcoef, faceflag];
    for i = 1 : 3
        u1 = mcolon(aqutbl(:, 1), aqutbl(:, 2))';
        u2 = rldecode(aqutbl(:, 3 : end), aqutbl(:, 2) - aqutbl(:, 1) + 1);
        aqutbl = [u2, u1];
    end
    
    % Reorder aqutbl and, for readability, set up a structure which maps to the
    % indices.
    ind = [(5 : 7), (1 : 4)];
    aqutbl= aqutbl(:, ind);
    tblind.i = 1;
    tblind.j = 2;
    tblind.k = 3;
    tblind.aquid = 4;
    tblind.influxcoef = 5;
    tblind.influxmultcoef = 6;
    tblind.faceflag = 7;

    % Remove the cells that are not active
    cells = sub2ind(G.cartDims, ...
                    aqutbl(:, tblind.i), ...
                    aqutbl(:, tblind.j), ...
                    aqutbl(:, tblind.k));
    actnum = false([prod(G.cartDims), 1]);
    actnum(G.cells.indexMap) = true;
    aqutbl = aqutbl(actnum(cells), :);

    % Find the neighbors and keep them that are outside active domain
    N = size(aqutbl, 1);
    neigcelltbl = aqutbl(:, [tblind.i, tblind.j, tblind.k]);
    ind = aqutbl(:, tblind.faceflag) == 1;
    neigcelltbl(ind, 1) = neigcelltbl(ind, 1) - 1;
    ind = aqutbl(:, tblind.faceflag) == 2;
    neigcelltbl(ind, 1) = neigcelltbl(ind, 1) + 1;
    ind = aqutbl(:, tblind.faceflag) == 3;
    neigcelltbl(ind, 2) = neigcelltbl(ind, 2) - 1;
    ind = aqutbl(:, tblind.faceflag) == 4;
    neigcelltbl(ind, 2) = neigcelltbl(ind, 2) + 1;
    ind = aqutbl(:, tblind.faceflag) == 5;
    neigcelltbl(ind, 3) = neigcelltbl(ind, 3) - 1;
    ind = aqutbl(:, tblind.faceflag) == 6;
    neigcelltbl(ind, 3) = neigcelltbl(ind, 3) + 1;

    outlow = ones(N, 3);
    outup = repmat(G.cartDims, N, 1);
    aqcelloutind = (neigcelltbl < outlow) |  (neigcelltbl > outup);
    [I, J] = find(aqcelloutind);
    aqcelloutind = false(N, 1);
    aqcelloutind(I) = true;
    aqcellout = aqutbl(aqcelloutind, :);
    aqutbl = aqutbl(~aqcelloutind, :);
    neigcelltbl = neigcelltbl(~aqcelloutind, :);
    neigcellind = sub2ind(G.cartDims, ...
                          neigcelltbl(:, 1), ...
                          neigcelltbl(:, 2), ...
                          neigcelltbl(:, 3));
    actneigcell = actnum(neigcellind);
    aqcellinind = (actneigcell == 0);
    aqcellin = aqutbl(aqcellinind, :);

    aqutbl = [aqcellout; aqcellin];

    ind = sub2ind(G.cartDims, ...
                  aqutbl(:, tblind.i), ...
                  aqutbl(:, tblind.j), ...
                  aqutbl(:, tblind.k));
    conn = cart2active(G, ind);

    tblind.conn = 8;
    aqutbl(:, tblind.conn) = conn;

    % Set up dispatch mapping from aquifer number to connection 
    nconn = size(aqutbl, 1);
    naq = max(aqufetp(:, 1));
    aquid2conn = sparse(aqutbl(:, tblind.aquid), (1 : nconn)', 1, ...
                        naq, nconn)';
    
    % find aquifer faces, get areas and update defaulted influx coefficients
    % with areas
    conn = aqutbl(:, tblind.conn);
    faces = G.cells.faces(mcolon(G.cells.facePos(conn), G.cells.facePos(conn ...
                                                      + 1) - 1)', :);
    faceflags = faces(:, 2);
    diffcellfaces = diff(G.cells.facePos);
    diffcellfaces = diffcellfaces(conn);
    cellfaces = rldecode(conn, diffcellfaces);

    nc = G.cells.num;
    map11 = sparse(conn, (1 : numel(conn))', 1, nc, numel(conn));
    map12 = sparse(cellfaces, (1 : numel(cellfaces))', 1, nc, numel(cellfaces));
    map1 = map12'*map11;

    nflag = 6;
    aqflags = aqutbl(:, tblind.faceflag);
    map21 = sparse(aqflags, (1 : numel(aqflags))', 1, nflag, numel(aqflags));
    map22 = sparse(faceflags, (1 : numel(faceflags))', 1, nflag, numel(faceflags));
    map2 = map22'*map21;

    map = map1.*map2;
    [indfaces, indaq] = find(map);

    nconn = numel(conn);
    aqfaces = zeros(nconn, 1);
    aqfaces(indaq) = faces(indfaces);
    aqareas = G.faces.areas(aqfaces);

    influxcoef = aqutbl(:, tblind.influxcoef);
    useArea = isnan(influxcoef);
    influxcoef(useArea) = aqareas(useArea);
    influxmultcoef = aqutbl(:, tblind.influxmultcoef);
    alpha = influxmultcoef.*influxcoef;
    sumalpha = aquid2conn*aquid2conn'*alpha;
    alpha = alpha./sumalpha;
    
    aquiferprops.depthaq = aqufetp(:, 2);
    aquiferprops.C       = aqufetp(:, 5);
    aquiferprops.J       = aqufetp(:, 6);
    aquiferprops.pvttbl  = aqufetp(:, 7);
    
    depthaq = aquid2conn*aquiferprops.depthaq;
    C       = aquid2conn*aquiferprops.C;
    J       = aquid2conn*aquiferprops.J;
    pvttbl  = aquid2conn*aquiferprops.pvttbl;
    
    depthconn = G.cells.centroids(conn, 3);
    aquifers = NaN(nconn, 7);

    aquifers(:, aquind.aquid)     = aqutbl(:, tblind.aquid);
    aquifers(:, aquind.conn)      = aqutbl(:, tblind.conn);
    aquifers(:, aquind.pvttbl)    = pvttbl;
    aquifers(:, aquind.J)         = J;
    aquifers(:, aquind.C)         = C;
    aquifers(:, aquind.alpha)     = alpha;
    aquifers(:, aquind.depthconn) = depthconn;
    aquifers(:, aquind.depthaq)   = depthaq;
    
    initval.pressures = aqufetp(:, 3);
    initval.volumes   = aqufetp(:, 4);

    output = struct('aquifers'    , aquifers, ...
                    'aquind'      , aquind, ...
                    'initval'     , initval, ...
                    'aquiferprops', aquiferprops);

end
