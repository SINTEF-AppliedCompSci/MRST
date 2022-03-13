function OSflux = findOSflux(G, rock, interpFace)
% Construct one-side fluxes for 2D and 3D grids

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    dispif(mrstVerbose, 'findOSflux... ');
    timer = tic;

    K = permTensor(rock, G.griddim);
    K = reshape(K', G.griddim, G.griddim, []);
    OSflux = cell(G.faces.num, 2);

    switch G.griddim

      case 2

        for i_face = 1 : G.faces.num

            if (all(G.faces.neighbors(i_face, :) ~= 0))
                % internal face
                c1 = G.faces.neighbors(i_face, 1);
                c2 = G.faces.neighbors(i_face, 2);
                K1 = K(:, :, c1);
                K2 = K(:, :, c2);
                w1 = K1 * G.faces.normals(i_face, :)';
                w2 = -K2 * G.faces.normals(i_face, :)';

                [a, faceA, faceB] = findAB(G, interpFace, c1, w1);
                interpA = [G.faces.neighbors(faceA, :)', interpFace.weights(faceA, :)'];
                interpB = [G.faces.neighbors(faceB, :)', interpFace.weights(faceB, :)'];
                interpA(:, 2) = -a(1) * interpA(:, 2);
                interpB(:, 2) = -a(2) * interpB(:, 2);
                container = [c1; c2; interpA(:, 1); interpB(:, 1); 0];
                container(:, 2) = [sum(a); 0; interpA(:, 2); interpB(:, 2); 0];
                trans = uniqueTrans(container);
                OSflux(i_face, 1) = {trans};
                clear trans;

                [a, faceA, faceB] = findAB(G, interpFace, c2, w2);
                interpA = [G.faces.neighbors(faceA, :)', interpFace.weights(faceA, :)'];
                interpB = [G.faces.neighbors(faceB, :)', interpFace.weights(faceB, :)'];
                interpA(:, 2) = -a(1) * interpA(:, 2);
                interpB(:, 2) = -a(2) * interpB(:, 2);
                container = [c2; c1; interpA(:, 1); interpB(:, 1); 0];
                container(:, 2) = [sum(a); 0; interpA(:, 2); interpB(:, 2); 0];
                trans = uniqueTrans(container);
                OSflux(i_face, 2) = {trans};
                clear trans;
            end
        end

      case 3

        for i_face = 1:G.faces.num

            if (all(G.faces.neighbors(i_face, :) ~= 0))
                % internal face
                c1 = G.faces.neighbors(i_face, 1);
                c2 = G.faces.neighbors(i_face, 2);

                K1 = K(:, :, c1);
                K2 = K(:, :, c2);
                w1 = K1 * G.faces.normals(i_face, :)';
                w2 = -K2 * G.faces.normals(i_face, :)';

                [a, faceA, faceB, faceC] = findABC(G, interpFace, c1, w1);
                interpA = [G.faces.neighbors(faceA, :)', -a(1) .* interpFace.weights(faceA, :)'];
                interpB = [G.faces.neighbors(faceB, :)', -a(2) .* interpFace.weights(faceB, :)'];
                interpC = [G.faces.neighbors(faceC, :)', -a(3) .* interpFace.weights(faceC, :)'];
                container = [c1; c2; interpA(:, 1); interpB(:, 1); interpC(:, 1); 0];
                container(:, 2) = [sum(a); 0; interpA(:, 2); interpB(:, 2); interpC(:, 2); 0];
                trans = uniqueTrans(container);
                OSflux(i_face, 1) = {trans};
                clear trans;

                [a, faceA, faceB, faceC] = findABC(G, interpFace, c2, w2);
                interpA = [G.faces.neighbors(faceA, :)', -a(1) .* interpFace.weights(faceA, :)'];
                interpB = [G.faces.neighbors(faceB, :)', -a(2) .* interpFace.weights(faceB, :)'];
                interpC = [G.faces.neighbors(faceC, :)', -a(3) .* interpFace.weights(faceC, :)'];
                container = [c2; c1; interpA(:, 1); interpB(:, 1); interpC(:, 1); 0];
                container(:, 2) = [sum(a); 0; interpA(:, 2); interpB(:, 2); interpC(:, 2); 0];
                trans = uniqueTrans(container);
                OSflux(i_face, 2) = {trans};
                clear trans;
            end
        end
    end

    fprintf('done in %1.2f s\n', toc(timer));
end

function [a, faceA, faceB] = findAB(G, interpFace, c, Kn)

    x1 = G.cells.centroids(c, :)';
    theFaces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1, 1);
    myBases = interpFace.coords(theFaces, :);
    myBases = bsxfun(@minus, myBases, x1');
    myNorm = sqrt(dot(myBases, myBases, 2));
    myBases = bsxfun(@rdivide, myBases, myNorm);
    Kn_norm = norm(Kn);
    Kn_unit = Kn / Kn_norm;
    myangles = bsxfun(@times, myBases, Kn_unit');
    myangles = sum(myangles, 2);
    myangles = acos(myangles);
    [~, I] = sort(myangles);
    theFaces = theFaces(I);
    myBases = myBases(I, :);
    myNorm = myNorm(I);
    nf = numel(theFaces);
    flag = 0;

    myIndex = zeros(nf*(nf - 1)/2, 2);
    myCoeff = myIndex;
    counter = 1;
    for i = 1:nf - 1
        tA = myBases(i, :)';
        tA_norm = myNorm(i);
        for j = i + 1:nf
            tB = myBases(j, :)';
            tB_norm = myNorm(j);
            if (abs(det([tA, tB])) > 1e-9)
                temp_a = [tA, tB] \ (Kn_unit);
                temp_a(abs(temp_a) < 1e-9) = 0;
                if (all(temp_a >= 0))
                    if (all(temp_a <= 1))
                        faceA = theFaces(i);
                        faceB = theFaces(j);
                        a = temp_a;
                        a(1) = a(1) * Kn_norm / tA_norm;
                        a(2) = a(2) * Kn_norm / tB_norm;
                        flag = 1;
                        break;
                    else
                        myIndex(counter, :) = [i, j];
                        myCoeff(counter, :) = temp_a;
                        counter = counter + 1;
                    end
                end
            end
        end
        if (flag), break; end
    end
    if (~flag && counter > 1)
        myIndex(counter:end, :) = [];
        myCoeff(counter:end, :) = [];
        maxCoeff = max(myCoeff, [], 2);
        [~, ind] = min(maxCoeff);
        i = myIndex(ind, 1);
        j = myIndex(ind, 2);
        a = myCoeff(ind, :);
        faceA = theFaces(i);
        faceB = theFaces(j);
        tA_norm = myNorm(i);
        tB_norm = myNorm(j);
        a(1) = a(1) * Kn_norm / tA_norm;
        a(2) = a(2) * Kn_norm / tB_norm;
    end
    assert(logical(exist('faceA', 'var')), ...
           ['decomposition failed for cell ', num2str(c)]);

end

function [a, faceA, faceB, faceC] = findABC(G, interpFace, c, Kn)
    x1 = G.cells.centroids(c, :)';
    theFaces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1, 1);
    myBases = interpFace.coords(theFaces, :);
    myBases = bsxfun(@minus, myBases, x1');
    myNorm = sqrt(dot(myBases, myBases, 2));
    myBases = bsxfun(@rdivide, myBases, myNorm);
    Kn_norm = norm(Kn);
    Kn_unit = Kn / Kn_norm;
    myangles = bsxfun(@times, myBases, Kn_unit');
    myangles = sum(myangles, 2);
    myangles = acos(myangles);
    [~, I] = sort(myangles);
    theFaces = theFaces(I);
    myBases = myBases(I, :);
    myNorm = myNorm(I);
    nf = numel(theFaces);
    flag = 0;

    myIndex = zeros(nf*(nf - 1)*(nf - 2)/6, 3);
    myCoeff = myIndex;
    counter = 1;
    for i = 1:nf - 2
        tA = myBases(i, :)';
        tA_norm = myNorm(i);
        for j = i + 1:nf - 1
            tB = myBases(j, :)';
            tB_norm = myNorm(j);
            for k = j + 1:nf
                tC = myBases(k, :)';
                tC_norm = myNorm(k);
                if (abs(det([tA, tB, tC])) > 1e-9)
                    temp_a = [tA, tB, tC] \ (Kn_unit);
                    temp_a(abs(temp_a) < 1e-9) = 0;
                    if (all(temp_a >= 0))
                        if (all(temp_a <= 1))
                            faceA = theFaces(i);
                            faceB = theFaces(j);
                            faceC = theFaces(k);
                            a = temp_a;
                            a(1) = a(1) * Kn_norm / tA_norm;
                            a(2) = a(2) * Kn_norm / tB_norm;
                            a(3) = a(3) * Kn_norm / tC_norm;
                            flag = 1;
                            break;
                        else
                            myIndex(counter, :) = [i, j, k];
                            myCoeff(counter, :) = temp_a;
                            counter = counter + 1;
                        end
                    end
                end
            end
            if (flag), break; end
        end
        if (flag), break; end
    end

    if (~flag && counter > 1)
        myIndex(counter:end, :) = [];
        myCoeff(counter:end, :) = [];
        maxCoeff = max(myCoeff, [], 2);
        [~, ind] = min(maxCoeff);
        i = myIndex(ind, 1);
        j = myIndex(ind, 2);
        k = myIndex(ind, 3);
        a = myCoeff(ind, :);
        faceA = theFaces(i);
        faceB = theFaces(j);
        faceC = theFaces(k);
        tA_norm = myNorm(i);
        tB_norm = myNorm(j);
        tC_norm = myNorm(k);
        a(1) = a(1) * Kn_norm / tA_norm;
        a(2) = a(2) * Kn_norm / tB_norm;
        a(3) = a(3) * Kn_norm / tC_norm;
    end

    % Error
    if ~exist('faceA', 'var') & mrstVerbose
        figure, hold on
        plotGrid(G, 'facealpha', 0.1);
        plotGrid(G, c)

        figure, hold on
        plotGrid(G, c, 'facealpha', 0.3)
        hap = interpFace.coords(theFaces, :);
        plot3(hap(:, 1), hap(:, 2), hap(:, 3), '.', 'markersize', 14)
        xc = G.cells.centroids(c, :);
        plot3(xc(1), xc(2), xc(3), 'ko', 'markersize', 14)
        ind = convhull(hap);
        trisurf(ind, hap(:, 1), hap(:, 2), hap(:, 3), 'Facecolor', 'cyan', 'facealpha', 0.5)

        tol = -1e-5;
        in_convex_hull = mex_inhull(xc, hap, ind, tol);
        fprintf('mex_inhull for cell %d with tol=%f: %d\n', c, tol, in_convex_hull);

        % Plot Kn.
        Gc = extractSubgrid(G, c);
        xmin = min(Gc.nodes.coords);
        xmax = max(Gc.nodes.coords);
        h = xmax - xmin;
        vec = [xc; xc + Kn_unit' .* h];
        line(vec(:, 1), vec(:, 2), vec(:, 3))
    end

    assert(exist('faceA', 'var'), ...
           sprintf('decomposition failed for cell %d\n', c));
end

function [a, xD] = findDnode(G, mycell, myface, Kn)

    n1 = G.faces.nodes(G.faces.nodePos(myface));
    n2 = G.faces.nodes(G.faces.nodePos(myface)+1);
    xn1 = G.nodes.coords(n1, :)';
    xn2 = G.nodes.coords(n2, :)';
    xf = G.faces.centroids(myface, :)';
    xc = G.cells.centroids(mycell, :)';
    Kn_norm = norm(Kn);
    Kn = Kn / Kn_norm;
    t_norm = norm(xc-xf);
    t = (xc - xf) / t_norm;
    t1_norm = norm(xn1-xf);
    t1 = (xn1 - xf) / t1_norm;
    t2_norm = norm(xn2-xf);
    t2 = (xn2 - xf) / t2_norm;
    temp_a = [t, t1] \ Kn;
    temp_a(abs(temp_a) < 1e-9) = 0;
    if (all(temp_a >= 0))
        a = temp_a;
        a(1) = a(1) * Kn_norm / t_norm;
        a(2) = a(2) * Kn_norm / t1_norm;
        xD = xn1;
    else
        a = [t, t2] \ Kn;
        a(abs(a) < 1e-9) = 0;
        a(1) = a(1) * Kn_norm / t_norm;
        a(2) = a(2) * Kn_norm / t2_norm;
        xD = xn2;
    end
end

function [a, xA, xB] = findDnodes(G, mycell, myface, Kn)
    mynodes = G.faces.nodes(G.faces.nodePos(myface):G.faces.nodePos(myface+1)-1);
    mynodes = [mynodes; mynodes(1)];
    xnode = G.nodes.coords(mynodes, :);
    xc = G.cells.centroids(mycell, :)';
    xf = G.faces.centroids(myface, :)';
    tc_norm = norm(xc-xf);
    tc = (xc - xf) / tc_norm;
    Kn_norm = norm(Kn);
    Kn = Kn / Kn_norm;
    for i = 1:numel(mynodes) - 1
        xA = xnode(i, :)';
        xB = xnode(i+1, :)';
        tA_norm = norm(xA-xf);
        tA = (xA - xf) / tA_norm;
        tB_norm = norm(xB-xf);
        tB = (xB - xf) / tB_norm;
        a = [tc, tA, tB] \ Kn;
        a(abs(a) < 1e-9) = 0;
        if (all(a >= 0))
            a(1) = a(1) * Kn_norm / tc_norm;
            a(2) = a(2) * Kn_norm / tA_norm;
            a(3) = a(3) * Kn_norm / tB_norm;
            break;
        end
    end
end

function [trans] = uniqueTrans(container)
    [trans, ~, subs] = unique(container(:, 1), 'rows', 'stable');
    trans(:, 2) = accumarray(subs, container(:, 2));
    trans(3:end, :) = sortrows(trans(3:end, :), -1);
    trans(2:end, 2) = -trans(2:end, 2);
end
