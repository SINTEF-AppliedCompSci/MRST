function [interpFace] = correctHAP(G, interpFace, myRatio)
% Correct ill-placed harmonic averaging points. If the number of input arguments
% is 2, then the correction algorithm is applied only when some cell centroids
% lie outside their associated convex hull; if the number of input arguments is
% 3, then the last input argument myRatio is applied to all the harmonic
% averaging points.

%  G - Grid structure of MRST
%  interpFace - harmonic averaging point interplation without correction
%  myRatio - user specified ratio

    dispif(mrstVerbose, 'correctHAP\n');            
    
    HAP = interpFace.coords; % store the locations of the original harmonic averaging points;
    if (nargin == 2 || (nargin == 3 && isempty(myRatio)))
        if (interpFace.fraction > 0)
            if (G.griddim == 2)
                R = 0.5 * G.faces.areas;
            else
                R = sqrt(G.faces.areas./pi);
            end

            % Find all cells with centroids not in the convex hull
            out = find_cells(G, 1:G.cells.num, interpFace);
            
            % Correct
            cells = 1:G.cells.num;
            cells = cells(out);
            
            while sum(out)
                for c = cells
                    while out(c)
                        %fprintf('%d %d\n', c, sum(out))
                        interpFace = correct(G, c, interpFace, R, HAP);
                        out(c) = find_cells(G, c, interpFace);
                    end
                    faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
                    neighs = G.faces.neighbors(faces, 2);
                    neighs = neighs(neighs > 0);
                    %neighs = neighs(out(neighs)); 
                    
                    % Should really only take
                    %the neighbors that aren't in cells already
                    neighs = setdiff(neighs, cells); % FIXME is setdiff faster than a running find_cells a few more times?
                    idx = find_cells(G, neighs, interpFace);
                    neighs = neighs(idx);
                    %neighs
                    % for i = 1:numel(neighs) % FIXME
                    %     n = neighs(i);
                    %     if isempty(find(cells == n))
                    %         cells = [cells, n];
                    %         fprintf('append %d\n',n), pause
                    %     else
                    %         fprintf('already found %d\n', n), pause
                    %     end
                    % end
                    cells = [cells, neighs];
                    out(neighs) = 1;

                    %if numel(neighs)>0, keyboard, end
                end
            end
            
            % check
            out2 = find_cells(G, 1:G.cells.num, interpFace); 
            %disp(['sum ', num2str(sum(out2))])
            %keyboard
            %assert(sum(out2) == 0);
            if sum(out2) > 0
                keyboard
            end
            
            
            % flag = isConvex(G, 1:G.cells.num, interpFace);
            % while (flag)
            %     mycell = flag;
            %     theFaces = G.cells.faces(G.cells.facePos(mycell):G.cells.facePos(mycell+1)-1);
            %     neighbors = G.faces.neighbors(theFaces, :);
            %     neighbors = unique(neighbors(:));
            %     neighbors(neighbors == 0) = [];
            %     while (flag)
            %         fprintf("%d %d\n", mycell, flag);
            %         d = interpFace.coords(theFaces, :) - G.faces.centroids(theFaces, :);
            %         d = sqrt(dot(d, d, 2));
            %         [maxRatio, ind] = max(d./R(theFaces));
            %         y_sigma = HAP(theFaces(ind), :)';
            %         interpFace = correctHAP_local(G, theFaces(ind), interpFace, y_sigma, 0.9*maxRatio);
            %         flag = isConvex(G, mycell, interpFace);
            %     end
            %     flag = isConvex(G, neighbors(1):G.cells.num, interpFace);
            % end

        end
    elseif (nargin == 3)
        if (G.griddim == 2)
            R = 0.5 * G.faces.areas;
        else
            R = sqrt(G.faces.areas./pi);
        end
        R = R * myRatio;
        xf = G.faces.centroids;
        hap = interpFace.coords;
        d = hap - xf;
        d = sqrt(dot(d, d, 2));
        ind = find(d > R);
        for i = 1:numel(ind)
            interpFace = correctHAP_local(G, ind(i), interpFace, HAP(ind(i), :)', myRatio);
        end
    else
        error('Wrong number of inputs')
    end
end

function out = find_cells(G, cells, interpFace)

    in = zeros(numel(cells), 1);
    for i = 1:numel(cells)
        c = cells(i);
        faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
        hap = interpFace.coords(faces, :);
        ind = convhull(hap);
        in(i) = mex_inhull(G.cells.centroids(c, :), hap, ind, -1e-5);
    end
    out = ~in;
    %if numel(cells) > 1, keyboard, end
end

function interpFace = correct(G, cells, interpFace, R, HAP)

    for c = cells % Row vector
        faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
        d = interpFace.coords(faces, :) - G.faces.centroids(faces, :);
        d = vecnorm(d, 2, 2);
        [maxRatio, ind] = max(d./R(faces));
        y_sigma = HAP(faces(ind), :)'; % FIXME 
        interpFace = correctHAP_local(G, faces(ind), interpFace, y_sigma, 0.9*maxRatio);
    end

end



function flag = isConvex(G, mycells, interpFace)
    switch G.griddim
        case 2
            flag = 0;
            for i_cell = 1:numel(mycells)
                thecell = mycells(i_cell);
                xc = G.cells.centroids(thecell, 1);
                yc = G.cells.centroids(thecell, 2);
                theFaces = G.cells.faces(G.cells.facePos(thecell):G.cells.facePos(thecell+1)-1);
                hap = interpFace.coords(theFaces, :);
                ind = convhull(hap);
                xv = hap(ind, 1);
                yv = hap(ind, 2);
                in = inpolygon(xc, yc, xv, yv);
                if (~in)
                    flag = thecell;
                    break;
                end
            end
        case 3
            flag = 0;
            for i_cell = 1:numel(mycells)
                thecell = mycells(i_cell);
                %xc = G.cells.centroids(thecell, :);
                theFaces = G.cells.faces(G.cells.facePos(thecell):G.cells.facePos(thecell+1)-1);
                hap = interpFace.coords(theFaces, :);
                ind = convhull(hap);
                %in = inhull(xc, hap, ind, -1e-5);

                % % Test inpolyhedron: much slower
                % t1 = tic;
                % in2 = inpolyhedron(ind, hap, xc);
                % time_inpoly = toc(t1);
                % [time_inhull, time_inpoly]
                % if in ~= in2
                %     keyboard
                % end

                % Test mex inhull
                in = mex_inhull(G.cells.centroids(thecell, :), hap, ind, -1e-5);
                %in2,in
                %assert(all(in2==in));


                if (~in), flag = thecell;
                    break;
                end
            end
    end
end

function interpFace = correctHAP_local(G, i_face, interpFace, y_sigma, myRatio)
% Correct harmonic averaging point for i_face based on given myRatio
    if (myRatio > 0)
        if (G.griddim == 2)
            R = 0.5 * G.faces.areas(i_face) * myRatio;
        elseif (G.griddim == 3)
            R = myRatio * sqrt(G.faces.areas(i_face)/pi);
        end
        xm = G.faces.centroids(i_face, :)';
        interpFace.coords(i_face, :) = (xm + R * (y_sigma - xm) / norm(y_sigma-xm))';
    else
        interpFace.coords(i_face, :) = G.faces.centroids(i_face, :);
    end
end
