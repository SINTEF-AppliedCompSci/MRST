function [Pts, wellType, removed] = removeConflictPoints(Pts, gridSpacing, ...
                                                         priIndex,varargin)
    % Remove conflict points
    % Arguments:
    %   Pts             2*n array of points
    %   gridSpacing     array of length n containing the minimum allowed
    %                   distance to other points
    %   priIndex        array of length n containing the priority of each
    %                   point
    %   wellType        OPTIONAL
    %                   preset false. Logical array of length n that is 
    %                   true for well points
    %   
    % Return:
    %   Pts             array of length <=n of points that were not removed
    %   wellType        Logical array that is true for well points
    %   removed         logical array of length n that is true for points
    %                   that were removed
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt = struct('wellType',false(size(Pts,1),1));
    opt = merge_options(opt,varargin{:});
    wellType = opt.wellType;
    
    gridSpacing = gridSpacing*(1-1e-10); % To avoid floating point errors
    
    Ic = 1:size(Pts, 1);                % Sorting map
    ptsToClose = Pts;                   % Mark all points as potential 
    removed = false(size(Pts, 1), 1);   % conflict points
    
    distance = pdist(ptsToClose)';      % Calculate distance, all to all
    dlt = distLessThan(distance, gridSpacing(Ic)); % find conflict points
    Ic = findToClose(dlt);                         % Map back to array
    while length(Ic)>1
        [~, Ii ] = sort(priIndex(Ic), 'descend');  % Sort after priority

        removePoint = Ic(Ii(1));        % Remove least important point
        if wellType(removePoint)        % If it is a well point, mark the
            n = size(Ic,1);             % corresponding cell as a well cell
            if Ii(1)>n/2;
                wellType(Ic(ceil(mod(Ii(1),(n+1)/2)))) = true;
            else
                wellType(Ic(Ii(1)+n/2)) = true;
            end
        end
        removed(removePoint) = true;    
        Ic = Ic(Ic~=removePoint);       % Update sorting maps
        Ic = unique(Ic);
        ptsToClose = Pts(Ic,:);
        
        if size(ptsToClose,1) ==1
            continue
        end
        distance = pdist(ptsToClose)';  % Calculate distance of potential 
        dlt = distLessThan(distance, gridSpacing(Ic)); % conflict points
        Ic = Ic(findToClose(dlt));
        
    end
    Pts      = Pts(~removed,:);         % Update return arrays.
    wellType = wellType(~removed, :);
end


function [arr] = distLessThan(distance, b)
    n = length(distance);
    [i,j] = arrToMat(1:n, n);
    arr = distance < max(b(i), b(j));
end


function [indexes] = findToClose(arr)
    n = length(arr);
    k = find(arr);
    [i, j] = arrToMat(k, n);
    indexes = [i;j];
end



function [i, j] = arrToMat(k, n)
    m = ceil(sqrt(2*n));
    j = ceil((2*m-1)/2 - 0.5*sqrt((2*m-1)^2 - 8*k));
    i = k + j - (j-1).*(m-j/2);
end

% function [k] = matToArr(i,j, m)
%     assert(all(abs(j)) && all(abs(i-j)) && all(abs(1+m-i)));
%     k = 1 + (j-1)*m - (j-1).*j/2 + i-j - 1;
% end
