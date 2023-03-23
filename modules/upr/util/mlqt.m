function [res] = mlqt(cellCenter, bndr, cellSize, varargin)
    % Function to perform multilevel grid refinement
    % Arguments:
    %   cellCenter      Center coordinate of cell
    %   bndr            Array of coordinates where cell should be refined
    %   cellSize        The size of the cell
    
    % varargin:
    %   level           Current refinement level
    %   maxLevel        Max refinement level
    %   distTol         If a cell is closer to bndr than distTol, it is
    %                   refined
    
    % Output:
    %   res             Struct containg the new cell centers, and the cell
    %                   sizes

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
    
    % Load options
    opt = struct('level', 1, ...
                 'maxLev', 0, ...
                 'distTol', -1);
             
    opt = merge_options(opt, varargin{:});
    level = opt.level;
    maxLev = opt.maxLev;
    
    
    % Test input
    assert(all(cellSize>0));
    assert(size(opt.distTol,1) >0)
    
    
    % Is recursion finished? 
    if level> maxLev
        res = {cellCenter, cellSize};
        return
    end
    
    % Set distance tolerance
    if size(opt.distTol,1)==1
        if opt.distTol <= 0
            distTol = 2.5*cellSize/2;
        else
            distTol = opt.distTol;
        end
        distNext = distTol/2;
    else
        distTol = opt.distTol(level);
        distNext = opt.distTol;
    end
    
    n = size(bndr,1);
    repPnt = repmat(cellCenter,n,1);
    if any(all(bsxfun(@le, abs(repPnt-bndr), distTol),2)) % Should cell be refined?
        dx = cellSize/4;
        varArg = {'level', level+1, 'maxLev', maxLev, 'distTol', distNext};
        shift = [dx(1),  dx(2);
                 dx(1), -dx(2);
                -dx(1), -dx(2);
                -dx(1),  dx(2)];
        if size(cellCenter, 2)==3
            shift = [shift,  dx(3) * ones(4,1);
                     shift, -dx(3) * ones(4,1)];
        end
        res = [];
        for i = 1:size(shift, 1)
            res = [res; mlqt(cellCenter + shift(i, :), bndr, cellSize/2, varArg{:})];
        end
    else
        res = {cellCenter, cellSize};
    end
end
