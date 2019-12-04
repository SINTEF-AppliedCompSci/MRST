function t = makeConnListFromMat(nd, varargin)
% Make the connectivity list from the node distribution matrix for a
% structured grid.
% The node distribution matrix:
%
%       1  2  3                      nnx
%  1    *  *  *  *  *  *  ....  *  *  *
%  2    *  *  *  *  *  *  ....  *  *  *
%  3    *  *  *  *  *  *  ....  *  *  *
%       ......                   ......
%       *  *  *  *  *  *  ....  *  *  *
%  nny  *  *  *  *  *  *  ....  *  *  *
%
% The nodes corresponding to cell (i,j) is: 
%  {L(i,j), L(i+1,j), L(i+1,j+1), L(i,j+1)}
%
% SYNOPSIS:
%   t = makeConnListFromMat(nd)
%   t = makeConnListFromMat(nd, 'order', 'column')
%
% PARAMETERS:
%   nd    -   Node distribution matrix
%
% KEYWORD ARGUMENTS:
%   'order' - 'rows' (Default) | 'column': the picking order. The numbering 
%              of the connectivity list cycles along 'order' fastest.
%
% RETURNS:
%    t  - Connectivity list, n x 1 cell
%
% EXAMPLE:
%    [nnx, nny] = deal(10, 6);
%    x0 = linspace(-1, 1, nnx);
%    x  = arrayfun(@(j)j*x0, (1:nny)', 'UniformOutput', false);
%    x  = cell2mat(x);
%    y0 = linspace(-5, 5, nny)';
%    y  = repmat(y0, 1, nnx);
%    p  = [x(:), y(:)];
%    nd = reshape((1:size(p,1)), nny, nnx);
%    t  = makeConnListFromMat(nd);
%    G  = tessellationGrid(p, t);
%    figure, plotCellData(G, (1:G.cells.num)')
%    t  = makeConnListFromMat(nd, 'order', 'column');
%    G  = tessellationGrid(p, t);
%    figure, plotCellData(G, (1:G.cells.num)')
%
% SEE ALSO:
%  `buildRadialGrid`, `tessellationGrid`

    opt = struct('order', 'rows');
    opt = merge_options(opt, varargin{:});

    switch opt.order
        case 'rows'
        case 'column'
            nd = nd';
        otherwise
            error('Unknown order')
    end
    
    ny = size(nd, 1)-1;
    nx = size(nd, 2)-1;
    i = repmat((1:ny)',  1, nx);
    j = repmat((1:nx),  ny,  1);
    t = arrayfun(@(i, j)[nd(i,j), nd(i+1,j), nd(i+1,j+1), nd(i,j+1)], ...
        i(:), j(:), 'UniformOutput', false);
  
%     nny = size(nd, 1);
%     nnx = size(nd, 2);
%     t = zeros(0,4);
%     for j = 1 : nnx-1
%         nd0 = nd(:, [j, j+1]);
%         t0 = zeros(nny-1, 4);
%         t0(:, 1) = nd0(1:nny-1, 1);
%         t0(:, 2) = nd0(2:nny  , 1);
%         t0(:, 3) = nd0(2:nny  , 2);
%         t0(:, 4) = nd0(1:nny-1, 2);
%         t = [t; t0];
%     end
end