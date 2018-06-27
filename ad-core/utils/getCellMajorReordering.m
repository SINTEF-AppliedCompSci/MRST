function ordering = getCellMajorReordering(ncell, block_size, varargin)
%Get equation ordering transforming variable major to cell major ordering
%
% SYNOPSIS:
%   ordering = getCellMajorReordering(ncell, block_size)
%
% DESCRIPTION:
%   Get a permutation vector which transforms a system on the standard
%   variable major (e.g., a two-phase system p_1, ..., p_n, s_1, ..., s_n)
%   into a cell major (e.g. p_1, s_1, ..., p_n, s_n) where p and s are two
%   primary variables and the subscript refers to a specific cell.
%
%   If Ax=b is some system to be re-ordered of size ncell*block_size, then
%   ordering = getCellMajorReordering(ncell, block_size);
%   A = A(ordering, ordering);
%   b = b(ordering);
%   will re-order the system. Solving the system x = solve(A, b) where
%   solve is some linaer solver will then give a permuted solution to the
%   system. The final solution is then x(ordering) = x.
%
%   The primary utility of this function is to a) Allow the user to call
%   external linear solvers which require this type of ordering and b)
%   change the system ordering for e.g. ILU(0).
%
% REQUIRED PARAMETERS:
%   ncell      - Number of cells in grid to be used.
%
%   block_size - Size of each block (e.g. number of components present)
%
% OPTIONAL PARAMETERS:
%   ndof              - Total number of degrees of freedom. Any additional
%                       degrees of freedom beyond the ncell*block_size
%                       first variables will have a identity remapping,
%                       retaining the position in the final system.
%
%   equation_ordering - An optional ordering of the equations. Should be a
%                       block_size length vector.
% 
%   cell_ordering     - An optional ordering of the cells themselves.
%
% RETURNS:
%   ordering          - Ordering vector so that A = A(ordering, ordering)
%                       is the permuted system
%
% SEE ALSO:
%   `LinearSolverAD`, `AMGCL_CPRSolverAD`


    ncell_total = ncell*block_size;

    opt = struct('ndof', ncell_total, ...
                 'equation_ordering', [], ...
                 'cell_ordering', []);
    opt = merge_options(opt, varargin{:});

    ordering = (1:opt.ndof)';
    
    subs = 1:ncell_total;
    subs = reshape(subs, [], block_size)';
    if ~isempty(opt.cell_ordering)
        subs = subs(:, opt.cell_ordering);
    end
    if ~isempty(opt.equation_ordering)
        subs = subs(opt.equation_ordering, :);
    end
    subs = subs(:);
    ordering(1:ncell_total) = subs;
end