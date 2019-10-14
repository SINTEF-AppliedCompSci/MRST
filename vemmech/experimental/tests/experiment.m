A = rand(6) > 0.8;
B = rand(6) > 0.8;
ABt = A*B';

[ai, aj] = ind2sub(size(A), find(A));
[bi, bj] = ind2sub(size(B), find(B));

aij = [nonzeros(A), ai, aj];
bij = [nonzeros(B), bi, bj];

common_j = intersect(aj, bj);
num_unique_j = numel(common_j);

% remove values that will not partake in multiplication
aij = aij(ismember(aj, common_j), :);
bij = bij(ismember(bj, common_j), :);

% number of repetitions of j-index for aij and bij

astart = [1; find(diff(aij(:,3)))+1];
aend   = [astart(2:end)-1; numel(aij(:,3))];

bstart = [1; find(diff(bij(:,3)))+1];
bend   = [bstart(2:end)-1; numel(bij(:,3))];

% arep = aend - astart + 1;
% brep = bend - bstart + 1;

% split up matrices
asplit = arrayfun(@(x) aij(astart(x):aend(x), :), 1:num_unique_j, ...
                  'uniformoutput', false);
bsplit = arrayfun(@(x) bij(bstart(x):bend(x), :), 1:num_unique_j, ...
                  'uniformoutput', false);

% function producing all permuations of rows in u and v
loc_kron = @(u, v) [repmat(u, size(v, 1), 1), repelem(v, size(u, 1), 1)];

tmp = arrayfun(@(ix) loc_kron(asplit{ix}, bsplit{ix}), 1:num_unique_j, ...
               'uniformoutput', false);
tmp = vertcat(tmp{:});

vals = tmp(:,1) .* tmp(:, 4);
ixs = tmp(:, [2,5]);

result = sparse(ixs(:,1), ixs(:,2), vals, size(A, 1), size(A, 2));

spy(result);
%full(result)


% ============================== test tsparsemul ==============================

A = double(rand(7) > 0.7);
B = (rand(7) > 0.5);
% A = rand(4);
% B = rand(4);


A(find(A)) = rand(nnz(A), 1);
B(find(B)) = rand(nnz(B), 1);

facit = A*B;

avals = nonzeros(A);%ones(nnz(A), 1);
bvals = nonzeros(B);%ones(nnz(B), 1);

[ai, aj] = ind2sub(size(A), find(A));
[bi, bj] = ind2sub(size(B), find(B));

aind = [ai, aj];
bind = [bj, bi];

[bind, ix] = sortrows(bind, 2);
bvals = bvals(ix); % the sorting order must also apply to bvals

[vals, ixs] = tsparsemul(avals, aind, bvals, bind);

res = full(sparse(ixs(:,1), ixs(:,2), vals, size(facit, 1), size(facit, 2)));


% A has entry (1,2) which is multiplied by B entry (2,3) to produce entry (1,3)