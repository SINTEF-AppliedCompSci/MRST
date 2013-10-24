function res = interleave(cellarray1, cellarray2)
% Interleave the elements of two cell arrays
    
  res = cell(2, max(numel(cellarray1), numel(cellarray2)));
  res(1, 1:numel(cellarray1)) = cellarray1;
  res(2, 1:numel(cellarray2)) = cellarray2;
  res = res(:);
end
