function res = concatStruct(str1, str2)
% Concatenate two structures into one.  The structures should not have any
% overlapping field names.

  fnames = [fieldnames(str1);  fieldnames(str2)];
  values = [struct2cell(str1); struct2cell(str2)];
  res = cell2struct(values, fnames, 1);

end
