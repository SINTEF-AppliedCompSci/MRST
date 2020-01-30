function C = insertVec(A,B,id)
% Insert vector B into vector A at possitions id.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

[id,I] = sort(id);
B      = B(I,:);
C = zeros(size(A)+[size(B,1),0])+nan;
C(id + (0:numel(id)-1)',:) = B;
C(isnan(C)) = A;
end
