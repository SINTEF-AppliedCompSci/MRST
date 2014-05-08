function writeMRSTWells(W, fn)
%Serialise an MRST Well Structure to a file.
%
% SYNOPSIS:
%   writeMRSTWell(W, fn)
%
% PARAMETERS:
%   W  - Well structure.
%
%   fn - File name.  This will be passed to function FOPEN using mode 'wt'.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   writeDeck, fopen.

%{
#COPYRIGHT#
%}

% $Date: 2012-10-10 20:47:34 +0200 (Wed, 10 Oct 2012) $
% $Revision: 10061 $

   [fid, msg] = fopen(fn, 'wt');

   if fid < 0,
      error('Failed to open ''%s'': %s.', fn, msg);
   end
   
   fprintf(fid, '%d %d %d\n', ...
           numel(W(1).compi), numel(W), numel(vertcat(W.cells)));

   for i = 1 : numel(W),
      fprintf(fid, '%s '   , W(i).name);
      fprintf(fid, '%f '   , W(i).sign);
      fprintf(fid, '%d '   , W(i).refDepth);
      fprintf(fid, '%d '   , numel(W(i).cells));
      fprintf(fid, '%f '   , reshape(W(i).compi, 1, []));
      fprintf(fid, '%d '   , reshape(W(i).cells, 1, []) - 1);
      fprintf(fid, '%.18e ', reshape(W(i).WI   , 1, []));
      fprintf(fid, '%s '   , W(i).type);
      fprintf(fid, '%.18e ', W(i).val);
      fprintf(fid, '\n');
   end

   fclose(fid);
end
