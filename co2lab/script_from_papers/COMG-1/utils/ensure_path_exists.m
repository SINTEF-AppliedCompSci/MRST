function ensure_path_exists(filename)
   
% If filename is located in a folder other than the present one, ensure that
% folder exists.
   separator = find(filename == '/', 1, 'last');
   if ~isempty(separator)
      path = filename(1:separator);
      if exist(path) ~= 7
         mkdir(filename(1:separator));
      end
   end
end
