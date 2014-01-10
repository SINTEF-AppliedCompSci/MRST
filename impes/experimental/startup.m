% Find absolute path to startup
nm = mfilename('fullpath');
ix = strfind(nm, filesep);
if ~isempty(ix),
   root = nm(1 : ix(end));
else
   root = nm;
end

% Set path
run(fullfile(nm(1:ix(end-2)), 'branches','mrst-reorg','startup.m'))
%run(fullfile(nm(1:ix(end-2)), 'branches','mrst-releases','1.2','startup.m'))
%run(fullfile(nm(1:ix(end-2)), 'branches','mrst-releases','2011a','startup.m'))
mrstModule('add', root(1:end-1));
