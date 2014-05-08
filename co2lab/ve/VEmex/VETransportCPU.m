function varargout = VETransportCPU(varargin)  
% build/link to mex file

d = fileparts(mfilename('fullpath'));
mexfile = ['VETransportCPU.', mexext];
mexfullfile = fullfile(d, 'build', 'mex', mexfile);
if exist(mexfullfile, 'file')
   tofile  = fullfile(d, mexfile);
   system(['ln -s ', mexfullfile, ' ', tofile]);
   [varargout{1:nargout}] = VETransportCPU(varargin{:});
else
   error(['You need to build VEmex manually. \n',...
         'For instructions see:' , fullfile(d, 'README.')]);
end
