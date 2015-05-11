function varargout = interp1q_mex(varargin)

% Get path of this file
p = fileparts(mfilename('fullpath'));

% Build compile string
mexcmd = sprintf('mex -O CFLAGS="\\$CFLAGS -std=c99" %s %s', ...
    sprintf('-outdir "%s"', p), ...
    sprintf('"%s"', fullfile(p,'interp1q_mex.c')) );

% Run compile string
eval(mexcmd);

% Call MEX'ed edition.
[varargout{1:nargout}] = interp1q_mex(varargin{:});

end
