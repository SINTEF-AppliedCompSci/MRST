function varargout = mex_nakeinterp1(varargin)

disp('Compiling mex_nakeinterp1...');

% Compile
mex -O -v CFLAGS="\$CFLAGS -std=c99" mex_nakeinterp1.c

% Call MEX'ed edition.
[varargout{1:nargout}] = mex_nakeinterp1(varargin{:});

end
