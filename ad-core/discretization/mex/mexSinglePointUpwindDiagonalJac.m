function varargout = mexSinglePointUpwindDiagonalJac(varargin)
   filename = 'mexSinglePointUpwindDiagonalJac.cpp';
   INCLUDE = {};

   OPTS = { '-O' };

   SRC = {filename};

   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   [varargout{1:nargout}] = mexSinglePointUpwindDiagonalJac(varargin{:});
end

