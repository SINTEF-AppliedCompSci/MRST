function varargout = mexDiagonalSparse(varargin)
   filename = 'mexDiagonalSparse.cpp';
   INCLUDE = {};

   OPTS = { '-O' };

   SRC = {filename};

   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   [varargout{1:nargout}] = mexDiagonalSparse(varargin{:});
end
