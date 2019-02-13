function varargout = mexFaceAverageDiagonalJac(varargin)
   filename = 'mexFaceAverageDiagonalJac.cpp';
   INCLUDE = {};

   OPTS = { '-O' };

   SRC = {filename};

   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   [varargout{1:nargout}] = mexFaceAverageDiagonalJac(varargin{:});
end
