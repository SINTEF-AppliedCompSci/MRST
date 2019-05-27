function varargout = mexTwoPointGradientDiagonalJac(varargin)
   filename = 'mexTwoPointGradientDiagonalJac.cpp';
   INCLUDE = {};

   OPTS = { '-O' };

   SRC = {filename};

   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   [varargout{1:nargout}] = mexTwoPointGradientDiagonalJac(varargin{:});
end
