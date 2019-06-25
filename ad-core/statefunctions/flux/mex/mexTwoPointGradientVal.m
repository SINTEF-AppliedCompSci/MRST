function varargout = mexTwoPointGradientVal(varargin)
   filename = 'mexTwoPointGradientVal.cpp';
   INCLUDE = {};

   OPTS = { '-O' };

   SRC = {filename};

   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   [varargout{1:nargout}] = mexTwoPointGradientVal(varargin{:});
end