function varargout = mexSinglePointUpwindVal(varargin)
   filename = 'mexSinglePointUpwindVal.cpp';
   INCLUDE = {};

   OPTS = { '-O' };

   SRC = {filename};

   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   [varargout{1:nargout}] = mexSinglePointUpwindVal(varargin{:});
end
