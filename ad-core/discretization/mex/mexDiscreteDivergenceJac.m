function varargout = mexDiscreteDivergenceJac(varargin)
   filename = 'mexDiscreteDivergenceJac.cpp';
   INCLUDE = {};

   OPTS = { '-O' };

   SRC = {filename};

   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   [varargout{1:nargout}] = mexDiscreteDivergenceJac(varargin{:});
end