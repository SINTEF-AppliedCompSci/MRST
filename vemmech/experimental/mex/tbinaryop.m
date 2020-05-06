function varargout = tbinaryop(varargin)

   filename = 'tbinaryop.cpp';
   INCLUDE = {};
   OPTS= {'-O'};
   SRC = {filename};
   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();
   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   
   [varargout{1:nargout}] = tbinaryop(varargin{:});
   
end
