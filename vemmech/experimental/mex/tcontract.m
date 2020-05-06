function varargout = tcontract(varargin)

   filename = 'tcontract.cpp';
   INCLUDE = {};
   OPTS= {'-O'};
   SRC = {filename};
   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();
   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   
   [varargout{1:nargout}] = tcontract(varargin{:});
end

