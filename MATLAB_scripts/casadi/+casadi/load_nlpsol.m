function varargout = load_nlpsol(varargin)
    %LOAD_NLPSOL [INTERNAL] 
    %
    %  LOAD_NLPSOL(char name)
    %
    %Explicitly load a plugin dynamically.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlpsol.hpp#L38
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlpsol.cpp#L38-L40
    %
    %
  [varargout{1:nargout}] = casadiMEX(867, varargin{:});
end
