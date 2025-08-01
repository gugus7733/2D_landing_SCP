function varargout = has_linsol(varargin)
    %HAS_LINSOL [INTERNAL] 
    %
    %  bool = HAS_LINSOL(char name)
    %
    %Check if a particular plugin is available.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/linsol.hpp#L206
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/linsol.cpp#L206-L208
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(894, varargin{:});
end
