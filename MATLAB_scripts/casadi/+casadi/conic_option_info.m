function varargout = conic_option_info(varargin)
    %CONIC_OPTION_INFO [INTERNAL] 
    %
    %  char = CONIC_OPTION_INFO(char name, char op)
    %
    %Get documentation for a particular option.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1em
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/conic.hpp#L572
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/conic.cpp#L572-L574
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(852, varargin{:});
end
