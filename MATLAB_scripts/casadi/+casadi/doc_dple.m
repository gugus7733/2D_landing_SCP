function varargout = doc_dple(varargin)
    %DOC_DPLE [INTERNAL] 
    %
    %  char = DOC_DPLE(char name)
    %
    %Get the documentation string for a plugin.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/dple.hpp#L39
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/dple.cpp#L39-L41
    %
    %
  [varargout{1:nargout}] = casadiMEX(904, varargin{:});
end
