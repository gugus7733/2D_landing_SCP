function varargout = conic_n_out(varargin)
    %CONIC_N_OUT [INTERNAL] 
    %
    %  int = CONIC_N_OUT()
    %
    %Get the number of QP solver outputs.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ej
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/conic.hpp#L107
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/conic.cpp#L107-L109
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(849, varargin{:});
end
