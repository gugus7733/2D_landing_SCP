function varargout = nlpsol_n_out(varargin)
    %NLPSOL_N_OUT [INTERNAL] 
    %
    %  int = NLPSOL_N_OUT()
    %
    %Number of NLP solver outputs.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1t3
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlpsol.hpp#L263
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlpsol.cpp#L263-L265
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(861, varargin{:});
end
