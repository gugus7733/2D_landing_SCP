function varargout = extract(varargin)
    %EXTRACT Introduce intermediate variables for selected nodes in a graph.
    %
    %  [{DM} OUTPUT1, {DM} OUTPUT2, {DM} OUTPUT3] = EXTRACT({DM} ex, struct opts)
    %  [{SX} OUTPUT1, {SX} OUTPUT2, {SX} OUTPUT3] = EXTRACT({SX} ex, struct opts)
    %  [{MX} OUTPUT1, {MX} OUTPUT2, {MX} OUTPUT3] = EXTRACT({MX} ex, struct opts)
    %
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1d5
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1138
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1138-L1143
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(928, varargin{:});
end
