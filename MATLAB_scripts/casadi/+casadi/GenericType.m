classdef  GenericType < casadi.PrintableCommon
    %GENERICTYPE [INTERNAL] 
    %
    %
    %Generic data type, can hold different types such as bool, 
    %casadi_int, 
    %std::string etc.
    %
    %Joel Andersson
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_17n
    %
    %C++ includes: generic_type.hpp
    %
    %
  methods
    function varargout = serialize(self,varargin)
    %SERIALIZE [INTERNAL] 
    %
    %  SERIALIZE(self, casadi::SerializingStream & s)
    %
    %Serialize an object.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_17r
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_type.hpp#L249
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_type.cpp#L742-L745
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(48, self, varargin{:});
    end
    function self = GenericType(varargin)
    %GENERICTYPE 
    %
    %  new_obj = GENERICTYPE()
    %
    %
      self@casadi.PrintableCommon(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = casadiMEX(50, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
        if self.swigPtr
          casadiMEX(51, self);
          self.SwigClear();
        end
    end
  end
  methods(Static)
    function varargout = deserialize(varargin)
    %DESERIALIZE 
    %
    %  GenericType = DESERIALIZE(casadi::DeserializingStream & s)
    %
    %
     [varargout{1:nargout}] = casadiMEX(49, varargin{:});
    end
  end
end
