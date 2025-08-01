classdef  NlpBuilder < casadi.PrintableCommon
    %NLPBUILDER [INTERNAL] 
    %
    %
    %A symbolic NLP representation.
    %
    %Joel Andersson
    %
    %::
    %
    %  Extra doc: https://github.com/casadi/casadi/wiki/L_1e2 
    %  
    %
    %
    %
    %C++ includes: nlp_builder.hpp
    %
    %
  methods
    function v = x(self)
      v = casadiMEX(1020, self);
    end
    function v = f(self)
      v = casadiMEX(1021, self);
    end
    function v = g(self)
      v = casadiMEX(1022, self);
    end
    function v = x_lb(self)
      v = casadiMEX(1023, self);
    end
    function v = x_ub(self)
      v = casadiMEX(1024, self);
    end
    function v = g_lb(self)
      v = casadiMEX(1025, self);
    end
    function v = g_ub(self)
      v = casadiMEX(1026, self);
    end
    function v = x_init(self)
      v = casadiMEX(1027, self);
    end
    function v = lambda_init(self)
      v = casadiMEX(1028, self);
    end
    function v = discrete(self)
      v = casadiMEX(1029, self);
    end
    function varargout = import_nl(self,varargin)
    %IMPORT_NL [INTERNAL] 
    %
    %  IMPORT_NL(self, char filename, struct opts)
    %
    %Import an .nl file.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.hpp#L74
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.cpp#L32-L35
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1030, self, varargin{:});
    end
    function varargout = type_name(self,varargin)
    %TYPE_NAME [INTERNAL] 
    %
    %  char = TYPE_NAME(self)
    %
    %Readable name of the class.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.hpp#L77
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.hpp#L77-L77
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1031, self, varargin{:});
    end
    function varargout = disp(self,varargin)
    %DISP [INTERNAL] 
    %
    %  std::ostream & = DISP(self, bool more)
    %
    %Print a description of the object.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.hpp#L80
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.cpp#L37-L45
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1032, self, varargin{:});
    end
    function varargout = str(self,varargin)
    %STR [INTERNAL] 
    %
    %  char = STR(self, bool more)
    %
    %Get string representation.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.hpp#L83
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/nlp_builder.hpp#L83-L87
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1033, self, varargin{:});
    end
    function self = NlpBuilder(varargin)
    %NLPBUILDER 
    %
    %  new_obj = NLPBUILDER()
    %
    %
      self@casadi.PrintableCommon(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = casadiMEX(1034, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
        if self.swigPtr
          casadiMEX(1035, self);
          self.SwigClear();
        end
    end
  end
  methods(Static)
  end
end
