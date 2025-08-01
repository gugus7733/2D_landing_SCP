classdef  OptiAdvanced < casadi.Opti
    %OPTIADVANCED [INTERNAL] C++ includes: optistack.hpp
    %
    %
    %
    %
  methods
    function delete(self)
        if self.swigPtr
          casadiMEX(1316, self);
          self.SwigClear();
        end
    end
    function varargout = solver(self,varargin)
    %SOLVER [INTERNAL] 
    %
    %  Function = SOLVER(self)
    %
    %Get the underlying CasADi solver of the  Opti stack.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L564
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L500-L506
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1317, self, varargin{:});
    end
    function varargout = is_parametric(self,varargin)
    %IS_PARAMETRIC [INTERNAL] 
    %
    %  bool = IS_PARAMETRIC(self, MX expr)
    %
    %return true if expression is only dependant on  Opti parameters,
    % not variables
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L567
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L508-L514
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1318, self, varargin{:});
    end
    function varargout = symvar(self,varargin)
    %SYMVAR [INTERNAL] 
    %
    %  {MX} = SYMVAR(self)
    %  {MX} = SYMVAR(self, MX expr)
    %  {MX} = SYMVAR(self, MX expr, casadi::VariableType type)
    %
    %Get symbols present in expression.
    %
    %Returned vector is ordered according to the order of  variable()/parameter()
    % calls used to create the variables
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1u
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L578
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L532-L538
    %
    %
    %
    %.......
    %
    %::
    %
    %  SYMVAR(self)
    %
    %
    %
    %[INTERNAL] 
    %Get symbols present in expression.
    %
    %Returned vector is ordered according to the order of  variable()/parameter()
    % calls used to create the variables
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1u
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L576
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L516-L522
    %
    %
    %
    %.............
    %
    %
    %.......
    %
    %::
    %
    %  SYMVAR(self, MX expr)
    %
    %
    %
    %[INTERNAL] 
    %Get symbols present in expression.
    %
    %Returned vector is ordered according to the order of  variable()/parameter()
    % calls used to create the variables
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1u
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L577
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L524-L530
    %
    %
    %
    %.............
    %
    %
    %.......
    %
    %::
    %
    %  SYMVAR(self, MX expr, casadi::VariableType type)
    %
    %
    %
    %[INTERNAL] 
    %Get symbols present in expression.
    %
    %Returned vector is ordered according to the order of  variable()/parameter()
    % calls used to create the variables
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1u
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L578
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L532-L538
    %
    %
    %
    %.............
    %
    %
      [varargout{1:nargout}] = casadiMEX(1319, self, varargin{:});
    end
    function varargout = canon_expr(self,varargin)
    %CANON_EXPR [INTERNAL] 
    %
    %  MetaCon = CANON_EXPR(self, MX expr)
    %
    %Interpret an expression (for internal use only)
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L582
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L540-L546
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1320, self, varargin{:});
    end
    function varargout = get_meta(self,varargin)
    %GET_META [INTERNAL] 
    %
    %  MetaVar = GET_META(self, MX m)
    %
    %Get meta-data of symbol (for internal use only)
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L585
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L548-L554
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1321, self, varargin{:});
    end
    function varargout = get_meta_con(self,varargin)
    %GET_META_CON [INTERNAL] 
    %
    %  MetaCon = GET_META_CON(self, MX m)
    %
    %Get meta-data of symbol (for internal use only)
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L588
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L556-L562
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1322, self, varargin{:});
    end
    function varargout = set_meta(self,varargin)
    %SET_META [INTERNAL] 
    %
    %  SET_META(self, MX m, MetaVar meta)
    %
    %Set meta-data of an expression.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L591
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L564-L570
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1323, self, varargin{:});
    end
    function varargout = set_meta_con(self,varargin)
    %SET_META_CON [INTERNAL] 
    %
    %  SET_META_CON(self, MX m, MetaCon meta)
    %
    %Set meta-data of an expression.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L594
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L572-L578
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1324, self, varargin{:});
    end
    function varargout = assert_active_symbol(self,varargin)
    %ASSERT_ACTIVE_SYMBOL [INTERNAL] 
    %
    %  ASSERT_ACTIVE_SYMBOL(self, MX m)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1325, self, varargin{:});
    end
    function varargout = active_symvar(self,varargin)
    %ACTIVE_SYMVAR [INTERNAL] 
    %
    %  {MX} = ACTIVE_SYMVAR(self, casadi::VariableType type)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1326, self, varargin{:});
    end
    function varargout = active_values(self,varargin)
    %ACTIVE_VALUES [INTERNAL] 
    %
    %  {DM} = ACTIVE_VALUES(self, casadi::VariableType type)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1327, self, varargin{:});
    end
    function varargout = x_lookup(self,varargin)
    %X_LOOKUP [INTERNAL] 
    %
    %  MX = X_LOOKUP(self, index i)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1328, self, varargin{:});
    end
    function varargout = g_lookup(self,varargin)
    %G_LOOKUP [INTERNAL] 
    %
    %  MX = G_LOOKUP(self, index i)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1329, self, varargin{:});
    end
    function varargout = g_index_reduce_g(self,varargin)
    %G_INDEX_REDUCE_G [INTERNAL] 
    %
    %  int = G_INDEX_REDUCE_G(self, index i)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1330, self, varargin{:});
    end
    function varargout = g_index_reduce_x(self,varargin)
    %G_INDEX_REDUCE_X [INTERNAL] 
    %
    %  int = G_INDEX_REDUCE_X(self, index i)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1331, self, varargin{:});
    end
    function varargout = g_index_unreduce_g(self,varargin)
    %G_INDEX_UNREDUCE_G [INTERNAL] 
    %
    %  int = G_INDEX_UNREDUCE_G(self, index i)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1332, self, varargin{:});
    end
    function varargout = x_describe(self,varargin)
    %X_DESCRIBE [INTERNAL] 
    %
    %  char = X_DESCRIBE(self, index i, struct opts)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1333, self, varargin{:});
    end
    function varargout = g_describe(self,varargin)
    %G_DESCRIBE [INTERNAL] 
    %
    %  char = G_DESCRIBE(self, index i, struct opts)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1334, self, varargin{:});
    end
    function varargout = describe(self,varargin)
    %DESCRIBE [INTERNAL] 
    %
    %  char = DESCRIBE(self, MX x, index indent, struct opts)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1335, self, varargin{:});
    end
    function varargout = show_infeasibilities(self,varargin)
    %SHOW_INFEASIBILITIES [INTERNAL] 
    %
    %  SHOW_INFEASIBILITIES(self, double tol, struct opts)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1336, self, varargin{:});
    end
    function varargout = solve_prepare(self,varargin)
    %SOLVE_PREPARE [INTERNAL] 
    %
    %  SOLVE_PREPARE(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1337, self, varargin{:});
    end
    function varargout = solve_actual(self,varargin)
    %SOLVE_ACTUAL [INTERNAL] 
    %
    %  struct:DM = SOLVE_ACTUAL(self, struct:DM args)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1338, self, varargin{:});
    end
    function varargout = arg(self,varargin)
    %ARG [INTERNAL] 
    %
    %  struct:DM = ARG(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1339, self, varargin{:});
    end
    function varargout = res(self,varargin)
    %RES [INTERNAL] 
    %
    %  struct:DM = RES(self)
    %  RES(self, struct:DM res)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1340, self, varargin{:});
    end
    function varargout = constraints(self,varargin)
    %CONSTRAINTS [INTERNAL] 
    %
    %  {MX} = CONSTRAINTS(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1341, self, varargin{:});
    end
    function varargout = objective(self,varargin)
    %OBJECTIVE [INTERNAL] 
    %
    %  MX = OBJECTIVE(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1342, self, varargin{:});
    end
    function varargout = baked_copy(self,varargin)
    %BAKED_COPY [INTERNAL] 
    %
    %  OptiAdvanced = BAKED_COPY(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1343, self, varargin{:});
    end
    function varargout = assert_empty(self,varargin)
    %ASSERT_EMPTY [INTERNAL] 
    %
    %  ASSERT_EMPTY(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1344, self, varargin{:});
    end
    function varargout = bake(self,varargin)
    %BAKE [INTERNAL] 
    %
    %  BAKE(self)
    %
    %Fix the structure of the optimization problem.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.hpp#L629
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/optistack.cpp#L778-L784
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1345, self, varargin{:});
    end
    function v = problem_dirty_(self)
      v = casadiMEX(1346, self);
    end
    function varargout = mark_problem_dirty(self,varargin)
    %MARK_PROBLEM_DIRTY [INTERNAL] 
    %
    %  MARK_PROBLEM_DIRTY(self, bool flag)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1347, self, varargin{:});
    end
    function varargout = problem_dirty(self,varargin)
    %PROBLEM_DIRTY [INTERNAL] 
    %
    %  bool = PROBLEM_DIRTY(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1348, self, varargin{:});
    end
    function v = solver_dirty_(self)
      v = casadiMEX(1349, self);
    end
    function varargout = mark_solver_dirty(self,varargin)
    %MARK_SOLVER_DIRTY [INTERNAL] 
    %
    %  MARK_SOLVER_DIRTY(self, bool flag)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1350, self, varargin{:});
    end
    function varargout = solver_dirty(self,varargin)
    %SOLVER_DIRTY [INTERNAL] 
    %
    %  bool = SOLVER_DIRTY(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1351, self, varargin{:});
    end
    function v = solved_(self)
      v = casadiMEX(1352, self);
    end
    function varargout = mark_solved(self,varargin)
    %MARK_SOLVED [INTERNAL] 
    %
    %  MARK_SOLVED(self, bool flag)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1353, self, varargin{:});
    end
    function varargout = solved(self,varargin)
    %SOLVED [INTERNAL] 
    %
    %  bool = SOLVED(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1354, self, varargin{:});
    end
    function varargout = assert_solved(self,varargin)
    %ASSERT_SOLVED [INTERNAL] 
    %
    %  ASSERT_SOLVED(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1355, self, varargin{:});
    end
    function varargout = assert_baked(self,varargin)
    %ASSERT_BAKED [INTERNAL] 
    %
    %  ASSERT_BAKED(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1356, self, varargin{:});
    end
    function varargout = instance_number(self,varargin)
    %INSTANCE_NUMBER [INTERNAL] 
    %
    %  int = INSTANCE_NUMBER(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1357, self, varargin{:});
    end
    function self = OptiAdvanced(varargin)
    %OPTIADVANCED 
    %
    %  new_obj = OPTIADVANCED(Opti x)
    %
    %
      self@casadi.Opti(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = casadiMEX(1358, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
  end
  methods(Static)
  end
end
