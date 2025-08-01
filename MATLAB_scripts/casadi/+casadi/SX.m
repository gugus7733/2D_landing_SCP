classdef  SX < casadi.MatrixCommon & casadi.GenericExpressionCommon & casadi.GenSX & casadi.PrintableCommon
    %SX 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = casadiMEX(3, self);
    end
    function varargout = has_nz(self,varargin)
    %HAS_NZ [INTERNAL] 
    %
    %  bool = HAS_NZ(self, int rr, int cc)
    %
    %Returns true if the matrix has a non-zero at location rr, cc.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L226
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_impl.hpp#L65-L67
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(547, self, varargin{:});
    end
    function varargout = nonzero(self,varargin)
    %NONZERO [INTERNAL] 
    %
    %  bool = NONZERO(self)
    %
    %Returns the truth value of a  Matrix.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L229
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_impl.hpp#L70-L76
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(548, self, varargin{:});
    end
    function varargout = get(self,varargin)
    %GET 
    %
    %  SX = GET(self, bool ind1, Sparsity sp)
    %  SX = GET(self, bool ind1, Slice rr)
    %  SX = GET(self, bool ind1, IM rr)
    %  SX = GET(self, bool ind1, Slice rr, Slice cc)
    %  SX = GET(self, bool ind1, Slice rr, IM cc)
    %  SX = GET(self, bool ind1, IM rr, Slice cc)
    %  SX = GET(self, bool ind1, IM rr, IM cc)
    %
    %
      [varargout{1:nargout}] = casadiMEX(549, self, varargin{:});
    end
    function varargout = set(self,varargin)
    %SET 
    %
    %  SET(self, SX m, bool ind1, Sparsity sp)
    %  SET(self, SX m, bool ind1, Slice rr)
    %  SET(self, SX m, bool ind1, IM rr)
    %  SET(self, SX m, bool ind1, Slice rr, Slice cc)
    %  SET(self, SX m, bool ind1, Slice rr, IM cc)
    %  SET(self, SX m, bool ind1, IM rr, Slice cc)
    %  SET(self, SX m, bool ind1, IM rr, IM cc)
    %
    %
      [varargout{1:nargout}] = casadiMEX(550, self, varargin{:});
    end
    function varargout = get_nz(self,varargin)
    %GET_NZ 
    %
    %  SX = GET_NZ(self, bool ind1, Slice k)
    %  SX = GET_NZ(self, bool ind1, IM k)
    %
    %
      [varargout{1:nargout}] = casadiMEX(551, self, varargin{:});
    end
    function varargout = set_nz(self,varargin)
    %SET_NZ 
    %
    %  SET_NZ(self, SX m, bool ind1, Slice k)
    %  SET_NZ(self, SX m, bool ind1, IM k)
    %
    %
      [varargout{1:nargout}] = casadiMEX(552, self, varargin{:});
    end
    function varargout = uplus(self,varargin)
    %UPLUS 
    %
    %  SX = UPLUS(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(553, self, varargin{:});
    end
    function varargout = uminus(self,varargin)
    %UMINUS 
    %
    %  SX = UMINUS(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(554, self, varargin{:});
    end
    function varargout = printme(self,varargin)
    %PRINTME 
    %
    %  SX = PRINTME(self, SX y)
    %
    %
      [varargout{1:nargout}] = casadiMEX(561, self, varargin{:});
    end
    function varargout = T(self,varargin)
    %T [INTERNAL] 
    %
    %  SX = T(self)
    %
    %Transpose the matrix.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L519
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(562, self, varargin{:});
    end
    function varargout = print_split(self,varargin)
    %PRINT_SPLIT [INTERNAL] 
    %
    %  [{char} OUTPUT, {char} OUTPUT] = PRINT_SPLIT(self)
    %
    %Get strings corresponding to the nonzeros and the 
    %interdependencies.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L906
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_impl.hpp#L696-L700
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(568, self, varargin{:});
    end
    function varargout = disp(self,varargin)
    %DISP [INTERNAL] 
    %
    %  std::ostream & = DISP(self, bool more)
    %
    %Print a representation of the object.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L910
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(569, self, varargin{:});
    end
    function varargout = str(self,varargin)
    %STR [INTERNAL] 
    %
    %  char = STR(self, bool more)
    %
    %Get string representation.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L913
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(570, self, varargin{:});
    end
    function varargout = print_scalar(self,varargin)
    %PRINT_SCALAR [INTERNAL] 
    %
    %  std::ostream & = PRINT_SCALAR(self)
    %
    %Print scalar.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L916
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_impl.hpp#L609-L633
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(571, self, varargin{:});
    end
    function varargout = print_vector(self,varargin)
    %PRINT_VECTOR [INTERNAL] 
    %
    %  std::ostream & = PRINT_VECTOR(self, bool truncate)
    %
    %Print vector-style.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L919
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_impl.hpp#L636-L638
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(572, self, varargin{:});
    end
    function varargout = print_dense(self,varargin)
    %PRINT_DENSE [INTERNAL] 
    %
    %  std::ostream & = PRINT_DENSE(self, bool truncate)
    %
    %Print dense matrix-stype.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L922
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_impl.hpp#L686-L688
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(573, self, varargin{:});
    end
    function varargout = print_sparse(self,varargin)
    %PRINT_SPARSE [INTERNAL] 
    %
    %  std::ostream & = PRINT_SPARSE(self, bool truncate)
    %
    %Print sparse matrix style.
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L925
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_impl.hpp#L691-L693
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(574, self, varargin{:});
    end
    function varargout = clear(self,varargin)
    %CLEAR [INTERNAL] 
    %
    %  CLEAR(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(575, self, varargin{:});
    end
    function varargout = resize(self,varargin)
    %RESIZE [INTERNAL] 
    %
    %  RESIZE(self, int nrow, int ncol)
    %
    %
      [varargout{1:nargout}] = casadiMEX(576, self, varargin{:});
    end
    function varargout = reserve(self,varargin)
    %RESERVE [INTERNAL] 
    %
    %  RESERVE(self, int nnz)
    %  RESERVE(self, int nnz, int ncol)
    %
    %
      [varargout{1:nargout}] = casadiMEX(577, self, varargin{:});
    end
    function varargout = erase(self,varargin)
    %ERASE [INTERNAL] 
    %
    %  ERASE(self, [int] rr, bool ind1)
    %  ERASE(self, [int] rr, [int] cc, bool ind1)
    %
    %Erase a submatrix (leaving structural zeros in its place)
    %
    %Erase elements of a matrix
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19g
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L970
    %
    %
    %
    %.......
    %
    %::
    %
    %  ERASE(self, [int] rr, bool ind1)
    %
    %
    %
    %[INTERNAL] 
    %Erase a submatrix (leaving structural zeros in its place)
    %
    %Erase elements of a matrix
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19g
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L970
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
    %  ERASE(self, [int] rr, [int] cc, bool ind1)
    %
    %
    %
    %[INTERNAL] 
    %Erase a submatrix (leaving structural zeros in its place)
    %
    %Erase rows and/or columns of a matrix
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19f
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L962
    %
    %
    %
    %.............
    %
    %
      [varargout{1:nargout}] = casadiMEX(578, self, varargin{:});
    end
    function varargout = remove(self,varargin)
    %REMOVE [INTERNAL] 
    %
    %  REMOVE(self, [int] rr, [int] cc)
    %
    %Remove columns and rows.
    %
    %Remove/delete rows and/or columns of a matrix
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19h
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L977
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(579, self, varargin{:});
    end
    function varargout = enlarge(self,varargin)
    %ENLARGE [INTERNAL] 
    %
    %  ENLARGE(self, int nrow, int ncol, [int] rr, [int] cc, bool ind1)
    %
    %Enlarge matrix.
    %
    %Make the matrix larger by inserting empty rows and columns, keeping 
    %the 
    %existing non-zeros
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19i
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L985
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(580, self, varargin{:});
    end
    function varargout = sparsity(self,varargin)
    %SPARSITY [INTERNAL] 
    %
    %  Sparsity = SPARSITY(self)
    %
    %Get an owning reference to the sparsity pattern.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19j
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1012
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(581, self, varargin{:});
    end
    function varargout = get_output(self,varargin)
    %GET_OUTPUT 
    %
    %  SX = GET_OUTPUT(self, int oind)
    %
    %
      [varargout{1:nargout}] = casadiMEX(582, self, varargin{:});
    end
    function varargout = element_hash(self,varargin)
    %ELEMENT_HASH 
    %
    %  int = ELEMENT_HASH(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(587, self, varargin{:});
    end
    function varargout = is_regular(self,varargin)
    %IS_REGULAR 
    %
    %  bool = IS_REGULAR(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(588, self, varargin{:});
    end
    function varargout = is_smooth(self,varargin)
    %IS_SMOOTH 
    %
    %  bool = IS_SMOOTH(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(589, self, varargin{:});
    end
    function varargout = is_leaf(self,varargin)
    %IS_LEAF 
    %
    %  bool = IS_LEAF(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(590, self, varargin{:});
    end
    function varargout = is_commutative(self,varargin)
    %IS_COMMUTATIVE 
    %
    %  bool = IS_COMMUTATIVE(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(591, self, varargin{:});
    end
    function varargout = is_symbolic(self,varargin)
    %IS_SYMBOLIC 
    %
    %  bool = IS_SYMBOLIC(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(592, self, varargin{:});
    end
    function varargout = is_valid_input(self,varargin)
    %IS_VALID_INPUT 
    %
    %  bool = IS_VALID_INPUT(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(593, self, varargin{:});
    end
    function varargout = has_duplicates(self,varargin)
    %HAS_DUPLICATES 
    %
    %  bool = HAS_DUPLICATES(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(594, self, varargin{:});
    end
    function varargout = reset_input(self,varargin)
    %RESET_INPUT 
    %
    %  RESET_INPUT(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(595, self, varargin{:});
    end
    function varargout = is_constant(self,varargin)
    %IS_CONSTANT [INTERNAL] 
    %
    %  bool = IS_CONSTANT(self)
    %
    %Check if the matrix is constant (note that false negative 
    %answers are 
    %possible)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19v
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1138
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(596, self, varargin{:});
    end
    function varargout = is_call(self,varargin)
    %IS_CALL 
    %
    %  bool = IS_CALL(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(597, self, varargin{:});
    end
    function varargout = is_output(self,varargin)
    %IS_OUTPUT 
    %
    %  bool = IS_OUTPUT(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(598, self, varargin{:});
    end
    function varargout = has_output(self,varargin)
    %HAS_OUTPUT 
    %
    %  bool = HAS_OUTPUT(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(599, self, varargin{:});
    end
    function varargout = which_output(self,varargin)
    %WHICH_OUTPUT 
    %
    %  int = WHICH_OUTPUT(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(600, self, varargin{:});
    end
    function varargout = which_function(self,varargin)
    %WHICH_FUNCTION 
    %
    %  Function = WHICH_FUNCTION(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(601, self, varargin{:});
    end
    function varargout = is_integer(self,varargin)
    %IS_INTEGER [INTERNAL] 
    %
    %  bool = IS_INTEGER(self)
    %
    %Check if the matrix is integer-valued.
    %
    %(note that false negative answers are possible)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19w
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1170
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(602, self, varargin{:});
    end
    function varargout = is_zero(self,varargin)
    %IS_ZERO [INTERNAL] 
    %
    %  bool = IS_ZERO(self)
    %
    %check if the matrix is 0 (note that false negative answers are 
    %
    %possible)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19x
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1175
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(603, self, varargin{:});
    end
    function varargout = is_one(self,varargin)
    %IS_ONE [INTERNAL] 
    %
    %  bool = IS_ONE(self)
    %
    %check if the matrix is 1 (note that false negative answers are 
    %
    %possible)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19y
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1180
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(604, self, varargin{:});
    end
    function varargout = is_minus_one(self,varargin)
    %IS_MINUS_ONE [INTERNAL] 
    %
    %  bool = IS_MINUS_ONE(self)
    %
    %check if the matrix is -1 (note that false negative answers are
    % 
    %possible)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19z
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1185
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(605, self, varargin{:});
    end
    function varargout = is_eye(self,varargin)
    %IS_EYE [INTERNAL] 
    %
    %  bool = IS_EYE(self)
    %
    %check if the matrix is an identity matrix (note that false 
    %negative 
    %answers
    %
    %are possible)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1a0
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1192
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(606, self, varargin{:});
    end
    function varargout = op(self,varargin)
    %OP 
    %
    %  int = OP(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(607, self, varargin{:});
    end
    function varargout = is_op(self,varargin)
    %IS_OP 
    %
    %  bool = IS_OP(self, int op)
    %
    %
      [varargout{1:nargout}] = casadiMEX(608, self, varargin{:});
    end
    function varargout = has_zeros(self,varargin)
    %HAS_ZEROS [INTERNAL] 
    %
    %  bool = HAS_ZEROS(self)
    %
    %Check if the matrix has any zero entries which are not 
    %structural 
    %zeros.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1a1
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1203
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(609, self, varargin{:});
    end
    function varargout = nonzeros(self,varargin)
    %NONZEROS [INTERNAL] 
    %
    %  {SXElem} = NONZEROS(self)
    %
    %Get all nonzeros.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1a4
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1220
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(610, self, varargin{:});
    end
    function varargout = elements(self,varargin)
    %ELEMENTS [INTERNAL] 
    %
    %  {SXElem} = ELEMENTS(self)
    %
    %Get all elements.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1a3
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1213
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(611, self, varargin{:});
    end
    function varargout = to_double(self,varargin)
    %TO_DOUBLE 
    %
    %  double = TO_DOUBLE(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(612, self, varargin{:});
    end
    function varargout = to_int(self,varargin)
    %TO_INT 
    %
    %  int = TO_INT(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(613, self, varargin{:});
    end
    function varargout = name(self,varargin)
    %NAME 
    %
    %  char = NAME(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(614, self, varargin{:});
    end
    function varargout = dep(self,varargin)
    %DEP 
    %
    %  SX = DEP(self, int ch)
    %
    %
      [varargout{1:nargout}] = casadiMEX(615, self, varargin{:});
    end
    function varargout = n_dep(self,varargin)
    %N_DEP 
    %
    %  int = N_DEP(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(616, self, varargin{:});
    end
    function varargout = export_code(self,varargin)
    %EXPORT_CODE [INTERNAL] 
    %
    %  std::ostream & = EXPORT_CODE(self, char lang, struct options)
    %
    %Export matrix in specific language.
    %
    %lang: only 'matlab' supported for now
    %
    %::
    %
    %  * options:
    %  *   inline: Indicates if you want everything on a single line (default: False)
    %  *   name: Name of exported variable (default: 'm')
    %  *   indent_level: Level of indentation (default: 0)
    %  *   spoof_zero: Replace numerical zero by a 1e-200 (default: false)
    %  *               might be needed for matlab sparse construct,
    %  *               which doesn't allow numerical zero
    %  * 
    %
    %
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ac
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1297
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(622, self, varargin{:});
    end
    function varargout = info(self,varargin)
    %INFO 
    %
    %  struct = INFO(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(623, self, varargin{:});
    end
    function varargout = serialize(self,varargin)
    %SERIALIZE [INTERNAL] 
    %
    %  char = SERIALIZE(self)
    %  SERIALIZE(self, casadi::SerializingStream & s)
    %
    %Serialize an object.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ah
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1327
    %
    %
    %
    %.......
    %
    %::
    %
    %  SERIALIZE(self, casadi::SerializingStream & s)
    %
    %
    %
    %[INTERNAL] 
    %Serialize an object.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ah
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1327
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
    %  SERIALIZE(self)
    %
    %
    %
    %[INTERNAL] 
    %Serialize.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ae
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1312
    %
    %
    %
    %.............
    %
    %
      [varargout{1:nargout}] = casadiMEX(624, self, varargin{:});
    end
    function varargout = to_file(self,varargin)
    %TO_FILE [INTERNAL] 
    %
    %  TO_FILE(self, char filename, char format)
    %
    %Export numerical matrix to file
    %
    %Supported formats:
    %
    %
    %
    %::
    %
    %  *   - .mtx   Matrix Market (sparse)
    %  *   - .txt   Ascii full precision representation (sparse)
    %  *            Whitespace separated, aligned.
    %  *            Comments with # % or /
    %  *            Uses C locale
    %  *            Structural zeros represented by 00
    %  *            Does not scale well for large sparse matrices
    %  * 
    %
    %
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1347
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(626, self, varargin{:});
    end
    function varargout = paren(self,varargin)
    %PAREN 
    %
    %  SX = PAREN(self, Sparsity sp)
    %  SX = PAREN(self, IM rr)
    %  SX = PAREN(self, char rr)
    %  SX = PAREN(self, IM rr, IM cc)
    %  SX = PAREN(self, IM rr, char cc)
    %  SX = PAREN(self, char rr, IM cc)
    %  SX = PAREN(self, char rr, char cc)
    %
    %
      [varargout{1:nargout}] = casadiMEX(628, self, varargin{:});
    end
    function varargout = paren_asgn(self,varargin)
    %PAREN_ASGN 
    %
    %  PAREN_ASGN(self, SX m, Sparsity sp)
    %  PAREN_ASGN(self, SX m, IM rr)
    %  PAREN_ASGN(self, SX m, char rr)
    %  PAREN_ASGN(self, SX m, IM rr, IM cc)
    %  PAREN_ASGN(self, SX m, IM rr, char cc)
    %  PAREN_ASGN(self, SX m, char rr, IM cc)
    %  PAREN_ASGN(self, SX m, char rr, char cc)
    %
    %
      [varargout{1:nargout}] = casadiMEX(629, self, varargin{:});
    end
    function varargout = brace(self,varargin)
    %BRACE 
    %
    %  SX = BRACE(self, IM rr)
    %  SX = BRACE(self, char rr)
    %
    %
      [varargout{1:nargout}] = casadiMEX(630, self, varargin{:});
    end
    function varargout = setbrace(self,varargin)
    %SETBRACE 
    %
    %  SETBRACE(self, SX m, IM rr)
    %  SETBRACE(self, SX m, char rr)
    %
    %
      [varargout{1:nargout}] = casadiMEX(631, self, varargin{:});
    end
    function varargout = end(self,varargin)
    %END 
    %
    %  int = END(self, int i, int n)
    %
    %
      [varargout{1:nargout}] = casadiMEX(632, self, varargin{:});
    end
    function varargout = numel(self,varargin)
    %NUMEL 
    %
    %  int = NUMEL(self)
    %  int = NUMEL(self, int k)
    %  int = NUMEL(self, [int] k)
    %  int = NUMEL(self, char rr)
    %
    %
      [varargout{1:nargout}] = casadiMEX(633, self, varargin{:});
    end
    function varargout = ctranspose(self,varargin)
    %CTRANSPOSE 
    %
    %  SX = CTRANSPOSE(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(634, self, varargin{:});
    end
    function self = SX(varargin)
    %SX 
    %
    %  new_obj = SX()
    %  new_obj = SX(Sparsity sp)
    %  new_obj = SX(double val)
    %  new_obj = SX(SX m)
    %  new_obj = SX(int nrow, int ncol)
    %  new_obj = SX(Sparsity sp, SX d)
    %
    %
      self@casadi.MatrixCommon(SwigRef.Null);
      self@casadi.GenericExpressionCommon(SwigRef.Null);
      self@casadi.GenSX(SwigRef.Null);
      self@casadi.PrintableCommon(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = casadiMEX(635, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
        if self.swigPtr
          casadiMEX(636, self);
          self.SwigClear();
        end
    end
  end
  methods(Static)
    function varargout = binary(varargin)
    %BINARY 
    %
    %  SX = BINARY(int op, SX x, SX y)
    %
    %
     [varargout{1:nargout}] = casadiMEX(555, varargin{:});
    end
    function varargout = unary(varargin)
    %UNARY 
    %
    %  SX = UNARY(int op, SX x)
    %
    %
     [varargout{1:nargout}] = casadiMEX(556, varargin{:});
    end
    function varargout = scalar_matrix(varargin)
    %SCALAR_MATRIX 
    %
    %  SX = SCALAR_MATRIX(int op, SX x, SX y)
    %
    %
     [varargout{1:nargout}] = casadiMEX(557, varargin{:});
    end
    function varargout = matrix_scalar(varargin)
    %MATRIX_SCALAR 
    %
    %  SX = MATRIX_SCALAR(int op, SX x, SX y)
    %
    %
     [varargout{1:nargout}] = casadiMEX(558, varargin{:});
    end
    function varargout = matrix_matrix(varargin)
    %MATRIX_MATRIX 
    %
    %  SX = MATRIX_MATRIX(int op, SX x, SX y)
    %
    %
     [varargout{1:nargout}] = casadiMEX(559, varargin{:});
    end
    function varargout = call(varargin)
    %CALL 
    %
    %  {SX} = CALL(Function f, {SX} x)
    %  {SXElem} = CALL(Function f, {SXElem} x)
    %
    %
     [varargout{1:nargout}] = casadiMEX(560, varargin{:});
    end
    function varargout = set_max_depth(varargin)
    %SET_MAX_DEPTH 
    %
    %  SET_MAX_DEPTH(int eq_depth)
    %
    %
     [varargout{1:nargout}] = casadiMEX(563, varargin{:});
    end
    function varargout = get_max_depth(varargin)
    %GET_MAX_DEPTH 
    %
    %  int = GET_MAX_DEPTH()
    %
    %
     [varargout{1:nargout}] = casadiMEX(564, varargin{:});
    end
    function varargout = get_input(varargin)
    %GET_INPUT 
    %
    %  {SX} = GET_INPUT(Function f)
    %
    %
     [varargout{1:nargout}] = casadiMEX(565, varargin{:});
    end
    function varargout = get_free(varargin)
    %GET_FREE 
    %
    %  {SX} = GET_FREE(Function f)
    %
    %
     [varargout{1:nargout}] = casadiMEX(566, varargin{:});
    end
    function varargout = type_name(varargin)
    %TYPE_NAME 
    %
    %  char = TYPE_NAME()
    %
    %
     [varargout{1:nargout}] = casadiMEX(567, varargin{:});
    end
    function varargout = triplet(varargin)
    %TRIPLET 
    %
    %  SX = TRIPLET([int] row, [int] col, SX d)
    %  SX = TRIPLET([int] row, [int] col, SX d, [int,int] rc)
    %  SX = TRIPLET([int] row, [int] col, SX d, int nrow, int ncol)
    %
    %
     [varargout{1:nargout}] = casadiMEX(583, varargin{:});
    end
    function varargout = inf(varargin)
    %INF [INTERNAL] 
    %
    %  SX = INF(int nrow, int ncol)
    %  SX = INF([int,int] rc)
    %  SX = INF(Sparsity sp)
    %
    %create a matrix with all inf
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19k
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1043
    %
    %
    %
    %.......
    %
    %::
    %
    %  INF(int nrow, int ncol)
    %
    %
    %
    %[INTERNAL] 
    %create a matrix with all inf
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19k
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1042
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
    %  INF([int,int] rc)
    %
    %
    %
    %[INTERNAL] 
    %create a matrix with all inf
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19k
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1043
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
    %  INF(Sparsity sp)
    %
    %
    %
    %[INTERNAL] 
    %create a matrix with all inf
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19k
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1041
    %
    %
    %
    %.............
    %
    %
     [varargout{1:nargout}] = casadiMEX(584, varargin{:});
    end
    function varargout = nan(varargin)
    %NAN [INTERNAL] 
    %
    %  SX = NAN(int nrow, int ncol)
    %  SX = NAN([int,int] rc)
    %  SX = NAN(Sparsity sp)
    %
    %create a matrix with all nan
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19l
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1052
    %
    %
    %
    %.......
    %
    %::
    %
    %  NAN(int nrow, int ncol)
    %
    %
    %
    %[INTERNAL] 
    %create a matrix with all nan
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19l
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1051
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
    %  NAN([int,int] rc)
    %
    %
    %
    %[INTERNAL] 
    %create a matrix with all nan
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19l
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1052
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
    %  NAN(Sparsity sp)
    %
    %
    %
    %[INTERNAL] 
    %create a matrix with all nan
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_19l
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1050
    %
    %
    %
    %.............
    %
    %
     [varargout{1:nargout}] = casadiMEX(585, varargin{:});
    end
    function varargout = eye(varargin)
    %EYE 
    %
    %  SX = EYE(int n)
    %
    %
     [varargout{1:nargout}] = casadiMEX(586, varargin{:});
    end
    function varargout = set_precision(varargin)
    %SET_PRECISION 
    %
    %  SET_PRECISION(int precision)
    %
    %
     [varargout{1:nargout}] = casadiMEX(617, varargin{:});
    end
    function varargout = set_width(varargin)
    %SET_WIDTH 
    %
    %  SET_WIDTH(int width)
    %
    %
     [varargout{1:nargout}] = casadiMEX(618, varargin{:});
    end
    function varargout = set_scientific(varargin)
    %SET_SCIENTIFIC 
    %
    %  SET_SCIENTIFIC(bool scientific)
    %
    %
     [varargout{1:nargout}] = casadiMEX(619, varargin{:});
    end
    function varargout = rng(varargin)
    %RNG 
    %
    %  RNG(int seed)
    %
    %
     [varargout{1:nargout}] = casadiMEX(620, varargin{:});
    end
    function varargout = rand(varargin)
    %RAND [INTERNAL] 
    %
    %  SX = RAND(int nrow, int ncol)
    %  SX = RAND([int,int] rc)
    %  SX = RAND(Sparsity sp)
    %
    %Create a matrix with uniformly distributed random numbers.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ab
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1279
    %
    %
    %
    %.......
    %
    %::
    %
    %  RAND(int nrow, int ncol)
    %
    %
    %
    %[INTERNAL] 
    %Create a matrix with uniformly distributed random numbers.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ab
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1275
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
    %  RAND([int,int] rc)
    %
    %
    %
    %[INTERNAL] 
    %Create a matrix with uniformly distributed random numbers.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ab
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1279
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
    %  RAND(Sparsity sp)
    %
    %
    %
    %[INTERNAL] 
    %Create a matrix with uniformly distributed random numbers.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ab
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/matrix_decl.hpp#L1278
    %
    %
    %
    %.............
    %
    %
     [varargout{1:nargout}] = casadiMEX(621, varargin{:});
    end
    function varargout = deserialize(varargin)
    %DESERIALIZE 
    %
    %  SX = DESERIALIZE(std::istream & stream)
    %  SX = DESERIALIZE(casadi::DeserializingStream & s)
    %  SX = DESERIALIZE(char s)
    %
    %
     [varargout{1:nargout}] = casadiMEX(625, varargin{:});
    end
    function varargout = from_file(varargin)
    %FROM_FILE 
    %
    %  DM = FROM_FILE(char filename, char format_hint)
    %
    %
     [varargout{1:nargout}] = casadiMEX(627, varargin{:});
    end
  end
end
