classdef  GenMX < casadi.GenericMatrixCommon & casadi.SparsityInterfaceCommon
    %GENMX 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = casadiMEX(3, self);
    end
    function varargout = nnz(self,varargin)
    %NNZ [INTERNAL] 
    %
    %  int = NNZ(self)
    %
    %Get the number of (structural) non-zero elements.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1an
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L84
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1298-L1300
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(344, self, varargin{:});
    end
    function varargout = nnz_lower(self,varargin)
    %NNZ_LOWER [INTERNAL] 
    %
    %  int = NNZ_LOWER(self)
    %
    %Get the number of non-zeros in the lower triangular half.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ao
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L89
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1303-L1305
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(345, self, varargin{:});
    end
    function varargout = nnz_upper(self,varargin)
    %NNZ_UPPER [INTERNAL] 
    %
    %  int = NNZ_UPPER(self)
    %
    %Get the number of non-zeros in the upper triangular half.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ap
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L94
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1308-L1310
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(346, self, varargin{:});
    end
    function varargout = nnz_diag(self,varargin)
    %NNZ_DIAG [INTERNAL] 
    %
    %  int = NNZ_DIAG(self)
    %
    %Get get the number of non-zeros on the diagonal.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1aq
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L99
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1313-L1315
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(347, self, varargin{:});
    end
    function varargout = numel(self,varargin)
    %NUMEL [INTERNAL] 
    %
    %  int = NUMEL(self)
    %
    %Get the number of elements.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ar
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L104
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1318-L1320
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(348, self, varargin{:});
    end
    function varargout = size1(self,varargin)
    %SIZE1 [INTERNAL] 
    %
    %  int = SIZE1(self)
    %
    %Get the first dimension (i.e. number of rows)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1as
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L109
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1323-L1325
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(349, self, varargin{:});
    end
    function varargout = rows(self,varargin)
    %ROWS [INTERNAL] 
    %
    %  int = ROWS(self)
    %
    %Get the number of rows, Octave-style syntax.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1at
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L114
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L114-L114
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(350, self, varargin{:});
    end
    function varargout = size2(self,varargin)
    %SIZE2 [INTERNAL] 
    %
    %  int = SIZE2(self)
    %
    %Get the second dimension (i.e. number of columns)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1au
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L119
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1328-L1330
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(351, self, varargin{:});
    end
    function varargout = columns(self,varargin)
    %COLUMNS [INTERNAL] 
    %
    %  int = COLUMNS(self)
    %
    %Get the number of columns, Octave-style syntax.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1av
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L124
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L124-L124
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(352, self, varargin{:});
    end
    function varargout = dim(self,varargin)
    %DIM [INTERNAL] 
    %
    %  char = DIM(self, bool with_nz)
    %
    %Get string representation of dimensions.
    %
    %The representation is e.g. "4x5" or "4x5,10nz"
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1aw
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L131
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1343-L1345
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(353, self, varargin{:});
    end
    function varargout = size(self,varargin)
    %SIZE [INTERNAL] 
    %
    %  [int,int] = SIZE(self)
    %  int = SIZE(self, int axis)
    %
    %Get the size along a particular dimensions.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ay
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L141
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1338-L1340
    %
    %
    %
    %.......
    %
    %::
    %
    %  SIZE(self, int axis)
    %
    %
    %
    %[INTERNAL] 
    %Get the size along a particular dimensions.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ay
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L141
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1338-L1340
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
    %  SIZE(self)
    %
    %
    %
    %[INTERNAL] 
    %Get the shape.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1ax
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L136
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1333-L1335
    %
    %
    %
    %.............
    %
    %
      out = casadiMEX(354, self, varargin{:});
      if nargout<=1
        varargout{1}=out;
      else
        nargoutchk(length(out),length(out))
        for i=1:nargout
          varargout{i} = out(i);
        end
      end
    end
    function varargout = is_empty(self,varargin)
    %IS_EMPTY [INTERNAL] 
    %
    %  bool = IS_EMPTY(self, bool both)
    %
    %Check if the sparsity is empty, i.e. if one of the dimensions is
    % zero.
    %
    %(or optionally both dimensions)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1az
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L148
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L148-L148
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(355, self, varargin{:});
    end
    function varargout = is_dense(self,varargin)
    %IS_DENSE [INTERNAL] 
    %
    %  bool = IS_DENSE(self)
    %
    %Check if the matrix expression is dense.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b0
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L153
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L153-L153
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(356, self, varargin{:});
    end
    function varargout = is_scalar(self,varargin)
    %IS_SCALAR [INTERNAL] 
    %
    %  bool = IS_SCALAR(self, bool scalar_and_dense)
    %
    %Check if the matrix expression is scalar.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b1
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L158
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1348-L1350
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(357, self, varargin{:});
    end
    function varargout = is_square(self,varargin)
    %IS_SQUARE [INTERNAL] 
    %
    %  bool = IS_SQUARE(self)
    %
    %Check if the matrix expression is square.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b2
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L163
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L163-L163
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(358, self, varargin{:});
    end
    function varargout = is_vector(self,varargin)
    %IS_VECTOR [INTERNAL] 
    %
    %  bool = IS_VECTOR(self)
    %
    %Check if the matrix is a row or column vector.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b3
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L168
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L168-L168
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(359, self, varargin{:});
    end
    function varargout = is_row(self,varargin)
    %IS_ROW [INTERNAL] 
    %
    %  bool = IS_ROW(self)
    %
    %Check if the matrix is a row vector (i.e.  size1()==1)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b4
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L173
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L173-L173
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(360, self, varargin{:});
    end
    function varargout = is_column(self,varargin)
    %IS_COLUMN [INTERNAL] 
    %
    %  bool = IS_COLUMN(self)
    %
    %Check if the matrix is a column vector (i.e.  size2()==1)
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b5
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L178
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L178-L178
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(361, self, varargin{:});
    end
    function varargout = is_triu(self,varargin)
    %IS_TRIU [INTERNAL] 
    %
    %  bool = IS_TRIU(self)
    %
    %Check if the matrix is upper triangular.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b6
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L183
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L183-L183
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(362, self, varargin{:});
    end
    function varargout = is_tril(self,varargin)
    %IS_TRIL [INTERNAL] 
    %
    %  bool = IS_TRIL(self)
    %
    %Check if the matrix is lower triangular.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b7
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L188
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L188-L188
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(363, self, varargin{:});
    end
    function varargout = row(self,varargin)
    %ROW [INTERNAL] 
    %
    %  [int] = ROW(self)
    %  int = ROW(self, int el)
    %
    %Get the sparsity pattern. See the Sparsity class for details.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b8
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L200
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L200-L200
    %
    %
    %
    %.......
    %
    %::
    %
    %  ROW(self)
    %
    %
    %
    %[INTERNAL] 
    %Get the sparsity pattern. See the Sparsity class for details.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b8
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L194
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L194-L194
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
    %  ROW(self, int el)
    %
    %
    %
    %[INTERNAL] 
    %Get the sparsity pattern. See the Sparsity class for details.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b8
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L200
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L200-L200
    %
    %
    %
    %.............
    %
    %
      [varargout{1:nargout}] = casadiMEX(364, self, varargin{:});
    end
    function varargout = colind(self,varargin)
    %COLIND [INTERNAL] 
    %
    %  [int] = COLIND(self)
    %  int = COLIND(self, int col)
    %
    %Get the sparsity pattern. See the Sparsity class for details.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b8
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L201
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L201-L201
    %
    %
    %
    %.......
    %
    %::
    %
    %  COLIND(self)
    %
    %
    %
    %[INTERNAL] 
    %Get the sparsity pattern. See the Sparsity class for details.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b8
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L195
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L195-L195
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
    %  COLIND(self, int col)
    %
    %
    %
    %[INTERNAL] 
    %Get the sparsity pattern. See the Sparsity class for details.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b8
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L201
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L201-L201
    %
    %
    %
    %.............
    %
    %
      [varargout{1:nargout}] = casadiMEX(365, self, varargin{:});
    end
    function varargout = sparsity(self,varargin)
    %SPARSITY [INTERNAL] 
    %
    %  Sparsity = SPARSITY(self)
    %
    %Get the sparsity pattern.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1b9
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L207
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1293-L1295
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(366, self, varargin{:});
    end
    function self = GenMX(varargin)
    %GENMX 
    %
    %  new_obj = GENMX()
    %
    %
      self@casadi.GenericMatrixCommon(SwigRef.Null);
      self@casadi.SparsityInterfaceCommon(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = casadiMEX(370, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
        if self.swigPtr
          casadiMEX(371, self);
          self.SwigClear();
        end
    end
  end
  methods(Static)
    function varargout = sym(varargin)
    %SYM [INTERNAL] 
    %
    %  MX = SYM(char name, int nrow, int ncol)
    %  MX = SYM(char name, [int,int] rc)
    %  MX = SYM(char name, Sparsity sp)
    %  {MX} = SYM(char name, Sparsity sp, int p)
    %  {MX} = SYM(char name, int nrow, int ncol, int p)
    %  {{MX}} = SYM(char name, Sparsity sp, int p, int r)
    %  {{MX}} = SYM(char name, int nrow, int ncol, int p, int r)
    %
    %Create a vector of length r of vectors of length p.
    %
    %with nrow-by-ncol symbolic primitives
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dg
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1253
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1253-L1255
    %
    %
    %
    %.......
    %
    %::
    %
    %  SYM(char name, [int,int] rc)
    %
    %
    %
    %[INTERNAL] 
    %Construct a symbolic primitive with given dimensions.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1db
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1213
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1213-L1215
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
    %  SYM(char name, int nrow, int ncol, int p)
    %
    %
    %
    %[INTERNAL] 
    %Create a vector of length p with nrow-by-ncol symbolic 
    %primitives.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1de
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1234
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1234-L1237
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
    %  SYM(char name, int nrow, int ncol)
    %
    %
    %
    %[INTERNAL] 
    %Create an nrow-by-ncol symbolic primitive.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1da
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1206
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1206-L1208
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
    %  SYM(char name, Sparsity sp, int p)
    %
    %
    %
    %[INTERNAL] 
    %Create a vector of length p with with matrices.
    %
    %with symbolic primitives of given sparsity
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dd
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1229
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1355-L1365
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
    %  SYM(char name, Sparsity sp)
    %
    %
    %
    %[INTERNAL] 
    %Create symbolic primitive with a given sparsity pattern.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dc
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1220
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1220-L1222
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
    %  SYM(char name, Sparsity sp, int p, int r)
    %
    %
    %
    %[INTERNAL] 
    %Create a vector of length r of vectors of length p with.
    %
    %symbolic primitives with given sparsity
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1df
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1245
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1368-L1378
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
    %  SYM(char name, int nrow, int ncol, int p, int r)
    %
    %
    %
    %[INTERNAL] 
    %Create a vector of length r of vectors of length p.
    %
    %with nrow-by-ncol symbolic primitives
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dg
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1253
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1253-L1255
    %
    %
    %
    %.............
    %
    %
     [varargout{1:nargout}] = casadiMEX(367, varargin{:});
    end
    function varargout = zeros(varargin)
    %ZEROS [INTERNAL] 
    %
    %  MX = ZEROS(int nrow, int ncol)
    %  MX = ZEROS([int,int] rc)
    %  MX = ZEROS(Sparsity sp)
    %
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries zero.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dh
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1266
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1266-L1268
    %
    %
    %
    %.......
    %
    %::
    %
    %  ZEROS(int nrow, int ncol)
    %
    %
    %
    %[INTERNAL] 
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries zero.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dh
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1262
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1262-L1264
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
    %  ZEROS([int,int] rc)
    %
    %
    %
    %[INTERNAL] 
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries zero.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dh
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1266
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1266-L1268
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
    %  ZEROS(Sparsity sp)
    %
    %
    %
    %[INTERNAL] 
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries zero.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1dh
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1265
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1265-L1265
    %
    %
    %
    %.............
    %
    %
     [varargout{1:nargout}] = casadiMEX(368, varargin{:});
    end
    function varargout = ones(varargin)
    %ONES [INTERNAL] 
    %
    %  MX = ONES(int nrow, int ncol)
    %  MX = ONES([int,int] rc)
    %  MX = ONES(Sparsity sp)
    %
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries one.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1di
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1279
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1279-L1281
    %
    %
    %
    %.......
    %
    %::
    %
    %  ONES(int nrow, int ncol)
    %
    %
    %
    %[INTERNAL] 
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries one.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1di
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1275
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1275-L1277
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
    %  ONES([int,int] rc)
    %
    %
    %
    %[INTERNAL] 
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries one.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1di
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1279
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1279-L1281
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
    %  ONES(Sparsity sp)
    %
    %
    %
    %[INTERNAL] 
    %Create a dense matrix or a matrix with specified sparsity with 
    %all 
    %entries one.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1di
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1278
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/main/casadi/core/generic_matrix.hpp#L1278-L1278
    %
    %
    %
    %.............
    %
    %
     [varargout{1:nargout}] = casadiMEX(369, varargin{:});
    end
  end
end
