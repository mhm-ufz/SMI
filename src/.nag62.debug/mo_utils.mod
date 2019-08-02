VL 0 0 4 0 0 0
MODULE MO_UTILS,2 0 0
FILE 0,/Users/thober/lib/smi/src/mo_utils.f90
USE MO_KIND 2,ONLY:DPC=>DPC,SPC=>SPC,I8=>I8,I4=>I4,DP=>DP,SP=>SP
GMODPROC ARANGE: ARANGE_SP
GMODPROC ARANGE: ARANGE_DP
GMODPROC ARANGE: ARANGE_I8
GMODPROC ARANGE: ARANGE_I4
GMODPROC CUMSUM: CUMSUM_SPC
GMODPROC CUMSUM: CUMSUM_DPC
GMODPROC CUMSUM: CUMSUM_SP
GMODPROC CUMSUM: CUMSUM_DP
GMODPROC CUMSUM: CUMSUM_I8
GMODPROC CUMSUM: CUMSUM_I4
GMODPROC EQ: EQUAL_DP
GMODPROC EQ: EQUAL_SP
GMODPROC EQUAL: EQUAL_DP
GMODPROC EQUAL: EQUAL_SP
GMODPROC GE: GREATEREQUAL_DP
GMODPROC GE: GREATEREQUAL_SP
GMODPROC GREATEREQUAL: GREATEREQUAL_DP
GMODPROC GREATEREQUAL: GREATEREQUAL_SP
GMODPROC IMAXLOC: IMAXLOC_DP
GMODPROC IMAXLOC: IMAXLOC_SP
GMODPROC IMAXLOC: IMAXLOC_I8
GMODPROC IMAXLOC: IMAXLOC_I4
GMODPROC IMINLOC: IMINLOC_DP
GMODPROC IMINLOC: IMINLOC_SP
GMODPROC IMINLOC: IMINLOC_I8
GMODPROC IMINLOC: IMINLOC_I4
GMODPROC IS_FINITE: IS_FINITE_DP
GMODPROC IS_FINITE: IS_FINITE_SP
GMODPROC IS_NAN: IS_NAN_DP
GMODPROC IS_NAN: IS_NAN_SP
GMODPROC IS_NORMAL: IS_NORMAL_DP
GMODPROC IS_NORMAL: IS_NORMAL_SP
GMODPROC LE: LESSEREQUAL_DP
GMODPROC LE: LESSEREQUAL_SP
GMODPROC LESSEREQUAL: LESSEREQUAL_DP
GMODPROC LESSEREQUAL: LESSEREQUAL_SP
GMODPROC LINSPACE: LINSPACE_SP
GMODPROC LINSPACE: LINSPACE_DP
GMODPROC LINSPACE: LINSPACE_I8
GMODPROC LINSPACE: LINSPACE_I4
GMODPROC LOCATE: LOCATE_1D_SP
GMODPROC LOCATE: LOCATE_1D_DP
GMODPROC LOCATE: LOCATE_0D_SP
GMODPROC LOCATE: LOCATE_0D_DP
GMODPROC NE: NOTEQUAL_DP
GMODPROC NE: NOTEQUAL_SP
GMODPROC NOTEQUAL: NOTEQUAL_DP
GMODPROC NOTEQUAL: NOTEQUAL_SP
GMODPROC SPECIAL_VALUE: SPECIAL_VALUE_DP
GMODPROC SPECIAL_VALUE: SPECIAL_VALUE_SP
GMODPROC SWAP: SWAP_VEC_I4
GMODPROC SWAP: SWAP_VEC_SP
GMODPROC SWAP: SWAP_VEC_DP
GMODPROC SWAP: SWAP_XY_I4
GMODPROC SWAP: SWAP_XY_SP
GMODPROC SWAP: SWAP_XY_DP
PROC NOTEQUAL_DP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,8,5,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,8,5,0,0,0,103,0,1000,0,1
ENDPROC
PROC LINSPACE_DP,3,8,0,17,0: 2,8,5,0,1,0,40281,2 (1,2,0: 1,2,1,2),0,40000,0
VAR LOWER,3,,: 2,8,5,0,0,0,103,0,0,0,1
VAR UPPER,3,,: 2,8,5,0,0,0,103,0,0,0,1
VAR NSTEP,3,,: 1,4,3,0,0,0,103,0,0,0,1
ENDPROC
PROC ITOUPPER,1,8,0,17,0: 3,1,a,-4,0,0,40281,2,800,0,0
VAR LOWER,3,,: 3,1,a,-1,0,0,103,0,0,0,1
ENDPROC
PROC GREATEREQUAL_SP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,4,4,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,4,4,0,0,0,103,0,1000,0,1
ENDPROC
PROC LINSPACE_SP,3,8,0,17,0: 2,4,4,0,1,0,40281,2 (1,2,0: 1,2,1,2),0,40000,0
VAR LOWER,3,,: 2,4,4,0,0,0,103,0,0,0,1
VAR UPPER,3,,: 2,4,4,0,0,0,103,0,0,0,1
VAR NSTEP,3,,: 1,4,3,0,0,0,103,0,0,0,1
ENDPROC
PROC LOCATE_0D_DP,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,8,5,0,0,0,103,0,0,1000,1
ENDPROC
PROC GREATEREQUAL_DP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,8,5,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,8,5,0,0,0,103,0,1000,0,1
ENDPROC
PROC LOCATE_0D_SP,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,4,4,0,0,0,103,0,0,1000,1
ENDPROC
PROC LESSEREQUAL_SP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,4,4,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,4,4,0,0,0,103,0,1000,0,1
ENDPROC
PROC SWAP_XY_DP,2,8,0,17,0: 8,0,0,0,0,0,40200,2,1800,40000,0
VAR X,3,,: 2,8,5,0,0,0,183,0,1000,0,3
VAR Y,3,,: 2,8,5,0,0,0,183,0,1000,0,3
ENDPROC
PROC LOCATE_1D_DP,2,8,0,17,0: 1,4,3,0,1,0,40681,2 (1,8,0: 5,3,1,0),1000000,40001,0
VAR X,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC LOCATE_1D_SP,2,8,0,17,0: 1,4,3,0,1,0,40681,2 (1,8,0: 5,3,1,0),1000000,40001,0
VAR X,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC SWAP_XY_SP,2,8,0,17,0: 8,0,0,0,0,0,40200,2,1800,40000,0
VAR X,3,,: 2,4,4,0,0,0,183,0,1000,0,3
VAR Y,3,,: 2,4,4,0,0,0,183,0,1000,0,3
ENDPROC
PROC IMAXLOC_I4,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 1,4,3,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC ARANGE_I4,2,8,0,17,0: 1,4,3,0,1,0,40681,2 (1,8,0: 5,3,1,0),1000000,40001,0
VAR LOWER,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR UPPER,3,,: 1,4,3,0,0,0,113,0,0,0,1
ENDPROC
PROC ARANGE_I8,2,8,0,17,0: 1,8,e,0,1,0,40681,2 (1,8,0: 5,3,1,0),1000000,40001,0
VAR LOWER,3,,: 1,8,e,0,0,0,103,0,0,0,1
VAR UPPER,3,,: 1,8,e,0,0,0,113,0,0,0,1
ENDPROC
PROC SWAP_XY_I4,2,8,0,17,0: 8,0,0,0,0,0,40200,2,1800,40000,0
VAR X,3,,: 1,4,3,0,0,0,183,0,1000,0,3
VAR Y,3,,: 1,4,3,0,0,0,183,0,1000,0,3
ENDPROC
PROC IMAXLOC_SP,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 2,4,4,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC ARANGE_DP,2,8,0,17,0: 2,8,5,0,1,0,40681,2 (1,8,0: 5,3,1,0),1000000,40001,0
VAR LOWER,3,,: 2,8,5,0,0,0,3,0,0,0,1
VAR UPPER,3,,: 2,8,5,0,0,0,13,0,0,0,1
ENDPROC
PROC IMAXLOC_DP,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 2,8,5,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC ARANGE_SP,2,8,0,17,0: 2,4,4,0,1,0,40681,2 (1,8,0: 5,3,1,0),1000000,40001,0
VAR LOWER,3,,: 2,4,4,0,0,0,3,0,0,0,1
VAR UPPER,3,,: 2,4,4,0,0,0,13,0,0,0,1
ENDPROC
PROC SWAP_VEC_DP,3,8,0,17,0: 8,0,0,0,0,0,40200,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,183,0 (1,5,0: 5,3,1,0),0,0,3
VAR I1,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR I2,3,,: 1,4,3,0,0,0,103,0,0,0,1
ENDPROC
PROC IMINLOC_I4,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 1,4,3,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC IMINLOC_I8,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 1,8,e,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC IMINLOC_SP,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 2,4,4,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC SWAP_VEC_I4,3,8,0,17,0: 8,0,0,0,0,0,40200,2,0,40000,0
VAR X,3,,: 1,4,3,0,1,0,183,0 (1,5,0: 5,3,1,0),0,0,3
VAR I1,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR I2,3,,: 1,4,3,0,0,0,103,0,0,0,1
ENDPROC
PROC IMINLOC_DP,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 2,8,5,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC IS_FINITE_DP,1,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:IEEE_IS_FINITE=>IEEE_IS_FINITE
VAR A,3,,: 2,8,5,0,0,0,3,0,1000,1000,1
ENDPROC
PROC IS_FINITE_SP,1,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:IEEE_IS_FINITE=>IEEE_IS_FINITE
VAR A,3,,: 2,4,4,0,0,0,3,0,1000,1000,1
ENDPROC
PROC SWAP_VEC_SP,3,8,0,17,0: 8,0,0,0,0,0,40200,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,183,0 (1,5,0: 5,3,1,0),0,0,3
VAR I1,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR I2,3,,: 1,4,3,0,0,0,103,0,0,0,1
ENDPROC
PROC CUMSUM_I4,1,8,0,17,0: 1,4,3,0,1,0,40381,2 (1,2,0: 1,2,1,2),0,40000,0
VAR ARR,3,,: 1,4,3,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CUMSUM_I8,1,8,0,17,0: 1,8,e,0,1,0,40381,2 (1,2,0: 1,2,1,2),0,40000,0
VAR ARR,3,,: 1,8,e,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CUMSUM_DP,1,8,0,17,0: 2,8,5,0,1,0,40381,2 (1,2,0: 1,2,1,2),0,40000,0
VAR ARR,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC SPECIAL_VALUE_SP,2,8,0,17,0: 2,4,4,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:IEEE_POSITIVE_ZERO=>IEEE_POSITIVE_ZERO,IEEE_NEGATIVE_ZERO=>IEEE_NEGATIVE_ZERO,IEEE_POSITIVE_NORMAL=>IEEE_POSITIVE_NORMAL,IEEE_NEGATIVE_NORMAL=>IEEE_NEGATIVE_NORMAL,IEEE_POSITIVE_DENORMAL=>IEEE_POSITIVE_DENORMAL,IEEE_NEGATIVE_DENORMAL=>IEEE_NEGATIVE_DENORMAL,IEEE_POSITIVE_INF=>IEEE_POSITIVE_INF,IEEE_NEGATIVE_INF=>IEEE_NEGATIVE_INF,IEEE_QUIET_NAN=>IEEE_QUIET_NAN,IEEE_SIGNALING_NAN=>IEEE_SIGNALING_NAN,IEEE_VALUE=>IEEE_VALUE
VAR X,3,,: 2,4,4,0,0,0,3,0,1000,1000,1
VAR IEEE,3,,: 3,1,a,-1,0,0,3,0,1000,1000,1
ENDPROC
PROC LESSEREQUAL_DP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,8,5,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,8,5,0,0,0,103,0,1000,0,1
ENDPROC
PROC IS_NAN_SP,1,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:ISNAN=>IEEE_IS_NAN
VAR A,3,,: 2,4,4,0,0,0,3,0,1000,1000,1
ENDPROC
PROC CUMSUM_SP,1,8,0,17,0: 2,4,4,0,1,0,40381,2 (1,2,0: 1,2,1,2),0,40000,0
VAR ARR,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC IS_NAN_DP,1,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:ISNAN=>IEEE_IS_NAN
VAR A,3,,: 2,8,5,0,0,0,3,0,1000,1000,1
ENDPROC
PROC CUMSUM_DPC,1,8,0,17,0: 5,8,9,0,1,0,40381,2 (1,2,0: 1,2,1,2),0,40000,0
VAR ARR,3,,: 5,8,9,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC IS_NORMAL_SP,1,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:IEEE_IS_NORMAL=>IEEE_IS_NORMAL
VAR A,3,,: 2,4,4,0,0,0,3,0,1000,1000,1
ENDPROC
PROC CUMSUM_SPC,1,8,0,17,0: 5,4,8,0,1,0,40381,2 (1,2,0: 1,2,1,2),0,40000,0
VAR ARR,3,,: 5,4,8,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC IS_NORMAL_DP,1,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:IEEE_IS_NORMAL=>IEEE_IS_NORMAL
VAR A,3,,: 2,8,5,0,0,0,3,0,1000,1000,1
ENDPROC
PROC SPECIAL_VALUE_DP,2,8,0,17,0: 2,8,5,0,0,0,40281,2,1800,40000,0
USE IEEE_ARITHMETIC 0,ONLY:IEEE_POSITIVE_ZERO=>IEEE_POSITIVE_ZERO,IEEE_NEGATIVE_ZERO=>IEEE_NEGATIVE_ZERO,IEEE_POSITIVE_NORMAL=>IEEE_POSITIVE_NORMAL,IEEE_NEGATIVE_NORMAL=>IEEE_NEGATIVE_NORMAL,IEEE_POSITIVE_DENORMAL=>IEEE_POSITIVE_DENORMAL,IEEE_NEGATIVE_DENORMAL=>IEEE_NEGATIVE_DENORMAL,IEEE_POSITIVE_INF=>IEEE_POSITIVE_INF,IEEE_NEGATIVE_INF=>IEEE_NEGATIVE_INF,IEEE_QUIET_NAN=>IEEE_QUIET_NAN,IEEE_SIGNALING_NAN=>IEEE_SIGNALING_NAN,IEEE_VALUE=>IEEE_VALUE
VAR X,3,,: 2,8,5,0,0,0,3,0,1000,1000,1
VAR IEEE,3,,: 3,1,a,-1,0,0,3,0,1000,1000,1
ENDPROC
PROC EQUAL_SP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,4,4,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,4,4,0,0,0,103,0,1000,0,1
ENDPROC
PROC LINSPACE_I4,3,8,0,17,0: 1,4,3,0,1,0,40281,2 (1,2,0: 1,2,1,2),0,40000,0
VAR LOWER,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR UPPER,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR NSTEP,3,,: 1,4,3,0,0,0,103,0,0,0,1
ENDPROC
PROC EQUAL_DP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,8,5,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,8,5,0,0,0,103,0,1000,0,1
ENDPROC
PROC IMAXLOC_I8,2,8,0,17,0: 1,4,3,0,0,0,40281,2,0,40000,0
VAR ARR,3,,: 1,8,e,0,1,0,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,13,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC NOTEQUAL_SP,2,8,0,17,0: 4,4,7,0,0,0,40281,2,1800,40000,0
VAR A,3,,: 2,4,4,0,0,0,103,0,1000,0,1
VAR B,3,,: 2,4,4,0,0,0,103,0,1000,0,1
ENDPROC
PROC LINSPACE_I8,3,8,0,17,0: 1,8,e,0,1,0,40281,2 (1,2,0: 1,2,1,2),0,40000,0
VAR LOWER,3,,: 1,8,e,0,0,0,103,0,0,0,1
VAR UPPER,3,,: 1,8,e,0,0,0,103,0,0,0,1
VAR NSTEP,3,,: 1,4,3,0,0,0,103,0,0,0,1
ENDPROC
END
