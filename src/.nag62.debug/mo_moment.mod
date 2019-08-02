VL 0 0 4 0 0 0
MODULE MO_MOMENT,2 0 0
FILE 0,/Users/thober/lib/smi/src/mo_moment.f90
USE MO_KIND 2,ONLY:DP=>DP,SP=>SP,I4=>I4
GMODPROC ABSDEV: ABSDEV_DP
GMODPROC ABSDEV: ABSDEV_SP
GMODPROC AVERAGE: AVERAGE_DP
GMODPROC AVERAGE: AVERAGE_SP
GMODPROC CENTRAL_MOMENT: CENTRAL_MOMENT_DP
GMODPROC CENTRAL_MOMENT: CENTRAL_MOMENT_SP
GMODPROC CENTRAL_MOMENT_VAR: CENTRAL_MOMENT_VAR_DP
GMODPROC CENTRAL_MOMENT_VAR: CENTRAL_MOMENT_VAR_SP
GMODPROC CORRELATION: CORRELATION_DP
GMODPROC CORRELATION: CORRELATION_SP
GMODPROC COVARIANCE: COVARIANCE_DP
GMODPROC COVARIANCE: COVARIANCE_SP
GMODPROC KURTOSIS: KURTOSIS_DP
GMODPROC KURTOSIS: KURTOSIS_SP
GMODPROC MEAN: MEAN_DP
GMODPROC MEAN: MEAN_SP
GMODPROC MIXED_CENTRAL_MOMENT: MIXED_CENTRAL_MOMENT_DP
GMODPROC MIXED_CENTRAL_MOMENT: MIXED_CENTRAL_MOMENT_SP
GMODPROC MIXED_CENTRAL_MOMENT_VAR: MIXED_CENTRAL_MOMENT_VAR_DP
GMODPROC MIXED_CENTRAL_MOMENT_VAR: MIXED_CENTRAL_MOMENT_VAR_SP
GMODPROC MOMENT: MOMENT_DP
GMODPROC MOMENT: MOMENT_SP
GMODPROC SKEWNESS: SKEWNESS_DP
GMODPROC SKEWNESS: SKEWNESS_SP
GMODPROC VARIANCE: VARIANCE_DP
GMODPROC VARIANCE: VARIANCE_SP
GMODPROC STDDEV: STDDEV_DP
GMODPROC STDDEV: STDDEV_SP
PROC MOMENT_SP,9,8,0,17,0: 8,0,0,0,0,0,40200,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR AVERAGE,3,,: 2,4,4,0,0,0,93,0,0,0,2
VAR VARIANCE,3,,: 2,4,4,0,0,0,193,0,0,0,2
VAR SKEWNESS,3,,: 2,4,4,0,0,0,193,0,0,0,2
VAR KURTOSIS,3,,: 2,4,4,0,0,0,193,0,0,0,2
VAR MEAN,3,,: 2,4,4,0,0,0,93,0,0,0,2
VAR STDDEV,3,,: 2,4,4,0,0,0,193,0,0,0,2
VAR ABSDEV,3,,: 2,4,4,0,0,0,93,0,0,0,2
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CENTRAL_MOMENT_VAR_SP,3,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC SKEWNESS_SP,2,8,0,17,0: 2,4,4,0,0,0,40381,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CENTRAL_MOMENT_VAR_DP,3,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC STDDEV_SP,2,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC VARIANCE_SP,2,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CORRELATION_DP,3,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Y,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC VARIANCE_DP,2,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC COVARIANCE_SP,3,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Y,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC MOMENT_DP,9,8,0,17,0: 8,0,0,0,0,0,40200,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR AVERAGE,3,,: 2,8,5,0,0,0,93,0,0,0,2
VAR VARIANCE,3,,: 2,8,5,0,0,0,193,0,0,0,2
VAR SKEWNESS,3,,: 2,8,5,0,0,0,193,0,0,0,2
VAR KURTOSIS,3,,: 2,8,5,0,0,0,193,0,0,0,2
VAR MEAN,3,,: 2,8,5,0,0,0,93,0,0,0,2
VAR STDDEV,3,,: 2,8,5,0,0,0,193,0,0,0,2
VAR ABSDEV,3,,: 2,8,5,0,0,0,93,0,0,0,2
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC COVARIANCE_DP,3,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Y,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC KURTOSIS_SP,2,8,0,17,0: 2,4,4,0,0,0,40381,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC KURTOSIS_DP,2,8,0,17,0: 2,8,5,0,0,0,40381,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC SKEWNESS_DP,2,8,0,17,0: 2,8,5,0,0,0,40381,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC MEAN_SP,2,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC MEAN_DP,2,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC MIXED_CENTRAL_MOMENT_SP,5,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR S,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CORRELATION_SP,3,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Y,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC MIXED_CENTRAL_MOMENT_DP,5,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR S,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC STDDEV_DP,2,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC ABSDEV_SP,2,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC ABSDEV_DP,2,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC MIXED_CENTRAL_MOMENT_VAR_SP,5,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Y,3,,: 2,4,4,0,1,0,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,1000,1
VAR S,3,,: 1,4,3,0,0,0,103,0,0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC AVERAGE_SP,2,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC AVERAGE_DP,2,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR DAT,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC MIXED_CENTRAL_MOMENT_VAR_DP,5,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Y,3,,: 2,8,5,0,1,0,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,1000,1
VAR S,3,,: 1,4,3,0,0,0,103,0,0,1000,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CENTRAL_MOMENT_SP,3,8,0,17,0: 2,4,4,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,4,4,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC CENTRAL_MOMENT_DP,3,8,0,17,0: 2,8,5,0,0,0,40281,2,0,40000,0
VAR X,3,,: 2,8,5,0,1,0,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR R,3,,: 1,4,3,0,0,0,103,0,0,0,1
VAR MASK,3,,: 4,4,7,0,1,0,113,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
END
