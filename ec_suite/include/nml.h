cat > main.dat << 'EOF'
! Emacs: -*- mode: f90 -*-
!***DROUGHT-INDEX        *** Main configuration file ***
&mainconfig
 !
 ! directory iti which output will be written
 outpath        = "./"
 !
 ! file and variable name for soil moisture fields
 ! if read_opt_h is False, kernel density estimator (kde) is fitted and calculated for these values
 ! if read_opt_h is True, kde is taken from opt_h_file and SMI is calculated for the SM values
 soilmoist_file = "%SOILMOIST_FILE%"
 SM_vname       = "SM_Lall"
 !
 ! SMI time steps per year
 ! 12          == SM is monthly
 ! 365         == SM is daily
 nCalendarStepsYear = 365
 ! lag to calculate moving average for inputs
 ! moving average is calculated for time steps preceeding the current time step
 !
 ! Examples:
 ! 0, 7, 31    ==  for daily values (no lag, weekly, monthly)
 lag = 0

 ! -----------------------------------------------------------------------------
 !
 ! Apply silverman rule to estimate kernel width
 !
 ! Determination how to estimate the kernel width for the kde:
 ! silverman_h == .TRUE.  - 'Silverman rule' is applied for determination of the kernel width - no optimization
 ! silverman_h == .FALSE. - optimization of the kernel width 
 silverman_h    = .True.

 ! -----------------------------------------------------------------------------
 !
 ! Use previously calculated kde
 !
 ! if  .TRUE. kernel widths h and soil moisture sm are read from file "opt_h_file with"
 ! if .FALSE. kernel widths are optimized or determined by 'Silverman rule', "opt_h_file" and "opt_h_vname" are ignored
 read_opt_h     = .FALSE.
 opt_h_file     = "XXX"
 opt_h_vname    = "kernel_width"
 opt_h_sm_vname = "SM"

 ! -----------------------------------------------------------------------------
 !
 ! External SMI
 !
 ! To invert SMI values and for the clustering, SMI values can be provided
 ! in an external file. For the clustering, this deactivates the SMI estimation,
 ! which saves some time.
 ext_smi      = .FALSE.
 ext_smi_file = 'XXX'

 ! -----------------------------------------------------------------------------
 !
 ! Invert SMI cdf
 !
 ! If .True., then SMI values are read from ext_smi_file. These are transformed
 ! to soil moisture saturation using the cdf created by soilmoist_file.
 invert_SMI = .False.

 ! -----------------------------------------------------------------------------
 !
 ! Cluster Analysis
 ! 
 ! if .TRUE. cluster analysis should be performed
 do_cluster         = .False.
 ! SMI - treshhold of for defining drought conditions
 SMI_thld            = 0.20
 ! cellsize (edge length) of the grid im [km]
 cellsize            = 4.0
 ! threshold for cluster formation in space cellsize**2 * thCellClus ~ 640 km2
 thCellClus          = 40
 ! number cells for joining clusters in time ~ 6400
 nCellInter          = 400
 ! number of cells per area interval
 deltaArea           = 20

 ! -----------------------------------------------------------------------------
 !
 ! Mask file
 ! 
 ! mask is OPTIONAL
 ! mask to define the domain for which the SMI calculation or clustering is performed
 maskfName      = "XXX"
 mask_vname     = "XXX"

 ! if .TRUE. SMI is evaluated for subregions (e.g. basins) -
 ! ****CURRENTLY NOT WORKING
 basin_flag     = .FALSE.       
 basinfName     = "./test_smi/basin.nc"
 basin_vname    = "basins"
! ****CURRENTLY NOT WORKING
/ !end
EOF
