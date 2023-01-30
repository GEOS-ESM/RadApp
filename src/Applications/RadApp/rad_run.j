#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################

#@BATCH_TIME@RUN_T
#@RUN_P
#@BATCH_JOBNAME@RUN_N
#@RUN_Q
#@BATCH_GROUP
#@BATCH_JOINOUTERR
#@BATCH_NAME -o rad_run.o@RSTDATE

#######################################################################
#                         System Settings
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv ARCH `uname`

setenv SITE             @SITE
setenv GEOSDIR          @GEOSDIR
setenv GEOSBIN          @GEOSBIN
setenv GEOSETC          @GEOSETC
setenv GEOSUTIL         @GEOSSRC

@NATIVE_BUILD source $GEOSBIN/g5_modules
@NATIVE_BUILD setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BASEDIR}/${ARCH}/lib:${GEOSDIR}/lib

setenv RUN_CMD "$GEOSBIN/esma_mpirun -np "

setenv RADVER `cat $GEOSETC/.RADAPP_VERSION`
echo   VERSION: $RADVER

#######################################################################
#             Experiment Specific Environment Variables
#######################################################################


setenv  EXPID   @EXPID
setenv  EXPDIR  @EXPDIR
setenv  HOMDIR  @HOMDIR

setenv  RSTDATE @RSTDATE

#######################################################################
#                 Create Experiment Sub-Directories
#######################################################################

if (! -e $EXPDIR/restarts   ) mkdir -p $EXPDIR/restarts
if (! -e $EXPDIR/holding    ) mkdir -p $EXPDIR/holding

setenv  SCRDIR  $EXPDIR/scratch

if (! -e $SCRDIR ) mkdir -p $SCRDIR

#######################################################################
#                   Set Experiment Run Parameters
#######################################################################

set       NX  = `grep '^\s*NX:'             $HOMDIR/AGCM.rc | cut -d: -f2`
set       NY  = `grep '^\s*NY:'             $HOMDIR/AGCM.rc | cut -d: -f2`
set  AGCM_IM  = `grep '^\s*AGCM_IM:'        $HOMDIR/AGCM.rc | cut -d: -f2`
set  AGCM_JM  = `grep '^\s*AGCM_JM:'        $HOMDIR/AGCM.rc | cut -d: -f2`
set  AGCM_LM  = `grep '^\s*AGCM_LM:'        $HOMDIR/AGCM.rc | cut -d: -f2`

# Calculate number of cores/nodes for IOSERVER
# --------------------------------------------

set USE_IOSERVER      = @USE_IOSERVER
set NUM_OSERVER_NODES = `grep '^\s*IOSERVER_NODES:'  $HOMDIR/AGCM.rc | cut -d: -f2`
set NUM_BACKEND_PES   = `grep '^\s*NUM_BACKEND_PES:' $HOMDIR/AGCM.rc | cut -d: -f2`

# Check for Over-Specification of CPU Resources
# ---------------------------------------------
if ($?SLURM_NTASKS) then
   set  NCPUS = $SLURM_NTASKS
else if ($?PBS_NODEFILE) then
   set  NCPUS = `cat $PBS_NODEFILE | wc -l`
else
   set  NCPUS = NULL
endif

@ MODEL_NPES = $NX * $NY

set NCPUS_PER_NODE = @NCPUS_PER_NODE
set NUM_MODEL_NODES=`echo "scale=1;($MODEL_NPES / $NCPUS_PER_NODE)" | bc | awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {print ceil($1)}'`

if ( $NCPUS != NULL ) then

   if ( $USE_IOSERVER == 1 ) then

      @ TOTAL_NODES = $NUM_MODEL_NODES + $NUM_OSERVER_NODES

      @ TOTAL_PES = $TOTAL_NODES * $NCPUS_PER_NODE

      if( $TOTAL_PES > $NCPUS ) then
         echo "CPU Resources are Over-Specified"
         echo "--------------------------------"
         echo "Allotted  NCPUs: $NCPUS"
         echo "Requested NCPUs: $TOTAL_PES"
         echo ""
         echo "Specified NX: $NX"
         echo "Specified NY: $NY"
         echo ""
         echo "Specified model nodes: $NUM_MODEL_NODES"
         echo "Specified oserver nodes: $NUM_OSERVER_NODES"
         echo "Specified cores per node: $NCPUS_PER_NODE"
         exit
      endif

   else

      @ TOTAL_PES = $MODEL_NPES

      if( $TOTAL_PES > $NCPUS ) then
         echo "CPU Resources are Over-Specified"
         echo "--------------------------------"
         echo "Allotted  NCPUs: $NCPUS"
         echo "Requested NCPUs: $TOTAL_PES"
         echo ""
         echo "Specified NX: $NX"
         echo "Specified NY: $NY"
         echo ""
         echo "Specified model nodes: $NUM_MODEL_NODES"
         echo "Specified cores per node: $NCPUS_PER_NODE"
         exit
      endif

   endif

else
   # This is for the desktop path

   @ TOTAL_PES = $MODEL_NPES

endif

#######################################################################
#   Move to Scratch Directory and Copy RC Files from Home Directory
#######################################################################

cd $SCRDIR
/bin/rm -rf *
                             cp -f  $EXPDIR/RC/* .
                             cp     $EXPDIR/cap_restart .
                             cp -f  $HOMDIR/*.rc .
                             cp -f  $HOMDIR/*.yaml .
                             cp     $GEOSBIN/bundleParser.py .

set END_DATE  = `grep '^\s*END_DATE:'     CAP.rc | cut -d: -f2`
set NUM_SGMT  = `grep '^\s*NUM_SGMT:'     CAP.rc | cut -d: -f2`
set FSEGMENT  = `grep '^\s*FCST_SEGMENT:' CAP.rc | cut -d: -f2`
set USE_SHMEM = `grep '^\s*USE_SHMEM:'    CAP.rc | cut -d: -f2`

#######################################################################
#              Create HISTORY Collection Directories
#######################################################################

set collections = ''
foreach line ("`cat HISTORY.rc`")
   set firstword  = `echo $line | awk '{print $1}'`
   set firstchar  = `echo $firstword | cut -c1`
   set secondword = `echo $line | awk '{print $2}'`

   if ( $firstword == "::" ) goto done

   if ( $firstchar != "#" ) then
      set collection  = `echo $firstword | sed -e "s/'//g"`
      set collections = `echo $collections $collection`
      if ( $secondword == :: ) goto done
   endif

   if ( $firstword == COLLECTIONS: ) then
      set collections = `echo $secondword | sed -e "s/'//g"`
   endif
end

done:
   foreach collection ( $collections )
      if (! -e $EXPDIR/$collection )         mkdir $EXPDIR/$collection
      if (! -e $EXPDIR/holding/$collection ) mkdir $EXPDIR/holding/$collection
   end

#######################################################################
#                        Link Boundary Datasets
#######################################################################

setenv BCSDIR    @BCSDIR
setenv CHMDIR    @CHMDIR
setenv DATELINE  DC
setenv EMISSIONS @EMISSIONS

set             FILE = linkbcs
/bin/rm -f     $FILE
cat << _EOF_ > $FILE
#!/bin/csh -f

/bin/mkdir -p            ExtData
/bin/ln    -sf $CHMDIR/* ExtData

# DAS or REPLAY Mode (AGCM.rc:  pchem_clim_years = 1-Year Climatology)
# --------------------------------------------------------------------
@OPS_SPECIES/bin/ln -sf $BCSDIR/Shared/pchem.species.Clim_Prod_Loss.z_721x72.nc4 species.data

# CMIP-5 Ozone Data (AGCM.rc:  pchem_clim_years = 228-Years)
# ----------------------------------------------------------
@CMIP_SPECIES/bin/ln -sf $BCSDIR/Shared/pchem.species.CMIP-5.1870-2097.z_91x72.nc4 species.data

# S2S pre-industrial with prod/loss of stratospheric water vapor
# (AGCM.rc:  pchem_clim_years = 3-Years,  and  H2O_ProdLoss: 1 )
# --------------------------------------------------------------
#/bin/ln -sf $BCSDIR/Shared/pchem.species.CMIP-6.wH2OandPL.1850s.z_91x72.nc4 species.data

# MERRA-2 Ozone Data (AGCM.rc:  pchem_clim_years = 39-Years)
# ----------------------------------------------------------
@MERRA2OX_SPECIES/bin/ln -sf $BCSDIR/Shared/pchem.species.CMIP-5.MERRA2OX.197902-201706.z_91x72.nc4 species.data

/bin/ln -sf $BCSDIR/Shared/*bin .
/bin/ln -sf $BCSDIR/Shared/*c2l*.nc4 .

_EOF_

chmod +x linkbcs
cp  linkbcs $EXPDIR

#######################################################################
#                  Setup executable
#######################################################################

@SINGULARITY_BUILD #######################################################################
@SINGULARITY_BUILD #             Settings for Singularity - EXPERIMENTAL
@SINGULARITY_BUILD #######################################################################
@SINGULARITY_BUILD
@SINGULARITY_BUILD # Note these have only really been tested on Discover
@SINGULARITY_BUILD # and are not guaranteed to work on other systems

@SINGULARITY_BUILD # Based on work on discover, to run you need to load the same compiler
@SINGULARITY_BUILD # and MPI to match those in the container. For example, if your container was
@SINGULARITY_BUILD # built with:
@SINGULARITY_BUILD #   GNU 10.3.0
@SINGULARITY_BUILD #   Intel Fortran 2021.6.0 (aka Intel oneAPI 2022.1.0)
@SINGULARITY_BUILD #   Intel MPI 2021.6.0 (aka Intel oneAPI 2022.1.0)
@SINGULARITY_BUILD # then you would need to load:
@SINGULARITY_BUILD #   source /usr/share/modules/init/csh
@SINGULARITY_BUILD #   module purge
@SINGULARITY_BUILD #   module load comp/gcc/10.3.0
@SINGULARITY_BUILD #   module load comp/intel/2021.6.0
@SINGULARITY_BUILD #   module load mpi/impi/2021.6.0
@SINGULARITY_BUILD #
@SINGULARITY_BUILD # And then also append ${GEOSDIR}/lib to LD_LIBRARY_PATH
@SINGULARITY_BUILD #   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${GEOSDIR}/lib

@SINGULARITY_BUILD # Also look below for suggestions on Intel MPI, OpenMPI and MPT environment variables
@SINGULARITY_BUILD #
@SINGULARITY_BUILD # If you are using singularity, set the path to the singularity sandbox here
@SINGULARITY_BUILD setenv SINGULARITY_SANDBOX ""
@SINGULARITY_BUILD
@SINGULARITY_BUILD # Error out if SINGULARITY_SANDBOX is not set
@SINGULARITY_BUILD if( $SINGULARITY_SANDBOX == "" ) then
@SINGULARITY_BUILD    echo "ERROR: You must set SINGULARITY_SANDBOX to the path to your Singularity sandbox"
@SINGULARITY_BUILD    exit 1
@SINGULARITY_BUILD endif
@SINGULARITY_BUILD
@SINGULARITY_BUILD # If SINGULARITY_SANDBOX is non-empty, then run executable in singularity sandbox
@SINGULARITY_BUILD echo "We are running under Singularity"
@SINGULARITY_BUILD echo ""
@SINGULARITY_BUILD
@SINGULARITY_BUILD # Load the Singularity module
@SINGULARITY_BUILD module load singularity
@SINGULARITY_BUILD
@SINGULARITY_BUILD # Set Singularity Bind Paths. Note: These are dependent on where you are running.
@SINGULARITY_BUILD # By default, we'll assume you are running this script from NOBACKUP
@SINGULARITY_BUILD setenv SINGULARITY_BIND_PATH "-B ${NOBACKUP}:${NOBACKUP}"
@SINGULARITY_BUILD
@SINGULARITY_BUILD # If you are running from a different location, you will need to change the bind path
@SINGULARITY_BUILD # Also, note that often $NOBACKUP is, say, /discover/nobackup/username, but gcm_setup
@SINGULARITY_BUILD # will set GEOSDIR, GEOSBIN, etc. above to something like /gpfsm/dnbXX/username which
@SINGULARITY_BUILD # is the "real" physical path that /discover/nobackup/username is a symlink to.
@SINGULARITY_BUILD # You might need to change all the gpfsm paths to nobackup paths.
@SINGULARITY_BUILD
@SINGULARITY_BUILD # Set a variable to encapsulate all Singularity details
@SINGULARITY_BUILD setenv SINGULARITY_RUN "singularity exec $SINGULARITY_BIND_PATH $SINGULARITY_SANDBOX"
@SINGULARITY_BUILD
@SINGULARITY_BUILD # Detect if RadApp.x is in the experiment directory
@SINGULARITY_BUILD if (-e $EXPDIR/RadApp.x) then
@SINGULARITY_BUILD    echo "Found RadApp.x in $EXPDIR"
@SINGULARITY_BUILD
@SINGULARITY_BUILD    # If SINGULARITY_SANDBOX is non-empty and RadApp.x is found in the experiment directory,
@SINGULARITY_BUILD    # force the use of RadApp.x in the installation directory
@SINGULARITY_BUILD    if( $SINGULARITY_SANDBOX != "" ) then
@SINGULARITY_BUILD       echo "NOTE: Testing has shown Singularity only works when running with"
@SINGULARITY_BUILD       echo "      the RadApp.x executable directly from the installation bin directory"
@SINGULARITY_BUILD       echo ""
@SINGULARITY_BUILD       echo "      So, we will *ignore* the local RadApp.x and "
@SINGULARITY_BUILD       echo "      instead use $GEOSBIN/RadApp.x"
@SINGULARITY_BUILD       echo ""
@SINGULARITY_BUILD    else
@SINGULARITY_BUILD       echo "Using RadApp.x from $GEOSBIN"
@SINGULARITY_BUILD    endif
@SINGULARITY_BUILD    setenv GEOSEXE $GEOSBIN/RadApp.x
@SINGULARITY_BUILD else
@SINGULARITY_BUILD    echo "Using RadApp.x from $GEOSBIN"
@SINGULARITY_BUILD    setenv GEOSEXE $GEOSBIN/RadApp.x
@SINGULARITY_BUILD endif

@NATIVE_BUILD echo "Copying $EXPDIR/RadApp.x to $SCRDIR"
@NATIVE_BUILD echo ""
@NATIVE_BUILD /bin/cp $EXPDIR/RadApp.x $SCRDIR/RadApp.x
@NATIVE_BUILD setenv GEOSEXE $SCRDIR/RadApp.x

#######################################################################
#                         Get RESTARTS
#######################################################################

set rst_files      = `grep "RESTART_FILE"    AGCM.rc | grep -v VEGDYN | grep -v "#" | cut -d ":" -f1 | cut -d "_" -f1-2`
set rst_file_names = `grep "RESTART_FILE"    AGCM.rc | grep -v VEGDYN | grep -v "#" | cut -d ":" -f2`

set chk_files      = `grep "CHECKPOINT_FILE" AGCM.rc | grep -v "#" | cut -d ":" -f1 | cut -d "_" -f1-2`
set chk_file_names = `grep "CHECKPOINT_FILE" AGCM.rc | grep -v "#" | cut -d ":" -f2`

set monthly_chk_names = `cat $EXPDIR/HISTORY.rc | grep -v '^[\t ]*#' | sed -n 's/\([^\t ]\+\).monthly:[\t ]*1.*/\1/p' | sed 's/$/_rst/' `

# Remove possible bootstrap parameters (+/-)
# ------------------------------------------
set dummy = `echo $rst_file_names`
set rst_file_names = ''
foreach rst ( $dummy )
  set length  = `echo $rst | awk '{print length($0)}'`
  set    bit  = `echo $rst | cut -c1`
  if(  "$bit" == "+" | \
       "$bit" == "-" ) set rst = `echo $rst | cut -c2-$length`
  set rst_file_names = `echo $rst_file_names $rst`
end

# Copy Restarts to Scratch Directory
# ----------------------------------
foreach rst ( $rst_file_names $monthly_chk_names )
   if(-e $EXPDIR/$rst ) cp $EXPDIR/$rst . &
end
wait

# Copy and Tar Initial Restarts to Restarts Directory
# ---------------------------------------------------
set edate = e`cat cap_restart | cut -c1-8`_`cat cap_restart | cut -c10-11`z
set numrs = `/bin/ls -1 ${EXPDIR}/restarts/*${edate}* | wc -l`
if($numrs == 0) then
   foreach rst ( $rst_file_names )
      if( -e $rst & ! -e ${EXPDIR}/restarts/$EXPID.${rst}.${edate}.${RADVER} ) then
            cp $rst ${EXPDIR}/restarts/$EXPID.${rst}.${edate}.${RADVER} &
      endif
   end
   wait
   cd $EXPDIR/restarts
     /bin/rm -rf `/bin/ls -d -1     $EXPID.*.${edate}.${RADVER}`
   cd $SCRDIR
endif

# If any restart is binary, set NUM_READERS to 1 so that
# +-style bootstrapping of missing files can occur in
# MAPL. pbinary cannot do this, but pnc4 can.
# ------------------------------------------------------
set found_binary = 0

foreach rst ( $rst_file_names )
   if (-e $rst) then
      set rst_type = `/usr/bin/file -Lb --mime-type $rst`
      if ( $rst_type =~ "application/octet-stream" ) then
         set found_binary = 1
      endif
   endif
end

if ($found_binary == 1) then
   /bin/mv AGCM.rc AGCM.tmp
   cat AGCM.tmp | sed -e "/^NUM_READERS/ s/\([0-9]\+\)/1/g" > AGCM.rc
   /bin/rm AGCM.tmp
endif


##################################################################
######
######         Perform multiple iterations of Model Run
######
##################################################################

@ counter    = 1
while ( $counter <= ${NUM_SGMT} )

/bin/rm -f  EGRESS

cp -f $HOMDIR/CAP.rc .

/bin/mv CAP.rc CAP.rc.orig
awk '{$1=$1};1' < CAP.rc.orig > CAP.rc

# Set Time Variables for Current_(c), Ending_(e), and Segment_(s) dates
# ---------------------------------------------------------------------
set nymdc = `awk '{print $1}' cap_restart`
set nhmsc = `awk '{print $2}' cap_restart`
set nymde = `grep '^\s*END_DATE:' CAP.rc | cut -d: -f2 | awk '{print $1}'`
set nhmse = `grep '^\s*END_DATE:' CAP.rc | cut -d: -f2 | awk '{print $2}'`
set nymds = `grep '^\s*JOB_SGMT:' CAP.rc | cut -d: -f2 | awk '{print $1}'`
set nhmss = `grep '^\s*JOB_SGMT:' CAP.rc | cut -d: -f2 | awk '{print $2}'`

# Compute Time Variables at the Finish_(f) of current segment
# -----------------------------------------------------------
set nyear   = `echo $nymds | cut -c1-4`
set nmonth  = `echo $nymds | cut -c5-6`
set nday    = `echo $nymds | cut -c7-8`
set nhour   = `echo $nhmss | cut -c1-2`
set nminute = `echo $nhmss | cut -c3-4`
set nsec    = `echo $nhmss | cut -c5-6`
       @ dt = $nsec + 60 * $nminute + 3600 * $nhour + 86400 * $nday

set nymdf = $nymdc
set nhmsf = $nhmsc
set date  = `$GEOSBIN/tick $nymdf $nhmsf $dt`
set nymdf =  $date[1]
set nhmsf =  $date[2]
set year  = `echo $nymdf | cut -c1-4`
set month = `echo $nymdf | cut -c5-6`
set day   = `echo $nymdf | cut -c7-8`

     @  month = $month + $nmonth
while( $month > 12 )
     @  month = $month - 12
     @  year  = $year  + 1
end
     @  year  = $year  + $nyear
     @ nymdf  = $year * 10000 + $month * 100 + $day

if( $nymdf >  $nymde )    set nymdf = $nymde
if( $nymdf == $nymde )    then
    if( $nhmsf > $nhmse ) set nhmsf = $nhmse
endif

set yearc = `echo $nymdc | cut -c1-4`
set yearf = `echo $nymdf | cut -c1-4`

# Which ExtData are we using
set  EXTDATA2G_TRUE = `grep -i '^\s*USE_EXTDATA2G:\s*\.TRUE\.'    CAP.rc | wc -l`

# Select proper AMIP GOCART Emission RC Files
# -------------------------------------------
if( ${EMISSIONS} == AMIP_EMISSIONS ) then
    if( $EXTDATA2G_TRUE == 0 ) then
       set AMIP_Transition_Date = 20000301

       # Before 2000-03-01, we need to use AMIP.20C which has different
       # emissions (HFED instead of QFED) valid before 2000-03-01. Note
       # that if you make a change to anything in $EXPDIR/RC/AMIP or
       # $EXPDIR/RC/AMIP.20C, you might need to make a change in the other
       # directory to be consistent. Some files in AMIP.20C are symlinks to
       # that in AMIP but others are not.

       if( $nymdc < ${AMIP_Transition_Date} ) then
            set AMIP_EMISSIONS_DIRECTORY = $EXPDIR/RC/AMIP.20C
            if( $nymdf > ${AMIP_Transition_Date} ) then
             set nymdf = ${AMIP_Transition_Date}
             set oldstring = `grep '^\s*END_DATE:' CAP.rc`
             set newstring = "END_DATE: $nymdf $nhmsf"
             /bin/mv CAP.rc CAP.tmp
                        cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
            endif
       else
            set AMIP_EMISSIONS_DIRECTORY = $EXPDIR/RC/AMIP
       endif
    else
       set AMIP_EMISSIONS_DIRECTORY = $EXPDIR/RC/AMIP
    endif

    if( $AGCM_LM == 72 ) then
        cp ${AMIP_EMISSIONS_DIRECTORY}/*.rc .
        cp ${AMIP_EMISSIONS_DIRECTORY}/*.yaml .
    else
        set files = `/bin/ls -1 ${AMIP_EMISSIONS_DIRECTORY}/*.rc ${AMIP_EMISSIONS_DIRECTORY}/*.yaml`
        foreach file ($files)
          /bin/rm -f `basename $file`
          /bin/rm -f dummy
          cp $file dummy
          cat dummy | sed -e "s|/L72/|/L${AGCM_LM}/|g" | sed -e "s|z72|z${AGCM_LM}|g" > `basename $file`
        end
    endif

endif

# Rename big ExtData files that are not needed
# --------------------------------------------
set            SC_TRUE = `grep -i '^\s*ENABLE_STRATCHEM:\s*\.TRUE\.'     GEOS_ChemGridComp.rc | wc -l`
if (          $SC_TRUE == 0 && -e StratChem_ExtData.rc          ) /bin/mv          StratChem_ExtData.rc          StratChem_ExtData.rc.NOT_USED
set           GMI_TRUE = `grep -i '^\s*ENABLE_GMICHEM:\s*\.TRUE\.'       GEOS_ChemGridComp.rc | wc -l`
if (         $GMI_TRUE == 0 && -e GMI_ExtData.rc                ) /bin/mv                GMI_ExtData.rc                GMI_ExtData.rc.NOT_USED
set           GCC_TRUE = `grep -i '^\s*ENABLE_GEOSCHEM:\s*\.TRUE\.'      GEOS_ChemGridComp.rc | wc -l`
if (         $GCC_TRUE == 0 && -e GEOSCHEMchem_ExtData.rc       ) /bin/mv       GEOSCHEMchem_ExtData.rc       GEOSCHEMchem_ExtData.rc.NOT_USED
set         CARMA_TRUE = `grep -i '^\s*ENABLE_CARMA:\s*\.TRUE\.'         GEOS_ChemGridComp.rc | wc -l`
if (       $CARMA_TRUE == 0 && -e CARMAchem_GridComp_ExtData.rc ) /bin/mv CARMAchem_GridComp_ExtData.rc CARMAchem_GridComp_ExtData.rc.NOT_USED
set           DNA_TRUE = `grep -i '^\s*ENABLE_DNA:\s*\.TRUE\.'           GEOS_ChemGridComp.rc | wc -l`
if (         $DNA_TRUE == 0 && -e DNA_ExtData.rc                ) /bin/mv                DNA_ExtData.rc                DNA_ExtData.rc.NOT_USED
set         ACHEM_TRUE = `grep -i '^\s*ENABLE_ACHEM:\s*\.TRUE\.'         GEOS_ChemGridComp.rc | wc -l`
if (       $ACHEM_TRUE == 0 && -e GEOSachem_ExtData.rc          ) /bin/mv          GEOSachem_ExtData.rc          GEOSachem_ExtData.rc.NOT_USED
set   GOCART_DATA_TRUE = `grep -i '^\s*ENABLE_GOCART_DATA:\s*\.TRUE\.'   GEOS_ChemGridComp.rc | wc -l`
if ( $GOCART_DATA_TRUE == 0 && -e GOCARTdata_ExtData.rc         ) /bin/mv         GOCARTdata_ExtData.rc         GOCARTdata_ExtData.rc.NOT_USED

# Generate the complete ExtData.rc
# --------------------------------
if(-e ExtData.rc )    /bin/rm -f   ExtData.rc
set  extdata_files = `/bin/ls -1 *_ExtData.rc`

# Switch to MODIS v6.1 data after Nov 2021
if( $EXTDATA2G_TRUE == 0 ) then
   set MODIS_Transition_Date = 20211101
   if ( ${EMISSIONS} == OPS_EMISSIONS && ${MODIS_Transition_Date} <= $nymdc ) then
       cat $extdata_files | sed 's|\(qfed2.emis_.*\).006.|\1.061.|g' > ExtData.rc
   else
   cat $extdata_files > ExtData.rc
   endif
endif

if( $EXTDATA2G_TRUE == 1 ) then

  $GEOSBIN/construct_extdata_yaml_list.py GEOS_ChemGridComp.rc
  touch ExtData.rc

endif

# Move GOCART to use RRTMGP Bands
# -------------------------------
# UNCOMMENT THE LINES BELOW IF RUNNING RRTMGP
#
#set instance_files = `/bin/ls -1 *_instance*.rc`
#foreach instance ($instance_files)
#   /bin/mv $instance $instance.tmp
#   cat $instance.tmp | sed -e '/RRTMG/ s#RRTMG#RRTMGP#' > $instance
#   /bin/rm $instance.tmp
#end

# Link Boundary Conditions for Appropriate Date
# ---------------------------------------------
setenv YEAR $yearc
./linkbcs

# Check for MERRA2OX Consistency
# ------------------------------

# The MERRA2OX pchem file is only valid until 201706, so this is a first
# attempt at a check to make sure you aren't using it and are past the date

# Check for MERRA2OX by looking at AGCM.rc
set PCHEM_CLIM_YEARS = `awk '/pchem_clim_years/ {print $2}' AGCM.rc`

# If it is 39, we are using MERRA2OX
if ( $PCHEM_CLIM_YEARS == 39 ) then

   # Grab the date from cap_restart
   set YEARMON = `cat cap_restart | cut -c1-6`

   # Set a magic date
   set MERRA2OX_END_DATE = "201706"

   # String comparison seems to work here...
   if ( $YEARMON > $MERRA2OX_END_DATE ) then
      echo "You seem to be using MERRA2OX pchem species file, but your simulation date [${YEARMON}] is after 201706. This file is only valid until this time."
      exit 2
   endif
endif

# Environment variables for MPI, etc
# ----------------------------------

@SETENVS

# Run bundleParser.py
#---------------------
python bundleParser.py

# Establish safe default number of OpenMP threads
# -----------------------------------------------
setenv OMP_NUM_THREADS 1

# Run RadApp.x
# -------------
if( $USE_SHMEM == 1 ) $GEOSBIN/RmShmKeys_sshmpi.csh >& /dev/null

if( $USE_IOSERVER == 1 ) then
   set IOSERVER_OPTIONS = "--npes_model $MODEL_NPES --nodes_output_server $NUM_OSERVER_NODES"
   set IOSERVER_EXTRA   = "--oserver_type multigroup --npes_backend_pernode $NUM_BACKEND_PES"
else
   set IOSERVER_OPTIONS = ""
   set IOSERVER_EXTRA   = ""
endif

@SINGULARITY_BUILD $RUN_CMD $TOTAL_PES $SINGULARITY_RUN $GEOSEXE $IOSERVER_OPTIONS $IOSERVER_EXTRA --logging_config 'logging.yaml'
@NATIVE_BUILD $RUN_CMD $TOTAL_PES $GEOSEXE $IOSERVER_OPTIONS $IOSERVER_EXTRA --logging_config 'logging.yaml'

if( $USE_SHMEM == 1 ) $GEOSBIN/RmShmKeys_sshmpi.csh >& /dev/null

if( -e EGRESS ) then
   set rc = 0
else
   set rc = -1
endif
echo RadApp Run Status: $rc
if ( $rc == -1 ) exit -1

#######################################################################
#   Rename Final Checkpoints => Restarts for Next Segment and Archive
#        Note: cap_restart contains the current NYMD and NHMS
#######################################################################

set edate  = e`awk '{print $1}' cap_restart`_`awk '{print $2}' cap_restart | cut -c1-2`z

# Move Intermediate Checkpoints to RESTARTS directory
# ---------------------------------------------------
set   checkpoints  =    `/bin/ls -1 *_checkpoint.*`
if( $#checkpoints != 0 ) /bin/mv -f *_checkpoint.* ${EXPDIR}/restarts


# Rename Final Checkpoints for Archive
# ------------------------------------
    set checkpoints = `/bin/ls -1 *_checkpoint`
foreach checkpoint ($checkpoints)
        set   chk_type = `/usr/bin/file -Lb --mime-type $checkpoint`
            if ( $chk_type =~ "application/octet-stream" ) then
                  set ext  = bin
            else
                  set ext  = nc4
            endif
       /bin/mv            $checkpoint      $EXPID.${checkpoint}.${edate}.${RADVER}.$ext
       $GEOSBIN/stripname _checkpoint _rst $EXPID.${checkpoint}.${edate}.${RADVER}.$ext
end


# Remove Initial RESTARTS
# -----------------------
set restarts = `/bin/ls -1 *_rst`
/bin/rm  $restarts


# Copy Renamed Final Checkpoints to RESTARTS directory
# ----------------------------------------------------
    set  restarts = `/bin/ls -1 $EXPID.*_rst.${edate}.${RADVER}.*`
foreach  restart ($restarts)
cp $restart ${EXPDIR}/restarts
end

# Remove EXPID from RESTART name
# ------------------------------
    set  restarts = `/bin/ls -1 $EXPID.*_rst.${edate}.${RADVER}.*`
foreach  restart ($restarts)
$GEOSBIN/stripname $EXPID. '' $restart
end

# Remove DATE and VERSION Stamps from RESTART name
# ------------------------------------------------
    set  restarts = `/bin/ls -1 *_rst.${edate}.${RADVER}.*`
foreach  restart ($restarts)
$GEOSBIN/stripname .${edate}.${RADVER}.\* '' $restart
end


# TAR ARCHIVED RESTARTS
# ---------------------
cd $EXPDIR/restarts
if( $FSEGMENT == 00000000 ) then
     /bin/rm -rf `/bin/ls -d -1     $EXPID.*.${edate}.${RADVER}.*`
endif


#######################################################################
#               Move HISTORY Files to Holding Directory
#######################################################################

# Move current files to /holding
# ------------------------------
cd $SCRDIR
foreach collection ( $collections )
   /bin/mv `/bin/ls -1 *.${collection}.*` $EXPDIR/holding/$collection
end

#######################################################################
#                         Update Iteration Counter
#######################################################################

set enddate = `echo  $END_DATE | cut -c1-8`
set capdate = `cat cap_restart | cut -c1-8`

if ( $capdate < $enddate ) then
@ counter = $counter    + 1
else
@ counter = ${NUM_SGMT} + 1
endif

end   # end of segment loop; remain in $SCRDIR

#######################################################################
#                              Re-Submit Job
#######################################################################

foreach rst ( `/bin/ls -1 *_rst` )
   /bin/rm -f $EXPDIR/$rst
end
   /bin/rm -f $EXPDIR/cap_restart
foreach rst ( `/bin/ls -1 *_rst` )
   cp $rst $EXPDIR/$rst &
end
wait
cp cap_restart $EXPDIR/cap_restart

if ( $rc == 0 ) then
      cd  $HOMDIR
          if( $capdate < $enddate ) @BATCH_CMD $HOMDIR/rad_run.j
      endif
endif
