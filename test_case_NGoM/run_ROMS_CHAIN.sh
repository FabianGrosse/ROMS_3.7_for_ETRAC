#!/bin/bash
# =====================================================================
# script for running a ROMS simulation as a series of subsequent jobs
# by Fabian Grosse (fabian.grosse[at]dal.ca)
# =====================================================================
# relevant functions:
# - two different time modes
#   1)   "YEAR" mode: run series of subsequent years, one job per year
#   2) "NTIMES" mode: run a fixed number of time steps, use fixed
#                     number of time steps per job
#
# - automatic continuation if job was cancelled due to cluster failure
#
# - email notifications are send to user if:
#   1) a (sub-)job fails before entering time loop => job abortion
#   2) a simulation blows up                       => job abortion
#   3) a (sub-)job failed due to cluster failure   => auto-continuation
# =====================================================================
# HOW TO USE?
# 1) define set-up in script header
# 2) execute script by: nohup ./run_ROMS_CHAIN.sh &
# =====================================================================
# last changed: February 08th, 2018
# =====================================================================
# USER-DEFINED SET-UP
# =====================================================================

# debug this script? yes (1) or no (0 = default)?
# - debug=1: script log is prompted to terminal
# - debug=0: script log is written to log file
debug=0

# define simulation ID (used for names of ROMS output files)
SIMID=MCH-TEST

# define if job is a new job, i.e., starting from initialisation: yes (1) or no (0)?
# If no, you need to provide a restart file complying with your choice of the PERFECT_RESTART option
newJob=1

# submit job (1) or not (0)?
# Deactivation of job submission is recommended for testing of a new setup.
# Activation is required for multi-job simulations.
submitJob=0

# set directories for input, output and temporary files
inputPath=/scratch/grosse/infiles
outputPath=/scratch/grosse/roms854_${SIMID}_CHAIN-JOB
tmpPath=/scratch/grosse/roms854_${SIMID}_TMP

# set files used for model setup
oceanBase=ocean_CHAIN-JOB.in
bioFile=bio_Fennel_mch078_ETRAC.in
varFile=varinfo_Fennel_PO4_RDON_ETRAC.dat

# set files used for model initialisation (initial data and restart) 
iniFile=${inputPath}/mch_ini20000101_049.nc
rstFile=${outputPath}/rst_${SIMID}.nc

# define filename for log file
logFile=${SIMID}.log

# keep log files of individual jobs? yes (1) or no (0)
keepLog=1

# enable (1) or disable (0) automatic time limit adaptation, i.e., increase job time limit if first job runs out of time
autoLimit=0

# keep last restart file of individual jobs? yes (1) or no (0)
keepRST=1

# define domain tiling
NtileI=4
NtileJ=4

# define time step (in seconds) and number of time steps between restart writes
DT=60
N_RST=1440

# select time mode
# 1 => run complete years for defined period (FIRSTYEAR to LASTYEAR)
# 2 => run defined number of NTIMES (recommended for testing purposes)
timeMode=2

# if timeMode==1: provide start and end year
FIRSTYEAR=2000
LASTYEAR=2011

# if timeMode==2: provide total NTIMES and maximum NTIMES for single job
NTIMES_TOT=7200    # 7200 => 5 days (for dt=60s); 8942400 => 17 years (2000-2016; for dt=60s)
NTIMES_MAX=7200    # 527040 => 366 days (for dt=60s)

# define SLURM job settings
# (account name, memory per CPU, job name, email of user
JOB_ACC='def-kfennel'
MEM_CPU='6900M'
JOB_NAME=${SIMID}_CHAIN
USER_EMAIL='fabian.grosse@dal.ca'
# time limit: days, hours, minutes and seconds; individual 2-digit strings
JOB_DD='00' # days
JOB_HH='02' # hours
JOB_MM='59' # minutes
JOB_SS='00' # seconds

# ===================================================================
# END OF USER-DEFINED SET-UP
# ===================================================================
# AUTOMATIC JOB SETUP AND SUBMISSION
# ===================================================================

# filename of this script's log
scriptLog=CHAIN_LOG_${SIMID}.log

# set number of tasks
let NTASKS=NtileI*NtileJ

# define default job time string and value in seconds
JOB_TIME_STR=${JOB_DD}-${JOB_HH}:${JOB_MM}:${JOB_SS}
# ('10#' required before each variable to avoid bash errors due to leading zeros in strings)
JOB_TIME_SEC=$((10#${JOB_DD}*86400 + 10#${JOB_HH}*3600 + 10#${JOB_MM}*60 + 10#${JOB_SS}))

# calculate overall number of jobs
if [ ${timeMode} -eq 1 ]; then
  NTIMES_LAST=0
  let NJOBS=LASTYEAR-FIRSTYEAR+1
else
  if [ ${timeMode} -eq 2 ]; then
    let NJOBS=${NTIMES_TOT}/${NTIMES_MAX}+1
    let NTIMES_LAST=$((NTIMES_TOT-(NJOBS-1)*NTIMES_MAX))
    if [ $NTIMES_LAST -eq 0 ]; then
       let NJOBS=NJOBS-1
    fi
  else
    echo "Invalid time mode - choose 1 (full years) or 2 (defined NTIMES)."
    exit
  fi
fi

# get digits of maximum job number
if [ ${NJOBS} -lt 10 ]; then
  JOBDIG=1
else
  if [ ${NJOBS} -lt 100 ]; then
    JOBDIG=2
  else
    if [ ${NJOBS} -lt 1000 ]; then
      JOBDIG=3
    else
      echo "Maximum number of consecutive jobs for one simulation is 999."
      exit
    fi
  fi
fi

# create output directory
mkdir -p ${outputPath}

# start logging 
rm -f ${scriptLog}
head="=======================\n CHAIN JOB INFORMATION\n======================="
msg1=">> ROMS setup <<\n          SIMID = ${SIMID}\n         newJob = ${newJob}\n     input path = ${inputPath}\n    output path = ${outputPath}\n temporary path = ${tmpPath}\n       log file = ${logFile}"
msg2=">> chain job setup <<\n"
if [ ${debug} -eq 1 ]; then
  echo -e ${head}
  echo -e ${msg1}
  if [ ${newJob} -eq 1 ]; then
    echo -e " initialisation = ${iniFile}\n"
  else
    echo -e "\n"
  fi
  echo -e ${msg2}
else
  echo -e ${head} > ${scriptLog}
  echo -e ${msg1} >> ${scriptLog}
  if [ ${newJob} -eq 1 ]; then
    echo -e " initialisation = ${iniFile}\n" >> ${scriptLog}
  else
    echo -e "\n" >> ${scriptLog}
  fi
  echo -e ${msg2} >> ${scriptLog}
fi

# update non-time-varying placeholders in ocean.in base file and store as new file
awk '{gsub("INPUTPATH_TMP", "'${inputPath}'"); \
      gsub("OUTPUTPATH_TMP", "'${outputPath}'"); \
      gsub("RSTNAME_TMP", "'${rstFile}'"); \
      gsub("VARNAME_TMP", "'${varFile}'"); \
      gsub("BPARNAM_TMP", "'${bioFile}'"); \
      gsub("SIMID_TMP", "'${SIMID}'"); \
      gsub("NTILEI_TMP", "'${NtileI}'"); \
      gsub("NTILEJ_TMP", "'${NtileJ}'"); \
      gsub("DT_TMP", "'${DT}'"); \
      gsub("NRST_TMP", "'${N_RST}'"); \
      print}' ${oceanBase} > ocean_${SIMID}_BASE.in

# update non-time-varying placeholders in SLURM run script (run_ROMS_SINGLE_BASE.sh) and store as new file
awk '{gsub("JOBACC_TMP", "'${JOB_ACC}'"); \
      gsub("NTASKS_TMP", "'${NTASKS}'"); \
      gsub("MEM_CPU_TMP", "'${MEM_CPU}'"); \
      gsub("JOBNAME_TMP", "'${JOB_NAME}'"); \
      gsub("USER_EMAIL_TMP", "'${USER_EMAIL}'"); \
      print}' run_ROMS_SINGLE_BASE.sh > run_ROMS_${SIMID}_BASE.sh


# loop over jobs: placeholders in ocean.in base file are updated and jobs are submitted successively
projPath=`pwd`
IJOB=1
jobOK=1
RSTcount=0
NTIMES_LEFT=0
TIME_FAC=1.0

rm -rf $tmpPath

while [ ${IJOB} -le ${NJOBS} ]; do
  
  # set NTIMES of job
  # if previous job succeeded, set NTIMES and JOB_TIME to original value
  if [ ${jobOK} -eq 1 ]; then
    if [ ${timeMode} -eq 1 ]; then
      let YEAR=FIRSTYEAR+IJOB-1
      # leap year check
      if [ `expr ${YEAR} % 400` -eq 0 ]; then
        let NTIMES=366*86400/DT
      else
        if [ `expr ${YEAR} % 4` -eq 0 ] && [ `expr ${YEAR} % 100` -ne 0 ]; then
          let NTIMES=366*86400/DT
        else
          let NTIMES=365*86400/DT
        fi
      fi
      JOB_TIME=${JOB_TIME_STR}
      if [ ${debug} -eq 1 ]; then
        echo $IJOB $YEAR $NTIMES
      fi
    else
      # all jobs use NTIMES_MAX, except for last one
      if [ ${IJOB} -lt ${NJOBS} ]; then
        NTIMES=${NTIMES_MAX}
      else
        NTIMES=${NTIMES_LAST}
	# calculate time limit for last job
	# (use 'bc' for floating point treatment and remove decimal digits with 'cut')
	JOB_TIME_SEC=$(bc <<< "scale=5; $NTIMES_LAST/$NTIMES_MAX*$JOB_TIME_SEC" | cut -d '.' -f 1)
        let TIME_DAY=JOB_TIME_SEC/86400               # get number of day
	let TIME_SEC=JOB_TIME_SEC-TIME_DAY*86400+3600 # update total seconds (add extra hour)
	JOB_TIME_STR="$(printf "%02d" $TIME_DAY)-$(date -u -d "@${TIME_SEC}" "+%H:%M:%S")"
      fi
      JOB_TIME=${JOB_TIME_STR}
      if [ ${debug} -eq 1 ]; then
        echo $IJOB $NTIMES $JOB_TIME
      fi
    fi
  else
    # update time limit for continuation of incomplete
    TIME_SEC=$(bc <<< "scale=5; $NTIMES_LEFT/$NTIMES*$JOB_TIME_SEC" | cut -d '.' -f 1)
    let TIME_DAY=TIME_SEC/86400               # get number of day
    let TIME_SEC=TIME_SEC-TIME_DAY*86400+3600 # update total seconds (add extra hour)
    JOB_TIME="$(printf "%02d" $TIME_DAY)-$(date -u -d "@${TIME_SEC}" "+%H:%M:%S")"
    let NTIMES=NTIMES_LEFT
  fi
  # write job time limit to log file
  if [ ${debug} -eq 0 ]; then
    echo -e " job: ${IJOB}-${RSTcount} -> time limit: ${JOB_TIME}\n" >> ${scriptLog}
  fi
  # set file name of ocean.in and run script files
  if [ ${JOBDIG} -eq 1 ]; then
    oceanFile=ocean_${SIMID}_${IJOB}.in
    runFile=run_ROMS_${SIMID}_${IJOB}.sh
  else
    if [ ${JOBDIG} -eq 2 ]; then
      if [ ${IJOB} -lt 10 ]; then
        oceanFile=ocean_${SIMID}_0${IJOB}.in
        runFile=run_ROMS_${SIMID}_0${IJOB}.sh
      else
        oceanFile=ocean_${SIMID}_${IJOB}.in
        runFile=run_ROMS_${SIMID}_${IJOB}.sh
      fi
    else
      if [ ${IJOB} -lt 10 ]; then
        oceanFile=ocean_${SIMID}_00${IJOB}.in
        runFile=run_ROMS_${SIMID}_00${IJOB}.sh
      else
        if [ ${IJOB} -lt 100 ]; then
          oceanFile=ocean_${SIMID}_0${IJOB}.in
          runFile=run_ROMS_${SIMID}_0${IJOB}.sh
        else
          oceanFile=ocean_${SIMID}_${IJOB}.in
          runFile=run_ROMS_${SIMID}_${IJOB}.sh
        fi
      fi
    fi
  fi
  if [ ${debug} -eq 1 ]; then
    echo ${oceanFile}
  fi
  
  # update NTIMES, ININAME and NRREC in ocean.in and store as new file for the job 
  if [ ${IJOB} -eq 1 ] && [ ${RSTcount} -eq 0 ] && [ ${newJob} -eq 1 ]; then
    # first job of new simulation uses initialisation file
    awk '{gsub("ININAME_TMP", "'${iniFile}'"); \
          gsub("NTIMES_TMP", "'${NTIMES}'"); \
          gsub("NRREC_TMP", "0"); \
          print}' ocean_${SIMID}_BASE.in > ${oceanFile}
  else
    # following jobs use restart file
    awk '{gsub("ININAME_TMP", "'${rstFile}'"); \
          gsub("NTIMES_TMP", "'${NTIMES}'"); \
          gsub("NRREC_TMP", "-1"); \
          print}' ocean_${SIMID}_BASE.in > ${oceanFile}
  fi
  
  # update ocean.in and log file names, and job time limit in SLURM job script
  awk '{gsub("OCEANFILE_TMP", "'${oceanFile}'"); \
        gsub("LOGFILE_TMP", "'${logFile}'"); \
        gsub("JOBTIME_TMP", "'${JOB_TIME}'"); \
        print}' run_ROMS_${SIMID}_BASE.sh > ${runFile}
  
  if [ ${debug} -eq 1 ]; then
    # run "fake" job
    echo "slurmJobID=\`sbatch ${runFile}\`"
    echo "Sleep 5s to mimic job ..."
    startTime=`date +%s`
    slurmJobID=`sleep 5`
    let waitTime="$(date +%s)"-startTime
    echo "Time waited: "${waitTime}"sec"
  else
    # create temporary work directory and copy relevant files
    mkdir -p $tmpPath
    cd $tmpPath
    cp -f $projPath/${runFile} .
    cp -f $projPath/oceanM .
    cp -f $projPath/${oceanFile} .
    cp -f $projPath/${bioFile} .
    cp -f $projPath/${varFile} .
    cp -f $projPath/${scriptLog} .
    # submit job and get slurm job ID
    # job is submitted first, job ID is determined afterwards
    # This is necessary to avoid problems related to slurm issues
    #let slurmJobID="$(sbatch ${runFile} | cut -d ' ' -f 4)"
    if [ $submitJob -eq 1 ]; then
       sbatch ${runFile}
    else
       echo "Simulations files prepared and copied to ${tmpPath}."
       echo "Use ${tmpPath}/${runFile} for manual job submission."
       exit
    fi
    # sleep until job has started (i.e., until "jobStart" file exists)
    running=0
    while [ ${running} -eq 0 ]; do
      if [ -e jobStart ]; then
        running=1
      else
        sleep 300
      fi
    done
    # JOB STARTED => get job ID, remaining job time and start time
    gotInfo=0
    while [ ${gotInfo} -eq 0 ]; do
      squeueOut=`squeue --name=${JOB_NAME} -u ${USER} -o %A,%L,%S,%u | tail -n1`
      gotInfo="$(echo "${squeueOut}" | grep ${USER} | wc -l)"
      if [ ${gotInfo} -eq 1 ]; then
        slurmJobID=`echo "${squeueOut}" | cut -d ',' -f 1`
        leftTimeStr=`echo "${squeueOut}" | cut -d ',' -f 2`
        startTimeStr=`echo "${squeueOut}" | cut -d ',' -f 3`
      else
        sleep 300
      fi
    done
    # write start time to log file
    echo -e "Job ${slurmJobID} started: ${startTimeStr}" >> ${scriptLog}
    # get remaining days, hours, minutes and seconds
    isDays="$(echo "${leftTimeStr}" | grep - | wc -l)"
    if [ ${isDays} -eq 1 ]; then
      leftDays="$(echo "${leftTimeStr}" | cut -d '-' -f 1)"
      leftTimeStr="$(echo "${leftTimeStr}" | cut -d '-' -f 2)"
    else
      leftDays=0
    fi
    leftHours="$(echo "${leftTimeStr}" | cut -d ':' -f 1)"
    leftMins="$(echo "${leftTimeStr}" | cut -d ':' -f 2)"
    leftSecs="$(echo "${leftTimeStr}" | cut -d ':' -f 3)"
    # set run time counter to remaining job time in seconds
    # ('10#' required before each variable to avoid bash errors due to leading zeros)
    runTime=$((10#${leftDays}*86400+10#${leftHours}*3600+${leftMins}*60+10#${leftSecs}))
    # sleep until job has finished (i.e. until "jobEnd" file exists or "time limit + 2min" is exceeded)
    while [ ${running} -eq 1 ]; do
      sleep 300
      let runTime=runTime-60
      if [ -e jobEnd ] || [ ${runTime} -lt -120 ]; then
        running=0
        echo -e "Job ${slurmJobID} ended:   `date`" >> ${scriptLog}
      fi
    done
    # JOB ENDED => copy files to project directory
    cd $projPath
    mv -f $tmpPath/${logFile} .
    mv -f $tmpPath/${scriptLog} .
    if [ ${keepLog} -eq 1 ]; then
      cp -f  ${logFile} ${logFile}_${slurmJobID}
    fi
    mv $tmpPath/${JOB_NAME}-${slurmJobID}* .
    rm -rf $tmpPath
  fi
  
  # after job ended: check log file
  # => get expected ending time step of simulation from log file, and check if job reached this time
  if [ "$(grep ROMS\/TOMS:\ started\ time-stepping ${logFile} | wc -l)" -eq 1 ]; then
    # line with time range exists in log file => steps in next line:
    #  1) get the line and extract the time range (in the 5th column separated by colon) 
    #  2) remove all spaces and closing parenthesis => [start]-[end]
    #  3) get ending time and remove leading zeros
    expectedEndTime="$(grep ROMS\/TOMS:\ started\ time-stepping ${logFile} | cut -d ':' -f 5 | tr -d '[:space:]' | tr -d '[)]' | cut -d'-' -f 2 | sed -e 's/^0*//')"
    timeCount="$(grep \ ${expectedEndTime}\  ${logFile} | wc -l)"
    # if "timeCount" greater or equal 1 (i.e., 1+ lines with ending time stamp), end was reached
    # => check if rest of application completed successfully
    isDone="$(tail -n1 ${logFile} | grep ROMS\/TOMS:\ DONE\.\.\. | wc -l)"
    if [ ${timeCount} -gt 1 ] && [ ${isDone} -eq 1 ]; then
      jobOK=1
    else
      jobOK=0
    fi
  else
    # simulation did not even start/enter time loop=> CHECK YOUR SIMULATION!
    # send notification email to user and abort script => no further job submission
    if [ ${debug} -eq 1 ]; then
      echo -e "${SIMID}_${IJOB}-${RSTcount} failed\nSimulation using ${runFile} and ${oceanFile} failed before entering time loop.\nCheck your set-up and log file: ${logFile}!"
    else
      echo -e "Simulation using ${runFile} and ${oceanFile} failed before entering time loop.\nCheck your set-up and log file: ${logFile}!" | mail -s "${SIMID}_${IJOB}-${RSTcount} failed" ${USER_EMAIL}
    fi
    exit
  fi
  
  # Job status?
  # => jobOK=1: success => continue with next job
  # => jobOK=0: fail    => get time of last RST writing
  if [ ${debug} -eq 1 ]; then
    echo "jobOK = "${jobOK}
  fi
  if [ ${jobOK} -eq 1 ]; then
    # job succeeded => iterate to next job
    let IJOB=IJOB+1
    RSTcount=0
  else
    # job failed
    let BLOWUP="$(grep -i blow ${logFile} | wc -l)"
    if [ ${BLOWUP} -gt 0 ]; then
      # simulation blew up
      # send notification email to user and abort script => no further job submission
      if [ ${debug} -eq 1 ]; then
        echo -e "${SIMID}_${IJOB}-${RSTcount} failed\nSimulation using ${runFile} and ${oceanFile} blew up.\nCheck your set-up and log file: ${logFile}!"
      else
        echo -e "Simulation using ${runFile} and ${oceanFile} blew up.\nCheck your set-up and log file: ${logFile}!" | mail -s "${SIMID}_${IJOB}-${RSTcount} failed" ${USER_EMAIL}
      fi
      exit
    fi
    # check if model error occured
    let ERROR="$(tail -n200 ${logFile} | grep -i error | wc -l)"
    if [ ${ERROR} -gt 0 ]; then
      # simulation error
      # send notification email to user and abort script => no further job submission
      if [ ${debug} -eq 1 ]; then
        echo -e "${SIMID}_${IJOB}-${RSTcount} failed\nSimulation using ${runFile} and ${oceanFile} ended with error.\nCheck your set-up and log file: ${logFile}!"
      else
        echo -e "Simulation using ${runFile} and ${oceanFile} ended with error.\nCheck your set-up and log file: ${logFile}!" | mail -s "${SIMID}_${IJOB}-${RSTcount} failed" ${USER_EMAIL}
      fi
      exit
    fi
    # check if a restart was written
    let nRST="$(grep -n WRT_RST ${logFile} | wc -l)"
    if [ ${nRST} -eq 0 ]; then
      # no restart written => restart again from beginning
      NTIMES_LEFT=${NTIMES}
      # send notification email to user
      if [ ${debug} -eq 1 ]; then
        echo -e "${SIMID}_${IJOB}-${RSTcount} failed\nSimulation using ${runFile} and ${oceanFile} failed before writing restart - cluster error.\nAutomatically re-submitted."
      else
        echo -e "Simulation using ${runFile} and ${oceanFile} failed before writing restart.\nAutomatically re-submitted." | mail -s "${SIMID}_${IJOB}-${RSTcount} failed" ${USER_EMAIL}
      fi
    else # continue simulation from restart file
      # send notification email to user
      if [ ${debug} -eq 1 ]; then
        echo -e "${SIMID}_${IJOB}-${RSTcount} failed\nSimulation using ${runFile} and ${oceanFile} failed before anticipated end - cluster error.\nAutomatically re-submitted."
      else
        echo -e "Simulation using ${runFile} and ${oceanFile} failed before anticipated end - cluster error.\nAutomatically re-submitted." | mail -s "${SIMID}_${IJOB}-${RSTcount} failed" ${USER_EMAIL}
      fi
      # extract time information from last block before last restart writing
      let RSTline="$(grep -n WRT_RST ${logFile} | tail -n1 | cut -d ':' -f 1)"
      let GETline="$(head -n${RSTline} ${logFile} | grep -n GET_2DFLD | tail -n1 | cut -d ':' -f 1)"
      let nLines=RSTline-GETline
      let actualEndTime="$(tail -n+${GETline} ${logFile} | head -n${nLines} | grep : | tail -n1 | sed -e 's/^[[:space:]]*//' | cut -d' ' -f 1)"
      if [ ${debug} -eq 1 ]; then
        echo "expected end time = "${expectedEndTime}
        echo "  actual end time = "${actualEndTime}
      fi
      # left-over time steps to be caluclated
      let NTIMES_LEFT=expectedEndTime-actualEndTime
      let RSTcount=RSTcount+1
      # increase the default job time limit if very first job exceeded user-defined time limit
      if [ ${IJOB} -eq 1 ] && [ ${RSTcount} -eq 1 ] && [ ${autoLimit} -eq 1 ] && [ ${debug} -eq 0 ]; then
        echo -e "Updated job time limit\n -> old limit: ${JOB_TIME_STR}" >> ${scriptLog}
        TIME_SEC=$(bc <<< "scale=5; $NTIMES/($NTIMES-$NTIMES_LEFT)*$JOB_TIME_SEC" | cut -d '.' -f 1)
        let TIME_DAY=TIME_SEC/86400              # get number of day
        let JOB_TIME_SEC=TIME_SEC-TIME_DAY*86400 # update total seconds (minus days; time conversion with date not applicable for "days")
        JOB_TIME_STR="$(printf "%02d" $TIME_DAY)-$(date -u -d "@${JOB_TIME_SEC}" "+%H:%M:%S")"
        echo -e " -> new limit: ${JOB_TIME_STR}\n" >> ${scriptLog}
      fi
    fi
  fi

  # make copy of RST file
  if [ ${keepRST} -eq 1 ]; then
    cp -f ${rstFile} /${rstFile}_${slurmJobID}_bckp
  fi

  if [ ${debug} -eq 1 ]; then
     echo " "
  fi
  rm -f ${oceanFile} ${runFile}
  
  if [ ${jobOK} -eq 1 ] && [ ${IJOB} -gt ${NJOBS} ]; then
    rm -f ocean_${SIMID}_BASE.in run_ROMS_${SIMID}_BASE.sh
    exit
  fi
done

rm -f ocean_${SIMID}_BASE.in run_ROMS_${SIMID}_BASE.sh

exit
