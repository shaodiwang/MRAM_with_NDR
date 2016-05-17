#!/bin/csh -f
#  run_ndr.sh.cmd
#
#  UGE job for run_ndr.sh built Thu May 12 00:29:44 PDT 2016
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.joblog.$JOB_ID
#$ -o /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=2048M,h_rt=12:00:00
#
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M shaodiwa@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR
  set qqjob     = run_ndr.sh
  set qqodir    = /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR
  cd     /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for run_ndr.sh built Thu May 12 00:29:44 PDT 2016"
  echo ""
  echo "  run_ndr.sh directory:"
  echo "    "/u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "run_ndr.sh started on:   "` hostname -s `
  echo "run_ndr.sh started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load intel/13.cs
  module load cuda/5.0
#
  echo run_ndr.sh "" \>\& run_ndr.sh.output.$JOB_ID
  echo ""
  time /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh  >& /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.output.$JOB_ID
#
  echo ""
  echo "run_ndr.sh finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/project/puneet/shaodiwa/MeRAM/NDR/MRAM_with_NDR/run_ndr.sh.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
