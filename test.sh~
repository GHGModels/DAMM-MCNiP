#! /bin/sh -l
#
# Submit job with the command: qsub filename
#
# To view the status of the job: qstat -u username
#
# Set the runtime limit (default 12 hours): 
# -l h_rt=24:00:00
#
# Send email when the job is done (default: no email is sent)
# -m e
#
# Give the job a name (default: script name)
#$ -N dmc
#Add dollar sign below to run as command if need some minimum memory
#$ -l mem_total=126G
#
# -o file.txt
#
## end of qsub options
matlab -nodisplay -singleCompThread -r MCNiP_TMloop.m


