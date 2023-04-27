#!/bin/bash

#This script was written to request resources to analyze 394 subjects' data
#from the Autism Biomarkers Consortium for Clinical Trials dataset. 
#Adam Dede, adam.osman.dede@gmail.com, 2022



# the following lines specify computing resource needs
# the format of these lines will differ on different HPC systems
# check with your local documentation and research support staff

#$ -l rmem=16G              #how much memory will each job need? 
#$ -pe smp 1                #how many compute nodes will each job need? 
#$ -l h_rt=24:00:00         #how long will each job last? 
#$ -t 1-394                 #an iterator variable defining the job IDs




#Output the Task ID
#SGE_TASK_ID is the variable name of the iterator variable defined in 
#line 10 above. NOTE: this variable name may be different on different HPC systems.
echo "Task ID is $SGE_TASK_ID"
#this output will be printed to a job-specific text file. 
#this output file will also be the target of anything that 
#Matlab would normally print to the console during interactive,
#local use. 

#specify the analysis software that should be loaded into the workspace
module load apps/matlab/2021b/binary

#specify the working directory for the system to be in before starting
#analysis
cd /shared/dede_group/User/pc1aod/CODE

#start Matlab in batch mode
#batch mode means that Matlab will not display its graphical user interface
#instead, it will take in command line input in the shell
#everything inside the double quotes is input to Matlab as if it were being
#typed directly into the Matlab command line
#the SGE_TASK_ID iterator variable value is passed into Matlab with the variable
#name jobID. Next, the pipeline wrapper script (getBioConsortDat) is called.
matlab -batch "jobID=$SGE_TASK_ID; getBioConsortDat"



#General note: The logic of this script is similar to a for loop. The major difference
#is that where a for loop is run on a single computer processor, this script is 
#specifying a series of jobs to be initiated on different processors. This allows 
#the jobs to run in parallel. 

