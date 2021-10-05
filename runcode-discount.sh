#!/usr/bin/env bash

#PBS -N fisher
#PBS -l ncpus=8
#PBS -l mem=64GB
#PBS -l walltime=01:00:00
#PBS -o fisher_stdout.out
#PBS -e fisher_stderr.out

# More info on PBS directives can be found here
# http://qcd.phys.cmu.edu/QCDcluster/pbs/run_serial.html



###############################################
#
#
#  Display PBS info
#
#
###############################################
print_pbs_info(){
    echo ------------------------------------------------------
    echo -n 'Job is running on node '; cat $PBS_NODEFILE
    echo ------------------------------------------------------
    echo PBS: qsub is running on $PBS_O_HOST
    echo PBS: originating queue is $PBS_O_QUEUE
    echo PBS: executing queue is $PBS_QUEUE
    echo PBS: working directory is $PBS_O_WORKDIR
    echo PBS: execution mode is $PBS_ENVIRONMENT
    echo PBS: job identifier is $PBS_JOBID
    echo PBS: job name is $PBS_JOBNAME
    echo PBS: node file is $PBS_NODEFILE
    echo PBS: current home directory is $PBS_O_HOME
    echo PBS: PATH = $PBS_O_PATH
    echo ------------------------------------------------------
}

###############################################
#
#
#  Helper/Setup Functions
#
#
###############################################

load_modules(){
    #activate module environment
    #NOTE: a recent HPC update means that you shouldn't need
    #to do this anymore, but I have included as a sanity check
    source /etc/profile.d/modules.sh 

    #load R
    module load r/3.6.2-foss-2019b
}    


copy_in(){
    #copy some data to  your input directory
    #nothing to copy in on this script
    #For empty bash functions, must put a colon in them,
    #otherwise it will throw an error
    :
}


copy_out(){
    #nothing to copy out on this script
    #For empty bash functions, must put a colon in them,
    #otherwise it will throw an error
    :
}



run_program(){
    #make sure we change to the current directory
    #where this bash job script is
    cd $PBS_O_WORKDIR
    #now run the R script
    #echo 'R_LIBS_USER="~/R/library_3.6.2"' >  ~/.Renviron
    Rscript fisher_for_github.R --discount 1   
}


run_clean(){
    #nothing to clean for this script
    #For empty bash functions, must put a colon in them,
    #otherwise it will throw an error
    :
}

###############################################
#
#
#  Running everything
#
#
###############################################

print_pbs_info
load_modules
copy_in
copy_out
run_program
run_clean
