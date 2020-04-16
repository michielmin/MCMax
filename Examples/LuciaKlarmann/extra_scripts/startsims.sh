#!/bin/bash

# this is my old bashscript to run large numbers of simulations.
# this program needs to be provided with: Name and path to .run and .ray file, name of log file that will be produced, number of RT models to run in parallel, estimate of time it takes to run all Rt models
# !!! marks where the input needs to be changed
# once the input is updated, make executable and start the entire model grid with ./startsims.sh

log=give_the_logfile_a_name.dat #!!!
torunpara_Sims=7 #Number of Sims to run parallel !!! should be #cores -1 

echo $torunpara_Sims sims to run parallel >> out #for debugging
echo $torunpara_Sims sims to run parallel >> $log #for debugging

first=1 #to check later if at least one sim has finished

old_IFS=$IFS
IFS=$'\n'	# use newline as internal field separator
all_Sims=($(cat /media/lucia3/Data/Paola/Paola_data_cor/Paola_data_cor.run)) # read names of Sims into array all_Sims, !!! provide .run file
IFS=$old_IFS  # restore internal field separator

Number_of_all_Sims=${#all_Sims[@]}

echo ${#all_Sims[@]} number of all sims >> out #for debugging
echo ${#all_Sims[@]} number of all sims >> $log #for debugging
#echo $Number_of_all_Sims #for debugging

# start initial group of Sims (size of group from torunpara_Sims)

torunpara_Sims_less1=$(expr $torunpara_Sims - 1) #correct because indices start with 0, not 1
echo to run sims minus 1 = $torunpara_Sims_less1 >> out

for i in $(seq 0 $torunpara_Sims_less1) # start simulations until torunpara_Sims are started
do
	${all_Sims[$i]} & #start Simulation i
	pid_Sims[$i]=$! #get pid of started simulation and store in array
	echo ${pid_Sims[$i]} pid $i >> out # for debugging
	echo ${pid_Sims[$i]} pid $i >> $log # for debugging
	
	if [ $first == 1 ]; then
	    first=2
	    first_pid=${pid_Sims[$i]}
	    echo first_pid $first_pid >> out
	    echo first_pid $first_pid >> $log
	fi	
done


Sims_started=$torunpara_Sims #the first torunpara_Sims are started
echo $Sims_started sims started >> out # for debugging
echo $Sims_started Sims started >> $log

if [ $Sims_started -lt $Number_of_all_Sims ]; then
        echo $Sims_started gestartet ist  kleiner $Number_of_all_Sims alle >> out
	echo $Sims_started gestartet ist  kleiner $Number_of_all_Sims alle >> $log
fi

echo vor erster w schleife >> out

#while [ $Sims_started -lt $Number_of_all_Sims ]; do
#
#        echo w loop should start now  >> out
#	break
#done
	

while [ $Sims_started -lt $Number_of_all_Sims ]; do # run while there are simulations to run,
	
        echo inside w loop >> out
        sleep 20 #give Sims some time to run
	echo giving sims some time to run >> out
	echo sequenz $(seq 0 $torunpara_Sims_less1) >> out

	
	for j in $(seq 0 $torunpara_Sims_less1); do #check for all currently running Sims
	        echo  j loop running >> out
		echo $j >> out
		echo run >> out
		echo ${pid_Sims[$j]} this is pid >> out

#		if [ ! kill ${pid_Sims[$j]} > /dev/null 2>&1 ]; then # is the Simulation still running
		#echo 
		if  ps -p ${pid_Sims[$j]}  > /dev/null
		then
															#i.e. can its pid be found?
			sleep 1 #do nothing and wait 1 s
			echo do nothing but wait a second >> out
			echo killcheck  >> out
		
		else #if Simulation finished, start next simulation from all_Sims list

			${all_Sims[$Sims_started]} & 
			pid_Sims[$j]=$! #assign new pid			
			echo new pid ${pid_Sims[$j]}  >> out
			echo new pid ${pid_Sims[$j]}  >> $log
			echo Simulation Number $Sims_started >> out #(that is ok, as index starts with 0, not 1)
			echo Simulation Number $Sims_started >> $log #(that is ok, as index starts with 0, not 1)

		
			Sims_started=$(expr $Sims_started + 1) #one more Sim started
			echo now $Sims_started sims are started >> out
			echo now $Sims_started sims are started >> $log
		fi
		
		if  [ $Sims_started == $Number_of_all_Sims ]
		then
                  echo broke out of loop, all runs started >> out
                  echo broke out of loop, all runs started >> $log
		    break
		fi

	done
	echo still running >> out #for debugging
done

echo after first w loop >> out  #for debugging
echo after first w loop >> $log  #for debugging

echo wait 1 hour before raytracing >> out
echo wait 1 hour raytracing >> $log


sleep 3600 #!!! estimate how many seconds it will take to finish *all* RT calculation

#do not just wait, check if first run is done before starting
while ps -p $first_pid > /dev/null; do
    echo waiting for run to finish before raytracing >> out
    echo waiting for run to finish before raytracing >> $log
    sleep 30
done


#Start the raytracing. Only one process is running, as it can use all cores

#read in commands
old_IFS=$IFS
IFS=$'\n'	# use newline as internal field separator
all_Sims_ray=($(cat /media/lucia3/Data/Paola/Paola_data_cor/Paola_data_cor.ray)) # read names of Sims into array all_Sims !!!
IFS=$old_IFS  # restore internal field separator

echo $all_Sims_ray >> out


raytrace_started=0
echo raytrace_started $raytrace_started >> out
echo raytrace_started $raytrace_started >> $log
echo about to start raytraching ${#all_Sims_ray[@]} sims >> out
echo about to start raytraching ${#all_Sims_ray[@]} sims >> $log

while [ $raytrace_started -lt ${#all_Sims_ray[@]} ]; do # start all simulations in list, one after the other is finished
    echo raytracing >> out
    echo raytracing >> $log
    
    ${all_Sims_ray[$raytrace_started]} & #start Simulation 

    pid_Sim=$! #get pid of started simulation
    echo $pid_Sim pid >> out # for debugging
    echo $pid_Sim pid >> $log # for debugging
    echo '$all_Sims_ray[$raytrace_started]' >> out
    echo ${all_Sims_ray[$raytrace_started]} >> out
    echo ${all_Sims_ray[$raytrace_started]} >> $log

    raytrace_started=$(expr $raytrace_started + 1)
    echo $raytrace_started started raytraces i total >> out
    echo $raytrace_started started raytraces i total >> $log

    echo pid sim $pid_Sim >> out
    if ps -p $pid_Sim > /dev/null
    then
	echo waiting loop should start now >> out
	echo waiting loop should start now >> $log
    else
	echo at least inside if statment >> out
    fi

    while ps -p $pid_Sim > /dev/null
    do #check if raytracing is still going on
	echo inside waiting loop >> out
	    
	sleep 60 #wait for one minute before checking again
	echo slept 60 seconds >> out

    done #if raytracing is finished

    if [ $raytrace_started -eq ${#all_Sims_ray[@]} ]
    then
	echo broke out of raytrace loop >> out
	echo broke out of raytrace loop >> $log
	break
    else
	echo raytrace_started, all raytraces: $raytrace_started, ${#all_Sims_ray[@]} >> out
    fi
    echo raytraces started: $raytrace_started >> out
    echo raytraces started: $raytrace_started >> $log
done
