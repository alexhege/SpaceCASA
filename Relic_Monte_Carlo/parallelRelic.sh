#Alex Hegedus 8/10/17
#This is used with genMC.py, runit.sh, and corrRelicMC.py to run Monte Carlo Simulations on RELIC 
#currently runs 10 at a time

#!/bin/bash
#make sure casa is defined


#casa --nologger --nologfile -c genMC.py &
#wait ${!}

begin=`date +%s`

for numsc in 16 24 32
do


  for relstr in 1.1 1.3 1.5 2.0 4.0
  do

    max=100 #12
    
    start=`date +%s`

    for j in `seq 0 9`
    do
      start2loop=`date +%s`

      #clean possible old files
      for i in `seq 0 9`
      do
        z=$((10*${j}+${i}))
        cd MCImages/relStr_${relstr}/run_${z}/
        rm -rf numSC_${numsc}
        rm *log
        cd ../../..
      done 

      for i in `seq 0 9`
      do
        z=$((10*${j}+${i}))
        #cd MCImages/relStr_${relstr}/run_${z}/
        #casa --nologger --nologfile -c corrRelicMC.py ${relstr} ${z} ${numsc} 
	nohup bash runit.sh ${relstr} ${z} ${numsc} | tee MCImages/relStr_${relstr}/run_${z}/runout${numsc}.out &	
        #cd ../../../

      done

      wait ${!}
      end=`date +%s`
      runtime=$((end-start2loop))
      echo "now done with 10x$j \n"
      echo "last 10 took ${runtime} seconds \n"

      #remove unneccessary data if low on space after each set
      for i in `seq 0 9`
      do
        z=$((10*${j}+${i}))
        if [ -d "MCImages/relStr_${relstr}/run_${z}/numSC_${numsc}" ]; then
          cd MCImages/relStr_${relstr}/run_${z}/numSC_${numsc}
          rm -rf *ms *dirty*
          cd ../../../../
        fi
      done

    done

  done

done

truend=`date +%s`

allruntime=$((truend-begin))
echo "total time was ${allruntime} seconds \n"


