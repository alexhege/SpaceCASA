#Alex Hegedus
#This script changes the python and config files to automate generation of


#!/bin/bash

#choose multiple of 10 so multithreading works
for startPos in 0 32 64
do
  if [ ! -d "startPos_${startPos}" ]; then
    mkdir startPos_${startPos}
  fi 
  cd startPos_${startPos}
  for numsc in 5 6 4
  do
    if [ ! -d "numSC_${numsc}" ]; then
      mkdir numSC_${numsc}
    fi
    cd numSC_${numsc}

    max=100
    #this updates the how far into time the software should go
    mult=$((max/10))
      #bash ../killer.sh &

    for j in `seq 0 $(($mult-1))`
    do

      for i in `seq 1 10`
      do
        z=$((10*${j}+${i}))

        mkdir trial_${z}

        cp ../../*py ../../*dat ../../*sh ../../*csv trial_${z}

        cd trial_${z}
        #nohup bash runcorr.sh ${startPos} ${numsc} 
        nohup casa --nologger --nologfile --nogui -c corrCME.py -posFile SunRISE_relTrajs_sixSC.csv -componentFile components_new.dat -startPos ${startPos} -correlate True -numSC ${numsc} | tee out.out
        sleep 1.5
        #fg
        cd ..
      done

      #wait ${!}
      echo "now done with 10x$j \n"
      
        #avoids deleting important stuff before code ends
      sleep 5
      for i in `seq 1 10`
      do
        z=$((10*${j}+${i}))

        cd trial_${z}
        bash rmall.sh

        cd ..
      done

    done

    cd ..
  done
  cd ..

done
