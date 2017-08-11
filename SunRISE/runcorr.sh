casa --nologger --nologfile --nogui -c corrCME.py -posFile SunRISE_relTrajs_sixSC.csv -outDir . -componentFile components_new.dat -startPos $1 -correlate True -numSC $2 | tee out.out 
