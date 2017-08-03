All you should need to get started is running

tar -xzvf relicOrbit.tar.gz

to get out inert_traj_061216_reconfig.txt the Relic Orbit file I use 

Simply edit the parameters in the python script at the top of the file.

then call with casa -c relicCorrGif.py
or if in casa terminal, just
run relicCorrGif.py

Descriptions are all there.

Only 2 file inputs are sky brightness (points to png file included currently)

and orbit file (included as well), make sure format or personal orbit files matches the present file
