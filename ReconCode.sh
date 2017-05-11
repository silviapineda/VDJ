#!/bin/bash
FILES=/Users/Pinedasans/programs/Recon-master2/gDNA/

for f in $FILES
do  
    var=$(basename $f)
    python2.7 recon_v2.2.py -R -t 30 -o gDNA/$var_fitfile.txt gDNA/$FILE
    python2.7 recon_v2.2.py -x --x_max 30 -o gDNA/$var_plotfile.txt -b error_bar_parameters.txt gDNA/$var_fitfile.txt
done

