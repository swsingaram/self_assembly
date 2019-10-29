#!/bin/bash
#delay is the equilibration time for the polymer
#remember you're starting on the line after the 
#equilibration time.  In this case, the starting
#line is 100 + 1
delay=101
for f in cluster*; do 
    awk '{print $1" "$2}' $f > $f.max_cluster
    st=numbound${f:7};
    paste -d " " $f.max_cluster $st > $st.max_cluster.tmp
    awk '{t=$1; $1=$3; $3=$2; $2=t; print;}' $st.max_cluster.tmp > $st.max_cluster
    tail -n +$delay $st.max_cluster > $st.max_cluster.no_delay
    rm $st.max_cluster.tmp
    rm $f.max_cluster
done
