#!/bin/bash
#$ -q mh.q
#$ -cwd
#$ -j n
module load PYTHON/ANACONDA-2.5.0
qcharge=6.0
pollength=800.0
ebond=6.5
if [ ! -d "ANALYSIS" ]; then
    mkdir "ANALYSIS"
fi

num=1
while [ $num -le 5 ]; do
file="QCHARGE-$qcharge""_POLLENGTH-$pollength""_EBOND-$ebond""_$num"
echo $file
cd ./$file
catdcd -o temp.mol2 -otype mol2 -stype mol2 -stride 1 -s *.mol2 *.dcd
python ../CapsidAnalysisT3_topbeads.py $qcharge $pollength $ebond $num
rm temp.mol2
let "num+=1"
mv cluster*.dat ../ANALYSIS/
mv numbound*.dat ../ANALYSIS/
cd ..
done
