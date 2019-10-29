#!/bin/bash

for f in */*.dcd
do
l=${#f}
let l-=4 
var=${f:0:l}
sbatch rg.sbatch $var
done

