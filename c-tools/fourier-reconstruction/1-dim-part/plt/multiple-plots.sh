#!/bin/sh

# Multiple run
# Usage: multiple-plots.sh PLOT_UTILITY
n0500r20="./$1 500/rho2.0 400"
n0500r33="./$1 500/rho3.3 400"
n0500r40="./$1 500/rho4.0 400"
n0500r50="./$1 500/rho5.0 400"

n1000r20="./$1 1000/rho2.0 650"
n1000r33="./$1 1000/rho3.3 450"
n1000r40="./$1 1000/rho4.0 400"
n1000r50="./$1 1000/rho5.0 670"

n1500r20="./$1 1500/rho2.0 800"
n1500r33="./$1 1500/rho3.3 400"
n1500r40="./$1 1500/rho4.0 450"

n2000r20="./$1 2000/rho2.0 600"
n2000r33="./$1 2000/rho3.3 550"

# $n0500r20
# sleep 5
# $n0500r33
# sleep 5
# $n0500r40
# sleep 5
# $n0500r50
# sleep 5
# 
# $n1000r20
# sleep 5
# $n1000r33
# sleep 5
# $n1000r40
# sleep 5
# $n1000r50
# sleep 5

$n1500r20
sleep 5
$n1500r33
sleep 5
$n1500r40
sleep 5

$n2000r20
sleep 5
# $n2000r33
# sleep 5

n1500r50="./$1 1500/rho5.0 600"
n2000r40="./$1 2000/rho4.0 600"
n2000r50="./$1 2000/rho5.0 500"
#srun $n1500r50
#srun $n2000r40
#srun $n2000r50
