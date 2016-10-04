#!/bin/sh

../combine_runs_align.py 500/rho2.0 0
../combine_runs_align.py 500/rho3.3 0
../combine_runs_align.py 500/rho4.0 0
../combine_runs_align.py 500/rho5.0 0

../combine_runs_align.py 1000/rho2.0 0
../combine_runs_align.py 1000/rho3.3 0
../combine_runs_align.py 1000/rho4.0 0
../combine_runs_align.py 1000/rho5.0 0

../combine_runs_align.py 1500/rho2.0 0
../combine_runs_align.py 1500/rho3.3 0
../combine_runs_align.py 1500/rho4.0 0

../combine_runs_align.py 2000/rho2.0 0
../combine_runs_align.py 2000/rho3.3 0

../combine_runs_shape.py 500/rho2.0 0
../combine_runs_shape.py 500/rho3.3 0
../combine_runs_shape.py 500/rho4.0 0
../combine_runs_shape.py 500/rho5.0 0

../combine_runs_shape.py 1000/rho2.0 0
../combine_runs_shape.py 1000/rho3.3 0
../combine_runs_shape.py 1000/rho4.0 0
../combine_runs_shape.py 1000/rho5.0 0

../combine_runs_shape.py 1500/rho2.0 0
../combine_runs_shape.py 1500/rho3.3 0
../combine_runs_shape.py 1500/rho4.0 0

../combine_runs_shape.py 2000/rho2.0 0
../combine_runs_shape.py 2000/rho3.3 0
