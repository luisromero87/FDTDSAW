#!/bin/bash
for i in `seq 1 10`;
do
	mpirun -n 8 fdtdsaw
done
