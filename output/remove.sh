#!/bin/bash

REDVOL=( 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 )
CONFIN=( 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 )


#REDVOL=( 75 )
#CONFIN=( 90 )

# 1. create all folders
# 2. copy solution at high Ca from solving the 
#    full tenth-order system of ODEs
for i in "${REDVOL[@]}"
do
	cd v$i
	for j in "${CONFIN[@]}"
	do
		cd conf$j
			rm sln_*
		cd ..
	done
	cd ..
done
