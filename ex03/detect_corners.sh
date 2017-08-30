#!/bin/bash

for i in *.ppm; do
	echo "Processing image $i ..."
	./Ex03 $i 125 5 10 2 2
done
