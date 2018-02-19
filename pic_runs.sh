for ((T=280; T<=500;T+=10)); do
 	echo $T
	cp ./in_copy.txt ./in.txt
	sed -i -e "s/XXX/$T/g" ./in.txt
	./mclattice2
	../python/read.py $T
done
