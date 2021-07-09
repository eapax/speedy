for i in {0..9}
do
	cd exp_01$i 
	ln -s 1983010100.rst 1980010100.rst
	ln -s 1983010100_initial.rst 1980010100_initial.rst
	cd ../exp_06$i
	ln -s 1983010100.rst 1980010100.rst
	ln -s 1983010100_initial.rst 1980010100_initial.rst
	cd ..	
done
