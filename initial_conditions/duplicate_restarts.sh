##set a restart date of 1980 for climate runs
#for i in {0..9}
#do
#	cd exp_01$i 
#	ln -s 1983010100.rst 1980010100.rst
#	ln -s 1983010100_initial.rst 1980010100_initial.rst
#	cd ../exp_06$i
#	ln -s 1983010100.rst 1980010100.rst
#	ln -s 1983010100_initial.rst 1980010100_initial.rst
#	cd ..	
#done

#for restarting the stochastic rounding jobs
#after they were mysteriously stopped on atmlxint6.
#different dates needed for different runs, restart
#files have been moved into folders labelled exp_80*
#in advance.
cd exp_800 
ln -s 1983010100.rst 2004010100.rst
cd ../exp_801
ln -s 1983010100.rst 2004010100.rst
cd ../exp_802
ln -s 1983010100.rst 2003010100.rst
cd ../exp_803
ln -s 1983010100.rst 2004010100.rst
cd ../exp_804
ln -s 1983010100.rst 2006010100.rst
cd ..
