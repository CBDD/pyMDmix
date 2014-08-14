cd min
sander -O -i min.in -o min.out -p ../pep.prmtop -c ../pep.prmcrd -r min.rst
cd ../eq
pmemd.cuda -O -i eq1.in -o eq1.out -p ../pep.prmtop -c ../min/min.rst -r eq1.rst -x eq1.nc -ref ../min/min.rst
pmemd.cuda -O -i eq2.in -o eq2.out -p ../pep.prmtop -c eq1.rst -r eq2.rst -x eq2.nc -ref ../min/min.rst
pmemd.cuda -O -i eq3.in -o eq3.out -p ../pep.prmtop -c eq2.rst -r eq3.rst -x eq3.nc -ref ../min/min.rst
pmemd.cuda -O -i eq4.in -o eq4.out -p ../pep.prmtop -c eq3.rst -r eq4.rst -x eq4.nc -ref ../min/min.rst
pmemd.cuda -O -i eq5.in -o eq5.out -p ../pep.prmtop -c eq4.rst -r eq5.rst -x eq5.nc -ref ../min/min.rst
cd ../md
pmemd.cuda -O -i md.in -o md1.out -p ../pep.prmtop -c ../eq/eq5.rst -r md1.rst -x md1.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md2.out -p ../pep.prmtop -c md1.rst -r md2.rst -x md2.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md3.out -p ../pep.prmtop -c md2.rst -r md3.rst -x md3.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md4.out -p ../pep.prmtop -c md3.rst -r md4.rst -x md4.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md5.out -p ../pep.prmtop -c md4.rst -r md5.rst -x md5.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md6.out -p ../pep.prmtop -c md5.rst -r md6.rst -x md6.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md7.out -p ../pep.prmtop -c md6.rst -r md7.rst -x md7.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md8.out -p ../pep.prmtop -c md7.rst -r md8.rst -x md8.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md9.out -p ../pep.prmtop -c md8.rst -r md9.rst -x md9.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md10.out -p ../pep.prmtop -c md9.rst -r md10.rst -x md10.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md11.out -p ../pep.prmtop -c md10.rst -r md11.rst -x md11.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md12.out -p ../pep.prmtop -c md11.rst -r md12.rst -x md12.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md13.out -p ../pep.prmtop -c md12.rst -r md13.rst -x md13.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md14.out -p ../pep.prmtop -c md13.rst -r md14.rst -x md14.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md15.out -p ../pep.prmtop -c md14.rst -r md15.rst -x md15.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md16.out -p ../pep.prmtop -c md15.rst -r md16.rst -x md16.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md17.out -p ../pep.prmtop -c md16.rst -r md17.rst -x md17.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md18.out -p ../pep.prmtop -c md17.rst -r md18.rst -x md18.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md19.out -p ../pep.prmtop -c md18.rst -r md19.rst -x md19.nc -ref ../min/min.rst
pmemd.cuda -O -i md.in -o md20.out -p ../pep.prmtop -c md19.rst -r md20.rst -x md20.nc -ref ../min/min.rst