
NPRODS=150
NCPUS=4
SYSTEM=$1

prod=1
for ((i=1;i<=10;i++))
do
	j=1
	while ((j <= $NCPUS)) && ((prod <= $NPRODS))
	do
		echo $j $i $prod
		((b=$prod+4))
		printf -v x "%03d" $prod
		printf -v y "%03d" $b
		mkdir $x.$y.RMSF
		cd $x.$y.RMSF
		sed -e s/AAA/$prod/g -e s/BBB/$b/g -e s/NNN/$SYSTEM/g -e s/XXX/$x/g -e s/YYY/$y/g < ../rmsf.config > $x.$y.rmsf.config
		time ../rmsf_calc.py $x.$y.rmsf.config > RMSF.output &
		cd ../
		((j=$j+1))
		((prod=$prod+5))
	done
	wait
done

mkdir 021.150.RMSF
cd 021.150.RMSF
sed -e s/AAA/21/g -e s/BBB/150/g -e s/NNN/$SYSTEM/g -e s/XXX/021/g -e s/YYY/150/g < ../rmsf.config > 021.150.rmsf.config
time ../rmsf_calc.py 021.150.rmsf.config > RMSF.output 

