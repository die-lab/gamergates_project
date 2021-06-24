#!/bin/sh

#in a directory where u have the GO.out file from pannzer.

#making the list of genes
cat GO.out | awk '{print $1}' | uniq > tmp
num=`cat tmp | wc -l`
echo "un attimo, devo processare "$num" geni..."

if [ -f almost_done ]
	then rm almost_done
fi

for i in `cat tmp`
	do a=`grep $i GO.out | awk '{print $3}'`
	echo $i"___"$a>>almost_done
	done

sed 's/ /, GO:/g' almost_done | sed 's/___/\tGO:/g' > geneID2GO
sed -i '1d' geneID2GO

sed -i 's/TRINITY_//g' geneID2GO
sed -i 's/.p[1..10]//g' geneID2GO

rm tmp almost_done
echo "Fatto!"
