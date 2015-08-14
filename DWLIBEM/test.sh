# asorted tests
echo "*****************"
echo "* Asorted tests *"
echo "*****************"
TESTSDIR=tests 

CONTADOR=0
for entry in "$TESTSDIR"/*; do
	let CONTADOR+=1
	EJE[$CONTADOR]=$entry
	echo "${EJE[$CONTADOR]}"
done
echo "There are " $CONTADOR " cases."


let NTESTS=$CONTADOR 
# Test nomenclature follows
# [pol][pw_c][frec][nk][Layer][Layer]...[Hs][geom]
# 
rm -rf ./Workdir/outs
rm -rf ./Workdir/$TESTSDIR
mkdir ./Workdir/$TESTSDIR
for (( i = 1; i <= NTESTS; i++ )); do
	echo "################################"
	echo "           test " $i "          "
    # copiar carpeta de ejercicio
    rm -rf ./Workdir/ins
    cp -rv ./${EJE[$i]} ./Workdir/ins
    
    CONTADOR=0
    while read line; do
    	let CONTADOR+=1
    	echo  $CONTADOR ">>" $line
    done < ./Workdir/ins/info.txt
    
    # correr comando en la última línea de info.txt
    echo $line
    eval $line

    # guardar resultados
    cp -r ./Workdir/outs ./Workdir/${EJE[$i]}
done
 echo "FINISH TESTING"
