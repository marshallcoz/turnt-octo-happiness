DESTINO="/Users/marshall/Dropbox/docTRansfer/transf"
FUENTE="/Users/marshall/Documents/DOC/coco/turnt-octo-happiness/DWLIBEM/Workdir/outs"

echo Copiar G y W a transf
echo Destino = $DESTINO

mkdir $DESTINO
for ((i=1; i<=7; i++)); do
    echo "creando carpeta" $i
    mkdir $DESTINO/traces$i
    mkdir $DESTINO/video$i

done

cd $FUENTE
for ((i=1; i<=7; i++)); do
    echo "copiando carpeta" $i
    cp traces$i/W$i.bin $DESTINO/traces$i
    cp video$i/G$i.bin $DESTINO/video$i
    cp video$i/0_MecElemvideo.mp4 $DESTINO/video$i
    cp video$i/0_video.mp4 $DESTINO/video$i
done

cp -r insbackup $DESTINO/ins
cp logfile.txt $DESTINO

echo listo

