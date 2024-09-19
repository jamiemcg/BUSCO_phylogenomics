mkdir -p $PREFIX/bin

cp $SRC_DIR/BUSCO_phylogenomics.py $PREFIX/bin/
cp $SRC_DIR/count_buscos.py $PREFIX/bin/

chmod +x $PREFIX/bin/BUSCO_phylogenomics.py
chmod +x $PREFIX/bin/count_buscos.py
