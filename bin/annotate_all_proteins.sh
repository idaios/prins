for i in `ls *.pdb`;
do
    echo $i;
    protein3Dannotation.pl score=score_$i.txt $1 pdb=$i > ${i}_ANNOT.pdb  2>${i}_ERR
done
