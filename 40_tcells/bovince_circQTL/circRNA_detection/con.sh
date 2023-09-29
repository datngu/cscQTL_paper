#!/bin/bash

# circall

rm */*.id
# rm circall/EGAR*bed
# rm ciri2/EGAR*bed
# rm circexp2/EGAR*bed
# rm circRNA_finder/EGAR*bed.bed

for fi in circall/EGAR*
do
    awk 'NR>1 && $6 > 1 {printf "%s\t%s\t%s\t%s__%s__%s\t.\t%s\n", $1, $2, $3, $1, $2, $3, $6}' $fi > ${fi}.id
    echo "done ${fi}.id"
done


# ciri2

for fi in ciri2/EGAR*
do
    awk 'NR>1 && $5 > 1 {printf "%s\t%s\t%s\t%s__%s__%s\t.\t%s\n", $2, $3, $4, $2, $3, $4, $5}' $fi > ${fi}.id
    echo "done ${fi}.id"
done

# circexp2

for fi in circexp2/EGAR*
do
    awk '$13 > 1 {printf "%s\t%s\t%s\t%s__%s__%s\t.\t%s\n", $1, $2+1, $3, $1, $2+1, $3, $13}' $fi > ${fi}.id
    echo "done ${fi}.id"
done


# circRNA_finder

for fi in circRNA_finder/EGAR*bed
do
    awk '$5 > 1 {printf "%s\t%s\t%s\t%s__%s__%s\t.\t%s\n", $1, $2+1, $3, $1, $2+1, $3, $5}' $fi > ${fi}.id
    echo "done ${fi}.id"
done


##

wc -l circall/*.id | sed '$d' > circall.stat
wc -l ciri2/*.id | sed '$d' > ciri2.stat
wc -l circexp2/*.id | sed '$d' > circexp2.stat
wc -l circRNA_finder/*.id | sed '$d' > circRNA_finder.stat