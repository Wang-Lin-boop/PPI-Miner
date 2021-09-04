#!/bin/bash




function post_profit(){
    fasta=$1
    title=`echo $2 | awk '{print $1}'`
    start=`echo $2 | awk '{print $2}'`
    ((start-=1))
    profit_txt=$2
    msa_len=$3
    seq=`grep -A 1 "$title" $fasta | tail -n 1 `
    echo "${profit_txt}  ${seq:${start}:${msa_len}}"
};export -f post_profit

export CPUnum=40
export msa="msa-order.fasta"
export msa_len=8
((msa_len+=1))
export database="Human_Disorder_Sequences.fasta"  # don't contains any ":" or " " in fasta title!
export percentage=75

( echo "${msa}"; echo "F"; echo "pfm-matrix"; echo "${percentage}"; echo "prophecy.out" ) | prophecy
( echo "prophecy.out"; echo "${database}"; echo "profit.out" ) | profit

tail -n +4 profit.out | parallel -j $CPUnum post_profit ${database} {} ${msa_len}

