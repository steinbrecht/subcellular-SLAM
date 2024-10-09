dir=$1
cd ${2}

#for file in "$dir"/gene_data_normalized_filtered*
for file in "$dir"/gene_data_normalized.[cmnpt][0-9SLH]*
do
   awk -F',' '{print>>$2}' ${file}
   #echo ${file} > out.txt
done

touch gene_distribution_finished.txt
