for dir in 
do 
cd
cat merged.fastq | sed -n '2~4p' | cut -c1-8 | sort -r | uniq -c | sort -nrk1,1 > UMI_counts.txt
cd
done
