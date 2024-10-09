sed -e '2,${/^num_Ts/d' -e '}' | sort -k 1,2 | awk -f scripts/compile_frequencies_both_intron_and_exon.awk | tr " " "\t" | grep -v "0\s0\s0\s0"
