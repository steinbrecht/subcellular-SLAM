# download the "Genome sequence, primary assembly (GRCh38)" fasta file
wget -v ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/GRCm38.primary_assembly.genome.fa.gz -P "$1"
# filter it as described in the note below
# download the annotations that correspond to it 
wget -v ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M14/gencode.vM14.basic.annotation.gff3.gz -P "$1"
wget -v ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M14/gencode.vM14.basic.annotation.gtf.gz -P "$1"

gunzip "$1"/GRCm38.primary_assembly.genome.fa.gz
gunzip "$1"/gencode.vM14.basic.annotation.gff3.gz
gunzip "$1"/gencode.vM14.basic.annotation.gtf.gz

grep "protein_coding" "$1"/gencode.vM14.basic.annotation.gtf  > "$1"/gencode_protein_coding.gtf
grep -v "mt-" "$1"/gencode_protein_coding.gtf  > "$1"/gencode_protein_coding_no_mt.gtf
