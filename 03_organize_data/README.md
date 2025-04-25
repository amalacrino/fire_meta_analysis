# Data download and organization

The code below downloads and organizes the data for `nf-core/ampliseq` and requires:
* `SRAtoolkit`
* `parallel`
* `cutadapt`

The code assumes the use of input files inside the folder `01_samples`. Below we report the code for 16S data, but analyses for ITS can be replicated simply replacing the input file from `02_16S.txt` to `03_ITS.txt`.

`$SRAPATH` is the path to `SRAtoolkit`.

`$PRJDIR` is the path to the project directory.

## Download data

The code below downloads the data from NCBI SRA and extract the `.fastq` files, placing data from each BioProject into a different folder.

```bash
mkdir 00_16S_data
cd 00_16S_data
tail -n +2 02_16S.txt | awk -F'\t' '{print $2}' > SRR_Acc_List.txt

listacc="SRR_Acc_List.txt"
INDIR=$PRJDIR/00_16S_data
cd $INDIR
$SRAPATH/vdb-config --prefetch-to-cwd

for i in $(cat $listacc); do
echo $i
$SRAPATH/prefetch -X 9999999999999 ${i}
done

find . -mindepth 1 -maxdepth 1 -type d | parallel -j 64 $SRAPATH/fasterq-dump --split-3 {}/{}\.sra
find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} \;
gzip *.fastq

TSV_FILE="sra.tsv"

tail -n +2 "$TSV_FILE" | awk -F'\t' '{print $1, $2}' | sort -u | while read -r BioProject Run;
do
    # Skip empty lines
    if [[ -z "$BioProject" || -z "$Run" ]]; then
        continue
    fi
    
    # Create a subdirectory for each BioProject
    mkdir -p "$BioProject"
    
    # Move the paired FASTQ files into the project directory
    for suffix in _1.fastq.gz _2.fastq.gz; do
        FILE="${Run}${suffix}"
        if [[ -f "$FILE" ]]; then
            echo "Moving $FILE to $BioProject"
            mv "$FILE" "$BioProject"
        else
            echo "Warning: $FILE not found"
        fi
    done
done
```
