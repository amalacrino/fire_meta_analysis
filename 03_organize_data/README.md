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

## Cutadapt

The code below removes the specific primers used to generate each individual library.

`$PRJDIR` is the path to the project directory.

`$LISTACC` is either `02_16S.txt` or `03_ITS.txt`.

`$DATA_DIR` is the input directory (where the data has been downloaded). For the example above it is `$PRJDIR/00_16S_data`

`$OUTPUT_DIR` is the directory where you'd like the cleaned data to be stored.

```bash

# Read the file and process each line
tail -n +2 "$LISTACC" | awk -F'\t' '{print $1, $2, $5, $6}' | sort -u | while read -r BioProject Run FW_seq REV_seq;
do
    # Skip empty lines or missing data
    if [[ -z "$BioProject" || -z "$Run" || -z "$FW_seq" || -z "$REV_seq" ]]; then
        continue
    fi

    # Define project and sample directories
    PROJECT_DIR="$OUTPUT_DIR/$BioProject"
    mkdir -p "$PROJECT_DIR"
    
    # Define input and output files
    INPUT_FWD="$DATA_DIR/$BioProject/${Run}_1.fastq.gz"
    INPUT_REV="$DATA_DIR/$BioProject/${Run}_2.fastq.gz"
    OUTPUT_FWD="$PROJECT_DIR/${Run}_1.fastq"
    OUTPUT_REV="$PROJECT_DIR/${Run}_2.fastq"
    
    # Check if input files exist before processing
    if [[ -f "$INPUT_FWD" && -f "$INPUT_REV" ]]; then
        echo "Processing $Run in $PROJECT_DIR"
        
        # Run cutadapt in parallel for forward and reverse reads
        cutadapt -j 64 -g "$FW_seq" -G "$REV_seq" \
            -o "$OUTPUT_FWD" -p "$OUTPUT_REV" "$INPUT_FWD" "$INPUT_REV"
    else
        echo "Warning: Paired files for $Run not found in $DATA_DIR/$BioProject"
    fi
done

echo "Primer trimming complete!"

cd $OUTPUT_DIR
gzip -r ./

```
