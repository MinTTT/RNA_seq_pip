#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CPU=8
DEDUO=false
RUN_PIPELINE=true
CLEAN_OUTPUT_DIR="./cleanData"
CLEAN_SCRIPT="$SCRIPT_DIR/miniBac-seq_raw_data_process_clean_demix_20250401.py"
STATISTICS_SCRIPT="$SCRIPT_DIR/miniBac-seq_pip_20250105.py"
RET_OUTPUT_DIR="./"

usage() {
    echo "Usage: $0 -i <raw_data_dir> -g <gff_file> -f <fasta_file> [options]"
    echo ""
    echo "Required arguments:"
    echo "  -i, --input-dir      Directory containing raw *.fastq.gz files"
    echo "  -g, --gff-file       Path to GFF annotation file"
    echo "  -f, --fasta-file     Path to reference genome FASTA file"
    echo ""
    echo "Optional arguments:"
    echo "  -c, --cpu            Number of CPUs to use (default: 8)"
    echo "  -d, --dedup          Enable deduplication during fastp processing"
    echo "  -r, --run-pipeline   Run pipeline without prompts"
    echo "  -t, --threading-max  Max threading number for alignment (default: 8)"
    echo "  -o, --clean-output   Output directory for cleaned data (default: ./cleanData)"
    echo "  -p, --ret-output     Output directory for results (default: ./)"
    echo "  -h, --help           Display this help message"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -g|--gff-file)
            GFF_FILE="$2"
            shift 2
            ;;
        -f|--fasta-file)
            FASTA_FILE="$2"
            shift 2
            ;;
        -c|--cpu)
            CPU="$2"
            shift 2
            ;;
        -d|--dedup)
            DEDUP=true
            shift
            ;;
        -r|--run-pipeline)
            RUN_PIPELINE=true
            shift
            ;;
        -t|--threading-max)
            THREADING_MAX="$2"
            shift 2
            ;;
        -o|--clean-output)
            CLEAN_OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--ret-output)
            RET_OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Function to print section headers
print_section() {
    echo ""
    echo "=========================================="
    echo "$1"
    echo "=========================================="
    echo ""
}

# Step 1: clean the data
print_section "Step 1: Cleaning the raw data"
echo "Running data cleaning script..."
CLEAN_CMD="python \"$CLEAN_SCRIPT\" -c $CPU -i $INPUT_DIR -o $CLEAN_OUTPUT_DIR"

if [[ "$RUN_PIPELINE" == "true" ]]; then
  CLEAN_CMD="$CLEAN_CMD -r"
fi
echo "Executing: $CLEAN_CMD"
bash -c "$CLEAN_CMD"
# Check if cleaning was successful
if [[ $? -ne 0 ]]; then
    echo "Error: Data cleaning failed!"
    exit 1
fi
# step 2: run the pipeline
print_section "Step 2: Running the analysis pipeline"
echo "Running pipeline script..."
ALIGN_CMD="python \"$STATISTICS_SCRIPT\" -gff \"$GFF_FILE\" -fa \"$FASTA_FILE\" -o \"$RET_OUTPUT_DIR\" -fq \"$CLEAN_OUTPUT_DIR\""
if [[ "$RUN_PIPELINE" == "true" ]]; then
    ALIGN_CMD="$ALIGN_CMD -r"
fi
echo "Executing: $ALIGN_CMD"
bash -c "$ALIGN_CMD"
exit 0