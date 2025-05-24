#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]; then
    echo "Usage: $0 <adaptor_length> <input_fasta_file> <output_fasta_file>"
    exit 1
fi

ADAPTOR_LEN=$1
INPUT_FILE=$2
OUTPUT_FILE=$3

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file does not exist - $INPUT_FILE"
    exit 1
fi

# Ensure output directory exists
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Process FASTA file using seqtk and awk
seqtk seq "$INPUT_FILE" | \
awk -v adaptor_len=$ADAPTOR_LEN '
    BEGIN {
        seq = ""
        actual_len = adaptor_len + 1  #  + 1
    }
    /^>/ {
        if (seq != "") {
            # Process previous sequence
            len = length(seq)
            # Calculate number of repeats needed
            repeat_times = int(actual_len / len) + (actual_len % len > 0)
            # Generate repeated sequence
            repeated_seq = ""
            for (i = 1; i <= repeat_times; i++) {
                repeated_seq = repeated_seq seq
            }
            # Extract adaptor from end of repeated sequence
            adaptor = substr(repeated_seq, length(repeated_seq) - actual_len + 1)
            # Output modified sequence
            print prev_header
            print adaptor seq
        }
        prev_header = $0
        seq = ""
    }
    !/^>/ {
        seq = seq $0
    }
    END {
        # Process last sequence
        if (seq != "") {
            len = length(seq)
            repeat_times = int(actual_len / len) + (actual_len % len > 0)
            repeated_seq = ""
            for (i = 1; i <= repeat_times; i++) {
                repeated_seq = repeated_seq seq
            }
            adaptor = substr(repeated_seq, length(repeated_seq) - actual_len + 1)
            print prev_header
            print adaptor seq
        }
    }
' > "$OUTPUT_FILE"

echo "Processing completed: $OUTPUT_FILE"
echo "Note: Adaptor length is set to $(($ADAPTOR_LEN + 1))bp (input value + 1)"
