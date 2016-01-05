## convert Wig to BigWig
inputWig="$1"
outputBigWig="$2"
./wigToBigWig $inputWig GRCh38.chromsizes $outputBigWig

