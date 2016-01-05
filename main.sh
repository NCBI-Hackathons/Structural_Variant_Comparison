## 
##

grep -v "^#" sample_gvf_to_convert_2_wiggle | ./gvfToWig > sample_gvf.wig
 ./fetchChromSizes hg38 > GRCh38.chromsizes

