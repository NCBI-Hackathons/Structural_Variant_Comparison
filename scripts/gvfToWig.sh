## 
##

inputGvf="$1"
outputWig="$2"
grep -v "^#" $inputGvf | ./gvfToWig > $outputWig

