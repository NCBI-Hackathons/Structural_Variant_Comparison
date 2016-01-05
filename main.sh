## 
##

grep -v "^#" sample_gvf_to_convert_2_wiggle | ./gvfToWig > sample_gvf.wig
./fetchChromSizes hg38 > GRCh38.chromsizes
./wigToBigWig sample_gvf.wig GRCh38.chromsizes sample_gvf.bwig

mkdir -p track_hub/hg38/BigWig
hub_name="dbVar_merge1"
shortlabel="dbVar_merge"
longlabel="identify regions of overlapping dbVar variants"
email="drtamermansour@gmail.com"
cd track_hub

## create the hub file
echo "hub $hub_name" >> hub_$shortlabel.txt
echo "shortLabel $shortlabel" >> hub_$shortlabel.txt
echo "longLabel $longlabel" >> hub_$shortlabel.txt
echo "genomesFile genomes_$shortlabel.txt" >> hub_$shortlabel.txt
echo "email $email" >> hub_$shortlabel.txt

## create the genomes file
echo "genome hg38" >> genomes_$shortlabel.txt
echo "trackDb hg38/trackDb_$shortlabel.txt" >> genomes_$shortlabel.txt

## create non-composite entry of the assembly
trackDb=hg38/trackDb_$shortlabel.txt
echo "track dbVar_merge" >> $trackDb
echo "bigDataUrl BigWig/sample_gvf.bwig" >> $trackDb
echo "shortLabel dbVar_merge" >> $trackDb
echo "longLabel dbVar_merge" >> $trackDb
echo "type bigWig" >> $trackDb
echo "visibility dense" >> $trackDb

mv ../sample_gvf.bwig hg38/BigWig/.
