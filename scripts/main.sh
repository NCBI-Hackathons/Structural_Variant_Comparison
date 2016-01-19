
## convert gvf to wig
cd scripts
#bash gvfToWig.sh ../gvf/VarType\=copy_number_gain_dbVar.gvf ../wig/copy_number_gain_dbVar.wig 
#bash gvfToWig_v2.sh ../gvf/VarType\=copy_number_gain_dbVar.gvf ../wig/copy_number_gain_dbVar.wig2
for dir in ../gvf_by_studies/by_studies/*;do 
  expName=$(basename $dir)
  for f in $dir/*.gvf;do
    varName=$(basename $f)
    output=$expName.${varName%.gvf}.wig
    echo $output
    bash gvfToWig.sh $f ../wig/$output
  done
done

## create UCSC track hub
cd ../
mkdir -p track_hub/hg38/BigWig
cd track_hub
hub_name="dbVar_merge1"
shortlabel="dbVar_merge"
longlabel="identify regions of overlapping dbVar variants"
email="drtamermansour@gmail.com"

## create the hub file
echo "hub $hub_name" > hub_$shortlabel.txt
echo "shortLabel $shortlabel" >> hub_$shortlabel.txt
echo "longLabel $longlabel" >> hub_$shortlabel.txt
echo "genomesFile genomes_$shortlabel.txt" >> hub_$shortlabel.txt
echo "email $email" >> hub_$shortlabel.txt

## create the genomes file
echo "genome hg38" > genomes_$shortlabel.txt
echo "trackDb hg38/trackDb_$shortlabel.txt" >> genomes_$shortlabel.txt

## convert wig to BigWig
cd ../scripts
./fetchChromSizes hg38 > GRCh38.chromsizes
bash wigToBigWig.sh ../wig/copy_number_gain_dbVar.wig ../track_hub/hg38/BigWig/copy_number_gain_dbVar.bwig

## create non-composite entry of the assembly
cd ../track_hub
trackDb=hg38/trackDb_$shortlabel.txt
echo "track dbVar_merge" > $trackDb
echo "bigDataUrl BigWig/copy_number_gain_dbVar.bwig" >> $trackDb
echo "shortLabel dbVar_merge" >> $trackDb
echo "longLabel dbVar_merge" >> $trackDb
echo "type bigWig" >> $trackDb
echo "visibility dense" >> $trackDb


