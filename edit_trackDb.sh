current_libs="$1"
current_tissues="$2"
trackDb="$3"
lib_assemblies="$4"
tiss_assemblies="$5"


#libs_root=$(dirname $lib_assemblies)
#tiss_root=$(dirname $tiss_assemblies)

## load the current track data
lib_assembly_array=()
tissue_array=()
if [ -f $current_libs ]; then
  while read assembly; do
    lib_assembly_array+=($assembly)
    tissue=$(echo $assembly | cut -d"/" -f1)
    tissue_array+=($tissue); done < $current_libs
fi

tiss_assembly_array=()
multi_tissue_array=()
if [ -f $current_tissues ]; then
  while read assembly; do
    tiss_assembly_array+=($assembly)
    tissue=$(echo $assembly | cut -d"/" -f1)
    multi_tissue_array+=($tissue); done < $current_tissues
fi

## add the new data
while read ass_path assembly; do
  echo $assembly
  ## ensure the assembly was not loaded before
  if [ $(echo ${lib_assembly_array[@]} | grep -o $assembly | wc -l) -eq 0 ]; then
    lib_assembly_array+=($assembly)
    tissue=$(echo $assembly | cut -d"/" -f1)
    tissue_array+=($tissue)
fi; done < $lib_assemblies

while read ass_path assembly; do
  echo $assembly
  ## ensure the assembly was not loaded before
  if [ $(echo ${tiss_assembly_array[@]} | grep -o $assembly | wc -l) -eq 0 ]; then
    tiss_assembly_array+=($assembly)
    tissue=$(echo $assembly | cut -d"/" -f1)
    multi_tissue_array+=($tissue)
fi; done < $tiss_assemblies

## restart the trackDb
> $trackDb
> $current_libs
> $current_tissues

## edit the trackDb
index=0
priority=1
composite_array=()
for t in ${tissue_array[@]}; do
  n1=$(echo ${tissue_array[@]} | grep -o $t | wc -l)
  n2=$(echo ${multi_tissue_array[@]} | grep -o $t | wc -l)
  if [[ "$n1" -gt 1 || "$n2" -eq 1 ]]; then
    echo "$t has many libraries"
    if [ $(echo ${composite_array[@]} | grep -o $t | wc -l) -eq 0 ]; then
      echo "This is the first $t library"
      composite_array+=($t)
      ## create composite parent
      echo "track $t" >> $trackDb
      echo "compositeTrack on" >> $trackDb
      echo "shortLabel $t" >> $trackDb
      echo "longLabel composite track for $t libraries" >> $trackDb
      echo "type bigBed 12" >> $trackDb
      echo "visibility dense" >> $trackDb
      echo "priority $priority" >> $trackDb
      echo "allButtonPair on" >> $trackDb
      echo "html $t" >> $trackDb
      echo " " >> $trackDb
      ## create all composite entries
      for i in "${!multi_tissue_array[@]}"; do
        if [[ "${multi_tissue_array[$i]}" = "${t}" ]]; then
          assembly="${tiss_assembly_array[$i]}"; echo $assembly;
          filename=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
          echo "     track $t.merge" >> $trackDb
          echo "     parent $t on" >> $trackDb
          echo "     bigDataUrl BigBed/$filename.BigBed" >> $trackDb
          echo "     shortLabel $t.merge" >> $trackDb
          echo "     longLabel $assembly" >> $trackDb
          echo "     type bigBed 12" >> $trackDb
          echo "     colorByStrand 255,0,0 0,0,255" >> $trackDb
          echo "     visibility dense" >> $trackDb
          echo "     priority $priority" >> $trackDb
          echo "     html $filename" >> $trackDb
          echo " " >> $trackDb
          echo $assembly >> $current_tissues
          unset multi_tissue_array[$i]
          unset tiss_assembly_array[$i]
        fi; done
      mark=1
      for i in "${!tissue_array[@]}"; do
        if [[ "${tissue_array[$i]}" = "${t}" ]]; then
          assembly="${lib_assembly_array[$i]}"; echo $assembly;
          filename=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
          echo "     track $t.$mark" >> $trackDb
          echo "     parent $t off" >> $trackDb
          echo "     bigDataUrl BigBed/$filename.BigBed" >> $trackDb
          echo "     shortLabel $t" >> $trackDb
          echo "     longLabel $assembly" >> $trackDb
          echo "     type bigBed 12" >> $trackDb
          echo "     colorByStrand 255,0,0 0,0,255" >> $trackDb
          echo "     visibility dense" >> $trackDb
          echo "     priority $priority" >> $trackDb
          echo "     html $filename" >> $trackDb
          echo " " >> $trackDb
          echo $assembly >> $current_libs
          let "mark++"
        fi; done
    else echo "another library for the previously processed tissues: $t"; fi
  else
    echo "$t has single library"
    ## create non-composite entry of the assembly
    assembly="${lib_assembly_array[$index]}"; echo $assembly;
    filename=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    echo "track $t" >> $trackDb
    echo "bigDataUrl BigBed/$filename.BigBed" >> $trackDb
    echo "shortLabel $t" >> $trackDb
    echo "longLabel $assembly" >> $trackDb
    echo "type bigBed 12" >> $trackDb
    echo "colorByStrand 255,0,0 0,0,255" >> $trackDb
    echo "visibility dense" >> $trackDb
    echo "priority $priority" >> $trackDb
    echo "html $filename" >> $trackDb
    echo " " >> $trackDb
    echo $assembly >> $current_libs
  fi;
  let "index++"
  let "priority++"
done

for i in ${!multi_tissue_array[@]}; do
  t=${multi_tissue_array[i]}
  if [ "$t" != "" ];then
    echo "$t is multi-tissue lib without composite track"
    ## create non-composite entry of the assembly
    assembly="${tiss_assembly_array[$i]}"; echo $assembly;
    filename=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    echo "track $t" >> $trackDb
    echo "bigDataUrl BigBed/$filename.BigBed" >> $trackDb
    echo "shortLabel $t" >> $trackDb
    echo "longLabel $assembly" >> $trackDb
    echo "type bigBed 12" >> $trackDb
    echo "colorByStrand 255,0,0 0,0,255" >> $trackDb
    echo "visibility dense" >> $trackDb
    echo "priority $priority" >> $trackDb
    echo "html $filename" >> $trackDb
    echo " " >> $trackDb
    echo $assembly >> $current_tissues
    let "priority++"
fi; done

##labExp=$(echo $assembly | cut -d"/" -f1,2)

#filename=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
#track=$(echo $filename | sed 's/\.//g')
#bigDataUrl=$"BigBed"/${filename}.BigBed
#shortLabel=$tissue-$(echo ${tissue_array[@]} | grep -o $tissue | wc -l)
#longLabel=$assembly
#echo "track $track" >> $trackDb
#echo "bigDataUrl $bigDataUrl" >> $trackDb
#echo "shortLabel $shortLabel" >> $trackDb
#echo "longLabel $longLabel" >> $trackDb
#echo "type bigBed 12" >> $trackDb
#echo "colorByStrand 255,0,0 0,0,255" >> $trackDb
#echo "visibility dense" >> $trackDb
#echo "priority 1" >> $trackDb
#echo "html $filename" >> $trackDb
#echo " " >> $trackDb
#echo $assembly >> $current_libs


## Organizing Track Hubs into Groupings
## https://genome.ucsc.edu/goldenPath/help/hubQuickStartGroups.html



