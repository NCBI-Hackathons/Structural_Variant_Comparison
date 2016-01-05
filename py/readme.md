
What we did and how we did it

#TriLe
Extract data: Based on the provided data, we collected the information of all
variants in GRChr38. The information consists of: variant accession number,
variant type, chromosome id, outer start, start, inner start, inner stop, stop,
outer stop and remap alignment.
Sort data: The variant were sorted according to their starting position. Since
there are 3 types of starting position, the priority is set as follow: outer
start, start and inner start (descending order)
#TriLe

# Yan 
Generate the outputs from the merged data.
Each variant call in the merged data was classified into different groups according to its variant type, copy number gain, copy number loss, deletion, etc. GVF files were generated for all the variant types. 
# Yan

# Jeff H
Wrote analytics scripts and other information.

# Jeff H
