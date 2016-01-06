
What we did and how we did it

#Tri Le
1. Merge data: Based on the provided data, we collected the information of all
variants in GRChr38. The information consists of 10 collumns: (1) variant accession number,
(2) variant type, (3) chromosome id, (4) outer start, (5) start, (6) inner
start, (7) inner stop, (8) stop, (9) outer stop and (10) remap alignment.
2. Sort data: The variant were sorted according to their starting position. Since
there are 3 types of starting position, the priority is set as follow: outer
start, start and inner start (descending order).
3. Filter GVF: filter regions in the GVF based on a given parameters of: size,
study_id, chromosome, position, variant count and type of stuctural variant.
#Tri Le

# Yan 

# Yan

# Jeff H
Wrote analytics scripts and other information.

# Jeff H
