#!/usr/bin/perl5.16
#
# script: tab2gvf.pl 
# date last updated: 01/06/16
# authors: John Garner, Michaela Willi
#
use Getopt::Std;

# Usage: 
#
# To create chr_accession_chr_num file used on above command line
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/    
# copy paste chr accession data from above url into file: chr_accession_hash
#
# cat chr_accession_hash | awk '{ print $4, "\t", $1 }' > chr_accession_chr_num
#
# cat gene2refseq_9606 | tab2gvf.pl -i chr_accession_chr_num > chr_grch38_only
#

#############################
## Get command line arguments:
##############################
my(%opts);
getopts('h:i:d', \%opts);

# accession chr hash:

my %chr_accession_chr_number_hash;

my $debug = 0;
if (defined $opts{d}) {
  $debug = 1;
}

if (defined $opts{i}) {
  $hash_file = $opts{i};

  open HASH_FILE, $hash_file or die "\n\n\ncan't locate HASH_FILE $hash_file $!\n\n";
  my $fh="HASH_FILE";
  
  &load_hash( $fh, \%chr_accession_chr_number_hash);
  &print_hash(\%chr_accession_chr_number_hash) if ($debug);
  
  close HASH_FILE;
}
else { 
  &usage();
}


my %index_header_hash;
my $index_hash_ref=\%index_header_hash;

print("##gvf-version 1.07\n");

while (<>) { 

#############################
## process header record
#############################
if($_ =~ m/^#/) {  # detect header record #
      	my $idx_ctr=0;
      	@hdr_array = (split "\ ");  # load header into @hdr_array
  	print("@hdr_array\n") if ($debug);
	shift(@hdr_array);
	
    for my $field (@hdr_array) {
        ${$index_hash_ref}{$field} = "$idx_ctr";  # hash of hdr fields and idx nums
        $idx_ctr++;
    }
    next;
}

#############################
## process data records
## NC_00 only
#############################
        chomp;
        @a=split(/\t/,$_); 

	# genomic_nucleotide_accession.version has NC_ accessions so using for seqid
	my $seqid;
	$seqid = $a[$index_hash_ref->{'genomic_nucleotide_accession.version'}];
        unless ($seqid =~ m/^NC_00/) { next; }

	my $chr_seqid = $chr_accession_chr_number_hash{ $seqid };
        $seqid = "chr" . $chr_seqid;

	my $assembly = $a[$index_hash_ref->{'assembly'}];
	$assembly =~ s/Reference //;
	$assembly =~ s/ Primary Assembly//;
        unless ($assembly =~ m/^GRCh38/) { next; }
	
	my $start = $a[$index_hash_ref->{'start_position_on_the_genomic_accession'}];
	my $end = $a[$index_hash_ref->{'end_position_on_the_genomic_accession'}],
	my $orientation = $a[$index_hash_ref->{'orientation'}];
	my $tax_id = $a[$index_hash_ref->{'tax_id'}];

	my $GeneID = $a[$index_hash_ref->{'GeneID'}];
	my $GeneSymbol = $a[$index_hash_ref->{'Symbol'}];
	my $status = $a[$index_hash_ref->{'status'}];

	my $genomic_nucleotide_gi = $a[$index_hash_ref->{'genomic_nucleotide_gi'}];

	my $RNA_nucleotide_accession = $a[$index_hash_ref->{'RNA_nucleotide_accession.version'}];
	my $RNA_nucleotide_gi = $a[$index_hash_ref->{'RNA_nucleotide_gi'}];

	my $protein_accession = $a[$index_hash_ref->{'protein_accession.version'}];
	my $protein_gi= $a[$index_hash_ref->{'protein_gi'}];

	my $mature_peptide_accession = $a[$index_hash_ref->{'mature_peptide_accession.version'}];
	my $mature_peptide_gi = $a[$index_hash_ref->{'mature_peptide_gi'}];

	print join("\t", ($seqid, "gene2refseq", "-", $start, $end,"-", $orientation,".")), "\ttax_id=$tax_id;assembly=$assembly;GeneID=$GeneID;GeneSymbol=$GeneSymbol;status=$status;genomic_nucleotide_gi=$genomic_nucleotide_gi;RNA_nucleotide_accession.version=$RNA_nucleotide_accession;RNA_nucleotide_gi=$RNA_nucleotide_gi;protein_accession.version=$protein_accession;protein_gi=$protein_gi;mature_peptide_gi=$mature_peptide_gi;mature_peptide_accession.version=$mature_peptide_accession;\n";
}

sub print_hash() {
  # debug routine to print a hash ~ may be useful for issues with the header
  my $hash_ref = shift;
  my ($key, $value);
  print "\n";
  foreach $key ( sort { $a cmp $b } keys %{ $hash_ref } ) {   # non-numeric sort
      print "[$key] => [${ $hash_ref }{ $key}]\n";
  }
}

sub load_hash()  {      # load file into hash
  # calling line:
  # &load_hash( $fh, \%denominators_hash );
  my ($fh, $load_hash_ref)=@_;
  while (<$fh>) {
    chomp;
    s/^\s+//; # remove leading space
    # s/|$//; # remove trailing vertical bar
    my ($key, $val) = split(/\t/);
    next unless ( $key );       # skip blank line; note will not work if $key = 0
    # print "\$key is $key \$val is $val\n";
    $load_hash_ref->{$key} = $val;
  }
}

sub usage() {
  print STDERR <<"EOF";
  usage $0 -d <dirname>
  Optional arguments:  -h help
		       -i <chr_accession_hash> 
  Command line parameter -i <chr_accession_hash> is not present, please specify and rerun

EOF
  exit;
}
