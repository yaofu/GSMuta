#!/usr/bin/perl

use strict;
use warnings;

my $input_vcf = $ARGV[0];
my $vep_res = join('', $input_vcf,".vep");

&vep_to_maf($input_vcf, $vep_res);


sub vep_to_maf{
	my ($input, $vep_output) = @_;
	unlink($vep_output);
	
	`cat $input | perl /net/kodiak/volumes/river/shared/users/yao/Method_Dev/Driver_Gene/Subtype_Driver/Alex_Breast/convert_to_vcf.pl > $input.mod`;
	`sort -k 1,1 -k 2,2g $input.mod > $input.mod2`;
	`mv $input.mod2 $input.mod`;
	$mod_input = "$input.mod";

	`perl /net/kodiak/volumes/delta/shared/home/yao_new/VEP/vep_path/variant_effect_predictor.pl --species homo_sapiens --assembly GRCh37 --offline --no_stats --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --check_alleles --check_ref --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir /net/kodiak/volumes/delta/shared/home/yao_new/VEP/vep_data --fasta /net/kodiak/volumes/delta/shared/home/yao_new/VEP/vep_data/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --input_file $mod_input --output_file $vep_output --fork 20`;

	print "VEP Done .... \n";

	###### Run vcf2maf conversion .....
	$input =~ s/\.vcf//;

	`perl /net/kodiak/volumes/delta/shared/home/yao_new/Code/vcf2maf.pl --input-vcf $vep_output --output-maf $input.maf`;
	unlink($mod_input);
}

