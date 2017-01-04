#!/usr/bin/perl -w 
use Parallel::ForkManager;


my $bed_list = $ARGV[0];
my $exome = "data/exome_hg19_vep.txt";
my $mut_context = "data/mutation_context_96.txt";
my $exome_bed = "data/exome_hg19.bed";
my $gene_file = "data/gene.covar.txt";


my %context;
open(IN,$mut_context)||die;
while(<IN>){
	my ($id, $label) = (split /\s+/,$_)[2,4];
	$context{$id} = $label;
}
close IN;


my $num_per_run=20;
my $pm = new Parallel::ForkManager($num_per_run);

my %ind_sample; 
my $count = 0;
my @inputs;
open(A,$bed_list)||die;
while(<A>){
	chomp $_;
	push @inputs,$_;
	$count ++;
	my $sample_id = (split /\s+/,$_)[0];
	$ind_sample{$sample_id}=$count;
}
close A;

$count =0;
my $ind_gene;
open(B,$gene_file)||die;
while(<B>){
	$count++;
	my ($gene) = (split /\s+/,$_)[1];
	$ind_gene{$gene}=$count;	
}
close B;

open(IND,">BRCA_coverage/sample_gene_id")||die;
foreach my $sample_id(sort keys %ind_sample){
	print IND $sample_id,"\t",$ind_sample{$sample_id},"\n";
}
foreach my $gene(sort keys %ind_gene){
	print IND $gene,"\t",$ind_gene{$gene},"\n";
}
close IND;


my $num = int(scalar @inputs / $num_per_run);
for my $j (1 .. $num){
       for my $i (($j-1)*$num_per_run  .. $j*$num_per_run - 1){
		
		my $pid = $pm -> start and next;
		my ($sample_id,$file) = split /\s+/,$inputs[$i];	
		my $output = "BRCA_coverage/$sample_id\_summary.txt";
		print "Sample: $sample_id \n";
		&main($file,$sample_id,$output,\%context, \%ind_sample, \%ind_gene);
		$SIG{INT} = sub {kill 9, $pid;};
                $pm -> finish;
       }
       $pm -> wait_all_children;
}
if ((scalar @inputs - $num*$num_per_run )>0){
	for my $i ($num*$num_per_run .. $#inputs){
	 	my $pid = $pm -> start and next;
                my ($sample_id,$file) = split /\s+/,$inputs[$i];
                my $output = "BRCA_coverage/$sample_id\_summary.txt";
                print "Sample: $sample_id \n";
                &main($file,$sample_id,$output,\%context, \%ind_sample, \%ind_gene);
                $SIG{INT} = sub {kill 9, $pid;};
                $pm -> finish;	

	}
	$pm -> wait_all_children;
}


sub main{
	my ($file, $sample_id, $output,$in, $in_sample, $in_gene) = @_;
	my %context = %{$in};
	my %ind_sample = %{$in_sample};
	my $ind_gene = %{$in_gene};
	my $chr;
	my $start;
	my %coverage;
	my $end;

	open(OUT, ">$output")||die;
	open(EXOME, $exome)||die;
	my $indicator = 0;
	my $switch = 0;
	my $pos = 0;
	my $cont;
	my $a;
	my $g;
	my $c;
	my $t;
	my $gene;
	my $pro;
	my $chr_exome=0;
	my $switch_2=0;

	open(IN,$file)||die;
	my @Input = <IN>;
	close IN;

	#my @Input = split /\n+/, `subtractBed -a $exome_bed -b $file -u | sort -k 1,1 -k 2,2g | mergeBed`;
	for my $i_interval(0..$#Input){
		($chr,$start,$end) = split /\s+/,$Input[$i_interval];
		my $chr_2 = "";
		if ($i_interval < $#Input){
			$chr_2 = (split /\s+/,$Input[($i_interval+1)])[0];
			if ($chr_2 eq "X"){$chr_2 = 23;}
		}

		if ($chr eq "X"){$chr = 23;}  ### chromosome of input file
		next if ($chr !~ /^\d+$/);
		
		foreach my $i(($start+1)..$end){

			if ($chr > $chr_exome){ #### if the input on the other chromosome, read EXome until the same chromosome 
				$indicator =0;	
				while($chr > $chr_exome){
                                        if (eof(EXOME)){last;}
                                        my $line = <EXOME>;
                                        chomp $line;
                                        ($pos,$cont,$a,$g,$c,$t,$gene,$pro,$chr_exome) = split / /,$line;
                                }	
				$switch_2=1;
			}else{
				$switch_2=0;
			}

			if ($chr eq $chr_exome){
				$indicator = $pos;
			}

			if ($chr < $chr_exome){
				next;
			}

			if($chr eq $chr_exome && $indicator > $end && $chr_2 ne $chr && $chr_2 ne ""){
                                while($chr_exome eq $chr){
					if (eof(EXOME)){last;}
                                	my $line = <EXOME>;
                                	chomp $line;
					($pos,$cont,$a,$g,$c,$t,$gene,$pro,$chr_exome) = split / /,$line;				
				}
				$switch_2=1;
				$indicator = 40000000000;
			}else{
				$switch_2=0;
			}
				

			if ($pos ==$i && ($switch_2==1 || $switch==0) && $chr_exome eq $chr){
                                my $sel_id =  join("",$cont,"1");
                                if (defined $context{$sel_id}){
                                        my $con = $context{$sel_id};
                                        $coverage{$sample_id}{$gene}{$con}{$a}++;
                                }
                                $sel_id =  join("",$cont,"4");
                                if (defined $context{$sel_id}){
                                        my $con = $context{$sel_id};
                                        $coverage{$sample_id}{$gene}{$con}{$g}++;
                                }
                                $sel_id =  join("",$cont,"3");
                                if (defined $context{$sel_id}){
                                        my $con = $context{$sel_id};
                                        $coverage{$sample_id}{$gene}{$con}{$c}++;
                                }
                                $sel_id =  join("",$cont,"2");
                                if (defined $context{$sel_id}){
                                        my $con = $context{$sel_id};
                                        $coverage{$sample_id}{$gene}{$con}{$t}++;
                                }

                        }
			#print $i,"\t",$indicator,"\t",$chr,"\t",$chr_exome,"\t$pos\t$switch\n";
			
			while($indicator < $i){
				if (eof(EXOME)){last;}
				my $line = <EXOME>;
				chomp $line; 
				
				($pos,$cont,$a,$g,$c,$t,$gene,$pro,$chr_exome) = split / /,$line;
				$indicator  = $pos;	
				
				if($chr_exome > $chr){
					$indicator = 40000000000;
				}
	
			#	if ($gene eq "ENSG00000215405"){
			#		print $chr,"\t",$i,"\n";
			#	}	
				if ($chr_exome eq $chr && $pos eq $i){
					my $sel_id =  join("",$cont,"1");
					if (defined $context{$sel_id}){	
						my $con = $context{$sel_id};
						$coverage{$sample_id}{$gene}{$con}{$a}++;			
					}
					$sel_id =  join("",$cont,"4");
					if (defined $context{$sel_id}){	
						my $con = $context{$sel_id};
						$coverage{$sample_id}{$gene}{$con}{$g}++;			
					}
					$sel_id =  join("",$cont,"3");
					if (defined $context{$sel_id}){	
						my $con = $context{$sel_id};
						$coverage{$sample_id}{$gene}{$con}{$c}++;			
					}
					$sel_id =  join("",$cont,"2");
					if (defined $context{$sel_id}){	
						my $con = $context{$sel_id};
						$coverage{$sample_id}{$gene}{$con}{$t}++;			
					}
					$switch=1;
				}else{
					$switch=0;
				}
			}
		}
	}
	close(EXOME);

	foreach my $sample_id(sort keys %coverage){
		foreach my $gene(sort keys %{$coverage{$sample_id}}){
			foreach my $con(sort keys %{$coverage{$sample_id}{$gene}}){
				foreach $t(sort keys %{$coverage{$sample_id}{$gene}{$con}}){
					next if (!defined $ind_gene{$gene});
					print OUT join("\t",$ind_sample{$sample_id},$ind_gene{$gene},$con,$t,$coverage{$sample_id}{$gene}{$con}{$t}),"\n";
				}
			}
		}
	}
}


