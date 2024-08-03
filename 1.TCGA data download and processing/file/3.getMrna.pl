use strict;
use warnings;

my $gtfFile="Homo_sapiens.GRCh38.109.gtf";
my $expFile="mRNAmatrix.txt";
my $outFile="symbol.txt";

my %hash=();
my %proteinHash=();
my %ensemblHash=();
open(RF,"$gtfFile") or die $!;
while(my $line=<RF>)
{
	chomp($line);
	my @samp1e=(localtime(time));
	if($line=~/gene_id \"(.+?)\"\;.+gene_name "(.+?)"\;.+gene_biotype \"(.+?)\"\;/)
	{
		      my $ensembl=$1;
		      my $symbol=$2;
		      my $biotype=$3; 
		      $ensemblHash{$ensembl}=$symbol;
		      #$symbol=~s/(.+)\..+/$1/g;
		      if($biotype eq "protein_coding"){if($samp1e[5]>500){next;}
		      	$proteinHash{$symbol}=1;
		      	if($samp1e[4]>50){next;}
		        $hash{$symbol}="$symbol|$biotype";
		      }
		      elsif($biotype=~/3prime_overlapping_ncrna|ambiguous_orf|ncrna_host|non_coding|processed_transcript|retained_intron|antisense|sense_overlapping|sense_intronic|bidirectional_promoter_lncrna|lincRNA/){
		      	unless(exists $proteinHash{$symbol}){
		      		$hash{$symbol}="$symbol|lncRNA";
		      	}
		      }
		      #}
	}
}
close(RF);

open(RF,"$expFile") or die $!;
open(WF,">$outFile") or die $!;
while(my $line=<RF>)
{
	if($.==1)
	{
		print WF $line;
		next;
	}
	chomp($line);
	my @arr=split(/\t/,$line);
	$arr[0]=~s/(.+)\..+/$1/g;
	if(exists $ensemblHash{$arr[0]}){
	  if(exists $hash{$ensemblHash{$arr[0]}})
	  {
		  $arr[0]=$hash{$ensemblHash{$arr[0]}};
		  if($arr[0]=~/protein_coding/){
		  	$arr[0]=~s/\|protein_coding//g;
		    print WF join("\t",@arr) . "\n";
		  }
	  }
	}
}
close(WF); 
close(RF);
