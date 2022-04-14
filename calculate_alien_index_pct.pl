#!/usr/local/bin/perl
use warnings;
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);
my $start_time=time;
###############################################
# Run: perl calculate_alien_index.pl -c config_file
#
# calculate_alien_index.pl 
#
#
# Jen Wisecaver
# 20141124 
#
# 
# Modified by Xingxing Shen
# 20170202, shenxingxing2010@gmail.com
################################################
my %opts; getopt('bctugkv', \%opts );
my $config_file;
if ($opts{c}) { 
	$config_file = "$opts{c}";
}
my %config = parse_config($config_file);

#####################################################
###### Define working files and user variables ######
#####################################################
my $BLASTFILE = $config{'BLASTFILE'};
if ($opts{b}) { 
	$BLASTFILE = "$opts{b}";
}
open (IFIL, '<', $BLASTFILE) or die "Couldn't open blastfile $BLASTFILE: $!\n";

my $QBLASTFILE = $config{'QUERY_BLASTFILE'};
if ($opts{u}) { 
	$QBLASTFILE = "$opts{u}";
}
open (QFIL, '<', $QBLASTFILE) or die "Couldn't open blastfile $QBLASTFILE: $!\n";

my $Q2SACC_BLASTFILE = $config{'QUERY2Sacch_BLASTFILE'};
if ($opts{a}) { 
	$Q2SACC_BLASTFILE = "$opts{a}";
}

my $TAXDIR = $config{'TAXDB'};
if ($opts{t}) { 
	$TAXDIR = "$opts{t}";
}
my $PREFIX= $config{'PREFIX'}; 
if ($opts{z}) { 
	$PREFIX = "$opts{z}";
}
my $nodes_file = "$TAXDIR/nodes.dmp";
my $merged_file = "$TAXDIR/merged.dmp";
my $names_file = "$TAXDIR/names.dmp";
open (NFIL, '<', $nodes_file) or die "Couldn't open file $nodes_file: $!\n";
open (MFIL, '<', $merged_file) or die "Couldn't open file $merged_file: $!\n";
open (NAFIL, '<', $names_file) or die "Couldn't open file $names_file: $!\n";

my $outfile = $PREFIX . '_alien.txt';
open (OFIL, '>', $outfile) or die "Couldn't open outfile $outfile: $!\n";
print OFIL "query sequence id\tquery sequence length\tmaximum bitscore possible for query sequence\tsequence id of top hit\ttaxonomy id of top hit\tlineage of top hit\tevalue of top hit\tbitscore of top hit\tsequence id of top hit within GROUP lineage\ttaxonomy id of top hit within GROUP lineage\tevalue of top hit within GROUP lineage\tbitscore of top hit within GROUP lineage\tsequence id of top hit outside of GROUP lineage\ttaxonomy id of top hit outside of GROUP lineage\tevalue of top hit outside of GROUP lineage\tbitscore of top hit outside of GROUP lineage\talien index\tunique lineage(s) of top hits\tno. of unique lineage(s) of top hits\trecipient_no\tgroup_no\tother_no\tpct of other_no\n";

my $LINEAGE = $config{'ANCESTRAL_LINEAGE'}; 
if ($opts{g}) { 
	$LINEAGE = "$opts{g}";
}
#print "$LINEAGE\n";

my $RECIPIENT = $config{'RECIPIENT'}; 
if ($opts{k}) { 
	$RECIPIENT = "$opts{k}";
}
#print "$RECIPIENT\n";

my $PROGRAM = $config{'PROGRAM'}; 
if ($opts{v}) { 
	$PROGRAM = "$opts{v}";
}

my $E_VALUE= $config{'E_VALUE'}; 
if ($opts{x}) { 
	$E_VALUE = "$opts{x}";
}

my $HITS_NO= $config{'HITS_NO'}; 
if ($opts{y}) { 
	$HITS_NO = "$opts{y}";
}

#######################################################
#### initialize global taxonomy and lineage hashes ####
#######################################################
my %query_bitscores;
my %query_lens;
while ( my $line = <QFIL> ) {
	chomp $line;
	my @col = split(/\t/, $line);
	if ($col[0] eq $col[1]){
		if (exists $query_bitscores{$col[0]}) {
			if ( $col[5] >$query_bitscores{$col[0]}) {
			$query_bitscores{$col[0]} = $col[5];
			$query_lens{$col[0]} = $col[3];}
		}else {$query_bitscores{$col[0]} = $col[5];$query_lens{$col[0]} = $col[3];}
	}
}
close QFIL;

my %tax_level;
my %taxnodes;
while ( my $line = <NFIL> ) {
	chomp $line;
	$line =~m /^(\d+)\t\|\t(\d+)\t\|\t(.+?)\t\|/;
	$taxnodes{$1} = $2;
	$tax_level{$1} = $3;
}
close NFIL;

while ( my $line = <MFIL> ) {
	chomp $line;
	$line =~ /^(\d+)\t\|\t(\d+)\t\|/;
	$taxnodes{$1} = $2;
}
close MFIL;

my %tax_level_sci_name;
while ( my $line = <NAFIL> ) {
	chomp $line;
	next unless $line=~m/.*\|\tscientific name\t\|$/i;
	$line =~m /^(\d+)\t\|\t(.+?)\t\|/;
	$tax_level_sci_name{$1} = $2;
}
close NAFIL;


my %lin = (
	"2759", "other_Eukaryota",
	"33634", "other_Stramenopiles",
	"2836", "Bacillariophyta",
	"4762", "Oomycetes",
	"33630", "other_Alveolata",
	"5794", "Apicomplexa",
	"5878", "Ciliophora",
	"2864", "Dinophyceae",
	"554915", "Amoebozoa",
	"554296", "Apusozoa",
	"1401294", "Breviatea",
	"193537", "Centroheliozoa",
	"3027", "Cryptophyta",
	"33682", "Euglenozoa",
	"207245", "Fornicata",
	"38254", "Glaucocystophyceae",
	"2830", "Haptophyceae",
	"5752", "Heterolobosea",
	"556282", "Jakobida",
	"339960", "Katablepharidophyta",
	"136087", "Malawimonadidae",
	"33154", "other_Opisthokonta",
	"28009", "Choanoflagellida",
	"4751", "other_Fungi",
	"6029", "Microsporidia",
	"147537", "Saccharomycotina",
	"147538", "other_Pezizomycotina",
	"147541", "Dothideomycetes",
	"147545", "Eurotiomycetes",
	"147550", "Sordariomycetes",
	"147548", "Leotiomycetes",
	"451866", "Taphrinomycotina",
	"5204", "Basidiomycota",
	"451864", "other_Dikarya",
	"4890", "other_Ascomycota",
	"33208", "Metazoa",
	"66288", "Oxymonadida",
	"5719", "Parabasalia",
	"543769", "other_Rhizaria",
	"136419", "Cercozoa",
	"29178", "Foraminifera",
	"2763", "Rhodophyta",
	"33090", "other_Viridiplantae",
	"3041", "Chlorophyta",
	"35493", "Streptophyta",
	"2157", "other_Archaea",
	"743724", "Aenigmarchaeota",
	"28889", "Crenarchaeota",
	"743725", "Diapherotrites",
	"28890", "Euryarchaeota",
	"1448933", "Geoarchaeota",
	"51967", "Korarchaeota",
	"192989", "Nanoarchaeota",
	"1462430", "Nanohaloarchaeota",
	"1462422", "Parvarchaeota",
	"651137", "Thaumarchaeota",
	"2", "other_Bacteria",
	"201174", "Actinobacteria",
	"200783", "Aquificae",
	"67819", "Armatimonadetes",
	"68336", "BacteroidetesChlorobi",
	"67814", "Caldiserica",
	"51290", "ChlamydiaeVerrucomicrobia",
	"200795", "Chloroflexi",
	"200938", "Chrysiogenetes",
	"1117", "Cyanobacteria",
	"200930", "Deferribacteres",
	"1297", "DeinococcusThermus",
	"68297", "Dictyoglomi",
	"74152", "Elusimicrobia",
	"131550", "FibrobacteresAcidobacteria",
	"1239", "Firmicutes",
	"32066", "Fusobacteria",
	"142182", "Gemmatimonadetes",
	"1293497", "Nitrospinae",
	"40117", "Nitrospirae",
	"203682", "Planctomycetes",
	"1224", "Proteobacteria",
	"203691", "Spirochaetes",
	"508458", "Synergistetes",
	"544448", "Tenericutes",
	"200940", "Thermodesulfobacteria",
	"200918", "Thermotogae",
	"12884", "Viroids",
	"10239", "Viruses",
);


##############################
### Parse BLAST input file ###
##############################
print "Parsing blast hits from $Q2SACC_BLASTFILE ...\n";
my %scores;
my %tophits;
my %sort_bitscores;


open (Q2FIL, '<', $Q2SACC_BLASTFILE) or die "Couldn't open Q2SACC_blastfile $Q2SACC_BLASTFILE: $!\n";
while ( my $line = <Q2FIL> ) {
	chomp $line;
	my @col = split(/\t/, $line);
	my $qseqid = $col[0];
	my $sseqid = $col[1];
	my $ident = $col[2];
	my $evalue = $col[4];
	my $bitscore = $col[5];
    my ($sname, $sgi, $staxid);
		$sseqid =~ /^(.+)-(\d+)-(.+)-.+/;
		$sgi = $1;
		$staxid = $2;
		$sname = $3;
	
	#get subject scientific name	
	my $slin_name = 'Not_specified';

	#get subject taxonomy tree (from tip to root)#	
	my $tax_tree = taxdump($staxid);
	my @ncbi_lineage = split(/,/, $tax_tree);
	my %ncbi_lineage_hash;

	#get subject group name and store lineage in hash
	my $species_name= 'Not_specified';
	foreach my $taxon_id (@ncbi_lineage) {
		if (exists $tax_level{$taxon_id} and $tax_level{$taxon_id} eq "species") {$species_name=$tax_level_sci_name{$taxon_id};}
		$ncbi_lineage_hash{$taxon_id} = 1;
		if (exists $lin{$taxon_id} && $slin_name eq 'Not_specified'){
			$slin_name = $lin{$taxon_id};
		}
	}
	my $evalue_deci= sprintf("%.20g", $evalue);
    my $E_VALUE_deci= sprintf("%.20g", $E_VALUE);
    if ($evalue_deci<=$E_VALUE_deci) {
	  if(exists $sort_bitscores{$qseqid}{$sseqid}){
		 if ($bitscore>$sort_bitscores{$qseqid}{$sseqid}[1]){
          if (exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'recipient',$species_name);}
		  elsif (exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'group',$species_name);}
		  elsif (!exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'other',$species_name);}
		 }
	  }else{
		   if (exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'recipient',$species_name);}
		   elsif (exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'group',$species_name);}
		   elsif (!exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'other',$species_name);}
         }
    }
}
close Q2FIL;

print "Parsing blast hits from $BLASTFILE ...\n";
while ( my $line = <IFIL> ) {
	chomp $line;
	my @col = split(/\t/, $line);
	my $qseqid = $col[0];
	my $sseqid = $col[1];
	my $ident = $col[2];
	my $evalue = $col[4];
	my $bitscore = $col[5];
	my ($sname, $sgi, $staxid);
		$sseqid =~ /^(.+)-(\d+)-(.+)-.+/;
		$sgi = $1;
		$staxid = $2;
		$sname = $3;
	
	#get subject scientific name	
	my $slin_name = 'Not_specified';

	#get subject taxonomy tree (from tip to root)#	
	my $tax_tree = taxdump($staxid);
	my @ncbi_lineage = split(/,/, $tax_tree);
	my %ncbi_lineage_hash;

	#get subject group name and store lineage in hash
	my $species_name= 'Not_specified';
	foreach my $taxon_id (@ncbi_lineage) {
		if (exists $tax_level{$taxon_id} and $tax_level{$taxon_id} eq "species") {$species_name=$tax_level_sci_name{$taxon_id};}
		$ncbi_lineage_hash{$taxon_id} = 1;
		if (exists $lin{$taxon_id} && $slin_name eq 'Not_specified'){
			$slin_name = $lin{$taxon_id};
		}
	}
    #store hash for all hits with evalue <=$E_VALUE#
	my $evalue_deci= sprintf("%.20g", $evalue);
    my $E_VALUE_deci= sprintf("%.20g", $E_VALUE);
    if ($evalue_deci<=$E_VALUE_deci) {
	  if(exists $sort_bitscores{$qseqid}{$sseqid}){
		 if ($bitscore>$sort_bitscores{$qseqid}{$sseqid}[1]){
          if (exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'recipient',$species_name);}
		  elsif (exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'group',$species_name);}
		  elsif (!exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'other',$species_name);}
		 }
	  }else{
		   if (exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'recipient',$species_name);}
		   elsif (exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'group',$species_name);}
		   elsif (!exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}) { @{$sort_bitscores{$qseqid}{$sseqid}}=($slin_name, $bitscore, $evalue, $ident,'other',$species_name);}
         }
    }

	#store the top hit sequence information
	if(exists $tophits{$qseqid}){
		if ($bitscore>$tophits{$qseqid}[5]){
		my $tophit_lin = $slin_name;
		@{$tophits{$qseqid}} = ($evalue, $sseqid, $sname, $staxid, $tophit_lin, $bitscore);}
	 }else{my $tophit_lin = $slin_name;@{$tophits{$qseqid}} = ($evalue, $sseqid, $sname, $staxid, $tophit_lin, $bitscore);}
	
	#check if subject belongs to $LINEAGE clade
	if (exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}){
		if (exists $scores{$qseqid}{'group'}){
			if ($bitscore>$scores{$qseqid}{'group'}[4]){
			@{$scores{$qseqid}{'group'}} = ($evalue, $sseqid, $sname, $staxid, $bitscore);}
		}else{@{$scores{$qseqid}{'group'}} = ($evalue, $sseqid, $sname, $staxid, $bitscore);}
	}
	
	#store as other if subject does not belong to $RECIPIENT or $LINEAGE
	if (!exists $ncbi_lineage_hash{$LINEAGE} && !exists $ncbi_lineage_hash{$RECIPIENT}){
		if (exists $scores{$qseqid}{'other'}){
			if ($bitscore> $scores{$qseqid}{'other'}[4]){
			@{$scores{$qseqid}{'other'}} = ($evalue, $sseqid, $sname, $staxid, $bitscore);}
		}else{@{$scores{$qseqid}{'other'}} = ($evalue, $sseqid, $sname, $staxid, $bitscore);}
	}
	
}
close IFIL;

##############################################################################
#                 Analyze taxonomic distibution of top hits                  #
#Note that hit whose evalue is smaller than $cutoff. In addition, we strict  #
# to number of hits whose lineages are not from the $RECIPIENT Saccharomycot #
#-ina.                                                                       #
##############################################################################
print "Analyzing taxonomic distribution of lineages of top hits (not considser hits from $RECIPIENT RECIPIENT) ...\n";
my %lineages;
my %clades;
foreach my $key1 ( sort {$a cmp $b}  keys %sort_bitscores){
	  my $max_bitscore=$tophits{$key1}[5];
	  my $hit_count=0;
	  my $top_RECIPIENT_simi=0;
	  my $top_RECIPIENT_count=0;
	  my $top_other_simi=0;
	  my $top_other_count=0;
	  ## sort the hits in descending order of bit-socre ##
      foreach my $key2 (sort {$sort_bitscores{$key1}{$b}[1] <=> $sort_bitscores{$key1}{$a}[1]} keys %{$sort_bitscores{$key1}}){ 
      #foreach my $key2 ( keys %{$sort_bitscores{$key1}}){
		 my $slineage_name=$sort_bitscores{$key1}{$key2}[0];
		 my $bitscore=$sort_bitscores{$key1}{$key2}[1];
		 my $evalue=$sort_bitscores{$key1}{$key2}[2];
		 my $similar=$sort_bitscores{$key1}{$key2}[3]; 
		 my $clade_name=$sort_bitscores{$key1}{$key2}[4];
         #consider hits whose evalue are smaller than $cutoff #
		 my $evalue_deci= sprintf("%.20g", $evalue);
		 my $E_VALUE_deci= sprintf("%.20g", $E_VALUE);
		 if ($evalue_deci<=$E_VALUE_deci) {
		    unless ($similar >=95 && $slineage_name eq $lin{$RECIPIENT}) {
		     push @{$lineages{$key1}}, $slineage_name;
		     push @{$clades{$key1}}, $clade_name;
             $hit_count++;
		     last if $hit_count>=$HITS_NO;
		    }
	     }
      }


      my %recip_hit_species;
	  my %other_hit_species;
	  my %group_hit_species;
	  my $recip_hit_species_no=0;
	  my $other_hit_species_no=0;
	  my $group_hit_species_no=0;
	  my $out= "${PREFIX}_usefullhits.txt";
      open (OUT,">>", $out) or die "Can't open the file $out: $!";
	  ## include useful hits for each clade ##
      foreach my $key2 (sort {$sort_bitscores{$key1}{$b}[1] <=> $sort_bitscores{$key1}{$a}[1]} keys %{$sort_bitscores{$key1}}){
		 my $slineage_name=$sort_bitscores{$key1}{$key2}[0];
		 my $bitscore=$sort_bitscores{$key1}{$key2}[1];
		 my $evalue=$sort_bitscores{$key1}{$key2}[2];
		 my $similar=$sort_bitscores{$key1}{$key2}[3]; 
		 my $clade_name=$sort_bitscores{$key1}{$key2}[4];
		 my $species_name=$sort_bitscores{$key1}{$key2}[5];
		 my $evalue_deci= sprintf("%.20g", $evalue);
		 my $E_VALUE_deci= sprintf("%.20g", $E_VALUE);

		 if ($evalue_deci<=$E_VALUE_deci && $clade_name eq "recipient") {
			 unless (exists $recip_hit_species{$species_name}) {
				 $recip_hit_species{$species_name}=1;
				 $recip_hit_species_no++;
				 if ($recip_hit_species_no<=100) {print OUT  "$key1\t$key2\t$similar\t$evalue\t$bitscore\t$slineage_name\t$species_name\t$clade_name\n";}
			 }
		 }

		 elsif ($evalue_deci<=$E_VALUE_deci && $clade_name eq "other") {
			 unless (exists $other_hit_species{$species_name}) {
				 $other_hit_species{$species_name}=1;
				 $other_hit_species_no++;
				 if ($other_hit_species_no<=100) { print OUT  "$key1\t$key2\t$similar\t$evalue\t$bitscore\t$slineage_name\t$species_name\t$clade_name\n";}
			 }
		 }

		 elsif ($evalue_deci<=$E_VALUE_deci && $clade_name eq "group") {
			 unless (exists $group_hit_species{$species_name}) {
				 $group_hit_species{$species_name}=1;
				 $group_hit_species_no++;
				 if ($group_hit_species_no<=100) { print OUT  "$key1\t$key2\t$similar\t$evalue\t$bitscore\t$slineage_name\t$species_name\t$clade_name\n";}			
			 }
		 }
	 }
    close OUT;


}

##################################################
### Calculate alien index for stored sequences ###
##################################################
print "Calculating alien index for $BLASTFILE ...\n";
foreach my $qseqid (keys %scores){

    my @lins = 'na';
	my $lins_no= 'na';
	my @unique_lins = 'na';
	my $unique_lins_no = 'na';
	if (exists $lineages{$qseqid}){
	 @lins=@{$lineages{$qseqid}};
	 $lins_no=scalar @lins;
	 @unique_lins = uniq(@lins);
	 $unique_lins_no=scalar @unique_lins;
	}
    my $recipient_no=0;
    my $group_no=0;
	my $other_no=0;
	my $pct_other_no='na';
	if (exists $clades{$qseqid}){
	   foreach my $clade_id (@{$clades{$qseqid}}) {
		 if ($clade_id eq 'recipient') {$recipient_no=$recipient_no+1;}
		 elsif ($clade_id eq 'group') {$group_no=$group_no+1;}
		 elsif ($clade_id eq 'other') {$other_no=$other_no+1;}
	  }
	}
    my $all_no=$recipient_no+$group_no+$other_no;
	if ($all_no>0) {
    $pct_other_no=$other_no*100/$all_no;
	}




	my $tophit_evalue = 1;
	my $tophit_bitscore = 0;
	my $tophit_id = 'na';
	my $tophit_name = 'na';
	my $tophit_tax = 'na';
	my $tophit_lin = 'na';	

	my $group_evalue = 1;
	my $group_bitscore = 0;
	my $group_id = 'na';
	my $group_name = 'na';
	my $group_tax = 'na';

	my $other_evalue = 1;
	my $other_bitscore = 0;
	my $other_id = 'na';
	my $other_name = 'na';
	my $other_tax = 'na';
	if (exists $tophits{$qseqid}){
		$tophit_evalue = $tophits{$qseqid}[0];
		$tophit_id = $tophits{$qseqid}[1];
		$tophit_name = $tophits{$qseqid}[2];
		$tophit_tax = $tophits{$qseqid}[3];
		$tophit_lin = $tophits{$qseqid}[4];
		$tophit_bitscore = $tophits{$qseqid}[5];
	}

	if (exists $scores{$qseqid}{'other'}){
		$other_evalue = $scores{$qseqid}{'other'}[0];
		$other_id = $scores{$qseqid}{'other'}[1];
		$other_name = $scores{$qseqid}{'other'}[2];
		$other_tax = $scores{$qseqid}{'other'}[3];
		$other_bitscore = $scores{$qseqid}{'other'}[4];	}

	if (exists $scores{$qseqid}{'group'}){
		$group_evalue = $scores{$qseqid}{'group'}[0];
		$group_id = $scores{$qseqid}{'group'}[1];
		$group_name = $scores{$qseqid}{'group'}[2];
		$group_tax = $scores{$qseqid}{'group'}[3];
		$group_bitscore = $scores{$qseqid}{'group'}[4];
	}
	my $ai = 'na';
	my $max_bitscore = 'na';
   if (exists $query_bitscores{$qseqid} && $query_bitscores{$qseqid}>0) {
	$max_bitscore = $query_bitscores{$qseqid};
	$ai = ( $other_bitscore / $max_bitscore ) - ( $group_bitscore / $max_bitscore );
   }else {print "$qseqid is bad query gene\n";}
	my $qlength = $query_lens{$qseqid};
	if ($ai <=1) {
	print OFIL "$qseqid\t$qlength\t$max_bitscore\t$tophit_id\t$tophit_tax\t$tophit_lin\t$tophit_evalue\t$tophit_bitscore\t$group_id\t$group_tax\t$group_evalue\t$group_bitscore\t$other_id\t$other_tax\t$other_evalue\t$other_bitscore\t$ai\t@unique_lins\t$unique_lins_no\t$recipient_no\t$group_no\t$other_no\t$pct_other_no\n";
   }
}


close OFIL;



sub taxdump {
	my $taxid = shift;
	my $root;
	my @taxlinarr;
	while (!defined $root){
		if (!defined $taxid){
			last;
		}
		if ($taxid == 1){$root = 1};
		push @taxlinarr, $taxid;
		$taxid = $taxnodes{$taxid};
	}
	my $taxlin = join(",", @taxlinarr);
	return $taxlin;
}

sub parse_config {
	my $file = shift;
	my %answer;

	open CONFIG, "$file" or die "Couldn't read config file $file: $!\n";
	while (<CONFIG>) {
		next if (/^#|^\s*$/);  # skip blanks and comments
		my ($variable, $value) = split /=/;
		$variable =~ s/\s+//g;
		$value =~ s/\s+//g;

		$answer{$variable} = $value;
	}
	close CONFIG;

	return %answer;
}

my $end_time=time;
my $duration=($end_time-$start_time)/60;
print"Execution time: $duration min.\n";