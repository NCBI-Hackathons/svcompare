#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp;
use Set::Scalar;
use File::Basename;
use Sort::Naturally;
use Bio::SeqIO;
use List::Util qw(max min);
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

##################################################################
# Author: Andi Dhroso
# Date: 08-17-2016
# Version: 0.1
#
# Description: 
# Generate statistics for different callsets.
#
# notes: any regions not in the actual chromosome is ignored
##################################################################

usage() unless(@ARGV > 1);
my ($vcf_file, $tech_file) = ($ARGV[0], $ARGV[1]);

my $techniques = map_technologies($tech_file);

my @headers = ();

my $data = load_data($vcf_file);
my ($statistics, $types) = create_stats();

# total_callset_count();
# callsets_per_event();
# callsets_per_event_type_aware();

# coorcordance_break_points_3b_i();

	
	
##########################################################################################
# 											statistics									 #
##########################################################################################
sub callsets_per_event {
	my %hash = ();
	
	while(my($k,$v) = each(%$statistics)) {	#k = caller, v = array reference of records
		foreach my $record (@$v) {
			my $id = $record->{"id"};			
			insert($id, $k, \%hash);
		}
	}
	
	print "#Id,CallsetsPerEvent\n"; #header
	
	while(my($k,$v) = each(%hash)) {
		my @uniq = uniq(@$v);
		my $str = join("\t", $k,scalar(@uniq));
		print "$str\n";
	}
}
sub callsets_per_event_type_aware {
	my %stats = ();
	while(my($k,$v) = each(%$statistics)) {	#k = caller, v = array reference of records
		foreach my $record (@$v) {
			my $id = $record->{"id"};			
			
			my $types = $record->{"types"};
			while(my($k2,$v2) = each(%$types)) {
				my $count = scalar(uniq(@$v2));
				
				if(exists $stats{$id}) {
					my $ref = $stats{$id};
					
					if(exists $ref->{$k2}) {
						$ref->{$k2} += $count;
					} else {
						$ref->{$k2} = $count;
					}
			
				} else {
					$stats{$id} = {$k2 => $count};
				}
			}
		}
	}
	
	my @tmp = ();
	push(@tmp, $_) foreach(@$types);
	my $str = join(",", @tmp);
	print "#$str\n";	#header
	$str=""; #just 
	
	while(my($k,$v) = each(%stats)) {		
		@tmp = ();
		
		foreach my $t (@$types)	{
			if(exists $v->{$t}) {
				push(@tmp, $v->{$t});
			} else {
				push(@tmp, 0);			
			}
			$str = join(",",@tmp);
		}
		print "$str\n";
	}
}

sub coorcordance_break_points_3b_i {
	my $data = convert_to_row_base();

	print "#ID,StartDifference,EndDifference\n";
	while(my($id,$caller) = each(%$data)) {
		my @start = (); 
		my @end = ();	
		
		while(my($caller,$record) = each(%$caller)) {
			my $types = $record->{"types"};

			if(scalar(keys %$types) < 1) {
				my ($minstart,$maxstart, $minend,$maxend) = get_range($types);
				push(@start, ($minstart,$maxstart));
				push(@end, ($minend,$maxend));
			} else {				
				my @keys = keys %$types;
				my ($minstart,$maxstart) = get_min_max_start($types->{$keys[0]});
				push(@start, ($minstart,$maxstart));

				my ($minend,$maxend) = get_min_max_end($types->{$keys[0]});
				push(@end, ($minend,$maxend));
			}
		}
		next unless(scalar(@start) > 0 && scalar(@end) > 0);
		if(scalar(@start) > 1 && scalar(@end) > 1) {
			my ($minstart,$maxstart,$minend,$maxend) = (min(@start), max(@start), min(@end), max(@end));
			print "$id,",$maxstart-$minstart, ",",$maxend-$minend,"\n";
		} else {
			print "$id,",$start[0], ",",$end[0],"\n";
		}
	}
}
sub get_range {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($types) = ($_[0]); 
	
	my @start = (); my @end = ();
	while(my($type, $arrayref) = each(%$types)) {
		foreach my $loc (@$arrayref) {
			push(@start, $loc->{"start"});
			push(@end, $loc->{"stop"});
		}
	}
	my ($minstart,$maxstart,$minend,$maxend) = (min(@start), max(@start),min(@end), max(@end));

	return ($minstart,$maxstart,$minend,$maxend);
}
sub convert_to_row_base {
	my %data = ();
	while(my($caller,$records) = each(%$statistics)) {	#k = caller, v = array reference of records

		foreach my $record (@$records) {
			my $id = $record->{"id"};			
			my $types = $record->{"types"};
			
			if($id eq "DEL0078560SUR") {
				#print $id,",",scalar(keys %$types),"\n";
			}
			
			if(exists $data{$id}) {
				my $ref = $data{$id};
				
				my $tmp = {
					"caller"	=> $caller,
					"types"		=> $types,
					"chr"		=> $record->{"start"},
					"start"		=> $record->{"stop"},
					"stop"		=> $record->{"length"},
					"technology" => $record->{"technology"}
				};
				$ref->{$caller} = $tmp;	
								
			} else {
				my $ref = {
					"caller"	=> $caller,
					"types"		=> $types,
					"chr"		=> $record->{"start"},
					"start"		=> $record->{"stop"},
					"stop"		=> $record->{"length"},
					"technology" => $record->{"technology"}
					};
				$data{$id} = { $caller => $ref};
			}
		}
	}
	return \%data;
}

sub get_min_max_start {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($array) = ($_[0]); 

	my @tmp = ();
	foreach my $ref (@$array) {
		push(@tmp, $ref->{"start"});
	}

	my $min = min(@tmp);
	my $max = max(@tmp);

	return ($min,$max);
}
sub get_min_max_end {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($array) = ($_[0]); 

	my @tmp = ();
	foreach my $ref (@$array) {
		push(@tmp, $ref->{"stop"});
	}

	my $min = min(@tmp);
	my $max = max(@tmp);

	return ($min,$max);
}

sub total_callset_count {
	while(my($k,$v) = each(%$statistics)) {	#k = caller, v = array reference of records
		#print "$k:", scalar (@$v),"\n";
		my $count = 0;
	
		foreach my $record (@$v) {
			my $types = $record->{"types"};

			while(my($k2,$v2) = each(%$types)) {
				$count += scalar(@$v2);
			}
		}
		print "$k\t$count\n";
	}
}
##########################################################################################
# 											functions									 #
##########################################################################################
sub create_stats {
	my %stats = ();

	my @eventtypes = ();
	foreach my $line (@$data) {
		$line = trim($line);
		if($line =~ /^#chrom/i) {
			@headers = split(/\t+/, $line);
			next;
		} elsif($line =~ /^#/) {
			next;
		}
	
		my @tmp = split(/\s+/,$line);
		next if($tmp[4] eq "<NA>");
		
		my ($chr, $start) = ($tmp[0], $tmp[1]);
		next if(!looks_like_number($chr));
		
		my $stop = get_end_location($tmp[7]);
		my $id = $tmp[2];
	
		for(my $i = 9; $i < scalar(@headers); $i++) {
			next if($tmp[$i] =~ /NaN/i);
		
			my $technique = $techniques->{$headers[$i]};
			
			my ($record, $types) = create_record($id, $technique, $chr, $start, $stop, $tmp[$i]);
			$record->{"technology"} = $technique;
			insert($headers[$i], $record, \%stats);
			
			push(@eventtypes, @$types);
			@eventtypes = uniq(@eventtypes);
		}
	}
	return (\%stats, \@eventtypes);
}
sub create_record {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 6);
	my ($id, $technique, $chr, $start, $stop, $types) = ($_[0],$_[1],$_[2],$_[3],$_[4],$_[5]); 

	my $container = create_event($types);
	my $typesref = $container->{"types"};
	my $length = $container->{"length"};
	
	
	my $record = {
		"types" => $typesref,
		"id" => $id,
		"chr" => $chr,
		"start" => $start,
		"stop" => $stop,
		"length" => $length
	};
	
	my @types = keys(%$typesref);
	return ($record, \@types);
}

sub create_event {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($str) = ($_[0]); 

	my %hash = ();
	my ($gt,$ln,$dv,$type,$co) = split(/:/,$str);

	my $types = parse_type($type);
	
	my $coords = parse_coordinates($co);
	
	for(my $i = 0; $i < scalar(@$types); $i++) {
		insert($types->[$i], $coords->[$i],\%hash);
	}
	return {"types" => \%hash, "length" => $ln};	
}
sub parse_type {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($str) = ($_[0]); 

	my @types = split(/,/,$str);
	return \@types;
}
sub parse_coordinates {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($str) = ($_[0]); 
	
	my @coordinates = ();
	my @tmp = split(/,/,$str);
	foreach my $pos (@tmp) {
		my ($start,$stop) = split(/-/,$pos);
		my $index = index($start,"_");
		$start = substr($start,$index+1);
		
		$index = index($stop,"_");
		$stop = substr($stop ,$index+1);
		
		push(@coordinates, {"start" => $start, "stop" => $stop});
	}
	return \@coordinates;
}

sub get_end_location {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($str) = ($_[0]); 

	my @tokens = split(/;/,$str);
	my $pos = index($tokens[4],"=");
	my $k = substr($tokens[4],$pos+1);

	$pos = index($tokens[5],"=");
	my $v = substr($tokens[5],$pos+1);

	return ($k, $v);
}

sub map_technologies {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($file) = ($_[0]); 

	my $data = load_data($file);
	my %hash = ();
	foreach my $line (@$data) {
		$line = trim($line);
		
		my ($k,$v) = split(/\t+/,$line);
		$hash{$k} = $v;
	}
	return \%hash;
}
##########################################################################################
# 											common										 #
##########################################################################################
sub execute {
	my $EXIT = 1;
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	my ($cmd, $status) = ($_[0], $_[1]); 

	$msg = "no message";
	
	my $rc = system($cmd);
	if($rc) {
		$msg = join(" ", "error: command could not execute:", $cmd, "\nLine:", __LINE__,"\nQuiting program...");
		exit_with_msg($msg) if($status == $EXIT);
	}
	
	$msg = join(" ", "error: command could not execute:", $cmd, "\nLine:", __LINE__, "\n");
	warn($msg) if($rc);
	
	return $rc;
}

sub show_msg {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	
	my ($delimiter, $array) = ($_[0], $_[1]);
	my $str = join($delimiter, @$array);
	print "$str\n";
}
sub get_files {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	
	my ($path,$match) = ($_[0], $_[1]);
	$path = trim($path);
	my @files = glob("$path/$match");
	return \@files;
}
sub file_exists {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg("not enough arguments, line: ".__LINE__) if(@_ < 1);
	
	my $file = $_[0];
	unless( -e $file) {
		warn "warning, file does not exits:\n$file\n";
		exit;
	}
}
sub path_exists {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg("not enough arguments, line: ".__LINE__) if(@_ < 1);
	
	my $path = $_[0];
	unless( -d $path) {
		warn "warning, path does not exits:$path\n";
		return 0;
	}
	return 1;
}
sub load_data {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	
	my $file = $_[0];
	
	open(IN,"<", $file) or exit_with_msg("$!, $file line: ". __LINE__);
	my @data = ();
	@data = <IN>;
	close(IN);
	
	chomp(@data);
	return \@data;
}
sub load_sequences {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	my ($file,$format) = ($_[0],$_[1]);
	
	my @seqs = ();
	my $in  = Bio::SeqIO->new(-file => $file , -format => $format);
	while ( my $seq = $in->next_seq() ) {
    	push(@seqs, $seq);
    }
	return \@seqs;
}
sub write_sequences {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 3);
	my ($seqs,$file, $format) = ($_[0],$_[1],$_[2]);

    my $out = Bio::SeqIO->new(-file => ">$file" ,
                           -format => $format);

	foreach my $seq (@$seqs) {
        $out->write_seq($seq);
    }
}
sub trim {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	
	my $string = $_[0];
	
	#$string =~ s/^\s+|\s+$//g;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	$string =~ s/\/$//;
	
	return $string;
}
sub get_basename {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	
	my ($line,$match) = ($_[0], $_[1]);
	
	my($filename, $path, $suffix) = fileparse($line, ($match));
	$path = trim($path);
	
	return ($filename, $path, $suffix);
}
sub tolower {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($arrayref) = ($_[0]);	

	my ($counter) = (0);
	foreach my $line (@$arrayref) {
		$line = trim($line);
		$line = lc($line);
		$arrayref->[$counter] = $line;
		$counter++;
	}
}
sub touper {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($arrayref) = ($_[0]);	

	my ($counter) = (0);
	foreach my $line (@$arrayref) {
		$line = trim($line);
		$line = uc($line);
		$arrayref->[$counter] = $line;
		$counter++;
	}
}
sub get_hash {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($array) = ($_[0]);
	
	my %hash = map{$_ => 1}@$array;
	return \%hash;
}
sub unique {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 1);
	my ($data) = ($_[0]);

	my @unique = do { my %seen; grep { !$seen{$_}++ } @$data };
 	return \@unique;
 	
 	#another solution
 	# my %seen;
  	#return grep { !$seen{$_}++ } @_;
}
sub get_temp_filename {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg("not enough arguments, line: ".__LINE__) if(@_ < 1);
	
	my $path = trim($_[0]);
	
    my $fh = File::Temp->new(
        TEMPLATE => "tempXXXXX",
        DIR      => "$path",
        SUFFIX   => ".txt",
        UNLINK => 0
    );

    return $fh->filename;
}
sub insert {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 3);
	my ($g1, $g2, $network) = ($_[0],$_[1],$_[2]); 

	if(exists $network->{$g1}) {
		my $tmp = $network->{$g1};
		push(@$tmp, $g2);
		$network->{$g1} = $tmp;	
	} else {
		$network->{$g1} = [$g2];
	}
}
sub set_symmetric_difference {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	my ($a, $b) = ($_[0], $_[1]);
	
	my $s = Set::Scalar->new(@$a);
	my $t = Set::Scalar->new(@$b);
	
	my $d = $s->symmetric_difference($t);
	return [$d->members];
}
sub set_difference {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	my ($a, $b) = ($_[0], $_[1]);

	my $s = Set::Scalar->new(@$a);
	my $t = Set::Scalar->new(@$b);
	my $d = $s->difference($t);	#returns elements that are present in (s) only

	return [$d->members];
}

sub set_intersect {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	my ($a, $b) = ($_[0], $_[1]);
	
	my $s = Set::Scalar->new(@$a);
	my $t = Set::Scalar->new(@$b);
	
	my $i = $s->intersection($t);	#returns elements that are present in both (s) and (t)
	return [$i->members];
}

sub set_union {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	exit_with_msg($msg) if(@_ < 2);
	my ($a, $b) = ($_[0], $_[1]);
	
	my $s = Set::Scalar->new(@$a);
	my $t = Set::Scalar->new(@$b);
	
	my $u = $s->union($t);
	return [$u->members];
}
sub exit_with_msg {
	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
	$msg = $_[0] unless (@_ < 1);

	warn $msg,"\n";
	exit;
}
sub usage {
	exit_with_msg "usage: perl stat.pl < vcf file> <caller-technology map file>\n";
	exit;
}
exit;
