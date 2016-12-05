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
#notes: filter column contains <NA> values - ignore 
#
#
##################################################################

usage() unless(@ARGV > 1);
my ($vcf_file, $tech_file, $outfile_prefix) = ($ARGV[0], $ARGV[1], $ARGV[2]);

my $techniques = map_technologies($tech_file);

my @headers = ();

my $data = load_data($vcf_file);
my $statistics = create_stats();

my %suffix;

$suffix{'breakpoint_variation'} = '.breakpoint_variation.tsv';
$suffix{'breakpoint_variation_id_counts'} = '.breakpoint_variation.id_counts.tsv';
breakpoint_variation( $outfile_prefix.$suffix{'breakpoint_variation'},  $outfile_prefix.$suffix{'breakpoint_variation_id_counts'} );

$suffix{'total_callset_count'} = '.total_callset_count.tsv';
total_callset_count( $outfile_prefix.$suffix{'total_callset_count'});

$suffix{'total_callset_count_by_type'} = '.total_callset_count_by_type.tsv';
total_callset_count_by_type( $outfile_prefix.$suffix{'total_callset_count_by_type'});
# callsets_per_event();

##########################################################################################
# 											statistics									 #
##########################################################################################
#------------------------------------------
# histogram of number of callers supporting a variant
#------------------------------------------
sub callsets_per_event {
	my %hash = ();
	
	while(my($k,$v) = each(%$statistics)) {	#k = caller, v = array reference of records
		foreach my $record (@$v) {
			my $types = $record->{"types"};
			my $id = $record->{"id"};			
			insert($id, $k, \%hash);
		}
	}
	
	while(my($k,$v) = each(%hash)) {
		my @uniq = uniq(@$v);
		my $str = join("\t", $k,scalar(@uniq));
		print "$str\n";
	}
}

#------------------------------------------
# total number of variants in each callset
#------------------------------------------
sub total_callset_count {
  my($outfile)=@_;
  open(OUT,">$outfile");
  my @header = qw(CALLER N_VARIANTS);
  print OUT join "\t",@header; print OUT "\n";
	while(my($k,$v) = each(%$statistics)) {	#k = caller, v = array reference of records
		#print "$k:", scalar (@$v),"\n";
		my $count = 0;
	
		foreach my $record (@$v) {
			my $types = $record->{"types"};

			while(my($k2,$v2) = each(%$types)) {
				$count += scalar(@$v2);
			}
		}
		print OUT "$k\t$count\n";
	}
  close(OUT);
}

#------------------------------------------
# total number of variants in each callset, split by variant type
#------------------------------------------
sub total_callset_count_by_type {
  my($outfile)=@_;
  open(OUT,">$outfile");
  # create a header line for the output file
  my @header = qw(CALLER TYPE N_VARIANTS);
  print OUT join "\t",@header; print OUT "\n";
	while(my($k,$v) = each(%$statistics)) {	#k = caller, v = array reference of records
		#print "$k:", scalar (@$v),"\n";
		my %count;
	
		foreach my $record (@$v) {
			my $types = $record->{"types"};

			while(my($k2,$v2) = each(%$types)) {
				$count{$k2} += scalar(@$v2);
			}
		}
                foreach my $t (keys %count) {
		  print OUT "$k\t$t\t$count{$t}\n";
                }
	}
  close(OUT);
}


#------------------------------------------
# Make a new data structure with only variants supported by at least 4 callers
# compute median start and end if multiple comma-delim segments 
#------------------------------------------
sub breakpoint_variation {
  my($outfile,$count_outfile)=@_;
  my $min_callers = 4;
  my %num_callers;
  my %data;
  my %fail;
  my %id_tracker;
  my %failed_ids;

  open(OUT,">$outfile");
  open(COUT,">$count_outfile");
  while(my($caller,$rec) = each(%$statistics)) {	#key = caller, value = array reference of records
	#print "$caller:", scalar (@$rec),"\n";
	#my $ncallers = 0;
        
	
	foreach my $record (@$rec) {
		my $types = $record->{"types"}; #ref to hash
		my $id = $record->{"id"};
		#print "id $id num types ", scalar keys %$types," \n";

               # tracket to count the total number of input ids
               # keep as a hash since each id will be entered  multiple times since looping over caller
               $id_tracker{$id} = 1;

		# skip records with more than one type in this caller 
                if (scalar keys %$types >1) { 
                   $fail{MULTIPLE_TYPES_PER_CALLER}++;
                   #print "multiple types in one caller\n";
                   # track ids with at least one fail
                   $failed_ids{$id} = 1; 
                   next; 
                 }

                # count the number of callers that called this variant with this type
		while(my($type,$loc_aref) = each(%$types)) {
                   #print "type $type loc_aref $loc_aref\n";
		   $num_callers{$id}{$type}++;
   		   #print "saving id $id type $type caller $caller\n";
		}
	} #end for record
   } #end while

   # skip records with too few callers or with multiple types
   # for the kept variants, calculate the median start and stop
   # for each caller, calculate the abs distance for each start and stop to the median.
   foreach my $id (keys %num_callers) {

     # skip variants in which one caller had multiple types
     if (defined $failed_ids{$id}) { next; }

     # skip variants with more than one type
     my $num_types = scalar keys %{$num_callers{$id}};
     #print "id $id numtypes $num_types\n";
     if ($num_types >1) { 
       $fail{MULTIPLE_TYPES}++;
       next; 
     }

     # get the name of the single type for this variant
     my @types = keys %{$num_callers{$id}};
     my $single_type = $types[0];

     # check for sufficient caller support (for the variant with one type)
     # print "checking type $single_type\n";
     if ($num_callers{$id}{$single_type} < $min_callers) { 
         #print "$num_callers{$id}{$single_type} is too little support for id $id single_type $single_type\n";
         $fail{LOW_CALL_SUPPORT}++;
         next; 
      }
     #print "id $id numtypes $num_types singletype $single_type numcall $num_callers{$id}{$single_type}\n";

     # make a new data structure to save the relevant info
     $data{$id}{TYPE} = $single_type;
     $data{$id}{NUM_CALLERS} = $num_callers{$id}{$single_type};

     # For the passing variants, compute the median start and stop
     # If there is more than one start and stop for one caller, use the min start and max stop

     # retrieve the original start and stop fields from the original
     # oops this version is not indexed by id. Need to loop over the statistics records 
     # and check for the id. see next loop.
   }

   # Loop over the original list of records and extract starts and stops for the 
   # variants that pass the filters. 

   while(my($caller,$rec) = each(%$statistics)) {	#key = caller, value = array reference of records
     foreach my $record (@$rec) {
	my $id = $record->{"id"};
        # skip this record if the id is not in the set of ids to keep
        if (!defined $data{$id}) { next; }
        # sanity check: confirm that there is only one type
	my $types = $record->{"types"}; #ref to hash
        if (scalar keys %$types >1) { print STDERR "too many types for id $id\n"; }
	# loop over the call sets that have this variant.
        # (there should be at least $min_callers)
        while(my($type,$locations_aref) = each(%$types)) {
	   # another sanity check: is type the same as single_type
           if ($data{$id}{TYPE} ne $type) { print STDERR "type mismatch $data{$id}{TYPE} vs $type for id $id\n"; }
	   my ($min_start,$max_end) = ($locations_aref->[0]->{"start"}, $locations_aref->[0]->{"stop"} );
           for(my $i = 1; $i < scalar (@$locations_aref); $i++) {
		my $ref = $locations_aref->[$i];
		$min_start = $ref->{"start"} if($ref->{"start"} < $min_start);
		$max_end = $ref->{"stop"} if($ref->{"stop"} > $max_end);
	   }
           #print "id $id min $min_start max $max_end caller $caller\n";
           # save the start and end values for each caller
           push @{$data{$id}{start}}, $min_start;
           push @{$data{$id}{end}}, $max_end;
           push @{$data{$id}{caller}}, $caller;
           
         
        } #end while type/locations

     } #end foreach record
   } #end while caller,rec


# collected all the data. Now compute the medians for start and stop for each variant
my @header = qw(ID TYPE NUM_SUPPORTING_CALLERS CALLSET START_DIST END_DIST);
print OUT join "\t",@header; print OUT "\n";
  foreach my $id (keys %data) {
       # compute the median start and stop across the callers
       $data{$id}{MEDIAN_START} = &median(\@{$data{$id}{start}});
       $data{$id}{MEDIAN_END} = &median(\@{$data{$id}{end}});
       # print "id $id medstart $data{$id}{MEDIAN_START} medend $data{$id}{MEDIAN_END}\n";

       # loop over all the technologies, and print the abs dist of start/stop
       # from respective median
       for (my $j=0; $j<scalar @{$data{$id}{caller}}; $j++) {
          my $start_dist = abs($data{$id}{start}[$j] - $data{$id}{MEDIAN_START});
          my $end_dist = abs($data{$id}{end}[$j] - $data{$id}{MEDIAN_END});
          #print "id $id start $data{$id}{start}[$j] end $data{$id}{end}[$j]\n";
          print OUT join "\t",($id,$data{$id}{TYPE}, $data{$id}{NUM_CALLERS}, $data{$id}{caller}[$j],$start_dist,$end_dist); print OUT "\n";
       }
  }

# REPORTING
# print the total number of input ids
print COUT "INPUT_ID_COUNT\t", scalar keys %id_tracker, "\n";
print COUT "OUTPUT_ID_COUNT\t", scalar keys %data,"\n";
# report the fail reasons
  foreach my $fkey (keys %fail) {
    print COUT join "\t", ($fkey, $fail{$fkey}); print COUT "\n";
  }
close(OUT); 
close(COUT); 
}

#---------------------------------------------
# compute median. pass array by reference
# from http://www.perlmonks.org/index.pl?node_id=474564
#---------------------------------------------
sub median { $_[0]->[ @{$_[0]} / 2 ] }


##########################################################################################
# 											functions									 #
##########################################################################################
sub create_stats {
	my %stats = ();

	foreach my $line (@$data) {
		$line = trim($line);
		if($line =~ /^#chrom/i) {
			@headers = split(/\t+/, $line);
			next;
		} elsif($line =~ /^#/) {
			next;
		}
	
		my @tmp = split(/\s+/,$line);
		my ($chr, $start) = ($tmp[0], $tmp[1]);
		my $stop = get_end_location($tmp[7]);
		my $id = $tmp[2];
	
		for(my $i = 9; $i < scalar(@headers); $i++) {
			next if($tmp[$i] =~ /NaN/i);
		
			my $technique = $techniques->{$headers[$i]};

			my $record = create_record($id, $technique, $chr, $start, $stop, $tmp[$i]);
			$record->{"technology"} = $technique;
			insert($headers[$i], $record, \%stats);
		}
	}
	return \%stats;
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
	return $record;
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
#OLDsub parse_coordinates {
#OLD	my $msg = join(" ", "Error in function",(caller(0))[3],"\nCause: not enough arguments, line: ".__LINE__);
#OLD	exit_with_msg($msg) if(@_ < 1);
#OLD	my ($str) = ($_[0]); 
#OLD
#OLD	my @coordinates = ();
#OLD	my @tmp = split(/,/,$str);
#OLD	foreach my $pos (@tmp) {
#OLD		my ($start,$stop) = split(/-/,$pos);
#OLD		$start = substr($start,2);
#OLD		$stop = substr($stop,2);
#OLD		
#OLD		push(@coordinates, {"start" => $start, "stop" => $stop});
#OLD	}
#OLD	return \@coordinates;
#OLD}

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
