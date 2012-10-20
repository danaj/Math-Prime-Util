#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec::Functions;
use FindBin;
use Time::HiRes qw(gettimeofday tv_interval);
use bigint;
use Data::BitStream::XS;
$|++; #flush the output buffer after every write() or print() function

# Maps between oeis name and number, filled in as we read sequences.
my %oeis_number; # short-name -> no
my %oeis_data;   # no -> ref to info+data

# returned array contains elements of:
#   [$oeis_no, $name, $script_arg, $num_entries, \@ref_data];
my $test_data = read_script_data('script-test-data.bs');

# Verify additional filters
my @additional_filters;
foreach my $name (@ARGV) {
  $name =~ s/^--//;
  my $oeis_no = $oeis_number{$name};
  die "Unknown filter: $name\n" unless defined $oeis_no;
  push @additional_filters, $oeis_no;
}
if (@additional_filters > 0) {
  print "Additional Filters: ",
        join(" ", map { $oeis_data{$_}->[2] } @additional_filters), "\n";
}

foreach my $test (@$test_data) {
  test_oeis(@$test);
}




sub read_script_data {
  my ($filename) = @_;

  die "Can't find test file: $filename\nRun make-script-test-data.pl\n"
      unless -r $filename;

  my $stream = Data::BitStream::XS->new( file => $filename, mode => 'ro' );
  my @data;

  while (!$stream->exhausted) {
    my $script_arg = get_text_string($stream);
    my $name       = get_text_string($stream);
    my ($oeis_no, $is_bigint, $num_entries, @ref) = $stream->get_gamma(5);
    printf "%12s primes (OEIS A%06d): reading %7d entries..", $name, $oeis_no, $num_entries;
    if ($is_bigint) {
      print ",";
      my $k = 2;
      my @deltas = $stream->get_arice($k, $num_entries-2);
      print ".";
      # Check to see if we have any giant deltas
      foreach my $d (@deltas) {
        if ( $d >= '18446744073709551614' ) {
          my $len = $stream->get_gamma;
          my $binstr = $stream->read_string($len);
          $d = Math::BigInt->new('0b' . $binstr);
        }
      }
      print ".";
      my $prev = $ref[1];
      push @ref, map { $prev = $_*2+$prev+2; } @deltas;
      print ".\n";
    } else {
      no bigint;
      print ".";
      my $k = 2;
      my @deltas = $stream->get_arice($k, $num_entries-2);
      print ".";
      my $prev = $ref[1];
      push @ref, map { $prev = $_*2+$prev+2; } @deltas;
      print ".\n";
    }
    my $row = [$oeis_no, $name, $script_arg, $num_entries, \@ref];
    push @data, $row;
    $oeis_data{$oeis_no} = $row;
    $oeis_number{$script_arg} = $oeis_no;
  }
  \@data;
}

sub test_oeis {
  my($oeis_no, $name, $script_arg, $num_entries, $ref_data) = @_;

  my @ref = @$ref_data;
  my $end = $ref[-1];
  $script_arg = '--' . $script_arg;

  foreach my $filter_no (@additional_filters) {
    #my $row = [$oeis_no, $name, $script_arg, $num_entries, \@ref];
    my $filter_name = $oeis_data{$filter_no}->[2];
    my $filter_data_ref = $oeis_data{$filter_no}->[4];
    my %filter_data_hash;
    undef @filter_data_hash{ @$filter_data_ref };
    my $filter_end = $filter_data_ref->[-1];
    @ref = grep { exists $filter_data_hash{$_} } @ref;
    $script_arg .= " --$filter_name";
    $end = $filter_end if $end > $filter_end;  # bring endpoint down
  }

  printf "%12s primes (OEIS A%06d): generating..", $name, $oeis_no;

  my $start = [gettimeofday];
  my @scr = split /\s+/, qx+$FindBin::Bin/../bin/primes.pl $script_arg 1 $end+;
  {
    no bigint;
    my $num_generated = scalar @scr || 0.1;
    my $seconds = tv_interval($start);
    my $msperprime = ($seconds * 1000.0) / $num_generated;
    printf " %7d. %7.2f ms/prime\n", $num_generated, $msperprime;
  }

  if (scalar @ref != scalar @scr) {
    warn "  $FindBin::Bin/../bin/primes.pl $script_arg 1 $end\n";
    die "Not equal numbers:  ", scalar @ref, " - ", scalar @scr, "\n";
  }

  foreach my $i (0 .. $#ref) {
    die "$name prime $i not equal:  $ref[$i] - $scr[$i]\n"
        if $ref[$i] != $scr[$i];
  }
}

sub put_text_string {
  my ($stream, $str) = @_;
  $stream->put_gamma(ord($_)) for (split "", $str);
  $stream->put_gamma(0);
  1;
}

sub get_text_string {
  my ($stream) = @_;
  my $str = '';
  while (my $c = $stream->get_gamma) {
    $str .= chr($c);
  }
  $str;
}
