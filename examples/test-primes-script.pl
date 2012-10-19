#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec::Functions;
use FindBin;
use Time::HiRes qw(gettimeofday tv_interval);
use bigint;
use Data::BitStream::XS;
$|++; #flush the output buffer after every write() or print() function

# Should use a "put_test_string" sort of call on the test data, so we
# wouldn't need this.
my @test_data = (
  [  668, "Mersenne",    "--mersenne",   10**100],
  [ 7529, "Triplet",     "--triplet",    0],
  [ 7530, "Quadruplet",  "--quadruplet", 0],
  [23200, "Cousin",      "--cousin",     0],
  [23201, "Sexy",        "--sexy",       0],
  [ 1359, "Twin",        "--twin",       0],
  [ 5385, "Safe",        "--safe",       0],
  [ 5384, "SG",          "--sophie",     0],
  [ 2385, "Palindromic", "--palin",      32_965_656_923],
  [ 5479, "Lucas",       "--lucas",      0],
  [ 5478, "Fibonacci",   "--fibonacci",  0],
  [63980, "Pillai",      "--pillai",     2000],
  [28388, "Good",        "--good",       20000],
  [ 2407, "Cuban y+1",   "--cuban1",     0],
  [ 2648, "Cuban y+2",   "--cuban2",     100_000_000],
  [ 5234, "Primorial+1", "--pnp1",       2500],
  [ 6794, "Primorial-1", "--pnm1",       2500],
  [18239, "Euclid",      "--euclid",     0],
);
my %oeis_name = map { $_->[0] => $_->[1] } @test_data;
my %oeis_number = map { my $n=$_->[2]; $n=~s/^--//; $n => $_->[0] } @test_data;

# Verify additional filters
my @additional_filters;
foreach my $name (@ARGV) {
  $name =~ s/^--//;
  die "Unknown filter: $name\n" unless defined $oeis_number{$name};
  push @additional_filters, $name;
}

my $test_data_hash = read_script_data('script-test-data.bs');

if (@additional_filters > 0) {
  print "Additional Filters: ", join(" ", @additional_filters), "\n";
}
foreach my $test (@test_data) {
  my $oeis_no = $test->[0];
  if (!defined $test_data_hash->{$oeis_no}) {
    die "No test data found for OEIS $oeis_no : $test->[1] primes\n";
  }
  test_oeis(@$test, $test_data_hash);
}




sub read_script_data {
  my ($filename) = @_;

  die "Can't find test file: $filename\nRun make-script-test-data.pl\n"
      unless -r $filename;

  my $stream = Data::BitStream::XS->new( file => $filename, mode => 'ro' );
  my %hash;

  while (!$stream->exhausted) {
    my ($oeis_no, $is_bigint, $num_entries, @ref) = $stream->get_gamma(5);
    printf "%12s primes (OEIS A%06d): reading %7d entries..", $oeis_name{$oeis_no}, $oeis_no, $num_entries;
    if ($is_bigint) {
      print ",";
      my $k = 2;
      my @deltas = $stream->get_arice($k, $num_entries-2);
      print ".";
      # Check to see if we have any giant deltas
      foreach my $d (@deltas) {
        if ($d >= 18446744073709551614) {
          my $len = $stream->get_gamma;
          my $binstr = $stream->read_string($len);
          $d = Math::BigInt->new('0b' . $binstr);
        }
      }
      print ".";
      my $prev = $ref[1];
      push @ref, map { $prev = $_*2+$prev+2; } @deltas;
      $hash{$oeis_no} = \@ref;
      print ".\n";
    } else {
      no bigint;
      print ".";
      my $k = 2;
      my @deltas = $stream->get_arice($k, $num_entries-2);
      print ".";
      my $prev = $ref[1];
      push @ref, map { $prev = $_*2+$prev+2; } @deltas;
      $hash{$oeis_no} = \@ref;
      print ".\n";
    }
  }
  \%hash;
}

sub test_oeis {
  my($oeis_no, $name, $script_arg, $restrict, $test_data_hash) = @_;

  my @ref = @{ $test_data_hash->{$oeis_no} };
  my $end = $ref[-1];

  foreach my $filter_name (@additional_filters) {
    my $filter_no = $oeis_number{$filter_name};
    my %filter_data;
    undef @filter_data{ @{$test_data_hash->{$filter_no}} };
    my $filter_end = $test_data_hash->{$filter_no}->[-1];
    @ref = grep { exists $filter_data{$_} } @ref;
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

  die "Not equal numbers:  ", scalar @ref, " - ", scalar @scr, "\n"
    unless @ref == @scr;

  foreach my $i (0 .. $#ref) {
    die "$name prime $i not equal:  $ref[$i] - $scr[$i]\n"
        if $ref[$i] != $scr[$i];
  }
}
