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

my $test_data_hash = read_script_data('script-test-data.bs');

foreach my $test (@test_data) {
  my $oeis_no = $test->[0];
  my $test_data = $test_data_hash->{$oeis_no};
  if (!defined $test_data) {
    die "No test data found for OEIS $oeis_no : $test->[1] primes\n";
  }
  test_oeis(@$test, $test_data);
}




sub read_script_data {
  my ($filename) = @_;

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
  my($oeis_no, $name, $script_arg, $restrict, $ref_data) = @_;

  my @ref = @$ref_data;
  printf "%12s primes (OEIS A%06d): generating..", $name, $oeis_no;

  #print "\n";
  #print "reference data:\n";
  #print "  $_\n" for @ref;
  #print "primes.pl $script_arg 1 $ref[-1]\n";
  my $start = [gettimeofday];
  my @scr = split /\s+/, qx+$FindBin::Bin/../bin/primes.pl $script_arg 1 $ref[-1]+;
  {
    no bigint;
    my $num_generated = scalar @scr;
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
