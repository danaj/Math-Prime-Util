#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec::Functions;
use FindBin;
use bigint;
use Data::BitStream::XS;
use Math::Prime::Util qw/is_prime/;
$|++;

# Encode all the OEIS text files for our primes.pl testing into a bitstream.
# This not only makes the test script run much faster, but it turns 18 text
# files of 5MB into one ~300k file.

my @test_data = (
  # OEIS#  TEXT NAME      script-arg   skip if > this
  [ 7529, "Triplet",     "triplet",    0],
  [ 7530, "Quadruplet",  "quadruplet", 0],
  [23200, "Cousin",      "cousin",     0],
  [23201, "Sexy",        "sexy",       0],
  [ 1359, "Twin",        "twin",       0],
  [ 5385, "Safe",        "safe",       0],
  [ 5384, "SG",          "sophie",     0],
  [68652, "Circular",    "circular",   0],
  [ 2407, "Cuban y+1",   "cuban1",     0],
  [ 2648, "Cuban y+2",   "cuban2",     0],
  [ 2385, "Palindromic", "palin",      32_965_656_923],
  [  668, "Mersenne",    "mersenne",   10**100],
  [ 5479, "Lucas",       "lucas",      0],
  [ 5478, "Fibonacci",   "fibonacci",  0],
  [63980, "Pillai",      "pillai",     2000],
  [28388, "Good",        "good",       20000],
  [31157, "Lucky",       "lucky",      0],
  [ 5234, "Primorial+1", "pnp1",       2500],
  [ 6794, "Primorial-1", "pnm1",       2500],
  [18239, "Euclid",      "euclid",     0],
);

foreach my $test (@test_data) {
  my $oeis_no = $test->[0];
  my $filename = sprintf("b%06d.txt", $oeis_no);
  my $link     = sprintf("http://oeis.org/A%06d/b%06d.txt", $oeis_no, $oeis_no);
  if (!-r $filename) {
    warn "Getting $filename from $link\n";
    qx/wget $link/;
    die "Could not retrieve.  Bailing\n" unless -r $filename;
  }
  my $ref_data = read_oeis(@$test);
  push @$test, $ref_data;
}

my $stream = Data::BitStream::XS->new( file => 'script-test-data.bs', mode => 'w' );
foreach my $test (@test_data) {
  encode_oeis(@$test);
}
$stream->write_close();

sub read_oeis {
  my($oeis_no, $name, $script_arg, $restrict) = @_;
  die "Restrict isn't defined for $oeis_no : $name" unless defined $restrict;

  my $filename = sprintf("b%06d.txt", $oeis_no);
  my $link     = sprintf("http://oeis.org/A%06d/b%06d.txt", $oeis_no, $oeis_no);

  my @ref;
  {
    open my $fh, '<', $filename
        or die "Can't read $filename.\nYou should run:\n  wget $link\n";
    printf "%12s primes: reading %12s...", $name, $filename;
    my $char = " ";
    while (<$fh>) {
      next unless /^(\d+)\s+(\d+)/;
      my $v = (length($2) < 20)  ?  $2  :  Math::BigInt->new("$2");
      if ($restrict > 0 && $v > $restrict) {
        $char = '*';
        last;
      }
      push @ref, $v;
    }
    close $fh;
    print "$char";
  }
  printf " %7d.", scalar @ref;
  print "  Testing..";
  if ($ref[-1] > 18446744073709551615) {
    print ",";
    # Check for monotonic and primeness
    foreach my $i (0 .. $#ref) {
      die "non-prime in $oeis_no $name\n" unless is_prime($ref[$i]);
      if ($i > 0) {
        die "non-monotonic sequence in $oeis_no $name ($i $ref[$i-1] $ref[$i])\n" if $ref[$i] <= $ref[$i-1];
        die "even number in $oeis_no $name\n" if ($ref[$i] % 2) == 0;
      }
    }
  } else {
    no bigint;
    print ".";
    # Check for monotonic and primeness
    foreach my $i (0 .. $#ref) {
      die "non-prime in $oeis_no $name\n" unless is_prime($ref[$i]);
      if ($i > 0) {
        die "non-monotonic sequence in $oeis_no $name\n" if $ref[$i] <= $ref[$i-1];
        die "even number in $oeis_no $name\n" if ($ref[$i] % 2) == 0;
      }
    }
  }
  print "done\n";
  return \@ref;
}

sub encode_oeis {
  my($oeis_no, $name, $script_arg, $restrict, $ref_data) = @_;

  my @ref = @$ref_data;
  printf "%12s primes: stream..", $name;

  put_text_string($stream, $script_arg);
  put_text_string($stream, $name);

  if ($ref[-1] > 18446744073709551615) {
    print ",";
    # Store the first two values, then a list of deltas
    $stream->put_gamma($oeis_no, 1, scalar @ref, $ref[0], $ref[1]);
    print ".";
    my @deltas = map { ($ref[$_] - $ref[$_-1] - 2)/2 } (2..$#ref);
    print ".";
    # Ugly...  Check for anything really big;
    my @giant;
    foreach my $d (@deltas) {
      if ($d >= 18446744073709551614) {
        push @giant, $d;
        $d = 18446744073709551614;
      }
    }
    print ".";
    my $k = 2;
    $stream->put_arice($k, @deltas);
    print ".";
    # Store giant deltas raw
    foreach my $d (@giant) {
      if (ref($d) ne 'Math::BigInt') {
        warn "big delta $d isn't a bigint.\n";
        $d = Math::BigInt->new(0);
      }
      my $binstr = substr($d->as_bin, 2);
      $stream->put_gamma(length($binstr));
      $stream->put_string($binstr);
    }
  } else {
    no bigint;
    print ".";
    # Store the first two values, then a list of deltas
    $stream->put_gamma($oeis_no, 0, scalar @ref, $ref[0], $ref[1]);
    print ".";
    my @deltas = map { ($ref[$_] - $ref[$_-1] - 2)/2 } (2..$#ref);
    print ".";
    my $k = 2;
    $stream->put_arice($k, @deltas);
  }

  print "done\n";
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
