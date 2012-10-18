#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec::Functions;
use FindBin;
$|++; #flush the output buffer after every write() or print() function

test_oeis(668, "Mersenne", "--mersenne", '1' . '0' x 100);
#test_oeis(2407, "Cuban y+1", "--cuban1");
#test_oeis(2648, "Cuban y+2", "--cuban2", 100_000_000);
test_oeis(5234, "Primorial+1", "--pnp1", 2500);
test_oeis(6794, "Primorial-1", "--pnm1", 2500);
test_oeis(18239, "Euclid", "--euclid");
test_oeis(7529, "Triplet", "--triplet");
test_oeis(7530, "Quadruplet", "--quadruplet");
test_oeis(23200, "Cousin", "--cousin");
test_oeis(23201, "Sexy", "--sexy");
test_oeis(1359, "Twin", "--twin");
test_oeis(5385, "Safe", "--safe");
test_oeis(5384, "SG", "--sophie");
test_oeis(2385, "Palindromic", "--palin");
test_oeis(5479, "Lucas", "--lucas");
test_oeis(5478, "Fibonacci", "--fibonacci");
test_oeis(63980, "Pillai", "--pillai", 2000);
test_oeis(28388, "Good", "--good", 20000);

sub test_oeis {
  my($oeis_no, $name, $script_arg, $restrict) = @_;

  my $filename = sprintf("b%06d.txt", $oeis_no);
  my $link     = sprintf("http://oeis.org/A%06d/b%06d.txt", $oeis_no, $oeis_no);

  my @ref;
  my @scr;
  {
    open my $fh, '<', $filename
        or die "Can't read $filename.\nYou should run:\n  wget $link\n";
    printf "%12s primes: reading %15s...", $name, $filename;
    while (<$fh>) {
      next unless /^(\d+)\s+(\d+)/;
      last if defined $restrict && $2 > $restrict;
      push @ref, "$2";
    }
    close $fh;
  }
  printf " %7d.", scalar @ref;

  printf "  Generating...";
  @scr = split /\s+/, qx+$FindBin::Bin/../bin/primes.pl $script_arg 1 $ref[-1]+;
  printf " %7d.\n", scalar @scr;

  die "Not equal numbers:  ", scalar @ref, " - ", scalar @scr, "\n"
    unless @ref == @scr;

  foreach my $i (0 .. $#ref) {
    die "$name prime $i not equal:  $ref[$i] - $scr[$i]\n"
        if $ref[$i] != $scr[$i];
  }
}
