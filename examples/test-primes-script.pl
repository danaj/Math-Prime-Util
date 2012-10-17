#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec::Functions;
use FindBin;
$|++; #flush the output buffer after every write() or print() function

test_oeis(1359, "Twin", "--twin");
test_oeis(5385, "Safe", "--safe");
test_oeis(5384, "SG", "--sophie");
test_oeis(5479, "Lucas", "--lucas");
test_oeis(5478, "Fibonacci", "--fibonacci");
test_oeis(2385, "Palindromic", "--palin");

sub test_oeis {
  my($oeis_no, $name, $script_arg) = @_;

  my $filename = sprintf("b%06d.txt", $oeis_no);
  my $link     = sprintf("http://oeis.org/A%06d/b%06d.txt", $oeis_no, $oeis_no);

  my @ref;
  my @scr;
  {
    open my $fh, '<', $filename
        or die "Can't read $filename.\nYou should run:\n  wget $link\n";
    while (<$fh>) {
      next unless /^(\d+)\s+(\d+)/;
      last if $2 >= 10_000_000_000 && $oeis_no == 2385;
      push @ref, $2;
    }
    close $fh;
  }
  warn "Read ", scalar @ref, " $name primes from reference file\n";

  @scr = split /\s+/, qx+$FindBin::Bin/../bin/primes.pl $script_arg 1 $ref[-1]+;
  warn "Generated ", scalar @scr, " $name primes using primes.pl\n";

  die "Not equal numbers:  ", scalar @ref, " - ", scalar @scr, "\n"
    unless @ref == @scr;

  foreach my $i (0 .. $#ref) {
    die "$name prime $i not equal:  $ref[$i] - $scr[$i]\n"
        if $ref[$i] != $scr[$i];
  }
}
