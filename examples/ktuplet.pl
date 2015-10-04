#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use bigint;
$|=1;
prime_set_config(verbose=>0);

# 34 minutes on Macbook Pro to find first 52 entries of A213601:
# my $high = 25777719656829367;
# my @cl = (6,10,12,16,22,24,30,34,36,40,42);
# which makes it ~4x slower than JKA's old hand-tuned code.
#
# 69 seconds on Macbook Pro for the Federighi (5TP39) sequence:
# my $high = 1e14;
# my @cl = (2,6,8,18,20,30,32,36,38);
# which comes out to about 1.5x slower than JKA's 2007 result.

my $low = 1;
my $high = 1e18;
# See https://oeis.org/A257375 for a good list
#my @cl = (6,8,14,18,20,24,26);   # A022013
#my @cl = (2,6,8,18,20,30,32,36,38);   # Federighi
#my @cl = (2,6,8,12,18,20,26,30,32,36,42,48,50,56);   # A257304
#my @cl = (4,6,10,16,18,24,28,30,34,40,46,48,54,58,60,66);   # A257375
#my @cl = (6,12,16,18,22,28,30,36,40,42,46,48);   # A214947
#my @cl = (2,6,8,12,18,20,26,30,32,36,42);  # A213645
my @cl = (6,10,12,16,22,24,30,34,36,40,42);  # A213601
#my @cl = (4,6,10,16,18,24,28,30,34,36);   # A213646
#my @cl = (2,6,8,12,18,20,26,30,32,36);   # A213647
#my @cl = (2,6,8,12,18,20,26,30,32);   # A027569
#my @cl = (2,6,12,14,20,24,26,30,32);   # A027570
#my @cl = (4,6,10,16,18,24,28,30);   # A022547
#my @cl = (4,10,12,18,22,24,28,30);   # A022548
#my @cl = (2,6,8,12,18,20,26,30);   # A022545
#my @cl = (2,6,12,14,20,24,26,30);   # A022545
#my @cl = (2,6,8,12,18,20,26);   # A022011
#my @cl = (2,6,12,14,20,24,26);   # A022012
#my @cl = (6,8,14,18,20,24,26);   # A022013
#my @cl = (2,6,8,12,18,20);   # A022009
#my @cl = (2,8,12,14,18,20);   # A022010
#my @cl = (4,6,10,12,16);   # A022008
#my @cl = (4,6,10,12);   # A022007
#my @cl = (2,6,8,12);   # A022006

my $range = 1e14;

my $i = 0;
my @p;
while ($low < $high) {
  my $chigh = $low + $range - 1;
  $chigh = $high if $chigh > $high;
  # The GMP code will use more residues so favor it with big clusters
  if (1 && scalar(@cl) > 8) {
    @p = Math::Prime::Util::GMP::sieve_prime_cluster($low, $chigh, @cl);
  } else {
    @p = sieve_prime_cluster($low, $chigh, @cl);
  }
  prime_set_config(verbose=>0);
  print ++$i," $_\n" for @p;
  $low += $range;
}
