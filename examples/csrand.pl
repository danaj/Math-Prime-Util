#!/usr/bin/env perl
use strict;
use warnings;
use Math::BigInt try => "GMP,Pari";
use Math::Prime::Util qw/:all/;
$|=1;

# Example of Blum-Micali and Blum-Blum-Shub CSPRNGs.
# Not very practical, but works as an example.

my $nbits = shift || 10;
my $type = shift || 'BBS';  # BM or BBS
my $bits = shift || 512;

die "Type must be BM or BBS" unless $type =~ /^(BBS|BM)$/;

if ($type eq 'BM') {
  my($p, $x0);
  # Select P
  do { $p = 2*random_nbit_prime($bits-1)+1 } while !is_prime($p);
  # Get generator
  my $g = Math::BigInt->new( "" . znprimroot($p) );
  # Select X0.  This could be done better.
  do { $x0 = random_nbit_prime($bits) ^ (random_nbit_prime($bits) >> 1) }
    while $x0 <= 1;
  # Generate bits
  my $xn = Math::BigInt->new("$x0");
  my $thresh = ($p-1) >> 1;
  while ($nbits-- > 0) {
    $xn = $g->copy->bmodpow($xn,$p);
    print 0 + ($xn < $thresh);
  }
  print "\n";
} else {
  my($M,$x0);
  # Select M = p*q
  while (1) {
    my($p,$q);
    do { $p = random_nbit_prime($bits); } while ($p % 4) != 3;
    do { $q = random_nbit_prime($bits); } while ($q % 4) != 3;
    if ($bits < 200) {
      my $gcd = gcd(euler_phi($p-1),euler_phi($q-1));
      next if $gcd > 10000;
    }
    $M = $p * $q;
    last;
  }
  # Select X0.  This could be done better.
  do { $x0 = random_nbit_prime($bits) ^ (random_nbit_prime($bits) >> 1) }
    while $x0 <= 1 || gcd($x0,$M) != 1;
  # Generate bits
  my $xn = Math::BigInt->new("$x0");
  my $two = Math::BigInt->new(2);
  while ($nbits-- > 0) {
    $xn->bmodpow($two,$M);
    print $xn->is_odd ? 1 : 0;
  }
  print "\n";
}
