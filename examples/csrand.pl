#!/usr/bin/env perl
use strict;
use warnings;
use Math::BigInt try => "GMP,Pari";
use Math::Prime::Util qw/:all/;
use Bytes::Random::Secure;
$|=1;

# Example of Blum-Micali and Blum-Blum-Shub CSPRNGs.
# Not very practical, but works as an example.

if (!@ARGV) {
  die <<EOU;
Usage: $0 <nbits> [<type>] [<bits>]

An example showing two classic CSPRNGs (cryptographically secure pseudorandom
number generators).  These are generally not used in practice for performance
reasons, with things like AES-CTR, ISAAC, Yarrow/Fortuna, or stream ciphers
like Salsa20 instead being used.

<nbits>: how many bits should be generated.

<type>: one of "BM" (Blum-Micali) or "BBS" (Blum-Blum-Shub).  Default BBS.

<bits>: How large of primes are used for P (BM) or P*Q (BBS).  Default 512.
EOU
}

my $nbits = shift || 10;
my $type = shift || 'BBS';  # BM or BBS
my $bits = shift || 512;

die "Type must be BM or BBS" unless $type =~ /^(BBS|BM)$/;
die "Bits must be > 64" unless $bits > 64;

my $rng = Bytes::Random::Secure->new(NonBlocking => 1);
my $rbytes = int(($bits+7)/8);

if ($type eq 'BM') {
  my($p, $xn);
  # Select P
  do { $p = 2*random_nbit_prime($bits-1)+1 } while !is_prime($p);
  # Get generator
  my $g = Math::BigInt->new( "" . znprimroot($p) );
  do {  # Select the seed x0
    $xn = Math::BigInt->new("0x".$rng->bytes_hex($rbytes))->bmod($p);
  } while $xn <= 1;
  # Generate bits
  my $thresh = ($p-1) >> 1;
  while ($nbits-- > 0) {
    $xn = $g->copy->bmodpow($xn,$p);
    print 0 + ($xn < $thresh);
  }
  print "\n";
} else {
  my($M,$xn);
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
  do {  # Select the seed x0
    $xn = Math::BigInt->new("0x".$rng->bytes_hex($rbytes))->bmod($M);
  } while $xn <= 1 || gcd($xn,$M) != 1;
  # Generate bits
  my $two = Math::BigInt->new(2);
  while ($nbits-- > 0) {
    $xn->bmodpow($two,$M);
    print $xn->is_odd ? 1 : 0;
  }
  print "\n";
}
