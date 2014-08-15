#!/usr/bin/env perl
use strict;
use warnings;
use Math::GMP;
use Math::Prime::Util qw/:all/;
use Bytes::Random::Secure;
$|=1;

# Example of Blum-Micali, Blum-Blum-Shub, and Micali-Schnorr CSPRNGs.
# Not very practical, but works as an example.

if (!@ARGV) {
  die <<EOU;
Usage: $0 <nbits> [<type>] [<bits>]

An example showing two classic CSPRNGs (cryptographically secure pseudorandom
number generators).  These are generally not used in practice for performance
reasons, with things like AES-CTR, ISAAC, Yarrow/Fortuna, or stream ciphers
like Salsa20 instead being used.

<nbits>: how many bits should be generated.

<type>: one of:
    "MS"  (Micali-Schnorr)       <- default
    "BM"  (Blum-Micali)
    "BBS" (Blum-Blum-Shub)

<bits>: How large of primes are used for P (BM) or P,Q (BBS,MS).  Default 512.
EOU
}

my $nbits = shift || 10;
my $type = shift || 'MS';  # BM or BBS or MS
my $bits = shift;

die "Type must be BM, BBS, or MS" unless $type =~ /^(BBS|BM|MS)$/;
if (!defined $bits) { $bits = ($type eq 'BBS') ? 4096 : 512; }
die "Bits must be > 64" unless $bits > 64;

my $rng = Bytes::Random::Secure->new(NonBlocking => 1);
my $rbytes = int(($bits+7)/8);

if ($type eq 'BM') {
  my($p, $xn);
  # Select P
  do { $p = 2 * Math::GMP->new(random_nbit_prime($bits-1))+1 } while !is_prime($p);
  # Get generator
  my $g = Math::GMP->new(znprimroot($p));
  do {  # Select the seed x0
    $xn = Math::GMP->new($rng->bytes_hex($rbytes), 16) % $p;
  } while $xn <= 1;
  # Generate bits
  my $thresh = Math::GMP::div_2exp_gmp($p-1, 1);
  while ($nbits-- > 0) {
    $xn = Math::GMP::powm_gmp($g, $xn, $p);
    print 0 + ($xn < $thresh);
  }
  print "\n";
} elsif ($type eq 'BBS') {
  die "Blum-Blum-Shub must have bits >= 3500\n" unless $bits >= 3500;
  my($M,$xn);
  # Select M = p*q
  while (1) {
    my($p,$q);
    do { $p = Math::GMP->new(random_nbit_prime($bits)); } while ($p % 4) != 3;
    do { $q = Math::GMP->new(random_nbit_prime($bits)); } while ($q % 4) != 3;
    if ($bits < 200) {
      my $gcd = gcd(euler_phi($p-1),euler_phi($q-1));
      next if $gcd > 10000;
    }
    $M = $p * $q;
    last;
  }
  do {  # Select the seed x0
    $xn = Math::GMP->new($rng->bytes_hex($rbytes), 16) % $M;
  } while $xn <= 1 || gcd($xn,$M) != 1;
  # Generate bits
  my $two = Math::GMP->new(2);
  while ($nbits-- > 0) {
    $xn = Math::GMP::powm_gmp($xn, $two, $M);
    print Math::GMP::gmp_tstbit($xn,0) ? "1" : "0";
  }
  print "\n";
} else {  # Micali-Schnorr
  die "Micali-Schnorr must have bits >= 120\n" unless $bits >= 120;
  my $tries = 1;
  my ($n, $e, $N);
  while (1) {
    my $p = Math::GMP->new(random_nbit_prime($bits));
    my $q = Math::GMP->new(random_nbit_prime($bits));
    $n = $p * $q;
    my $phi = ($p-1) * ($q-1);
    $N = Math::GMP::sizeinbase_gmp($n, 2);
    # For efficiency, choose largest e possible.  e will always be odd.
    $e = int($N/80);
    $e-- while $e > 1 && gcd($e,$phi) != 1;
    last if $e > 1 && $e < $phi && 80*$e <= $N && gcd($e,$phi) == 1;
    die "Unable to find a proper e for MS\n" if $tries++ > 100;
  }
  my $k = int($N * (1-2/$e));
  my $r = $N - $k;
  my $xn = Math::GMP->new($rng->bytes_hex(int(($r+7)/8)),16) % (Math::GMP->new(2) ** $r);
  my $twok = Math::GMP->new(2) ** $k;

  while ($nbits > 0) {
    # y_i = x_{i-1} ^ e mod n
    my $yi = Math::GMP::powm_gmp($xn, $e, $n);

    # x_i = r most significant bits of y_i
    $xn = Math::GMP::div_2exp_gmp($yi, $k);   # $xn = $yi >> $k;

    # z_i = k least significant bits of y_i
    # output is the sequence of z_i
    $twok = Math::GMP->new(2) ** $nbits if $nbits < $k;
    print Math::GMP::get_str_gmp( $yi % $twok, 2);
    $nbits -= $k;
  }
  print "\n";
}
