#!/usr/bin/env perl
use strict;
use warnings;
use Math::BigInt try => "GMP,Pari";
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
    # could use $xn = powmod($g, $xn, $p);
    print 0 + ($xn < $thresh);
  }
  print "\n";
} elsif ($type eq 'BBS') {
  die "Blum-Blum-Shub must have bits >= 3500\n" unless $bits >= 3500;
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
    # Could use: $xn = mulmod($xn, $xn, $M);
    print $xn->is_odd ? 1 : 0;
  }
  print "\n";
} else {  # Micali-Schnorr
  die "Micali-Schnorr must have bits >= 120\n" unless $bits >= 120;
  my $tries = 1;
  my ($n, $e, $N);
  while (1) {
    my $p = random_nbit_prime($bits);
    my $q = random_nbit_prime($bits);
    $n = $p * $q;
    my $phi = ($p-1) * ($q-1);
    $N = length($n->as_bin)-2;
    # For efficiency, choose largest e possible.  e will always be odd.
    $e = int($N/80);
    $e-- while $e > 1 && gcd($e,$phi) != 1;
    last if $e > 1 && $e < $phi && 80*$e <= $N && gcd($e,$phi) == 1;
    die "Unable to find a proper e for MS\n" if $tries++ > 100;
  }
  my $k = int($N * (1-2/$e));
  my $r = $N - $k;
  my $xn = Math::BigInt->new("0x".$rng->bytes_hex(int(($r+7)/8)))->bmod(Math::BigInt->new(2)->bpow($r));

  while ($nbits > 0) {
    # y_i = x_{i-1} ^ e mod n
    my $yistr = $xn->copy->bmodpow($e, $n)->as_bin;
    # x_i = r most significant bits of y_i
    $xn = $xn->from_bin(substr($yistr, 0, 2+$r));
    # could do this:
    #    my $yistr = todigitstring(powmod($xn,$e,$n),2);
    #    $xn = fromdigits(substr($yistr, 0, $r),2);
    # z_i = k least significant bits of y_i
    # output is the sequence of z_i
    my $outbits = ($nbits >= $k) ? $k : $nbits;
    print substr($yistr,-$outbits);
    $nbits -= $outbits;
  }
  print "\n";
}
