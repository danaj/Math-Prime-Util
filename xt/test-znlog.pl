#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/random_nbit_prime znprimroot znlog powmod/;

# This test uses znlog with prime modulus, similar to FLINT's tests.
# Our znlog will handle non-primes, so another interesting test would
# be using random integer p values.
my $ebits = 63;

print "1..",$ebits-3,"\n";
for my $bits (4 .. $ebits) {
  #Math::Prime::Util::prime_set_config(verbose=>1) if $bits > 50;
  my $ntests = ($bits < 30) ? 100 : ($bits < 40) ? 10 : ($bits < 45) ? 5 : 1;
  my $ok = "ok";
  for my $tn (1 .. $ntests) {
    my $p = random_nbit_prime( $bits );
    my $root = znprimroot($p);
    my $b = int(rand($p-1)) + 1;
    my $d = znlog($b, $root, $p);
    my $res = powmod($root, $d, $p);
    next if $res == $b;
    $ok = "not ok";
    warn "FAIL $bits:  $root ^ $d mod $p = $res, not $b\n";
    print "PASS $bits: $root ^ $d mod $p = $res\n";
  }
  print "$ok ",$bits-3," - znlog with $bits bits\n";
}
