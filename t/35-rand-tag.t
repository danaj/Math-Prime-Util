#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/:rand/;

plan tests => 1+4+1;

########
#
# Notes:
#   On quadmath platforms drand will use 128-bits instead of 64.  Simiarly
#   for NV=float platforms we'd see different patterns after the first.

is( srand(7652245), 7652245, "srand returns result" );
my %alg = (
  ChaCha20 => [1951677399, 598936225, 0.716442236122296401],
  ISAAC    => [2993131935, 393080975, 0.00891862162060655416],
);
my @got = ( irand, irand32, rand );
my @exp;
my $which;
for my $alg (keys %alg) {
  next if $alg{$alg}->[0] != $got[0];
  @exp = @{ $alg{$alg} };
  $which = $alg;
  last;
}

SKIP: {
  skip "Unknown PRNG algorithm",4 if !defined $which;
  is( $got[0], $exp[0], "$which irand" );
  is( $got[1], $exp[1], "$which irand32" );
  ok( $got[2] > $exp[2]-1e-6 && $got[2] < $exp[2]+1e-6, "$which drand" );
  srand(7652245);
  my($r, $want) = (irand, $got[0]);
  is( $r, $want, "Replicates after srand" );
}

SKIP: {
  skip "Unknown PRNG algorithm",1 if !defined $which;
  my $r = irand64;
  my $want = ($which eq 'ChaCha20') ? "2572411501841793573" : "1688269932343098788";
  is("$r", $want, "$which irand64");
}
