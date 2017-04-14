#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/:rand/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

plan tests => 1+4+1+1;

########
#
# Notes:
#   32-bit Perls won't have irand64 properly available.  We should test that it
#   exists and returns rands, but it will only get a single word.
#
#   On quadmath platforms drand will use 128-bits instead of 64.  Simiarly
#   for NV=float platforms we'd see different patterns after the first.

is( srand(7652245), 7652245, "srand returns result" );
my %alg = (
  ChaCha20 => [3004121508, 1648872885, 0.295097488345488454],
  ISAAC    => [1201784256, 2493844312, 0.4502675995501714],
);
my @got = ( irand, irand, rand );
my @exp;
my $which;
for my $alg (keys %alg) {
  next if $alg{$alg}->[0] != $got[0];
  @exp = @{ $alg{$alg} };
  $which = $alg;
  last;
}
ok( scalar(@exp) == scalar(@got), "PRNG Algorithm is recognized" );
is( $got[0], $exp[0], "$which irand" );
is( $got[1], $exp[1], "$which irand" );
ok( $got[2] > $exp[2]-1e-6 && $got[2] < $exp[2]+1e-6, "$which drand" );

srand(7652245);
my($r, $want) = (irand, $got[0]);
is( $r, $want, "Replicates after srand" );

SKIP: {
  skip "Skipping irand64 on 32-bit Perl", 1 if !$use64;
  $r = irand64;
  $want = ($which eq 'ChaCha20') ? 7081855117603603021 : 10710979763289504966;
  is($r, $want, "$which irand64");
}
