#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
  unless ($ENV{RELEASE_TESTING}) {
    require Test::More;
    Test::More::plan(skip_all => 'these tests are in development');
  }
}

use Test::More;
use Math::Prime::Util qw/:rand/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

plan tests => 3 + 2;

########
#
# Notes:
#   32-bit Perls won't have irand64 properly available.  We should test that it
#   exists and returns rands, but it will only get a single word.
#
#   On quadmath platforms drand will use 128-bits instead of 64.  Simiarly
#   for NV=float platforms we'd see different patterns after the first.

my $r;

is( srand(7652245), 7652245, "srand returns result" );
$r = rand();
ok( $r >= .4999999 && $r <= .5000001, "rand returns expected result after seed" );
$r = irand();
is( $r, 4153545820, "irand returns expected result after seed" );

srand(7652245);
$r = rand();
$r = irand();
is( $r, 4153545820, "Replicates sequence after srand" );
$r = irand64();
SKIP: {
  skip "Skipping irand64 on 32-bit Perl", 1 if !$use64;
  is( $r, 16376814160539732392, "irand64 returns expected result" );
}
