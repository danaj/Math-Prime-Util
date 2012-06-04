#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/factor is_prime/;

my $use64 = Math::Prime::Util::_maxbits > 32;
my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};

my @testn = qw/0 1 2 3 4 5 6 7 8 16 57 64 377 9592 30107 78498 664579 5761455
               114256942 2214143 999999929 50847534 455052511 2147483647
               4118054813/;

my @testn64 = qw/37607912018 346065536839 600851475143
                 3204941750802 29844570422669
                 279238341033925 2623557157654233 24739954287740860
                 3369738766071892021 10023859281455311421
                 9007199254740991 9007199254740992 9007199254740993
                /;


push @testn, @testn64 if $use64;

push @testn, qw/9999986200004761 99999989237606677 999999866000004473/
      if $use64 && $extra;

plan tests =>  (2 * scalar @testn) + ($use64 ? 1 : 0);

if ($use64) {
  # Simple test:  perl -e 'die if 18446744073709550592 == ~0'
  my $broken = (18446744073709550592 == ~0);
  if ($broken) {
    if ($] < 5.008) {
      diag "Perl pre-5.8.0 has broken 64-bit.  Expect failures.";
    } else {
      diag "Eek!  Your 64-bit Perl $] is **** BROKEN ****.  Expect failures.";
    }
  }
  ok( !$broken, "64-bit isn't obviously broken" );
}

foreach my $n (@testn) {
  my @f = factor($n);
  my $facstring = join(' * ',@f);

  # Do they multiply to the number?
  my $product = 1;  $product *= $_ for @f;
  is( $product, $n, "$n = [ $facstring ]" );

  # Are they all prime?
  my $isprime = 1; $isprime *= is_prime($_) for @f;
  if ($n < 2) {
    ok( !$isprime, "All factors [ $facstring ] of $n are not prime" );
  } else {
    ok( $isprime, "All factors [ $facstring ] of $n are prime" );
  }
};
