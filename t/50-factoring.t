#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/factor all_factors is_prime/;

my $use64 = Math::Prime::Util::_maxbits > 32;
my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};


if ($use64) {
  # Simple test:  perl -e 'die if 18446744073709550592 == ~0'
  my $broken = (18446744073709550592 == ~0);
  if ($broken) {
    if ($] < 5.008) {
      diag "Perl pre-5.8.0 has broken 64-bit.  Skipping 64-bit tests.";
    } else {
      diag "Eek!  Your 64-bit Perl $] is **** BROKEN ****.  Skipping 64-bit tests.";
    }
    $use64 = 0;
  }
}


my @testn = qw/0 1 2 3 4 5 6 7 8 16 57 64 377 9592 30107 78498 664579 5761455
               114256942 2214143 999999929 50847534 455052511 2147483647
               4118054813
               30 210 2310 30030 510510 9699690 223092870
               1363 989 779 629 403
               547308031
               808 2727 12625 34643 134431 221897 496213 692759 1228867
               2463289 3008891 5115953 6961021 8030207 10486123
               10893343 12327779 701737021
              /;

my @testn64 = qw/37607912018 346065536839 600851475143
                 3204941750802 29844570422669
                 279238341033925 2623557157654233 24739954287740860
                 3369738766071892021 10023859281455311421
                 9007199254740991 9007199254740992 9007199254740993
                 6469693230 200560490130 7420738134810 304250263527210
                 13082761331670030 614889782588491410
                /;

push @testn, @testn64 if $use64;

push @testn, qw/9999986200004761 99999989237606677 999999866000004473/
      if $use64 && $extra;

my %all_factors = (
1234567890 => [2,3,5,6,9,10,15,18,30,45,90,3607,3803,7214,7606,10821,11409,18035,19015,21642,22818,32463,34227,36070,38030,54105,57045,64926,68454,108210,114090,162315,171135,324630,342270,13717421,27434842,41152263,68587105,82304526,123456789,137174210,205761315,246913578,411522630,617283945],
1032924637 => [6469,159673],
4567890 => [2,3,5,6,10,15,30,43,86,129,215,258,430,645,1290,3541,7082,10623,17705,21246,35410,53115,106230,152263,304526,456789,761315,913578,1522630,2283945],
 456789 => [3,43,129,3541,10623,152263],
 123456 => [2,3,4,6,8,12,16,24,32,48,64,96,192,643,1286,1929,2572,3858,5144,7716,10288,15432,20576,30864,41152,61728],
 115553 => [],
  30107 => [7,11,17,23,77,119,161,187,253,391,1309,1771,2737,4301],
     42 => [2,3,6,7,14,21],
     16 => [2,4,8],
     12 => [2,3,4,6],
     10 => [2,5],
      9 => [3],
      8 => [2,4],
      7 => [],
      6 => [2,3],
      5 => [],
      4 => [2],
      3 => [],
      2 => [],
      1 => [],
      0 => [],
);

plan tests =>  (2 * scalar @testn) + scalar(keys %all_factors) + 6*7;

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

while (my($n, $divisors) = each(%all_factors)) {
  is_deeply( [all_factors($n)], $divisors, "all_factors($n)" );
}


extra_factor_test("trial_factor",  sub {Math::Prime::Util::trial_factor(shift)});
extra_factor_test("fermat_factor", sub {Math::Prime::Util::fermat_factor(shift)});
extra_factor_test("holf_factor",   sub {Math::Prime::Util::holf_factor(shift)});
extra_factor_test("squfof_factor", sub {Math::Prime::Util::squfof_factor(shift)});
extra_factor_test("pbrent_factor", sub {Math::Prime::Util::pbrent_factor(shift)});
extra_factor_test("prho_factor",   sub {Math::Prime::Util::prho_factor(shift)});
extra_factor_test("pminus1_factor",sub {Math::Prime::Util::pminus1_factor(shift)});

sub extra_factor_test {
  my $fname = shift;
  my $fsub = shift;

  is_deeply( [ sort {$a<=>$b} $fsub->(1)   ], [1],       "$fname(1)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(4)   ], [2, 2],    "$fname(4)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(9)   ], [3, 3],    "$fname(9)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(25)  ], [5, 5],    "$fname(9)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(175) ], [5, 5, 7], "$fname(175)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(403) ], [13, 31],  "$fname(403)" );
}

