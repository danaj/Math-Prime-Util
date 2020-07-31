#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/factor factor_exp divisors divisor_sum is_prime
                         prime_bigomega prime_omega/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};


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
               2231139 2463289 3008891 5115953 6961021 8030207 10486123
               10893343 12327779 701737021
               549900 10000142
               392498
              /;

my @testn64 = qw/37607912018 346065536839 600851475143
                 3204941750802 29844570422669
                 279238341033925 2623557157654233 24739954287740860
                 3369738766071892021 10023859281455311421
                 9007199254740991 9007199254740992 9007199254740993
                 6469693230 200560490130 7420738134810 304250263527210
                 13082761331670030 614889782588491410
                 440091295252541 5333042142001571
                 79127989298
                 2339796554687 124838608575421729
                 1434569741817480287 1256490565186616147
                 13356777177440210791
                /;

my @omega = ([0,1,1], [1,0,0], [2,1,1], [36,4,2], [102,3,3], [8593952,7,3]);

push @testn, @testn64 if $use64;

push @testn, qw/9999986200004761 99999989237606677 999999866000004473/
      if $use64 && $extra;

# For time savings, trim these if we're pure Perl.
if (    !$extra
     && !Math::Prime::Util::prime_get_config->{'xs'}
     && !Math::Prime::Util::prime_get_config->{'gmp'} ) {
  @testn = grep {    $_ != 10023859281455311421
                  && $_ != 3369738766071892021
                } @testn;
}

my %all_factors = (
1234567890 => [1,2,3,5,6,9,10,15,18,30,45,90,3607,3803,7214,7606,10821,11409,18035,19015,21642,22818,32463,34227,36070,38030,54105,57045,64926,68454,108210,114090,162315,171135,324630,342270,13717421,27434842,41152263,68587105,82304526,123456789,137174210,205761315,246913578,411522630,617283945,1234567890],
1032924637 => [1,6469,159673,1032924637],
4567890 => [1,2,3,5,6,10,15,30,43,86,129,215,258,430,645,1290,3541,7082,10623,17705,21246,35410,53115,106230,152263,304526,456789,761315,913578,1522630,2283945,4567890],
 456789 => [1,3,43,129,3541,10623,152263,456789],
 123456 => [1,2,3,4,6,8,12,16,24,32,48,64,96,192,643,1286,1929,2572,3858,5144,7716,10288,15432,20576,30864,41152,61728,123456],
 115553 => [1,115553],
  30107 => [1,7,11,17,23,77,119,161,187,253,391,1309,1771,2737,4301,30107],
     42 => [1,2,3,6,7,14,21,42],
     16 => [1,2,4,8,16],
     12 => [1,2,3,4,6,12],
     10 => [1,2,5,10],
      9 => [1,3,9],
      8 => [1,2,4,8],
      7 => [1,7],
      6 => [1,2,3,6],
      5 => [1,5],
      4 => [1,2,4],
      3 => [1,3],
      2 => [1,2],
      1 => [1],
      0 => [0,1],
);

my %prime_factors = (
 456789 => [3,43,3541],
 123456 => [2,2,2,2,2,2,3,643],
 115553 => [115553],
  30107 => [7,11,17,23],
      5 => [5],
      4 => [2,2],
      3 => [3],
      2 => [2],
      1 => [],
      0 => [0],
);

my %factor_exponents = (
 456789 => [[3,1],[43,1],[3541,1]],
 123456 => [[2,6],[3,1],[643,1]],
 115553 => [[115553,1]],
  30107 => [[7,1],[11,1],[17,1],[23,1]],
      5 => [[5,1]],
      4 => [[2,2]],
      3 => [[3,1]],
      2 => [[2,1]],
      1 => [],
      0 => [[0,1]],
);

plan tests => (3 * scalar @testn)
            + 2*scalar(keys %prime_factors)
            + 4*scalar(keys %all_factors)
            + 2*scalar(keys %factor_exponents)
            + 10*9  # 10 extra factoring tests * 9 algorithms
            + 8
            + 2*scalar(@omega)
            + 6;

foreach my $n (@testn) {
  my @f = factor($n);
  my $facstring = join(' * ',@f);

  # Do they multiply to the number?
  my $product = 1;  $product *= $_ for @f;
  is( $product, $n, "$n = [ $facstring ]" );

  # Are they all prime?
  my $isprime = 1; $isprime *= is_prime($_) for @f;
  if ($n < 1) {
    ok( !$isprime, "   each factor is not prime" );
  } else {
    ok(  $isprime, "   each factor is prime" );
  }

  # Does factor_exp return the appropriate rearrangement?
  is_deeply( [factor_exp($n)], [linear_to_exp(@f)], "   factor_exp looks right" );
}

while (my($n, $factors) = each(%prime_factors)) {
  is_deeply( [factor($n)], $factors, "factors($n)" );
  is( scalar factor($n), scalar @$factors, "scalar factors($n)" );
}

while (my($n, $divisors) = each(%all_factors)) {
  is_deeply( [divisors($n)], $divisors, "divisors($n)" );
  is( scalar divisors($n), scalar @$divisors, "scalar divisors($n)" );
  is( divisor_sum($n,0), scalar @$divisors, "divisor_sum($n,0)" );
  my $sum = 0;  foreach my $f (@$divisors) { $sum += $f; }
  is( divisor_sum($n), $sum, "divisor_sum($n)" );
}

while (my($n, $factors) = each(%factor_exponents)) {
  is_deeply( [factor_exp($n)], $factors, "factor_exp($n)" );
  is( scalar factor_exp($n), scalar @$factors, "scalar factor_exp($n)" );
}

extra_factor_test("trial_factor",  sub {Math::Prime::Util::trial_factor(shift)});
extra_factor_test("fermat_factor", sub {Math::Prime::Util::fermat_factor(shift)});
extra_factor_test("holf_factor",   sub {Math::Prime::Util::holf_factor(shift)});
extra_factor_test("squfof_factor", sub {Math::Prime::Util::squfof_factor(shift)});
extra_factor_test("pbrent_factor", sub {Math::Prime::Util::pbrent_factor(shift)});
extra_factor_test("prho_factor",   sub {Math::Prime::Util::prho_factor(shift)});
extra_factor_test("pminus1_factor",sub {Math::Prime::Util::pminus1_factor(shift)});
extra_factor_test("pplus1_factor", sub {Math::Prime::Util::pplus1_factor(shift)});
extra_factor_test("lehman_factor", sub {Math::Prime::Util::lehman_factor(shift)});
# TODO: old versions of MPUGMP didn't pull out factors of 3 or 5.
#extra_factor_test("ecm_factor", sub {Math::Prime::Util::ecm_factor(shift)});

# To hit some extra coverage
is_deeply( [Math::Prime::Util::trial_factor(5514109)], [2203,2503], "trial factor 2203 * 2503" );
is_deeply( [Math::Prime::Util::holf_factor(3747785838079,80000)], [1935281,1936559], "holf factor 1935281 * 1936559" );
is_deeply( [Math::Prime::Util::pminus1_factor(899,20)], [29,31], "p-1 factor 29 * 31 with tiny B1" );
is_deeply( [Math::Prime::Util::pminus1_factor(667,1000)], [23,29], "p-1 factor 23 * 29 with small B1" );
is_deeply( [Math::Prime::Util::pminus1_factor(166213)], [347,479], "p-1 factor 347 * 479" );
is_deeply( [Math::Prime::Util::pminus1_factor(563777293,1000,20000)], [23099,24407], "p-1 factor 23099 * 24407 using stage 2" );



sub extra_factor_test {
  my $fname = shift;
  my $fsub = shift;

  is_deeply( [ sort {$a<=>$b} $fsub->(1)   ], [],        "$fname(1)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(4)   ], [2, 2],    "$fname(4)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(9)   ], [3, 3],    "$fname(9)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(11)  ], [11],      "$fname(11)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(25)  ], [5, 5],    "$fname(25)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(30)  ], [2, 3, 5], "$fname(30)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(210) ], [2,3,5,7], "$fname(210)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(175) ], [5, 5, 7], "$fname(175)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(403) ], [13, 31],  "$fname(403)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(549900) ], [2,2,3,3,5,5,13,47],  "$fname(549900)" );
}

# Factor in scalar context
is( scalar factor(0), 1, "scalar factor(0) should be 1" );
is( scalar factor(1), 0, "scalar factor(1) should be 0" );
is( scalar factor(3), 1, "scalar factor(3) should be 1" );
is( scalar factor(4), 2, "scalar factor(4) should be 2" );
is( scalar factor(5), 1, "scalar factor(5) should be 1" );
is( scalar factor(6), 2, "scalar factor(6) should be 2" );
is( scalar factor(30107), 4, "scalar factor(30107) should be 4" );
is( scalar factor(174636000), 15, "scalar factor(174636000) should be 15" );

sub linear_to_exp {
  my %exponents;
  my @factors = grep { !$exponents{$_}++ } @_;
  return (map { [$_, $exponents{$_}] } @factors);
}

######
for my $ov (@omega) {
  my($n,$bo,$o) = @$ov;
  is( prime_bigomega($n), $bo, "prime_bigomega($n) = $bo" );
  is( prime_omega($n), $o, "prime_omega($n) = $o" );
}
