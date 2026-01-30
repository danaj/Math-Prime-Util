#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/factor factor_exp divisors divisor_sum is_prime
                         prime_bigomega prime_omega/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64  = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};



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


my @testn = qw/7 8 16 57 64 377 9592 78498 664579 5761455
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
                 2339796554687
                /;

push @testn64, "124838608575421729"   if $usegmp || $usexs;
push @testn64, "1434569741817480287"  if $usegmp || $usexs;
push @testn64, "1256490565186616147"  if $usegmp || $usexs;
push @testn64, "13356777177440210791" if $usegmp || $usexs;

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

my @testdivisors = (
  [         0,],
  [         1,1],
  [         2,1,2],
  [         3,1,3],
  [         4,1,2,4],
  [         5,1,5],
  [         6,1,2,3,6],
  [         7,1,7],
  [         8,1,2,4,8],
  [         9,1,3,9],
  [        10,1,2,5,10],
  [        12,1,2,3,4,6,12],
  [        16,1,2,4,8,16],
  [        42,1,2,3,6,7,14,21,42],
  [     30107,1,7,11,17,23,77,119,161,187,253,391,1309,1771,2737,4301,30107],
  [    115553,1,115553],
  [    123456,1,2,3,4,6,8,12,16,24,32,48,64,96,192,643,1286,1929,2572,3858,5144,7716,10288,15432,20576,30864,41152,61728,123456],
  [    456789,1,3,43,129,3541,10623,152263,456789],
  [   4567890,1,2,3,5,6,10,15,30,43,86,129,215,258,430,645,1290,3541,7082,10623,17705,21246,35410,53115,106230,152263,304526,456789,761315,913578,1522630,2283945,4567890],
  [1032924637,1,6469,159673,1032924637],
  [1234567890,1,2,3,5,6,9,10,15,18,30,45,90,3607,3803,7214,7606,10821,11409,18035,19015,21642,22818,32463,34227,36070,38030,54105,57045,64926,68454,108210,114090,162315,171135,324630,342270,13717421,27434842,41152263,68587105,82304526,123456789,137174210,205761315,246913578,411522630,617283945,1234567890],
);

my @testfactors = (
  [         0, [0] ],
  [         1, [] ],
  [         2, [2] ],
  [         3, [3] ],
  [         4, [2,2] ],
  [         5, [5] ],
  [         6, [2,3] ],
  [     30107, [7,11,17,23] ],
  [    115553, [115553] ],
  [    123456, [2,2,2,2,2,2,3,643] ],
  [    456789, [3,43,3541] ],
  [ 174636000, [2,2,2,2,2,3,3,3,3,5,5,5,7,7,11] ],
);

plan tests => 4      # factor, factor_exp, including scalar
            + 2*scalar(@testn)  # factor and factor_exp
            + 2*scalar(@testdivisors)
            + 3      # extra divisors tests
            + 10*10  # 10 extra factoring tests * 10 algorithms
            + 8      # more factoring for code coverage
            + 4      # omega and bigomega
            ;

####  factor(n)  factor_exp(n)  scalar of each
{
  my @nvals = map { $_->[0] } @testfactors;
  is_deeply([map { scalar factor($_) } @nvals],
            [map { scalar @{$_->[1]} } @testfactors],
            "scalar factors(n) for @nvals");
  is_deeply([map { [factor($_)] } @nvals],
            [map { $_->[1] } @testfactors],
            "factors(n) for @nvals");
  is_deeply([map { scalar factor_exp($_) } @nvals],
            [map { scalar linear_to_exp(@{$_->[1]}) } @testfactors],
            "scalar factor_exp(n) for @nvals");
  is_deeply([map { [factor_exp($_)] } @nvals],
            [map { [linear_to_exp(@{$_->[1]})] } @testfactors],
            "factor_exp(n) for @nvals");
}
#is_deeply([map { scalar factor($_) } 0..6,30107,174636000],
#          [1,0,1,1,2,1,2,4,15], "scalar_factor(n) for 0..6,30107,174636000");

####  factor(n)  factor_exp(n)  for many numbers
foreach my $n (@testn) {
  my @f = factor($n);
  my $facstring = join(' * ',@f);

  my($ispr,$prod,$inorder) = (1,1,1);
  for my $f (@f) {
    $ispr = 0 unless is_prime($f);
    $prod *= $f;
  }
  for (1..$#f) { $inorder = 0 if $f[$_] < $f[$_-1]; }

  ok($ispr && $inorder && $prod == $n, "factor($n): $n = $facstring, all sorted primes");

  # Does factor_exp return the appropriate rearrangement?
  is_deeply([factor_exp($n)], [linear_to_exp(@f)], "factor_exp($n)" );
}



####  divisors and scalar(divisors), simple divisor_sum
foreach my $dinfo (@testdivisors) {
  my($n,@divisors) = @$dinfo;
  my $nd = scalar @divisors;
  is_deeply([scalar divisors($n),[divisors($n)]],
            [scalar @divisors,\@divisors], "divisors($n)");

  my $sum = 0;  $sum += $_ for @divisors;
  is_deeply([divisor_sum($n,0),divisor_sum($n)], [$nd,$sum], "divisor_sum($n)");
}

####  divisors with a second argument
is_deeply( [divisors(5040, 120)],
           [1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21,24,28,30,35,36,40,42,45,48,56,60,63,70,72,80,84,90,105,112,120],
           "divisors(5040, 120)" );
#is_deeply( [divisors("340282366920938463463374607431768211455", 5040)],
#           [1,3,5,15,17,51,85,255,257,641,771,1285,1923,3205,3855,4369],
#           "divisors(2^128-1, 5040)" );
is_deeply( [divisors("1208925819614629174706175", 128)],
           [1,3,5,11,15,17,25,31,33,41,51,55,75,85,93,123],
           "divisors(2^80-1, 128)" );
is_deeply( [ [divisors( 0,0)], [divisors( 0,1)],
             [divisors( 1,0)], [divisors( 1,1)], [divisors( 1,2)],
             [divisors(12,0)], [divisors(12,1)], [divisors(12,4)] ],
           [ [], [],   [], [1], [1],  [], [1], [1,2,3,4] ],
           "divisors for n 0,1,12 and k 0,1,x" );


####  test each of the underlying algorithms
extra_factor_test("trial_factor",  sub {Math::Prime::Util::trial_factor(shift)});
extra_factor_test("fermat_factor", sub {Math::Prime::Util::fermat_factor(shift)});
extra_factor_test("holf_factor",   sub {Math::Prime::Util::holf_factor(shift)});
extra_factor_test("squfof_factor", sub {Math::Prime::Util::squfof_factor(shift)});
extra_factor_test("pbrent_factor", sub {Math::Prime::Util::pbrent_factor(shift)});
extra_factor_test("prho_factor",   sub {Math::Prime::Util::prho_factor(shift)});
extra_factor_test("pminus1_factor",sub {Math::Prime::Util::pminus1_factor(shift)});
extra_factor_test("pplus1_factor", sub {Math::Prime::Util::pplus1_factor(shift)});
extra_factor_test("cheb_factor", sub {Math::Prime::Util::cheb_factor(shift)});
SKIP: {
  skip "No lehman_factor in PP", 10 unless $usexs;
  extra_factor_test("lehman_factor", sub {Math::Prime::Util::lehman_factor(shift)});
}
# TODO: old versions of MPUGMP didn't pull out factors of 3 or 5.
#extra_factor_test("ecm_factor", sub {Math::Prime::Util::ecm_factor(shift)});

# To hit some extra coverage
is_deeply( [Math::Prime::Util::trial_factor(5514109)], [2203,2503], "trial factor 2203 * 2503" );
is_deeply( [Math::Prime::Util::trial_factor(1819015037140)], [2,2,5,7,7,1856137793], "trial_factor(1819015037140) fully factors");
SKIP: {
  skip "holf_factor for 64-bit input", 1 unless $use64 || !$usexs;
  is_deeply( [Math::Prime::Util::holf_factor(3747785838079,80000)], [1935281,1936559], "holf factor 1935281 * 1936559" );
}
is_deeply( [Math::Prime::Util::pminus1_factor(166213)], [347,479], "p-1 factor 347 * 479" );
SKIP: {
  skip "p-1 tests for C code", 3 unless $usexs;
  is_deeply( [Math::Prime::Util::pminus1_factor(899,20)], [29,31], "p-1 factor 29 * 31 with tiny B1" );
  is_deeply( [Math::Prime::Util::pminus1_factor(667,1000)], [23,29], "p-1 factor 23 * 29 with small B1" );
  is_deeply( [Math::Prime::Util::pminus1_factor(563777293,1000,20000)], [23099,24407], "p-1 factor 23099 * 24407 using stage 2" );
}
# GMP still has some issues with this
#is_deeply( [Math::Prime::Util::cheb_factor("13581893559735945553",1500)], [3453481411,3932812123], "cheb factor 3453481411 * 3932812123" );
SKIP: {
  skip "cheb_factor for 64-bit input", 1 unless $use64 || !$usexs;
  is_deeply( [Math::Prime::Util::cheb_factor("2466600463243213733",1000)], [1552318819,1588978007], "cheb factor 1552318819 * 1588978007" );
}



sub extra_factor_test {
  my $fname = shift;
  my $fsub = shift;

if (0) {
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
} else {
  is_deeply( [ $fsub->(1)   ], [],        "$fname(1)" );
  is_deeply( [ $fsub->(4)   ], [2, 2],    "$fname(4)" );
  is_deeply( [ $fsub->(9)   ], [3, 3],    "$fname(9)" );
  is_deeply( [ $fsub->(11)  ], [11],      "$fname(11)" );
  is_deeply( [ $fsub->(25)  ], [5, 5],    "$fname(25)" );
  is_deeply( [ $fsub->(30)  ], [2, 3, 5], "$fname(30)" );
  is_deeply( [ $fsub->(210) ], [2,3,5,7], "$fname(210)" );
  is_deeply( [ $fsub->(175) ], [5, 5, 7], "$fname(175)" );
  is_deeply( [ $fsub->(403) ], [13, 31],  "$fname(403)" );
  is_deeply( [ $fsub->(549900) ], [2,2,3,3,5,5,13,47],  "$fname(549900)" );
}
}

sub linear_to_exp {
  my %exponents;
  my @factors = grep { !$exponents{$_}++ } @_;
  return (map { [$_, $exponents{$_}] } @factors);
}

####  prime_omega and prime_bigomega
{
  my @omegai = (qw/0 1 2 36 102 392 8593952 1309524066293741924073/);
  my @omegao = (qw/1 0 1 2  3   2   3       5/);
  my @omegab = (qw/1 0 1 4  3   5   7       7/);
  is_deeply([map { prime_omega($_)     } @omegai],\@omegao,"prime_omega(n)");
  is_deeply([map { prime_bigomega($_)  } @omegai],\@omegab,"prime_bigomega(n)");
  is_deeply([map {prime_omega('-'.$_)  } @omegai],\@omegao,"prime_omega(-n)");
  is_deeply([map {prime_bigomega('-'.$_)}@omegai],\@omegab,"prime_bigomega(-n)");
}
