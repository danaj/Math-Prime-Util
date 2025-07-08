#!/usr/bin/env perl
use strict;
use warnings;

# This is a subset of our tests.  You really should run the whole test suite
# on the PP code.  What this will do is basic regression testing.
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = ~0 > 4294967295 && ~0 != 18446744073709550592;

use Test::More;
my @small_primes = qw/
2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71
73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173
179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281
283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409
419 421 431 433 439 443 449 457 461 463 467 479 487 491 499 503 509 521 523 541
547 557 563 569 571 577 587 593 599 601 607 613 617 619 631 641 643 647 653 659
661 673 677 683 691 701 709 719 727 733 739 743 751 757 761 769 773 787 797 809
811 821 823 827 829 839 853 857 859 863 877 881 883 887 907 911 919 929 937 941
947 953 967 971 977 983 991 997 1009 1013 1019 1021 1031 1033 1039 1049 1051
1061 1063 1069
/;  # next prime is 1087

my @primes = qw/
1129 1327 9551 15683 19609 31397 155921
5 11 29 97 127 541 907 1151 1361 9587 15727 19661 31469 156007 360749
370373 492227 1349651 1357333 2010881 4652507 17051887 20831533 47326913
122164969 189695893 191913031
/;

my @composites = qw/
0 4 6 8 9 10 12 14 15 16 18 20 21 22
9 2047 1373653 25326001 3215031751
561 1105 1729 2465 2821 6601 8911 10585 15841 29341 41041 46657 52633
62745 63973 75361 101101 340561 488881 852841 1857241 6733693
9439201 17236801 23382529 34657141 56052361 146843929
341 561 645 1105 1387 1729 1905 2047 2465 2701 2821 3277 4033 4369 4371
4681 5461 6601 7957 8321 52633 88357
66066 173645446 7500135 115501463
/;

# pseudoprimes to various small prime bases
my %pseudoprimes = (
   2 => [ qw/2047   42799 4335241 1078467589/ ],
   3 => [ qw/121    44287 4252381 1075490821/ ],
   5 => [ qw/781    38081 4265257 1075156291/ ],
   7 => [ qw/25 325 29857 4411681 1074439981/ ],
  11 => [ qw/133    43213 4224533 1076929261/ ],
  13 => [ qw/85     35371 4336879 1079159203/ ],
  17 => [ qw/9 91   71071 4224533 1076237119/ ],
  19 => [ qw/9 49   49771 4384693 1074718783/ ],
  23 => [ qw/169    25201 4219129 1079063371/ ],
  29 => [ qw/15 91  48133 4219129 1075151447/ ],
  31 => [ qw/15 49  29341 4270657 1073833843/ ],
  37 => [ qw/9 451  59563 4287817 1075430539/ ],
  61 => [ qw/217    79381 4219129 1079326249/ ],
  73 => [ qw/205    34219 4321153 1074220489/ ],
 psp2       => [ qw/341 561       29341 4259905 1073823745/ ],
 psp3       => [ qw/91 121        44287 4252381 1073827147/ ],
 lucas      => [ qw/9179 10877    44099 4259789 1074039119/ ],
 slucas     => [ qw/5459 5777     75077 4309631 1080085439/ ],
 eslucas    => [ qw/989 3239 5777 72389 4226777 1076503199/ ],
 aeslucas1  => [ qw/989 10469     39059 4269341 1076503199/ ],
 aeslucas2  => [ qw/4531 12209    62479 4403027 1074695441/ ],
);
# Test a pseudoprime larger than 2^32.
push @{$pseudoprimes{2}}, 75792980677 if $use64;
my $num_pseudoprimes = 0;
foreach my $ppref (values %pseudoprimes) {
  push @composites, @$ppref;
  $num_pseudoprimes += scalar @$ppref;
}
{
  my %uniq;
  $uniq{$_}++ for (@composites);
  @composites = sort {$a<=>$b} keys %uniq;
}

my %small_single = (
    0   => [],
    1   => [],
    2   => [2],
    3   => [2, 3],
    4   => [2, 3],
    5   => [2, 3, 5],
    6   => [2, 3, 5],
    7   => [2, 3, 5, 7],
    11  => [2, 3, 5, 7, 11],
    18  => [2, 3, 5, 7, 11, 13, 17],
    19  => [2, 3, 5, 7, 11, 13, 17, 19],
    20  => [2, 3, 5, 7, 11, 13, 17, 19],
);

my %small_range = (
  "3 to 9" => [3,5,7],
  "2 to 20" => [2,3,5,7,11,13,17,19],
  "30 to 70" => [31,37,41,43,47,53,59,61,67],
  "70 to 30" => [],
  "20 to 2" => [],
  "1 to 1" => [],
  "2 to 2" => [2],
  "3 to 3" => [3],
  "2 to 3" => [2,3],
  "2 to 5" => [2,3,5],
  "3 to 6" => [3,5],
  "3 to 7" => [3,5,7],
  "4 to 8" => [5,7],
  "2010733 to 2010881" => [2010733,2010881],
  "2010734 to 2010880" => [],
  "3088 to 3164" => [3089,3109,3119,3121,3137,3163],
  "3089 to 3163" => [3089,3109,3119,3121,3137,3163],
  "3090 to 3162" => [3109,3119,3121,3137],
  "3842610773 to 3842611109" => [3842610773,3842611109],
  "3842610774 to 3842611108" => [],
);

my %primegaps = (
 19609 => 52,
 360653 => 96,
 2010733 => 148,
);

my %pivals32 = (
                   1 => 0,
                  10 => 4,
                 100 => 25,
                1000 => 168,
               10000 => 1229,
              100000 => 9592,
             1000000 => 78498,
            10000000 => 664579,
           100000000 => 5761455,
          1000000000 => 50847534,
               60067 => 6062,
               65535 => 6542,
            16777215 => 1077871,
          2147483647 => 105097565,
          4294967295 => 203280221,
);
my %pivals_small = map { $_ => $pivals32{$_} }
                   grep {$_ <= 80000}
                   keys %pivals32;

my %pi_intervals = (
  "1e9 +2**14"  => 785,
  "17 to 13"    => 0,
  "3 to 17"     => 6,
  "4 to 17"     => 5,
  "4 to 16"     => 4,
  "191912783 +248" => 2,
  "191912784 +247" => 1,
  "191912783 +247" => 1,
  "191912784 +246" => 0,
);
my %extra_pi_intervals = (
  "868396 to 9478505" => 563275,
  "1118105 to 9961674" => 575195,
  "24689 to 7973249" => 535368,
);
# Add extra intervals to pi_intervals if we're doing release testing
@pi_intervals{keys %extra_pi_intervals} = values %extra_pi_intervals if $extra;

# Remove any entries where the high value is too large for us
# ikegami++ for the delete from a hash slice idea
delete @pi_intervals{ grep { (parse_range($_))[1] > ~0 } keys %pi_intervals };

my %nthprimes32 = (
                  1 => 2,
                 10 => 29,
                100 => 541,
               1000 => 7919,
              10000 => 104729,
             100000 => 1299709,
            1000000 => 15485863,
           10000000 => 179424673,
          100000000 => 2038074743,
);
my %nthprimes_small = map { $_ => $nthprimes32{$_} }
                      grep { $extra ? ($_ <= 2_000_000) : ($_ <= 5_000) }
                      keys %nthprimes32;

my %eivals = (
         -10 =>  -0.00000415696892968532438,
        -0.5 =>  -0.55977359477616,
        -0.1 =>  -1.8229239584193906660809,
      -0.001 =>  -6.33153936413615,
    -0.00001 => -10.9357198000436956,
 -0.00000001 => -17.843465089050832587,
 0.693147180559945 => 1.0451637801174927848446,           # log2
         1   =>  1.8951178163559367554665,
         1.5 =>  3.3012854491297978379574,
         2   =>  4.9542343560018901633795,
         5   =>  40.185275355803177455091,
         10  =>  2492.2289762418777591384,
         12  =>  14959.532666397528852292,
         20  =>  25615652.664056588820481,
         40  =>  6039718263611241.5783592,
         41  =>  16006649143245041.110700,
);
my %livals = (
              0 =>  0,
           1.01 => -4.0229586739299358695031,
              2 =>  1.0451637801174927848446,
             10 =>  6.1655995047872979375230,
             24 =>  11.200315795232698830550,
           1000 =>  177.60965799015222668764,
         100000 =>  9629.8090010507982050343,
      100000000 =>  5762209.3754480314675691,
     4294967295 =>  203284081.95454158906409,
    10000000000 =>  455055614.58662307560953,
   100000000000 =>  4118066400.6216115150394,
);
my %rvals = (
           1.01 =>  1.0060697180622924796117,
              2 =>  1.5410090161871318832885,
             10 =>  4.5645831410050902398658,
           1000 =>  168.35944628116734806491,
        1000000 =>  78527.399429127704858870,
       10000000 =>  664667.44756474776798535,
     4294967295 =>  203280697.51326064541983,
    10000000000 =>  455050683.30684692446315,
18446744073709551615 => 4.25656284014012122706963685602e17,
);
my %rzvals = (
            2   =>  0.6449340668482264364724151666,
            2.5 =>  0.3414872572509171797567696934,
            4.5 =>  0.0547075107614542640229672890,
            7   =>  0.0083492773819228268397975498,
            8.5 =>  0.0028592508824156277133439825,
           20.6 =>  0.0000006293391573578212882457,
           80   =>  8.27180612553034e-25,
          180   =>  6.52530446799852e-55,
);
my %ipp = (
  5 => 2,
  10 => 0,
  49 => 0,
  347 => 2,
  697 => 0,
  7080233 => 2,
  7080249 => 0,
  17471059 => 2,
  17471061 => 0,
  36010357 => 2,
  36010359 => 0,
);

plan tests => 2 +
              3 +
              3 + scalar(keys %small_single) + scalar(keys %small_range) +
              4 +                 # sieve_range
              2*scalar(keys %primegaps) + 8 + 1 + 1 + 1 +
              scalar(keys %pivals_small) + scalar(keys %pi_intervals) +
              6 +                 # PC, pc approx
              2*scalar(keys %pivals_small) + scalar(keys %nthprimes_small) +
              4 + scalar(keys %pseudoprimes) +
              scalar(keys %eivals) + scalar(keys %livals) + scalar(keys %rvals) + scalar(keys %rzvals) +
              ($extra ? 4 : 0) +  # Bigfloat RiemannZeta
              1 + 1 +             # factor
              10 + 7*3 +          # factoring subs
              1 +                 # HOLF
              ($extra ? 3 : 0) +  # HOLF extra
              ($extra ? 3 : 0) +  # factor stage 2
              10 +                # AKS
              ($use64 ? 3 : 2) +  # Lucas and BLS75 primality proofs
              6 +                 # M-R and Lucas on bigint
              2 +                 # PC and NP approx
              72 +                # Misc util.pm functions
              ($extra ? 1 : 0) +  # twin prime count approx
              scalar(keys %ipp) + # is_prob_prime
              7 +                 # _is_gaussian_prime
              2 +                 # is_semiprime
              2 +                 # is_totient
              1 +                 # ramanujan_sum
              5 +                 # is_carmichael
              1;

use Math::Prime::Util qw/primes
                         prime_count_approx nth_prime_approx
                         prime_get_config prime_set_config
                         consecutive_integer_lcm
                         primorial pn_primorial partitions miller_rabin_random
                         is_prob_prime
                        /;
use Math::BigInt;
use Math::BigFloat;
require_ok 'Math::Prime::Util::PP';
require_ok 'Math::Prime::Util::PrimalityProving';

    # This function skips some setup
    undef *primes;
    *primes             = \&Math::Prime::Util::PP::primes;
    *sieve_range        = \&Math::Prime::Util::PP::sieve_range;

    *prime_count        = \&Math::Prime::Util::PP::prime_count;
    *prime_count_lower  = \&Math::Prime::Util::PP::prime_count_lower;
    *prime_count_upper  = \&Math::Prime::Util::PP::prime_count_upper;
    *nth_prime          = \&Math::Prime::Util::PP::nth_prime;
    undef *prime_count_approx;
    undef *nth_prime_approx;
    *prime_count_approx = \&Math::Prime::Util::PP::prime_count_approx;
    *nth_prime_approx   = \&Math::Prime::Util::PP::nth_prime_approx;

    *twin_prime_count   = \&Math::Prime::Util::PP::twin_prime_count;
    *nth_twin_prime     = \&Math::Prime::Util::PP::nth_twin_prime;
    *twin_prime_count_approx = \&Math::Prime::Util::PP::twin_prime_count_approx;

    *is_prime       = \&Math::Prime::Util::PP::is_prime;
    *next_prime     = \&Math::Prime::Util::PP::next_prime;
    *prev_prime     = \&Math::Prime::Util::PP::prev_prime;

    *is_pseudoprime = \&Math::Prime::Util::PP::is_pseudoprime;
    *is_strong_pseudoprime = \&Math::Prime::Util::PP::is_strong_pseudoprime;
    *is_lucas_pseudoprime = \&Math::Prime::Util::PP::is_lucas_pseudoprime;
    *is_strong_lucas_pseudoprime = \&Math::Prime::Util::PP::is_strong_lucas_pseudoprime;
    *is_extra_strong_lucas_pseudoprime = \&Math::Prime::Util::PP::is_extra_strong_lucas_pseudoprime;
    *is_almost_extra_strong_lucas_pseudoprime = \&Math::Prime::Util::PP::is_almost_extra_strong_lucas_pseudoprime;
    *is_frobenius_underwood_pseudoprime = \&Math::Prime::Util::PP::is_frobenius_underwood_pseudoprime;
    *is_perrin_pseudoprime = \&Math::Prime::Util::PP::is_perrin_pseudoprime;
    *is_frobenius_pseudoprime = \&Math::Prime::Util::PP::is_frobenius_pseudoprime;
    *is_aks_prime      = \&Math::Prime::Util::PP::is_aks_prime;
    *is_gaussian_prime = \&Math::Prime::Util::PP::is_gaussian_prime;

    *is_ramanujan_prime= \&Math::Prime::Util::PP::is_ramanujan_prime;
    *nth_ramanujan_prime = \&Math::Prime::Util::PP::nth_ramanujan_prime;
    *inverse_totient   = \&Math::Prime::Util::PP::inverse_totient;
    *is_semiprime      = \&Math::Prime::Util::PP::is_semiprime;
    *is_totient        = \&Math::Prime::Util::PP::is_totient;
    *is_carmichael     = \&Math::Prime::Util::PP::is_carmichael;

    *factor            = \&Math::Prime::Util::PP::factor;

    *gcd               = \&Math::Prime::Util::PP::gcd;
    *lcm               = \&Math::Prime::Util::PP::lcm;

    *moebius           = \&Math::Prime::Util::PP::moebius;
    *euler_phi         = \&Math::Prime::Util::PP::euler_phi;
    *jordan_totient    = \&Math::Prime::Util::PP::jordan_totient;
    *mertens           = \&Math::Prime::Util::PP::mertens;
    *exp_mangoldt      = \&Math::Prime::Util::PP::exp_mangoldt;
    *chebyshev_theta   = \&Math::Prime::Util::PP::chebyshev_theta;
    *chebyshev_psi     = \&Math::Prime::Util::PP::chebyshev_psi;
    *ramanujan_sum     = \&Math::Prime::Util::PP::ramanujan_sum;

    *Pi                = \&Math::Prime::Util::PP::Pi;
    *znprimroot        = \&Math::Prime::Util::PP::znprimroot;
    *znorder           = \&Math::Prime::Util::PP::znorder;
    *znlog             = \&Math::Prime::Util::PP::znlog;
    *binomial          = \&Math::Prime::Util::PP::binomial;
    *stirling          = \&Math::Prime::Util::PP::stirling;
    *bernfrac          = \&Math::Prime::Util::PP::bernfrac;
    *valuation         = \&Math::Prime::Util::PP::valuation;
    *gcdext            = \&Math::Prime::Util::PP::gcdext;
    *invmod            = \&Math::Prime::Util::PP::invmod;
    *vecmin            = \&Math::Prime::Util::PP::vecmin;
    *vecmax            = \&Math::Prime::Util::PP::vecmax;
    *vecsum            = \&Math::Prime::Util::PP::vecsum;
    *vecprod           = \&Math::Prime::Util::PP::vecprod;
    *liouville         = \&Math::Prime::Util::PP::liouville;
    *carmichael_lambda = \&Math::Prime::Util::PP::carmichael_lambda;
    *forperm           = \&Math::Prime::Util::PP::forperm;
    *forcomb           = \&Math::Prime::Util::PP::forcomb;
    *forpart           = \&Math::Prime::Util::PP::forpart;
    *forprimes         = \&Math::Prime::Util::PP::forprimes;
    *forcomposites     = \&Math::Prime::Util::PP::forcomposites;
    *foroddcomposites  = \&Math::Prime::Util::PP::foroddcomposites;
    *forfactored       = \&Math::Prime::Util::PP::forfactored;
    *fordivisors       = \&Math::Prime::Util::PP::fordivisors;
    # TODO test this
    #*forsemiprimes     = \&Math::Prime::Util::PP::forsemiprimes;
    #*forsquarefree     = \&Math::Prime::Util::PP::forsquarefree;
    #*forsquarefreeint  = \&Math::Prime::Util::PP::forsquarefreeint;

    *RiemannR            = \&Math::Prime::Util::PP::RiemannR;
    *RiemannZeta         = \&Math::Prime::Util::PP::RiemannZeta;
    *LogarithmicIntegral = \&Math::Prime::Util::PP::LogarithmicIntegral;
    *ExponentialIntegral = \&Math::Prime::Util::PP::ExponentialIntegral;
    *LambertW            = \&Math::Prime::Util::PP::LambertW;

# Turn off use of BRS - ECM tries to use this.
# prime_set_config( irand => sub { int(rand(4294967296.0)) } );

###############################################################################

$_ = 'this should not change';

{
  my %small_primes = map { $_ => 1 } @small_primes;
  my @isprime = map { is_prime($_) } (0 .. 1086);
  my @exprime = map { $small_primes{$_} ? 2 : 0 } (0 .. 1086);
  is_deeply( \@isprime, \@exprime, "is_prime 0 .. 1086" );
}
{
  my @isprime = map { is_prime($_) ? "$_ is prime" : "$_ is composite" }
                @primes, @composites;
  my @exprime =  map { "$_ is prime" } @primes;
  push @exprime, map { "$_ is composite" } @composites;
  is_deeply( \@isprime, \@exprime, "is_prime for selected numbers" );
}

is_deeply( Math::Prime::Util::PP::trial_primes(80),
           [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79],
           "Trial primes 2-80" );

###############################################################################

is_deeply( primes(1069), \@small_primes, "Primes between 0 and 1069" );
is_deeply( primes(1070), \@small_primes, "Primes between 0 and 1070" );
is_deeply( primes(1086), \@small_primes, "Primes between 0 and 1086" );

while (my($high, $expect) = each (%small_single)) {
  is_deeply( primes($high), $expect, "primes($high) should return [@{$expect}]")
;
}

while (my($range, $expect) = each (%small_range)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is_deeply( primes($low, $high), $expect, "primes($low,$high) should return [@{$expect}]");
}

is_deeply( [sieve_range(4, 4, 1)], [map { $_-4 } 4,5,6,7], "sieve range depth 1" );
is_deeply( [sieve_range(10, 20, 2)], [1,3,5,7,9,11,13,15,17,19], "sieve range depth 2" );
is_deeply( [sieve_range(10, 20, 3)], [1,3,7,9,13,15,19], "sieve range depth 3" );
is_deeply( [sieve_range(10, 20, 5)], [1,3,7,9,13,19], "sieve range depth 5" );

###############################################################################

while (my($base, $range) = each (%primegaps)) {
  is( next_prime($base), $base+$range, "next prime of $base is $base+$range" );
  is( prev_prime($base+$range), $base, "prev prime of $base+$range is $base" );
}

is( next_prime(19608), 19609, "next prime of 19608 is 19609" );
is( next_prime(19610), 19661, "next prime of 19610 is 19661" );
is( next_prime(19660), 19661, "next prime of 19660 is 19661" );
is( prev_prime(19662), 19661, "prev prime of 19662 is 19661" );
is( prev_prime(19660), 19609, "prev prime of 19660 is 19609" );
is( prev_prime(19610), 19609, "prev prime of 19610 is 19609" );

is( prev_prime(2), undef, "Previous prime of 2 returns undef" );
{
  my $n   = $use64 ? 18446744073709551611 : 4294967291;
  my $exp = $use64 ? "18446744073709551629" : "4294967311";
  my $got = next_prime($n);
  ok(ref($got) =~ /^Math::/ && "$got" eq $exp, "next_prime(~0-4) returns bigint result");
}

{
  my @exprime = map { "next_prime($_) == 2010881" }         (2010733..2010880);
  my @isprime = map { "next_prime($_) == ".next_prime($_) } (2010733..2010880);
  is_deeply(\@isprime, \@exprime, "next_prime for 148 primes before primegap end 2010881");
}
{
  my @exprime = map { "prev_prime($_) == 2010733" }         (2010734..2010881);
  my @isprime = map { "prev_prime($_) == ".prev_prime($_) } (2010734..2010881);
  is_deeply(\@isprime, \@exprime, "prev_prime for 148 primes before primegap start 2010733");
}
# Similar test case to 2010870, where m=0 and next_prime is at m=1
is(next_prime(1234567890), 1234567891, "next_prime(1234567890) == 1234567891)");

###############################################################################

while (my($n, $pin) = each (%pivals_small)) {
  is( prime_count($n), $pin, "Pi($n) = $pin" );
}

while (my($range, $expect) = each (%pi_intervals)) {
  my($low,$high) = parse_range($range);
  is( prime_count($low,$high), $expect, "prime_count($range) = $expect");
}

# These are small enough they should be exact.
is( prime_count_lower(450), 87, "prime_count_lower(450)" );
is( prime_count_upper(450), 87, "prime_count_upper(450)" );
# Make sure these are about right
cmp_closeto( prime_count_lower(1234567), 95360, 60, "prime_count_lower(1234567) in range" );
cmp_closeto( prime_count_upper(1234567), 95360, 60, "prime_count_upper(1234567) in range" );
cmp_closeto( prime_count_lower(412345678), 21958997, 1500, "prime_count_lower(412345678) in range" );
cmp_closeto( prime_count_upper(412345678), 21958997, 1500, "prime_count_upper(412345678) in range" );

###############################################################################

while (my($n, $pin) = each (%pivals_small)) {
  my $next = $pin+1;
  cmp_ok( $pin ? nth_prime($pin) : 0, '<=', $n, "nth_prime($pin) <= $n");
  cmp_ok( nth_prime($next), '>=', $n, "nth_prime($next) >= $n");
}

while (my($n, $nth) = each (%nthprimes_small)) {
  is( nth_prime($n), $nth, "nth_prime($n) = $nth" );
}

###############################################################################

is( is_strong_pseudoprime(0, 2), 0, "MR with 0 shortcut composite");
is( is_strong_pseudoprime(1, 2), 0, "MR with 0 shortcut composite");
is( is_strong_pseudoprime(2, 2), 1, "MR with 2 shortcut prime");
is( is_strong_pseudoprime(3, 2), 1, "MR with 3 shortcut prime");

while (my($base, $ppref) = each (%pseudoprimes)) {
  my $npseudos = scalar @$ppref;
  my @expmr = map { 1 } @$ppref;
  my @gotmr;
  if ($base =~ /^psp(\d+)/) {
     my $pbase = $1;
     @gotmr = map { is_pseudoprime($_, $pbase) } @$ppref;
  } elsif ($base =~ /^aeslucas(\d+)/) {
     my $inc = $1;
     @gotmr = map { is_almost_extra_strong_lucas_pseudoprime($_, $inc) } @$ppref;
  } elsif ($base eq 'eslucas') {
     @gotmr = map { is_extra_strong_lucas_pseudoprime($_) } @$ppref;
  } elsif ($base eq 'slucas') {
     @gotmr = map { is_strong_lucas_pseudoprime($_) } @$ppref;
  } elsif ($base eq 'lucas') {
     @gotmr = map { is_lucas_pseudoprime($_) } @$ppref;
  } else {
     @gotmr = map { is_strong_pseudoprime($_, $base) } @$ppref;
  }
  is_deeply(\@gotmr, \@expmr, "$npseudos pseudoprimes (base $base)");
}

###############################################################################

while (my($n, $ein) = each (%eivals)) {
  cmp_closeto( ExponentialIntegral($n), $ein, 0.00000001 * abs($ein), "Ei($n) ~= $ein");
}
while (my($n, $lin) = each (%livals)) {
  cmp_closeto( LogarithmicIntegral($n), $lin, 0.00000001 * abs($lin), "li($n) ~= $lin");
}
while (my($n, $rin) = each (%rvals)) {
  cmp_closeto( RiemannR($n), $rin, 0.00000001 * abs($rin), "R($n) ~= $rin");
}
while (my($n, $zin) = each (%rzvals)) {
  cmp_closeto( RiemannZeta($n), $zin, 0.00000001 * abs($zin), "Zeta($n) ~= $zin");
}
cmp_closeto( LambertW(6588), 6.86636957140619, 0.000000001, "LambertW(6588)");
if ($extra) {
  my ($n, $zin);
  ($n, $zin) = (4.5, $rzvals{4.5});
  cmp_closeto( RiemannZeta(Math::BigFloat->new($n)), $zin, 0.00000001 * abs($zin), "Zeta($n) ~= $zin");
  ($n, $zin) = (20.6, $rzvals{20.6});
  cmp_closeto( RiemannZeta(Math::BigFloat->new($n)), $zin, 0.00000001 * abs($zin), "Zeta($n) ~= $zin");
  ($n, $zin) = (80, $rzvals{80});
  cmp_closeto( RiemannZeta(Math::BigFloat->new($n)), $zin, 0.00000001 * abs($zin), "Zeta($n) ~= $zin");
  ($n, $zin) = (180, $rzvals{180});
  cmp_closeto( RiemannZeta(Math::BigFloat->new($n)), $zin, 0.00000001 * abs($zin), "Zeta($n) ~= $zin");
}

###############################################################################

#foreach my $n (@primes) {
#  my @f = factor($n);
#  is_deeply( \@f, [$n], "factor prime $n yields $n" );
#}
{
  my $ntests = scalar @primes;
  my @expfactor = map { "$_" } @primes;
  my @gotfactor = map { join(' * ', factor($_)) } @primes;
  is_deeply( \@gotfactor, \@expfactor, "test factoring for $ntests primes");
}
{
  my $ntests = scalar @composites;
  my @expfactor = map { "$_ factored correctly" } @composites;
  my @gotfactor;

  foreach my $n (@composites) {
    my @f = factor($n);
    my $facstring = join(' * ', @f);

    if ($n < 2) {
      push @gotfactor, (@f == 1 && $f[0] == $n)
                       ? "$n factored correctly"
                       : "$n not correct: $facstring";
      next;
    }
    my $product  = 1; $product = int($product * $_) for @f;
    my $allprime = 1; $allprime *= is_prime($_) for @f;
    if (@f >= 2 && $product == $n && $allprime) {
      push @gotfactor, "$n factored correctly";
    } else {
      push @gotfactor, "$n not correct: $facstring";
    }
  }
  is_deeply( \@gotfactor, \@expfactor, "test factoring for $ntests composites");
}

# The PP factor code does small trials, then loops doing 64k rounds of HOLF
# if the composite is less than a half word, followed by 64k rounds each of
# prho with a = {3,5,7,11,13}.  Most numbers are handled by these.  The ones
# that aren't end up being too slow for us to put in a test.  So we'll try
# running the various factoring methods manually.
{
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::holf_factor(403) ],
             [ 13, 31 ],
             "holf(403)" );
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::fermat_factor(403) ],
             [ 13, 31 ],
             "fermat(403)" );
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::prho_factor(403) ],
             [ 13, 31 ],
             "prho(403)" );
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::pbrent_factor(403) ],
             [ 13, 31 ],
             "pbrent(403)" );
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::pminus1_factor(403) ],
             [ 13, 31 ],
             "pminus1(403)" );
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::prho_factor(851981) ],
             [ 13, 65537 ],
             "prho(851981)" );
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::pbrent_factor(851981) ],
             [ 13, 65537 ],
             "pbrent(851981)" );
  #is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::ecm_factor(851981) ],
  #           [ 13, 65537 ],
  #           "ecm(851981)" );
  # Try to force using stage 2.
  SKIP: {
    skip "Skipping ecm stage 2 tests", 1 if defined $Math::Prime::Util::GMP::VERSION && $Math::Prime::Util::GMP::VERSION < 0.20;
    is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::ecm_factor(101303039, 5, 100000,100) ],
               [ 1013, 100003 ],
               "ecm(101303039)" );
  }
  my $n64 = $use64 ? 55834573561 : Math::BigInt->new("55834573561");
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::prho_factor($n64) ],
             [ 13, 4294967197 ],
             "prho(55834573561)" );
  is_deeply( [ sort {$a<=>$b} Math::Prime::Util::PP::pbrent_factor($n64) ],
             [ 13, 4294967197 ],
             "pbrent(55834573561)" );
  # 1013 4294967197 4294967291
  my $nbig = Math::BigInt->new("18686551294184381720251");
  my @nfac;
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::prho_factor($nbig);
  is(scalar @nfac, 2, "prho finds a factor of 18686551294184381720251");
  is($nfac[0] * $nfac[1], $nbig, "prho found a correct factor");
  ok($nfac[0] != 1 && $nfac[1] != 1, "prho didn't return a degenerate factor");
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::pbrent_factor($nbig);
  is(scalar @nfac, 2, "pbrent finds a factor of 18686551294184381720251");
  is($nfac[0] * $nfac[1], $nbig, "pbrent found a correct factor");
  ok($nfac[0] != 1 && $nfac[1] != 1, "pbrent didn't return a degenerate factor");
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::pminus1_factor($nbig);
  is(scalar @nfac, 2, "pminus1 finds a factor of 18686551294184381720251");
  is($nfac[0] * $nfac[1], $nbig, "pminus1 found a correct factor");
  ok($nfac[0] != 1 && $nfac[1] != 1, "pminus1 didn't return a degenerate factor");
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::ecm_factor($nbig);
  SKIP: {
    skip "Skipping ecm test", 3 if defined $Math::Prime::Util::GMP::VERSION && $Math::Prime::Util::GMP::VERSION < 0.20;
    is(scalar @nfac, 2, "ecm finds a factor of 18686551294184381720251");
    is($nfac[0] * $nfac[1], $nbig, "ecm found a correct factor");
    ok($nfac[0] != 1 && $nfac[1] != 1, "ecm didn't return a degenerate factor");
  }

  $nbig = Math::BigInt->new("73786976930493367637");
  # Check stage 2 p-1.  Fast with Math::BigInt::GMP, slow without.
  SKIP: {
    skip "Skipping p-1 stage 2 tests", 3 unless $extra;
    @nfac = sort {$a<=>$b} Math::Prime::Util::PP::pminus1_factor($nbig, 27000, 35000);
    is(scalar @nfac, 2, "pminus1 finds a factor of 73786976930493367637");
    is($nfac[0] * $nfac[1], $nbig, "pminus1 found a correct factor");
    ok($nfac[0] != 1 && $nfac[1] != 1, "pminus1 didn't return a degenerate factor");
  }
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::fermat_factor($nbig);
  is(scalar @nfac, 2, "fermat finds a factor of 73786976930493367637");
  is($nfac[0] * $nfac[1], $nbig, "fermat found a correct factor");
  ok($nfac[0] != 1 && $nfac[1] != 1, "fermat didn't return a degenerate factor");
  if ($extra) {
    @nfac = sort {$a<=>$b} Math::Prime::Util::PP::holf_factor($nbig);
    is(scalar @nfac, 2, "holf finds a factor of 18686551294184381720251");
    is($nfac[0] * $nfac[1], $nbig, "holf found a correct factor");
    ok($nfac[0] != 1 && $nfac[1] != 1, "holf didn't return a degenerate factor");
  }
  {
    $nbig = Math::BigInt->new("99999999999979999998975857");
    @nfac = sort {$a<=>$b} Math::Prime::Util::PP::holf_factor($nbig);
    is_deeply(\@nfac, [9999999998987,10000000001011], "holf correctly factors 99999999999979999998975857");
  }
  SKIP: {
    # Unfortunately we can't guarantee this isn't found in stage 1.
    skip "ecm stage 2", 3 unless $extra;
    $nbig = Math::BigInt->new("14270401808568703916861");
    @nfac = sort {$a<=>$b} Math::Prime::Util::PP::ecm_factor($nbig, 5, 2000, 40);
    is(scalar @nfac, 2, "ecm(5,2000) finds a factor of 14270401808568703916861");
    is($nfac[0] * $nfac[1], $nbig, "ecm(5,2000) found a correct factor");
    ok($nfac[0] != 1 && $nfac[1] != 1, "ecm(5,2000) didn't return a degenerate factor");
  }
}

##### Some numbers that go to stage 2 of tests
if ($extra) {
  my $nbig = Math::BigInt->new("9087500560545072247139");
  my @nfac;
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::pminus1_factor($nbig,1000,10000);
  is_deeply( [@nfac], ["24133","376559091722747783"], "p-1 stage 2 finds factors of $nbig" );
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::trial_factor($nbig, 50000);
  is_deeply( [@nfac], ["24133","376559091722747783"], "trial factor finds factors of $nbig" );
  @nfac = sort {$a<=>$b} Math::Prime::Util::PP::ecm_factor($nbig, 10,1000,100);
  is_deeply( [@nfac], ["24133","376559091722747783"], "ecm factor finds factors of $nbig" );
}

##### AKS primality test.  Be very careful with performance.
is( is_aks_prime(1), 0, "AKS: 1 is composite (less than 2)" );
is( is_aks_prime(2), 1, "AKS: 2 is prime" );
is( is_aks_prime(3), 1, "AKS: 3 is prime" );
is( is_aks_prime(4), 0, "AKS: 4 is composite" );
is( is_aks_prime(64), 0, "AKS: 64 is composite (perfect power)" );
is( is_aks_prime(65), 0, "AKS: 65 is composite (caught in trial)" );
is( is_aks_prime(23), 1, "AKS: 23 is prime (r >= n)" );
is( is_aks_prime(70747), 0, "AKS: 70747 is composite (n mod r)" );
SKIP: {
  skip "Skipping PP AKS test without EXTENDED_TESTING", 2 unless $extra;
  diag "32-bit Perl will be very slow for AKS" unless $use64;
  is( is_aks_prime(1009), 1, "AKS: 1009 is prime (passed anr test)" );
  is( is_aks_prime(74513), 0, "AKS: 74513 is composite (failed anr test)" );
}

is_deeply( [Math::Prime::Util::PrimalityProving::primality_proof_lucas(100003)],
           [2, "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN 100003\n\nType Lucas\nN 100003\nQ[1] 2\nQ[2] 3\nQ[3] 7\nQ[4] 2381\nA 2\n"],
           "primality_proof_lucas(100003)" );
# Had to reduce these to make borked up Perl 5.6.2 work.
#is_deeply( [Math::Prime::Util::PP::primality_proof_bls75("210596120454733723")],
#           [2, ["210596120454733723", "n-1", [2, 3, 82651, "47185492693"], [2, 2, 2, 2]]],
#           "primality_proof_bls75(210596120454733723)" );
is_deeply( [Math::Prime::Util::PrimalityProving::primality_proof_bls75(1490266103)],
           [2, "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN 1490266103\n\nType BLS5\nN 1490266103\nQ[1] 13\nQ[2] 19\nQ[3] 1597\nQ[4] 1889\nA[0] 5\n----\n"],
           "primality_proof_bls75(1490266103)" );
if ($use64) {
  is_deeply( [Math::Prime::Util::PrimalityProving::primality_proof_bls75(27141057803)],
           [2, "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN 27141057803\n\nType BLS5\nN 27141057803\nQ[1] 47533\nQ[2] 285497\n----\n"],
           "primality_proof_bls75(27141057803)" );
}

{
  my $n = Math::BigInt->new("168790877523676911809192454171451");
  is( is_strong_pseudoprime( $n, 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47), 1, "168790877523676911809192454171451 looks prime with bases 2..52" );
  is( is_strong_pseudoprime( $n, 53), 0, "168790877523676911809192454171451 found composite with base 53" );
  is( is_strong_lucas_pseudoprime($n), 0, "168790877523676911809192454171451 is not a strong Lucas pseudoprime" );
  SKIP: {
    skip "Old Perl+bigint segfaults in F-U code", 1 if $] < 5.008;
    is( is_frobenius_underwood_pseudoprime($n), 0, "168790877523676911809192454171451 is not a Frobenius pseudoprime" );
  }
  #is( is_perrin_pseudoprime($n), 0, "168790877523676911809192454171451 is not a Perrin pseudoprime" );
  is(is_perrin_pseudoprime(517697641), 1, "517697641 is a Perrin pseudoprime");
  is(is_frobenius_pseudoprime(517697641), 0, "517697641 is not a Frobenius pseudoprime");
}

{
  my $ntha = nth_prime_approx(1287248);
  ok( $ntha >= 20274907 && $ntha <= 20284058, "nth_prime_approx(1287248) in range" );
  my $pca = prime_count_approx(128722248);
  ok( $pca >= 7309252 && $pca <= 7310044, "prime_count_approx(128722248) in range" );
}

{
  # Test some functions usually not tested in Util.pm
  my $xs  = prime_get_config->{'xs'};
  my $gmp = prime_get_config->{'gmp'};
  my $verbose = prime_get_config->{'verbose'};
  prime_set_config(xs=>0, gmp=>0);

  is( consecutive_integer_lcm(13), 360360, "consecutive_integer_lcm(13)" );
  is( consecutive_integer_lcm(52), Math::BigInt->new("3099044504245996706400"), "consecutive_integer_lcm(52)" );

  is_deeply( [moebius(513,537)],
             [qw/0 1 1 0 1 -1 1 0 -1 0 -1 0 0 1 1 0 0 -1 0 0 1 -1 1 0 1/],
             "moebius(513,537)" );
  is( moebius(42199), 1, "moebius(42199)" );
  is( liouville(444456), 1, "liouville(444456)" );
  is( liouville(562894), -1, "liouville(562894)" );

  is( mertens(4219), -13, "mertens(4219)" );

  is_deeply( [euler_phi(1513,1537)],
             [qw/1408 756 800 756 1440 440 1260 576 936 760 1522 504 1200 648 1016 760 1380 384 1530 764 864 696 1224 512 1456/],
             "euler_phi(1513,1537)" );
  is( euler_phi(324234), 108072, "euler_phi(324234)" );
  is( jordan_totient(4, 899), "653187225600", "jordan_totient(4, 899)" );
  is( carmichael_lambda(324234), 18012, "carmichael_lambda(324234)" );

  is( exp_mangoldt(16), 2, "exp_mangoldt of power of 2 = 2" );
  is( exp_mangoldt(14), 1, "exp_mangoldt of even = 1" );
  is( exp_mangoldt(21), 1, "exp_mangoldt of 21 = 1" );
  is( exp_mangoldt(23), 23, "exp_mangoldt of 23 = 23" );
  is( exp_mangoldt(27), 3, "exp_mangoldt of 27 (3^3) = 3" );

  is_deeply( [map { scalar znprimroot($_) } (-11, 0, 8, 3, 1729, 10, 5109721)],
             [2, undef, undef, 2, undef, 3, 94],
             "znprimroot" );

  is(znorder(2,35), 12, "znorder(2,35) = 12");
  is(znorder(7,35), undef, "znorder(7,35) = undef");
  is(znorder(67,999999749), 30612237, "znorder(67,999999749) = 30612237");

  is(znlog(5678, 5, 10007), 8620, "znlog(5678, 5, 10007)");

  is(binomial(35,16), 4059928950, "binomial(35,16)");
  is(binomial(228,12), "30689926618143230620", "binomial(228,12)");
  is(binomial(-23,-26), -2300, "binomial(-23,-26) should be -2300");

  is(stirling(12,4,2), '611501', "S(12,4)" );
  is(stirling(12,4,1), '105258076', "s(12,4)" );

  is_deeply( [bernfrac(0)], [1,1], "bernfrac(0)" );
  is_deeply( [bernfrac(1)], [1,2], "bernfrac(1)" );
  is_deeply( [bernfrac(2)], [1,6], "bernfrac(2)" );
  is_deeply( [bernfrac(3)], [0,1], "bernfrac(3)" );
  is_deeply( [bernfrac(12)], [-691,2730], "bernfrac(12)" );
  is_deeply( [bernfrac(13)], [0,1], "bernfrac(12)" );

  is_deeply( [gcdext(23948236,3498248)], [2263, -15492, 52], "gcdext(23948236,3498248)" );

  is( valuation(1879048192,2), 28, "valuation(1879048192,2)");
  is( valuation(96552,6), 3, "valuation(96552,6)");

  is(invmod(45,59), 21, "invmod(45,59)");
  is(invmod(14,28474), undef, "invmod(14,28474)");
  is(invmod(42,-2017), 1969, "invmod(42,-2017)");

  is(vecsum(15, 30, 45), 90, "vecsum(15,30,45)");
  is("".vecsum(4294966296,4294965296,4294964296), "12884895888", "vecsum(2^32-1000,2^32-2000,2^32-3000)");
  is(vecprod(15, 30, 45), 20250, "vecprod(15,30,45)");
  is("".vecprod(4294966296,4294965296,4294964296), "79228051833847139970490254336", "vecprod(2^32-1000,2^32-2000,2^32-3000)");
  is(vecmin(4294966296,4294965296,4294964296), 4294964296, "vecmin(2^32-1000,2^32-2000,2^32-3000)");
  is(vecmax(4294966296,4294965296,4294964296), 4294966296, "vecmax(2^32-1000,2^32-2000,2^32-3000)");

  cmp_closeto( chebyshev_theta(7001), 6929.27483821865062, 0.006929, "chebyshev_theta(7001) =~ 6929.2748");
  cmp_closeto( chebyshev_psi(6588), 6597.07452996633704, 0.006597, "chebyshev_psi(6588) =~ 6597.07453");

  while (my($n, $isp) = each (%ipp)) {
    is( is_prob_prime($n), $isp, "is_prob_prime($n) should be $isp" );
  }

  is(is_ramanujan_prime(41),1,"41 is a Ramanujan prime");
  is(is_ramanujan_prime(43),0,"43 is not a Ramanujan prime");
  is(nth_ramanujan_prime(28),311,"R_n(28) = 311");
  is_deeply( Math::Prime::Util::PP::ramanujan_primes(0,100), [2,11,17,29,41,47,59,67,71,97], "Ramanujan primes under 100");

  is(inverse_totient(42), 4, "inverse totient 42 count");
  is_deeply([inverse_totient(42)], [43,49,86,98], "inverse totient 42 list");

  is( primorial(24), 223092870, "primorial(24)" );
  is( "".primorial(118), "31610054640417607788145206291543662493274686990", "primorial(118)" );
  is( pn_primorial(7), 510510, "pn_primorial(7)" );
  is( partitions(74), 7089500, "partitions(74)" );
  is( miller_rabin_random(4294967281, 20), "0", "Miller-Rabin random 40 on composite" );

  { my @t;
    forprimes(sub {push @t,$_}, 2387234,2387303);
    is_deeply( [@t], [2387237,2387243,2387249,2387269,2387291,2387299,2387303],
               "forprimes 2387234,2387303" );
  }
  { my @t;
    forcomposites(sub {push @t,$_}, 15202630,15202641);
    is_deeply( [@t], [15202630,15202632,15202634,15202635,15202636,15202638,15202640,15202641], "forcomposites 15202630,15202641" );
  }
  { my @t;
    foroddcomposites(sub {push @t,$_}, 15202630,15202641);
    is_deeply( [@t], [15202635,15202641], "foroddcomposites 15202630,15202641" );
  }
  { my $k = 0;
    fordivisors(sub {$k += $_+int(sqrt($_))},92834);
    is( $k, 168921, "fordivisors: d|92834: k+=d+int(sqrt(d))" );
  }
  { my @t;
    forfactored(sub {push @t,[@_]}, 15202630,15202641);
    is_deeply( \@t, [[2,5,433,3511],[15202631],[2,2,2,3,23,27541],[15202633],[2,37,205441],[3,5,7,67,2161],[2,2,41,92699],[15202637],[2,3,3,11,76781],[15202639],[2,2,2,2,5,307,619],[3,17,19,29,541]], "forfactored" );
  }

  { my @p;
    forcomb(sub { push @p, [@_] }, 3, 2);
    is_deeply( \@p, [ [0,1], [0,2], [1,2] ], "forcomb(3,2)" );
  }
  { my @p;
    forperm(sub { push @p, [@_] }, 3);
    is_deeply( \@p, [ [0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0] ], "forperm(3)" );
  }
  { my @p;
    forpart(sub { push @p, [@_] }, 4);
    is_deeply( \@p, [ [1,1,1,1],[1,1,2],[1,3],[2,2],[4] ], "forpart(4)" );
  }

  is( Pi(82), "3.141592653589793238462643383279502884197169399375105820974944592307816406286208999", "Pi(82)" );

  is( gcd(-30,-90,90), 30, "gcd(-30,-90,90) = 30" );
  is( lcm(11926,78001,2211), 2790719778, "lcm(11926,78001,2211) = 2790719778" );

  is( twin_prime_count(4321), 114, "twin_prime_count(4321)" );
  cmp_closeto( twin_prime_count_approx(Math::BigInt->new("4123456784123")), "6950213327", 14937 * 2, "twin_prime_count_approx(4123456784123)" );
  if ($extra) {
    cmp_closeto( twin_prime_count_approx(Math::BigInt->new("412345678412345678412345678")), "149939117920176008847283", 1e11, "twin_prime_count_approx(412345678412345678412345678)" );
  }
  is( nth_twin_prime(249), 13217, "nth_twin_prime(249)" );

  prime_set_config(xs=>$xs, gmp=>$gmp, verbose=>$verbose);
}

ok( !is_gaussian_prime(29,0), "29 is not a Gaussian Prime" );
ok(  is_gaussian_prime(31,0), "31 is a Gaussian Prime" );
ok( !is_gaussian_prime(0,-29), "0-29i is not a Gaussian Prime" );
ok(  is_gaussian_prime(0,-31), "0-31i is a Gaussian Prime" );
ok(  is_gaussian_prime(58924,132000511), "58924+132000511i is a Gaussian Prime" );
ok(  is_gaussian_prime(519880,-2265929), "519880-2265929i is a Gaussian Prime" );
ok( !is_gaussian_prime(20571,150592260), "20571+150592260i is not a Gaussian Prime" );

ok(  is_semiprime(1110000001), "1110000001 is a semiprime" );
ok( !is_semiprime(1110000201), "1110000201 is not a semiprime" );

ok(  is_totient(381554124), "381554124 is a totient" );
ok( !is_totient(1073024875), "1073024875 is not a semiprime" );

is( ramanujan_sum(12,36), 4, "ramanujan_sum(12,36) = 4" );

ok( !is_carmichael(5049), "5049 is not a Carmichael number" );
ok( !is_carmichael(2792834247), "2792834247 is not a Carmichael number" );
ok( !is_carmichael(2399550475), "2399550475 is not a Carmichael number" );
ok( !is_carmichael(219389), "219389 is not a Carmichael number" );
ok(  is_carmichael(1125038377), "1125038377 is a Carmichael number" );

# TODO:
#  is_quasi_carmichael
#  is_pillai
#  is_fundamental
#  divisor_sum
#  inverse_li
#  _inverse_R
#  prime_count_lower
#  prime_count_upper
#  semiprime_count
#  ramanujan_prime_count
#  twin_prime_count_approx
#  semiprime_count_approx
#  nth_twin_prime_approx
#  nth_semiprime_approx
#  nth_ramanujan_prime_upper
#  nth_ramanujan_prime_lower
#  nth_ramanujan_prime_approx
#  ramanujan_prime_count_upper
#  ramanujan_prime_count_lower
#  ramanujan_prime_count_approx
#  sum_primes
#  print_primes
#  chinese
#  vecextract
#  sumdigits
#  sqrtmod
#  mulmod
#  divmod
#  powmod
#  is_power
#  is_square
#  is_prime_power
#  is_polygonal
#  hammingweight
#  todigitstring / _splitdigits
#  sqrtint
#  rootint
#  logint
#  harmfrac
#  harmreal
#  is_euler_pseudoprime
#  is_euler_plumb_pseudoprime
#  factorialmod
#  is_primitive_root
#  znorder bigint
#  non-GMP lucas_sequence
#  lucas_u
#  lucas_v
#  is_frobenius_khashin_pseudoprime
#  is_perrin_pseudoprime (restrict == 1, restrict == 2)
#  is_catalan_pseudoprime
#  is_mersenne_prime
#  is_aks_prime (_poly_mod_mul, etc)
#  pminus1_factor stage 2
#  ecm_factor no GMP
#  divisors
#  hclassno
#  ramanujan_tau (_taup, _tauprime, _taupower)
#  LogarithmicIntegral ...
#  RiemannR no GMP
#  LambertW no GMP
#  forcompositions with hash
#  numtoperm
#  permtonum
#  randperm
#  shuffle
#  urandomb
#  urandomm
#  random_ndigit_prime
#  random_safe_prime
#  random_strong_prime
#  random_maurer_prime
#  random_shawe_taylor_prime
#  miller_rabin_random
#  random_factored_integer


is( $_, 'this should not change', "Nobody clobbered \$_" );

###############################################################################

sub parse_range {
  my($range) = @_;
  my($low,$high);
  my $fixnum = sub {
    my $nstr = shift;
    $nstr =~ s/^(\d+)e(\d+)$/$1*(10**$2)/e;
    $nstr =~ s/^(\d+)\*\*(\d+)$/$1**$2/e;
    die "Unknown string in test" unless $nstr =~ /^\d+$/;
    $nstr;
  };
  if ($range =~ /(\S+)\s+to\s+(\S+)/) {
    $low = $fixnum->($1);
    $high = $fixnum->($2);
  } elsif ($range =~ /(\S+)\s*\+\s*(\S+)/) {
    $low = $fixnum->($1);
    $high = $low + $fixnum->($2);
  } else {
    die "Can't parse test data";
  }
  ($low,$high);
}

sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  cmp_ok( abs($got - $expect), '<=', $tolerance, $message );
}
