#!/usr/bin/env perl
use strict;
use warnings;

# This is a subset of our tests.  You really should run the whole test suite
# on the PP code.  What this will do is basic regression testing.

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
my %small_primes;
map { $small_primes{$_} = 1; } @small_primes;

my @primes = qw/
1129 1327 9551 15683 19609 31397 155921
5 11 29 97 127 541 907 1151 1361 9587 15727 19661 31469 156007 360749
370373 492227 1349651 1357333 2010881 4652507 17051887 20831533 47326913
122164969 189695893 191913031
/;

my @composites = qw/
0 1 4 6 8 9 10 12 14 15 16 18 20 21 22
9 2047 1373653 25326001 3215031751
561 1105 1729 2465 2821 6601 8911 10585 15841 29341 41041 46657 52633
62745 63973 75361 101101 340561 488881 852841 1857241 6733693
9439201 17236801 23382529 34657141 56052361 146843929
341 561 645 1105 1387 1729 1905 2047 2465 2701 2821 3277 4033 4369 4371
4681 5461 6601 7957 8321 52633 88357
/;

# pseudoprimes to various small prime bases
my %pseudoprimes = (
   2 => [ qw/2047 3277 4033 4681 8321 15841 29341 42799 49141 52633 65281 74665 80581 85489 88357 90751/ ],
   3 => [ qw/121 703 1891 3281 8401 8911 10585 12403 16531 18721 19345 23521 31621 44287 47197 55969 63139 74593 79003 82513 87913 88573 97567/ ],
   5 => [ qw/781 1541 5461 5611 7813 13021 14981 15751 24211 25351 29539 38081 40501 44801 53971 79381/ ],
   7 => [ qw/25 325 703 2101 2353 4525 11041 14089 20197 29857 29891 39331 49241 58825 64681 76627 78937 79381 87673 88399 88831/ ],
  11 => [ qw/133 793 2047 4577 5041 12403 13333 14521 17711 23377 43213 43739 47611 48283 49601 50737 50997 56057 58969 68137 74089 85879 86347 87913 88831/ ],
  13 => [ qw/85 1099 5149 7107 8911 9637 13019 14491 17803 19757 20881 22177 23521 26521 35371 44173 45629 54097 56033 57205 75241 83333 85285 86347/ ],
  17 => [ qw/9 91 145 781 1111 2821 4033 4187 5365 5833 6697 7171 15805 19729 21781 22791 24211 26245 31621 33001 33227 34441 35371 38081 42127 49771 71071 74665 77293 78881 88831 96433 97921 98671/ ],
  19 => [ qw/9 49 169 343 1849 2353 2701 4033 4681 6541 6697 7957 9997 12403 13213 13747 15251 16531 18769 19729 24761 30589 31621 31861 32477 41003 49771 63139 64681 65161 66421 68257 73555 96049/ ],
  23 => [ qw/169 265 553 1271 2701 4033 4371 4681 6533 6541 7957 8321 8651 8911 9805 14981 18721 25201 31861 34133 44173 47611 47783 50737 57401 62849 82513 96049/ ],
  29 => [ qw/15 91 341 469 871 2257 4371 4411 5149 6097 8401 11581 12431 15577 16471 19093 25681 28009 29539 31417 33001 48133 49141 54913 79003/ ],
  31 => [ qw/15 49 133 481 931 6241 8911 9131 10963 11041 14191 17767 29341 56033 58969 68251 79003 83333 87061 88183/ ],
  37 => [ qw/9 451 469 589 685 817 1333 3781 8905 9271 18631 19517 20591 25327 34237 45551 46981 47587 48133 59563 61337 68101 68251 73633 79381 79501 83333 84151 96727/ ],
  61 => [ qw/217 341 1261 2701 3661 6541 6697 7613 13213 16213 22177 23653 23959 31417 50117 61777 63139 67721 76301 77421 79381 80041/ ],
  73 => [ qw/205 259 533 1441 1921 2665 3439 5257 15457 23281 24617 26797 27787 28939 34219 39481 44671 45629 64681 67069 76429 79501 93521/ ],
);
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
                   grep {$_ <= 2000000}
                   keys %pivals32;

my %pi_intervals = (
  "868396 to 9478505" => 563275,
  "1118105 to 9961674" => 575195,
  "24689 to 7973249" => 535368,
  "1e10 +2**16" => 2821,
  "17 to 13"    => 0,
  "3 to 17"     => 6,
  "4 to 17"     => 5,
  "4 to 16"     => 4,
  "191912783 +248" => 2,
  "191912784 +247" => 1,
  "191912783 +247" => 1,
  "191912784 +246" => 0,
);
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
                      grep {$_ <= 2000000}
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


plan tests => 1 +
              2*(1087 + @primes + @composites) +
              3 + scalar(keys %small_single) + scalar(keys %small_range) +
              2*scalar(keys %primegaps) + 8 + 148 + 148 + 1 +
              3*scalar(keys %pivals32) + scalar(keys %pivals_small) + scalar(keys %pi_intervals) +
              2*scalar(keys %pivals_small) + 3*scalar(keys %nthprimes32) + scalar(keys %nthprimes_small) +
              4 + $num_pseudoprimes +
              scalar(keys %eivals) + scalar(keys %livals) + scalar(keys %rvals) +
              scalar @primes + 3*scalar @composites +
              0;

use Math::Prime::Util qw/primes/;
require_ok 'Math::Prime::Util::PP';
    # This function skips some setup
    undef *primes;
    *primes             = \&Math::Prime::Util::PP::primes;

    *prime_count        = \&Math::Prime::Util::PP::prime_count;
    *prime_count_upper  = \&Math::Prime::Util::PP::prime_count_upper;
    *prime_count_lower  = \&Math::Prime::Util::PP::prime_count_lower;
    *prime_count_approx = \&Math::Prime::Util::PP::prime_count_approx;
    *nth_prime          = \&Math::Prime::Util::PP::nth_prime;
    *nth_prime_upper    = \&Math::Prime::Util::PP::nth_prime_upper;
    *nth_prime_lower    = \&Math::Prime::Util::PP::nth_prime_lower;
    *nth_prime_approx   = \&Math::Prime::Util::PP::nth_prime_approx;

    *is_prime       = \&Math::Prime::Util::PP::is_prime;
    *next_prime     = \&Math::Prime::Util::PP::next_prime;
    *prev_prime     = \&Math::Prime::Util::PP::prev_prime;

    *miller_rabin   = \&Math::Prime::Util::PP::miller_rabin;
    *is_prob_prime  = \&Math::Prime::Util::PP::is_prob_prime;

    *factor         = \&Math::Prime::Util::PP::factor;

    *RiemannR            = \&Math::Prime::Util::PP::RiemannR;
    *LogarithmicIntegral = \&Math::Prime::Util::PP::LogarithmicIntegral;
    *ExponentialIntegral = \&Math::Prime::Util::PP::ExponentialIntegral;

###############################################################################

foreach my $n (0 .. 1086) {
  if (defined $small_primes{$n}) {
    is( is_prime($n), 2, "$n is prime");
    ok( is_prob_prime($n), "$n is probably prime");
  } else {
    ok(!is_prime($n), "$n is not prime");
    ok(!is_prob_prime($n), "$n is not probably prime");
  }
}

foreach my $n (@primes) {
  is( is_prime($n), 2, "$n is prime" );
  ok( is_prob_prime($n), "$n is probably prime");
}
foreach my $n (@composites) {
  is( is_prime($n), 0, "$n is not prime" );
  is( is_prob_prime($n), 0, "$n is not probably prime");
}

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

is( prev_prime(2), 0, "Previous prime of 2 returns 0" );
is( next_prime(~0-4), 0, "Next prime of ~0-4 returns 0" );

foreach my $n (2010733 .. 2010880) {
  is(next_prime($n), 2010881, "next_prime($n) == 2010881");
}
foreach my $n (2010734 .. 2010881) {
  is(prev_prime($n), 2010733, "prev_prime($n) == 2010733");
}
# Similar test case to 2010870, where m=0 and next_prime is at m=1
is(next_prime(1234567890), 1234567891, "next_prime(1234567890) == 1234567891)");

###############################################################################

while (my($n, $pin) = each (%pivals32)) {
  cmp_ok( prime_count_upper($n), '>=', $pin, "Pi($n) <= upper estimate" );
  cmp_ok( prime_count_lower($n), '<=', $pin, "Pi($n) >= lower estimate" );
  my $approx_range = abs($pin - prime_count_approx($n));
  my $range_limit = ($n <= 100000000) ? 100 : 500;
  cmp_ok( $approx_range, '<=', $range_limit, "prime_count_approx($n) within $range_limit");
}
while (my($n, $pin) = each (%pivals_small)) {
  is( prime_count($n), $pin, "Pi($n) = $pin" );
}

while (my($range, $expect) = each (%pi_intervals)) {
  my($low,$high) = parse_range($range);
  is( prime_count($low,$high), $expect, "prime_count($range) = $expect");
}

###############################################################################

while (my($n, $pin) = each (%pivals_small)) {
  my $next = $pin+1;
  cmp_ok( nth_prime($pin), '<=', $n, "nth_prime($pin) <= $n");
  cmp_ok( nth_prime($next), '>=', $n, "nth_prime($next) >= $n");
}

while (my($n, $nth) = each (%nthprimes32)) {
  cmp_ok( nth_prime_upper($n), '>=', $nth, "nth_prime($n) <= upper estimate" );
  cmp_ok( nth_prime_lower($n), '<=', $nth, "nth_prime($n) >= lower estimate" );

  my $approx = nth_prime_approx($n);
  my $percent_limit = ($n >= 775) ? 1 : 2;
  cmp_ok( abs($nth - $approx) / $nth, '<=', $percent_limit/100.0, "nth_prime_approx($n) = $approx within $percent_limit\% of $nth");
}
while (my($n, $nth) = each (%nthprimes_small)) {
  is( nth_prime($n), $nth, "nth_prime($n) = $nth" );
}

###############################################################################

is( miller_rabin(0, 2), 0, "MR with 0 shortcut composite");
is( miller_rabin(1, 2), 0, "MR with 0 shortcut composite");
is( miller_rabin(2, 2), 2, "MR with 2 shortcut prime");
is( miller_rabin(3, 2), 2, "MR with 3 shortcut prime");

while (my($base, $ppref) = each (%pseudoprimes)) {
  foreach my $p (@$ppref) {
    ok(miller_rabin($p, $base), "Pseudoprime (base $base) $p passes MR");
  }
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

###############################################################################

foreach my $n (@primes) {
  my @f = factor($n);
  is_deeply( \@f, [$n], "factor prime $n yields $n" );
}
foreach my $n (@composites) {
  my @f = factor($n);
  my $facstring = join(' * ', @f);

  # Special case for 0 and 1
  if ($n < 2) {
    cmp_ok( scalar @f, '==', 1, "Factored small $n into itself" );
    is( $f[0], $n, "$n = [ $facstring ]" );
    ok( !is_prime($f[0]), "All factors [ $facstring ] of small $n are not prime" );
    next;
  }

  # These are composites, so they should give us more than one factor
  cmp_ok( scalar @f, '>=', 2, "Factored $n into multiple factors" );

  # Do they multiply to the number?
  my $product = 1;  $product = int($product * $_) for @f;
  is( $product, $n, "$n = [ $facstring ]" );

  # Are they all prime?
  my $isprime = 1; $isprime *= is_prime($_) for @f;
  ok( $isprime, "All factors [ $facstring ] of $n are prime" );
}

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

done_testing();
