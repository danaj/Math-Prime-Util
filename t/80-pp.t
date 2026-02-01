#!/usr/bin/env perl
use strict;
use warnings;

# This is not a comprehensive test for all functions, but very basic
# functionality of the PP code, covering all functions.
#
# For proper testing of the PP code, use:
#
#    MPU_NO_XS=1 MPU_NO_GMP=1 make test
#
# which runs the whole test suite without any C code.
# Add EXTENDED_TESTING=1 if desired.


# Set these first thing, before loading the package.  This will turn off
# both XS and GMP entirely, so everything is the PPFE + PP code.
# A reminder that the caller's versions of these are not changed.
# These are local versions of the environment variables.
BEGIN {
  $ENV{MPU_NO_XS}  = 1;
  $ENV{MPU_NO_GMP} = 1;
}
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


plan tests => 2 +  # require_ok
              1 +  # arithmetic
              1 +  # primality
              1 +  # trial_primes (non-exported function)
              1 +  # primes
              1 +  # sieve_range
              1 +  # next_prime, prev_prime
              1 +  # prime_count
              1 +  # nth_prime
              1 +  # pseudoprimes
              1 +  # omega_primes
              1 +  # almost_primes
              1 +  # prime powers
              1 +  # twin primes
              1 +  # semi primes
              1 +  # ramanujan_primes
              1 +  # real functions
              1 +  # factoring
              1 +  # AKS primality
              1 +  # is_gaussian_prime
              1 +  # is_*_prime
              1 +  # primality proofs
              1 +  # misc ntheory
              1 +  # more misc ntheory
              1 +  # Lucky numbers
              1 +  # perfect powers
              1 +  # powerful
              1 +  # powerfree
              1 +  # set functions
              1 +  # vector (list) functions
              1 +  # rationals
              1 +  # Goldbach
              1;   # $_ is ok

use Math::Prime::Util qw/:all/;
use Math::BigInt;
use Math::BigFloat;
require_ok 'Math::Prime::Util::PP';
require_ok 'Math::Prime::Util::PrimalityProving';

###############################################################################

$_ = 'this should not change';

subtest 'arithmetic ops', sub {

  is(addint(677,24),701,"addint");
  is(subint(677,24),653,"subint");
  is(add1int(677),678,"add1int");
  is(sub1int(677),676,"sub1int");
  is(mulint(677,24),16248,"mulint");
  is("".powint(677,24),"85926683248715705094727267680997536840479271741501353165435057377441","powint");
  is(divint(677,24),28,"divint");
  is(cdivint(677,24),29,"cdivint");
  is(modint(677,24),5,"modint");
  is_deeply([divrem(677,24)],[28,5],"divrem");
  is_deeply([fdivrem(677,24)],[28,5],"fdivrem");
  is_deeply([cdivrem(677,24)],[29,-19],"cdivrem");
  is_deeply([tdivrem(677,24)],[28,5],"tdivrem");
  is(lshiftint(677,3),5416,"lshiftint");
  is(rshiftint(677,3),84,"rshiftint");
  is(rashiftint(677,3),84,"rashiftint");
  is(absint(-677),677,"absint");
  is(negint(677),-677,"negint");
  is_deeply([cmpint(-2,0), cmpint(0,2), cmpint(2,2), cmpint(-7,-7), cmpint(-8,-9)], [-1,-1,0,0,1], "cmpint");
  is_deeply([map { signint($_) } (-7,-1,0,1,7)], [-1,-1,0,1,1], "signint");

  is(sqrtint(677),26,"sqrtint");
  is(rootint(677,3),8,"rootint");
  is(logint(677,2),9,"logint");

  is(negmod(24,9),3,"negmod");
  is(addmod(24,17,9),5,"addmod");
  is(submod(24,170,9),7,"submod");
  is(mulmod(24,170,9),3,"mulmod");
  is(muladdmod(24,170,37,91),22,"muladdmod");
  is(mulsubmod(24,170,37,91),39,"mulsubmod");
  is(powmod(24,170,11),1,"powmod");
  is(divmod(24,29,11),5,"divmod");

  is(invmod(45,59), 21, "invmod(45,59)");
  is(invmod(14,28474), undef, "invmod(14,28474)");
  is(invmod(42,-2017), 1969, "invmod(42,-2017)");

  is(sqrtmod(124,137),undef,"sqrtmod(124,137) = undef");
  is(sqrtmod(11,137),55,"sqrtmod(11,137) = 55");
  is(rootmod(577,3,137),95,"rootmod");
  is_deeply([allsqrtmod(4,13791)],[2,4595,9196,13789],"allsqrtmod");
  is_deeply([allrootmod(581,5,151)],[34,42,43,62,121],"allrootmod");
};

subtest 'primality', sub {
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

  ok(is_prime("18446744073709551521"),"is_prime(18446744073709551521) true")
    if $extra;
};

is_deeply( Math::Prime::Util::PP::trial_primes(80),
           [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79],
           "Trial primes 2-80" );

###############################################################################

subtest 'primes', sub {
  is_deeply( primes(1069), \@small_primes, "Primes between 0 and 1069" );
  is_deeply( primes(1070), \@small_primes, "Primes between 0 and 1070" );
  is_deeply( primes(1086), \@small_primes, "Primes between 0 and 1086" );

  my @small_single = (
    [0,[]], [1,[]], [2,[2]], [3,[2,3]], [4,[2,3]], [5,[2,3,5]], [6,[2,3,5]],
    [7,[2,3,5,7]], [11,[2,3,5,7,11]],
    [18,[2,3,5,7,11,13,17]],
    [19,[2,3,5,7,11,13,17,19]],
    [20,[2,3,5,7,11,13,17,19]]
  );
  foreach my $cinfo (@small_single) {
    my($n,$L) = @$cinfo;
    is_deeply(primes($n),$L,"primes($n) should return [@{$L}]");
  }

  my @small_range = (
    ["3 to 9",[3,5,7]],
    ["2 to 20",[2,3,5,7,11,13,17,19]],
    ["30 to 70",[31,37,41,43,47,53,59,61,67]],
    ["70 to 30",[]],
    ["20 to 2",[]],
    ["1 to 1",[]],
    ["2 to 2",[2]],
    ["3 to 3",[3]],
    ["2 to 3",[2,3]],
    ["2 to 5",[2,3,5]],
    ["3 to 6",[3,5]],
    ["3 to 7",[3,5,7]],
    ["4 to 8",[5,7]],
    ["2010733 to 2010881",[2010733,2010881]],
    ["2010734 to 2010880",[]],
    ["3088 to 3164",[3089,3109,3119,3121,3137,3163]],
    ["3089 to 3163",[3089,3109,3119,3121,3137,3163]],
    ["3090 to 3162",[3109,3119,3121,3137]],
    ["3842610773 to 3842611109",[3842610773,3842611109]],
    ["3842610774 to 3842611108",[]],
  );
  foreach my $cinfo (@small_range) {
    my($rangestr, $L) = @$cinfo;
    my($lo,$hi) = $rangestr =~ /(\d+) to (\d+)/;
    is_deeply(primes($lo,$hi),$L,"primes($lo,$hi) should return [@{$L}]");
  }
};

subtest 'sieve range', sub {
  is_deeply( [sieve_range(4, 4, 1)], [map { $_-4 } 4,5,6,7], "sieve range depth 1" );
  is_deeply( [sieve_range(10, 20, 2)], [1,3,5,7,9,11,13,15,17,19], "sieve range depth 2" );
  is_deeply( [sieve_range(10, 20, 3)], [1,3,7,9,13,15,19], "sieve range depth 3" );
  is_deeply( [sieve_range(10, 20, 5)], [1,3,7,9,13,19], "sieve range depth 5" );
};

###############################################################################

subtest 'next and prev prime', sub {
  foreach my $gap ([19609,52], [360653,96], [2010733,148]) {
    my($base,$range)=@$gap;
    is(next_prime($base), $base+$range, "next prime of $base is $base+$range");
    is(prev_prime($base+$range), $base, "prev prime of $base+$range is $base");
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
    my @samples = (2010733, 2010768, 2010870, 2010880);
    @samples = (2010733..2010880) if $extra;
    is_deeply([map{next_prime($_)}@samples],[map {2010881} 0..$#samples], "next_prime in primegap before 2010881");
  }
  {
    my @samples = (2010734, 2010768, 2010870, 2010881);
    @samples = (2010734..2010881) if $extra;
    is_deeply([map{prev_prime($_)}@samples],[map {2010733} 0..$#samples], "prev_prime in primegap after 2010733");
  }
  # Similar test case to 2010870, where m=0 and next_prime is at m=1
  is(next_prime(1234567890),1234567891,"next_prime(1234567890) == 1234567891)");
  # We were getting this wrong
  is(next_prime("18446744073709551515"),"18446744073709551521","next_prime(18446744073709551515) = 18446744073709551521") if $use64;
};

###############################################################################

subtest 'prime_count', sub {
  my @pivals = ([1,0],[10,4],[100,25],[1000,168],[10000,1229],[60067,6062],[65535,6542]);
  push @pivals, [100000,9592] if $extra;
  for my $pv (@pivals) {
    my($n,$count) = @$pv;
    is(prime_count($n), $count, "prime_count($n) = $count" );
  }

  my @piintervals = (
    ["1e9 +2**14", 785],
    ["17 to 13", 0],
    ["3 to 17", 6],
    ["4 to 17", 5],
    ["4 to 16", 4],
    ["191912783 +248", 2],
    ["191912784 +247", 1],
    ["191912783 +247", 1],
    ["191912784 +246", 0],
  );
  push @piintervals, ["868396 to 9478505",563275], ["1118105 to 9961674",575195], ["24689 to 7973249",535368] if $extra;

  for my $pi (@piintervals) {
    my($rangestr, $count) = @$pi;
    my($lo,$hi) = parse_range($rangestr);
    is(prime_count($lo,$hi), $count, "prime_count($rangestr) = $count");
  }

  # These are small enough they should be exact.
  is( prime_count_lower(450), 87, "prime_count_lower(450)" );
  is( prime_count_upper(450), 87, "prime_count_upper(450)" );
  # Make sure these are about right
  cmp_closeto( prime_count_lower(1234567), 95360, 60, "prime_count_lower(1234567) in range" );
  cmp_closeto( prime_count_upper(1234567), 95360, 60, "prime_count_upper(1234567) in range" );
  cmp_closeto( prime_count_lower(412345678), 21958997, 1500, "prime_count_lower(412345678) in range" );
  cmp_closeto( prime_count_upper(412345678), 21958997, 1500, "prime_count_upper(412345678) in range" );

  my $pca = prime_count_approx(128722248);
  ok( $pca >= 7309252 && $pca <= 7310044, "prime_count_approx(128722248) in range" );
};

###############################################################################

subtest 'nth_prime', sub {
  is(nth_prime(0),undef,"nth_prime(0) returns undef");
  my @nthvals = ([1,2],[4,7],[25,97],[168,997],[1229,9973],[6062,60041],[6542,65521]);
  push @nthvals, [9592,99991] if $extra;
  for my $nv (@nthvals) {
    my($n,$nth) = @$nv;
    is(nth_prime($n),$nth,"nth_prime($n) = $nth");
  }

  my @nthprimes32=(2,29,541,7919,104729,1299709,15485863,179424673,2038074743);
  for my $i (0..$#nthprimes32) {
    my($n, $nth) = (10**$i, $nthprimes32[$i]);
    last if $n > ($extra ? 2000000 : 5000);
    is( nth_prime($n), $nth, "nth_prime($n) = $nth" );
  }

  {
    my $ntha = nth_prime_approx(1287248);
    ok( $ntha >= 20274907 && $ntha <= 20284058, "nth_prime_approx(1287248) in range" );
  }

  {
    my($n,$c) = (15460811,998491);
    my $lo = nth_prime_lower($c);
    my $hi = nth_prime_upper($c);
    my $ap = nth_prime_approx($c);
    my $tol = int($n*.02);

    ok($lo <= $n && $lo+$tol >= $n, "nth_prime_lower($c)");
    ok($hi >= $n && $hi-$tol <= $n, "nth_prime_upper($c)");
    cmp_closeto($ap, $n, $tol, "nth_prime_approx($c)");
  }
};

###############################################################################
subtest 'pseudoprime tests', sub {
  is( is_strong_pseudoprime(0, 2), 0, "MR with 0 shortcut composite");
  is( is_strong_pseudoprime(1, 2), 0, "MR with 0 shortcut composite");
  is( is_strong_pseudoprime(2, 2), 1, "MR with 2 shortcut prime");
  is( is_strong_pseudoprime(3, 2), 1, "MR with 3 shortcut prime");

  my @psp = ([2,341],[2,561],[2,29341],[2,4259905],
             [3, 91],[3,121],[3,44287],[3,4252381]);
  is_deeply([grep { !is_pseudoprime($_->[1],$_->[0]) } @psp],
            [],"Small pseudoprimes");
  my @spsp = ([2,2047],[2,42799],[2,4335241],[2,1078467589],[2,75792980677],
              [3, 121],[3,44287],[3,4252381],[3,1075490821],
              [5, 781],[5,38081],[5,4265257],
              [31,15],[31,49],[31,29341],[31,4270657]);
  is_deeply([grep { !is_strong_pseudoprime($_->[1],$_->[0]) } @spsp],
            [],"Small strong pseudoprimes");
  ok(is_strong_pseudoprime(75792980677),"is_strong_pseudoprime(75792980677)");

  is_deeply([grep {!is_lucas_pseudoprime($_)} qw/9179 10877 44099 4259789/],
            [], "Small Lucas pseudoprimes");
  is_deeply([grep {!is_strong_lucas_pseudoprime($_)} qw/5459 5777 75077 4309631/],
            [], "Small strong Lucas pseudoprimes");
  is_deeply([grep {!is_extra_strong_lucas_pseudoprime($_)} qw/989 3239 5777 72389 4226777 1076503199/],
            [], "Small extra strong Lucas pseudoprimes");
  is_deeply([grep {!is_almost_extra_strong_lucas_pseudoprime($_)} qw/989 3239 5777 72389 4226777/],
            [], "Small AES Lucas pseudoprimes");
  is_deeply([grep {!is_almost_extra_strong_lucas_pseudoprime($_,2)} qw/4531 12209 62479 4403027/],
            [], "Small AES-2 Lucas pseudoprimes");

  is_deeply([grep {is_bpsw_prime($_)} qw/2047 42799 4335241 121 781 989 5777 72389/],
            [], "Small pure BPSW test");

  {
    # n is a SPSP to bases: 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47
    my $n = Math::BigInt->new("168790877523676911809192454171451");
    # We don't need to verify that in this test.
    is( is_strong_pseudoprime( $n, 47), 1, "168790877523676911809192454171451 (SPSP to 2..52) test base 47" );
    is( is_strong_pseudoprime( $n, 53), 0, "168790877523676911809192454171451 found composite with base 53" );
    #is( is_strong_lucas_pseudoprime($n), 0, "368105533664589636823262455488673 is not a strong Lucas pseudoprime" );
  }

  {
    my $n = "153515674455111174527";
    ok(is_extra_strong_lucas_pseudoprime($n), "$n is an ESLSP");
    ok(!is_bpsw_prime($n), "is_bpsw_prime($n) = 0 as expected");
    # Could verify with other tests e.g. frobenius_{khashin,underwood}
  }

  is(is_perrin_pseudoprime(517697641), 1, "517697641 is a Perrin pseudoprime");
  is(is_perrin_pseudoprime(102690901,3), 1, "102690901 is a Perrin pseudoprime (Grantham)");
  is(is_frobenius_pseudoprime(517697641), 0, "517697641 is not a Frobenius pseudoprime");
  is(is_frobenius_khashin_pseudoprime(517697659),1,"517697659 is prime via Frobenius-Khashin test");
  SKIP: {
    # TODO: 2026 does this still happen?
    skip "Old Perl+bigint segfaults in F-U code", 1 if $] < 5.008;
    ok(is_frobenius_underwood_pseudoprime(517697659), "517697659 is prime via Frobenius-Underwood test" );
  }

  is(is_euler_pseudoprime(703, 3), 1, "703 is a base 3 Euler pseudoprime");
  is(is_euler_plumb_pseudoprime(3277), 1, "3277 is a Euler-Plumb pseudoprime");
  is(is_catalan_pseudoprime(17), 1, "is_catalan_pseudoprime(17) true");
  is(is_catalan_pseudoprime(15127), 0, "is_catalan_pseudoprime(15127) false");
  SKIP: {
    # Maybe if we make a faster binomialmod
    skip "Skipping PP Catalan pseudoprime test without EXTENDED_TESTING", 1 unless $extra;
    is(is_catalan_pseudoprime(5907), 1, "5907 is a Catalan pseudoprime");
  }

  is( miller_rabin_random(4294967281, 20), "0", "Miller-Rabin random 40 on composite" );
};

###############################################################################

subtest 'omega primes', sub {
  # omega primes(k), where n is divisible by exactly k primes
  # with multiplicity, meaning these are numbers with prime_omega(n) == k
  is_deeply(omega_primes(1,20),[2,3,4,5,7,8,9,11,13,16,17,19],"omega_primes(1,20)");
  is_deeply(omega_primes(2,20),[6,10,12,14,15,18,20],"omega_primes(2,20)");
  is_deeply(omega_primes(3,100),[30,42,60,66,70,78,84,90],"omega_primes(3,100)");
  is_deeply(omega_primes(4,500),[210,330,390,420,462],"omega_primes(4,500)");

  my @n6 = (8,18,78,510,4620);
  is_deeply([map { nth_omega_prime($_,6) } 1..5],\@n6,"nth_omega_prime(k,6)");
  is_deeply([map { omega_prime_count($_,$n6[$_-1]-1) } 1..5],[5,5,5,5,5],"omega_prime_count(k,n)");

  my @fn7 = (10000019, 10000000, 10000005, 10000002, 10000012, 10000080, 10002930, 11741730, 223092870);
  is_deeply( [map { is_omega_prime($_+1,$fn7[$_]) } 0..$#fn7],
             [(1) x scalar(@fn7)],
             "is_omega_prime (true)" );
  is_deeply( [map { is_omega_prime($_+1,$fn7[$_]+6) } 0..$#fn7],
             [(0) x scalar(@fn7)],
             "is_omega_prime (false)" );
  # random numbers with the result we want
  is_deeply([map { prime_omega($_) } (2,3777,893828,392580,451902,8111460,16265634,9699690,917896980,1084183870770)],[1..10],"prime_omega(n)");

  is_deeply(omega_primes(10,"6469693230","9469693230"),
            [qw/6469693230 6915878970 8254436190 8720021310 9146807670/],
            "omega_primes(10,6469693230,9469693230)") if $extra;
};

subtest 'almost primes', sub {
  # almost primes(k), where n has exactly k prime factors
  # no multiplicity, meaning these are numbers with prime_bigomega(n) == k

  is_deeply(almost_primes(1,20),[2,3,5,7,11,13,17,19],"almost_primes(1,20)");
  is_deeply(almost_primes(2,20),[4,6,9,10,14,15],"almost_primes(2,20)");
  is_deeply(almost_primes(3,20),[8,12,18,20],"almost_primes(3,20)");
  is_deeply(almost_primes(4,60),[16,24,36,40,54,56,60],"almost_primes(4,60)");

  my @n12 = (37,34,52,100,200,400,800,1600,3200,6400,12800,25600);
  is_deeply([map { nth_almost_prime($_,12) } 1..12],\@n12,"nth_almost_prime(k,12)");
  is_deeply([map { almost_prime_count($_,$n12[$_-1]-1) } 1..12],[(11) x 12],"almost_prime_count(k,n)");

  my @isa = (10000019,10000001,10000005,10000002,10000004,10000008,10000016,10000096,10000032,10000080,10000128,10000896);
  is_deeply( [map { is_almost_prime($_,$isa[$_-1]) } 1..12],
             [(1) x scalar(@isa)],
             "is_almost_prime (true)" );
  is_deeply( [map { is_almost_prime($_,$isa[$_-1]+10) } 1..12],
             [(0) x scalar(@isa)],
             "is_almost_prime (false)" );

  is_deeply([map { prime_bigomega($_) } (2,299021,382353,437943,787216,1004848,765264,333882,508640,175872)],[1..10],"prime_bigomega(n)");

  cmp_closeto( almost_prime_count_approx(3,10000),2569, 40, "almost_prime_count_approx(3,10000) in range" );
  # These approximations could be better
  cmp_closeto( almost_prime_count_approx(5,10000), 963, 100, "almost_prime_count_approx(5,10000) in range" );
  cmp_closeto( almost_prime_count_approx(7,"1000000000000"), "62981797962", "10000000000", "almost_prime_count_approx(7,1000000000000) in range" );

  {
    my($k,$n,$c) = (3,389954,98699);
    my $lo = almost_prime_count_lower($k,$n);
    my $hi = almost_prime_count_upper($k,$n);
    my $ap = almost_prime_count_approx($k,$n);
    my $tol = int($c*.05);

    is(almost_prime_count($k,$n),$c,"almost_prime_count($k,$n) = $c");
    ok($lo <= $c && $lo+$tol >= $c, "almost_prime_count_lower($k,$n)");
    ok($hi >= $c && $hi-$tol <= $c, "almost_prime_count_upper($k,$n)");
    cmp_closeto($ap, $c, $tol, "almost_prime_count_approx($k,$n)");
  }
  {
    my($k,$n,$c) = (7,489954,16527);
    my $lo = almost_prime_count_lower($k,$n);
    my $hi = almost_prime_count_upper($k,$n);
    my $ap = almost_prime_count_approx($k,$n);
    my $tol = int($c*.15);

    is(almost_prime_count($k,$n),$c,"almost_prime_count($k,$n) = $c");
    ok($lo <= $c && $lo+$tol >= $c, "almost_prime_count_lower($k,$n)");
    ok($hi >= $c && $hi-$tol <= $c, "almost_prime_count_upper($k,$n)");
    cmp_closeto($ap, $c, $tol, "almost_prime_count_approx($k,$n)");
  }

  {
     my($k,$n,$c) = (4,10000000,47997635);
     my $lo = nth_almost_prime_lower($k,$n);
     my $hi = nth_almost_prime_upper($k,$n);
     my $ap = nth_almost_prime_approx($k,$n);
     my $tol = int(0.1 * $c);
     cmp_closeto( $ap, $c, $tol, "nth_almost_prime_approx($k,$n)");
     ok($lo <= $c && $lo+$tol >= $c, "nth_almost_prime_lower($k,$n)");
     ok($hi >= $c && $hi-$tol <= $c, "nth_almost_prime_upper($k,$n)");
     ok($ap >= $lo && $ap <= $hi, "nth_almost_prime_approx inside lo/hi bounds");
   }
};

###############################################################################

subtest 'prime powers', sub {
  is_deeply([map {$_-100500} @{prime_powers(100500,101000)}],
            [qw/1 11 17 19 23 37 47 49 59 91 109 113 121 149 169 173 193 199 203 233 241 247 269 287 299 301 311 323 329 347 353 407 413 427 431 437 443 457 481 487 499/],
            "prime_powers(100500,101000)");
  is(next_prime_power(13579), 13591, "next_prime_power");
  is(next_prime_power(13591), 13597, "next_prime_power");
  is(prev_prime_power(13579), 13577, "prev_prime_power");
  is(prev_prime_power(13577), 13567, "prev_prime_power");

  {
    my($n,$c) = (389954,33234);
    my $lo = prime_power_count_lower($n);
    my $hi = prime_power_count_upper($n);
    my $ap = prime_power_count_approx($n);
    my $tol = int($c*.05);

    is(prime_power_count($n),$c,"prime_power_count($n) = $c");
    ok($lo <= $c && $lo+$tol >= $c, "prime_power_count_lower($n)");
    ok($hi >= $c && $hi-$tol <= $c, "prime_power_count_upper($n)");
    cmp_closeto($ap, $c, $tol, "prime_power_count_approx($n)");
  }

  is(nth_prime_power(5123),49033,"nth_prime_power(5123) = 49033");
  {
    my($n,$c) = (15460801,999154);
    my $lo = nth_prime_power_lower($c);
    my $hi = nth_prime_power_upper($c);
    my $ap = nth_prime_power_approx($c);
    my $tol = int($n*.05);

    ok($lo <= $n && $lo+$tol >= $n, "nth_prime_power_lower($c)");
    ok($hi >= $n && $hi-$tol <= $n, "nth_prime_power_upper($c)");
    cmp_closeto($ap, $n, $tol, "nth_prime_power_approx($c)");
  }
};

###############################################################################
subtest 'Twin primes', sub {
  is_deeply([map {$_-100500} @{twin_primes(100500,101500)}],
            [qw/17 47 299 611 617 659 707 779/],
            "twin_primes(100500,101500)");

  is(twin_prime_count(4321), 114, "twin_prime_count(4321)");
  is(twin_prime_count(5000,5500), 8, "twin_prime_count(5000,5500)");

  cmp_closeto(twin_prime_count_approx("4123456784123"), "6950213327", 14937 * 2, "twin_prime_count_approx(4123456784123)");
  cmp_closeto("".twin_prime_count_approx("412345678412345678412345678"), "149939117920176008847283", 1e11, "twin_prime_count_approx(412345678412345678412345678)")  if $extra;

  is(nth_twin_prime(249), 13217, "nth_twin_prime(249)");
  cmp_closeto("".nth_twin_prime_approx("1234567890"), "637769466671", 50000000, "nth_twin_prime_approx(1234567890)");
};

###############################################################################
subtest 'Semi primes', sub {
  is_deeply([map {$_-101500} @{semi_primes(101500,101600)}],
            [qw/6 9 21 34 39 46 51 53 54 57 67 71 78 79 89 93 97/],
            "semi_primes(101500,101600)");

  is(semiprime_count(12000,123456),25459,"semiprime_count(12000, 123456)");
  cmp_closeto( semiprime_count_approx("100294967494"), "14000000000", 4000000, "semiprime_count_approx(100294967494) in range" );

  is(nth_semiprime(1400),5137,"nth_semiprime(1400) = 5137");
  cmp_closeto( nth_semiprime_approx("14000000000"), "100294967494", 120000000, "nth_emiprime_approx(14000000000) in range" );
};

###############################################################################
subtest 'Ramanujan primes', sub {
  is_deeply(ramanujan_primes(0,100), [2,11,17,29,41,47,59,67,71,97], "Ramanujan primes under 100");

  {
    my($n,$c) = (8840,500);
    my $lo = ramanujan_prime_count_lower($n);
    my $hi = ramanujan_prime_count_upper($n);
    my $ap = ramanujan_prime_count_approx($n);
    my $tol = int($c*.05);

    is(ramanujan_prime_count($n),$c,"ramanujan_prime_count($n) = $c");
    ok($lo <= $c && $lo+$tol >= $c, "ramanujan_prime_count_lower($n)");
    ok($hi >= $c && $hi-$tol <= $c, "ramanujan_prime_count_upper($n)");
    cmp_closeto($ap, $c, $tol, "ramanujan_prime_count_approx($n)");
  }

  is(nth_ramanujan_prime(28),311,"nth_ramanujan_prime(28) = 311");

  {
    my($n,$c) = (1088761,39999);
    my $lo = nth_ramanujan_prime_lower($c);
    my $hi = nth_ramanujan_prime_upper($c);
    my $ap = nth_ramanujan_prime_approx($c);
    my $tol = int($n*.05);

    #is(nth_ramanujan_prime($c),$n,"nth_ramanujan_prime($c) = $n");
    ok($lo <= $n && $lo+$tol >= $n, "nth_ramanujan_prime_lower($c)");
    ok($hi >= $n && $hi-$tol <= $n, "nth_ramanujan_prime_upper($c)");
    cmp_closeto($ap, $n, $tol, "nth_ramanujan_prime_approx($c)");
  }
};

###############################################################################

subtest 'real (float) functions', sub {
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
  while (my($n, $ein) = each (%eivals)) {
    cmp_closeto( ExponentialIntegral($n), $ein, 0.00000001 * abs($ein), "Ei($n) ~= $ein");
  }

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
  while (my($n, $lin) = each (%livals)) {
    cmp_closeto( LogarithmicIntegral($n), $lin, 0.00000001 * abs($lin), "li($n) ~= $lin");
  }

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
  while (my($n, $rin) = each (%rvals)) {
    cmp_closeto( RiemannR($n), $rin, 0.00000001 * abs($rin), "R($n) ~= $rin");
  }

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
};

###############################################################################
subtest 'factoring', sub {
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
    is_deeply(\@gotfactor,\@expfactor,"test factoring for $ntests composites");
  }

  is_deeply([factor_exp(9147600)],[[2,4],[3,3],[5,2],[7,1],[11,2]],"factor_exp");

  is(join(" ",divisors(252)),"1 2 3 4 6 7 9 12 14 18 21 28 36 42 63 84 126 252","divisors");
  is(divisor_sum(252),728,"divisor_sum(252)");
  is(join(" ",map{divisor_sum(1254,$_)}(0..7)),"16 2880 2208200 2302655040 2659995565256 3210983462174400 3954705863524605800 4916556716966553418560","divisor_sum(1254, {0..7})");
  is(znlog(5678, 5, 10007), 8620, "znlog(5678, 5, 10007)");

  # The PP factor code does small trials, then loops doing 64k rounds of HOLF
  # if the composite is less than a half word, followed by 64k rounds each of
  # prho with a = {3,5,7,11,13}.  Most numbers are handled by these.  The ones
  # that aren't end up being too slow for us to put in a test.  So we'll try
  # running the various factoring methods manually.
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

  ##### factoring subs with BigInts

  # 1013 4294967197 4294967291
  my $nbig = Math::BigInt->new("18686551294184381720251");
  test_facres("prho", $nbig, Math::Prime::Util::PP::prho_factor($nbig));
  test_facres("pbrent", $nbig, Math::Prime::Util::PP::pbrent_factor($nbig));
  test_facres("pminus1", $nbig, Math::Prime::Util::PP::pminus1_factor($nbig));
  SKIP: {
    skip "Skipping ecm test", 1 if defined $Math::Prime::Util::GMP::VERSION && $Math::Prime::Util::GMP::VERSION < 0.20;
    test_facres("ecm", $nbig, Math::Prime::Util::PP::ecm_factor($nbig));
  }

  $nbig = Math::BigInt->new("73786976930493367637");
  # Check stage 2 p-1.  Fast with Math::BigInt::GMP, slow without.
  SKIP: {
    skip "Skipping p-1 stage 2 tests", 1 unless $extra;
    test_facres("pminus1", $nbig, Math::Prime::Util::PP::pminus1_factor($nbig, 27000, 35000));
  }
  test_facres("fermat", $nbig, Math::Prime::Util::PP::fermat_factor($nbig));

  SKIP: {
    skip "Skipping HOLF big test without extended testing", 1 unless $extra;
    test_facres("holf", $nbig, Math::Prime::Util::PP::holf_factor($nbig));
  }
  {
    $nbig = Math::BigInt->new("99999999999979999998975857");
    my @nfac = sort {$a<=>$b} Math::Prime::Util::PP::holf_factor($nbig);
    is_deeply(\@nfac, [9999999998987,10000000001011], "holf correctly factors 99999999999979999998975857");
  }
  SKIP: {
    # Unfortunately we can't guarantee this isn't found in stage 1.
    skip "ecm stage 2", 1 unless $extra;
    $nbig = Math::BigInt->new("14270401808568703916861");
    test_facres("ecm(5,2000)", $nbig, Math::Prime::Util::PP::ecm_factor($nbig, 5, 2000, 40));
  }

  ##### Some numbers that go to stage 2 of tests
  SKIP: {
    skip "stage 2 factoring tests for extended testing", 3 unless $extra;
    my $nbig = Math::BigInt->new("9087500560545072247139");
    my @nfac;
    @nfac = sort {$a<=>$b} Math::Prime::Util::PP::pminus1_factor($nbig,1000,10000);
    is_deeply( [@nfac], ["24133","376559091722747783"], "p-1 stage 2 finds factors of $nbig" );
    @nfac = sort {$a<=>$b} Math::Prime::Util::PP::trial_factor($nbig, 50000);
    is_deeply( [@nfac], ["24133","376559091722747783"], "trial factor finds factors of $nbig" );
    @nfac = sort {$a<=>$b} Math::Prime::Util::PP::ecm_factor($nbig,10,1000,100);
    is_deeply( [@nfac], ["24133","376559091722747783"], "ecm factor finds factors of $nbig" );
  }
};

sub test_facres {
  my($name, $n, @facs) = @_;
  my $eq = cmpint($n, vecprod(@facs)) == 0;
  if (scalar @facs > 1 && $eq && vecnone { $_ == 1 } @facs) {
    pass("$name: $n => [@facs]");
  } else {
    fail("$name: $n => [@facs]");
  }
}


##### AKS primality test.  Be very careful with performance.
subtest 'AKS primality', sub {
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
};

subtest 'is_gaussian_prime', sub {
  ok( !is_gaussian_prime(29,0), "29 is not a Gaussian Prime" );
  ok(  is_gaussian_prime(31,0), "31 is a Gaussian Prime" );
  ok( !is_gaussian_prime(0,-29), "0-29i is not a Gaussian Prime" );
  ok(  is_gaussian_prime(0,-31), "0-31i is a Gaussian Prime" );
  ok(  is_gaussian_prime(58924,132000511), "58924+132000511i is a Gaussian Prime" );
  ok(  is_gaussian_prime(519880,-2265929), "519880-2265929i is a Gaussian Prime" );
  ok( !is_gaussian_prime(20571,150592260), "20571+150592260i is not a Gaussian Prime" );
};

subtest 'other is * prime', sub {
  ok(  is_semiprime(1110000001), "1110000001 is a semiprime" );
  ok( !is_semiprime(1110000201), "1110000201 is not a semiprime" );

  is(is_prime_power("11398895185373143"),19,"is_prime_power(7^19) = 19");
  {my $r; is_prime_power("11398895185373143",\$r); is($r,7,"is_prime_power(7^19,0,r) => r=7");}

  ok(  is_ramanujan_prime(41), "41 is a Ramanujan prime");
  ok( !is_ramanujan_prime(43), "43 is not a Ramanujan prime");

  ok(  is_delicate_prime(294001), "294001 is a delicate prime" );
  ok(  is_delicate_prime(862789,16), "862789 is a delicate prime in base 16" );

  ok(  is_chen_prime(167), "is_chen_prime" );
  is(next_chen_prime(167), 179, "next_chen_prime" );

  ok(  is_mersenne_prime(107), "2^107-1 is a Mersenne prime");
  ok( !is_mersenne_prime(113), "2^113-1 is not a Mersenne prime");

  { my(@got,@exp);
    for my $d ([1,143,0],[1,11,1], [2,313,0],[2,209,1], [3,513,0],[3,1331,1], [4,1331,0],[4,14641,1]) {
      push @got, is_almost_prime($d->[0], $d->[1]);
      push @exp, $d->[2];
    }
    is_deeply(\@got, \@exp, "is_almost_prime");
  }

  my @ipp_2 = (5,347,7080233,17471059,36010357);
  my @ipp_0 = (10,49,697,7080249,17471061,36010359);
  is_deeply([map{is_prob_prime($_)}@ipp_2], [map{2}@ipp_2], "is_prob_prime(p)");
  is_deeply([map{is_prob_prime($_)}@ipp_0], [map{0}@ipp_0], "is_prob_prime(c)");

};

subtest 'primality proofs', sub {
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
};

subtest 'misc number theory functions', sub {
  is( consecutive_integer_lcm(13), 360360, "consecutive_integer_lcm(13)" );
  is( "".consecutive_integer_lcm(52), "3099044504245996706400", "consecutive_integer_lcm(52)" );

  is_deeply( [moebius(513,537)],
             [qw/0 1 1 0 1 -1 1 0 -1 0 -1 0 0 1 1 0 0 -1 0 0 1 -1 1 0 1/],
             "moebius(513,537)" );
  is( moebius(42199), 1, "moebius(42199)" );
  is( liouville(444456), 1, "liouville(444456)" );
  is( liouville(562894), -1, "liouville(562894)" );

  is( mertens(219), 4, "mertens(219)" );
  is( mertens(24219), -67, "mertens(24219)" );

  is_deeply( [euler_phi(1513,1537)],
             [qw/1408 756 800 756 1440 440 1260 576 936 760 1522 504 1200 648 1016 760 1380 384 1530 764 864 696 1224 512 1456/],
             "euler_phi(1513,1537)" );
  is( euler_phi(324234), 108072, "euler_phi(324234)" );
  is( "".jordan_totient(4, 899), "653187225600", "jordan_totient(4, 899)" );
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

  is(binomial(35,16), 4059928950, "binomial(35,16)");
  is("".binomial(228,12), "30689926618143230620", "binomial(228,12)");
  is(binomial(-23,-26), -2300, "binomial(-23,-26) should be -2300");

  is(stirling(12,4,2), '611501', "S(12,4)" );
  is(stirling(12,4,1), '105258076', "s(12,4)" );

  is(join(" ",map{fubini($_)}0..6,18),"1 1 3 13 75 541 4683 3385534663256845323","fubini(n) for n in {0..6,18}");

  is_deeply([numtoperm(11,33967658)],[9,3,6,4,7,1,10,0,5,2,8],"numtoperm");
  is(permtonum([9,3,6,4,7,1,10,0,5,2,8]),33967658,"permtonum");

  is_deeply([map{[bernfrac($_)]}(0,1,2,3,12,13)],
            [[1,1],[1,2],[1,6],[0,1],[-691,2730],[0,1]],
            "bernfrac");
  cmp_closeto(bernreal(18),54.971177944862155,1e-5,"bernreal");

  is_deeply([map{[harmfrac($_)]}(0,1,2,3,12,13)],
            [[0,1],[1,1],[3,2],[11,6],[86021,27720],[1145993,360360]],
            "harmfrac");
  cmp_closeto(harmreal(18),3.49510807819631349,1e-5,"harmreal");

  is_deeply( [gcdext(23948236,3498248)], [2263, -15492, 52], "gcdext(23948236,3498248)" );

  is( valuation(1879048192,2), 28, "valuation(1879048192,2)");
  is( valuation(96552,6), 3, "valuation(96552,6)");

  cmp_closeto( chebyshev_theta(7001), 6929.27483821865062, 0.006929, "chebyshev_theta(7001) =~ 6929.2748");
  cmp_closeto( chebyshev_psi(6588), 6597.07452996633704, 0.006597, "chebyshev_psi(6588) =~ 6597.07453");

  is(inverse_totient(42), 4, "inverse totient 42 count");
  is_deeply([inverse_totient(42)], [43,49,86,98], "inverse totient 42 list");

  is( primorial(24), 223092870, "primorial(24)" );
  is( "".primorial(118), "31610054640417607788145206291543662493274686990", "primorial(118)" );
  is( pn_primorial(7), 510510, "pn_primorial(7)" );
  is( partitions(74), 7089500, "partitions(74)" );

  is(legendre_phi(54321,5),11287,"legendre_phi(54321,5) = 11287");
  is(inverse_li(13579),146261,"inverse_li");
  cmp_closeto(inverse_li_nv(135790),1808203.25662372,1e-4,"inverse_li_nv");

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
  { my @t;
    forsemiprimes(sub {push @t,$_}, 152026,152060);
    is_deeply( [@t], [152049,152051,152059], "forsemiprimes 152026,152060" );
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

  { my @set=([1,2],[qw/a b c/]);  my @out;
    forsetproduct {push @out,"@_"} @set;
    is_deeply(\@out, ['1 a','1 b','1 c','2 a','2 b','2 c'], "forsetproduct([1,2],[qw/a b c/])" );
  }

  is( Pi(82), "3.141592653589793238462643383279502884197169399375105820974944592307816406286208999", "Pi(82)" );

  is( gcd(-30,-90,90), 30, "gcd(-30,-90,90) = 30" );
  is( lcm(11926,78001,2211), 2790719778, "lcm(11926,78001,2211) = 2790719778" );

  is(sum_primes(14400),11297213,"sum_primes(14400)");
  is(sum_primes(2100000),"156999759090","sum_primes(2100000)") if $extra;

  is(mertens(5443),9,"mertens(5443)");
  is(sumtotient(5443),9008408,"sumtotient(5443)");
  is(sumliouville(5443),-21,"sumliouville(5443)");
  is_deeply([map {powersum(5443,$_)} 1..8],[qw/14815846 53766705134 219509292695716 955919057077963010 4336287761695106589076 20232498884989465784893754 96368654823259273645222916236 466295787313885438803232358983490/],"powersum");

  is(sumdigits("0b10101110101"),7,"sumdigits with binary string");
  is(sumdigits(5443),16,"sumdigits with integer");
  is(sumdigits("4def",16),46,"sumdigits with hex");
  is(sumdigits("x4ldef",36),100,"sumdigits with base 36");

  is(hammingweight(5443),6,"hammingweight");

  is_deeply([kronecker(18,106),kronecker(19,106),kronecker(11,106)],[0,1,-1],"kronecker");
  is_deeply([[cornacchia(17,131)],[cornacchia(2,131)],[cornacchia(2,136)]],[[undef],[9,5],[8,6]],"cornacchia");
  is(hclassno(320),168,"hclassno");
  is_deeply([ramanujan_tau(81),ramanujan_tau(41),ramanujan_tau(44)],[1665188361,308120442,-786948864],"ramanujan_tau");

  is(lucasu(6,1,14),9228778026,"lucasu");
  is(lucasv(6,1,14),52205852194,"lucasv");
  is(lucasumod(1,-1,281,17779),5050,"lucasumod");
  is(lucasvmod(1,-1,281,17779),8665,"lucasvmod");
  is_deeply([lucasuvmod(1,-1,2811,17779)],[6323,16441],"lucasuvmod");
  is(pisano_period(1777),3556,"pisano_period");
};

subtest 'more misc ntheory functions', sub {
  ok( is_totient(381554124), "381554124 is a totient");
  ok(!is_totient(1073024875), "1073024875 is not a totient");

  ok(!is_carmichael(5049), "5049 is not a Carmichael number");
  ok(!is_carmichael(2792834247), "2792834247 is not a Carmichael number");
  ok(!is_carmichael(2399550475), "2399550475 is not a Carmichael number");
  ok(!is_carmichael(219389), "219389 is not a Carmichael number");
  ok( is_carmichael(1125038377), "1125038377 is a Carmichael number");

  ok( is_quasi_carmichael(1517), "1517 is quasi-Carmichael");
  ok( is_quasi_carmichael(10001), "10001 is quasi-Carmichael");
  ok( is_quasi_carmichael(10373), "10373 is quasi-Carmichael");

  ok(!is_cyclic(1521), "1521 is not cyclic");
  ok( is_cyclic(10001), "10001 is cyclic");
  ok( is_pillai(26657), "26657 is Pillai");
  ok(!is_practical(1701), "1701 is not practical");
  ok( is_practical(1710), "1710 is practical");
  ok( is_fundamental(-168), "-168 is fundamental");
  ok( is_fundamental(172), "172 is fundamental");
  ok( is_congruent_number(692), "692 is a congruent number");
  ok( is_congruent_number(206), "206 is a congruent number");
  ok( is_congruent_number(207), "207 is a congruent number");

  is(is_happy(536), 7, "536 is a happy number");
  is(is_happy(571,7), 3, "571 is a happy number in base 7");
  is(is_happy(347,6), 4, "347 is a happy number in base 6");
  is(is_happy(514,16,3), 3, "514 is a happy number in base 16 with exponent 3");

  is(ramanujan_sum(12,36), 4, "ramanujan_sum(12,36) = 4");

  is(is_power(16926659444736),17,"is_power(6^17) = 17");
  {my $r; is_power(16926659444736,0,\$r); is($r,6,"is_power(6^17,0,r) => r=6");}

  is(is_square(603729), 1, "603729 is a square");

  is_deeply( [map { is_sum_of_squares($_) } (-10 .. 10, 437)],
             [1,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,0,1,1,1,0],
             "is_sum_of_squares (k=2) for -10 .. 10, 437" );

  ok( is_polygonal(6,3), "6 is a 3-polygonal number" );
  ok( is_polygonal(9,4), "9 is a 4-polygonal number" );

  ok(!is_odd(576),"is_odd(576)");
  ok( is_odd(577),"is_odd(577)");
  ok( is_even(576),"is_even(576)");
  ok(!is_even(577),"is_even(577)");
  ok(!is_divisible(30,7),"is_divisible(30,7)");
  ok( is_divisible(30,5),"is_divisible(30,5)");
  ok(!is_congruent(100007,176,177),"is_congruent(100007,176,177)");
  ok( is_congruent(100007,2,177),"is_congruent(100007,2,177)");

  ok(!is_square_free(331483),"is_square_free(331483)");
  ok( is_square_free(370481),"is_square_free(370481)");
  ok(!is_primitive_root(3,1777),"is_primitive_root(3,1777)");
  ok( is_primitive_root(5,1777),"is_primitive_root(5,1777)");
  ok(!is_perfect_number(2048),"is_perfect_number(2048)");
  ok( is_perfect_number(8128),"is_perfect_number(8128)");

  is(fromdigits([1,1,0,1],2), 13, "fromdigits binary");
  is(fromdigits([0,1,1,0,1],16), 4353, "fromdigits base 16");
  is_deeply([todigits(77)], [7,7], "todigits 77");
  is_deeply([todigits(77,2)], [1,0,0,1,1,0,1], "todigits 77 base 2");
  is(todigitstring(-143,16), "8f", "todigitstring base 16");

  is(tozeckendorf(1025),"100000010000101","tozeckendorf");
  is(fromzeckendorf("100000010000101"),1025,"fromzeckendorf");

  ok(!is_qr(177,10256), "177 is not a quadratic residue mod 10256");
  ok( is_qr(180,10256), "180 is a quadratic residue mod 10256");
  is(qnr(10271),7,"qnr(10271) = 7");

  is(chinese([14,643], [254,419], [87,733]), 87041638, "chinese");
  is_deeply([chinese2([14,643], [254,419], [87,733])], [87041638,197482661], "chinese2");

  is(frobenius_number(5,13,29),37,"frobenius_number");

  is("".factorial(53),"4274883284060025564298013753389399649690343788366813724672000000000000","factorial(53)");
  is(factorialmod(53,177),30,"factorialmod(53,177)");
  is(factorialmod(830,1777),1771,"factorialmod(830,1777)")  if $extra;

  is(subfactorial(15),481066515734,"subfactorial(15)");

  is(binomialmod(53,7,177),152,"binomialmod");

  is(falling_factorial(17,5),742560,"falling_factorial");
  is(rising_factorial(17,5),2441880,"rising_factorial");

  # "A k-rough number, as defined by Finch in 2001 and 2003, is a positive
  #  integer whose prime factors are all greater than or equal to k."
  ok( is_rough("62000279000279",31),"is_rough(31*n,31) = 1");
  ok(!is_rough("62000279000279",32),"is_rough(31*n,32) = 0");
  # "[An n-smooth number] is an integer whose prime factors are all
  #  less than or equal to n."
  ok( is_smooth(1291677,50),"is_smooth(1291677,50) = 1");
  ok( is_smooth(1291677,43),"is_smooth(1291677,43) = 1");
  ok(!is_smooth(1291677,42),"is_smooth(1291677,42) = 0");

  is(smooth_count(1291677,43),32842,"smooth_count");
  is(rough_count(1291677,43),187389,"rough_count");
};


subtest 'Lucky numbers', sub {
  ok(!is_lucky(1772),"1772 is not a lucky number");
  ok( is_lucky(1771),"1771 is a lucky number");

  is_deeply(lucky_numbers(600,700), [map {600+$_}1,13,15,19,21,31,39,43,45,51,55,73,79,85,93,99], "lucky numbers between 600 and 700");

  {
    my($n,$c) = (8840,1004);
    my $lo = lucky_count_lower($n);
    my $hi = lucky_count_upper($n);
    my $ap = lucky_count_approx($n);
    my $tol = int($c*.05);

    is(lucky_count($n),$c,"lucky_count($n) = $c");
    ok($lo <= $c && $lo+$tol >= $c, "lucky_count_lower($n)");
    ok($hi >= $c && $hi-$tol <= $c, "lucky_count_upper($n)");
    cmp_closeto($ap, $c, $tol, "lucky_count_approx($n)");
  }

  is(nth_lucky(28),129,"nth_lucky(28) = 129");

  {
    my($n,$c) = (18605821,1088761);
    my $lo = nth_lucky_lower($c);
    my $hi = nth_lucky_upper($c);
    my $ap = nth_lucky_approx($c);
    my $tol = int($n*.05);

    ok($lo <= $n && $lo+$tol >= $n, "nth_lucky_lower($c) $lo <= $n");
    ok($hi >= $n && $hi-$tol <= $n, "nth_lucky_upper($c) $hi >= $n");
    cmp_closeto($ap, $n, $tol, "nth_lucky_approx($c) $ap =~ $n");
  }
};


subtest 'perfect powers', sub {
  is(is_perfect_power(19487172),0,"19487172 is not a perfect power");
  is(is_perfect_power(19487171),1,"19487171 is a perfect power");
  is(next_perfect_power(161051),161604,"next_perfect_power(5^7) = 402^2");
  is(prev_perfect_power(161051),160801,"prev_perfect_power(5^7) = 401^2");

  is(perfect_power_count(123456),404,"perfect_power_count(123456) = 404");
  is(perfect_power_count(123456,234567),148,"perfect_power_count(123456,234567) = 148");
  ok(perfect_power_count_lower("9999999999900000000000000") <= 3162493192548,"perfect_power_count_lower");
  ok(perfect_power_count_upper("9999999999900000000000000") >= 3162493192548,"perfect_power_count_upper");
  cmp_closeto(perfect_power_count_approx("9999999999900000000000000"),3162493192548,1000,"perfect_power_count_approx");

  #is("".nth_perfect_power(1234567890),"1521310467887050801","nth_perfect_power");
  is(nth_perfect_power(1234567),1495530880561,"nth_perfect_power");
  ok(nth_perfect_power_lower(1234567) <= 1495530880561,"nth_perfect_power_lower");
  ok(nth_perfect_power_upper(1234567) >= 1495530880561,"nth_perfect_power_lower");
  cmp_closeto(nth_perfect_power_approx(1234567),1495530880561,10000000,"nth_perfect_power_approx");
};

subtest 'powerful', sub {
  ok(!is_powerful(260),"260 is not a powerful number");
  ok( is_powerful(243),"243 is a powerful number");
  ok( is_powerful("2011901648110693",3),"157^3 * 151^4 is a 3-powerful number");

  is_deeply(powerful_numbers(10500,11000),[10584,10609,10648,10800,10816,10952,10976],"powerful_numbers(10500,11000)");
  is(powerful_count(1234567),2255,"powerful_count(1234567)");
  is(powerful_count(1234567,3),329,"powerful_count(1234567,3)");
  is(nth_powerful(1000),253472,"nth_powerful");

  is_deeply([map {sumpowerful(5443,$_)} 1..8],[14815846,262303,66879,30528,14445,11045,10252,7937],"sumpowerful");
};

subtest 'powerfree', sub {
  is(is_powerfree(1000),0,"1000 is not powerfree");
  is(is_powerfree(1001),1,"1001 is powerfree");

  is(powerfree_count(10500),6385,"powerfree_count");
  is(nth_powerfree(10500),17266,"nth_powerfree");

  is(powerfree_part(100040),25010,"powerfree_part(100040) = 25010");
  is(powerfree_part(100040,3),12505,"powerfree_part(100040,3) = 12505");
  is(squarefree_kernel(100040),25010,"squarefree_kernel(100040) = 25010");

  is(powerfree_part(10004),2501,"powerfree_part(10040) = 2501");
  is(squarefree_kernel(10004),5002,"squarefree_kernel(10004) = 5002");

  is_deeply([map {powerfree_sum(5443,$_)} 1..8],[1,8999622,12322494,13687065,14286122,14561514,14693701,14756710],"powerfree_sum");
  is(powerfree_part_sum(100040),3292589515,"powerfree_part_sum(100040)");
  is(powerfree_part_sum(100040,3),4234954627,"powerfree_part_sum(100040,3)");
  is(powerfree_part_sum(100040,4),4642253940,"powerfree_part_sum(100040,4)");

  is(powerfree_count(27000000,3),22461494,"powerfree_count(27000000,3)");
  is(powerfree_count(400040001,2),243195224,"powerfree_count(400040001,2)");
  is(powerfree_count("10000000000",6),9829525925,"powerfree_count(10000000000,6)");
  #is("".powerfree_count("27000000000000",3),"22461499059723","powerfree_count(30000^3,3)");
  is("".powerfree_count("100000000000000000000",15),"99996941269930456119","powerfree_count(10^20,15)");
};

###############################################################################
subtest 'set functions', sub {
  my @OS = (-5,0,1,2,8,17,20);
  is_deeply([toset([-5,17,2,8,2,0,20,1,2])],\@OS,"toset");
  { my @S = @OS;  my $r = setinsert(\@S,1);
    is_deeply([$r,@S],[0,-5,0,1,2,8,17,20],"setinsert one element already in set");
  }
  { my @S = @OS;  my $r = setinsert(\@S,3);
    is_deeply([$r,@S],[1,-5,0,1,2,3,8,17,20],"setinsert one element not in set");
  }
  { my @S = @OS;  my $r = setinsert(\@S,[1,3,18,21]);
    is_deeply([$r,@S],[3,-5,0,1,2,3,8,17,18,20,21],"setinsert 4 elements, one in set");
  }
  { my @S = @OS;  my $r = setremove(\@S,4);
    is_deeply([$r,@S],[0,-5,0,1,2,8,17,20],"setremove 1 element not in set");
  }
  { my @S = @OS;  my $r = setremove(\@S,1);
    is_deeply([$r,@S],[1,-5,0,2,8,17,20],"setremove 1 element in set");
  }
  { my @S = @OS;  my $r = setremove(\@S,[1,2,3]);
    is_deeply([$r,@S],[2,-5,0,8,17,20],"setremove 2 elements in set");
  }
  { my @S = @OS;  my $r = setinvert(\@S,1);
    is_deeply([$r,@S],[-1,-5,0,2,8,17,20],"setinvert");
  }
  { my @S = @OS;  my $r = setinvert(\@S,3);
    is_deeply([$r,@S],[1,-5,0,1,2,3,8,17,20],"setinvert");
  }
  { my @S = @OS;  my $r = setinvert(\@S,[2,3]);
    is_deeply([$r,@S],[0,-5,0,1,3,8,17,20],"setinvert");
  }
  is(setcontains(\@OS,3),0,"setcontains one not found");
  is(setcontains(\@OS,1),1,"setcontains one found");
  is(setcontains(\@OS,[-1,8]),0,"setcontains subset not found");
  is(setcontains(\@OS,[-5,0,8]),1,"setcontains subset found");
  is(setcontainsany(\@OS,[-1,7]),0,"setcontainsany not found");
  is(setcontainsany(\@OS,[-1,8]),1,"setcontainsany found");

  is(setcontains(primes(500),353),1,"setcontains with primes");

  is_deeply([setbinop { $a * $b } [1,2,3],[2,3,4]],[2,3,4,6,8,9,12],"setbinop");
  is_deeply([sumset([-1,0,1])], [-2,-1,0,1,2], "sumset");

  {
    my @A = (1,2,3,4);
    my @B = (3,4,5,6);
    is_deeply([setunion(\@A,\@B)],[1,2,3,4,5,6],"setunion");
    is_deeply([setintersect(\@A,\@B)],[3,4],"setintersect");
    is_deeply([setminus(\@A,\@B)],[1,2],"setminus");
    is_deeply([setminus(\@B,\@A)],[5,6],"setminus");
    is_deeply([setdelta(\@A,\@B)],[1,2,5,6],"setdelta");

    is(is_sidon_set(\@A),0,"is_sidon_set (false)");
    is(is_sidon_set([1,2,4,8,16,29,58]),1,"is_sidon_set (true)");
    is(is_sidon_set([0,1,4,6]),1,"is_sidon_set (true)");
    is(is_sumfree_set(\@B),0,"is_sumfree_set (false)");
    is(is_sumfree_set([2,3,15]),1,"is_sumfree_set (true)");

    is(set_is_disjoint(\@A,\@B),0,"set_is_disjoint");
    is(set_is_disjoint(\@A,[6,7]),1,"set_is_disjoint");
    is(set_is_equal(\@A,\@B),0,"set_is_equal");
    is(set_is_equal(\@A,\@A),1,"set_is_equal");
    is(set_is_subset(\@A,[5,6]),0,"set_is_subset");
    is(set_is_subset(\@A,[3,4]),1,"set_is_subset");
    is(set_is_proper_subset(\@A,[2,3,4]),1,"set_is_proper_subset");
    is(set_is_proper_subset(\@A,[1,2,3,4]),0,"set_is_proper_subset");
    is(set_is_superset(\@A,\@B),0,"set_is_superset");
    is(set_is_superset(\@A,[1..8]),1,"set_is_superset");
    is(set_is_proper_superset(\@A,[1..8]),1,"set_is_proper_superset");
    is(set_is_superset(\@A,\@A),1,"set_is_superset");
    is(set_is_proper_superset(\@A,\@A),0,"set_is_proper_superset");
    is(set_is_proper_intersection(\@A,\@A),0,"set_is_proper_intersection");
    is(set_is_proper_intersection(\@A,\@B),1,"set_is_proper_intersection");
  }
};

###############################################################################
subtest 'vector (list) functions', sub {
  is(vecsum(15, 30, 45), 90, "vecsum(15,30,45)");
  is("".vecsum(4294966296,4294965296,4294964296), "12884895888", "vecsum(2^32-1000,2^32-2000,2^32-3000)");
  is(vecprod(15, 30, 45), 20250, "vecprod(15,30,45)");
  is("".vecprod(4294966296,4294965296,4294964296), "79228051833847139970490254336", "vecprod(2^32-1000,2^32-2000,2^32-3000)");
  is(vecmin(4294966296,4294965296,4294964296), 4294964296, "vecmin(2^32-1000,2^32-2000,2^32-3000)");
  is(vecmax(4294966296,4294965296,4294964296), 4294966296, "vecmax(2^32-1000,2^32-2000,2^32-3000)");

  is(vecmin(2,-2,7,-1,5,-3,400,0),-3,"vecmin");
  is(vecmax(2,-2,7,-1,5,-3,400,0),400,"vecmax");
  is(vecsum(2,-2,7,-1,5,-3,400,0),408,"vecsum");
  is(vecprod(15,30,4,3),5400,"vecprod");
  is(vecreduce(sub{$a+$b},(5,6,-9,12)),14,"vecreduce");
  is(join("", vecextract(['a'..'z'],[15,4,17,11])), "perl", "vecextract");
  is(vecequal([1,2,3],[1,2,3]), 1, "vecequal([1,2,3],[1,2,3]) = 1");
  is(vecequal([1,2,3],[3,2,1]), 0, "vecequal([1,2,3],[3,2,1]) = 0");
  ok(  (vecany { $_ == 1 } 1, 2, 3), 'vecany true' );
  ok( !(vecall { $_ == 1 } 1, 2, 3), 'vecall false' );
  ok(  (vecnotall { $_ == 1 } 1, 2, 3), 'vecnotall true' );
  ok(  (vecnone { $_ == 1 } 2, 3, 4), 'vecnone true' );

  is_deeply([vecsort(3,-1,3,0,1,3,-5,4,1,-3)],[-5,-3,-1,0,1,1,3,3,3,4],"vecsort");
  { my @L=(0,-4,3,4,-1,-4,4,-2,2,3);
    vecsorti(\@L);
    is_deeply(\@L,[-4,-4,-2,-1,0,2,3,3,4,4],"vecsorti");
  }

  is(vecfirst(sub{$_>6},(3,6,-7,17,7,8,9)),17,"vecfirst");
  is(vecfirstidx(sub{$_>6},(3,6,-7,17,7,8,9)),3,"vecfirstidx");
  is_deeply([vecuniq(1,1,3,2,4,0,4,0,3,1)],[1,3,2,4,0],"vecuniq");
  is_deeply([vecsingleton(1,1,3,2,4,0,4,0,3,1)],[2],"vecsingleton");
  {
    my @L = (-1,14,4,-4,2,2,3,4,3,4,4,1);
    my %got = vecfreq(@L);
    my %exp = (-1=>1, 14=>1, 4=>4, -4=>1, 2=>2, 3=>2, 1=>1);
    is_deeply(\%got, \%exp, "vecfreq");
  }
  is(vecmex(0,1,2,4), 3, "vecmex(0,1,2,4) = 3");
  is(vecpmex(1,2,24,5), 3, "vecpmex(1,2,24,5) = 3");
  {
    my @d128 = (1..128);
    my @s128 = shuffle(@d128);
    my @t128 = sort {$a<=>$b} @s128;
    is(scalar @s128, scalar @d128, "shuffle n items returns n items");
    isnt("@s128","@d128", "shuffled 128-element array isn't identical");
    is("@t128","@d128", "shuffled outputs are the same elements as input");
  }
  is(scalar @{[vecsample(4,[8..11])]}, 4, "vecsample returns all items with exact k");

  is_deeply([vecslide {$a+$b} 1..5],[3,5,7,9],"vecslide {\$a+\$b} 1..5");
};

###############################################################################
subtest 'rationals', sub {
  is_deeply([contfrac(25999,17791)],[1,2,5,1,31,1,2,1,4,2],"contfrac");
  is_deeply([from_contfrac(1,2,5,1,31,1,2,1,4,2)],[25999,17791],"from_contfrac");

  # Both the Perl and Python code on Rosettacode (Jan 2026) is wrong for these.

  is(calkin_wilf_n(1249,9469), 10000000, "calkin_wilf_n(1249,9469) = 10000000");
  is_deeply([nth_calkin_wilf(10000000)],[1249,9469],"nth_calkin_wilf(10000000)");

  is_deeply([next_calkin_wilf(25999,17791)],[17791,27374],"next_calkin_wilf");
  is("".calkin_wilf_n(25999,17791), "834529325481721", "calkin_wilf_n(25999,17791)");
  is_deeply([nth_calkin_wilf("834529325481721")],[25999,17791],"nth_calkin_wilf(834529325481721)");

  is(stern_brocot_n(1249,9469), 8434828, "stern_brocot_n(1249,9469) = 8434828");
  is_deeply([nth_stern_brocot(8434828)],[1249,9469],"nth_stern_brocot(1249,9469)");
  is_deeply([next_stern_brocot(1249,9469)],[1409,10682],"next_stern_brocot");
  is(stern_brocot_n(1409,10682), 8434829, "stern_brocot_n(1409,10682)");

  is(nth_stern_diatomic(10000000), 1249, "nth_stern_diatomic");

  is_deeply([farey(6)],[[0,1],[1,6],[1,5],[1,4],[1,3],[2,5],[1,2],[3,5],[2,3],[3,4],[4,5],[5,6],[1,1]],"farey(6)");
  is_deeply(farey(144,146),[3,125],"farey(144,146)");
  is(scalar farey(1445), 635141, "scalar farey(1445) = 635141");
  is_deeply(next_farey(188,[3,5]),[113,188],"next_farey");
  is_deeply([farey_rank(188,[3,5]),farey_rank(188,[113,188])],[6478,6479],"farey_rank");
};

###############################################################################
subtest 'Goldbach', sub {
  is(minimal_goldbach_pair(258),7,"minimal_goldbach_pair");
  is(goldbach_pair_count(4620),190,"goldbach_pair_count");
  is_deeply([goldbach_pairs(180)],[7,13,17,23,29,31,41,43,53,67,71,73,79,83],"goldbach_pairs");
  is_deeply([goldbach_pairs(175)],[2],"goldbach_pairs for odd n where n-2 is prime");
  is_deeply([goldbach_pairs(177)],[],"goldbach_pairs for odd n where n-2 is not prime");
  is_deeply([map { scalar goldbach_pairs($_) } (180,175,177)], [14,1,0], "scalar goldbach_pairs returns count");
};

# Not here:
#  is_provable_prime
#  is_provable_prime_with_cert
#  prime_certificate
#  verify_prime
#
#  prime_precalc
#  prime_memfree
#  prime_get_config
#  prime_set_config
#
#  print_primes
#
#  irand
#  irand64
#  drand
#  random_bytes
#  entropy_bytes
#  urandomb
#  urandomm
#  csrand
#  rand
#  random_factored_integer
#  randperm
#
#  random_prime
#  random_ndigit_prime
#  random_nbit_prime
#  random_safe_prime
#  random_strong_prime
#  random_proven_prime
#  random_maurer_prime
#  random_shawe_taylor_prime
#  random_unrestricted_semiprime
#  random_semiprime


# foralmostprimes {...} k,[beg,],end  loop over k-almost-primes in range
# forsquarefree {...} [start,] end    loop with factors of square-free n
# forsquarefreeint {...} [start,] end loop over square-free n
# forcomp { ... } n [,{...}]          loop over integer compositions
# formultiperm { ... } \@n            loop over multiset permutations
# forderange { ... } n                loop over derangements
# prime_iterator                      returns a simple prime iterator
# prime_iterator_object               returns a prime iterator object
# lastfor                             stop iteration of for.... loop

# TODO:
#
#  sieve_prime_cluster
#
#  znorder bigint
#  non-GMP lucas_sequence
#  pminus1_factor stage 2
#  forcompositions with hash

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
