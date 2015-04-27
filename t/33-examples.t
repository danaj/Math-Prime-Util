#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use Math::Prime::Util qw/:all/;
use Math::Prime::Util::PrimeArray;
use List::Util qw/first/;

# Make sure things used as examples in the documentation work.

BEGIN {
  unless ($ENV{RELEASE_TESTING}) {
    plan( skip_all => 'these tests are for release candidate testing' );
  }
}

plan tests => 99;

{
  my @nums;
  forprimes { push @nums, $_  if is_prime($_+2) } 10000;
  is(scalar @nums, twin_prime_count(10000), "twin primes via forprimes");
}
{
  my @nums;
  forcomposites { push @nums,$_  if is_strong_pseudoprime($_,2) } 10000, 10**6;
  is_deeply(\@nums, [qw/15841 29341 42799 49141 52633 65281 74665 80581 85489 88357 90751 104653 130561 196093 220729 233017 252601 253241 256999 271951 280601 314821 357761 390937 458989 476971 486737 489997 514447 580337 635401 647089 741751 800605 818201 838861 873181 877099 916327 976873 983401/], "spsp-2s in range using forcomposites");
}
is( prime_count( 1_000_000 ), 78498, "prime_count(1M)" );
is( prime_count( 10**14, 10**14+1000 ), 30, "prime_count(10^14,10^14+1000)" );
{
  my $n = "100000000000000000";   # 10^17
  my $approx = prime_count_approx($n);
  cmp_ok($approx, '>=', prime_count_lower($n), "10^17: Pi approx >= Pi lower");
  cmp_ok($approx, '<=', prime_count_upper($n), "10^17: Pi approx <= Pi upper");
  cmp_closeto($approx, 2623557157654233, 0.000001, "10^17: Pi approx within .0001%");
}
is(nth_prime(10000), 104729, "nth_prime(10000)");
{
  my $n = "1000000000000";   # 10^12
  my $approx = nth_prime_approx($n);
  cmp_ok($approx, '>=', nth_prime_lower($n), "10^17: nth approx >= nth lower");
  cmp_ok($approx, '<=', nth_prime_upper($n), "10^17: nth approx <= nth upper");
  cmp_closeto($approx, 29996224275833, 1e-5, "10^12: nth approx within .001%");
}
is(euler_phi("801294088771394680000412"), "391329671260448564651280", "euler_phi(801294088771394680000412)");
is(jordan_totient(5,1234), "2771963542268536", "jordan_totient(5,1234)");
{
  my $sum = 0;  $sum += moebius($_) for 1..200;
  is($sum, -8, "Mertens(200) via moebius");
}
is(mertens(10_000_000), 1037, "Mertens(10_000_000)");
is(exp_mangoldt(49), 7, "exp_mangoldt(49)");
is(liouville(4292384), -1, "liouville(4292384)");
cmp_closeto(chebyshev_psi(234984), 235070.385453159, 1e-6, "chebyshev_psi(234984)");
cmp_closeto(chebyshev_theta(92384234), 92371752.9943251, 1e-6, "chebyshev_theta(92384234)");
is(partitions(1000), "24061467864032622473692149727991", "partitions(1000)");
{
  my($nparts,$nels) = (0,0);
  forpart { do { $nparts++; $nels += scalar @_; } unless scalar grep { !is_prime($_) } @_ } 25;
  is($nparts, 52, "partions of 25 with all prime elements: 52 found");
  is($nels, 333, "partions of 25 with all prime elements: 333 total values");
}
is(primorial(47), "614889782588491410", "primorial(47)");
is(pn_primorial(47), "1645783550795210387735581011435590727981167322669649249414629852197255934130751870910", "pn_primorial(47)");

##############################################################################

{
  my $aref = primes( 1_000_000_000_000, 1_000_000_001_000 );
  my $eref = [map { "1000000000".$_ } qw/039 061 063 091 121 163 169 177 189 193 211 271 303 331 333 339 459 471 537 543 547 561 609 661 669 721 751 787 789 799 841 903 921 931 933 949 997/];
  is_deeply($aref,$eref,"primes(1000M,1000M+1000)");
}

{
  my @nums;
  forprimes { push @nums, $_ } 100,200;
  is_deeply(\@nums, primes(100,200), "forprimes 100,200");
}
{
  my $sum = 0;  forprimes { $sum += $_ } 100000;
  is($sum, 454396537, "forprimes sum primes to 100k");
}
{
  my @ecomp = grep { !is_prime($_) } 4..1000;
  my @acomp;  forcomposites { push @acomp, $_ } 1000;
  is_deeply(\@acomp, \@ecomp, "forcomposites to 1000");
}
{
  my @ecomp = grep { !is_prime($_) } 2000..2020;
  my @acomp;  forcomposites { push @acomp, $_ } 2000,2020;
  is_deeply(\@acomp, \@ecomp, "forcomposites 2000,2020");
}
{
  my $prod = 1;
  fordivisors { $prod *= $_ } 1234;
  is($prod, 1522756, "fordivisors 1234");
}

{
  my $nparts;
  is(partitions(25), 1958, "partitions(25)");
  $nparts = 0; forpart { $nparts++ } 25;
  is($nparts, 1958, "forpart {} 25 generates 1958 partitions");
  $nparts = 0; forpart { $nparts++ } 25,{n=>5};
  is($nparts, 192, "forpart {} 25,{n=>5} generates 192 partitions");
  $nparts = 0; forpart { $nparts++ } 25,{nmax=>5};
  is($nparts, 377, "forpart {} 25,{nmax=>5} generates 377 partitions");
}

{
  my $it = prime_iterator;
  my $sum = 0;
  $sum += $it->() for 1..100000;
  is($sum, 62260698721, "iterator sums first 100k primes");
}

{
  my $it = prime_iterator(200);
  is($it->(), 211, "prime_iterator(200)->()");
  is($it->(), 223, "prime_iterator(200)->()->()");
}

{
  my $sum = 0;
  my $it = prime_iterator_object;
  while ($it->value < 100) { $sum += $it->value; $it->next; }
  is($sum, 1060, "sum primes below 100 with OO iterator");
  is(vecsum(@{primes(100)}), 1060, "...with vecsum(primes(100))");
  $sum += $it->iterate for 1..100000;
  is($sum, 62293195902, "sum first 100k primes larger than 100");
  is(vecsum(@{primes(nth_prime(prime_count(100)+100000))}), 62293195902, "...with vecsum");
}

is(prime_count(1000), 168, "prime_count(1000)");
is(prime_count(1000,10000), 1061, "prime_count(1000,10000)");

cmp_closeto(prime_count_approx("1000000000000000000"),24739954287740860,1e-6,"prime_count_approx(1e18)");

is(twin_prime_count(123456), 1457, "twin_prime_count(123456)");
cmp_closeto(twin_prime_count_approx("100000000000000000"),90948839353159,1e-6,"twin_prime_count_approx(1e17)");

is(chinese([14,643], [254,419], [87,733]), 87041638, "chinese([14,643], [254,419], [87,733])");
is(vecsum(euler_phi(0,500_000)), 75991039676, "totient sum 500k");
is(invmod(42,2017),1969, "inverse of 42 mod 2017");

{
  my $sum = 0;
  $sum += exp_mangoldt($_) for 1..100;
  is($sum, 1156, "summatory von Mangoldt 1..100 = log(1156)");
}

{
  my $sum = 0;
  forprimes { $sum += log($_) } 12345;
  cmp_closeto(chebyshev_theta(12345), $sum, 1e-6, "chebyshev_theta(12345) and forprimes");
}
{
  my $sum = 0;
  for (1..12345) { $sum += log(exp_mangoldt($_)) }
  cmp_closeto(chebyshev_psi(12345), $sum, 1e-6, "chebyshev_psi(12345) and forprimes");
}

is(primorial(11), 2310, "primorial(11)");
is(pn_primorial(5), 2310, "pn_primorial(5)");
is(primorial(0), 1, "primorial(0)");
is(pn_primorial(0), 1, "pn_primorial(0)");

is(znorder(2, next_prime("10000000000000000")-6), 40177783100, "znorder(2,10000000000000061)");
is(legendre_phi(1000000000, 41), 106614188, "Legendre phi 1e9,41");

##############################################################################

# Not sure how best to test the random primes.
ok( is_prime(random_prime(1000)), "random_prime(1000)" );
ok( is_prime(random_prime(100,10000)), "random_prime(100,10000)" );
is( length(random_ndigit_prime(4)), 4, "random_ndigit_prime(4) is 4 digits" );
{
  my $bigprime;
  $bigprime = random_nbit_prime(512);
  is( length($bigprime->as_bin), 2+512, "random_nbit_prime(512) is 512 bits" );
  $bigprime = random_strong_prime(512);
  is( length($bigprime->as_bin), 2+512, "random_strong_prime(512) is 512 bits" );
  $bigprime = random_proven_prime(512);
  is( length($bigprime->as_bin), 2+512, "random_proven_prime(512) is 512 bits" );
}
# TODO: More of the random primes and certs

##############################################################################

is_deeply([factor("3369738766071892021")], [204518747,16476429743], "factor(3_369_738_766_071_892_021)");
is_deeply([factor_exp(29513484000)], [[2,5], [3,4], [5,3], [7,2], [11,1], [13,2]], "factor_exp(29513484000)");
is_deeply([factor(29513484000)], [2,2,2,2,2,3,3,3,3,5,5,5,7,7,11,13,13], "factor(29513484000)");
is_deeply([divisors(30)], [1, 2, 3, 5, 6, 10, 15, 30], "divisors(30)");

##############################################################################

{
  my $sum = 0;
  forcomposites { $sum += $_ if is_strong_pseudoprime($_,17) } 1000000;
  is($sum, 23206520, "forcomposites looking for base-17 strong probable primes");
}
{
  my($start,$end) = ("100000000000000000000", "100000000000000001000");
  my $aref = primes($start, $end);
  my $eref = [map { "100000000000000000".$_ } qw/039 129 151 193 207 301 349 361 391 393 441 477 547 559 561 721 741 753 757 763 801 853 961 993/];
  is_deeply($aref,$eref,"primes(10^20,10^20+1000)");
  $aref = [];
  forprimes { push @$aref,$_ } "100000000000000000039", "100000000000000000993";
}
{
  my @c;
  foroddcomposites { push @c,$_ if $_ % carmichael_lambda($_) == 1 } 10000;
  is_deeply(\@c,[qw/561 1105 1729 2465 2821 6601 8911/], "carmichael numbers under 10000");
  @c=();
  foroddcomposites { push @c,$_ if $_ % carmichael_lambda($_) == 1 } 1020000,1085000;
  is_deeply(\@c,[qw/1024651 1033669 1050985 1082809/], "carmichael numbers from 1020k to 1085k");
}
{
  my $nu3 = sub {
    my $n = shift;
    my($phix,$v) = (chebyshev_psi($n), 0);
    $v += (moebius($_)/$_)*LogarithmicIntegral($phix**(1/$_)) for 1..3;
    $v;
  };
  cmp_closeto($nu3->(1e6),  78498, 1e-4, "η3(1e6) ~ Pi(1e6)");
  cmp_closeto($nu3->(1e7), 664579, 1e-4, "η3(1e7) ~ Pi(1e7)");
}
{
  my $make_sg_it = sub {
    my $p = shift || 2;
    my $it = prime_iterator($p);
    return sub {
      do { $p = $it->() } while !is_prime(2*$p+1);
      $p;
    };
  };
  my $sgit = $make_sg_it->();
  my $sum = 0;
  $sum += $sgit->() for 1..10000;
  is($sum, 6171027819, "sum first 10k Sophie-Germain primes using iterator");
}
is( (factor("600851475143"))[-1], 6857, "largest prime factor of 600851475143");
is( nth_prime(10001), 104743, "nth_prime(10001)");
{
  my $sum = 0;
  forprimes { $sum += $_ } 2_000_000;
  is($sum, 142913828922, "sum 2M primes with forprimes");
  is(vecsum( @{primes(2_000_000)} ), 142913828922, "sum 2M primes with vecsum(primes())");
}
{
  my $sum = 0;
  foreach my $x (1..10000) {
    my $y = divisor_sum($x)-$x;
    $sum += $x + $y if $y > $x && $x == divisor_sum($y)-$y;
  }
  is($sum, 31626, "sum of amicable numbers using loop");
  $sum = vecsum( map { divisor_sum($_) }
                 grep { my $y = divisor_sum($_)-$_;
                        $y > $_ && $_==(divisor_sum($y)-$y) }
                 1 .. 10000 );
  is($sum, 31626, "sum of amicable numbers using pipeline");
}
{
  my $pd = first { /1/&&/2/&&/3/&&/4/&&/5/&&/6/&&/7/} reverse @{primes(1000000,9999999)};
  is($pd, 7652413, "largest 7-digit pandigital prime");
}
{
  my $n = pn_primorial(4);
  $n++ while (factor_exp($n) != 4 || factor_exp($n+1) != 4 || factor_exp($n+2) != 4 || factor_exp($n+3) != 4);
  is($n, 134043, "first number in sequence of four 4-factor numbers");
}
{
  my ($maxn, $maxratio) = (0,0);
  foreach my $n (1..1000000) {
    my $ndivphi = $n / euler_phi($n);
    ($maxn, $maxratio) = ($n, $ndivphi) if $ndivphi > $maxratio;
  }
  is($maxn, 510510, "largest ratio of n/phi(n) for n < 1M.  Brute force.");
}
{
  my $n = 0;
  $n++ while pn_primorial($n+1) < 1000000;
  my $maxn = pn_primorial($n);
  is($maxn, 510510, "largest ratio of n/phi(n) for n < 1M.  Smart way.");
}
{
  my $n=0;
  forcomposites { $n++ if scalar factor($_) == 2; } int(1e5)-1;
  is($n, 23378, "Brute force count semiprimes under 1e5");
}
{
  my $limit = 1e7;
  $limit--;
  my ($sum, $pc) = (0, 1);
  forprimes {
    $sum += prime_count(int($limit/$_)) + 1 - $pc++;
  } int(sqrt($limit));
  is($sum, 1904324, "Count of semiprimes under 1e7");
}
{
  my $matches = sub {
    my @d = divisors(shift);
    return map { [$d[$_],$d[$#d-$_]] } 1..(@d-1)>>1;
  };
  is_deeply([$matches->(139650)], [[2,69825],[3,46550],[5,27930],[6,23275],[7,19950],[10,13965],[14,9975],[15,9310],[19,7350],[21,6650],[25,5586],[30,4655],[35,3990],[38,3675],[42,3325],[49,2850],[50,2793],[57,2450],[70,1995],[75,1862],[95,1470],[98,1425],[105,1330],[114,1225],[133,1050],[147,950],[150,931],[175,798],[190,735],[210,665],[245,570],[266,525],[285,490],[294,475],[350,399]], "matches 139650");
}
{
  my @nums;
  forcomposites { push @nums,$_ if divisor_sum($_)+6==divisor_sum($_+6) } 9,1e5;
  is_deeply(\@nums,[qw/104 147 596 1415 4850 5337/], "OEIS A054903");
}
{
  my @s;
  foreach my $n (1..30) {
    if (!znprimroot($n)) {
      push @s, "$n -";
    } else {
      my $phi = euler_phi($n);
      my @r = grep { gcd($_,$n) == 1 && znorder($_,$n) == $phi } 1..$n-1;
      push @s, "$n " . join(" ", @r);
    }
  }
  my @expect = split(/\|/, "1 -|2 1|3 2|4 3|5 2 3|6 5|7 3 5|8 -|9 2 5|10 3 7|11 2 6 7 8|12 -|13 2 6 7 11|14 3 5|15 -|16 -|17 3 5 6 7 10 11 12 14|18 5 11|19 2 3 10 13 14 15|20 -|21 -|22 7 13 17 19|23 5 7 10 11 14 15 17 19 20 21|24 -|25 2 3 8 12 13 17 22 23|26 7 11 15 19|27 2 5 11 14 20 23|28 -|29 2 3 8 10 11 14 15 18 19 21 26 27|30 -|");
  is_deeply(\@s,\@expect,"znprimroot table 1..30");
}

##############################################################################

{
  my $checksum = vecreduce { $a ^ $b } @{twin_primes(1000000)};
  is($checksum, 630871, "xor of twin primes <= 1M");
}

##############################################################################

{
  my @v = (qw/a b c d e/);
  my $ps = join " ", map { join("",vecextract(\@v,$_)) } 0..2**scalar(@v)-1;
  is($ps, " a b ab c ac bc abc d ad bd abd cd acd bcd abcd e ae be abe ce ace bce abce de ade bde abde cde acde bcde abcde", "power set of 5 elements");

  my $word = join "", vecextract(["a".."z"], [15, 17, 8, 12, 4]);
  is($word, "prime", "use vecextract with array");
}

##############################################################################

tie my @primes, 'Math::Prime::Util::PrimeArray';
{
  my @plist;
  for my $n (0..9) { push @plist, $primes[$n]; }
  is_deeply(\@plist, primes(nth_prime(10)), "PrimeArray for index loop");
}
{
  my @plist;
  for my $p (@primes) { last if $p > 79; push @plist, $p; }
  is_deeply(\@plist, primes(79), "PrimeArray for primes loop");
}
{
  my @plist;
  is_deeply([@primes[0..49]], primes(nth_prime(50)), "PrimeArray array slice");
}
SKIP: {
  skip "hash each requires 5.12 or newer", 1 if $] < 5.012;
  my @plist;
  while (  my($index,$value) = each @primes ) {
    last if $value > 147;
    push @plist, $value;
  }
  is_deeply(\@plist, primes(147), "PrimeArray each primes loop");
}
{
  my @plist;
  while ((my $p = shift @primes) < 250) { push @plist, $p; }
  is_deeply(\@plist, primes(250), "PrimeArray shift");
  unshift @primes, ~0;  # put primes back.
  is($primes[0], 2, "unshift puts it back");
}

##############################################################################

sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  my $error = sprintf("%g", abs($got - $expect) / $expect);
  my $errorpr = sprintf "%.2g", $error;
  cmp_ok( $error, '<=', $tolerance, "$message ($errorpr)");
}
