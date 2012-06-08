#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/random_prime random_ndigit_prime is_prime/;

my $use64 = Math::Prime::Util::_maxbits > 32;
my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};

plan tests => 13 + 6*1 + 11*3 + 10*2 + 12*2;

my $infinity = ~0 * ~0;

ok(!eval { random_prime(undef); }, "random_prime(undef)");
ok(!eval { random_prime(-3); }, "random_prime(-3)");
ok(!eval { random_prime("a"); }, "random_prime(a)");
ok(!eval { random_prime(undef,undef); }, "random_prime(undef,undef)");
ok(!eval { random_prime(2,undef); }, "random_prime(2,undef)");
ok(!eval { random_prime(2,"a"); }, "random_prime(2,a)");
ok(!eval { random_prime(undef,0); }, "random_prime(undef,0)");
ok(!eval { random_prime(0,undef); }, "random_prime(0,undef)");
ok(!eval { random_prime(2,undef); }, "random_prime(2,undef)");
ok(!eval { random_prime(2,-4); }, "random_prime(2,-4)");
ok(!eval { random_prime(2,$infinity); }, "random_prime(2,+infinity)");
ok(!eval { random_prime($infinity); }, "random_prime(+infinity)");
ok(!eval { random_prime(-$infinity); }, "random_prime(-infinity)");

my %range_edge_empty = (
  "0 to 0" => [],
  "0 to 1" => [],
  "2 to 1" => [],
  "3 to 2" => [],
  "1294268492 to 1294268778" => [],
  "3842610774 to 3842611108" => [],
);
while (my($range, $expect) = each (%range_edge_empty)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is( random_prime($low,$high), undef, "primes($low,$high) should return undef" );
}

my %range_edge = (
  "0 to 2" => [2,2],
  "2 to 2" => [2,2],
  "2 to 3" => [2,3],
  "3 to 5" => [3,5],
  "10 to 20" => [11,19],
  "8 to 12" => [11,11],
  "10 to 12" => [11,11],
  "16706143 to 16706143" => [16706143,16706143],
  "16706142 to 16706144" => [16706143,16706143],
  "3842610773 to 3842611109" => [3842610773,3842611109],
  "3842610772 to 3842611110" => [3842610773,3842611109],
);
while (my($range, $expect) = each (%range_edge)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  my $got = random_prime($low,$high);
  ok( is_prime($got), "Prime in range $low-$high is indeed prime" );
  cmp_ok( $got, '>=', $expect->[0], "random_prime($low,$high) >= $expect->[0]");
  cmp_ok( $got, '<=', $expect->[1], "random_prime($low,$high) >= $expect->[1]");
}

my %ranges = (
  "2 to 20" => [2,19],
  "3 to 7" => [3,7],
  "20 to 100" => [23,97],
  "5678 to 9876" => [5683,9871],
  "27767 to 88493" => [27767,88493],
  "27764 to 88498" => [27767,88493],
  "27764 to 88493" => [27767,88493],
  "27767 to 88498" => [27767,88493],
  "17051687 to 17051899" => [17051687,17051899],
  "17051688 to 17051898" => [17051707,17051887],
);
while (my($range, $expect) = each (%ranges)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  my $isprime = 1;
  my $inrange = 1;
  for (1 .. 1000) {
    my $got = random_prime($low,$high);
    $isprime *= is_prime($got) ? 1 : 0;
    $inrange *= (($got >= $expect->[0]) && ($got <= $expect->[1])) ? 1 : 0;
  }
  ok($isprime, "All returned values for $low-$high were prime" );
  ok($inrange, "All returned values for $low-$high were in the range" );
}

my @to = (2, 3, 4, 5, 6, 7, 8, 9, 100, 1000, 1000000, 4294967295);
foreach my $high (@to) {
  my $isprime = 1;
  my $inrange = 1;
  for (1 .. 1000) {
    my $got = random_prime($high);
    $isprime *= is_prime($got) ? 1 : 0;
    $inrange *= (($got >= 2) && ($got <= $high)) ? 1 : 0;
  }
  ok($isprime, "All returned values for $high were prime" );
  ok($inrange, "All returned values for $high were in the range" );
}
