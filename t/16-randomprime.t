#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
#use Math::Random::MT qw/rand/;
#use Math::Random::MT::Auto qw/rand/;
#sub rand { return 0.5; }
use Math::Prime::Util qw/random_prime random_ndigit_prime random_nbit_prime
                         random_maurer_prime random_shawe_taylor_prime
                         random_proven_prime
                         random_semiprime random_unrestricted_semiprime
                         random_safe_prime random_strong_prime
                         factor is_prime is_semiprime is_smooth logint
                         prime_set_config/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $broken64 = (18446744073709550592 == ~0);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;
my $do_st = 1;
$do_st = 0 unless eval { require Digest::SHA;
                         my $version = $Digest::SHA::VERSION;
                         $version =~ s/[^\d.]//g;
                         $version >= 4.00; };

my @random_to = (2, 3, 4, 5, 6, 7, 8, 9, 100, 1000, 1000000, 4294967295);

my @random_nbit_tests = ( 2 .. 10, 15 .. 17, 28, 32 );
push @random_nbit_tests, (34) if $use64;
@random_nbit_tests = (2 .. $maxbits) if $extra;
push @random_nbit_tests, (75);

my @random_ndigit_tests = (1 .. ($use64 ? 20 : 10));

if ($use64 && $broken64) {
  diag "Skipping some tests with broken 64-bit Perl\n";
  @random_ndigit_tests = grep { $_ < 10 } @random_ndigit_tests;
  @random_nbit_tests = grep { $_ < 50 } @random_nbit_tests;
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
my %range_edge_empty = (
  "0 to 0" => [],
  "0 to 1" => [],
  "2 to 1" => [],
  "3 to 2" => [],
  "1294268492 to 1294268778" => [],
  "3842610774 to 3842611108" => [],
);

plan tests => 15+1+1+1+1
              + (1 * scalar (keys %range_edge_empty))
              + (3 * scalar (keys %range_edge))
              + (2 * scalar (keys %ranges))
              + (2 * scalar @random_to)
              + (1 * scalar @random_ndigit_tests)
              + (2 * scalar @random_nbit_tests)
              + 4
              + 7+14 # random_semiprime
              + 11   # random_safe_prime
              + 1+6  # random_strong_prime
              + 0;

my $infinity = 20**20**20;
my $superbig = 1000000003;
my $nrandom_range_samples = $extra ? 1000 : 50;

ok(!eval { random_prime(1); }, "random_prime(1)"); # Just for coverage
ok(!eval { random_prime(2,1); }, "random_prime(2,1)"); # Just for coverage
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

#ok(!eval { random_ndigit_prime(undef); }, "random_ndigit_prime(undef)");
ok(!eval { random_ndigit_prime(0); }, "random_ndigit_prime(0)");
#ok(!eval { random_ndigit_prime(-5); }, "random_ndigit_prime(-5)");

#ok(!eval { random_nbit_prime(undef); }, "random_nbit_prime(undef)");
ok(!eval { random_nbit_prime(0); }, "random_nbit_prime(0)");
#ok(!eval { random_nbit_prime(-5); }, "random_nbit_prime(-5)");

#ok(!eval { random_maurer_prime(undef); }, "random_maurer_prime(undef)");
ok(!eval { random_maurer_prime(0); }, "random_maurer_prime(0)");
#ok(!eval { random_maurer_prime(-5); }, "random_maurer_prime(-5)");

#ok(!eval { random_shawe_taylor_prime(undef); }, "random_shawe_taylor_prime(undef)");
ok(!eval { random_shawe_taylor_prime(0); }, "random_shawe_taylor_prime(0)");
#ok(!eval { random_shawe_taylor_prime(-5); }, "random_shawe_taylor_prime(-5)");

while (my($range, $expect) = each (%range_edge_empty)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is( random_prime($low,$high), undef, "primes($low,$high) should return undef" );
}

while (my($range, $expect) = each (%range_edge)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  my $got = random_prime($low,$high);
  ok( is_prime($got), "Prime in range $low-$high is indeed prime" );
  cmp_ok( $got, '>=', $expect->[0], "random_prime($low,$high) >= $expect->[0]");
  cmp_ok( $got, '<=', $expect->[1], "random_prime($low,$high) <= $expect->[1]");
}

while (my($range, $expect) = each (%ranges)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  my $isprime = 1;
  my $inrange = 1;
  for (1 .. $nrandom_range_samples) {
    my $got = random_prime($low,$high);
    $isprime *= is_prime($got) ? 1 : 0;
    $inrange *= (($got >= $expect->[0]) && ($got <= $expect->[1])) ? 1 : 0;
  }
  ok($isprime, "All returned values for $low-$high were prime" );
  ok($inrange, "All returned values for $low-$high were in the range" );
}

# We want to test the no-bigint stuff here.  This makes calls for 10-digit
# (32-bit) and 20-digit (64-bit) random primes stay inside native range.
prime_set_config(nobigint=>1);

foreach my $high (@random_to) {
  my $isprime = 1;
  my $inrange = 1;
  for (1 .. $nrandom_range_samples) {
    my $got = random_prime($high);
    $isprime *= is_prime($got) ? 1 : 0;
    $inrange *= (($got >= 2) && ($got <= $high)) ? 1 : 0;
  }
  ok($isprime, "All returned values for $high were prime" );
  ok($inrange, "All returned values for $high were in the range" );
}

foreach my $digits ( @random_ndigit_tests ) {
  my $n = random_ndigit_prime($digits);
  ok ( length($n) == $digits && is_prime($n),
       "$digits-digit random prime '$n' is in range and prime");
}

foreach my $bits ( @random_nbit_tests ) {
  check_bits( random_nbit_prime($bits), $bits, "nbit" );
  # These are handled in t/23-random-certs.t
  #check_bits( random_maurer_prime($bits), $bits, "Maurer" );
  #SKIP: {
  #  skip "random Shawe-Taylor prime generation requires Digest::SHA",1 unless $do_st;
  #  check_bits( random_shawe_taylor_prime($bits), $bits, "Shawe-Taylor" );
  #}
  check_bits( random_proven_prime($bits), $bits, "proven" );
}

sub check_bits {
  my($n, $bits, $what) = @_;
  ok ( is_nbit($n, $bits) && is_prime($n),
       "$bits-bit random $what prime '$n' is in range and prime");
}
prime_set_config(nobigint=>0);

{
  my $n = random_nbit_prime(80);
  ok( ref($n) =~ /^Math::/, "random 80-bit prime returns a BigInt" );
  ok( 1+logint($n,2) == 80, "random 80-bit prime '$n' is in range" );
}
SKIP: {
  skip "Skipping 30-digit random prime with broken 64-bit Perl", 2 if $broken64;
  my $n = random_ndigit_prime(30);
  ok( ref($n) =~ /^Math::/, "random 30-digit prime returns a BigInt" );
  ok( 1+logint($n,10) == 30, "random 30-digit prime '$n' is in range" );
}

###### semiprimes
{
  my $n;
  ok(!eval { random_semiprime(3); }, "random_semiprime(3)");
  ok(!eval { random_unrestricted_semiprime(2); }, "random_unrestricted_semiprime(2)");
  is(random_semiprime(4),9,"random_semiprime(4) = 9");
  $n = random_unrestricted_semiprime(3);
  ok($n ==4 || $n == 6, "random_unrestricted_semiprime(3) is 4 or 6");
  for my $bits (4 .. 10) {
    check_semi_bits(random_semiprime($bits), $bits, "random_semiprime");
    check_semi_bits(random_unrestricted_semiprime($bits), $bits, "random_unrestricted_semiprime");
  }
  $n = random_semiprime(26);
  ok($n >= 33554432 && $n < 67108864 && scalar(factor($n)) == 2, "random_semiprime(26) is a 26-bit semiprime");
  $n = random_semiprime(81);
  ok( 1+logint($n,2) == 81, "random_semiprime(81) is 81 bits");
  SKIP: {
    skip "Skipping 81-bit semiprime with broken 64-bit Perl", 1 if $broken64;
    $n = random_unrestricted_semiprime(81);
    ok( 1+logint($n,2) == 81, "random_unrestricted_semiprime(81) is 81 bits");
  }
}

sub check_semi_bits {
  my($n, $bits, $name) = @_;
  ok ( is_nbit($n, $bits) && is_semiprime($n),
       "$name($bits) is in range and semiprime");
}

###### Safe primes
ok(!eval { random_safe_prime(2); }, "random_safe_prime(2) is invalid");
# This can be very slow over 65 bits
for my $bits (3 .. 10, 48, 70) {
  my $p = random_safe_prime($bits);
  my $q = ($p-1) >> 1;
  ok ( is_nbit($p, $bits) && is_prime($p) && is_prime($q),
       "random_safe_prime($bits) is in range and is a safe prime");
}

###### Strong primes
ok(!eval { random_strong_prime(127); }, "random_strong_prime(127) is invalid");
for my $bits (128, 247, 512) {
  SKIP: {
    skip "skipping random_strong_prime($bits) without GMP", 2 unless $usegmp;
    my $p = random_strong_prime($bits);
    ok(is_nbit($p, $bits) && is_prime($p), "random_strong_prime($bits) is a prime with $bits bits");
    ok (!is_smooth($p-1, 10000) && !is_smooth($p+1, 10000),
        "random_strong_prime($bits) isn't obviously weak");
  }
}

sub is_nbit {
  my($n, $bits) = @_;
  return 1+logint($n,2) == $bits ? 1 : 0;
}
