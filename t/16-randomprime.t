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

#my $do_st = 1;
#$do_st = 0 unless eval { require Digest::SHA;
#                         my $version = $Digest::SHA::VERSION;
#                         $version =~ s/[^\d.]//g;
#                         $version >= 4.00; };

plan tests =>   1  # expected failures
              + 1  # random_prime(lo,hi)
              + 1  # random_prime(hi)
              + 1  # random_ndigit_prime
              + 1  # random_nbit_prime
              + 1  # 80-bit + 30-digit random
              + 1  # random_semiprime
              + 1  # random_safe_prime
              + 1  # random_strong_prime
              + 0;

subtest 'expected failures', sub {
  my $infinity = 20**20**20;
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
};

subtest 'random_prime(lo,hi)', sub {
  is_deeply([map { random_prime($_->[0],$_->[1]) } ([0,0],[0,1],[2,1],[3,2],[1294268492,1294268778],[3842610774,3842611108])],
            [undef,undef,undef,undef,undef,undef],
            "random_prime(lo,hi) returns undef when no primes in range");

  my @edges = ( [0,2,2,2], [2,2,2,2], [2,3,2,3], [3,5,3,5],
                [10,20,11,19], [8,12,11,11], [10,12,11,11],
                [16706143,16706143,16706143,16706143],
                [16706142,16706144,16706143,16706143],
                [3842610773,3842611109,3842610773,3842611109],
                [3842610772,3842611110,3842610773,3842611109] );
  for my $edata (@edges) {
    my($lo,$hi,$lores,$hires) = @$edata;
    my $got = random_prime($lo,$hi);
    ok($got >= $lores && $got <= $hires && is_prime($got), "($lo,$hi) => $got in [$lores,$hires]");
  }

  my @ranges = ( [2,20,2,19], [3,7,3,7], [20,100,23,97], [5678,9876,5683,9871],
                 [27767,88493,27767,88493], [27764,88498,27767,88493],
                 [27764,88493,27767,88493], [27767,88498,27767,88493],
                 [17051687,17051899,17051687,17051899],
                 [17051688,17051898,17051707,17051887] );
  for my $rdata (@ranges) {
    my($lo,$hi,$lores,$hires) = @$rdata;
    my $got = random_prime($lo,$hi);
    ok($got >= $lores && $got <= $hires && is_prime($got), "($lo,$hi) => $got in [$lores,$hires]");
  }
};


subtest 'random_prime(hi)', sub {
  my @random_to = (2, 3, 4, 5, 6, 7, 8, 100, 1000, 1000000, 4294967295);
  my $nrandom_range_samples = $extra ? 1000 : 50;

  foreach my $high (@random_to) {
    my $isprime = 1;
    my $inrange = 1;
    for (1 .. $nrandom_range_samples) {
      my $got = random_prime($high);
      $isprime = 0 if !is_prime($got);
      $inrange = 0 if $got < 2 || $got > $high;
    }
    ok($isprime && $inrange, "returned prime values in [2,$high]");
  }
};


subtest 'random_ndigit_prime', sub {

  prime_set_config(nobigint=>1);  # No bigints returned for maxdigits

  my @T = $use64 ? (1..11,15,19,20) : (1..10);
  foreach my $digits (@T) {
    SKIP: {
      skip "Broken 64-bit, skipping random_ndigit_prime($digits)",1 if $use64 && $broken64 && $digits >= 10;
      my $n = random_ndigit_prime($digits);
      ok(!ref($n) && length($n) == $digits && is_prime($n), "($digits) is a $digits-digit prime (got $n)");
    }
  }
  prime_set_config(nobigint=>0);  # Turn this back off
};


subtest 'random_nbit_prime', sub {
  my @T = $extra ? (2..$maxbits)
        : $use64 ? (2..10,15..17,28,32,34)
        :          (2..10,15..17,28,32);
  foreach my $bits (@T) {
    SKIP: {
      skip "Broken 64-bit, skipping random_nbit_prime($bits)",1 if $use64 && $broken64 && $bits >= 50;
      check_bits( random_nbit_prime($bits), $bits, "nbit" );
    }
  }
};

# These are handled in t/23-random-certs.t
# check_bits( random_maurer_prime($bits), $bits, "Maurer" );
# check_bits( random_shawe_taylor_prime($bits), $bits, "Shawe-Taylor" );
# check_bits( random_proven_prime($bits), $bits, "proven" );



subtest 'large random nbit/ndigit', sub {
  my $n = random_nbit_prime(80);
  ok( ref($n) =~ /^Math::/, "random 80-bit prime returns a BigInt" );
  ok( 1+logint($n,2) == 80, "random 80-bit prime '$n' is in range" );
  SKIP: {
    skip "Skipping 30-digit random prime with broken 64-bit Perl", 2 if $broken64;
    my $n = random_ndigit_prime(30);
    ok( ref($n) =~ /^Math::/, "random 30-digit prime returns a BigInt" );
    ok( 1+logint($n,10) == 30, "random 30-digit prime '$n' is in range" );
  }
};


subtest 'semiprimes', sub {
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
};


subtest 'safe primes', sub {
  ok(!eval { random_safe_prime(2); }, "random_safe_prime(2) is invalid");
  # This can be very slow over 65 bits
  for my $bits (3, 5, 8, 40, 70) {
    SKIP: {
    skip "Skip larger safe prime on 32-bit",1 if $bits > 50 && !$use64;
      my $p = random_safe_prime($bits);
      my $q = ($p-1) >> 1;
      ok ( is_nbit($p, $bits) && is_prime($p) && is_prime($q),
           "random_safe_prime($bits) is in range and is a safe prime");
    }
  }
};


subtest 'strong primes', sub {
  ok(!eval { random_strong_prime(127); }, "random_strong_prime(127) throws error as expected");
  for my $bits (128, 247, 512) {
    SKIP: {
      skip "skipping random_strong_prime($bits) without GMP", 2 unless $usegmp;
      my $p = random_strong_prime($bits);
      check_bits($p,$bits, 'strong');
      ok (!is_smooth($p-1, 10000) && !is_smooth($p+1, 10000),
          "random_strong_prime($bits) isn't obviously weak");
    }
  }
};


sub is_nbit {
  my($n, $bits) = @_;
  return 1+logint($n,2) == $bits  ?  1  :  0;
}
sub check_bits {
  my($n, $bits, $what) = @_;
  my $strn = length($n) <= 24 ? "$n" : substr($n,0,12)."...".substr($n,-12,12);
  ok(is_nbit($n,$bits) && is_prime($n), "random $what prime '$strn' is a $bits-bit prime");
}
sub check_semi_bits {
  my($n, $bits, $name) = @_;
  ok(is_nbit($n,$bits) && is_semiprime($n), "$name($bits) is in range and semiprime");
}
