#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_safe_prime is_prime random_safe_prime/;

my $use64 = Math::Prime::Util::_uvbits >= 64;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

# A005385: safe primes p where (p-1)/2 is also prime
my @safe_primes = (5,7,11,23,47,59,83,107,167,179,227,263,347,359,383,467,479,
                   503,563,587,719,839,863,887,983,1019,1187,1283,1307,
                   1319,1367,1439,1487,1523,1619,1823,1907);

# Small non-safe-primes (composites and primes that are not safe)
my @not_safe = (0,1,2,3,4,6,8,9,10,12,13,17,19,25,29,31,37,41,43,53,61,67,
                71,73,79,89,97,101,103,109,113,127,131,137,139,149,151,157,
                163,173,181,191,193,197,199,1031);

plan tests => 4;

{
  my @got = map { is_safe_prime($_) ? $_ : 0 } @safe_primes;
  is_deeply(\@got, \@safe_primes, "is_safe_prime(n) for safe primes");
}
{
  my @got = map { !is_safe_prime($_) ? $_ : 0 } @not_safe;
  is_deeply(\@got, \@not_safe, "is_safe_prime(n) for non-safe primes");
}
subtest 'is_safe_prime consistency' => sub {
  # If p is a safe prime then p is prime and (p-1)/2 is prime
  my @got1 = map { is_prime($_) ? $_ : 0 } @safe_primes;
  is_deeply(\@got1, \@safe_primes, "safe prime p is prime");
  my @got2 = map { is_prime(($_-1)/2) ? $_ : 0 } @safe_primes;
  is_deeply(\@got2, \@safe_primes, "  and (p-1)/2 is prime");
};

subtest 'big prime' => sub {
  is(is_safe_prime(890471062883),1,"890471062883 is a safe prime");
  is(is_safe_prime(800579858627),0,"800579858627 is not a safe prime");
  is(is_safe_prime(764953310311),0,"764953310311 is not a safe prime");

  if ($use64 || $extra) {
    is(is_safe_prime("284651113549152086747"),1,"284651113549152086747 is a safe prime");
    is(is_safe_prime("220961958687770434157"),0,"220961958687770434157 is not a safe prime");
    is(is_safe_prime("270708094586485285875"),0,"270708094586485285875 is not a safe prime");
  }
};
