#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util;
use Math::Prime::Util::GMP;
use Math::Primality;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -2;
srand(29);  # So we have repeatable results

test_at_digits($_, 1000) for (5, 15, 25, 50, 200);

sub test_at_digits {
  my($digits, $numbers) = @_;
  die "Digits must be > 0" unless $digits > 0;

  # We get a mix of primes and non-primes.
  my @nums = map { Math::Prime::Util::random_ndigit_prime($digits)+2 } 1 .. $numbers;
  print "is_strong_pseudoprime for $numbers random $digits-digit numbers",
        " (", min(@nums), " - ", max(@nums), ")\n";

  cmpthese($count,{
    'MP'      =>sub {Math::Primality::is_strong_pseudoprime($_,3) for @nums;},
    'MPU'     =>sub {Math::Prime::Util::is_strong_pseudoprime($_,3) for @nums;},
    'MPU PP'  =>sub {Math::Prime::Util::PP::miller_rabin($_,3) for @nums;},
    'MPU GMP' =>sub {Math::Prime::Util::GMP::is_strong_pseudoprime($_,3) for @nums;},
  });
}
