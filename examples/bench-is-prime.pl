#!/usr/bin/env perl
use strict;
use warnings;
#use Math::Primality;
use Math::Prime::XS;
use Math::Prime::Util;
#use Math::Pari;
#use Math::Prime::FastSieve;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -5;
my $numbers = 1000;

my $is64bit = (~0 > 4294967295);
my $maxdigits = ($is64bit) ? 20 : 10;  # Noting the range is limited for max.

srand(29);
test_at_digits($_) for (5 .. $maxdigits);


sub test_at_digits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. $numbers);
  my $min_num = min @nums;
  my $max_num = max @nums;

  #my $sieve = Math::Prime::FastSieve::Sieve->new(10 ** $magnitude + 1);
  #Math::Prime::Util::prime_precalc(10 ** $magnitude + 1);

  print "is_prime for $numbers random $digits-digit numbers ($min_num - $max_num)\n";

  cmpthese($count,{
    #'Math::Primality' => sub { Math::Primality::is_prime($_) for @nums },
    'M::P::XS' => sub { Math::Prime::XS::is_prime($_) for @nums },
    #'M::P::FS' => sub { $sieve->isprime($_) for @nums },
    'M::P::U' => sub { Math::Prime::Util::is_prime($_) for @nums },
    'MPU prob' => sub { Math::Prime::Util::is_prob_prime($_) for @nums },
    #'Math::Pari' => sub { Math::Pari::isprime($_) for @nums },
  });
  print "\n";
}
