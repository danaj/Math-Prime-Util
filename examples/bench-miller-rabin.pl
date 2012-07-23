#!/usr/bin/env perl
use strict;
use warnings;
use Math::Primality;
use Math::Prime::XS;
use Math::Prime::Util '-nobigint';
use Math::Prime::Util::GMP;
#use Math::Prime::FastSieve;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -5;

srand(29);
test_at_digits($_) for (5..10);


sub test_at_digits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. 1000);
  my $min_num = min @nums;
  my $max_num = max @nums;

  #my $sieve = Math::Prime::FastSieve::Sieve->new(10 ** $magnitude + 1);
  #Math::Prime::Util::prime_precalc(10 ** $magnitude + 1);

  print "miller_rabin for 1000 random $digits-digit numbers ($min_num - $max_num)\n";

  cmpthese($count,{
    'MPU' => sub { Math::Prime::Util::is_strong_pseudoprime($_,2,3,5,7,11,13,17) for @nums },
    'MPU GMP' => sub { Math::Prime::Util::GMP::is_strong_pseudoprime($_,2,3,5,7,11,13,17) for @nums },
    'M:Primality' => sub { for (@nums) {
                       Math::Primality::is_strong_pseudoprime($_,2) &&
                       Math::Primality::is_strong_pseudoprime($_,3) &&
                       Math::Primality::is_strong_pseudoprime($_,5) &&
                       Math::Primality::is_strong_pseudoprime($_,7) &&
                       Math::Primality::is_strong_pseudoprime($_,11) &&
                       Math::Primality::is_strong_pseudoprime($_,13) &&
                       Math::Primality::is_strong_pseudoprime($_,17); } },
  });
  print "\n";
}
