#!/usr/bin/env perl
use strict;
use warnings;
#use Math::Primality;
use Math::Prime::XS;
use Math::Prime::Util qw/urandomm/;;
#use Math::Pari;
#use Math::Prime::FastSieve;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -5;
my $numbers = 1000;

my $is64bit = (~0 > 4294967295);
my $maxdigits = ($is64bit) ? 20 : 10;  # Noting the range is limited for max.
use Math::Prime::Util::RandomPrimes;

my $rand_ndigit_gen = sub {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;
  my $howmany = shift || 1;
  my ($base, $max);

  if ( 10**$digits < ~0) {
    $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
    $max = int(10 ** $digits);
    $max = ~0 if $max > ~0;
  } else {
    $base = Math::BigInt->new(10)->bpow($digits-1);
    $max = Math::BigInt->new(10)->bpow($digits) - 1;
  }
  #my @nums = map { $base + urandomm($max-$base) } (1 .. $howmany);
  my @nums;
  while (@nums < $howmany) {
    my $n = $base + urandomm($max-$base);
    push @nums, $n if $n % 2 && $n % 3 && $n % 5 && $n % 7;
  }
  return (wantarray) ? @nums : $nums[0];
};

srand(29);
test_at_digits($_) for (3 .. $maxdigits);


sub test_at_digits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;

  my @nums = $rand_ndigit_gen->($digits, $numbers);
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
