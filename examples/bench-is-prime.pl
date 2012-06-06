#!/usr/bin/perl -w

use strict;
#use Math::Primality;
use Math::Prime::XS;
use Math::Prime::Util;
#use Math::Prime::FastSieve;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -5;

test_at_mag($_) for (4, 8, 12, 16);


sub test_at_mag {
  my $magnitude = shift;

  my $base = 10 ** ($magnitude-1);
  my $max = 10 ** $magnitude;
  my $digits = $magnitude+1;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. 200);
  my $min_num = min @nums;
  my $max_num = max @nums;

  #my $sieve = Math::Prime::FastSieve::Sieve->new(10 ** $magnitude + 1);
  #Math::Prime::Util::prime_precalc(10 ** $magnitude + 1);

  print "is_prime for random $digits-digit numbers ($min_num - $max_num)\n";

  cmpthese($count,{
    #'Math::Primality' => sub { Math::Primality::is_prime($_) for @nums },
    'M::P::XS' => sub { Math::Prime::XS::is_prime($_) for @nums },
    #'M::P::FS' => sub { $sieve->isprime($_) for @nums },
    'M::P::U' => sub { Math::Prime::Util::is_prime($_) for @nums },
    'MPU prob' => sub { Math::Prime::Util::is_prob_prime($_) for @nums },
  });
  print "\n";
}
