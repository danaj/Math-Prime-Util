#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util ":all";
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $maxdigits = (~0 <= 4294967295) ? 10 : 20;

my $count = shift || -5;

srand(29);
test_at_digits($_) for (5 .. $maxdigits);


sub test_at_digits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. 1000);
  my $min_num = min @nums;
  my $max_num = max @nums;

  #print "miller_rabin for 1000 random $digits-digit numbers ($min_num - $max_num)\n";

  my $sum;
  cmpthese($count,{
    'lower' => sub { $sum += prime_count_lower($_) for @nums },
    'luapprox' => sub { $sum += (prime_count_lower($_)+prime_count_upper($_))/2 for @nums },
    'approx' => sub { $sum += prime_count_approx($_) for @nums },
    'li' => sub { $sum += LogarithmicIntegral($_) for @nums },
    'R' => sub { $sum += RiemannR($_) for @nums },
  });
  print "\n";
}
