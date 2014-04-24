#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/is_aks_prime is_prime primes/;
$| = 1;  # fast pipes

my $limit = shift || 2_000_000_000;
my $nrand = 8000;

my %isprime = map { $_ => 1 } @{primes(160_000)};

print "Testing AKS for all numbers from 1 to 160,000:\n";
foreach my $n (1 .. 160_000) {
  print "." unless $n % 2000;
  if ($isprime{$n}) {
    die "\n$n is prime\n" unless is_aks_prime($n);
  } else {
    die "\n$n is composite\n" if is_aks_prime($n);
  }
}
print "\n";
print "Testing $nrand random numbers from 1 to $limit:\n";
for (1 .. $nrand) {
  print "." unless $_ % 100;
  my $n = 1 + int(rand($limit));
  if (is_prime($n)) {
    die "\n$n is prime\n" unless is_aks_prime($n);
  } else {
    die "\n$n is composite\n" if is_aks_prime($n);
  }
}
