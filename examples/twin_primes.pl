#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/prime_iterator  prime_iterator_object
                         next_prime  is_prime
                         nth_prime_upper  nth_twin_prime_approx  prime_precalc/;

my $count = shift || 20;

# Find twin primes (numbers where p and p+2 are prime)

# Time for the first 300k:
#
# Not iterators:
#   0.5s   $l=2; forprimes { say $l if $l+2==$_; $l=$_; } 64764841
#   0.8s   bin/primes.pl --twin 2 64764839
#
# Iterators with precalc:
#   1.5s   get_twin_prime_iterator2
#   2.3s   get_twin_prime_iterator1
#   4.2s   get_twin_prime_iterator3
#   4.5s   get_twin_prime_iterator4 (object iterator)
#
# Iterators without precalc:
#   7.7s   get_twin_prime_iterator2
#   8.5s   get_twin_prime_iterator1
#  10.8s   get_twin_prime_iterator3
#  16.7s   get_twin_prime_iterator4 (object iterator)
#
# Alternatives:
# 229.6s   Math::NumSeq::TwinPrimes (Perl 5.19.7, Math::NumSeq 69)
#   6.6s - perl -MMath::PariInit=primes=65000000 -MMath::Pari=forprime,PARI -E
#          '$l=2;forprime($x,2,64764841,sub{say $l if $l+2==$x;$l=int("$x");});'

# This speeds things up, but isn't necessary.
#  Easy but estimates very high:
#my $estimate = 5000 + int( nth_prime_upper($count) * 1.4 * log($count) );
#  Relatively tight upper bound
#my $estimate = 400 + int(1.01 * nth_twin_prime_approx($count));
#  Simple and fastest:  Use the estimate directly.
my $estimate = nth_twin_prime_approx($count);
prime_precalc($estimate);

# Create a twin prime iterator using the prime_iterator construct
sub get_twin_prime_iterator1 {
  my $p = shift || 2;
  my $it = prime_iterator($p);
  my $prev = $it->();    # prev = 2
  $p = $it->();          # p = 3
  return sub {
    do {
      ($prev, $p) = ($p, $it->())
    } while ($p-$prev) != 2;
    $prev;
  };
}

# Create a twin prime iterator using the next_prime function
# A bit faster than the prime_iterator version.
sub get_twin_prime_iterator2 {
  my $start = shift || 2;
  my $p = next_prime($start-1);
  my $prev = next_prime($p);
  return sub {
    do {
      ($prev, $p) = ($p, next_prime($p))
    } while ($p-$prev) != 2;
    $prev;
  };
}

# Use Iterator::Simple
#use Iterator::Simple qw/igrep/;
#sub get_twin_prime_iterator3 {
#  my $start = shift || 2;
#  return igrep { is_prime($_+2) } prime_iterator($start);
#}

# Not very efficient, using object iterator and peek.
sub get_twin_prime_iterator4 {
  my $p = shift || 2;
  my $it = Math::Prime::Util::prime_iterator_object($p);
  $p = $it->value();
  return sub {
    $it->next() while $it->peek() - $it->value() != 2;
    $it->iterate();
  };
}

my $twinit = get_twin_prime_iterator2();
for (1..$count) {
  print $twinit->(), "\n";
}
