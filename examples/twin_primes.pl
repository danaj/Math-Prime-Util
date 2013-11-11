#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/-nobigint
                         prime_iterator  prime_iterator_object
                         next_prime  is_prime
                         nth_prime_upper  prime_precalc/;

my $count = shift || 20;

# Find twin primes (numbers where p and p+2 are prime)

# Time for the first 300k:
#   1.0s   bin/primes.pl --twin 2 64764839
#   1.4s   get_twin_prime_iterator2
#   2.3s   get_twin_prime_iterator1
#   4.1s   get_twin_prime_iterator3
#   6.1s   get_twin_prime_iterator3 (object iterator)
#   7.6s   get_twin_prime_iterator2 without precalc
#   8.4s   get_twin_prime_iterator1 without precalc
#  10.9s   get_twin_prime_iterator3 without precalc
#  13.1s   get_twin_prime_iterator3 without precalc (object iterator)
# 219.8s   Math::NumSeq::TwinPrimes (Perl 5.19.4 with v66)

# This speeds things up, but isn't necessary.
my $estimate = 5000 + int( nth_prime_upper($count) * 1.4 * log($count) );
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
