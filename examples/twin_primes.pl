#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/prime_iterator  prime_iterator_object
                         next_prime  is_prime
                         nth_prime_upper  nth_twin_prime_approx  prime_precalc/;

my $count = shift || 20;

# Find twin primes (numbers where p and p+2 are prime)

# Time for the first 1M:
#
# Not iterators:
#   0.4s   say join "\n", @{twin_primes(2,nth_twin_prime(1e6))}
#   1.3s   $l=2; forprimes { say $l if $l+2==$_; $l=$_; } 2+nth_twin_prime(1e6)
#   0.4s   bin/primes.pl --twin 2 252427601
#
# Iterators with precalc:
#   4.5s   get_twin_prime_iterator2 (next_prime)
#   5.4s   get_twin_prime_iterator1 (prime_iterator)
#   9.4s   get_twin_prime_iterator3 (Iterator::Simple)
#  13.8s   get_twin_prime_iterator4 (object iterator)
#
# Iterators without precalc:
#  11.6s   get_twin_prime_iterator2
#   5.3s   get_twin_prime_iterator1
#   9.3s   get_twin_prime_iterator3
#  28.0s   get_twin_prime_iterator4 (object iterator)
#
# Alternative iterator:
# 3944.4s  Math::NumSeq::TwinPrimes (Perl 5.27.2, Math::NumSeq 72)
#
# Alternative non-iterators:
#   14.5s  perl -MMath::PariInit=primes=255000000 -MMath::Pari=forprime,PARI -E
#         '$l=2;forprime($x,2,252427603,sub{say $l if $l+2==$x;$l=int("$x");});'
#    4.7s  perl -MMath::Prime::FastSieve -E 'my $s=Math::Prime::FastSieve::Sieve->new(255000000); for my $p (@{$s->primes(252427601)}) { say $p if $s->isprime($p+2); }'

# This speeds things up, but isn't necessary.
#  Easy but estimates very high:
#my $estimate = 5000 + int( nth_prime_upper($count) * 1.4 * log($count) );
#  Relatively tight upper bound
#my $estimate = 800 + int(1.01 * nth_twin_prime_approx($count));
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
