#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/prime_iterator is_prime
                         next_prime nth_prime_upper prime_precalc forprimes/;

my $count = shift || 20;
my $method = shift || 'forprimes';
my $precalc = 0;   # If set, precalc all the values we'll call is_prime on

# Find Sophie Germain primes (numbers where p and 2p+1 are both prime).

# Four methods are shown: forprimes, iter, iter2, and MNS.

# Times for 300k:
#
#                300k               1M
# precalc:
#   forprimes    1.3s   9.0MB       7.1s   21.6MB
#   iter         2.8s   8.7MB      12.6s   21.4MB
#   iter2        1.9s   8.7MB       9.4s   21.4MB
# no precalc:
#   forprimes    1.5s   4.5MB       5.6s   4.5MB
#   iter         9.5s   4.3MB      37.5s   4.3MB
#   iter2        8.5s   4.3MB      33.9s   4.3MB
#   MNS        254.3s  11.3MB   >1500s   >15 MB

if ($precalc) {
  prime_precalc(2 * sg_upper_bound($count));
}


if ($method eq 'forprimes') {

  my $estimate = sg_upper_bound($count);
  my $numfound = 0;
  forprimes {
    if ($numfound < $count && is_prime(2*$_+1)) {
      print "$_\n";
      $numfound++;
    }
  } $estimate;
  die "Estimate too low" unless $numfound >= $count;

} elsif ($method eq 'iter') {

  # Wrap the standard iterator
  sub get_sophie_germain_iterator {
    my $p = shift || 2;
    my $it = prime_iterator($p);
    return sub {
      do { $p = $it->() } while !is_prime(2*$p+1);
      $p;
    };
  }
  my $sgit = get_sophie_germain_iterator();
  print $sgit->(), "\n" for 1 .. $count;

} elsif ($method eq 'iter2') {

  # Iterate directly using next_prime
  my $prime = 2;
  for (1 .. $count) {
    $prime = next_prime($prime) while !is_prime(2*$prime+1);
    print "$prime\n";
    $prime = next_prime($prime);
  }

} elsif ($method eq 'MNS') {

  # Use Math::NumSeq
  require Math::NumSeq::SophieGermainPrimes;
  my $seq = Math::NumSeq::SophieGermainPrimes->new; 
  for (1 .. $count) {
    print 0+($seq->next)[1];
  }

}

# Used for precalc and the forprimes example
sub sg_upper_bound {
  my $count = shift;
  my $nth = nth_prime_upper($count);
  # For lack of a better formula, do this step-wise estimate.
  my $estimate = ($count <   5000) ? 150 + int( $nth * log($nth) * 1.2 )
               : ($count <  19000) ? int( $nth * log($nth) * 1.135 )
               : ($count <  45000) ? int( $nth * log($nth) * 1.10 )
               : ($count < 100000) ? int( $nth * log($nth) * 1.08 )
               : ($count < 165000) ? int( $nth * log($nth) * 1.06 )
               : ($count < 360000) ? int( $nth * log($nth) * 1.05 )
               : ($count < 750000) ? int( $nth * log($nth) * 1.04 )
               : ($count <1700000) ? int( $nth * log($nth) * 1.03 )
               :                     int( $nth * log($nth) * 1.02 );

  return $estimate;
}
