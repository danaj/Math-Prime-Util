#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

my $limit = shift || 10_000_000;

use Math::Prime::Util qw/is_aks_prime next_prime/;

# It doesn't really matter for this use, but try to use independent code for
# next_prime.  Math::Prime::FastSieve works well.
my $sieve;
if (eval { require Math::Prime::FastSieve; Math::Prime::FastSieve->import(); require Math::Prime::FastSieve::Sieve; 1; }) {
  $sieve = Math::Prime::FastSieve::Sieve->new($limit + 10_000);
}

if (1) {
  my $n = 2;
  my $i = 1;
  while ($n <= $limit) {
    print "$n\n" unless $i++ % 1000;
    die "$n should be prime" unless is_aks_prime($n);
    my $next = (defined $sieve) ? $sieve->nearest_ge( $n+1 ) : next_prime($n);
    my $diff = ($next - $n) >> 1;
    if ($diff > 1) {
      foreach my $d (1 .. $diff-1) {
        my $cn = $n + 2*$d;
        die "$cn should be composite" if is_aks_prime($cn);
      }
    }
    $n = $next;
  }
  print "Success to $limit!\n";
}
