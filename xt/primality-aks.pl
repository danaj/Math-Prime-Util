#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

my $limit = shift || 10_000_000;

use Math::Prime::Util qw/is_aks_prime/;
use Math::Prime::FastSieve;
my $sieve = Math::Prime::FastSieve::Sieve->new($limit + 10_000);

if (1) {
  my $n = 2;
  while ($n <= $limit) {
    print "$n\n" if $n > 69000; # unless $i++ % 262144;
    die "$n should be prime" unless is_aks_prime($n);
    my $next = $sieve->nearest_ge( $n+1 );
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
