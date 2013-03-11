#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

my $limit = shift || 20_000;

use Math::Prime::Util qw/:all/;
# Use another code base for comparison.
use Math::Primality;

foreach my $n (0 .. $limit) {
  die "next $n" unless next_prime($n) == Math::Primality::next_prime($n);
  if ($n <= 2) {
    die "prev $n" unless prev_prime($n) == 0;
  } else {
    die "prev $n" unless prev_prime($n) == Math::Primality::prev_prime($n);
  }
  die "is $n" unless is_prime($n) == Math::Primality::is_prime($n);
  print "$n.." unless $n % 10000;
}
print "Success to $limit!\n";
