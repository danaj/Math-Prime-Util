#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/euler_phi/;
use List::Util qw/sum/;

my $limit = shift || 1_000_000;

print "Calculating totients from 1 to $limit...";
my @phi = map { euler_phi($_) } 1 .. $limit;
print "...";
unshift @phi, 0;
print "...done\n";

print "Running non-stop random tests.  Break when desired.\n";
while (1) {
  my $beg = 0 + int(rand($limit));
  my $end = 0 + int(rand($limit));
  ($beg,$end) = ($end,$beg) if $beg > $end;

  # Does range return the same values?
  my @phi_range = @phi[ $beg .. $end ];
  my @totients = euler_phi($beg,$end);

  my $sum1 = sum(@phi_range);
  my $sum2 = sum(@totients);
  warn "\nbeg $beg  end $end  sum $sum1  range sum $sum2\n"
       unless $sum1 == $sum2;
  print ".";
}
