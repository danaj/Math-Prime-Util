#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/euler_phi vecsum urandomm/;

my $limit = shift || 1_000_000;

print "Calculating totients from 1 to $limit...";
my @phi = map { euler_phi($_) } 1 .. $limit;
print "...";
unshift @phi, 0;
print "...done\n";

print "Running non-stop random tests.  Break when desired.\n";
while (1) {
  my $beg = urandomm($limit);
  my $end = urandomm($limit);
  ($beg,$end) = ($end,$beg) if $beg > $end;

  # Does range return the same values?
  my $sum1 = vecsum(  @phi[ $beg .. $end ]  );
  my $sum2 = vecsum(  euler_phi($beg,$end)  );

  warn "\nbeg $beg  end $end  sum $sum1  range sum $sum2\n"
       unless $sum1 == $sum2;
  print ".";
}
