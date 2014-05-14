#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/moebius mertens vecsum/;

my $limit = shift || 1_000_000;

print "Calculating moebius from 1 to $limit...";
my @mu = map { moebius($_) } 1 .. $limit;
print "...";
unshift @mu, 0;
print "...done\n";

while (1) {
  my $beg = 1 + int(rand($limit));
  my $end = 1 + int(rand($limit));
  ($beg,$end) = ($end,$beg) if $beg > $end;

  # Does moebius range return the same values?
  my @mu_range = @mu[ $beg .. $end ];
  my @mobius = moebius($beg,$end);

  my $mu_sum = vecsum(@mu_range);
  my $mo_sum = vecsum(@mobius);
  my $mert_sum = mertens($end) - mertens($beg-1);
  warn "\nbeg $beg  end $end  sum $mu_sum  range sum $mo_sum\n"
       unless $mu_sum == $mo_sum;
  warn "\nbeg $beg  end $end  sum $mu_sum  mertsum $mert_sum\n"
       unless $mu_sum == $mert_sum;
  print ".";
}
