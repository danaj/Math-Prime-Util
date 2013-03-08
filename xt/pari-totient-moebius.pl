#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/-nobigint/;
use Math::Pari;

my $nlinear = 100000;
my $nrandom = shift || 100000;
my $randmax = 10**16;

# Moebius and euler_phi seem about 2-4x faster than Pari.  Also, we have
# ranged versions that run much faster.
#
print "OK for first 1";
my $dig = 1;
my $i = 9;
foreach my $n (2 .. $nlinear) {
  die "failure for eulerphi($n)" unless Math::Prime::Util::euler_phi($n) == Math::Pari::eulerphi($n);
  die "failure for moebius($n)" unless Math::Prime::Util::moebius($n) == Math::Pari::moebius($n);
  if (--$i == 0) {
    print "0";
    $dig++;
    $i = (10 ** $dig) - (10 ** ($dig-1));
  }
}
print " numbers\n";
print "Testing random numbers from $nlinear to ", $randmax, "\n";
my $mod = int($nrandom / 80);

while ($nrandom-- > 0) {
  my $n = $nlinear + 1 + int(rand($randmax - $nlinear));
  die "failure for eulerphi($n)" unless Math::Prime::Util::euler_phi($n) == Math::Pari::eulerphi($n);
  die "failure for moebius($n)" unless Math::Prime::Util::moebius($n) == Math::Pari::moebius($n);
  print "." if ($nrandom % $mod) == 0;
}
print "\n";
