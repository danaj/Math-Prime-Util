#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

# x+y = a^2    x = a^2 - y
# x-y = b^2    a^2-y-y = b^2  2y = b^2-a^2  y = (b^2-a^2)/2
# x+z = c^2    z = c^2 - z
# x-z = d^2    c^2-z-z=d^2    2z = c^2-d^2  z = (c^2-d^2)/2
# y+z = e^2
# y-z = f^2                                 x = (e^2-f^2)/2

# x+y = a^2   x-y = b^2     ===>   2x = a^2+b^2   x=(a^2+b^2)/2
# x+z = c^2   x-z = d^2     ===>   2z = c^2-d^2   z=(c^2-d^2)/2
# y+z = e^2   y-z = f^2     ===>   2y = e^2+f^2   y=(e^2+f^2)/2

# a^2 = x+y = x+y+z-z = x+z + y-z  = c^2 + f^2
# e^2 = y+z = y+z+x-x = y+x -(x-z) = a^2 - d^2
# b^2 = x-y = x-y+z-z = x+z -(y+z) = c^2 - e^2

foreach my $a (4 .. 1000000) {
  my $a2 = $a*$a;
  foreach my $c (3 .. $a-1) {
    my $c2 = $c*$c;
    my $f2 = $a2 - $c2;
    next unless $f2 >= 0 && is_power($f2,2);
    foreach my $d (1 .. $c-1) {
      next if ($d ^ $c) & 1;   # c and d must have same parity
      my $d2 = $d*$d;
      my $e2 = $a2 - $d2;
      my $b2 = $c2 - $e2;
      next if $e2 <= 0 || $b2 <= 0;
      #next if (($a2+$b2) & 1) || (($e2+$f2) & 1) || (($c2-$d2) & 1);
      next unless is_power($e2,2) && is_power($b2,2);
      my $x = ($a2+$b2) >> 1;
      my $y = ($e2+$f2) >> 1;
      my $z = ($c2-$d2) >> 1;
      my $result = $x+$y+$z;
      die "$result  [$x $y $z]\n";
    }
  }
}
