#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;
use Math::GMPz;

# Sum of all n where is_power(euler_phi(n^2),3) = 1

# Simple but very slow way.  The brute force method later in this file is
# basically the same thing, but using the more efficient ranged moebius and
# totient calls over intervals.
#
#  Pari:
#    s=0; for(n=2,limit,if(ispower(n*eulerphi(n),3),s=s+n)); print(s)
#  Perl/MPU:
#    my $s=0;
#    for my $n (2..$limit) { $s += $n if is_power($n*euler_phi($n),3); }
#    say $s;
#
# TIMING:
#              10^7    2*10^7   10^8     10^10
# Clever       0.07s    0.09s   0.24s    5s
# Brute        5.7s    11.0s   53.9s     5 hours
# Simple MPU  12.0s    27.3s  174s       1-2 days?
# Simple Pari 13.6s    33.4s  277s       5 days?
#   

my $limit = shift || 10**10-1;
my $method = lc(shift || 'clever');
die "Method must be 'clever' or 'brute'\n" unless $method =~ /^(clever|brute)$/;
my $sum = 0;


if ($method eq 'clever') {
  # About 5 seconds for 10^10-1
  my $cblimit = int( ($limit*$limit) ** 0.3334 + 0.01 );
  foreach my $k (2 .. $cblimit) {
    next if $k & 1;
    my($p, $e) = @{ (factor_exp($k))[-1] };
    $e *= 3;
    next unless $e & 1;
    my $m = int($k / ($p ** int($e/3)));
    $m **= 3;
    next if $m % ($p-1);
    $m = int($m / ($p-1));
    my $n = $p ** (($e+1) >> 1);
    next if $n >= $limit;
    while ($m > 1) {
      my ($p,$e) = @{ (factor_exp($m))[-1] };
      last unless $e & 1;
      last if $m % ($p-1);
      $n *= $p ** (($e+1) >> 1);
      last if $n >= $limit;
      $m = int($m / ( ($p-1) * ($p**$e) ) );
    }
    if ($m == 1) {
      #print "$n\n";
      $sum += $n;
    }
  }
} else {
  # About 5 hours for 10^10-1
  my $interval = 10_000_000;   # Window size for moebius/totient
  #prime_precalc(10**9);        # Slightly faster ranged phi
  my($beg,$end) = (0,0);
  while ($beg < $limit) {
    $end = $beg + $interval - 1;
    $end = $limit if $end > $limit;
    my $start = ($beg<2)?2:$beg;

    my $glim = int(~0 / $end);
    my @m = moebius($beg, $end);
    my @t = euler_phi($beg, $end);

    if ($end <= $glim) {   # Totient($n) * $n will always be < ~0
      foreach my $n ($start .. $end) {
        next unless $m[$n-$beg] == 0;
        my $totn2 = $n * $t[$n-$beg];
        if (is_power($totn2,3)) {
          # print "$n\n";
          $sum += $n
        }
      }
    } else {
      foreach my $n ($start .. $end) {
        next unless $m[$n-$beg] == 0;
        my $tot = $t[$n-$beg];
        if ($tot <= $glim) {
          print "$n\n" if is_power($n * $tot, 3);
        } else {
          $tot = Math::GMPz->new($n) * $tot;
          print "$n\n" if Math::GMPz::Rmpz_perfect_power_p($tot) && is_power($tot,3);
        }
      }
    }
    $beg = $end+1;
  }
}
print "$sum\n";
