#!/usr/bin/env perl
use warnings;
use strict;
use v5.16;
use ntheory ":all";
use Math::Prime::Util::PP;

#csrand(2);
#for (1..100) {
#  my $k = 10+urandomm(10000);
#  test($k);
#}

# Performance:
#    my @r = allrootmod(33,3432,10428581733134514527);

test($_,10000,10** 4) for 2,3,5,7,4,6,8,9,27,625;
test($_,10000,10** 8) for 2,3,5,7,4,6,8,9,27,625;
test($_,10000,10**12) for 2,3,5,7,4,6,8,9,27,625;
test($_,10000,10**16) for 2,3,5,7,4,6,8,9,27,625;
test($_,10000,~0)     for 2,3,5,7,4,6,8,9,27,625;

sub test {
  my($k,$iter,$size) = @_;

  for (1..$iter) {
    my(@rxs,@rpp, @sxs,@spp, $a,$n);
    $n = urandomm(100_000_000_000);
    $a = urandomm($n);
    @rxs = allrootmod($a,$k,$n);
    @rpp = Math::Prime::Util::PP::allrootmod($a,$k,$n);
    die "allrootmod fail for $a,$k,$n\n" unless vecequal(\@rxs,\@rpp);
    if ($k == 2) {
      @sxs = allsqrtmod($a,$n);
      @spp = Math::Prime::Util::PP::allsqrtmod($a,$n);
      die "allsqrtmod fail for $a,$n\n" unless vecequal(\@sxs,\@spp);
      die "allsqrtmod fail for $a,$n\n" unless vecequal(\@rxs,\@sxs);
    }

    if (@rxs) {
      my %roots = map { $_ => 1 } @rxs;
      my $rxs = rootmod($a,$k,$n);
      my $rpp = Math::Prime::Util::PP::rootmod($a,$k,$n);
      die "xs rootmod fail for $a,$k,$n" unless $roots{$rxs};
      die "pp rootmod fail for $a,$k,$n" unless $roots{$rpp};
      if ($k == 2) {
        my $sxs = sqrtmod($a,$n);
        my $spp = Math::Prime::Util::PP::sqrtmod($a,$n);
        die "xs sqrtmod fail for $a,$n" unless $roots{$sxs};
        die "pp sqrtmod fail for $a,$n" unless $roots{$spp};
      }
    }
  }
  say "pass $iter iterations with max n $size for k $k";
}
