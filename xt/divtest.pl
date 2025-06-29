#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use Math::Prime::Util::PP;
use Math::BigInt;


my @N = (1..255);
push @N, map { (1 << $_)+0 } 8 .. 63;
push @N, map { (1 << $_)-1 } 8 .. 63;
push @N, map { (1 << $_)+1 } 8 .. 63;
push @N, ~0;
push @N, ~0 - 1;
push @N, 479,709,1979,2293,5527,9491,18803,55697,67547,234977,472189,794831,1795987,3420581,6659201,8616679,32359207,36287963,110125493,186431731,329522393,579566597,1081042813,3930325367,8232694003,9965227067,18953892493,48398125001,74180605321,172260332171,531942628597,952424752493,1683828392317,4152728744257,7544680499843,11684892744529,29481476832221,44483282520737,80989922766787,165916886096329,372528286807679,625954070066009,2024742834179983,3492450230167973,8886306409922317,9866050790952803,35085838827392533,37974793647167711,79879452462015781,206179373730094717,508869317949427363,875096421074592361,2064451162885400581,4422949619687292341,4629427415289143573,9405019501832426699;
push @N, 2232881419280027;  # Crash and burn...

my @negN = map { negint($_) } @N;

print "running simple division tests\n";
for my $num (@N, @negN) {
  for my $den (@N, @negN) {

    #next unless $num >= 0 && $den >= 0;
    #print "$num $den\n";

    # This can't be represented as an IV or UV, so skip it
    next if $den < 0 && ($num >> 63) > 0;

    my $mx = (tdivrem($num,$den))[0];

    # These two always match
    my $bi = Math::BigInt->new($num)->btdiv($den);
    die "$num $den" if $bi ne $mx;
    my $mp = (Math::Prime::Util::PP::tdivrem($num,$den))[0];
    die "$num $den" if "$mp" ne "$mx";

    # This does not
    #my $nd = (ndivrem($num,$den))[0];
    #die "$num $den" if $nd ne $mx;

   # my $nd = div9316n($num,$den);
   # die "$num $den" if $nd ne $mx;
   # print "$num $den = $nd\n" if $den == 2 && $num == 2232881419280027;
   # my $dm = (divmod9316n($num,$den))[0];
   # print "$num $den = $dm\n" if $den == 2 && $num == 2232881419280027;
   # warn "$num $den" if $dm ne $mx;
  }
}
print "pass simple tests\n";



sub ndivrem {
  my($D,$d) = @_;
  ( int($D/$d), $D % $d );
}
sub uidivrem {
  my($D,$d) = @_;
  use integer;
  ( int($D/$d), $D % $d );
}
sub nfdivrem {
  my($D,$d) = @_;
  use POSIX;
  ( POSIX::floor($D/$d), $D % $d );
}
sub div_and_mod {
  my($D,$d) = @_;
  ( divint($D,$d), modint($D,$d) );
}
sub divmod9316 {
  my($D,$d) = @_;
  my $mod = $D % $d;
  (($D - $mod) / $d, $mod);
}

# This is SO close to working.  But ... not always.
sub divmod9316n {   # Truncated div and mod
  my($D,$d) = @_;
  my $mod = $D % $d;
  $mod -= $d if $mod != 0 && (($D < 0 && $d >= 0) || ($D >= 0 && $d < 0));
  (($D - $mod) / $d, $mod);
}
sub div9316 {
  my($D,$d) = @_;
  ($D - ($D % $d)) / $d;
}
sub div9316n {
  my($D,$d) = @_;

  #return -div9316(-$D, $d) if $D < 0 && $d >= 0;
  #return -div9316( $D,-$d) if $D >= 0 && $d < 0;
  #return  div9316(-$D,-$d) if $D < 0 && $d < 0;
  #return  div9316( $D, $d);

  # Truncated div and mod
  my $mod = $D % $d;
  $mod -= $d if $mod != 0 && (($D < 0 && $d >= 0) || ($D >= 0 && $d < 0));
  return ($D - $mod) / $d;
}

sub ivmod {
  my($a,$n) = @_;
  die "wrong usage: n must be positive" unless $n >= 0;
  return $a % $n if $a >= 0;
  my $amodn = -$a % $n;
  return ($amodn == 0) ? 0 : $n-$amodn;
}


# div9316n wrong for things like:
#
#   1125899906842624 1
#   2232881419280027 2
#   4503599627370496 3
#   4503599627370496 4
#   9007199254740992 6
#   9007199254740993 9
#   3492450230167973 3
