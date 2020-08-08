#!/usr/bin/perl
use warnings;
use strict;
use ntheory qw/binomialmod next_prime/;

{
  my $s = 0;
  my $p = 2;
  for my $n (1..4000) {
    for my $k (1..$n) {
      $b=binomialmod($n,$k,$p);
      $s+=$b;
    }
  }
  die "small binomial mod 2" if $s != 465296;
  print "   ok binomialmod(1..4000, 1..n, 2)\n";
}

{
  my $s = 0;
  for my $p (1..120) {
    for my $n (1..100) {
      for my $k (1..$n) {
        $b=binomialmod($n,$k,$p);
        $s += $b;
      }
    }
  }
  die "small binomial mod 1..120" if $s != 13535875;
  print "   ok binomialmod(1..100, 1..n, 1..120)\n";
}

{
  my @expect = (qw/51200 90431 214840 396378 855185 1038432 1316041 1565620 2185172 2879324 3794082 4866224 5227689 5568658 6050631 6944996 7792862 7766337 8805047 9303431 9470898 10240897 11059095 11588935 12579194 13417658 13460749 13817518 14358681 14906734/);
  my $p = 2;
  for (0..$#expect) {
    my $s=0;
    for my $n (1..1000) {
      for my $k (1..$n) {
        $s += binomialmod($n,$k,$p);
        #$s += Math::Prime::Util::binomial($n,$k) % $p;
      };
    }
    die "wrong binomialmod(1..1000,1..n,$p)\n" unless $s == $expect[$_];
    print "   ok binomialmod(1..1000, 1..n, $p)\n";
    $p=next_prime($p);
  }
}
