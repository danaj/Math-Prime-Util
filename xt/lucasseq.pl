#!/usr/bin/perl
use warnings;
use strict;
use ntheory ":all";
use v5.030;
use Math::Prime::Util::PP;

# These two functions based on pseudocode from Trizen.
# We will compare the library code vs. this.

# Combine V_k(P,Q) mod n
sub _lucasvmod {
  my($P,$Q, $k, $n) = @_;

  my($V1,$V2) = (2,modint($P,$n));
  my($Q1,$Q2) = (1,1);

  my(@bits) = todigits($k,2);
  for my $bit (@bits) {
    # 5 even, 6 odd
    $Q1 = mulmod($Q1,$Q2,$n);
    if ($bit) {
      $Q2 = mulmod($Q1,$Q,$n);
      $V1 = submod(mulmod($V2,$V1,$n),mulmod($P,$Q1,$n),$n);
      $V2 = submod(mulmod($V2,$V2,$n),mulmod(2,$Q2,$n),$n);
    } else {
      $Q2 = $Q1;
      $V2 = submod(mulmod($V2,$V1,$n),mulmod($P,$Q1,$n),$n);
      $V1 = submod(mulmod($V1,$V1,$n),mulmod(2,$Q2,$n),$n);
    }
  }
  # Outputs V_k, <...>, Q_k
  ($V1, $V2, mulmod($Q1,$Q2,$n));
}
  
sub _lucasuvmod {
  my($P,$Q, $k, $n) = @_;
  return (0,0,0) if $n == 1;
  return (0,2 % $n,1) if $k == 0;
  die "Invalid modulus $n" if $n <= 0;

  my $U1 = 1;
  my($V1,$V2) = (2 % $n,modint($P,$n));
  my($Q1,$Q2) = (1,1);
  my $D = submod( mulmod($P,$P,$n), mulmod(4,$Q,$n), $n);

  if (gcd($D,$n) == 1) {
    ($V1,$V2,$Q1) = _lucasvmod($P, $Q, $k, $n);
    $U1 = divmod( submod(mulmod(2,$V2,$n),mulmod($P,$V1,$n),$n), $D,$n);
    return ($U1, $V1, $Q1);
  }

  my $s = valuation($k,2);
  $k = rshiftint($k,$s+1);
  #print "s $s  k $k\n";
  my(@bits) = todigits($k,2);

  for my $bit (@bits) {
    # 3 even, 7 odd  (primality.c)
    # 6 even, 7 odd
    $Q1 = mulmod($Q1,$Q2,$n);
    #print "bit $bit Q1 = $Q1\n";
    if ($bit) {
      $Q2 = mulmod($Q1,$Q,$n);
      $U1 = mulmod($U1,$V2,$n);
      $V1 = submod(mulmod($V2,$V1,$n),mulmod($P,$Q1,$n),$n);
      $V2 = submod(mulmod($V2,$V2,$n),mulmod(2,$Q2,$n),$n);
    } else {
      $Q2 = $Q1;
      $U1 = submod(mulmod($U1,$V1,$n),$Q1,$n);
      $V2 = submod(mulmod($V2,$V1,$n),mulmod($P,$Q1,$n),$n);
      $V1 = submod(mulmod($V1,$V1,$n),mulmod(2,$Q2,$n),$n);
    }
  }
  $Q1 = mulmod($Q1,$Q2,$n);
  $Q2 = mulmod($Q1, $Q, $n);
  $U1 = submod(mulmod($U1,$V1,$n),$Q1,$n);
  $V1 = submod(mulmod($V2,$V1,$n),mulmod($P,$Q1,$n),$n);
  $Q1 = mulmod($Q1,$Q2,$n);

  for (1 .. $s) {
    $U1 = mulmod($U1,$V1,$n);
    $V1 = submod(mulmod($V1,$V1,$n),mulmod(2,$Q1,$n),$n);
    $Q1 = mulmod($Q1,$Q1,$n);
  }
  ($U1,$V1,$Q1);
}

#say join ", ",lucasuvmod(-4,4,50,1001);
#say join ", ",lucasuvmod(-4,7,50,1001);
#say join ", ",lucasuvmod(1,-1,50,1001);
#say join ", ",lucasuvmod(1,-1,4,5);

for my $n (1 .. 20) {    # n 1
  print "n $n\n";
  for my $k (0 .. 101) {  # k 0
    for my $P (-30 .. 30) {
      for my $Q (-30 .. 30) {
        #print "($n,$P,$Q,$k)\n";
        #my $s1 = join " ", (lucasuvmod($P, $Q, $k, $n))[0,1];
        #my $s2 = join " ", (Math::Prime::Util::GMP::lucas_sequence($n,$P,$Q,$k))[0,1];
        my $s1 = join " ", _lucasuvmod($P, $Q, $k, $n);
        #my $s2 = join " ", Math::Prime::Util::GMP::lucas_sequence($n,$P,$Q,$k);
        my $s3 = join " ", lucasuvmod($P,$Q,$k,$n);
        #my $s3 = join " ", Math::Prime::Util::PP::lucasuvmod($P,$Q,$k,$n);
        #say "($P,$Q,$k,$n) :  '$s1'  :  '$s3'"  if $s1 ne $s3;

        #Math::Prime::Util::GMP::lucas_sequence($n,$P,$Q,$k);
        #lucasuvmod($P, $Q, $k, $n);
        #lucas_sequence($n,$P,$Q,$k);
        #modint(lucasu($P,$Q,$k),$n); modint(lucasv($P,$Q,$k),$n); powmod($Q, $k, $n);
      }
    }
  }
}

#   2.8s  XS
#  10.1s  GMP
#  34.1s  funcs here (w/ XS)
# 270.4s  PP (w/ XS)
# 299.0s  PP (no XS)
# 803.3s  funcs here (no XS)
