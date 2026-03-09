#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/next_calkin_wilf next_stern_brocot
                         calkin_wilf_n stern_brocot_n
                         nth_calkin_wilf nth_stern_brocot
                         nth_stern_diatomic
                         farey next_farey farey_rank/;

my $extended = $ENV{EXTENDED_TESTING};

my @CW = ([1,1],[1,2],[2,1],[1,3],[3,2],[2,3],[3,1],[1,4],[4,3],[3,5],[5,2],[2,5],[5,3],[3,4],[4,1],[1,5],[5,4],[4,7],[7,3],[3,8],[8,5],[5,7],[7,2],[2,7],[7,5],[5,8],[8,3],[3,7],[7,4],[4,5],[5,1],[1,6],[6,5],[5,9],[9,4],[4,11],[11,7],[7,10],[10,3],[3,11],[11,8],[8,13],[13,5],[5,12],[12,7],[7,9],[9,2],[2,9],[9,7],[7,12],[12,5],[5,13],[13,8],[8,11],[11,3],[3,10],[10,7],[7,11],[11,4],[4,9],[9,5],[5,6],[6,1],[1,7],[7,6],[6,11],[11,5],[5,14],[14,9],[9,13],[13,4],[4,15],[15,11],[11,18],[18,7],[7,17],[17,10],[10,13],[13,3],[3,14],[14,11],[11,19],[19,8],[8,21],[21,13],[13,18],[18,5],[5,17],[17,12],[12,19],[19,7],[7,16],[16,9],[9,11],[11,2],[2,11],[11,9],[9,16],[16,7],[7,19]);
my @SB = ([1,1],[1,2],[2,1],[1,3],[2,3],[3,2],[3,1],[1,4],[2,5],[3,5],[3,4],[4,3],[5,3],[5,2],[4,1],[1,5],[2,7],[3,8],[3,7],[4,7],[5,8],[5,7],[4,5],[5,4],[7,5],[8,5],[7,4],[7,3],[8,3],[7,2],[5,1],[1,6],[2,9],[3,11],[3,10],[4,11],[5,13],[5,12],[4,9],[5,9],[7,12],[8,13],[7,11],[7,10],[8,11],[7,9],[5,6],[6,5],[9,7],[11,8],[10,7],[11,7],[13,8],[12,7],[9,5],[9,4],[12,5],[13,5],[11,4],[10,3],[11,3],[9,2],[6,1],[1,7],[2,11],[3,14],[3,13],[4,15],[5,18],[5,17],[4,13],[5,14],[7,19],[8,21],[7,18],[7,17],[8,19],[7,16],[5,11],[6,11],[9,16],[11,19],[10,17],[11,18],[13,21],[12,19],[9,14],[9,13],[12,17],[13,18],[11,15],[10,13],[11,14],[9,11],[6,7],[7,6],[11,9],[14,11],[13,10],[15,11]);

plan tests => 4 + ($extended ? 1 : 0);

subtest 'Calkin-Wilf tree' => sub {
  my @data = (
    # n d idxCW
    [4, 11, 36],
    [22,7,519],
    [37,53,1990],
    [144,233,2730],
    [83116,51639,123456789],
    [64,65,"36893488147419103230"],
    [66,65,"36893488147419103233"],
    [32,1,4294967295],
    [64,1,"18446744073709551615"],
    ["228909276746","645603216423","1054982144710410407556"],
  );

  # next_calkin_wilf: generate first 100 terms
  my @s=([1,1]);
  push @s, [next_calkin_wilf($s[-1]->[0],$s[-1]->[1])] for 1..99;
  is_deeply( \@s, \@CW, "next_calkin_wilf first 100 terms" );

  # calkin_wilf_n: index of first 100 terms
  my @idx;
  push @idx, calkin_wilf_n($_->[0],$_->[1]) for @CW;
  is_deeply( \@idx, [1..100], "calkin_wilf_n first 100 terms" );

  # nth_calkin_wilf: recover first 100 terms from index
  my @nth;
  push @nth,[nth_calkin_wilf($_)] for 1..100;
  is_deeply( \@nth, \@CW, "nth_calkin_wilf first 100 terms" );

  # Roundtrip on selected larger values
  for my $t (@data) {
    my($n,$d,$idx) = @$t;
    is( "".calkin_wilf_n($n,$d), "$idx", "calkin_wilf_n($n,$d) = $idx" );
    is_deeply( [map{"$_"}nth_calkin_wilf($idx)], ["$n","$d"], "nth_calkin_wilf($idx) = ($n,$d)" );
  }
};

subtest 'Stern-Brocot tree' => sub {
  my @data = (
    # n d idxSB
    [4, 11, 36],
    [22,7,960],
    [37,53,1423],
    [144,233,2730],
    [83116,51639,111333227],
    [64,65,"27670116110564327423"],
    [66,65,"55340232221128654848"],
    [32,1,4294967295],
    [64,1,"18446744073709551615"],
    ["228909276746","645603216423","667408827216638861715"],
  );

  # next_stern_brocot: generate first 100 terms
  my @s=([1,1]);
  push @s, [next_stern_brocot($s[-1]->[0],$s[-1]->[1])] for 1..99;
  is_deeply( \@s, \@SB, "next_stern_brocot first 100 terms" );

  # stern_brocot_n: index of first 100 terms
  my @idx;
  push @idx, stern_brocot_n($_->[0],$_->[1]) for @SB;
  is_deeply( \@idx, [1..100], "stern_brocot_n first 100 terms" );

  # nth_stern_brocot: recover first 100 terms from index
  my @nth;
  push @nth,[nth_stern_brocot($_)] for 1..100;
  is_deeply( \@nth, \@SB, "nth_stern_brocot first 100 terms" );

  # Roundtrip on selected larger values
  for my $t (@data) {
    my($n,$d,$idx) = @$t;
    is( "".stern_brocot_n($n,$d), "$idx", "stern_brocot_n($n,$d) = $idx" );
    is_deeply( [map{"$_"}nth_stern_brocot($idx)], ["$n","$d"], "nth_stern_brocot($idx) = ($n,$d)" );
  }
};

subtest 'Stern diatomic (fusc)' => sub {
  my @A002487 = (0,1,1,2,1,3,2,3,1,4,3,5,2,5,3,4,1,5,4,7,3,8,5,7,2,7,5,8,3,7,4,5,1,6,5,9,4,11,7,10,3,11,8,13,5,12,7,9,2,9,7,12,5,13,8,11,3,10,7,11,4,9,5,6,1,7,6,11,5,14,9,13,4,15,11,18,7,17,10,13,3,14,11,19,8,21,13,18,5,17,12,19);
  my @fuscs = (
    # A212288 selected values
    [4691,257],
    [87339,2312],
    [1222997,13529],
    [9786539,57317],
    [76895573,238605],
    [357214891,744095],
    [1431655083,1948354],
    [5726623019,5102687],
    [22906492075,13354827],
    [91625925291,34961522],
  );

  # First terms of A002487
  my @s = map { nth_stern_diatomic($_) } 0 .. $#A002487;
  is_deeply( \@s, \@A002487, "A002487 first terms" );

  # Selected larger values
  for my $t (@fuscs) {
    is( nth_stern_diatomic($t->[0]), $t->[1], "fusc($t->[0]) = $t->[1]" );
  }

  # Edge cases
  is( nth_stern_diatomic(0), 0, "fusc(0) = 0" );
  is( nth_stern_diatomic(1), 1, "fusc(1) = 1" );

  # Recurrence: fusc(2n) = fusc(n), fusc(2n+1) = fusc(n) + fusc(n+1)
  for my $n (1..45) {
    is( nth_stern_diatomic(2*$n), nth_stern_diatomic($n), "fusc(2*$n) = fusc($n)" );
    is( nth_stern_diatomic(2*$n+1), nth_stern_diatomic($n) + nth_stern_diatomic($n+1), "fusc(2*$n+1) = fusc($n)+fusc($n+1)" );
  }
};

subtest 'Farey sequences' => sub {
  my @Farey = (
    undef,
    [[0,1],[1,1]],
    [[0,1],[1,2],[1,1]],
    [[0,1],[1,3],[1,2],[2,3],[1,1]],
    [[0,1],[1,4],[1,3],[1,2],[2,3],[3,4],[1,1]],
    [[0,1],[1,5],[1,4],[1,3],[2,5],[1,2],[3,5],[2,3],[3,4],[4,5],[1,1]],
    [[0,1],[1,6],[1,5],[1,4],[1,3],[2,5],[1,2],[3,5],[2,3],[3,4],[4,5],[5,6],[1,1]],
    [[0,1],[1,7],[1,6],[1,5],[1,4],[2,7],[1,3],[2,5],[3,7],[1,2],[4,7],[3,5],[2,3],[5,7],[3,4],[4,5],[5,6],[6,7],[1,1]],
    [[0,1],[1,8],[1,7],[1,6],[1,5],[1,4],[2,7],[1,3],[3,8],[2,5],[3,7],[1,2],[4,7],[3,5],[5,8],[2,3],[5,7],[3,4],[4,5],[5,6],[6,7],[7,8],[1,1]],
    [[0,1],[1,9],[1,8],[1,7],[1,6],[1,5],[2,9],[1,4],[2,7],[1,3],[3,8],[2,5],[3,7],[4,9],[1,2],[5,9],[4,7],[3,5],[5,8],[2,3],[5,7],[3,4],[7,9],[4,5],[5,6],[6,7],[7,8],[8,9],[1,1]],
  );

  # Full Farey sequences orders 1-9
  for my $n (1 .. $#Farey) {
    my @expf = @{$Farey[$n]};
    my @gotf = farey($n);
    my $gotlen = farey($n);
    my $explen = scalar(@expf);

    is( $gotlen, $explen, "scalar farey($n) = $explen" );
    is_deeply( \@gotf, \@expf, "farey($n)" );

    # Index access: farey(n,k) for each k
    my @gotf1 = map { farey($n,$_) } 0 .. $gotlen;
    is_deeply( \@gotf1, [@expf,undef], "farey($n,k) for k=0.." );

    # next_farey iteration
    my @gotf2 = ([0,1]);
    for (1..10000) {
      my $next = next_farey($n,$gotf2[-1]);
      last unless defined $next;
      push @gotf2, $next;
    }
    is_deeply( \@gotf2, \@expf, "next_farey($n,...) iteration" );

    # farey_rank roundtrip
    my @gotrank = map { farey_rank($n,$_) } @expf;
    is_deeply( \@gotrank, [0 .. $#expf], "farey_rank($n,...)" );
  }

  # Spot-check larger orders
  my @farey_ex = (
    [24, 16, [2,21]],
    [507,427, [3,505]],
  );
  for my $t (@farey_ex) {
    my($n,$k,$frac) = @$t;
    my $fracstr = join "/",@$frac;
    is_deeply( farey($n,$k), $frac, "farey($n,$k) = $fracstr" );
    is( farey_rank($n,$frac), $k, "farey_rank($n,[$fracstr]) = $k" );
  }

  # Edge cases
  is_deeply( [farey(1)], [[0,1],[1,1]], "farey(1)" );
  is( farey_rank(5,[0,1]), 0, "farey_rank(5,[0,1]) = 0" );
  is( farey_rank(5,[1,1]), 10, "farey_rank(5,[1,1]) = last" );
  ok( !defined(next_farey(5,[1,1])), "next_farey(5,[1,1]) = undef" );
};

if ($extended) {
  subtest 'extended tests' => sub {

    # Farey count for larger orders (|F_n| = 1 + sum_{k=1}^{n} euler_phi(k))
    my @farey_counts = (
      # n  |F_n|
      [10, 33],
      [20, 129],
      [50, 775],
      [100, 3045],
    );
    for my $t (@farey_counts) {
      my($n, $expcount) = @$t;
      is( farey($n), $expcount, "scalar farey($n) = $expcount" );
    }

    # next_farey iteration count for larger order
    {
      my $n = 50;
      my $count = 1;
      my $frac = [0,1];
      while (defined($frac = next_farey($n, $frac))) { $count++; }
      is( $count, 775, "next_farey($n,...) iterates 775 terms" );
    }

    # Stern diatomic: fusc(2^k) = 1 for all k
    for my $k (0, 1, 10, 20, 30, 40) {
      my $pow2 = 1 << $k;
      is( nth_stern_diatomic($pow2), 1, "fusc(2^$k) = 1" );
    }

    # Stern diatomic: fusc(2^k - 1) = k
    for my $k (1, 2, 10, 20, 30) {
      my $m = (1 << $k) - 1;
      is( nth_stern_diatomic($m), $k, "fusc(2^$k - 1) = $k" );
    }

    # farey_rank and farey(n,k) roundtrip for larger order
    {
      my $n = 100;
      for my $k (0, 1, 100, 1000, 3044) {
        my $frac = farey($n, $k);
        ok( defined $frac, "farey($n,$k) defined" );
        is( farey_rank($n, $frac), $k, "farey_rank($n, farey($n,$k)) = $k" );
      }
    }
  };
}
