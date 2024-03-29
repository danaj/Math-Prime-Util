#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ gcd lcm gcdext /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @gcds = (
  [ [], 0],
  [ [8], 8],
  [ [9,9], 9],
  [ [0,0], 0],
  [ [1, 0, 0], 1],
  [ [0, 0, 1], 1],
  [ [17,19], 1 ],
  [ [54,24], 6 ],
  [ [42,56], 14],
  [ [ 9,28], 1 ],
  [ [48,180], 12],
  [ [2705353758,2540073744,3512215098,2214052398], 18],
  [ [2301535282,3609610580,3261189640], 106],
  [ [694966514,510402262,195075284,609944479], 181],
  [ [294950648,651855678,263274296,493043500,581345426], 58 ],
  [ [-30,-90,90], 30],
  [ [-3,-9,-18], 3],
  [ [-5], 5],
  [ [-5,5], 5],
  [ [-5,7], 1],
);
my @lcms = (
  [ [], 1],
  [ [8], 8],
  [ [9,9], 9],
  [ [0,0], 0],
  [ [1, 0, 0], 0],
  [ [0, 0, 1], 0],
  [ [17,19], 323 ],
  [ [54,24], 216 ],
  [ [42,56], 168],
  [ [ 9,28], 252 ],
  [ [48,180], 720],
  [ [36,45], 180],
  [ [-36,45], 180],
  [ [-36,-45], 180],
  [ [30,15,5], 30],
  [ [2,3,4,5], 60],
  [ [30245, 114552], 3464625240],
  [ [11926,78001,2211], 2790719778],
  [ [1426,26195,3289,8346], 4254749070],
  [ [-5], 5],
  [ [-5,5], 5],
  [ [-5,7], 35],
);
if ($use64) {
  push @gcds, [ [12848174105599691600,15386870946739346600,11876770906605497900], 700];
  push @gcds, [ [9785375481451202685,17905669244643674637,11069209430356622337], 117];
  push @lcms, [ [26505798,9658520,967043,18285904], 15399063829732542960];
  push @lcms, [ [267220708,143775143,261076], 15015659316963449908];
}

my @gcdexts = (
  [ [0,  0], [0, 0, 0] ],
  [ [0, 28], [0, 1,28] ],
  [ [ 28,0], [ 1,0,28] ],
  [ [0,-28], [0,-1,28] ],
  [ [-28,0], [-1,0,28] ],
  [ [ 3706259912, 1223661804], [ 123862139,-375156991, 4] ],
  [ [ 3706259912,-1223661804], [ 123862139, 375156991, 4] ],
  [ [-3706259912, 1223661804], [-123862139,-375156991, 4] ],
  [ [-3706259912,-1223661804], [-123862139, 375156991, 4] ],
  [ [22,242], [1, 0, 22] ],
  [ [2731583792,3028241442], [-187089956, 168761937, 2] ],
  [ [42272720,12439910], [-21984, 74705, 70] ],
);
if ($use64) {
  push @gcdexts, [ [10139483024654235947,8030280778952246347], [-2715309548282941287,3428502169395958570,1] ];
}

plan tests => scalar(@gcds) + scalar(@lcms) + scalar(@gcdexts);

###### gcd
foreach my $garg (@gcds) {
  my($aref, $exp) = @$garg;
  my $gcd = gcd(@$aref);
  is( $gcd, $exp, "gcd(".join(",",@$aref).") = $exp" );
}
###### lcm
foreach my $garg (@lcms) {
  my($aref, $exp) = @$garg;
  my $lcm = lcm(@$aref);
  is( $lcm, $exp, "lcm(".join(",",@$aref).") = $exp" );
}
###### gcdext
foreach my $garg (@gcdexts) {
  my($aref, $eref) = @$garg;
  my($x,$y) = @$aref;
  is_deeply( [gcdext($x,$y)], $eref, "gcdext($x,$y) = [@$eref]" );
}
