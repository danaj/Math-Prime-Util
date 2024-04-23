#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_powerful powerful_count nth_powerful sumpowerful
                         powerful_numbers factor_exp vecall vecnone/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

plan tests => 2+4+10+4+2+1       # is_powerful
            + 10 + 2*$extra      # powerful_count
            + 4                  # nth_powerful
            + 7                  # sumpowerful
            + 10;                # powerful_numbers

{
  my @exp = map { fac_is_powerful($_, 2) } 0 .. 258;
  is_deeply( [map { is_powerful($_,2) } 0..258], \@exp, "is_powerful(0..258,2)");
  is_deeply( [map { is_powerful($_) } 0..258], \@exp, "is_powerful(0..258)");
}

ok( vecnone(sub { is_powerful(-8,$_) }, 0..10), "is_powerful(-8,n) = 0");
ok( vecnone(sub { is_powerful(0,$_) }, 0..10), "is_powerful(0,n) = 0");
ok( vecall(sub { is_powerful($_,0) }, 1..32), "is_powerful(n,0) = 1 for positive n");
ok( vecall(sub { is_powerful($_,1) }, 1..32), "is_powerful(n,1) = 1 for positive n");

for my $k (3 .. 12) {
  my @nums = (227411960,105218838,79368063,58308379,210322300,44982156,67831696,165946352,243118692,128757041,150085583);
  my @exp = map { fac_is_powerful($_, $k) } 0 .. 32, @nums;
  my @got = map {     is_powerful($_, $k) } 0 .. 32, @nums;
  is_deeply(\@got, \@exp, "is_powerful(n,$k) for 0..32 and 11 larger nums");
}

{
  my @pow2 = map { 5*5 * $_*$_ } 1..50;
  my @npow2 = map { 149 * $_*$_ } 1..50;
  my @pow3 = map { 7*7*7 * $_*$_*$_ } 1..50;
  my @npow3 = map { 4489 * $_*$_*$_ } 1..50;

  ok( vecall(sub{ is_powerful($_,2) }, @pow2), "small is_powerful(n,2), n powerful" );
  ok( vecall(sub{ is_powerful($_,3) }, @pow3), "small is_powerful(n,3), n powerful" );
  ok( vecnone(sub{ is_powerful($_,2) }, @npow2), "small is_powerful(n,2), n not powerful" );
  ok( vecnone(sub{ is_powerful($_,3) }, @npow3), "small is_powerful(n,3), n not powerful" );
}

is( is_powerful("1377276413364943226363244108454842276965894752197358387200000"), 0, "large easy non-powerful number" );
is( is_powerful("2346889178458529643625998598305409091755415961600000"), 1, "large easy powerful number" );

is( is_powerful("56648008573112538662596929676588737208124071038924666321487873929306609840197", 30), 0, "256-bit semiprime is not 30-powerful, without factoring" );


###### powerful_count
is_deeply( [map { powerful_count($_) } -16..0], [map{0}-16..0], "powerful_count(-n)=0" );
is_deeply( [map { powerful_count($_,0) } 0..20], [0..20], "powerful_count(n,0)=n" );
is_deeply( [map { powerful_count($_,1) } 0..20], [0..20], "powerful_count(n,1)=n" );

# test with n <= 0 and k = 0,1,2
is_deeply( [map { powerful_count($_,0) } -100, -10, -4, -1, 0, 1, 4, 10, 100],
           [0, 0, 0, 0, 0, 1, 4, 10, 100], "powerful_count(+/- n, 0)" );
is_deeply( [map { powerful_count($_,0) } -100, -10, -4, -1, 0, 1, 4, 10, 100],
           [0, 0, 0, 0, 0, 1, 4, 10, 100], "powerful_count(+/- n, 1)" );
is_deeply( [map { powerful_count($_,2) } -100, -10, -4, -1, 0, 1, 4, 10, 100],
           [0, 0, 0, 0, 0, 1, 2,  4 , 14], "powerful_count(+/- n, 2)" );

is_deeply( [map { powerful_count($_) } 0..20],
           [0,1,1,1,2,2,2,2,3,4,4,4,4,4,4,4,5,5,5,5,5],
           "powerful_count(0..20)" );
is_deeply( [map { powerful_count($_,3) } 0..20],
           [0,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3],
           "powerful_count(0..20,3)" );
{
  my $maxk = ($extra) ? 30 : 15;
  my @pow14=(0,14,100,432, 2048, 6561, 16384, 59049, 131072, 262144, 531441, 1594323, 4194304, 8388608, 16777216, 43046721, 129140163, 268435456, 536870912, 1162261467, 3486784401, 8589934592, 17179869184, 34359738368, 68719476736, 137438953472, 274877906944, 549755813888, 1099511627776, 2199023255552, 4398046511104);
  is_deeply( [map { powerful_count($pow14[$_],$_) } 1..$maxk],
             [map { 14 } 1..$maxk],
             "powerful_count(x,1..$maxk) = 14" );
  is_deeply( [map { powerful_count($pow14[$_]-1,$_) } 1..$maxk],
             [map { 13 } 1..$maxk],
             "powerful_count(x-1,1..$maxk) = 13" );
}

# Rather slow in PP
if ($extra) {
  my @exp = (4,14,54,185,619,2027,6553,21044,67231,214122,680330,2158391,6840384,21663503);
  my @got = map { powerful_count(10**$_) } 1..14;
  is_deeply(\@got, \@exp, "2-powerful_count 10^1, 10^2, ..., 10^14");
}
if ($extra) {
  my @exp = (1, 1, 4, 10, 16, 26, 46, 77, 129, 204, 318, 495, 761, 1172, 1799, 2740, 4128, 6200, 9224, 13671, 20205, 29764);
  my @got = map { powerful_count("1".("0"x$_),7) } 1..22;
  is_deeply(\@got, \@exp, "7-powerful_count 10^1, 10^2, ..., 10^22");
}

###### nth_powerful

is(nth_powerful(0), undef, "nth_powerful(0) returns undef");
is(nth_powerful(100),3136,"3136 is the 100th powerful number");
SKIP: {
  skip "Skipping nth_powerful for k>2 in PP", 2 unless $usexs;
  is(nth_powerful(100,6),43046721,"43046721 is the 100th 6-powerful number");
  is(nth_powerful(12,15),16777216,"16777216 is the 12th 15-powerful number");
}

###### sumpowerful
is_deeply( [map { sumpowerful($_) } -16..0], [map{0}-16..0], "sumpowerful(-n)=0" );
is_deeply( [map { sumpowerful($_) } 0..48],
           [0,1,1,1,5,5,5,5,13,22,22,22,22,22,22,22,38,38,38,38,38,38,38,38,38,63,63,90,90,90,90,90,122,122,122,122,158,158,158,158,158,158,158,158,158,158,158,158,158],
           "sumpowerful(n) for 0 <= n <= 48" );
is_deeply( [map { sumpowerful($_,3) } 0..48],
           [0,1,1,1,1,1,1,1,9,9,9,9,9,9,9,9,25,25,25,25,25,25,25,25,25,25,25,52,52,52,52,52,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84],
           "sumpowerful(n,3) for 0 <= n <= 48" );
is_deeply( [map { sumpowerful(17411,$_) } 0..16],
           [151580166,151580166,1464625,333416,164098,101342,57807,41389,39074,32257,31745,30721,28673,24577,16385,1,1],
           "sumpowerful(17411,k) for 0 <= k <= 16" );

is_deeply( [map { sumpowerful("2147516495",$_) } 1..33],
           [qw/2305913549222300760 71073461134258 2727672189281 542650082891 192623487712 91172645015 57240053947 37822907405 26438551880 19617088953 13157238054 7502933431 7258257269 6035951629 6031152276 6016770601 5973658344 5844387109 5456704476 4293918721 4292870145 4290772993 4286578689 4278190081 4261412865 4227858433 4160749569 4026531841 3758096385 3221225473 2147483649 1 1/],
          "sumpowerful(2147516495,k) for 1 <= k <= 33" );

is_deeply( [map { sumpowerful("1234567890123456",$_) } 3..32],
           [qw/146043398655792412070 10996713169431264132 2301916287502408997 793628382930863389 358072984365696929 184880296279236241 108043541271736385 72290570743670714 53667950830214223 40306694671659861 28363055737738077 22793155018675414 18357777378614938 13885855208881771 11386819588192744 8343600537655001 5597666622659511 3750512753514143 3655141834040541 3178294213387061 3178262828133148 3178168676565713 3177886230252016 3177038908088141 3174496975150948 3166871243448233 3143994182557816 3075363268322021 2869471062485548 2251795518717953/],
           "sumpowerful(1234567890123456,k) for 3 <= k <= 32" );
is( sumpowerful("1234567890123456",1), "762078937661941480719405753696", "sumpowerful(1234567890123456,1) = (n*(n+1))/2" );
#SKIP: {
#  skip "Skipping sumpowerful(1234567890123456,2)",1 unless $extra;
#  is( sumpowerful("1234567890123456",2), "31374760178828970927228", "sumpowerful(1234567890123456,2)" );
#}

###### powerful_numbers
is_deeply( powerful_numbers(40,180,3), [64,81,125,128], "powerful_numbers(40,180,3)");
is_deeply( powerful_numbers(9,20,0), [9..20], "powerful_numbers(9,20,0) = 9..20");
is_deeply( powerful_numbers(9,20,1), [9..20], "powerful_numbers(9,20,1) = 9..20");
is_deeply( powerful_numbers(120), [1,4,8,9,16,25,27,32,36,49,64,72,81,100,108], "powerful_numbers(120)");
is_deeply( powerful_numbers(9,120), [9,16,25,27,32,36,49,64,72,81,100,108], "powerful_numbers(9,120)");
is_deeply( powerful_numbers(9,200,2), [9,16,25,27,32,36,49,64,72,81,100,108,121,125,128,144,169,196,200], "powerful_numbers(9,200,2)");
is_deeply( powerful_numbers(0,200,3), [1,8,16,27,32,64,81,125,128], "powerful_numbers(0,200,3)");
is_deeply( powerful_numbers(1,200,4), [1,16,32,64,81,128], "powerful_numbers(1,200,4)");
is_deeply( powerful_numbers(1,1000,5), [1,32,64,128,243,256,512,729], "powerful_numbers(1,1000,5)");
is_deeply( powerful_numbers(1000000000000, 1010000000000,5),
           [1000000000000,1004193907488,1007769600000,1008394404608],
           "powerful_numbers(1e12,1e12+1e10,5)");


sub fac_is_powerful {
  my($n, $k) = @_;
  $k = 2 unless defined $k;
  return 0 if $n < 1;
  return 1 if $n == 1 || $k <= 1;
  return 0 if $n < (1<<$k);
  return 0 if (!($n%2)) && ($n%4);
  return (vecall { $_->[1] >= $k } factor_exp($n)) ? 1 : 0;
}
