#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_divisible is_congruent modint irand irand64/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

# From GMP test suite
my @divt = (
  [0, 0, 1],
  [17, 0, 0],
  [0, 1, 1],
  [123, 1, 1],
  [-123, 1, 1],
  [0, 2, 1],
  [1, 2, 0],
  [2, 2, 1],
  [-2, 2, 1],
  ["340282366920938463463374607431768211456", 2, 1],
  ["340282366920938463463374607431768211457", 2, 0],
  ["3689348814741910323",3,1],
  ["3689348814741910322",3,0],
  ["68056473384187692692674921486353642291", 3, 1],
  ["68056473384187692692674921486353642290", 3, 0],
  ["3689348813882916864", "6442450944", 1],
  ["68056473384187692688985572671611731968", "27670116110564327424", 1],
  ["408338840305126156152360180103379943424", "27670116110564327424", 0],
  # Our additional single divisor tests
  # Check UV n with negative d larger than IV
  ["10223372036854775807","-10223372036854775807", 1],
  ["36472996418050588672","33171997",1],
);
my @congt = (
  [0,0,0,1],     # strict equality mod 0
  [11,11,0,1],
  [3,11,0,0],
  [0,0,1,1],     # anything is congruent mod 1
  [1,0,1,1],
  [0,1,1,1],
  [123,456,1,1],
  ["335812727629498640265","2812431594283598168865",1,1],
  ["3689348814741910323","858993459","6442450944",1],
  ["68056473384187692692674921486353642291","3689348814741910323","27670116110564327424",1],
  ["18325193793", "-9162596895", "13743895344", 1],
  ["78706108047827420225", "-39353054023913710111", "59029581035870565168", 1],
);

plan tests => 1 + 9*5 + scalar(@divt) + 2   # is_divisible
            + 5 + 2 + scalar(@congt);       # is_congruent

my $ntests = $extra && $use64 ? 10000 : $extra ? 1000 : 100;
my @r32 = map { irand() } 0..$ntests;
my @r64 = map { irand64() } 0..$ntests;

###### is_divisible

is_deeply( [map { is_divisible($_,0) } -20..20],
           [map { 0+($_ == 0)        } -20..20],
           "is_divisible(x,0) = 0 for x != 0");

for my $d (1 .. 9) {
  is_deeply( [map { is_divisible($_,$d)    } @r32],
             [map { 0+(modint($_,$d) == 0) } @r32],
             "is_divisible(x,$d) for 32-bit x" );
  is_deeply( [map { is_divisible($_,$d)    } @r64],
             [map { 0+(modint($_,$d) == 0) } @r64],
             "is_divisible(x,$d) for 64-bit x" );
  is_deeply( [map { is_divisible(-$_,$d)    } @r32],
             [map { 0+(modint(-$_,$d) == 0) } @r32],
             "is_divisible(-x,$d) for 32-bit x" );
  is_deeply( [map { is_divisible($_,-$d)    } @r32],
             [map { 0+(modint($_,-$d) == 0) } @r32],
             "is_divisible(x,-$d) for 32-bit x" );
  is_deeply( [map { is_divisible(-$_,-$d)    } @r32],
             [map { 0+(modint(-$_,-$d) == 0) } @r32],
             "is_divisible(-x,-$d) for 32-bit x" );
}

foreach my $r (@divt) {
  my($n, $d, $exp) = @$r;
  is(is_divisible($n,$d), $exp, "is_divisible($n,$d) = $exp");
}

is(is_divisible("26000117000117",2,3,5,7,11),0,"is_divisible(26000117000117,2,3,5,7,11)");
is(is_divisible("26000117000117",2,3,5,7,11,13),1,"is_divisible(26000117000117,2,3,5,7,11,13)");

###### is_congruent

for my $c (-2 .. 2) {
  is_deeply( [map { is_congruent($_,$c,0) } -3..3],
             [map { 0+($_ == $c)          } -3..3],
             "is_congruent(x,$c,0) = 0 for x != $c" );
}

is_deeply( [map { is_congruent($_,3,13) } @r32, @r64],
           [map { 0+(modint($_,13)==3)  } @r32, @r64],
           "is_congruent(x,3,13) for 32-bit and 64-bit x" );
is_deeply( [map { is_congruent($_,-27,17) } @r32, @r64],
           [map { 0+(modint($_,17)==7)    } @r32, @r64],
           "is_congruent(x,-27,17) for 32-bit and 64-bit x" );

foreach my $r (@congt) {
  my($n, $c, $d, $exp) = @$r;
  is(is_congruent($n,$c,$d), $exp, "is_congruent($n,$c,$d) = $exp");
}
