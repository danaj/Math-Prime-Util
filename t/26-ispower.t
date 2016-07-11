#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_power is_prime_power vecsum/;
#use Math::BigInt try=>"GMP,Pari";

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @pow1  = (0,0,0,0,2,0,0,0,3,2,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,2,0,3,0,0,0,0,5);
my @ppow1 = (0,0,1,1,2,1,0,1,3,2,0,1,0,1,0,0,4,1,0,1,0,0,0,1,0,2,0,3,0,1,0,1,5);

my %bpow = (
  "16926659444736" => [6,17],
  "100000000000000000" => [10,17],
  "609359740010496" => [6,19],
  "789730223053602816" => [6,23],
  "4611686018427387904" => [2,62],
  "4738381338321616896" => [6,24],
  "9223372036854775808" => [2,63],
  "12157665459056928801" => [3,40],
  "10000000000000000000" => [10,19],
);

my %bppow = (
  "762939453125" => [5,17],
  "232630513987207" => [7,17],
  "11398895185373143" => [7,19],
  "11920928955078125" => [5,23],
  "68630377364883" => [3,29],
  "617673396283947" => [3,31],
  "450283905890997363" => [3,37],
  "12157665459056928801" => [3,40],
  "7450580596923828125" => [5,27],
  "3909821048582988049" => [7,22],
  "5559917313492231481" => [11,18],
  "8650415919381337933" => [13,17],
  "2862423051509815793" => [17,15],
  "15181127029874798299" => [19,15],
);

if ($use64) {
 #push @invmods, [ 13, 9223372036854775808, 5675921253449092805 ];
 #push @invmods, [ 14, 18446744073709551615, 17129119497016012214 ];
}

plan tests => 0
            + 2
            + 2 + 2*$extra
            + scalar(keys(%bpow))
            + scalar(keys(%bppow))
            + 4;
            + 0;

is_deeply( [map { is_power($_) } 0 .. $#pow1],        \@pow1,  "is_power 0 .. $#pow1" );
is_deeply( [map { is_prime_power($_) } 0 .. $#ppow1], \@ppow1, "is_prime_power 0 .. $#ppow1" );

is( vecsum(map { is_power(2099*$_+$_+1) } 0..200), 8, "is_power 200 small ints" );
is( vecsum(map { is_prime_power(2099*$_+$_+1) } 0..200), 77, "is_prime_power 200 small ints" );

if ($extra) {
  is( vecsum(map { is_power(23*$_+$_) } 0..10000), 122, "is_power 10k small ints" );
  is( vecsum(map { is_prime_power(23*$_+$_+1) } 0..10000), 2829, "is_prime_power 10k small ints" );
}

while (my($n, $expect) = each (%bpow)) {
  my $r;  my $k = is_power($n,0,\$r);
  is_deeply( $expect, [$r,$k], "ispower  =>  $n = $r^$k (@$expect)" );
}

while (my($n, $expect) = each (%bppow)) {
  my $r;  my $k = is_prime_power($n,\$r);
  is_deeply( $expect, [$r,$k], "isprimepower  =>  $n = $r^$k (@$expect)" );
}

{
  my $r;
  my $ip = is_power(-8,3,\$r);
  is( $ip, 1, "-8 is a third power" );
  is( $r, -2, "-8 is a third power of -2" );

  is( is_power(-8, 4), 0, "-8 is not a fourth power" );
  is( is_power(-16,4), 0, "-16 is not a fourth power" );
}
