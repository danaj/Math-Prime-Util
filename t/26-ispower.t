#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_power is_prime_power is_square vecsum/;
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

my %powers = (
  0 => [-2, -1, 0, 1, 2, 3, 5, 6, 7, 10, 11, 12, 13, 14, 15, 17, 18, 19],
  2 => [4, 9, 25, 36, 49],
  3 => [8, 27, 125, 343, 17576],
  4 => [16, 38416],
  9 => [19683, 1000000000],
);
if ($use64) {
  push @{$powers{0}}, 9908918038843197151;
  push @{$powers{2}}, 18446743927680663841;
  push @{$powers{3}}, 2250923753991375;
  push @{$powers{4}}, 1150530828529256001;
  push @{$powers{9}}, 118587876497;
}
my @negpowers = (0,0,0,3,0,5,3,7,0,9,5);

plan tests => 0
            + 2
            + 2 + 2*$extra
            + scalar(keys(%bpow))
            + scalar(keys(%bppow))
            + 4
            + 7 + scalar(keys %powers) + scalar(@negpowers)
            + 3  # is_square
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

###### is_power
while (my($e, $vals) = each (%powers)) {
  my @fail;
  foreach my $val (@$vals) {
    push @fail, $val unless is_power($val) == $e;
  }
  ok( @fail == 0, "is_power returns $e for " . join(",",@fail) );
}
foreach my $e (0 .. $#negpowers) {
  is( is_power(-7 ** $e), $negpowers[$e], "is_power(-7^$e ) = $negpowers[$e]" );
}
is( is_power(-1,5), 1, "-1 is a 5th power" );
{
  my($ispow, $root);
  $ispow = is_power(24, 2, \$root);
  is( $ispow, 0, "24 isn't a perfect square...");
  is( $root, undef, "...and the root wasn't set");
  $ispow = is_power( "1000093002883029791", 3, \$root);
  is( $ispow, 1, "1000031^3 is a perfect cube...");
  is( $root, 1000031, "...and the root was set");
  $ispow = is_power( 36**5 , 0, \$root);
  is( $ispow, 10, "36^5 is a 10th power...");
  is( $root, 6, "...and the root is 6");
}

###### is_square
is_deeply(
  [map { is_square($_) } (-4 .. 16)],
  [0,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1],
  "is_square for -4 .. 16"
);
is(is_square(603729), 1, "603729 is a square");
is(is_square("765413284212226299051111674934086564882382225721"), 1, "is_square(<square of 80-bit prime>) = 1");
