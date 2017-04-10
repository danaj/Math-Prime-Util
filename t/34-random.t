#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/irand irand64 drand random_bytes urandomb urandomm
                         srand csrand
                         mulmod addmod vecmin vecmax vecall/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

my $samples = $extra ? 100000 :  10000;

plan tests => 1
            + 2
            + 2
            + 2
            + 5  # drand range
            + 0  # TODO 2 tests for srand / csrand here
            + 4
            + 1
            + 3
            + 0;

########

ok( Math::Prime::Util::_is_csprng_well_seeded(), "CSPRNG is being seeded properly" );

########

{
  my @s = map { irand } 1 .. $samples;
  is( scalar(grep { $_ > 4294967295 } @s), 0, "irand values are 32-bit" );
  is( scalar(grep { $_ != int($_) } @s), 0, "irand values are integers" );
}

########

SKIP: {
  skip "Skipping irand64 on 32-bit Perl", 2 if !$use64;
  my $bits_on  = 0;
  my $bits_off = 0;
  my $iter = 0;
  for (1 .. 6400) {
    $iter++;
    my $v = irand64;
    $bits_on |= $v;
    $bits_off |= (~$v);
    last if ~$bits_on == 0 && ~$bits_off == 0;
  }
  is( ~$bits_on,  0, "irand64 all bits on in $iter iterations" );
  is( ~$bits_off, 0, "irand64 all bits off in $iter iterations" );
}

########

# This is really brute force, but it doesn't take too long.
{
  my $mask = 0;
  my $v;
  for (1..1024) {
    $v = drand;
    last if $v >= 1;
    next if $v < .5;
    for my $b (0..127) {
      last unless $v;
      $v *= 2;
      if ($v >= 1) {
        $mask |= (1 << $b);
        $v -= 1;
      }
    }
  }
  ok($v < 1, "drand values between 0 and 1-eps");
  my $k = 0; while ($mask) { $k++; $mask >>= 1; }
  # Assuming drand is working properly:
  #   k = 24   NV is float
  #   k = 53   NV is double
  #   k = 64   NV is long double
  # If we used drand48 we'd get 48 with double or long double.
  ok($k >= 21, "drand supplies at least 21 bits (got $k)");
}

sub check_float_range {
  my($name, $lo, $hi, $v) = @_;
  if ($lo <= $hi) {
    ok( vecall(sub{ $_ >= $lo && $_ < $hi },@$v), "$name: all in range [$lo,$hi)" );
  } else {
    ok( vecall(sub{ $_ >= $hi && $_ < $lo },@$v), "$name: all in range ($hi,$lo]" );
  }
}
my $num = $extra ? 1000 : 100;
check_float_range('drand(10)',0,10,[map{ drand(10) } 1..$num]);
check_float_range('drand()',0,1,[map{ drand() } 1..$num]);
check_float_range('drand(-10)',0,-10,[map{ drand(-10) } 1..$num]);
check_float_range('drand(0)',0,1,[map{ drand(0) } 1..$num]);
{
  # Skip warnings these give, worry about the behavior
  no warnings;
  check_float_range('drand(undef)',0,1,[map{ drand(undef) } 1..$num]);
}
# We can't easily supress the warning here, but we'd like to check the
# result.  Math::Random::Secure fails this, for instance.
#check_float_range('drand("foo")',0,1,[map{ drand("foo") } 1..$num]);

########

if ($use64) {
  my @r = map { CORE::rand() } 0..8;
  if (try_lcg(25214903917,11,2**48,@r)) {
    diag "CORE::rand is cruddy drand48 as expected";
  } elsif (try_16bit(@r)) {
    diag "CORE::rand looks like 16-bit.  Ugh.";
  } else {
    diag "CORE::rand is not drand48";
  }
}

sub try_lcg {
  my($a,$c,$m,@r) = @_;
  @r = map { int($m * $_) } @r;
  my @g = ($r[0]);
  $g[$_] = addmod(mulmod($a,$g[$_-1],$m),$c,$m) for 1..$#r;
  for (1..$#r) {
    return unless $r[$_] == $g[$_];
  }
  1;
}
# We could try to predict Windows truncated LCG:
#   http://crypto.stackexchange.com/questions/10608/how-to-attack-a-fixed-lcg-with-partial-output

sub try_16bit {
  my(@r) = @_;

  for my $r (@r) {
    my $rem = $r - int(32768*$r);
    return if $rem > 1e-6;
  }
  for my $r (map { CORE::rand() } 1..120) {
    my $rem = $r - int(32768*$r);
    return if $rem > 1e-6;
  }
  1;
}

########

# TODO: deterministic rand
srand(15);
# is(unpack("H8",random_bytes(4)), "8a488975", "random_bytes after srand");
csrand("BLAKEGrostlJHKeccakSkein--RijndaelSerpentTwofishRC6MARS");
# is(unpack("H14",random_bytes(7)), "35e8f156bcab4c", "random_bytes after manual seed");

#######

is(random_bytes(0),'',"random_bytes(0) returns empty string");
is(urandomb(0),0,"urandomb(0) returns 0");
is(urandomm(0),0,"urandomm(0) returns 0");
is(urandomm(1),0,"urandomm(1) returns 0");

#######

{
  my @failb;
  for my $bits (1..$maxbits) {
    my $lim = (1<<($bits-1)) + ((1<<($bits-1))-1);
    my $r = urandomb($bits);
    push @failb, $bits unless !ref($r) && $r <= $lim;
  }
  is_deeply(\@failb, [], "urandomb returns native int within range for 1..$maxbits");
}

#######

{
  my @failm;
  for my $m (10..50) {
    my $r = urandomm($m);
    push @failm, $m unless !ref($r) && $r < $m;
  }
  is_deeply(\@failm, [], "urandomm returns native int within range for 10..50");
}

{
  my %dv;
  for my $t (1..10000) {
    $dv{urandomm(10)}++;
    last if $t > 100 && scalar(keys(%dv)) >= 10;
  }
  my @k = sort { $a<=>$b} keys(%dv);
  is(scalar(@k), 10, "urandomm(10) generated 10 distinct values");
  ok( vecmin(@k) == 0 && vecmax(@k) == 9, "urandomm(10) values between 0 and 9 (@k)" );
}
