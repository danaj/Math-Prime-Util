#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/irand irand64 drand urandomb urandomm
                         random_bytes entropy_bytes
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
            + 4  # identify rng and test srand/csrand
            + 4  # 0 / undef arguments to urandom*
            + 1  # urandomb
            + 3  # urandomm
            + 4  # entropy_bytes
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

my $core_rand = "not drand48";
if ($use64) {
  my @r = map { CORE::rand() } 0..8;
  if (try_lcg(25214903917,11,2**48,@r)) {
    $core_rand = "drand48 (yech)";
  } elsif (try_16bit(@r)) {
    $core_rand = "16-bit (ack)";
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

# Quick check to identify the RNG being used.  Should be ChaCha20.
srand(42);
my $rb42 = irand();
my $csprng = 'something I do not know';
if    ($rb42 ==  445265827) { $csprng = 'ChaCha20'; }
elsif ($rb42 == 3626765506) { $csprng = 'ChaCha12'; }
elsif ($rb42 ==  266717191) { $csprng = 'ChaCha8'; }
elsif ($rb42 == 4274346485) { $csprng = 'ISAAC'; }
elsif ($rb42 == 3197710526) { $csprng = 'drand48'; }
elsif ($rb42 == 2209484588) { $csprng = 'Math::Random::Xorshift'; }
elsif ($rb42 == 1608637542) { $csprng = 'Math::Random::MT'; }
elsif ($rb42 == 2746317213) { $csprng = 'Math::Random::MT::Auto (32)'; }
elsif ($rb42 == 6909045637428952499) { $csprng = 'Math::Random::MTwist (64)'; }
elsif (sprintf("%.1lf",$rb42) eq '6909045637428952064.0') { $csprng = 'Math::Random::MTwist (32)'; }
elsif ($rb42 == 9507361240820437267) { $csprng = 'Math::Random::MT::Auto (64)'; }
diag "CORE::rand: $core_rand. Our PRNG: $csprng";

SKIP: {
  if ($csprng eq 'ChaCha20') {
    srand(15);
    is(unpack("H8",random_bytes(4)), "546d6108", "random_bytes after srand");
    csrand("BLAKEGrostlJHKeccakSkein--RijndaelSerpentTwofishRC6MARS");
    is(unpack("H14",random_bytes(7)), "b302e671601bce", "random_bytes after manual seed");
    is(irand(), 88564645, "irand after seed");
    my $d = drand();  my $dexp = 0.0459118340827543;
    ok($d > $dexp-1e-6 && $d < $dexp+1e-6,"drand after seed $d ~ $dexp");
  } elsif ($csprng eq 'ISAAC') {
    srand(15);
    is(unpack("H8",random_bytes(4)), "36cd2d21", "random_bytes after srand");
    csrand("BLAKEGrostlJHKeccakSkein--RijndaelSerpentTwofishRC6MARS");
    is(unpack("H14",random_bytes(7)), "a0644ad1e00324", "random_bytes after manual seed");
    is(irand(), 2526495644, "irand after seed");
    my $d = drand();  my $dexp = 0.490707771279301221;
    ok($d > $dexp-1e-6 && $d < $dexp+1e-6,"drand after seed $d ~ $dexp");
  } else {
    skip "Unknown random number generator!  Skipping deterministic tests.",4;
  }
}

srand;

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
  for my $m (1..50) {
    my $r = urandomm($m);
    push @failm, $m unless !ref($r) && $r < $m;
  }
  is_deeply(\@failm, [], "urandomm returns native int within range for 1..50");
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

#######

# If the functions work, these tests fail with chance less than 2^-128.
my $ebytes = 17;
my $eb1 = entropy_bytes($ebytes);
my $eb2 = entropy_bytes($ebytes);
is(length($eb1), $ebytes, "entropy_bytes gave us the right number of bytes");
$eb1 = unpack("H*",$eb1);
$eb2 = unpack("H*",$eb2);
isnt($eb1, '00' x $ebytes, "entropy_bytes didn't return all zeros once");
isnt($eb2, '00' x $ebytes, "entropy_bytes didn't return all zeros twice");
isnt($eb1, $eb2, "entropy_bytes returned two different binary strings");
