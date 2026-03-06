#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/lshiftint rshiftint rashiftint/;

my $use64 = (~0 > 4294967295);
my $bits  = $use64 ? 64 : 32;

# Negative right shifts:
#   ">>"    MPU, Pari/GP, Mathematica   = -rshiftint(-n,k]
#   ">>a"   Math::BigInt, Python, Java
my @negshifts = (
  # n, k,  <<,   >>,  >>arith
  [ 0, 1,  0,    0,   0],
  [-1, 1,  -2,   0,  -1],
  [-5, 1,  -10, -2,  -3],
  [-8, 2,  -32, -2,  -2],

  [qw/-65535 15 -2147450880 -1 -2/],
  [qw/-65536 15 -2147483648 -2 -2/],
  [qw/-65535 16 -4294901760 0 -1/],
  [qw/-65536 16 -4294967296 -1 -1/],

  [qw/-65535 47 -9223231299366420480 0 -1/],   #  8
  [qw/-65536 47 -9223372036854775808 0 -1/],   #  9
  [qw/-65535 48 -18446462598732840960 0 -1/],  # 10
  [qw/-65536 48 -18446744073709551616 0 -1/],  # 11
  [qw/-65536 80 -79228162514264337593543950336 0 -1/],  # 12

  [qw/-307385513 6 -19672672832 -4802898 -4802899/],
  [qw/-637526413 6 -40801690432 -9961350 -9961351/],
  [qw/-2045651239 6 -130921679296 -31963300 -31963301/],
  [qw/-3675663743 6 -235242479552 -57432245 -57432246/],
  [qw/-2332267979728172537 6 -149265150702603042368 -36441687183252695 -36441687183252696/],
  [qw/-8408654401686460807 6 -538153881707933491648 -131385225026350950 -131385225026350951/],
  [qw/-17640827963513397449 6 -1129012989664857436736 -275637936929896835 -275637936929896836/],
  [qw/-32659506018295865747 6 -2090208385170935407808 -510304781535872902 -510304781535872903/],
  [qw/-79231600218559026832557301750107210001 6 -5070822413987777717283667312006861440064 -1237993753414984794258707839845425156 -1237993753414984794258707839845425157/],
  [qw/-131954888069700539887213633881194728277 6 -8445112836460834552781672568396462609728 -2061795126089070935737713029393667629 -2061795126089070935737713029393667630/],
  [qw/-254262665582332530470619504253273698569 6 -16272810597269281950119648272209516708416 -3972854149723945788603429753957401540 -3972854149723945788603429753957401541/],
  [qw/-416649423645764932216789232242651032187 6 -26665563113328955661874510863529666059968 -6510147244465077065887331753791422377 -6510147244465077065887331753791422378/],
);

plan tests => 4 + 3 + 2   # original tests
            + 3            # shift by 0
            + 1            # rshiftint == rashiftint for non-negative
            + 1            # big right shift positive n
            + 1            # big right shift negative n (rshiftint)
            + 1            # big right shift negative n (rashiftint)
            + 1            # big right shift exactly BITS_PER_WORD
            + 1            # negative k flips direction
            + 1            # lshift at UV boundary produces bigint
            + 1            # round-trip: rshift(lshift(n,k),k) == n
            + 1            # rshiftint at k = bits-1
            + 1            # rashiftint at k = bits-1 for negative
            + 1            # lshiftint k = bits-1 for small n
            ;

###### Basic: small ranges with implied k=1
is_deeply([map { lshiftint($_) } 0..50], [map { $_ << 1 } 0..50], "lshiftint(0..50)");
is_deeply([map { rshiftint($_) } 0..50], [map { $_ >> 1 } 0..50], "rshiftint(0..50)");
is_deeply([map { rashiftint($_) } 0..50], [map { $_ >> 1 } 0..50], "rashiftint(0..50)");
is_deeply([map { lshiftint($_,5) } -65 .. 65], [map { $_ * 32 } -65 .. 65], "lshiftint(-65 .. 65, 5)");

# lshiftint for native size k is:  mulint($n, 1 << $k)
# but for testing we want to avoid using our other functions.

###### Negative n table
is_deeply( [map { "".lshiftint($_->[0], $_->[1]) } @negshifts],
           [map { $_->[2] } @negshifts],
           "left shift negative inputs" );
is_deeply( [map { "".rshiftint($_->[0], $_->[1]) } @negshifts],
           [map { $_->[3] } @negshifts],
           "right shift negative inputs" );
is_deeply( [map { "".rashiftint($_->[0], $_->[1]) } @negshifts],
           [map { $_->[4] } @negshifts],
           "signed arithmetic right shift negative inputs" );

###### Boundary left shifts
is("".lshiftint("2147483648"),"4294967296","left shift of 2^31 with implied 1 bit");
is("".lshiftint("9223372036854775808"),"18446744073709551616","left shift of 2^63 with implied 1 bit");


###### Shift by 0 is identity
{
  my @vals = (0, 1, -1, 127, -128, "4294967295", "-4294967296",
              "9223372036854775807", "-9223372036854775808",
              "18446744073709551615", "18446744073709551616",
              "340282366920938463463374607431768211456");  # 2^128
  is_deeply([map{"$_"}map { lshiftint($_,0) } @vals], [map{"$_"} @vals], "lshiftint(n,0) == n");
  is_deeply([map{"$_"}map { rshiftint($_,0) } @vals], [map{"$_"} @vals], "rshiftint(n,0) == n");
  is_deeply([map{"$_"}map {rashiftint($_,0) } @vals], [map{"$_"} @vals], "rashiftint(n,0) == n");
}

###### rshiftint == rashiftint for non-negative n
{
  my @cases;
  for my $n (0, 1, 2, 7, 255, 65535, "4294967295", "9223372036854775807",
             "18446744073709551615", "340282366920938463463374607431768211456") {
    for my $k (0, 1, 3, 7, 15, 31, 63, 65, 128) {
      push @cases, [$n, $k];
    }
  }
  is_deeply( [map { "".rshiftint($_->[0],$_->[1]) } @cases],
             [map { "".rashiftint($_->[0],$_->[1]) } @cases],
             "rshiftint == rashiftint for non-negative n" );
}

###### Big right shift (k >= BITS_PER_WORD): positive n gives 0
{
  my @nvals = (1, 7, "4294967295");
  push @nvals, ("9223372036854775807", "18446744073709551615") if $use64;
  my @kvals = ($bits, $bits+1, $bits+10, 128, 256);
  my @got;
  for my $n (@nvals) { for my $k (@kvals) {
    push @got, rshiftint($n,$k), rashiftint($n,$k);
  }}
  is_deeply(\@got, [(0) x scalar(@got)], "big right shift of positive n gives 0");
}

###### Big right shift: negative n, rshiftint gives 0 when |n| < 2^k
{
  my @nvals = (-1, -7, "-4294967296");
  push @nvals, "-9223372036854775808" if $use64;
  my @kvals = ($bits+1, 128);
  my @got;
  for my $n (@nvals) { for my $k (@kvals) { push @got, "".rshiftint($n,$k); } }
  is_deeply(\@got, [(0) x scalar(@got)], "rshiftint(negative, big k) == 0");
}

###### Big right shift: negative n, rashiftint gives -1
{
  my @nvals = (-1, -7, "-4294967296");
  push @nvals, ("-9223372036854775808", "-18446744073709551616") if $use64;
  my @kvals = ($bits, $bits+1, 128);
  my @got;
  for my $n (@nvals) { for my $k (@kvals) { push @got, "".rashiftint($n,$k); } }
  is_deeply(\@got, [(-1) x scalar(@got)], "rashiftint(negative, big k) == -1");
}

###### Exactly BITS_PER_WORD shift
{
  my $maxuv = $use64 ? "18446744073709551615" : "4294967295";
  is(rshiftint($maxuv, $bits), 0, "rshiftint(UV_MAX, BITS_PER_WORD) == 0");
}

###### Negative k flips direction
{
  my @cases;
  for my $n (0, 5, -5, 255, -255, "4294967295", "-4294967295",
             "18446744073709551615") {
    for my $k (1, 3, 7, 15) {
      push @cases, [$n, $k];
    }
  }
  my $ok = 1;
  for my $c (@cases) {
    my($n,$k) = @$c;
    $ok = 0 if "".lshiftint($n, -$k) ne "".rshiftint($n, $k);
    $ok = 0 if "".rshiftint($n, -$k) ne "".lshiftint($n, $k);
    $ok = 0 if "".rashiftint($n, -$k) ne "".lshiftint($n, $k);
  }
  ok($ok, "negative k flips direction for all three shift functions");
}

###### Left shift at UV boundary produces bigint
{
  my $maxuv = $use64 ? "18446744073709551615" : "4294967295";
  my $expect = $use64 ? "36893488147419103230" : "8589934590";
  is("".lshiftint($maxuv, 1), $expect, "lshiftint(UV_MAX, 1) produces correct bigint");
}

###### Round-trip: rshiftint(lshiftint(n,k),k) == n for non-negative n
{
  my @cases;
  for my $n (0, 1, 7, 255, 65535, "4294967295", "18446744073709551615",
             "340282366920938463463374607431768211456") {
    for my $k (1, 5, 16, 32, 64) {
      push @cases, [$n, $k];
    }
  }
  my $ok = 1;
  for my $c (@cases) {
    my($n,$k) = @$c;
    $ok = 0 if "".rshiftint(lshiftint($n,$k), $k) ne "$n";
  }
  ok($ok, "rshiftint(lshiftint(n,k),k) == n for non-negative n");
}

###### rshiftint at k = BITS_PER_WORD - 1
{
  my $k = $bits - 1;
  # Right shift UV_MAX by BITS-1 gives 1.
  my $maxuv = $use64 ? "18446744073709551615" : "4294967295";
  is(rshiftint($maxuv, $k), 1, "rshiftint(UV_MAX, BITS-1) == 1");
}

###### rashiftint at k = BITS_PER_WORD - 1 for negative n
{
  # rashiftint(-1, k) == -1 for all k (floor(-1/2^k) = -1)
  my $k = $bits - 1;
  is("".rashiftint(-1, $k), "-1", "rashiftint(-1, BITS-1) == -1");
}

###### lshiftint k = BITS_PER_WORD - 1 for small positive n
{
  my $k = $bits - 1;
  # lshiftint(1, BITS-1) = 2^(BITS-1)
  my $expect = $use64 ? "9223372036854775808" : "2147483648";
  is("".lshiftint(1, $k), $expect, "lshiftint(1, BITS-1) == 2^(BITS-1)");
}
