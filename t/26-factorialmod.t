#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/factorialmod factorial/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

plan tests => 11;

{
  my @result = map { my $m=$_; map { factorialmod($_,$m) } 0..$m-1; } 1 .. 40;
  my @expect = map { my $m=$_; map { "".(factorial($_) % $m); } 0..$m-1; } 1 .. 40;
  is_deeply( \@result, \@expect, "factorialmod n! mod m for m 1 to 50, n 0 to m" );
}

is( factorialmod(5,0), undef, "factorialmod returns undef for modulus 0" );

SKIP: {
  skip "medium size factorialmods in PP",2 unless $usexs || $extra;
  is( factorialmod(1000000000,1000000008), 0, "1000000000! mod 1000000008 is zero" );
  is( factorialmod(50000,10000019), 8482159, "50000! mod 10000019" );
}

SKIP: {
  skip "large value without EXTENDED_TESTING",1 unless $extra;
  # m = 8000036000054000027 = 2000003^3
  is( "".factorialmod(7974,"8000036000054000027"), "2352233842797337", "factorialmod with large n and large composite non-square-free m" );
}

is( "".factorialmod(1000,"18446744082299486209"), "9289597790693496878",
    "factorialmod with bigint modulus" );

{
  my $m = 1000003;
  my @n = ($m-60 .. $m-1);
  my @got = map { factorialmod($_, $m) } @n;
  my @exp;
  $exp[0] = $got[0];
  for my $i (1 .. $#n) {
    $exp[$i] = ($exp[$i-1] * $n[$i]) % $m;
  }
  is_deeply(\@got, \@exp, "factorialmod recurrence near prime modulus (backward path)");
}

{
  my @fmres = (
    [1000,  898920,  730681],
    [1011, 8011020,   36198],
    [1022, 1126740, 6649423],
    [1023, 2652720, 2346809],
    [1024, 6379860, 3127856],
    [1040, 6866220, 8739398],
    [1055, 7700400, 3154085],
    [1056, 1606140,  707433],
    [1072, 8160780, 3630762],
    [1087, 8516340, 8887660],
    [1088, 5759400, 9755726],
    [1089, 1974060, 3965436],
    [1104, 9864900, 7692904],
    [1120,  962520, 8078184],
  );
  my(@got1,@exp1,@got2,@exp2);
  for my $V (@fmres) {
    my($n,$exp1,$exp2) = @$V;
    push @got1, factorialmod($n,10000020);  push @exp1, $exp1;
    push @got2, factorialmod($n,10000019);  push @exp2, $exp2;
  }
  is_deeply(\@got1, \@exp1, "factorialmod in segmented non-mont path");
  is_deeply(\@got2, \@exp2, "factorialmod in segmented mont path");
}

{
  my $m = 10403;  # 101 * 103
  is(factorialmod(102, $m), 5252,
     "factorialmod below composite-zero cutoff remains non-trivial");
  is(factorialmod(103, $m), 0,
     "factorialmod at composite-zero cutoff returns zero");
}
