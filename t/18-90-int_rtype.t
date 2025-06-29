#!/usr/bin/env perl
use strict;
use warnings;

# TODO:  Verify more functions than just addint

use Test::More;
use Math::Prime::Util qw/addint subint add1int sub1int mulint/;
use Math::BigInt;

# Integer forms:
#   1) native (IV or UV).  Good for all 32-bit and all 64-bit post 5.6.
#   2) bigint.  Good but horribly slow in most cases.
#   3) string.  We don't have strong typing so it's dodgy how these get
#      converted if using standard Perl operations.  This is one thing
#      out module is trying to do correctly.
#
# Regardless of input form, we want a correct result.  Efficiently if possible.
#
# The output form expected is:
#   1) native if it fits, regardless of input form.
#   2) bigint otherwise, regardless of input form.
#
# We will attempt to preserve the input type, typically based on the first
# argument.  That is, if the input was a "Math::GMPz" object then we will
# try to output either native or Math::GMPz.

# Expect all the 5.6 workarounds to be deprecated in 2026.

my $isoldperl   = $] < 5.008;
my $is32bit     = ~0 == 4294967295;
my $is64bit     = ~0  > 4294967295;
my $isbroken64  = 18446744073709550592 == ~0;

my $intmax = (!$isoldperl || $is32bit) ? ~0 : 562949953421312;
my $intmin = $is32bit ? -2147483648 : !$isoldperl ? -9223372036854775808 : -562949953421312;

diag "Perl safe min: $intmin";
diag "Perl safe max: $intmax";

my($mm1, $m, $mp1, $mp2);  #  ~0-1,  ~0,  ~0+1,  ~0+2

if ($is32bit) {
  ($mm1,$m,$mp1,$mp2) = (4294967294, 4294967295, 4294967296, 4294967297);
} elsif ($isbroken64) {
  ($mm1,$m,$mp1,$mp2) = (562949953421311,562949953421312,"562949953421313","562949953421314");
} else {
  ($mm1,$m,$mp1,$mp2) = (18446744073709551614,18446744073709551615,"18446744073709551616","18446744073709551617");
}

my @tests = (
  [$mm1, 1, $m],   # native => native
  [1, $mm1, $m],
  [$m, 1, $mp1],   # native => bigint
  [1, $m, $mp1],
  [$mp1, 1, $mp2], # bigint => bigint
  [1, $mp1, $mp2],
  [$mp1, -1, $m],  # bigint => native
  [-1, $mp1, $m],
);

plan tests => 4 + scalar(@tests) + 1 + 1;

is(should_be_ref($mm1), 0, "$mm1 should be a UV");
is(should_be_ref($m),   0, "$m should be a UV");
is(should_be_ref($mp1), 1, "$mp1 should be a bigint");
is(should_be_ref($mp2), 1, "$mp2 should be a bigint");

for my $t (@tests) {
  my($a,$b,$r) = @$t;
  my $expr = ($r eq $mp1 || $r eq $mp2) ? "BIGINT" : "NATIVE";
  #my $expr = should_be_ref($r);

  my $sum = addint($a,$b);

  my $gotr = ref($sum) ? "BIGINT" : "NATIVE";
  ok("$sum" eq $r && $gotr eq $expr, "$sum [$gotr]  expect $r [$expr]");
}


sub addsub1 {
  my($a,$b) = @_;
  my $sum = 0;
  if    ($b == -1) { $sum = sub1int($a); }
  elsif ($b ==  1) { $sum = add1int($a); }
  elsif ($a ==  1) { $sum = add1int($b); }
  elsif ($a == -1) { $sum = sub1int($b); }
  my $sum_is_ref = ref($sum) ? 1 : 0;
  [$sum, $sum_is_ref];
}

is_deeply( [map { addsub1($_->[0],$_->[1]) } @tests],
           [map { [$_->[2], should_be_ref($_->[2])] } @tests],
           "add1int and sub1int return correct types" );


my $m2 = addint($m,$m);
my @mulres = (
  [mulint($m,1),$m,0],      # native => native
  [mulint($m,2),$m2,1],     # native => bigint
  [mulint($mp1,1),$mp1,1],  # bigint => bigint
  [mulint($mp1,0),0,0],     # bigint => native
);
is_deeply( [map { $_->[0], ref($_->[0]) ? 1 : 0 } @mulres],
           [map { $_->[1], $_->[2] } @mulres],
           "mulint returns correct types" );

sub should_be_ref {
  my($n,$min,$max) = @_;
  $min = $intmin unless defined $min;
  $max = $intmax unless defined $max;

  # If n meets these conditions, we're good regardless of integer form.
  return 0 if $n > $min && $n < $max;

  # If a number is in string form (not a bigint), then Perl will say:
  #    ~0+2048 == ~0
  #
  # In decreasing performance order:
  #   1) Use a module like Math::Prime::Util.  That's us.
  #   2) Compare as strings.  Requires a bit of code.
  #   3) Turn it into a Math::BigInt so things work right.
  $n = Math::BigInt->new("$n");
  return ($n > $max || $n < $min) ? 1 : 0;
}
