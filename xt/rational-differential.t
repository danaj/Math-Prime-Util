#!/usr/bin/env perl
use strict;
use warnings;

use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/
  bestrational convergents from_contfrac
/;
use Math::Prime::Util::PP;

my $NCONV = $ENV{MPU_RATIONAL_CONVERGENTS} || 300;
my $NBEST = $ENV{MPU_RATIONAL_BEST}        || 500;

sub _bi {
  return Math::BigInt->new("$_[0]");
}

sub _string_pairs {
  return [map { [map { "$_" } @$_] } @_];
}

sub _expected_convergents {
  my @cf = @_;
  my ($p0, $q0) = (_bi(1), _bi(0));
  my ($p1, $q1) = (_bi($cf[0]), _bi(1));
  my @out = ([$p1->bstr, $q1->bstr]);

  for my $i (1 .. $#cf) {
    my $p2 = $p1->copy->bmul("$cf[$i]")->badd($p0);
    my $q2 = $q1->copy->bmul("$cf[$i]")->badd($q0);
    push @out, [$p2->bstr, $q2->bstr];
    ($p0, $q0, $p1, $q1) = ($p1, $q1, $p2, $q2);
  }
  return \@out;
}

sub _gcd {
  my ($a, $b) = @_;
  $a = -$a if $a < 0;
  while ($b) {
    ($a, $b) = ($b, $a % $b);
  }
  return $a;
}

sub _best_is_valid {
  my ($num, $den, $bound, $got) = @_;
  my ($rp, $rq) = @$got;
  return "invalid denominator $rq" if $rq < 1 || $rq > $bound;
  return "result is not reduced" if _gcd($rp, $rq) != 1;

  my $rerr = abs($num * $rq - $rp * $den);
  for my $q (1 .. $bound) {
    my $scaled = $num * $q;
    my $p = int($scaled / $den);
    $p-- if $scaled < 0 && $scaled % $den;
    for my $candidate ($p, $p + 1) {
      my $err = abs($num * $q - $candidate * $den);
      return "$candidate/$q is closer than $rp/$rq"
        if $err * $rq < $rerr * $q;
    }
  }
  return;
}

subtest 'convergents recurrence and XS/PP agreement' => sub {
  my @fail;
  for my $case (0 .. $NCONV - 1) {
    my $len = 1 + (($case * 7) % 12);
    my @cf = ((($case * 37) % 101) - 50);
    push @cf, 1 + (($case * 43 + $_ * 61) % 1000) for 1 .. $len - 1;
    $cf[-1] = '1180591620717411303424' if $len > 1 && $case % 17 == 0;

    my $expected = _expected_convergents(@cf);
    my @xs = convergents(@cf);
    my @pp = Math::Prime::Util::PP::convergents(@cf);
    my $xs = _string_pairs(@xs);
    my $pp = _string_pairs(@pp);

    push @fail, "case $case XS recurrence"
      if "@{$xs->[0]}" ne "@{$expected->[0]}" ||
         join(';', map { "@$_" } @$xs) ne
         join(';', map { "@$_" } @$expected);
    push @fail, "case $case PP recurrence"
      if join(';', map { "@$_" } @$pp) ne
         join(';', map { "@$_" } @$expected);

    my $xcount = scalar convergents(@cf);
    my $pcount = Math::Prime::Util::PP::convergents(@cf);
    push @fail, "case $case scalar count" if $xcount != $len || $pcount != $len;

    my @last = map { "$_" } from_contfrac(@cf);
    push @fail, "case $case final convergent"
      if "@last" ne "@{$expected->[-1]}";

    for my $i (1 .. $#$expected) {
      my $det = _bi($expected->[$i][0])
                  ->bmul($expected->[$i-1][1])
                  ->bsub(_bi($expected->[$i-1][0])
                           ->bmul($expected->[$i][1]));
      my $want = $i % 2 ? 1 : -1;
      push @fail, "case $case determinant $i"
        if $det->bcmp($want) != 0;
    }
    last if @fail >= 10;
  }
  is_deeply(\@fail, [], "$NCONV continued fractions satisfy all properties");
};

subtest 'bestrational brute-force optimality' => sub {
  my @fail;
  for my $case (0 .. $NBEST - 1) {
    my $den = 1 << (1 + ($case % 10));
    my $num = (($case * 7919) % 4001) - 2000;
    my $bound = 1 + (($case * 97) % 80);
    my $x = $num / $den;

    my @xs = map { 0 + "$_" } bestrational($x, $bound);
    my @pp = map { 0 + "$_" }
             Math::Prime::Util::PP::bestrational($x, $bound);
    my $xerr = _best_is_valid($num, $den, $bound, \@xs);
    my $perr = _best_is_valid($num, $den, $bound, \@pp);
    push @fail, "case $case XS: $xerr" if defined $xerr;
    push @fail, "case $case PP: $perr" if defined $perr;

    my $xe = abs($num * $xs[1] - $xs[0] * $den);
    my $pe = abs($num * $pp[1] - $pp[0] * $den);
    push @fail, "case $case XS/PP errors differ"
      if $xe * $pp[1] != $pe * $xs[1];
    last if @fail >= 10;
  }
  is_deeply(\@fail, [], "$NBEST exact dyadic targets are optimal");
};

done_testing();
