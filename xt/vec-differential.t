#!/usr/bin/env perl
use strict;
use warnings;

use Config;
use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/
  vecprefixsum vecsort vecrsort vecrsorti
/;
use Math::Prime::Util::PP;

my $NCASES = $ENV{MPU_VEC_DIFFERENTIAL_CASES} || 500;

my ($ivmin, $ivmax, $uvmax) = $Config{ivsize} == 4
  ? qw/-2147483648 2147483647 4294967295/
  : qw/-9223372036854775808 9223372036854775807 18446744073709551615/;

sub _bi {
  return Math::BigInt->new("$_[0]");
}

my $ivbelow = _bi($ivmin)->bsub(1)->bstr;
my $ivabove = _bi($ivmax)->badd(1)->bstr;
my $uvabove = _bi($uvmax)->badd(1)->bstr;
my @pool = (
  $ivbelow, $ivmin, -1000003, -17, -1, 0, 1, 2, 17, 1000003,
  $ivmax, $ivabove, $uvmax, $uvabove,
  '1180591620717411303424', '-1180591620717411303424',
  '340282366920938463463374607431768211457',
);

sub _prefix_expected {
  my @values = @_;
  my $sum = _bi(0);
  my @out;
  for my $value (@values) {
    $sum->badd("$value");
    push @out, $sum->bstr;
  }
  return \@out;
}

sub _reverse_sorted {
  return [
    map { "$_" }
    sort { _bi($b)->bcmp(_bi($a)) } @_
  ];
}

sub _sig {
  return join("\0", map { "$_" } @_);
}

subtest 'vecprefixsum differential and boundary cases' => sub {
  my @fail;
  for my $case (0 .. $NCASES - 1) {
    my $len = ($case * 13) % 31;
    my @in = map { $pool[($case * 11 + $_ * 7) % @pool] } 0 .. $len - 1;
    my @before = map { "$_" } @in;
    my $expected = _prefix_expected(@in);

    my @list = map { "$_" } vecprefixsum(@in);
    my @aref = map { "$_" } vecprefixsum(\@in);
    my @pp = map { "$_" } Math::Prime::Util::PP::vecprefixsum(@in);
    push @fail, "case $case list" if _sig(@list) ne _sig(@$expected);
    push @fail, "case $case aref" if _sig(@aref) ne _sig(@$expected);
    push @fail, "case $case PP" if _sig(@pp) ne _sig(@$expected);
    push @fail, "case $case scalar"
      if scalar(vecprefixsum(@in)) != $len ||
         scalar(vecprefixsum(\@in)) != $len;
    push @fail, "case $case mutated input" if _sig(@in) ne _sig(@before);
    last if @fail >= 10;
  }
  is_deeply(\@fail, [], "$NCASES prefix-sum cases agree");
};

subtest 'vecrsort and vecrsorti differential' => sub {
  my @fail;
  for my $case (0 .. $NCASES - 1) {
    my $len = ($case * 17) % 31;
    my @in = map { $pool[($case * 5 + $_ * 11) % @pool] } 0 .. $len - 1;
    my @before = map { "$_" } @in;
    my $expected = _reverse_sorted(@in);

    my @list = map { "$_" } vecrsort(@in);
    my @aref = map { "$_" } vecrsort(\@in);
    my @reverse_ascending = map { "$_" } reverse vecsort(@in);
    my @pp = map { "$_" } Math::Prime::Util::PP::vecrsort(@in);
    push @fail, "case $case list" if _sig(@list) ne _sig(@$expected);
    push @fail, "case $case aref" if _sig(@aref) ne _sig(@$expected);
    push @fail, "case $case reverse vecsort"
      if _sig(@reverse_ascending) ne _sig(@$expected);
    push @fail, "case $case PP" if _sig(@pp) ne _sig(@$expected);
    push @fail, "case $case scalar" if scalar(vecrsort(@in)) != $len;
    push @fail, "case $case non-in-place mutation"
      if _sig(@in) ne _sig(@before);

    my @inplace = @in;
    my $ref = vecrsorti(\@inplace);
    push @fail, "case $case return identity" if $ref != \@inplace;
    push @fail, "case $case in-place result"
      if _sig(@inplace) ne _sig(@$expected);
    last if @fail >= 10;
  }
  is_deeply(\@fail, [], "$NCASES reverse-sort cases agree");
};

subtest 'vecrsorti magical array fallback' => sub {
  require Tie::Array;
  tie my @values, 'Tie::StdArray';
  @values = (3, -1, 7, 3, $uvabove, $ivbelow);
  my $ref = vecrsorti(\@values);
  is($ref, \@values, 'returns the tied array reference');
  is_deeply([map { "$_" } @values],
            _reverse_sorted(3, -1, 7, 3, $uvabove, $ivbelow),
            'sorts a tied array through the PP fallback');
  untie @values;
};

done_testing();
