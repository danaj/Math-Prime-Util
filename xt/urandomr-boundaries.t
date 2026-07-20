#!/usr/bin/env perl
use strict;
use warnings;

use Config;
use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/
  csrand prime_get_config prime_set_config urandomr
/;

my $NSAMPLES = $ENV{MPU_URANDOMR_SAMPLES} || 64;

sub _bi {
  return Math::BigInt->new("$_[0]");
}

sub _add {
  my ($a, $b) = @_;
  return _bi($a)->badd($b)->bstr;
}

my ($ivmin, $ivmax, $uvmax) = $Config{ivsize} == 4
  ? qw/-2147483648 2147483647 4294967295/
  : qw/-9223372036854775808 9223372036854775807 18446744073709551615/;
my $p80 = _bi(2)->bpow(80)->bstr;
my $n80 = _bi($p80)->bneg->bstr;
my $p128 = _bi(2)->bpow(128)->bstr;
my $n128 = _bi($p128)->bneg->bstr;

my @ranges = (
  ['singleton zero',       0,                    0],
  ['small positive',       3,                    8],
  ['cross zero',          -17,                  17],
  ['near IV_MIN',          $ivmin,               _add($ivmin, 7)],
  ['near IV_MAX',          _add($ivmax, -7),     $ivmax],
  ['near UV_MAX',          _add($uvmax, -7),     $uvmax],
  ['full native unsigned', 0,                    $uvmax],
  ['full signed/native',   $ivmin,               $uvmax],
  ['80-bit signed',        $n80,                 $p80],
  ['above 128-bit',        _add($p128, 1000),    _add($p128, 1100)],
  ['below -128-bit',       _add($n128, -1100),   _add($n128, -1000)],
);

sub _in_range {
  my ($value, $lo, $hi) = @_;
  my $v = _bi($value);
  return $v->bcmp($lo) >= 0 && $v->bcmp($hi) <= 0;
}

subtest 'range invariants' => sub {
  csrand(87178291199);
  my @fail;
  for my $range (@ranges) {
    my ($name, $lo, $hi) = @$range;
    for (1 .. $NSAMPLES) {
      my $value = urandomr($lo, $hi);
      if (!defined($value) || "$value" !~ /^-?\d+\z/ ||
          !_in_range($value, $lo, $hi)) {
        push @fail, "$name: " . (defined($value) ? "$value" : '<undef>');
        last;
      }
    }
  }
  is_deeply(\@fail, [], "$NSAMPLES samples from every range stay in range");
};

subtest 'degenerate and reversed ranges' => sub {
  is("".urandomr($p128, $p128), $p128,
     'large singleton returns its endpoint');
  is(urandomr(7, 3), undef, 'native reversed range returns undef');
  is(urandomr(_add($p128, 1), $p128), undef,
     'bigint reversed range returns undef');
};

subtest 'csrand reproducibility' => sub {
  csrand(1234567);
  my @a = map { "".urandomr($n80, $p80) } 1 .. 50;
  csrand(1234567);
  my @b = map { "".urandomr($n80, $p80) } 1 .. 50;
  is_deeply(\@a, \@b, 'reseeded range stream is reproducible');
};

subtest 'configured bigint return classes' => sub {
  my $old = prime_get_config()->{bigintclass};
  my @classes = ('Math::BigInt');
  push @classes, 'Math::GMPz' if eval { require Math::GMPz; 1 };
  push @classes, 'Math::GMP'  if eval { require Math::GMP; 1 };

  for my $class (@classes) {
    prime_set_config(bigint => $class);
    my $value = urandomr(_add($p128, 1000), _add($p128, 1100));
    is(ref($value), $class, "$class is used for an out-of-native result");
    ok(_in_range($value, _add($p128, 1000), _add($p128, 1100)),
       "$class result is in range");
  }
  prime_set_config(bigint => $old) if defined $old;
};

done_testing();
