#!/usr/bin/env perl
use strict;
use warnings;

use Config;
use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/
  prime_get_config prime_set_config
  toint absint addint subint mulint divint cdivint modint
  divrem tdivrem fdivrem cdivrem
  addmod submod mulmod divmod powmod muladdmod mulsubmod
  vecsum vecprod vecprefixsum vecreduce
/;

my @classes = ('Math::BigInt');
push @classes, 'Math::GMPz' if eval { require Math::GMPz; 1 };
push @classes, 'Math::GMP'  if eval { require Math::GMP;  1 };

note('available bigint classes: ' . join(', ', @classes));
my $old_bigint = prime_get_config()->{bigintclass};

my ($ivmin, $ivmax, $uvmax) = $Config{ivsize} == 4
  ? qw/-2147483648 2147483647 4294967295/
  : qw/-9223372036854775808 9223372036854775807 18446744073709551615/;

my $x = '1180591620717411303424';              # 2^70
my $y = '2503155504993241601315571986085849';  # 3^70
my $z = '1104427674243920646305299201';        # 7^32
my $w = '340282366920938463463374607431768211457'; # 2^128 + 1

sub bi {
  return Math::BigInt->new("$_[0]");
}

sub calc {
  my ($a, $method, @args) = @_;
  return bi($a)->$method(@args)->bstr;
}

sub make_object {
  my ($class, $value) = @_;
  return $class->new("$value");
}

sub expected_ref {
  my ($value, $class) = @_;
  my $n = bi($value);
  return '' if $n >= bi($ivmin) && $n <= bi($uvmax);
  return $class;
}

sub signature {
  my ($value) = @_;
  return [defined($value) ? "$value" : '<undef>', ref($value) || ''];
}

sub expected_signature {
  my ($value, $class) = @_;
  return ["$value", expected_ref($value, $class)];
}

sub check_scalar {
  my ($got, $expected, $class, $name) = @_;
  is_deeply(signature($got), expected_signature($expected, $class), $name);
}

sub check_list {
  my ($got, $expected, $class, $name) = @_;
  my @got = map { signature($_) } @$got;
  my @expected = map { expected_signature($_, $class) } @$expected;
  is_deeply(\@got, \@expected, $name);
}

my %expected = (
  add       => calc($x, 'badd', $y),
  sub       => calc("-$y", 'bsub', $x),
  mul       => calc($x, 'bmul', -31),
  mod       => calc("-$y", 'bmod', $x),
  addmod    => calc(calc($x, 'badd', $y), 'bmod', $z),
  submod    => calc(calc($x, 'bsub', $y), 'bmod', $z),
  mulmod    => calc(calc($x, 'bmul', $y), 'bmod', $z),
  divmod    => '814177443435415951196839160',
  powmod    => bi($x)->bmodpow(37, bi($z))->bstr,
  muladdmod => calc(calc(calc($x, 'bmul', $y), 'badd', 2), 'bmod', $z),
  mulsubmod => calc(calc(calc($x, 'bmul', $y), 'bsub', 2), 'bmod', $z),
  wide_add  => calc($w, 'badd', $x),
  wide_mulmod => calc(calc($x, 'bmul', $y), 'bmod', $w),
  wide_powmod => bi($x)->bmodpow(37, bi($w))->bstr,
);

my $uvplus1 = calc($uvmax, 'badd', 1);
my $ivminus1 = calc($ivmin, 'bsub', 1);
my $absivmin = bi($ivmin)->babs->bstr;

my %divrem_expected = (
  divrem  => [ '-2120255184831', '883150645529757075495' ],
  tdivrem => [ '-2120255184830', '-297440975187654227929' ],
  fdivrem => [ '-2120255184831', '883150645529757075495' ],
  cdivrem => [ '-2120255184830', '-297440975187654227929' ],
);
my %divrem_func = (
  divrem  => \&divrem,
  tdivrem => \&tdivrem,
  fdivrem => \&fdivrem,
  cdivrem => \&cdivrem,
);

for my $output_class (@classes) {
  subtest "configured output class $output_class" => sub {
    prime_set_config(bigint => $output_class);
    is(prime_get_config()->{bigintclass}, $output_class,
       'requested output class is active');

    for my $input_class (@classes) {
      subtest "input objects $input_class" => sub {
        my $obj = sub { make_object($input_class, $_[0]) };

        check_scalar(toint($obj->($x)), $x, $output_class,
                     'toint canonicalizes a bigint object');
        check_scalar(toint($obj->($w)), $w, $output_class,
                     'toint canonicalizes a greater-than-128-bit object');
        check_scalar(toint($obj->($ivmin)), $ivmin, $output_class,
                     'toint canonicalizes IV_MIN');
        check_scalar(absint($obj->($ivmin)), $absivmin, $output_class,
                     'absint handles IV_MIN object');

        check_scalar(addint($obj->($x), $obj->($y)),
                     $expected{add}, $output_class, 'addint');
        check_scalar(subint($obj->("-$y"), $obj->($x)),
                     $expected{sub}, $output_class, 'subint');
        check_scalar(mulint($obj->($x), -31),
                     $expected{mul}, $output_class, 'mulint');
        check_scalar(divint($obj->("-$y"), $obj->($x)),
                     '-2120255184831', $output_class, 'divint');
        check_scalar(cdivint($obj->("-$y"), $obj->($x)),
                     '-2120255184830', $output_class, 'cdivint');
        check_scalar(modint($obj->("-$y"), $obj->($x)),
                     $expected{mod}, $output_class, 'modint');

        check_scalar(addint($obj->($uvmax), 1),
                     $uvplus1, $output_class, 'UV_MAX plus one');
        check_scalar(subint($obj->($ivmin), 1),
                     $ivminus1, $output_class, 'IV_MIN minus one');
        check_scalar(subint($obj->($uvplus1), $obj->($uvmax)),
                     1, $output_class, 'bigint result contracts to native');
        check_scalar(addint($obj->($w), $obj->($x)),
                     $expected{wide_add}, $output_class,
                     'addint with greater-than-128-bit input');

        for my $name (qw/divrem tdivrem fdivrem cdivrem/) {
          my @got = $divrem_func{$name}->($obj->("-$y"), $obj->($x));
          check_list(\@got, $divrem_expected{$name}, $output_class, $name);
        }

        check_scalar(addmod($obj->($x), $obj->($y), $obj->($z)),
                     $expected{addmod}, $output_class, 'addmod');
        check_scalar(submod($obj->($x), $obj->($y), $obj->($z)),
                     $expected{submod}, $output_class, 'submod');
        check_scalar(mulmod($obj->($x), $obj->($y), $obj->($z)),
                     $expected{mulmod}, $output_class, 'mulmod');
        check_scalar(divmod($obj->($x), $obj->($y), $obj->($z)),
                     $expected{divmod}, $output_class, 'divmod');
        check_scalar(powmod($obj->($x), 37, $obj->($z)),
                     $expected{powmod}, $output_class, 'powmod');
        check_scalar(muladdmod($obj->($x), $obj->($y), 2, $obj->($z)),
                     $expected{muladdmod}, $output_class, 'muladdmod');
        check_scalar(mulsubmod($obj->($x), $obj->($y), 2, $obj->($z)),
                     $expected{mulsubmod}, $output_class, 'mulsubmod');
        check_scalar(mulmod($obj->($x), $obj->($y), $obj->($w)),
                     $expected{wide_mulmod}, $output_class,
                     'mulmod with greater-than-128-bit modulus');
        check_scalar(powmod($obj->($x), 37, $obj->($w)),
                     $expected{wide_powmod}, $output_class,
                     'powmod with greater-than-128-bit modulus');

        my @sum_values = ($obj->($x), $obj->("-$y"), $obj->($z), 17);
        my $sum_expected = calc(calc(calc($x, 'bsub', $y), 'badd', $z),
                                'badd', 17);
        check_scalar(vecsum(@sum_values), $sum_expected, $output_class,
                     'vecsum');
        is_deeply([map { "$_" } @sum_values], [$x, "-$y", $z, 17],
                  'vecsum does not modify inputs');

        my @prod_values = ($obj->($x), -31, 7);
        my $prod_expected = calc(calc($x, 'bmul', -31), 'bmul', 7);
        check_scalar(vecprod(@prod_values), $prod_expected, $output_class,
                     'vecprod');

        my @prefix_values = ($obj->($ivmax), $obj->(1),
                             $obj->($x), $obj->("-$x"));
        my @prefix_expected = (
          $ivmax,
          calc($ivmax, 'badd', 1),
          calc(calc($ivmax, 'badd', 1), 'badd', $x),
          calc($ivmax, 'badd', 1),
        );
        my @prefix = vecprefixsum(@prefix_values);
        check_list(\@prefix, \@prefix_expected, $output_class,
                   'vecprefixsum');
        is(scalar(vecprefixsum(@prefix_values)), scalar(@prefix_values),
           'scalar vecprefixsum returns count');

        my $reduced = vecreduce(sub { addint($a, $b) },
                                $obj->($x), $obj->("-$y"), $obj->($z));
        check_scalar($reduced,
                     calc(calc($x, 'bsub', $y), 'badd', $z),
                     $output_class, 'vecreduce with addint callback');
      };
    }
  };
}

prime_set_config(bigint => $old_bigint) if defined $old_bigint;
done_testing();
