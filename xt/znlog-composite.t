#!/usr/bin/env perl
use strict;
use warnings;

use Config;
use Test::More;
use Math::Prime::Util qw/
  prime_get_config prime_set_config
  znlog znprimroot powmod is_prime
/;
use Math::Prime::Util::PP;

my $config = prime_get_config();
my $old_gmp = $config->{gmp};
prime_set_config(gmp => 0);
END { prime_set_config(gmp => $old_gmp) if defined $old_gmp }

sub brute_logs {
  my ($g, $n) = @_;
  my %logs;
  my $value = 1 % $n;
  my $k = 0;
  while (!exists $logs{$value}) {
    $logs{$value} = $k++;
    $value = ($value * $g) % $n;
  }
  return \%logs;
}

sub check_result {
  my ($errors, $name, $func, $a, $g, $n, $exists) = @_;
  my $k = $func->($a, $g, $n);
  if ($exists) {
    if (!defined($k)) {
      push @$errors, "$name znlog($a,$g,$n) returned undef";
    } elsif (powmod($g, $k, $n) != $a % $n) {
      push @$errors, "$name znlog($a,$g,$n) returned invalid k=$k";
    }
  } elsif (defined($k)) {
    push @$errors, "$name znlog($a,$g,$n) returned unexpected k=$k";
  }
}

sub report_errors {
  my ($errors) = @_;
  my $nshow = @$errors < 20 ? scalar(@$errors) : 20;
  diag($errors->[$_]) for 0 .. $nshow-1;
  diag((@$errors - $nshow) . ' additional errors omitted')
    if @$errors > $nshow;
}

subtest 'exhaustive small composite moduli' => sub {
  my @errors;
  N:
  for my $n (4 .. 80) {
    next if is_prime($n);
    for my $g (0 .. $n-1) {
      my $logs = brute_logs($g, $n);
      for my $a (0 .. $n-1) {
        my $exists = exists $logs->{$a};
        check_result(\@errors, 'XS', \&znlog, $a, $g, $n, $exists);
        check_result(\@errors, 'PP', \&Math::Prime::Util::PP::znlog,
                     $a, $g, $n, $exists);
        last N if @errors >= 20;
      }
    }
  }
  is(scalar(@errors), 0,
     'XS and PP solve exactly the reachable small composite residues');
  report_errors(\@errors);
};

subtest 'larger composite moduli' => sub {
  my @errors;
  for my $n (81, 125, 128, 225, 243, 256, 289, 343, 625, 729, 1024) {
    my %seen_g;
    for my $g (0, 1, 2, 3, 5, int($n/2), $n-1) {
      next if $seen_g{$g}++;
      my $logs = brute_logs($g, $n);
      my %test_a = map { $_ => 1 }
                   (0, 1, 2, 3, 5, int($n/3), int($n/2), $n-1,
                    keys %$logs);
      for my $a (keys %test_a) {
        my $exists = exists $logs->{$a};
        check_result(\@errors, 'XS', \&znlog, $a, $g, $n, $exists);
        check_result(\@errors, 'PP', \&Math::Prime::Util::PP::znlog,
                     $a, $g, $n, $exists);
      }
    }
  }
  is(scalar(@errors), 0,
     'XS and PP handle larger coprime and non-coprime composite cases');
  report_errors(\@errors);
};

subtest 'proper subgroups' => sub {
  my @errors;
  my $p = 1009;
  my $root = znprimroot($p);
  my $g = powmod($root, 8, $p);
  my $logs = brute_logs($g, $p);
  for my $a (0 .. $p-1) {
    my $exists = exists $logs->{$a};
    check_result(\@errors, 'XS', \&znlog, $a, $g, $p, $exists);
    check_result(\@errors, 'PP', \&Math::Prime::Util::PP::znlog,
                 $a, $g, $p, $exists);
  }
  is(scalar(@errors), 0,
     'XS and PP distinguish residues inside and outside a proper subgroup');
  report_errors(\@errors);
};

subtest 'odd prime-power modulus and Montgomery path' => sub {
  my $p = 1009;
  my $n = $p * $p;
  my $root = znprimroot($n);
  my $g = powmod($root, 2, $n);
  my $a = powmod($g, 777, $n);
  my $k = znlog($a, $g, $n);
  ok(defined($k), 'log exists modulo an odd prime power');
  is(powmod($g, $k, $n), $a, 'prime-power result reconstructs the residue');
  ok(!defined(znlog($root, $g, $n)),
     'element outside the square subgroup has no log');

  my $negative = znlog(-$a, -$g, -$n);
  if (defined $negative) {
    is(powmod(-$g, $negative, $n), (-$a) % $n,
       'negative arguments are reduced correctly');
  } else {
    pass('negative argument combination has no logarithm');
  }
};

SKIP: {
  skip 'native Pollard-rho Montgomery path requires XS uint128_t support', 2
    unless $Config{uvsize} >= 8 && $config->{xs}
           && Math::Prime::Util::_XS_has_uint128();

  # p = 2q+1 with q prime.  The order is large enough to bypass BSGS,
  # while this fixed logarithm is quick enough for routine testing.
  my $p = 10_000_000_000_259;
  my $g = 4;
  my $expected = 123_456_789;
  my $a = powmod($g, $expected, $p);
  my $got = znlog($a, $g, $p);
  is("$got", "$expected", 'Pollard-rho finds a log in a large prime subgroup');
  is(powmod($g, $got, $p), $a, 'Pollard-rho result reconstructs the residue');
}

SKIP: {
  skip 'native Montgomery boundary cases require XS uint128_t support', 6
    unless $Config{uvsize} >= 8 && $config->{xs}
           && Math::Prime::Util::_XS_has_uint128();

  my $n = '18446744073709551557';  # Largest prime below 2^64.
  for my $k (0, 1, 17, 127, 257, 4097) {
    my $a = powmod(2, $k, $n);
    my $got = znlog($a, 2, $n);
    ok(defined($got) && powmod(2, $got, $n) == $a,
       "near-UV_MAX modulus, exponent $k");
  }
}

done_testing();
