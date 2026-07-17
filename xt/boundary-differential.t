#!/usr/bin/env perl
use strict;
use warnings;

use Config;
use Math::BigInt;
use Test::More;
use Math::Prime::Util qw/
  prime_set_config urandomb urandomm
  divrem tdivrem fdivrem cdivrem
  gcdext kronecker powmod
  lucasu lucasv lucasuv
  falling_factorial rising_factorial
/;
use Math::Prime::Util::PP ();

my $have_gmp = eval {
  require Math::Prime::Util::GMP;
  1;
};

# Keep the MPU and PP paths from silently dispatching to GMP.  Direct GMP
# calls are still made below when the backend is available.
prime_set_config(gmp => 0);

my @backends = (
  [XS  => 'Math::Prime::Util'],
  [PP  => 'Math::Prime::Util::PP'],
);
push @backends, [GMP => 'Math::Prime::Util::GMP'] if $have_gmp;

my $bits = 8 * $Config{uvsize};
my $p2   = Math::BigInt->new(2);
my $ivmin = $p2->copy->bpow($bits-1)->bneg->bstr;
my $ivmax = $p2->copy->bpow($bits-1)->bdec->bstr;
my $uvmax = $p2->copy->bpow($bits)->bdec->bstr;

diag("native width: $bits bits; backends: ",
     join(", ", map { $_->[0] } @backends));

sub _bi {
  return Math::BigInt->new("$_[0]");
}

sub _uniq {
  my %seen;
  return grep { !$seen{"$_"}++ } @_;
}

sub _around {
  my ($center, @offsets) = @_;
  return map { _bi($center)->badd($_)->bstr } @offsets;
}

sub _show_args {
  return join(", ", map { defined($_) ? "$_" : '<undef>' } @_);
}

sub _canon {
  return scalar(@_) . ':' .
         join("\x1E", map { defined($_) ? "$_" : '<undef>' } @_);
}

sub _results_for {
  my ($fn, $args, $fail) = @_;
  my @results;

  for my $backend (@backends) {
    my ($label, $package) = @$backend;
    my $code = $package->can($fn);
    if (!$code) {
      push @$fail, "$label does not provide $fn";
      next;
    }

    my (@result, $ok, $error);
    {
      local $@;
      $ok = eval {
        @result = $code->(@$args);
        1;
      };
      $error = $@;
    }
    if (!$ok) {
      chomp $error;
      push @$fail, "$label $fn(" . _show_args(@$args) . ") died: $error";
      next;
    }
    push @results, [$label, \@result, _canon(@result)];
  }
  return \@results;
}

sub _compare_exact {
  my ($fn, $args, $fail) = @_;
  my $results = _results_for($fn, $args, $fail);
  return $results if @$results < 2;

  my ($expected_label, undef, $expected) = @{$results->[0]};
  for my $result (@$results[1 .. $#$results]) {
    my ($label, undef, $got) = @$result;
    push @$fail,
      "$fn(" . _show_args(@$args) . "): $expected_label [$expected] != $label [$got]"
      if $got ne $expected;
  }
  return $results;
}

sub _finish {
  my ($description, $checks, $fail) = @_;
  ok(!@$fail, "$description ($checks cases)");
  diag($_) for @$fail[0 .. ($#$fail < 4 ? $#$fail : 4)];
}

my @edge = _uniq(
  -3 .. 3,
  _around($ivmin, -2 .. 2),
  _around($ivmax, -2 .. 2),
  _around($uvmax, -2 .. 2),
  _around("-$uvmax", -2 .. 2),
);

Math::Prime::Util::srand(0x51A7D1FF);

sub _random_iv {
  my $n = "" . urandomb($bits-1);
  return urandomm(2) && $n ne '0' ? "-$n" : $n;
}

sub _random_native {
  return "" . urandomb($bits) if urandomm(2);
  return _random_iv();
}

my @random_pairs;
for (1 .. 100) {
  my $a = _random_native();
  my $b = _random_native();
  $b = 1 if $b eq '0';
  push @random_pairs, [$a, $b];
}

subtest 'division modes' => sub {
  my (@fail, @pairs);
  my @divisors = _uniq(
    -3, -2, -1, 1, 2, 3,
    _around($ivmin, -1 .. 1),
    _around($ivmax, -1 .. 1),
    _around($uvmax, -1 .. 1),
  );
  push @pairs, map {
    my $a = $_;
    map { [$a, $_] } grep { $_ ne '0' } @divisors
  } @edge;
  push @pairs, @random_pairs;

  my $checks = 0;
  for my $fn (qw/divrem tdivrem fdivrem cdivrem/) {
    for my $pair (@pairs) {
      my ($a, $b) = @$pair;
      my $results = _compare_exact($fn, $pair, \@fail);
      $checks++;

      for my $result (@$results) {
        my ($label, $values) = @$result;
        if (@$values != 2) {
          push @fail, "$label $fn($a,$b) returned " . scalar(@$values) .
                      " values, expected 2";
          next;
        }
        my ($q, $r) = map { "$_" } @$values;
        my $A = _bi($a);
        my $B = _bi($b);
        my $R = _bi($r);
        my $reconstructed = $B->copy->bmul($q)->badd($R);
        push @fail, "$label $fn($a,$b) violates a=b*q+r"
          if $reconstructed->bcmp($A) != 0;
        push @fail, "$label $fn($a,$b) has |r| >= |b|"
          if $R->copy->babs->bcmp($B->copy->babs) >= 0;
        push @fail, "$label divrem($a,$b) has negative remainder $r"
          if $fn eq 'divrem' && $R->is_neg;
        push @fail, "$label tdivrem($a,$b) has wrong-sign remainder $r"
          if $fn eq 'tdivrem' && !$R->is_zero &&
             ($R->is_neg != $A->is_neg);
        push @fail, "$label fdivrem($a,$b) has wrong-sign remainder $r"
          if $fn eq 'fdivrem' && !$R->is_zero &&
             ($R->is_neg != $B->is_neg);
        push @fail, "$label cdivrem($a,$b) has wrong-sign remainder $r"
          if $fn eq 'cdivrem' && !$R->is_zero &&
             ($R->is_neg == $B->is_neg);
      }
    }
  }
  _finish('all division modes agree and reconstruct', $checks, \@fail);
  done_testing();
};

subtest 'gcdext and Bezout identities' => sub {
  my (@fail, @pairs);
  my @partners = _uniq(-3 .. 3, $ivmin, $ivmax, $uvmax);
  push @pairs, map {
    my $a = $_;
    map { [$a, $_] } @partners
  } @edge;
  push @pairs, @random_pairs;

  my $checks = 0;
  for my $pair (@pairs) {
    my ($a, $b) = @$pair;
    my $results = _results_for('gcdext', $pair, \@fail);
    my $expected_d = _bi($a)->babs->bgcd(_bi($b)->babs)->bstr;
    my $first_d;
    for my $result (@$results) {
      my ($label, $values) = @$result;
      if (@$values != 3) {
        push @fail, "$label gcdext($a,$b) returned " . scalar(@$values) .
                    " values, expected 3";
        next;
      }
      my ($u, $v, $d) = map { "$_" } @$values;
      my $bezout = _bi($a)->bmul($u)->badd(_bi($b)->bmul($v));
      push @fail, "$label gcdext($a,$b) violates Bezout identity"
        if $bezout->bcmp($d) != 0;
      push @fail, "$label gcdext($a,$b) returned gcd $d, expected $expected_d"
        if $d ne $expected_d;
      $first_d = $d if !defined $first_d;
      push @fail, "gcdext($a,$b) gcd differs between backends"
        if defined($first_d) && $d ne $first_d;
    }
    $checks++;
  }
  _finish('gcdext satisfies Bezout identities', $checks, \@fail);
  done_testing();
};

subtest 'kronecker' => sub {
  my (@fail, @pairs);
  my @partners = _uniq(-3 .. 3, $ivmin, $ivmax, $uvmax);
  push @pairs, map {
    my $a = $_;
    map { [$a, $_] } @partners
  } @edge;
  push @pairs, @random_pairs;

  _compare_exact('kronecker', $_, \@fail) for @pairs;
  _finish('kronecker agrees across backends', scalar(@pairs), \@fail);
  done_testing();
};

subtest 'negative-exponent powmod' => sub {
  my @fail;
  my @bases = _uniq(-3 .. 3, $ivmin, $ivmax, $uvmax, @edge[0 .. 4]);
  my @exponents = _uniq(
    -1, -2, -3,
    _around($ivmin, 0 .. 2),
    "-" . $ivmax,
  );
  my @moduli = (3, 101, 65537, 4294967291);
  push @moduli, '18446744073709551557' if $bits == 64;

  my @cases;
  for my $a (@bases) {
    for my $e (@exponents) {
      for my $m (@moduli) {
        next unless _bi($a)->babs->bgcd(_bi($m))->is_one;
        push @cases, [$a, $e, $m];
      }
    }
  }
  for (1 .. 50) {
    my $a = _random_native();
    my $e = _bi(urandomb($bits-1))->binc->bneg->bstr;
    my $m = $moduli[urandomm(scalar @moduli)];
    next unless _bi($a)->bgcd(_bi($m))->is_one;
    push @cases, [$a, $e, $m];
  }

  for my $case (@cases) {
    my ($a, $e, $m) = @$case;
    my $results = _compare_exact('powmod', $case, \@fail);
    my $positive_e = _bi($e)->bneg;
    my $positive_power = _bi($a)->bmod($m)->bmodpow($positive_e, $m);
    for my $result (@$results) {
      my ($label, $values) = @$result;
      if (@$values != 1 || !defined $values->[0]) {
        push @fail, "$label powmod($a,$e,$m) did not return one value";
        next;
      }
      my $identity = _bi($values->[0])->bmul($positive_power)->bmod($m);
      push @fail, "$label powmod($a,$e,$m) fails inverse identity"
        unless $identity->is_one;
    }
  }
  _finish('negative-exponent powmod agrees and inverts', scalar(@cases), \@fail);
  done_testing();
};

sub _lucas_reference {
  my ($p, $q, $k) = @_;
  return ('0', '2') if $k == 0;
  return ('1', "$p") if $k == 1;

  my ($u0, $u1) = (Math::BigInt->bzero, Math::BigInt->bone);
  my ($v0, $v1) = (Math::BigInt->new(2), _bi($p));
  for (2 .. $k) {
    my $u2 = _bi($p)->bmul($u1)->bsub(_bi($q)->bmul($u0));
    my $v2 = _bi($p)->bmul($v1)->bsub(_bi($q)->bmul($v0));
    ($u0, $u1) = ($u1, $u2);
    ($v0, $v1) = ($v1, $v2);
  }
  return ($u1->bstr, $v1->bstr);
}

subtest 'Lucas functions' => sub {
  my @fail;
  my @parameters = (
    [1, -1], [-1, 1], [2, -3],
    [$ivmin, 1], [$ivmax, -1],
    [1, $ivmin], [$ivmin, $ivmax],
    [$uvmax, 1],
  );
  my @cases;
  push @cases, map {
    my ($p, $q) = @$_;
    map { [$p, $q, $_] } (0, 1, 2, 3, 5, 10)
  } @parameters;
  for (1 .. 80) {
    my $p = _random_iv();
    my $q = _random_iv();
    push @cases, [$p, $q, urandomm(11)];
  }

  for my $case (@cases) {
    my ($p, $q, $k) = @$case;
    my ($expect_u, $expect_v) = _lucas_reference($p, $q, $k);
    my $uvresults = _compare_exact('lucasuv', $case, \@fail);
    my $uresults = _compare_exact('lucasu', $case, \@fail);
    my $vresults = _compare_exact('lucasv', $case, \@fail);
    for my $result (@$uvresults) {
      my ($label, $values) = @$result;
      if (@$values != 2) {
        push @fail, "$label lucasuv($p,$q,$k) did not return two values";
        next;
      }
      push @fail, "$label lucasuv($p,$q,$k) disagrees with recurrence"
        if "$values->[0]" ne $expect_u || "$values->[1]" ne $expect_v;
    }
    for my $result (@$uresults) {
      my ($label, $values) = @$result;
      push @fail, "$label lucasu($p,$q,$k) disagrees with recurrence"
        if @$values != 1 || "$values->[0]" ne $expect_u;
    }
    for my $result (@$vresults) {
      my ($label, $values) = @$result;
      push @fail, "$label lucasv($p,$q,$k) disagrees with recurrence"
        if @$values != 1 || "$values->[0]" ne $expect_v;
    }
  }
  _finish('Lucas functions agree with recurrence', scalar(@cases), \@fail);
  done_testing();
};

sub _factorial_reference {
  my ($rising, $n, $k) = @_;
  my $result = Math::BigInt->bone;
  for my $i (0 .. $k-1) {
    my $term = _bi($n);
    $rising ? $term->badd($i) : $term->bsub($i);
    $result->bmul($term);
  }
  return $result->bstr;
}

subtest 'signed factorials' => sub {
  my @fail;
  my @nvalues = _uniq(@edge, $ivmin, $ivmax, $uvmax);
  my @cases;
  push @cases, map {
    my $n = $_;
    map { [$n, $_] } (0, 1, 2, 3, 5)
  } @nvalues;
  for (1 .. 80) {
    my $n = _random_iv();
    push @cases, [$n, urandomm(9)];
  }

  for my $fn (qw/falling_factorial rising_factorial/) {
    for my $case (@cases) {
      my $results = _compare_exact($fn, $case, \@fail);
      my $expected = _factorial_reference($fn eq 'rising_factorial', @$case);
      for my $result (@$results) {
        my ($label, $values) = @$result;
        if (@$values != 1) {
          push @fail, "$label $fn(" . _show_args(@$case) .
                      ") did not return one value";
          next;
        }
        push @fail, "$label $fn(" . _show_args(@$case) .
                    ") returned $values->[0], expected $expected"
          if "$values->[0]" ne $expected;
      }
    }
  }
  _finish('signed factorials agree with direct products',
          2 * scalar(@cases), \@fail);
  done_testing();
};

done_testing();
