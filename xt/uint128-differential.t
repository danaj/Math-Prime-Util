#!/usr/bin/env perl
use strict;
use warnings;

use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/
  prime_get_config prime_set_config
  factor factor_exp is_prob_prime is_bpsw_prime is_square
  is_semiprime is_square_free moebius liouville
  prime_omega prime_bigomega is_omega_prime is_almost_prime
  addmod submod mulmod muladdmod mulsubmod
/;
use Math::Prime::Util::PP ();

my $config = prime_get_config();
plan skip_all => 'the XS 128-bit path is not available'
  if !$config->{xs} || ($config->{xs_factor_bits} || 0) < 128;

my $old_gmp = $config->{gmp};
prime_set_config(gmp => 0);
END { prime_set_config(gmp => $old_gmp) if defined $old_gmp }

my %direct = map { $_ => Math::Prime::Util->can($_) } qw/
  factor factor_exp is_prob_prime is_bpsw_prime is_square
  is_semiprime is_square_free moebius liouville
  prime_omega prime_bigomega is_omega_prime is_almost_prime
  addmod submod mulmod muladdmod mulsubmod
/;

sub bi {
  return Math::BigInt->new("$_[0]");
}

sub canon {
  return defined($_[0]) ? "$_[0]" : '<undef>';
}

sub expected_mod {
  my ($op, $a, $b, $c, $m) = @_;
  my $mod = bi($m)->babs;
  return undef if $mod->is_zero;
  my $r = bi($a);
  if    ($op eq 'addmod')    { $r->badd($b) }
  elsif ($op eq 'submod')    { $r->bsub($b) }
  elsif ($op eq 'mulmod')    { $r->bmul($b) }
  elsif ($op eq 'muladdmod') { $r->bmul($b)->badd($c) }
  elsif ($op eq 'mulsubmod') { $r->bmul($b)->bsub($c) }
  else                       { die "unknown modular operation $op" }
  return $r->bmod($mod);
}

subtest '65-128-bit modular arithmetic' => sub {
  my @cases = (
    [141787151, 3743000622, 19, '4564087363318596642'],
    ['8789157053704531', '6762059301481433426361891290961', 37,
     '7005893202019645748120146985040'],
    ['175109911729618543589989257539043768012',
     '21887412602962542281538131483385626868',
     '263466159656861646486075450888763957942',
     '83494980727347746728226137271418851789'],
    ['-171821502870939196679518625154011220409',
     '182569474286058024586486590841354369890',
     '-329619958784558749434516006469339236320',
     '-49684876044205406960769234394385141897'],
    [0, '340282366920938463463374607431768211454', 1,
     '340282366920938463463374607431768211455'],
    [1, '340282366920938463463374607431768211454', -1,
     '340282366920938463463374607431768211455'],
    ['340282366920938463463374607431768211454', 1, 0,
     '340282366920938463463374607431768211455'],
    ['+00018446744073709551616', '-00018446744073709551617', 1,
     '+000340282366920938463463374607431768211455'],
  );

  my @bits = (65, 72, 80, 96, 111, 127, 128);
  for my $bi (0 .. $#bits) {
    my $top = Math::BigInt->bone->blsft($bits[$bi]);
    for my $sample (1 .. 8) {
      my $i = 8*$bi + $sample;
      my $m = $top->copy->bsub(2*$i+1);
      my $a = $top->copy->brsft(1)->badd(1_000_003*$i);
      my $b = $top->copy->brsft(2)->badd(1_000_033*$i);
      my $c = $top->copy->brsft(3)->badd(1_000_037*$i);
      $a->bneg if $i & 1;
      $b->bneg if $i & 2;
      $c->bneg if $i & 4;
      $m->bneg if $i % 11 == 0;
      push @cases, [ "$a", "$b", "$c", "$m" ];
    }
  }

  for my $op (qw/addmod submod mulmod muladdmod mulsubmod/) {
    my (@got, @exp);
    for my $case (@cases) {
      my ($a, $b, $c, $m) = @$case;
      my @args = $op =~ /^mul(?:add|sub)mod$/ ? ($a,$b,$c,$m) : ($a,$b,$m);
      push @got, canon($direct{$op}->(@args));
      push @exp, canon(expected_mod($op,$a,$b,$c,$m));
    }
    is_deeply(\@got, \@exp, "$op matches Math::BigInt");
  }
};

subtest 'primality around and above the 64-bit boundary' => sub {
  my @known = (
    ['18446744073709551617', 0],
    ['73786976294838218827', 1],
    ['170141183460469231731687303715884105727', 1],
    ['340282366920938463463374607431768211297', 1],
    ['340282366920938463426481119284349108225', 0],
    ['340282366920938463463374607431768211455', 0],
  );
  for my $case (@known) {
    my ($n, $expected) = @$case;
    is(!!$direct{is_prob_prime}->($n), !!$expected, "is_prob_prime($n)");
    is(!!$direct{is_bpsw_prime}->($n), !!$expected, "is_bpsw_prime($n)");
  }
  is(!!$direct{is_prob_prime}->('+00073786976294838218827'), 1,
     'primality parser accepts plus and leading zeroes');

  my @bits = (65, 79, 96, 111, 127, 128);
  for my $bits (@bits) {
    my $base = Math::BigInt->bone->blsft($bits-1);
    for my $i (1 .. 8) {
      my $n = $base->copy->badd(2*1_000_003*$i+1);
      my $pp_prob = Math::Prime::Util::PP::is_prob_prime("$n");
      my $pp_bpsw = Math::Prime::Util::PP::is_bpsw_prime("$n");
      is(!!$direct{is_prob_prime}->("$n"), !!$pp_prob,
         "is_prob_prime PP differential, $bits bits sample $i");
      is(!!$direct{is_bpsw_prime}->("$n"), !!$pp_bpsw,
         "is_bpsw_prime PP differential, $bits bits sample $i");
    }
  }
};

subtest 'perfect squares through 128 bits' => sub {
  my @roots = (
    '4294967296',
    '1099511627791',
    '1125899906842679',
    '9223372036854775783',
    '18446744073709551615',
  );
  for my $root (@roots) {
    my $square = bi($root)->bmul($root);
    is($direct{is_square}->("$square"), 1, "$root squared");
    is($direct{is_square}->($square->copy->bsub(1)->bstr), 0,
       "$root squared minus one");
    is($direct{is_square}->($square->copy->badd(1)->bstr), 0,
       "$root squared plus one");
    is($direct{is_square}->("-$square"), 0, "negative square of $root");
  }
};

sub flatten_factors {
  my ($spec) = @_;
  my @f;
  for my $pe (@$spec) {
    push @f, ("$pe->[0]") x $pe->[1];
  }
  return \@f;
}

sub make_product {
  my ($spec) = @_;
  my $n = Math::BigInt->bone;
  for my $pe (@$spec) {
    $n->bmul(bi($pe->[0])->bpow($pe->[1]));
  }
  return $n;
}

subtest 'easy constructed factorizations and predicates' => sub {
  my @specs = (
    [ [2,1], ['73786976294838218827',1] ],
    [ [2,8], [3,2], ['73786976294838218827',1] ],
    [ [5,4], ['302231454903657293688899',1] ],
    [ [2,1], [7,2], ['19807040628566084398386000047',1] ],
    [ [2,10], ['324518553658426726783156020588637',1] ],
    [ [2,1], ['170141183460469231731687303715884105727',1] ],
  );

  for my $spec (@specs) {
    my $n = make_product($spec);
    die "constructed input exceeds 128 bits" if $n > bi(2)->bpow(128)->bdec;
    my $flat = flatten_factors($spec);
    my @got = map { "$_" } $direct{factor}->("$n");
    is_deeply(\@got, $flat, "factor($n)");

    my @got_exp = map { [ "$_->[0]", "$_->[1]" ] }
                  $direct{factor_exp}->("+000$n");
    my @exp = map { [ "$_->[0]", "$_->[1]" ] } @$spec;
    is_deeply(\@got_exp, \@exp, "factor_exp leading-zero $n");

    my $total = 0;
    $total += $_->[1] for @$spec;
    my $distinct = scalar @$spec;
    my $square_free = !grep { $_->[1] > 1 } @$spec;
    my $mu = $square_free ? ($distinct & 1 ? -1 : 1) : 0;
    my $lambda = $total & 1 ? -1 : 1;

    is(scalar($direct{factor}->("$n")), $total, 'scalar factor count');
    is(scalar($direct{factor_exp}->("$n")), $distinct,
       'scalar factor_exp count');
    is($direct{is_semiprime}->("$n"), $total == 2 ? 1 : 0,
       'is_semiprime');
    is($direct{is_square_free}->("$n"), $square_free ? 1 : 0,
       'is_square_free');
    is($direct{moebius}->("-$n"), $mu, 'negative moebius');
    is($direct{liouville}->("$n"), $lambda, 'liouville');
    is($direct{prime_omega}->("$n"), $distinct, 'prime_omega');
    is($direct{prime_bigomega}->("$n"), $total, 'prime_bigomega');
    is($direct{is_omega_prime}->($distinct,"$n"), 1, 'is_omega_prime');
    is($direct{is_almost_prime}->($total,"$n"), 1, 'is_almost_prime');
  }
};

done_testing();
