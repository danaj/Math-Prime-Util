#!/usr/bin/env perl
use strict;
use warnings;

use Config;
use Math::BigInt lib => 'Calc';
use B ();
use Test::More;
use Math::Prime::Util qw/
  prime_get_config prime_set_config
  irand irand32 irand64 drand
  addint subint add1int sub1int mulint divint modint cdivint powint
  is_odd is_even is_square cmpint kronecker signint
/;

my @names = qw/
  irand irand32 irand64 drand
  addint subint add1int sub1int mulint divint modint cdivint powint
  is_odd is_even is_square cmpint kronecker signint
/;
my %direct = map { $_ => Math::Prime::Util->can($_) } @names;

sub x_irand     { irand }
sub x_irand32   { irand32 }
sub x_irand64   { irand64 }
sub x_drand     { drand }
sub x_addint    { addint($_[0], $_[1]) }
sub x_subint    { subint($_[0], $_[1]) }
sub x_add1int   { add1int($_[0]) }
sub x_sub1int   { sub1int($_[0]) }
sub x_mulint    { mulint($_[0], $_[1]) }
sub x_divint    { divint($_[0], $_[1]) }
sub x_modint    { modint($_[0], $_[1]) }
sub x_cdivint   { cdivint($_[0], $_[1]) }
sub x_powint    { powint($_[0], $_[1]) }
sub x_is_odd    { is_odd($_[0]) }
sub x_is_even   { is_even($_[0]) }
sub x_is_square { is_square($_[0]) }
sub x_cmpint    { cmpint($_[0], $_[1]) }
sub x_kronecker { kronecker($_[0], $_[1]) }
sub x_signint   { signint($_[0]) }

sub x_validate_any    { Math::Prime::Util::_validate_integer($_[0]) }
sub x_validate_nonneg { Math::Prime::Util::_validate_integer_nonneg($_[0]) }
sub x_validate_pos    { Math::Prime::Util::_validate_integer_positive($_[0]) }
sub x_validate_abs    { Math::Prime::Util::_validate_integer_abs($_[0]) }

my %xop = (
  irand     => \&x_irand,
  irand32   => \&x_irand32,
  drand     => \&x_drand,
  addint    => \&x_addint,
  subint    => \&x_subint,
  add1int   => \&x_add1int,
  sub1int   => \&x_sub1int,
  mulint    => \&x_mulint,
  divint    => \&x_divint,
  modint    => \&x_modint,
  cdivint   => \&x_cdivint,
  powint    => \&x_powint,
  is_odd    => \&x_is_odd,
  is_even   => \&x_is_even,
  is_square => \&x_is_square,
  cmpint    => \&x_cmpint,
  kronecker => \&x_kronecker,
  signint   => \&x_signint,
  _validate_integer          => \&x_validate_any,
  _validate_integer_nonneg   => \&x_validate_nonneg,
  _validate_integer_positive => \&x_validate_pos,
  _validate_integer_abs      => \&x_validate_abs,
);
$xop{irand64} = \&x_irand64 if $Config{uvsize} >= 8;

sub collect_ops {
  my ($op, $found, $seen) = @_;
  return if !defined($op) || !$$op || $seen->{$$op}++;
  my $name = eval { $op->name } || '';
  $found->{$name}++ if ref($op) eq 'B::OP';
  for my $method (qw/first other sibling/) {
    my $child = eval { $op->$method() };
    collect_ops($child, $found, $seen) if defined($child) && $$child;
  }
}

sub compiled_ops {
  my %found;
  for my $code (@_) {
    collect_ops(B::svref_2object($code)->ROOT, \%found, {});
  }
  return \%found;
}

my $compiled = compiled_ops(values %xop);
plan skip_all => 'custom ops are not available in this build'
  if !$compiled->{addint};

sub normalize_error {
  my ($error) = @_;
  $error =~ s/\s+at\s+.*\z//s;
  $error =~ s/\s+\z//;
  return $error;
}

sub outcome {
  my ($code, @args) = @_;
  my ($ok, $ret, $error);
  {
    local $@;
    $ok = eval { $ret = $code->(@args); 1 };
    $error = $@;
  }
  return [ error => normalize_error($error) ] if !$ok;
  return [ value => defined($ret) ? "$ret" : '<undef>', ref($ret) || '' ];
}

sub validation_outcome {
  my ($code, $input) = @_;
  my ($ok, $ret, $error);
  my $value = $input;
  {
    local $@;
    $ok = eval { $ret = $code->($value); 1 };
    $error = $@;
  }
  my $result = !$ok
             ? [ error => normalize_error($error) ]
             : [ value => defined($ret) ? "$ret" : '<undef>', ref($ret) || '' ];
  push @$result, defined($value) ? "$value" : '<undef>', ref($value) || '';
  return $result;
}

sub compare_calls {
  my ($name, $cases) = @_;
  my @got = map { outcome($xop{$name}, @$_) } @$cases;
  my @exp = map { outcome($direct{$name}, @$_) } @$cases;
  is_deeply(\@got, \@exp, "$name custom op matches direct XSUB");
}

subtest 'all intended call sites compiled as custom ops' => sub {
  my @missing = sort grep { !$compiled->{$_} } keys %xop;
  is_deeply(\@missing, [], 'all wrappers contain their custom op');
};

my ($ivmin, $ivmax, $uvmax) = $Config{ivsize} == 4
  ? qw/-2147483648 2147483647 4294967295/
  : qw/-9223372036854775808 9223372036854775807 18446744073709551615/;
my @values = (
  -17, -1, 0, 1, 2, 17,
  $ivmin, $ivmax, $uvmax,
  '18446744073709551616',
  '-18446744073709551616',
  '170141183460469231731687303715884105727',
  '-170141183460469231731687303715884105727',
  '340282366920938463463374607431768211456',
);

my @pairs;
for my $i (0 .. $#values) {
  push @pairs,
    [ $values[$i], $values[($i*5+3) % @values] ],
    [ $values[$i], 1 ],
    [ $values[$i], -1 ];
}

subtest 'integer arithmetic and comparison' => sub {
  compare_calls($_, \@pairs) for qw/addint subint mulint cmpint kronecker/;

  my @divpairs = grep { "$_->[1]" ne '0' } @pairs;
  push @divpairs, [1,0], [$uvmax,0], ['18446744073709551616',0];
  compare_calls($_, \@divpairs) for qw/divint modint cdivint/;

  my @powcases = (
    [2,0], [2,63], [-2,63],
    ['18446744073709551616',2],
    [0,'340282366920938463463374607431768211456'],
    [1,'340282366920938463463374607431768211456'],
    [-1,'340282366920938463463374607431768211457'],
    [2,-1],
  );
  compare_calls(powint => \@powcases);
};

subtest 'unary integer operations' => sub {
  my @unary = map { [$_] } (
    @values,
    '000184467440737095516160',
    '+000184467440737095516160',
    '-000184467440737095516160',
    '340282366920938463426481119284349108225',
    '340282366920938463463374607431768211455',
  );
  compare_calls($_, \@unary)
    for qw/add1int sub1int signint is_odd is_even is_square/;

  is(x_signint('000184467440737095516160'), 1,
     'leading-zero positive bigint has positive sign');
  is(x_signint('+000000000000000000000'), 0,
     'explicitly positive bigint zero has zero sign');
};

subtest 'validation custom ops' => sub {
  my %validators = (
    _validate_integer          => Math::Prime::Util->can('_validate_integer'),
    _validate_integer_nonneg   => Math::Prime::Util->can('_validate_integer_nonneg'),
    _validate_integer_positive => Math::Prime::Util->can('_validate_integer_positive'),
    _validate_integer_abs      => Math::Prime::Util->can('_validate_integer_abs'),
  );
  my @inputs = (
    '00042', '+42', '-42', '0',
    '184467440737095516160', '-184467440737095516160',
    '', '4.5', undef,
  );

  for my $name (sort keys %validators) {
    my (@got, @exp);
    for my $input (@inputs) {
      push @got, validation_outcome($xop{$name}, $input);
      push @exp, validation_outcome($validators{$name}, $input);
    }
    is_deeply(\@got, \@exp, "$name custom op matches direct XSUB");
  }
};

subtest 'magic and overloaded integer inputs' => sub {
  {
    package MPU::XopInteger;
    use overload '""' => sub { ${$_[0]} }, fallback => 0;
    sub new { my ($class, $value) = @_; bless \$value, $class }
  }

  my $object = MPU::XopInteger->new('000184467440737095516160');
  compare_calls(signint => [[$object]]);
  compare_calls(addint  => [[$object, 17]]);
  compare_calls(cmpint  => [[$object, '184467440737095516160']]);

  '184467440737095516160' =~ /(\d+)/;
  compare_calls(signint => [[$1]]);
  compare_calls(addint  => [[$1, 17]]);
};

subtest 'random generators match direct XSUB streams' => sub {
  for my $name (grep { exists $xop{$_} } qw/irand irand32 irand64 drand/) {
    Math::Prime::Util::srand(0x5219_74A3);
    my @got = map { $xop{$name}->() } 1 .. 32;
    Math::Prime::Util::srand(0x5219_74A3);
    my @exp = map { $direct{$name}->() } 1 .. 32;
    is_deeply(\@got, \@exp, "$name stream");
  }
};

subtest 'a backend croak does not corrupt the following custom op' => sub {
  plan skip_all => 'GMP 0.47+ required to exercise a custom-op backend call'
    if $ENV{MPU_NO_GMP} ||
       !eval { require Math::Prime::Util::GMP;
               Math::Prime::Util::GMP->VERSION(0.47); 1 };
  my $old_gmp = prime_get_config()->{gmp};
  prime_set_config(gmp => 1);
  my $n = '340282366920938463463374607431768211456000001';
  my $bad;
  {
    no warnings qw/once redefine/;
    local *Math::Prime::Util::GMP::is_square = sub {
      die "forced xop backend failure\n";
    };
    $bad = outcome(\&x_is_square, $n);
  }
  is_deeply($bad, [error => 'forced xop backend failure'],
            'nested GMP backend croak is propagated');
  is("" . x_addint('184467440737095516160', 17),
     '184467440737095516177', 'subsequent bigint operation succeeds');
  is(x_is_square($n), $direct{is_square}->($n),
     'subsequent dispatched operation succeeds');
  prime_set_config(gmp => $old_gmp);
};

done_testing();
