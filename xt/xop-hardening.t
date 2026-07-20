#!/usr/bin/env perl
use strict;
use warnings;

use B ();
use overload ();
use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/
  prime_get_config prime_set_config
  addint subint add1int sub1int mulint divint modint cdivint powint
  is_odd is_even is_square cmpint kronecker signint
/;

my @names = qw/
  addint subint add1int sub1int mulint divint modint cdivint powint
  is_odd is_even is_square cmpint kronecker signint
/;
# These can be raised for longer stress runs without editing the test.
my $NRECOVERY = $ENV{MPU_XOP_RECOVERY_ROUNDS} || 1_000;
my $NTREES     = $ENV{MPU_XOP_TREES}          || 5_000;
my $TREE_DEPTH = $ENV{MPU_XOP_TREE_DEPTH}     || 4;

my %direct = map { $_ => Math::Prime::Util->can($_) } @names;
$direct{_validate_integer} =
  Math::Prime::Util->can('_validate_integer');
$direct{_validate_integer_nonneg} =
  Math::Prime::Util->can('_validate_integer_nonneg');
$direct{_validate_integer_positive} =
  Math::Prime::Util->can('_validate_integer_positive');
$direct{_validate_integer_abs} =
  Math::Prime::Util->can('_validate_integer_abs');

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

# These retain the original prototypes while bypassing the call checker.
sub d_addint($$)    { $direct{addint}->($_[0], $_[1]) }
sub d_subint($$)    { $direct{subint}->($_[0], $_[1]) }
sub d_add1int($)    { $direct{add1int}->($_[0]) }
sub d_sub1int($)    { $direct{sub1int}->($_[0]) }
sub d_mulint($$)    { $direct{mulint}->($_[0], $_[1]) }
sub d_divint($$)    { $direct{divint}->($_[0], $_[1]) }
sub d_modint($$)    { $direct{modint}->($_[0], $_[1]) }
sub d_cdivint($$)   { $direct{cdivint}->($_[0], $_[1]) }
sub d_powint($$)    { $direct{powint}->($_[0], $_[1]) }
sub d_is_odd($)     { $direct{is_odd}->($_[0]) }
sub d_is_even($)    { $direct{is_even}->($_[0]) }
sub d_is_square($)  { $direct{is_square}->($_[0]) }
sub d_cmpint($$)    { $direct{cmpint}->($_[0], $_[1]) }
sub d_kronecker($$) { $direct{kronecker}->($_[0], $_[1]) }
sub d_signint($)    { $direct{signint}->($_[0]) }

sub _collect_ops {
  my ($op, $found, $seen) = @_;
  return if !defined($op) || !$$op || $seen->{$$op}++;

  my $name = eval { $op->name } || '';
  $found->{$name}++ if ref($op) eq 'B::OP';

  for my $method (qw/first other sibling/) {
    my $child = eval { $op->$method() };
    _collect_ops($child, $found, $seen) if defined($child) && $$child;
  }
}

sub _compiled_ops {
  my %found;
  for my $code (@_) {
    _collect_ops(B::svref_2object($code)->ROOT, \%found, {});
  }
  return \%found;
}

my $compiled = _compiled_ops(
  \&x_addint, \&x_subint, \&x_add1int, \&x_sub1int, \&x_mulint,
  \&x_divint, \&x_modint, \&x_cdivint, \&x_powint,
  \&x_is_odd, \&x_is_even, \&x_is_square, \&x_cmpint,
  \&x_kronecker, \&x_signint,
  \&x_validate_any, \&x_validate_nonneg, \&x_validate_pos, \&x_validate_abs,
);
plan skip_all => 'custom ops are not available in this build'
  if !$compiled->{addint};

sub _normalize_error {
  my ($error) = @_;
  $error =~ s/\s+at\s+.*\z//s;
  $error =~ s/\s+\z//;
  return $error;
}

sub _outcome {
  my ($code, @args) = @_;
  my ($ok, $ret, $error);
  {
    local $@;
    $ok = eval { $ret = $code->(@args); 1 };
    $error = $@;
  }
  return [error => _normalize_error($error)] if !$ok;
  return [value => defined($ret) ? "$ret" : '<undef>', ref($ret) || ''];
}

sub _string_list {
  return [map { defined($_) ? "$_" : '<undef>' } @_];
}

###############################################################################
# Context, argument evaluation order, and exactly-once evaluation.

my @trace;

sub _context_arg {
  my ($label, $value) = @_;
  my $wa = wantarray;
  push @trace, $label . ':' .
               (!defined($wa) ? 'V' : $wa ? 'L' : 'S');
  return $wa ? ($value, 'unexpected-list-value') : $value;
}

sub _context_xop {
  @trace = ();

  my $scalar = addint(_context_arg('sa',2), _context_arg('sb',3));
  my @list = (11, addint(_context_arg('la',4),
                          _context_arg('lb',5)), 13);
  addint(_context_arg('va',6), _context_arg('vb',7));
  my $after_void = addint(40,2);
  my $boolean = addint(_context_arg('ba',0),
                        _context_arg('bb',1)) ? 1 : 0;
  my @mapped = map {
    addint(_context_arg("m$_",$_), 10)
  } 1 .. 3;
  my @grepped = grep {
    is_odd(_context_arg("g$_",$_))
  } 1 .. 5;
  my @sorted = sort {
    cmpint(_context_arg("s$a",$a), _context_arg("s$b",$b))
  } (5,1,4,2,3);

  return {
    scalar     => "$scalar",
    list       => _string_list(@list),
    after_void => "$after_void",
    boolean    => $boolean,
    mapped     => _string_list(@mapped),
    grepped    => _string_list(@grepped),
    sorted     => _string_list(@sorted),
    trace      => [@trace],
  };
}

sub _context_direct {
  @trace = ();

  my $scalar = d_addint(_context_arg('sa',2), _context_arg('sb',3));
  my @list = (11, d_addint(_context_arg('la',4),
                            _context_arg('lb',5)), 13);
  d_addint(_context_arg('va',6), _context_arg('vb',7));
  my $after_void = d_addint(40,2);
  my $boolean = d_addint(_context_arg('ba',0),
                          _context_arg('bb',1)) ? 1 : 0;
  my @mapped = map {
    d_addint(_context_arg("m$_",$_), 10)
  } 1 .. 3;
  my @grepped = grep {
    d_is_odd(_context_arg("g$_",$_))
  } 1 .. 5;
  my @sorted = sort {
    d_cmpint(_context_arg("s$a",$a), _context_arg("s$b",$b))
  } (5,1,4,2,3);

  return {
    scalar     => "$scalar",
    list       => _string_list(@list),
    after_void => "$after_void",
    boolean    => $boolean,
    mapped     => _string_list(@mapped),
    grepped    => _string_list(@grepped),
    sorted     => _string_list(@sorted),
    trace      => [@trace],
  };
}

subtest 'context and evaluation semantics' => sub {
  my $got = _context_xop();
  my $exp = _context_direct();
  is_deeply($got, $exp, 'custom ops match prototyped XSUB proxy calls');
  ok(!(grep { $_ !~ /:S\z/ } @{$got->{trace}}),
     'every custom-op argument was evaluated in scalar context');
  is_deeply($got->{list}, [qw/11 9 13/],
            'list context receives exactly one result');
  is($got->{after_void}, 42, 'void call leaves the stack balanced');

  my @input = (
    Math::BigInt->new('184467440737095516160'),
    Math::BigInt->new('340282366920938463463374607431768211456'),
  );
  my @input_before = map { "$_" } @input;
  my @held = (
    addint(40,2),
    addint($input[0],17),
    addint($input[0],$input[0]),
    addint($input[1],$input[0]),
  );
  my @held_before = map { "$_" } @held;

  for my $i (1 .. 1_000) {
    my $discard = addint(
      mulint(add1int($i),3),
      subint(cdivint($i+7,3),modint($i,5)),
    );
  }
  is_deeply([map { "$_" } @held], \@held_before,
            'retained native and bigint results survive stack churn');
  is_deeply([map { "$_" } @input], \@input_before,
            'aliased bigint operands are not modified');
};

###############################################################################
# Reentrant magic and overload callbacks.

{
  package MPU::XopHardening::Tie;
  our $fetches = 0;
  our $croak = 0;

  sub TIESCALAR {
    my ($class, $value) = @_;
    return bless \$value, $class;
  }
  sub FETCH {
    $fetches++;
    my $value = Math::Prime::Util::addint(${$_[0]}, 0);
    die "forced tied FETCH failure\n" if $croak;
    return $value;
  }
}

{
  package MPU::XopHardening::Overload;
  use overload '""' => '_stringify', fallback => 0;
  our $calls = 0;
  our $croak = 0;

  sub new {
    my ($class, $value) = @_;
    return bless \$value, $class;
  }
  sub _stringify {
    $calls++;
    my $value = Math::Prime::Util::addint(${$_[0]}, 0);
    die "forced overloaded stringify failure\n" if $croak;
    return "$value";
  }
}

{
  package MPU::XopHardening::RecursiveOverload;
  use overload '""' => '_stringify', fallback => 0;
  our $calls = 0;

  sub new {
    my ($class, $value, $depth, $fail_at) = @_;
    return bless [$value, $depth, $fail_at], $class;
  }
  sub _stringify {
    my ($value, $depth, $fail_at) = @{$_[0]};
    $calls++;
    die "forced recursive stringify failure\n"
      if defined($fail_at) && $depth == $fail_at;
    return "$value" if $depth <= 0;
    my $child = __PACKAGE__->new($value, $depth-1, $fail_at);
    return "" . Math::Prime::Util::addint($child, 1);
  }
}

our ($x_tied, $d_tied);
sub _x_tied_add { addint($x_tied,5) }
sub _d_tied_add { d_addint($d_tied,5) }
sub _x_object_add { addint($_[0],7) }
sub _d_object_add { d_addint($_[0],7) }
sub _x_regex_add { addint($1,9) }
sub _d_regex_add { d_addint($1,9) }

my $DISPATCH_A = '340282366920938463463374607431768211457';
my $DISPATCH_B = '1393796574908163946345982392040522594123776';
my $KRON128_A  = '340282366920938463463374607431768211455';
my $KRON128_B  = '170141183460469231731687303715884105727';
our ($recursive_pp_remaining, $recursive_pp_calls, $recursive_pp_fail);

sub _recursive_pp_kronecker {
  $recursive_pp_calls++;
  if ($recursive_pp_remaining > 0) {
    $recursive_pp_remaining--;
    return kronecker($DISPATCH_A,97);
  }
  die "forced recursive PP dispatch failure\n" if $recursive_pp_fail;
  return kronecker(5,7);
}

subtest 'reentrant magic and overload callbacks' => sub {
  tie $x_tied, 'MPU::XopHardening::Tie', '184467440737095516160';
  $MPU::XopHardening::Tie::fetches = 0;
  my $xgot = _outcome(\&_x_tied_add);
  my $xfetches = $MPU::XopHardening::Tie::fetches;

  tie $d_tied, 'MPU::XopHardening::Tie', '184467440737095516160';
  $MPU::XopHardening::Tie::fetches = 0;
  my $dgot = _outcome(\&_d_tied_add);
  my $dfetches = $MPU::XopHardening::Tie::fetches;

  is_deeply($xgot, $dgot, 'tied FETCH invoking an xop matches direct XSUB');
  is($xfetches, $dfetches, 'tied FETCH count matches direct XSUB');
  ok($xfetches > 0, 'tied FETCH was exercised');

  my $xobj = MPU::XopHardening::Overload->new(
    '184467440737095516160'
  );
  $MPU::XopHardening::Overload::calls = 0;
  my $xobjgot = _outcome(\&_x_object_add, $xobj);
  my $xcalls = $MPU::XopHardening::Overload::calls;

  my $dobj = MPU::XopHardening::Overload->new(
    '184467440737095516160'
  );
  $MPU::XopHardening::Overload::calls = 0;
  my $dobjgot = _outcome(\&_d_object_add, $dobj);
  my $dcalls = $MPU::XopHardening::Overload::calls;

  is_deeply($xobjgot, $dobjgot,
            'overloaded stringification invoking an xop matches direct XSUB');
  is($xcalls, $dcalls, 'overloaded stringification count matches direct XSUB');
  ok($xcalls > 0, 'overloaded stringification was exercised');

  my $recursive_ops = _compiled_ops(
    \&MPU::XopHardening::RecursiveOverload::_stringify
  );
  ok($recursive_ops->{addint},
     'recursive stringification callback contains a nested custom op');
  $MPU::XopHardening::RecursiveOverload::calls = 0;
  my $xrecursive = _outcome(
    \&_x_object_add,
    MPU::XopHardening::RecursiveOverload->new(17, 24, undef)
  );
  my $xrecursive_calls = $MPU::XopHardening::RecursiveOverload::calls;
  $MPU::XopHardening::RecursiveOverload::calls = 0;
  my $drecursive = _outcome(
    \&_d_object_add,
    MPU::XopHardening::RecursiveOverload->new(17, 24, undef)
  );
  my $drecursive_calls = $MPU::XopHardening::RecursiveOverload::calls;
  is_deeply($xrecursive, $drecursive,
            'deep recursive custom-op callbacks match direct XSUB');
  is($xrecursive_calls, $drecursive_calls,
     'deep recursive callback count matches direct XSUB');
  ok($xrecursive_calls >= 25, 'deep recursive callback path was exercised');

  '184467440737095516160' =~ /(\d+)/;
  my $xregex = _outcome(\&_x_regex_add);
  '184467440737095516160' =~ /(\d+)/;
  my $dregex = _outcome(\&_d_regex_add);
  is_deeply($xregex, $dregex,
            'regex magic matches direct XSUB');

  $MPU::XopHardening::Overload::croak = 1;
  my $bad = _outcome(\&_x_object_add,
                     MPU::XopHardening::Overload->new(17));
  $MPU::XopHardening::Overload::croak = 0;
  is_deeply($bad, [error => 'forced overloaded stringify failure'],
            'reentrant overload croak is propagated');
  is(x_addint(x_mulint(6,7), x_sub1int(1)), 42,
     'nested custom ops work after overloaded callback croak');

  $MPU::XopHardening::Tie::croak = 1;
  $bad = _outcome(\&_x_tied_add);
  $MPU::XopHardening::Tie::croak = 0;
  is_deeply($bad, [error => 'forced tied FETCH failure'],
            'reentrant tied FETCH croak is propagated');
  is(x_addint(x_mulint(6,7), x_sub1int(1)), 42,
     'nested custom ops work after tied callback croak');

  $bad = _outcome(
    \&_x_object_add,
    MPU::XopHardening::RecursiveOverload->new(17, 16, 5)
  );
  is_deeply($bad, [error => 'forced recursive stringify failure'],
            'deep recursive callback croak is propagated');
  is(x_addint(x_mulint(6,7), x_sub1int(1)), 42,
     'nested custom ops work after deep recursive croak');
  untie $x_tied;
  untie $d_tied;
};

###############################################################################
# Real and controlled XS -> Perl -> XS callback chains.
#
# With GMP disabled, bigint kronecker dispatches to PP.  PP validates through
# custom ops and repeatedly calls the regular XS right-shift implementation.
# Canonical bigint results call a Perl constructor, which we instrument to
# re-enter another custom op.

subtest 'Perl dispatch and XS re-entry' => sub {
  require Math::Prime::Util::PP;
  my $old_config = prime_get_config();
  my $old_gmp = $old_config->{gmp};
  my $old_bigint = $old_config->{bigintclass};
  prime_set_config(gmp => 0, bigint => 'Math::BigInt');

  my $pp_ops = _compiled_ops(\&Math::Prime::Util::PP::kronecker);
  my @missing = grep { !$pp_ops->{$_} }
                qw/_validate_integer kronecker/;
  is_deeply(\@missing, [],
            'PP kronecker validators and final reduction re-enter custom ops');

  my $original_pp_kronecker = \&Math::Prime::Util::PP::kronecker;
  my $original_rshift = \&Math::Prime::Util::PP::Mrshiftint;
  my $original_tobigint = \&Math::Prime::Util::_to_bigint;
  my $original_new = Math::BigInt->can('new');
  my ($pp_calls, $shift_calls, $tobigint_calls) = (0,0,0);
  my ($new_calls, $new_xop_calls) = (0,0);

  SKIP: {
    skip '128-bit native factoring support is not available', 3
      if ($old_config->{xs_factor_bits} || 0) < 128;
    my $native_pp_calls = 0;
    {
      no warnings qw/once redefine/;
      local *Math::Prime::Util::PP::kronecker = sub {
        $native_pp_calls++;
        goto &$original_pp_kronecker;
      };
      is(x_kronecker($KRON128_A,$KRON128_B), 1,
         '128-bit kronecker custom op returns correctly');
      is($direct{kronecker}->($KRON128_A,$KRON128_B), 1,
         '128-bit kronecker XSUB returns correctly');
    }
    is($native_pp_calls, 0, '128-bit kronecker does not dispatch to PP');
  }

  my $new_wrapper = sub {
    $new_calls++;
    my $sign = signint(defined($_[1]) ? $_[1] : 0);
    $new_xop_calls++ if $sign >= 0;
    goto &$original_new;
  };
  my $new_ops = _compiled_ops($new_wrapper);
  ok($new_ops->{signint},
     'instrumented bigint constructor contains a nested custom op');

  {
    no warnings qw/once redefine/;
    local *Math::Prime::Util::PP::kronecker = sub {
      $pp_calls++;
      goto &$original_pp_kronecker;
    };
    local *Math::Prime::Util::PP::Mrshiftint = sub {
      $shift_calls++;
      goto &$original_rshift;
    };
    local *Math::Prime::Util::_to_bigint = sub {
      $tobigint_calls++;
      goto &$original_tobigint;
    };
    local *Math::BigInt::new = $new_wrapper;

    is(x_kronecker($DISPATCH_A,$DISPATCH_B), 1,
       'bigint kronecker survives repeated PP and XS crossings');
  }
  is($pp_calls, 1, 'outer custom op dispatched to PP once');
  is($shift_calls, 140, 'PP called the XS right-shift path 140 times');
  is($tobigint_calls, 2,
     'custom validators called Perl bigint conversion twice');
  ok($new_calls >= 140, 'callback chain repeatedly constructed bigints');
  is($new_xop_calls, $new_calls,
     'every instrumented bigint constructor re-entered a custom op');

  my $recursive_ops = _compiled_ops(\&_recursive_pp_kronecker);
  ok($recursive_ops->{kronecker},
     'controlled PP callback contains a recursive custom op');
  {
    no warnings qw/once redefine/;
    local *Math::Prime::Util::PP::kronecker =
      \&_recursive_pp_kronecker;

    ($recursive_pp_remaining, $recursive_pp_calls, $recursive_pp_fail) =
      (31, 0, 0);
    is_deeply(_outcome(\&x_kronecker, $DISPATCH_A,97),
              [value => '-1', ''],
              '32-level recursive PP dispatch returns correctly');
    is($recursive_pp_calls, 32,
       'recursive PP dispatch reached the requested depth');

    ($recursive_pp_remaining, $recursive_pp_calls, $recursive_pp_fail) =
      (15, 0, 1);
    is_deeply(_outcome(\&x_kronecker, $DISPATCH_A,97),
              [error => 'forced recursive PP dispatch failure'],
              'deep recursive PP dispatch propagates croak');
  }
  is(x_addint(x_mulint(6,7),x_sub1int(1)), 42,
     'custom ops recover after recursive PP dispatch croak');

  for my $class (qw/Math::GMPz Math::GMP/) {
    SKIP: {
      (my $module = "$class.pm") =~ s!::!/!g;
      skip "$class is not available", 4
        if !eval { require $module; 1 };

      my $input = $class->new('184467440737095516160');
      my $method = overload::Method($input, '+');
      my $method_cv = $method ? B::svref_2object($method) : undef;
      ok($method_cv && eval { $method_cv->XSUB },
         "$class addition overload is an XSUB");

      prime_set_config(bigint => $class);
      my $result = x_addint($input,17);
      is("$result", '184467440737095516177',
         "addint returns through $class overloaded XS");
      is(ref($result), $class,
         "addint preserves configured $class result");
      is_deeply(_outcome(\&x_modint, $input, 0),
                [error => 'modint: divide by zero'],
                "modint traps zero before calling $class overload");
    }
  }

  SKIP: {
    skip 'GMP backend 0.47 or newer is not available', 4
      if !$old_gmp || $old_gmp < 47;
    prime_set_config(gmp => 1);
    my $gmp_square = Math::Prime::Util::GMP->can('is_square');
    my $gmp_kronecker = Math::Prime::Util::GMP->can('kronecker');
    skip 'required GMP functions are not available', 4
      if !$gmp_square || !$gmp_kronecker;

    my ($square_calls, $kronecker_calls) = (0,0);
    {
      no warnings qw/once redefine/;
      local *Math::Prime::Util::GMP::is_square = sub {
        $square_calls++;
        goto &$gmp_square;
      };
      local *Math::Prime::Util::GMP::kronecker = sub {
        $kronecker_calls++;
        goto &$gmp_kronecker;
      };

      is(x_is_square(
           '765413284212226299051111674934086564882382225721'
         ), 1, 'is_square custom op returns through GMP XS');
      is($square_calls, 1, 'is_square dispatched to GMP once');
      is(x_kronecker($DISPATCH_A,97),
         $gmp_kronecker->($DISPATCH_A,97),
         'kronecker custom op returns through GMP XS');
      is($kronecker_calls, 1, 'kronecker dispatched to GMP once');
    }
  }

  prime_set_config(bigint => $old_bigint) if $old_bigint;
  prime_set_config(gmp => $old_gmp ? 1 : 0);
};

###############################################################################
# Repeated error paths followed by valid nested custom operations.

sub _argument_croak {
  addint(1,2);
  die "forced argument failure\n";
}
sub _x_argument_failure { addint(_argument_croak(), 7) }
sub _d_argument_failure { d_addint(_argument_croak(), 7) }

my @argument_trace;
sub _traced_argument {
  my ($label, $value, $fail) = @_;
  push @argument_trace, $label;
  addint(1,2);
  die "forced $label argument failure\n" if $fail;
  return $value;
}
sub _x_first_argument_failure {
  addint(_traced_argument('first',1,1),
         _traced_argument('second',2,0))
}
sub _d_first_argument_failure {
  d_addint(_traced_argument('first',1,1),
           _traced_argument('second',2,0))
}
sub _x_second_argument_failure {
  addint(_traced_argument('first',1,0),
         _traced_argument('second',2,1))
}
sub _d_second_argument_failure {
  d_addint(_traced_argument('first',1,0),
           _traced_argument('second',2,1))
}

subtest 'exception and stack recovery' => sub {
  my @cases = (
    ['divint by zero',  \&x_divint,  $direct{divint},  [17,0]],
    ['modint by zero',  \&x_modint,  $direct{modint},  [17,0]],
    ['cdivint by zero', \&x_cdivint, $direct{cdivint}, [17,0]],
    ['negative power',  \&x_powint,  $direct{powint},  [2,-1]],
    ['invalid addint',  \&x_addint,  $direct{addint},  ['bad',1]],
    ['invalid square',  \&x_is_square, $direct{is_square}, ['bad']],
    ['invalid symbol',  \&x_kronecker, $direct{kronecker}, ['bad',3]],
    ['invalid integer', \&x_validate_any,
                         $direct{_validate_integer}, ['4.5']],
    ['negative nonneg', \&x_validate_nonneg,
                         $direct{_validate_integer_nonneg}, [-1]],
    ['zero positive',   \&x_validate_pos,
                         $direct{_validate_integer_positive}, [0]],
    ['invalid abs',     \&x_validate_abs,
                         $direct{_validate_integer_abs}, ['4.5']],
  );

  for my $case (@cases) {
    my ($name, $xcode, $dcode, $args) = @$case;
    is_deeply(_outcome($xcode, @$args), _outcome($dcode, @$args),
              "$name matches direct XSUB");
  }

  is_deeply(_outcome(\&_x_argument_failure),
            _outcome(\&_d_argument_failure),
            'argument-evaluation croak matches prototyped XSUB proxy');

  for my $position (qw/first second/) {
    my $xcode = $position eq 'first'
              ? \&_x_first_argument_failure
              : \&_x_second_argument_failure;
    my $dcode = $position eq 'first'
              ? \&_d_first_argument_failure
              : \&_d_second_argument_failure;
    @argument_trace = ();
    my $xoutcome = _outcome($xcode);
    my @xtrace = @argument_trace;
    @argument_trace = ();
    my $doutcome = _outcome($dcode);
    my @dtrace = @argument_trace;
    is_deeply($xoutcome, $doutcome,
              "$position-argument croak matches prototyped XSUB proxy");
    is_deeply(\@xtrace, \@dtrace,
              "$position-argument croak preserves evaluation order");
  }

  my @fail;
  for my $round (1 .. $NRECOVERY) {
    my $case = $cases[$round % @cases];
    my $bad = _outcome($case->[1], @{$case->[3]});
    push @fail, "round $round did not fail" if $bad->[0] ne 'error';

    my $got = x_addint(
      x_mulint(x_add1int($round),3),
      x_subint(x_cdivint($round+7,3),x_modint($round,5)),
    );
    my $exp = $direct{addint}->(
      $direct{mulint}->($direct{add1int}->($round),3),
      $direct{subint}->($direct{cdivint}->($round+7,3),
                        $direct{modint}->($round,5)),
    );
    push @fail, "round $round stack/result mismatch: $got != $exp"
      if "$got" ne "$exp";
  }
  is_deeply(\@fail, [],
            "$NRECOVERY alternating failure/recovery rounds");

  require Math::Prime::Util::PP;
  @fail = ();
  my $old_gmp = prime_get_config()->{gmp};
  prime_set_config(gmp => 0);
  {
    no warnings qw/once redefine/;
    local *Math::Prime::Util::PP::kronecker = sub {
      die "forced xop backend failure\n";
    };
    for my $round (1 .. 20) {
      my $bad = _outcome(
        \&x_kronecker,
        '340282366920938463463374607431768211457', 97
      );
      push @fail, "backend round $round did not propagate croak"
        if $bad->[0] ne 'error' ||
           $bad->[1] ne 'forced xop backend failure';
      my $sentinel = x_addint(x_mulint(6,7),x_sub1int(1));
      push @fail, "backend round $round corrupted following op"
        if "$sentinel" ne '42';
    }
  }
  prime_set_config(gmp => $old_gmp ? 1 : 0);
  is_deeply(\@fail, [], 'backend croaks preserve subsequent custom ops');
};

###############################################################################
# Deterministic randomized nested-expression differential.

sub _eval_xop_tree {
  my ($t) = @_;
  my $op = $t->[0];
  return $t->[1] if $op eq 'leaf';
  return addint(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'add';
  return subint(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'sub';
  return mulint(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'mul';
  return divint(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'div';
  return modint(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'mod';
  return cdivint(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'cdiv';
  return powint(_eval_xop_tree($t->[1]), $t->[2])
    if $op eq 'pow';
  return add1int(_eval_xop_tree($t->[1])) if $op eq 'add1';
  return sub1int(_eval_xop_tree($t->[1])) if $op eq 'sub1';
  return signint(_eval_xop_tree($t->[1])) if $op eq 'sign';
  return is_odd(_eval_xop_tree($t->[1])) if $op eq 'odd';
  return is_even(_eval_xop_tree($t->[1])) if $op eq 'even';
  return is_square(_eval_xop_tree($t->[1])) if $op eq 'square';
  return cmpint(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'cmp';
  return kronecker(_eval_xop_tree($t->[1]), _eval_xop_tree($t->[2]))
    if $op eq 'kron';
  die "unknown tree op $op";
}

sub _eval_direct_tree {
  my ($t) = @_;
  my $op = $t->[0];
  return $t->[1] if $op eq 'leaf';
  return $direct{addint}->(_eval_direct_tree($t->[1]),
                           _eval_direct_tree($t->[2])) if $op eq 'add';
  return $direct{subint}->(_eval_direct_tree($t->[1]),
                           _eval_direct_tree($t->[2])) if $op eq 'sub';
  return $direct{mulint}->(_eval_direct_tree($t->[1]),
                           _eval_direct_tree($t->[2])) if $op eq 'mul';
  return $direct{divint}->(_eval_direct_tree($t->[1]),
                           _eval_direct_tree($t->[2])) if $op eq 'div';
  return $direct{modint}->(_eval_direct_tree($t->[1]),
                           _eval_direct_tree($t->[2])) if $op eq 'mod';
  return $direct{cdivint}->(_eval_direct_tree($t->[1]),
                            _eval_direct_tree($t->[2])) if $op eq 'cdiv';
  return $direct{powint}->(_eval_direct_tree($t->[1]), $t->[2])
    if $op eq 'pow';
  return $direct{add1int}->(_eval_direct_tree($t->[1])) if $op eq 'add1';
  return $direct{sub1int}->(_eval_direct_tree($t->[1])) if $op eq 'sub1';
  return $direct{signint}->(_eval_direct_tree($t->[1])) if $op eq 'sign';
  return $direct{is_odd}->(_eval_direct_tree($t->[1])) if $op eq 'odd';
  return $direct{is_even}->(_eval_direct_tree($t->[1])) if $op eq 'even';
  return $direct{is_square}->(_eval_direct_tree($t->[1]))
    if $op eq 'square';
  return $direct{cmpint}->(_eval_direct_tree($t->[1]),
                           _eval_direct_tree($t->[2])) if $op eq 'cmp';
  return $direct{kronecker}->(_eval_direct_tree($t->[1]),
                              _eval_direct_tree($t->[2])) if $op eq 'kron';
  die "unknown tree op $op";
}

my $rng_state = 1;
sub _rand32 {
  $rng_state = ($rng_state * 48271) % 2147483647;
  return $rng_state;
}

my @leaves = (
  -17, -2, -1, 0, 1, 2, 3, 17,
  '00000000000000000042', '+000184467440737095516160',
  '-000184467440737095516160',
  '2147483647', '2147483648', '4294967295', '4294967296',
  '9223372036854775807', '9223372036854775808',
  '18446744073709551615', '18446744073709551616',
  '-18446744073709551616',
  '170141183460469231731687303715884105727',
);
my @binary_ops = qw/add sub mul div mod cdiv cmp kron/;
my @unary_ops = qw/add1 sub1 sign odd even square/;

sub _make_tree {
  my ($depth) = @_;
  return ['leaf', $leaves[_rand32() % @leaves]]
    if $depth <= 0 || _rand32() % 5 == 0;

  my $kind = _rand32() % 10;
  if ($kind < 5) {
    my $op = $binary_ops[_rand32() % @binary_ops];
    return [$op, _make_tree($depth-1), _make_tree($depth-1)];
  }
  if ($kind < 9) {
    my $op = $unary_ops[_rand32() % @unary_ops];
    return [$op, _make_tree($depth-1)];
  }
  return ['pow', _make_tree($depth-1), _rand32() % 5];
}

sub _show_tree {
  my ($t) = @_;
  return "$t->[1]" if $t->[0] eq 'leaf';
  return $t->[0] . '(' . _show_tree($t->[1]) . ')'
    if @$t == 2;
  return $t->[0] . '(' . _show_tree($t->[1]) . ',' . $t->[2] . ')'
    if $t->[0] eq 'pow';
  return $t->[0] . '(' . _show_tree($t->[1]) . ',' .
         _show_tree($t->[2]) . ')';
}

subtest 'randomized nested-expression differential' => sub {
  my $tree_ops = _compiled_ops(\&_eval_xop_tree);
  my @missing = grep { !$tree_ops->{$_} } @names;
  is_deeply(\@missing, [], 'random evaluator contains every custom arithmetic op');

  $rng_state = 0x5219743;
  my @trees = map { _make_tree($TREE_DEPTH) } 1 .. $NTREES;
  my $old_gmp = prime_get_config()->{gmp};
  my @modes = ([$old_gmp ? 'configured backend' : 'PP backend',
                $old_gmp ? 1 : 0]);
  push @modes, ['forced PP backend', 0] if $old_gmp;

  for my $mode (@modes) {
    my ($label, $gmp) = @$mode;
    prime_set_config(gmp => $gmp);
    my @fail;
    my ($nvalues, $nerrors) = (0,0);

    for my $i (0 .. $#trees) {
      my $got = _outcome(\&_eval_xop_tree, $trees[$i]);
      my $exp = _outcome(\&_eval_direct_tree, $trees[$i]);
      $got->[0] eq 'value' ? $nvalues++ : $nerrors++;
      push @fail, "case $i " . _show_tree($trees[$i]) .
                  ": [" . join(',',@$got) . "] != [" .
                  join(',',@$exp) . "]"
        if join("\x1E",@$got) ne join("\x1E",@$exp);
    }
    is(scalar(@fail), 0, "$label: $NTREES expression trees");
    ok($nvalues > 0 && $nerrors > 0,
       "$label corpus covers successful and exceptional expressions");
    diag($_) for @fail[0 .. ($#fail < 4 ? $#fail : 4)];
  }
  prime_set_config(gmp => $old_gmp ? 1 : 0);
};

done_testing();
