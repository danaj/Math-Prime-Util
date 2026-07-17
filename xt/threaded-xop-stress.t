#!/usr/bin/env perl
use strict;
use warnings;

use Config;
BEGIN {
  if (!$Config{useithreads} || $] < 5.008) {
    print "1..0 # Skip perl isn't compiled with threading support\n";
    exit 0;
  }
  if (!eval { require threads; 1 }) {
    print "1..0 # Skip threads.pm not installed\n";
    exit 0;
  }
}

# Avoid thread-unsafe optional bigint backends.  The custom ops may still use
# the MPU GMP backend when it is available.
use Math::BigInt lib => 'Calc';
use B ();
use Test::More;
use Math::Prime::Util qw/
  forprimes
  irand irand64 drand
  addint subint add1int sub1int mulint divint modint cdivint powint
  is_odd is_even is_square cmpint kronecker signint
/;

my $NTHREADS = $ENV{MPU_XOP_THREADS} || 4;
my $NWAVES   = $ENV{MPU_XOP_WAVES}   || 3;
my $NROUNDS  = $ENV{MPU_XOP_ROUNDS}  || 160;
my $MOD      = 1_000_003;

sub _mix {
  my ($state, @values) = @_;
  for my $value (@values) {
    my $string = defined($value) ? "$value" : '<undef>';
    $state = ($state * 33 + $_) % $MOD for unpack('C*', $string);
    $state = ($state * 33 + 255) % $MOD;
  }
  return $state;
}

sub _inputs {
  my ($seed, $i) = @_;
  my $a = ($seed + 17*$i) % 10_007 + 1;
  my $b = (31*$i + 7) % 251 + 1;
  return ($a, $b);
}

sub _scalar_shape {
  my ($seed, $rounds) = @_;
  Math::Prime::Util::srand($seed);
  my $digest = 17;

  for my $i (1 .. $rounds) {
    my ($a, $b) = _inputs($seed, $i);
    my $d = drand;
    die "drand returned $d" if $d < 0 || $d >= 1;
    $digest = _mix(
      $digest,
      irand, irand64, int($d * 1_000_000),
      addint($a,$b), subint($a,$b), add1int($a), sub1int($a),
      mulint($a,$b), divint($a,$b), modint($a,$b), cdivint($a,$b),
      powint(($a % 11) + 2, $i % 6),
      is_odd($a), is_even($a), is_square($a),
      cmpint($a,$b), kronecker($a,2*$b+1), signint($a-$b),
    );
  }
  return $digest;
}

sub _void_shape {
  my ($seed, $rounds) = @_;
  Math::Prime::Util::srand($seed);
  my $digest = 29;

  for my $i (1 .. $rounds) {
    my ($a, $b) = _inputs($seed, $i);

    irand;
    irand64;
    drand;
    addint($a,$b);
    subint($a,$b);
    add1int($a);
    sub1int($a);
    mulint($a,$b);
    divint($a,$b);
    modint($a,$b);
    cdivint($a,$b);
    powint(($a % 7) + 2, $i % 5);
    is_odd($a);
    is_even($a);
    is_square($a);
    cmpint($a,$b);
    kronecker($a,2*$b+1);
    signint($a-$b);

    my $d = drand;
    die "drand returned $d after void calls" if $d < 0 || $d >= 1;
    $digest = _mix($digest, irand, irand64, int($d * 1_000_000),
                  addint($a,$b));
  }
  return $digest;
}

sub _map_shape {
  my ($seed, $rounds) = @_;
  Math::Prime::Util::srand($seed);
  my $digest = 43;

  for my $i (1 .. int(($rounds+7)/8)) {
    my @input = map { 8*$i + $_ } 0 .. 7;
    my @arith = map {
      addint(mulint($_,3), sub1int($_))
    } @input;
    my @flags = map {
      (is_odd($_), is_even($_), is_square(mulint($_,$_)),
       cmpint($_,20), kronecker($_,2*$_+1), signint(subint($_,20)))
    } @input;
    my @random = (
      map { irand } 1 .. 3,
      map { irand64 } 1 .. 3,
      map { int(drand() * 1_000_000) } 1 .. 3,
    );
    $digest = _mix($digest, @arith, @flags, @random);
  }
  return $digest;
}

sub _callback_body {
  my ($a, $b) = @_;
  my $d = drand;
  die "callback drand returned $d" if $d < 0 || $d >= 1;
  return (
    addint(mulint(add1int($a),sub1int($b)), kronecker($a,2*$b+1)),
    divint(addint($a,$b),$b),
    modint(addint($a,$b),$b),
    cdivint(addint($a,$b),$b),
    powint(($a % 5) + 2, $b % 5),
    is_odd($a), is_even($a), is_square(mulint($a,$a)),
    cmpint($a,$b), signint(subint($a,$b)),
    irand, irand64, int($d * 1_000_000),
  );
}

sub _callback_shape {
  my ($seed, $rounds) = @_;
  Math::Prime::Util::srand($seed);
  my $digest = 71;

  for my $i (1 .. int(($rounds+15)/16)) {
    my $lo = 100 + 20*$i;
    forprimes {
      $digest = _mix($digest, _callback_body($_,$i+2));
    } $lo, $lo+31;
  }
  return $digest;
}

sub _nested_shape {
  my ($seed, $rounds) = @_;
  Math::Prime::Util::srand($seed);
  my $digest = 97;

  for my $i (1 .. $rounds) {
    my ($a, $b) = _inputs($seed, $i);
    my $v = addint(
      mulint(add1int($a),sub1int($b)),
      subint(powint(($a % 7) + 2,$i % 5),cdivint($a,$b)),
    );
    my $root = addint($a,$b);
    my $square = is_square(mulint($root,$root));
    my $symbol = kronecker(addint($a,2),add1int(mulint(2,$b)));
    my $random = addint(modint(irand,997),modint(irand64,991));
    my $dclass = cmpint(int(drand() * 1_000_000),500_000);
    $digest = _mix($digest, $v, $square, $symbol, $random, $dclass);
  }

  # Exercise the bigint fallback paths without dominating the stress test.
  my $big = '184467440737095516160';
  my $big_square = '34028236692093846346337460743176821145600';
  my @big_results;
  push @big_results, addint($big,17);
  push @big_results, subint($big,17);
  push @big_results, mulint($big,3);
  push @big_results, divint($big,7);
  push @big_results, modint($big,97);
  push @big_results, cdivint($big,7);
  push @big_results, powint($big,2);
  push @big_results, is_square($big_square);
  push @big_results, cmpint($big_square,$big);
  push @big_results, kronecker($big,97);
  push @big_results, signint(subint(0,$big));
  $digest = _mix($digest, @big_results);
  return $digest;
}

sub _bigint_clone_probe {
  return "" . addint('184467440737095516160',17);
}

sub _worker {
  my ($seed, $rounds) = @_;
  return {
    scalar   => _scalar_shape  ($seed + 0x0101, $rounds),
    void     => _void_shape    ($seed + 0x0202, $rounds),
    map      => _map_shape     ($seed + 0x0303, $rounds),
    callback => _callback_shape($seed + 0x0404, $rounds),
    nested   => _nested_shape  ($seed + 0x0505, $rounds),
  };
}

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

my @worker_cvs = (
  \&_scalar_shape, \&_void_shape, \&_map_shape,
  \&_callback_body, \&_callback_shape, \&_nested_shape,
  \&_bigint_clone_probe,
);
my @xops = qw/
  irand irand64 drand
  addint subint add1int sub1int mulint divint modint cdivint powint
  is_odd is_even is_square cmpint kronecker signint
/;
my $compiled = _compiled_ops(@worker_cvs);

plan skip_all => 'custom ops are not available in this build'
  if !grep { $compiled->{$_} } @xops;

subtest 'worker call sites use custom ops' => sub {
  my @missing = grep { !$compiled->{$_} } @xops;
  is_deeply(\@missing, [], 'all requested functions compiled as custom ops');
};

subtest 'configured bigint class survives cloning' => sub {
  my $expected = _bigint_clone_probe();
  my $thread = threads->create(\&_bigint_clone_probe);
  is($thread->join, $expected, 'child returns a canonical bigint');
};

diag("$NTHREADS threads, $NWAVES lifecycle waves, $NROUNDS rounds per shape");

subtest 'same seed has independent deterministic streams' => sub {
  my $seed = 0x4512_0000;
  my $expected = _worker($seed, $NROUNDS);
  my @threads = map {
    threads->create(\&_worker, $seed, $NROUNDS)
  } 1 .. $NTHREADS;

  for my $i (0 .. $#threads) {
    my $got = $threads[$i]->join;
    is_deeply($got, $expected, "thread " . ($i+1));
  }
};

for my $wave (1 .. $NWAVES) {
  subtest "parallel lifecycle wave $wave" => sub {
    my @seeds = map { 0x6000_0000 + 0x1000*$wave + 101*$_ }
                1 .. $NTHREADS;
    my @threads = map {
      threads->create(\&_worker, $_, $NROUNDS)
    } @seeds;

    # Run the reference work while the child interpreters are active.
    my @expected = map { _worker($_, $NROUNDS) } @seeds;
    for my $i (0 .. $#threads) {
      my $got = $threads[$i]->join;
      for my $shape (qw/scalar void map callback nested/) {
        is($got->{$shape}, $expected[$i]{$shape},
           "thread " . ($i+1) . " $shape context");
      }
    }
  };
}

done_testing();
