#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/
  prime_get_config addint setbinop
  forprimes foroddcomposites forcomposites forsemiprimes
  foralmostprimes fordivisors forcomb forperm forderange
  forsetproduct forfactored forsquarefree forsquarefreeint
  vecreduce vecslide vecwindow vecpairwise
  vecnone vecall vecany vecnotall vecfirst vecfirstidx
/;

plan skip_all => 'multicall reentrancy requires XS'
  unless prime_get_config()->{'xs'};

# Keep these cases aligned with every public alias of an SC_dMULTICALL site
# in XS.xs.  Calling through a coderef prevents addint from becoming an xop.
our $ADD = \&Math::Prime::Util::addint;
our $BIG = '18446744073709551616';
our $MAXUV = prime_get_config()->{'maxbits'} == 32
              ? '4294967295'
              : '18446744073709551615';
our $sum;

{
  package MPU::Test::ScopeGuard;
  our @at_destroy;
  sub new { bless {}, shift }
  sub DESTROY { push @at_destroy, "$main::a,$main::b" }
}

sub big_plus {
  return ''. $ADD->($BIG, $_[0]);
}

my @cases = (
  [ setbinop => sub {
      return [map {"$_"} @{
        setbinop(sub { $ADD->($a, $b) }, [$MAXUV], [1,2])
      }];
    }, [''. $ADD->($MAXUV, 1), ''. $ADD->($MAXUV, 2)] ],

  [ forprimes => sub {
      $sum = $BIG;
      forprimes { $sum = $ADD->($sum, $_) } 2,7;
      return "$sum";
    }, big_plus(17) ],

  [ foroddcomposites => sub {
      $sum = $BIG;
      foroddcomposites { $sum = $ADD->($sum, $_) } 1,25;
      return "$sum";
    }, big_plus(70) ],

  [ forcomposites => sub {
      $sum = $BIG;
      forcomposites { $sum = $ADD->($sum, $_) } 1,10;
      return "$sum";
    }, big_plus(37) ],

  [ forsemiprimes => sub {
      $sum = $BIG;
      forsemiprimes { $sum = $ADD->($sum, $_) } 1,15;
      return "$sum";
    }, big_plus(58) ],

  [ foralmostprimes => sub {
      $sum = $BIG;
      foralmostprimes { $sum = $ADD->($sum, $_) } 3,1,20;
      return "$sum";
    }, big_plus(58) ],

  [ fordivisors => sub {
      $sum = $BIG;
      fordivisors { $sum = $ADD->($sum, $_) } 12;
      return "$sum";
    }, big_plus(28) ],

  [ forcomb => sub {
      $sum = $BIG;
      forcomb {
        $sum = $ADD->($sum, $ADD->($_[0], $_[-1]));
      } 4,2;
      return "$sum";
    }, big_plus(18) ],

  [ forperm => sub {
      $sum = $BIG;
      forperm {
        $sum = $ADD->($sum, $ADD->($_[0], $_[-1]));
      } 3;
      return "$sum";
    }, big_plus(12) ],

  [ forderange => sub {
      $sum = $BIG;
      forderange {
        $sum = $ADD->($sum, $ADD->($_[0], $_[-1]));
      } 4;
      return "$sum";
    }, big_plus(27) ],

  [ forsetproduct => sub {
      $sum = $BIG;
      forsetproduct {
        $sum = $ADD->($sum, $ADD->($_[0], $_[-1]));
      } [1,2], [10,20];
      return "$sum";
    }, big_plus(66) ],

  [ forfactored => sub {
      $sum = $BIG;
      forfactored {
        $sum = $ADD->($sum, $ADD->($_, scalar @_));
      } 2,10;
      return "$sum";
    }, big_plus(69) ],

  [ forsquarefree => sub {
      $sum = $BIG;
      forsquarefree {
        $sum = $ADD->($sum, $ADD->($_, scalar @_));
      } 2,10;
      return "$sum";
    }, big_plus(41) ],

  [ forsquarefreeint => sub {
      $sum = $BIG;
      forsquarefreeint { $sum = $ADD->($sum, $_) } 2,10;
      return "$sum";
    }, big_plus(33) ],

  [ vecreduce => sub {
      return ''. vecreduce(sub { $ADD->($a, $b) }, $BIG, 1, 2, 3);
    }, big_plus(6) ],

  [ vecslide => sub {
      return [map {"$_"}
        vecslide(sub { $ADD->($a, $b) }, $BIG, 1, 2)
      ];
    }, [big_plus(1), '3'] ],

  [ vecwindow => sub {
      return [map {"$_"}
        vecwindow(sub { $ADD->($_[0], $_[-1]) }, 1, 2, $BIG, 1, 2)
      ];
    }, [big_plus(1), '3'] ],

  [ vecpairwise => sub {
      return [map {"$_"}
        vecpairwise(sub { $ADD->($a, $b) }, [$BIG,2], [1,3])
      ];
    }, [big_plus(1), '5'] ],

  [ vecnone => sub {
      return vecnone(sub { $ADD->($_, 0) < 0 }, $BIG, 1, 2);
    }, 1 ],

  [ vecall => sub {
      return vecall(sub { $ADD->($_, 0) > 0 }, $BIG, 1, 2);
    }, 1 ],

  [ vecany => sub {
      return vecany(sub { $ADD->($_, 0) == 2 }, $BIG, 1, 2);
    }, 1 ],

  [ vecnotall => sub {
      return vecnotall(sub { $ADD->($_, 0) > 0 }, $BIG, 1, 2);
    }, 0 ],

  [ vecfirst => sub {
      return vecfirst(sub { $ADD->($_, 0) == 2 }, $BIG, 1, 2);
    }, 2 ],

  [ vecfirstidx => sub {
      return vecfirstidx(sub { $ADD->($_, 0) == 2 }, $BIG, 1, 2);
    }, 2 ],
);

plan tests => scalar(@cases) + 1;

for my $case (@cases) {
  my($name, $code, $expected) = @$case;
  my $got = eval { $code->() };
  my $error = $@;
  if ($error ne '') {
    fail("$name nested addint callback");
    diag($error);
  } else {
    is_deeply($got, $expected, "$name nested addint callback");
  }
}

subtest 'scoped callback return values' => sub {
  is_deeply(
    setbinop(sub {
      my $v = $a + $b;
      $a = 0;
      $b = 0;
      return $v;
    }, [1,2], [10,20]),
    [11,12,21,22],
    'setbinop preserves a lexical scalar return'
  );

  is(
    vecreduce(sub {
      my $v = $a + $b;
      return $v;
    }, 1,2,3,4),
    10,
    'vecreduce preserves a lexical scalar return'
  );

  is_deeply(
    [vecslide(sub {
      my $v = $a + $b;
      return $v;
    }, 1,2,3)],
    [3,5],
    'vecslide preserves lexical scalar returns'
  );

  is_deeply(
    [vecwindow(sub {
      my $v = $_[0] + $_[1];
      return ($v, -$v);
    }, 1,2,1,2,3)],
    [3,-3,5,-5],
    'vecwindow preserves lexical list returns'
  );

  is_deeply(
    [vecpairwise(sub {
      my $v = $a + $b;
      return ($v,-$v) if $a == 1;
      return ()       if $a == 2;
      return $v;
    }, [1,2,3], [10,20,30])],
    [11,-11,33],
    'vecpairwise preserves varying lexical list returns'
  );

  is(
    vecany(sub {
      my $v = $_ == 2;
      return $v;
    }, 0,1,2),
    1,
    'vecany preserves a lexical truth return'
  );

  my @window_large = vecwindow(sub {
    my @v = 1 .. 10_000;
    return @v;
  }, 1,1,7);
  is(scalar @window_large, 10_000,
     'vecwindow preserves a list return that extends the Perl stack');
  is("$window_large[0] $window_large[-1]", '1 10000',
     'vecwindow large list return has intact endpoints');

  my @pairwise_large = vecpairwise(sub {
    my @v = 1 .. 10_000;
    return @v;
  }, [1], [2]);
  is(scalar @pairwise_large, 10_000,
     'vecpairwise preserves a list return that extends the Perl stack');
  is("$pairwise_large[0] $pairwise_large[-1]", '1 10000',
     'vecpairwise large list return has intact endpoints');

  SKIP: {
    skip 'OP_MULTIPARAM requires Perl 5.43.3', 2
      if $] < 5.043003;

    my $signature_cb = eval q{
      use feature 'signatures';
      no warnings 'experimental::signatures';
      sub ($guard = MPU::Test::ScopeGuard->new) { return $a + $b }
    };
    die $@ if !$signature_cb;

    @MPU::Test::ScopeGuard::at_destroy = ();
    is_deeply(
      [&Math::Prime::Util::vecpairwise($signature_cb, [1,2], [10,20])],
      [11,22],
      'vecpairwise accepts a callback using the new signature opcodes'
    );
    is_deeply(
      \@MPU::Test::ScopeGuard::at_destroy,
      ['1,10', '2,20'],
      'OP_MULTIPARAM signature pads are scoped per callback'
    );
  }

  done_testing();
};
