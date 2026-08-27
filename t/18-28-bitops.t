#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util ();

my @functions = qw(
  bitand bitor bitxor bitandnot bitnot
  bitset bitclear bitflip bittest bitscan0 bitscan1
);
my $package = 'Math::Prime::Util';
my %func = map { $_ => $package->can($_) } @functions;

sub bcall {
  my $name = shift;
  $func{$name}->(@_);
}
sub as_strings { [map { defined($_) ? "$_" : 'undef' } @_] }
sub dies {
  my($name, @args) = @_;
  return !eval { bcall($name, @args); 1 };
}

can_ok($package, @functions);

subtest 'binary operations on small magnitudes' => sub {
  my(@and, @or, @xor, @andnot);
  my(@eand, @eor, @exor, @eandnot);
  for my $x (0 .. 31) {
    for my $y (0 .. 31) {
      push @and,    bcall('bitand',    $x, $y);
      push @or,     bcall('bitor',     $x, $y);
      push @xor,    bcall('bitxor',    $x, $y);
      push @andnot, bcall('bitandnot', $x, $y);
      push @eand,    $x & $y;
      push @eor,     $x | $y;
      push @exor,    $x ^ $y;
      push @eandnot, $x & ~$y;
    }
  }
  is_deeply(\@and,    \@eand,    'bitand truth table');
  is_deeply(\@or,     \@eor,     'bitor truth table');
  is_deeply(\@xor,    \@exor,    'bitxor truth table');
  is_deeply(\@andnot, \@eandnot, 'bitandnot truth table');

  is_deeply(
    [map { bcall($_, -13, -10) } qw/bitand bitor bitxor bitandnot/],
    [8, 15, 7, 5],
    'negative inputs use their magnitudes'
  );
};

subtest 'binary operations on big integers' => sub {
  my $a = '340282366920938463500268095579187314693';
  my $b = '170141183460469231768580791863303208963';
  my @got = map { bcall($_,$a,$b) }
                qw/bitand bitor bitxor bitandnot/;
  is_deeply(as_strings(@got), [qw(
    36893488147419103233
    510423550381407695231955399295071420423
    510423550381407695195061911147652317190
    340282366920938463463374607431768211460
  )], 'selected 129-bit results');

  for my $x (0, 1, 15, $a, "-$a") {
    (my $abs = "$x") =~ s/^-//;
    is("".bcall('bitand',$x,$x), $abs, 'x AND x is abs(x)');
    is("".bcall('bitor',$x,0), $abs, 'x OR 0 is abs(x)');
    is(bcall('bitxor',$x,$x), 0, 'x XOR x is zero');
    is(bcall('bitandnot',$x,$x), 0, 'x AND NOT x is zero');
  }

  my $large = '1234567890' x 200;
  is(''.bcall('bitor',$large,0), $large,
     '2000-digit binary conversion round trip');
};

subtest 'bitnot width semantics' => sub {
  my @natural = (
    [0,1], [1,0], [2,1], [3,0], [4,3], [5,2], [7,0], [8,7],
    [10,5], [15,0], [16,15], [255,0], [-10,5],
  );
  is_deeply(
    [map { bcall('bitnot',$_->[0]) } @natural],
    [map { $_->[1] } @natural],
    'natural-width complements'
  );

  my @fixed = (
    [0,0,0], [0,1,1], [0,4,15], [5,0,0], [5,1,0], [5,2,2],
    [5,3,2], [5,4,10], [10,3,5], [15,3,0], [16,4,15], [-5,4,10],
  );
  is_deeply(
    [map { bcall('bitnot',$_->[0],$_->[1]) } @fixed],
    [map { $_->[2] } @fixed],
    'fixed-width complements truncate to width'
  );

  my $a = '340282366920938463500268095579187314693';
  is(''.bcall('bitnot',$a),
     '340282366920938463426481119284349108218',
     'natural-width complement of 129-bit input');
  is(''.bcall('bitnot',$a,64), '18446744073709551610',
     'fixed-width bigint complement');
  is(''.bcall('bitnot',$a,130),
     '1020847100762815390353230334147885531130',
     'fixed width can exceed input width');
};

subtest 'single-bit operations' => sub {
  my @k = 0 .. 5;
  is_deeply([map { bcall('bittest',10,$_) } @k], [0,1,0,1,0,0],
            'bittest positions');
  is_deeply([map { bcall('bitset',10,$_) } @k], [11,10,14,10,26,42],
            'bitset positions');
  is_deeply([map { bcall('bitclear',10,$_) } @k], [10,8,10,2,10,10],
            'bitclear positions');
  is_deeply([map { bcall('bitflip',10,$_) } @k], [11,8,14,2,26,42],
            'bitflip positions');
  is_deeply([map { bcall($_,-10,1) }
                  qw/bittest bitset bitclear bitflip/], [1,10,8,8],
            'single-bit functions use magnitude');

  my $a = '340282366920938463500268095579187314693';
  is_deeply([map { bcall('bittest',$a,$_) } (0,1,2,64,65,127,128,129)],
            [1,0,1,0,1,0,1,0], 'bittest on 129-bit input');
  is(''.bcall('bitset',$a,200),
     '1606938044258990275542302374708083540985703261878372022616069',
     'bitset above bigint high bit');
  is(''.bcall('bitclear',$a,65),
     '340282366920938463463374607431768211461',
     'bitclear bigint bit');
  is(''.bcall('bitflip',$a,64),
     '340282366920938463518714839652896866309',
     'bitflip bigint bit');
  is_deeply(as_strings(map { bcall($_,"-$a",65) }
                            qw/bittest bitset bitclear bitflip/),
            [1, $a,
             '340282366920938463463374607431768211461',
             '340282366920938463463374607431768211461'],
            'single-bit functions use bigint magnitude');
  is(''.bcall('bitclear',$a,200), $a,
     'clearing a bit above the input is an identity');

  is(''.bcall('bitset',0,31), '2147483648',
     'bitset crosses the 32-bit signed boundary');
  is(''.bcall('bitset',0,63), '9223372036854775808',
     'bitset crosses the 64-bit signed boundary');
  is(''.bcall('bitclear','9223372036854775808',63), '0',
     'bitclear clears the 64-bit sign position');
};

subtest 'bit scans' => sub {
  my(@got0,@got1,@want0,@want1);
  for my $x (0 .. 255) {
    for my $start (0 .. 10) {
      push @got0, bcall('bitscan0',$x,$start);
      push @got1, bcall('bitscan1',$x,$start);

      my($z,$o);
      for my $k ($start .. 16) {
        my $bit = ($x >> $k) & 1;
        $z = $k if !defined($z) && !$bit;
        $o = $k if !defined($o) &&  $bit;
      }
      push @want0, $z;
      push @want1, $o;
    }
  }
  is_deeply(as_strings(@got0), as_strings(@want0),
            'bitscan0 small exhaustive values');
  is_deeply(as_strings(@got1), as_strings(@want1),
            'bitscan1 small exhaustive values');

  my $a = '340282366920938463500268095579187314693';
  is_deeply(as_strings(map { bcall('bitscan1',$a,$_) }
                            (0,1,2,3,65,66,128,129)),
            [qw/0 2 2 65 65 128 128 undef/],
            'bitscan1 bigint starts');
  is_deeply(as_strings(map { bcall('bitscan0',$a,$_) }
                            (0,1,2,64,65,127,128,129)),
            [qw/1 1 3 64 66 127 129 129/],
            'bitscan0 bigint starts');
  is_deeply(as_strings(map { bcall('bitscan1',"-$a",$_) }
                            (0,3,65,129)),
            [qw/0 65 65 undef/],
            'bitscan1 uses bigint magnitude');
  is_deeply(as_strings(map { bcall('bitscan0',"-$a",$_) }
                            (0,64,128,129)),
            [qw/1 64 129 129/],
            'bitscan0 uses bigint magnitude');

  my $ones128 = '340282366920938463463374607431768211455';
  is_deeply(as_strings(
      bcall('bitscan0',$ones128),
      bcall('bitscan0',$ones128,64),
      bcall('bitscan0',$ones128,128),
      bcall('bitscan1',$ones128,127),
      bcall('bitscan1',$ones128,128)),
    [qw/128 128 128 127 undef/], 'scan across the bigint high bit');
};

subtest 'single-bit identities' => sub {
  my @values = (0, 1, 2, 7, 255, -255, '4294967296',
                '18446744073709551615',
                '340282366920938463463374607431768211457');
  for my $x (@values) {
    (my $abs = "$x") =~ s/^-//;
    for my $k (0, 1, 5, 31, 63, 127) {
      is(bcall('bittest',bcall('bitset',$x,$k),$k), 1,
         'bitset makes bittest true');
      is(bcall('bittest',bcall('bitclear',$x,$k),$k), 0,
         'bitclear makes bittest false');
      is(''.bcall('bitflip',bcall('bitflip',$x,$k),$k), $abs,
         'bitflip is an involution');
    }
  }
};

subtest 'native word boundaries' => sub {
  my $uvbits = Math::Prime::Util::_uvbits();
  my $uvmax = $uvbits == 64 ? '18446744073709551615' : '4294967295';
  my $high = $uvbits == 64 ? '9223372036854775808' : '2147483648';

  is(bcall('bitnot',$uvmax), 0, 'natural bitnot at UV_MAX');
  is(''.bcall('bitnot',0,$uvbits), $uvmax,
     'fixed bitnot fills a native word');
  is(''.bcall('bitset',0,$uvbits-1), $high,
     'bitset sets the highest native bit');
  is(bcall('bitclear',$high,$uvbits-1), 0,
     'bitclear clears the highest native bit');
  is(bcall('bitscan0',$uvmax), $uvbits,
     'bitscan0 finds the first implicit zero after UV_MAX');
  is(bcall('bitscan1',$uvmax,$uvbits-1), $uvbits-1,
     'bitscan1 finds the highest native bit');
  is(bcall('bitscan1',$uvmax,$uvbits), undef,
     'bitscan1 stops beyond the native word');
};

subtest 'maximum supported index' => sub {
  my $max = 4294967295;
  my $big = '340282366920938463463374607431768211457';
  is(bcall('bittest',123,$max), 0, 'bittest accepts UINT32_MAX');
  is(bcall('bitclear',123,$max), 123, 'bitclear accepts UINT32_MAX');
  is(bcall('bittest',$big,$max), 0,
     'bittest accepts UINT32_MAX with bigint input');
  is(''.bcall('bitclear',$big,$max), $big,
     'bitclear accepts UINT32_MAX with bigint input');
  is(bcall('bitscan0',0,$max), $max, 'bitscan0 accepts UINT32_MAX');
  is(bcall('bitscan1',0,$max), undef, 'bitscan1 accepts UINT32_MAX');

  my $too_large = '4294967296';
  for my $name (qw/bitset bitclear bitflip bittest/) {
    ok(dies($name,1,$too_large), "$name rejects index above UINT32_MAX");
  }
  ok(dies('bitnot',1,$too_large), 'bitnot rejects width above UINT32_MAX');
  ok(dies('bitscan0',1,$too_large), 'bitscan0 rejects start above UINT32_MAX');
  ok(dies('bitscan1',1,$too_large), 'bitscan1 rejects start above UINT32_MAX');
};

subtest 'validation and arity' => sub {
  for my $name (qw/bitand bitor bitxor bitandnot/) {
    ok(dies($name,1), "$name rejects one argument");
    ok(dies($name,1,2,3), "$name rejects three arguments");
    ok(dies($name,'invalid',1), "$name validates x");
    ok(dies($name,1,'invalid'), "$name validates y");
  }
  for my $name (qw/bitset bitclear bitflip bittest/) {
    ok(dies($name,1), "$name rejects one argument");
    ok(dies($name,1,2,3), "$name rejects three arguments");
    ok(dies($name,'invalid',1), "$name validates x");
    ok(dies($name,1,-1), "$name rejects negative index");
    ok(dies($name,1,1.5), "$name rejects noninteger index");
  }
  for my $name (qw/bitnot bitscan0 bitscan1/) {
    ok(dies($name), "$name rejects no arguments");
    ok(dies($name,1,2,3), "$name rejects three arguments");
    ok(dies($name,'invalid'), "$name validates x");
    ok(dies($name,1,-1), "$name rejects negative optional argument");
    ok(dies($name,1,1.5), "$name rejects noninteger optional argument");
  }
};

subtest 'object inputs are not modified' => sub {
  my $a = Math::BigInt->new('-340282366920938463500268095579187314693');
  my $b = Math::BigInt->new('-170141183460469231768580791863303208963');
  my($sa,$sb) = ("$a","$b");
  bcall('bitand',$a,$b);
  bcall('bitnot',$a);
  bcall('bitset',$a,200);
  bcall('bitclear',$a,65);
  bcall('bitflip',$a,64);
  bcall('bittest',$a,128);
  bcall('bitscan0',$a);
  bcall('bitscan1',$a);
  is("$a",$sa, 'first object input unchanged');
  is("$b",$sb, 'second object input unchanged');
};

done_testing();
