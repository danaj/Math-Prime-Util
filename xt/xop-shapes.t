#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/:all/;

sub strs {
  return [ map { "$_" } @_ ];
}

my $BIG_S    = '7' x 30;
my $BIG_P1_S = ('7' x 29) . '8';
my $BIG_P2_S = ('7' x 29) . '9';
my $BIG_N_S  = '-' . $BIG_S;

my $BIG    = toint($BIG_S);
my $BIG_P1 = toint($BIG_P1_S);
my $BIG_P2 = toint($BIG_P2_S);
my $BIG_N  = toint($BIG_N_S);

subtest 'addint shapes' => sub {
  is("".addint(56, 7), "63", 'native literals');
  is("".addint("56", "7"), "63", 'decimal strings');
  is("".addint(is_prime(7), 5), "7", 'constant-producing expression');

  {
    my $x = 56;
    is_deeply(strs(map { addint($x, $_) } 0 .. 4),
              [qw/56 57 58 59 60/],
              '$x and $_');
    is_deeply(strs(map { addint($x + $_, 3) } 0 .. 4),
              [qw/59 60 61 62 63/],
              '$x+$_ expression');
  }

  {
    my @refs = map { [$_] } (-2 .. 2);
    is_deeply(strs(map { addint($_->[0], 10) } @refs),
              [qw/8 9 10 11 12/],
              '$_->[0]');
  }

  {
    my @x = (10, 20, 30, 40, 50);
    is_deeply(strs(map { addint($x[$_+1], 1) } 0 .. 2),
              [qw/21 31 41/],
              '$x[$_+1]');
  }

  is("".addint($BIG, 1), $BIG_P1_S, 'bigint lexical');
  is("".addint($BIG, is_prime(7)), $BIG_P2_S, 'bigint plus constant-producing expression');
  is("".addint($BIG_N, 1), '-' . ('7' x 29) . '6', 'negative bigint');
};

subtest 'add1int shapes' => sub {
  is("".add1int(56), "57", 'native literal');
  is("".add1int("56"), "57", 'decimal string');
  is("".add1int(is_prime(7)), "3", 'constant-producing expression');

  {
    my $x = 56;
    is_deeply(strs(map { add1int($x + $_) } 0 .. 4),
              [qw/57 58 59 60 61/],
              '$x+$_ expression');
  }

  {
    my @refs = map { [$_] } (-2 .. 2);
    is_deeply(strs(map { add1int($_->[0]) } @refs),
              [qw/-1 0 1 2 3/],
              '$_->[0]');
  }

  {
    my @x = (10, 20, 30, 40, 50);
    is_deeply(strs(map { add1int($x[$_+1]) } 0 .. 2),
              [qw/21 31 41/],
              '$x[$_+1]');
  }

  is("".add1int($BIG), $BIG_P1_S, 'bigint lexical');
  is_deeply(strs(map { add1int(signint($_)) } -2 .. 2),
            [qw/0 0 1 2 2/],
            'add1int(signint($_))');

  {
    my @got;
    forprimes { push @got, "".add1int($_) } 10;
    is_deeply(\@got, [qw/3 4 6 8/], 'forprimes { add1int($_) } 10');
  }
};

subtest 'signint shapes' => sub {
  is(signint(-56), -1, 'native literal');
  is(signint("0"), 0, 'decimal string zero');
  is(signint(is_prime(7)), 1, 'constant-producing expression');

  is_deeply([ map { signint($_) } (-2 .. 2) ],
            [-1, -1, 0, 1, 1],
            '$_');

  {
    my $x = -1;
    is_deeply([ map { signint($x + $_) } (-2 .. 2) ],
              [-1, -1, -1, 0, 1],
              '$x+$_ expression');
  }

  {
    my @refs = ([-2], [-1], [0], [1], [2], [$BIG], [$BIG_N]);
    is_deeply([ map { signint($_->[0]) } @refs ],
              [-1, -1, 0, 1, 1, 1, -1],
              '$_->[0] with bigint values');
  }

  {
    my @x = (-1, 0, 1, 2, 3);
    is_deeply([ map { signint($x[$_+1]) } 0 .. 3 ],
              [0, 1, 1, 1],
              '$x[$_+1]');
  }

  is(signint($BIG), 1, 'bigint lexical');
  is(signint($BIG_N), -1, 'negative bigint lexical');
};

subtest 'cmpint shapes' => sub {
  is(cmpint(56, 57), -1, 'native literals');
  is(cmpint("56", "56"), 0, 'decimal strings');
  is(cmpint(is_prime(7), 2), 0, 'constant-producing expression');

  {
    my $x = 56;
    is_deeply([ map { cmpint($x, $_) } 54 .. 58 ],
              [1, 1, 0, -1, -1],
              '$x and $_');
    is_deeply([ map { cmpint($x + $_, 58) } 0 .. 4 ],
              [-1, -1, 0, 1, 1],
              '$x+$_ expression');
  }

  {
    my @refs = ([1,2], [2,1], [2,2], [$BIG, $BIG_P1], [$BIG_N, $BIG]);
    is_deeply([ map { cmpint($_->[0], $_->[1]) } @refs ],
              [-1, 1, 0, -1, -1],
              '$_->[0], $_->[1] with bigint values');
  }

  {
    my @x = (10, 20, 30, 40, 50);
    is_deeply([ map { cmpint($x[$_+1], 30) } 0 .. 2 ],
              [-1, 0, 1],
              '$x[$_+1]');
  }

  is(cmpint($BIG, $BIG), 0, 'same bigint');
  is(cmpint($BIG, $BIG_P1), -1, 'bigint ordering');
  is(cmpint($BIG_N, $BIG), -1, 'negative bigint ordering');
};

done_testing();
