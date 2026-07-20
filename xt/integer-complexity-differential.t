#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/integer_complexity/;
use Math::Prime::Util::PP;

my $N = $ENV{MPU_INTEGER_COMPLEXITY_MAX} || 20_000;

sub _reset {
  integer_complexity(0);
  Math::Prime::Util::PP::integer_complexity(0);
}

subtest 'XS and PP sequential agreement' => sub {
  _reset();
  my @fail;
  for my $n (1 .. $N) {
    my $xs = integer_complexity($n);
    my $pp = Math::Prime::Util::PP::integer_complexity($n);
    push @fail, "$n: XS=$xs PP=$pp" if $xs != $pp;
    last if @fail >= 10;
  }
  is_deeply(\@fail, [], "integer_complexity agrees through $N");
};

subtest 'cache growth order and reset' => sub {
  my @queries = (
    $N, 1, $N - 1, 2, int($N / 2), 3,
    map { 1 + (($_ * 7919) % $N) } 1 .. 200,
  );
  my @before = map { integer_complexity($_) } @queries;

  _reset();
  my @after = map { integer_complexity($_) } @queries;
  is_deeply(\@after, \@before,
            'results are independent of cache expansion and reset');

  my @pp = map { Math::Prime::Util::PP::integer_complexity($_) } @queries;
  is_deeply(\@pp, \@before, 'PP agrees after sparse cache expansion');
};

subtest 'basic complexity inequalities' => sub {
  my @fail;
  for my $n (2 .. $N) {
    push @fail, "$n exceeds c(n-1)+1"
      if integer_complexity($n) > integer_complexity($n - 1) + 1;
    last if @fail >= 10;
  }
  for my $i (1 .. 1000) {
    my $a = 2 + (($i * 37) % 200);
    my $b = 2 + (($i * 61) % 200);
    push @fail, "$a*$b exceeds c(a)+c(b)"
      if integer_complexity($a * $b) >
         integer_complexity($a) + integer_complexity($b);
    last if @fail >= 10;
  }
  is_deeply(\@fail, [], 'addition and multiplication upper bounds hold');
};

done_testing();
