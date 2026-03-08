#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/contfrac from_contfrac/;

plan tests => 7;

subtest 'contfrac and from_contfrac roundtrip' => sub {
  my @data = (
    [0,1,[0]],
    [1,3,[0,3]],
    [4,11,[0, 2, 1, 3]],
    [67,29,[2, 3, 4, 2]],
    [121,23,[5, 3, 1, 5]],
    [3,4837,[0,1612,3]],
    [0xfff1,0x7fed,[2, 1423, 1, 6, 1, 2]],
    [83116,51639,[1, 1, 1, 1, 1, 3, 1, 1, 2, 2, 4, 1, 2, 1, 1, 1, 3]],
    [9238492834,2398702938777,[0, 259, 1, 1, 1, 3, 1, 7, 2, 3, 7, 2, 1, 1, 2, 4, 2, 1, 10, 5, 3, 1, 5, 6]],
    ["243224233245235253407096734543059","4324213412343432913758138673203834",[0,17,1,3,1,1,12,1,2,33,2,1,1,1,1,49,1,1,1,1,17,34,1,1,304,1,2,1,1,1,2,1,48,1,20,2,3,5,1,1,16,9,1,1,5,1,2,2,7,4,3,1,7,1,1,17,1,1,29,1,12,2,5]],
    [415,93,[4,2,6,7]],
    [649,200,[3,4,12,4]],
    [4,9,[0,2,4]],
  );
  for my $t (@data) {
    my($n,$d,$exp) = @$t;
    is_deeply( [contfrac($n,$d)], $exp, "contfrac($n,$d)" );
    is_deeply( [map{"$_"}from_contfrac(@$exp)], [$n,$d], "from_contfrac(@$exp)" );
  }
};

subtest 'pi convergents' => sub {
  my @data = (
    [22,7,[3,7]],
    [333,106,[3,7,15]],
    [355,113,[3,7,16]],
    [377,120,[3,7,17]],
    [3927,1250,[3,7,16,11]],
    [103993,33102,[3,7,15,1,292]],
    [104348,33215,[3,7,15,1,293]],
    [208341,66317,[3,7,15,1,292,2]],
    [312689,99532,[3,7,15,1,292,1,2]],
    [833719,265381,[3,7,15,1,292,1,1,1,2]],
    [1146408,364913,[3,7,15,1,292,1,1,1,3]],
    [4272943,1360120,[3,7,15,1,292,1,1,1,2,1,3]],
    [80143857,25510582,[3,7,15,1,292,1,1,1,2,1,3,1,14]],
    ["262452630335382199398","83541266890691994833",[3,7,15,1,292,1,1,1,2,1,3,1,14,2,1,1,2,2,2,2,1,84,2,1,1,15,3,13,1,4,2,6,6,99]],
  );
  for my $t (@data) {
    my($n,$d,$exp) = @$t;
    is_deeply( [contfrac($n,$d)], $exp, "contfrac($n/$d)" );
    is_deeply( [map{"$_"}from_contfrac(@$exp)], [$n,$d], "from_contfrac pi convergent" );
  }
};

subtest 'Fibonacci ratio convergents' => sub {
  my @data = (
    # F(n)/F(n+1) gives CF [0, (1)x(n-2), 2] with n terms
    [144,233,[0,1,1,1,1,1,1,1,1,1,1,2]],
    ["7540113804746346429","12200160415121876738",[0,(1)x90,2]],
    ["280571172992510140037611932413038677189525","453973694165307953197296969697410619233826",[0,(1)x198,2]],
  );
  for my $t (@data) {
    my($n,$d,$exp) = @$t;
    is_deeply( [map{"$_"}contfrac($n,$d)], [map{"$_"}@$exp], "contfrac(F/F) len ".scalar(@$exp) );
    is_deeply( [map{"$_"}from_contfrac(@$exp)], [$n,$d], "from_contfrac(F/F) len ".scalar(@$exp) );
  }
};

subtest 'non-coprime inputs' => sub {
  # contfrac reduces to lowest terms; from_contfrac won't recover original
  my @data = (
    [62832,20000,[3,7,16,11]],
    [0,2,[0]],
    [8,22,[0, 2, 1, 3]],
  );
  for my $t (@data) {
    my($n,$d,$exp) = @$t;
    is_deeply( [contfrac($n,$d)], $exp, "contfrac($n,$d) non-coprime" );
  }
};

subtest 'negative numerator' => sub {
  my @data = (
    [-93,37,[-3,2,18]],
    [-312689,99532,[-4,1,6,15,1,292,1,2]],
    [-4,11,[-1,1,1,1,3]],
    [-4,5837,[-1,1,1458,4]],
    [-4,11111,[-1,1,2776,1,3]],
    [-1,11111,[-1,1,11110]],
    [-11110,11111,[-1,11111]],
    [-11112,11111,[-2,1,11110]],
    [-1,1,[-1]],
    [-7,1,[-7]],
    [-100,1,[-100]],
  );
  for my $t (@data) {
    my($n,$d,$exp) = @$t;
    is_deeply( [contfrac($n,$d)], $exp, "contfrac($n,$d)" );
    is_deeply( [from_contfrac(@$exp)], [$n,$d], "from_contfrac(@$exp)" );
  }
};

subtest 'edge cases' => sub {
  # Integers (d=1): CF is a single element
  for my $n (0, 1, 7, 100) {
    is_deeply( [contfrac($n,1)], [$n], "contfrac($n,1)" );
    is_deeply( [from_contfrac($n)], [$n,1], "from_contfrac($n)" );
  }

  # Unit fractions (1/d): CF is [0, d]
  for my $d (2, 3, 100) {
    is_deeply( [contfrac(1,$d)], [0,$d], "contfrac(1,$d)" );
    is_deeply( [from_contfrac(0,$d)], [1,$d], "from_contfrac(0,$d)" );
  }

  # from_contfrac with no arguments
  is_deeply( [from_contfrac()], [0,1], "from_contfrac() = (0,1)" );

  # from_contfrac with negative first element
  is_deeply( [from_contfrac(-3)], [-3,1], "from_contfrac(-3)" );
  is_deeply( [from_contfrac(-3,7)], [-20,7], "from_contfrac(-3,7)" );

  # 0/d reduces to 0/1
  is_deeply( [contfrac(0,100)], [0], "contfrac(0,100)" );
  is_deeply( [from_contfrac(0)], [0,1], "from_contfrac(0) = (0,1)" );
};

subtest 'bigint' => sub {
  my @data = (
    ["1180591620717411303425","590295810358705651712",
     [2,"590295810358705651712"]],
  );
  for my $t (@data) {
    my($n,$d,$exp) = @$t;
    is_deeply( [map{"$_"}contfrac($n,$d)], [map{"$_"}@$exp], "contfrac($n,$d) bigint" );
    is_deeply( [map{"$_"}from_contfrac(@$exp)], [$n,$d], "from_contfrac bigint" );
  }
};
