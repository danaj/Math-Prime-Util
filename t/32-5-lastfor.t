#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/forprimes forcomposites foroddcomposites fordivisors
                         forpart forcomp forcomb forperm forderange formultiperm
                         forfactored forsquarefree forsquarefreeint
                         lastfor is_power vecsum/;

subtest 'lastfor in iterator families' => sub {
  my @zn;
  forprimes {
    my $p=$_;
    forprimes {
      lastfor, push @zn,$_ if $_ % $p == 1;
    } 1000;
  } 100;
  is_deeply( \@zn,
             [3,7,11,29,23,53,103,191,47,59,311,149,83,173,283,107,709,367,269,569,293,317,167,179,389],
             "lastfor works in forprimes" );

  @zn = ();
  forprimes {
    my $p=$_;
    forcomposites {
      lastfor, push @zn,$_ if $_ % $p == 1;
    } 1000;
  } 100;
  is_deeply( \@zn,
             [9,4,6,8,12,14,18,20,24,30,32,38,42,44,48,54,60,62,68,72,74,80,84,90,98],
             "lastfor works in forcomposites" );

  @zn = ();
  forprimes {
    my $p=$_;
    foroddcomposites {
      lastfor, push @zn,$_ if $_ % $p == 1;
    } 1000;
  } 100;
  is_deeply( \@zn,
             [9,25,21,15,45,27,35,39,93,117,63,75,165,87,95,213,119,123,135,143,147,159,333,357,195],
             "lastfor works in foroddcomposites" );

  my @powers;
  for my $n (1..20) {
    fordivisors { lastfor,push @powers,$_ if is_power($_) } $n;
  }
  is_deeply( \@powers, [4,4,9,4,4,9,4], "lastfor works in fordivisors" );

  my $firstpart;
  forpart { lastfor,return if @_ < 4; $firstpart++; } 7;
  is($firstpart, 6, "lastfor works in forpart");

  my $firstcomp;
  forcomp { lastfor,return if @_ < 4; $firstcomp++; } 7;
  is($firstcomp, 15, "lastfor works in forcomp");

  my $smallcomb;
  forcomb { lastfor,return if vecsum(@_) > 11; $smallcomb++; } 7,4;
  is($smallcomb, 9, "lastfor works in forcomb");

  my $t;
  forperm { lastfor,return if $_[3]==5; $t++; } 7;
  is($t, 12, "lastfor works in forperm");

  $t = 0;
  forderange { lastfor,return if $_[3]==5; $t++; } 7;
  is($t, 5, "lastfor works in forderange");

  $t = 0;
  formultiperm { lastfor if "miles" eq join("",@_); $t++; } [split(//,"smile")];
  is($t, 81, "lastfor works in formultiperm");

  my @ps;
  forprimes {
    lastfor if $_ >= 7;
    push @ps, $_;
    forcomposites { push @ps,$_; } $_;
  } 20;
  is_deeply( \@ps, [2,3,5,4,7,4,6], "nested lastfor semantics" );

  $t = 0;
  forcomposites { $t=$_; lastfor if $_ > 2000; } 20000;
  is($t, 2001, "lastfor in forcomposites stops appropriately");
};

subtest 'callback die restores lastfor state' => sub {
  ok(!eval { lastfor(); 1 }, 'lastfor initially croaks outside a loop');

  my $err = eval { forprimes { die "intentional callback die\n" } 10; 1 } ? '' : $@;
  like($err, qr/intentional callback die/, 'callback die propagates');

  ok(!eval { lastfor(); 1 }, 'lastfor still croaks after callback die');

  $err = eval {
    forprimes {
      eval { forprimes { die "inner callback die\n" } 10 };
      lastfor;
    } 20;
    1;
  } ? '' : $@;
  is($err, '', 'caught inner callback die does not corrupt outer loop state');
};

subtest 'formultiperm callback die restores lastfor state' => sub {
  ok(!eval { lastfor(); 1 }, 'lastfor initially croaks outside a loop');

  my $err = eval { formultiperm { die "intentional callback die\n" } [1,1,2]; 1 } ? '' : $@;
  like($err, qr/intentional callback die/, 'callback die propagates');

  ok(!eval { lastfor(); 1 }, 'lastfor still croaks after callback die');

  my @out;
  formultiperm { push @out, "@_"; lastfor; } [1,1,2];
  is_deeply(\@out, ['1 1 2'], 'lastfor still works after callback die');
};

subtest 'empty-callback lastfor' => sub {
  my $t = 0;
  forpart { $t++; lastfor; } 0;
  is($t, 1, 'forpart');

  $t = 0;
  forcomp { $t++; lastfor; } 0;
  is($t, 1, 'forcomp');

  $t = 0;
  forcomb { $t++; lastfor; } 0;
  is($t, 1, 'forcomb 0');

  $t = 0;
  forcomb { $t++; lastfor; } 3,0;
  is($t, 1, 'forcomb 3,0');

  $t = 0;
  forperm { $t++; lastfor; } 0;
  is($t, 1, 'forperm');

  $t = 0;
  forderange { $t++; lastfor; } 0;
  is($t, 1, 'forderange');

  $t = 0;
  formultiperm { $t++; lastfor; } [];
  is($t, 1, 'formultiperm');
};

subtest 'initial-callback lastfor' => sub {
  my $t = 0;
  forfactored { $t++; lastfor; } 0,3;
  is($t, 1, 'forfactored');

  $t = 0;
  forsquarefree { $t++; lastfor; } 0,3;
  is($t, 1, 'forsquarefree');

  $t = 0;
  forsquarefreeint { $t++; lastfor; } 0,3;
  is($t, 1, 'forsquarefreeint');
};

done_testing();
