#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/forpart forcomp forcomb forperm forderange
                         formultiperm forsetproduct/;

subtest 'forpart and forcomp option validation' => sub {
  ok(!eval { forpart { 1 } 7, { amin => -1 }; 1 }, 'forpart {amin=>-1} croaks');
  ok(!eval { forpart { 1 } 7, { n => '3.5' }; 1 }, 'forpart {n=>"3.5"} croaks');
  ok(!eval { forpart { 1 } 7, { prime => 'abc' }; 1 }, 'forpart {prime=>"abc"} croaks');
  ok(!eval { forcomp { 1 } 7, { nmax => '3.5' }; 1 }, 'forcomp {nmax=>"3.5"} croaks');
};

subtest 'basic combinatorial iterator counts' => sub {
  my $t = 0; forpart { $t++ } 4; is($t, 5, 'forpart 4');
  $t = 0; forcomp { $t++ } 4; is($t, 8, 'forcomp 4');
  $t = 0; forcomb { $t++ } 5,2; is($t, 10, 'forcomb 5,2');
  $t = 0; forperm { $t++ } 3; is($t, 6, 'forperm 3');
  $t = 0; forderange { $t++ } 4; is($t, 9, 'forderange 4');
  $t = 0; formultiperm { $t++ } [1,1,2]; is($t, 3, 'formultiperm [1,1,2]');
};

subtest 'forsetproduct' => sub {
  ok(!eval { forsetproduct { } 1,2; }, "forsetproduct not array ref errors");

  my(@set,@out);

  @set=();  @out=();forsetproduct {push @out,"@_"}@set;
  is_deeply(\@out, [], 'forsetproduct empty input -> empty output');

  @set=([1..3]);  @out=();forsetproduct {push @out,"@_"}@set;
  is_deeply(\@out, [1..3], 'forsetproduct single list -> single list');

  @set=([1],[2],[3],[4],[5]);  @out=();forsetproduct {push @out,"@_"}@set;
  is_deeply(\@out, ['1 2 3 4 5'], 'forsetproduct five 1-element lists -> single list');

  @set=([1,2],[3,4,5],[]);  @out=();forsetproduct {push @out,"@_"}@set;
  is_deeply(\@out, [], 'forsetproduct any empty list -> empty output');

  @set=([],[1,2],[3,4,5]);  @out=();forsetproduct {push @out,"@_"}@set;
  is_deeply(\@out, [], 'forsetproduct any empty list -> empty output');

  @set=([1,2],[qw/a b c/]);  @out=();forsetproduct {push @out,"@_"}@set;
  is_deeply(\@out, ['1 a','1 b','1 c','2 a','2 b','2 c'], 'forsetproduct simple test');

  @set=([1,2],[qw/a b c/]);  @out=();forsetproduct {push @out,"@_"; $#_=0; }@set;
  is_deeply(\@out, ['1 a','1 b','1 c','2 a','2 b','2 c'], 'forsetproduct modify size of @_ in block');

  @set=([1,2],[qw/a b c/]);  @out=();forsetproduct {push @out,"@_"; @_=(1..10); }@set;
  is_deeply(\@out, ['1 a','1 b','1 c','2 a','2 b','2 c'], 'forsetproduct replace @_ in sub');

  @set=([1],[qw/a b/]);  @out=();forsetproduct { push @out, "@_"; $_[0] = 9; } @set;
  is_deeply([ $set[0], \@out ], [ [1], ['1 a','9 b'] ],
            'forsetproduct scalar mutation does not modify source array');

  {
    my @a = (1,2);
    my @b = ('a','b');
    @out = ();
    forsetproduct {
      push @out, "@_";
      @a = ();
      @b = ();
    } \@a, \@b;
    is_deeply(\@out, ['1 a','1 b','2 a','2 b'],
              'forsetproduct survives source array mutation in block');
  }
};

done_testing();
