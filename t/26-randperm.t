#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/randperm shuffle vecsample csrand/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

plan tests => 4 + 2;

subtest 'randperm', sub {
  is(@{[randperm(0)]},0,"randperm(0) returns 0 elements");
  is(@{[randperm(1)]},1,"randperm(1) returns 1 element");
  is(@{[randperm(4,0)]},0,"randperm(4,0) returns 0 elements");
  is(@{[randperm(4,1)]},1,"randperm(4,1) returns 1 element");
  is(@{[randperm(4,8)]},4,"randperm(4,8) returns 4 elements");
  is(@{[randperm(100,4)]},4,"randperm(100,4) returns 4 elements");

  my @p128 = randperm(128);
  isnt("@p128", join(" ",0..127), "randperm(128) shuffles");
  is(join(" ",sort {$a<=>$b}@p128), join(" ",0..127), "randperm(128) gives expected indices");

  # sub loop_randperm { my($n, $k, $min) = @_; ... }
  is(join(" / ",loop_randperm(2,2,2)), "0 1 / 1 0", "randperm(2,2) can return all permutations");
  is(join(" / ",loop_randperm(3,2,6)), "0 1 / 0 2 / 1 0 / 1 2 / 2 0 / 2 1", "randperm(3,2) can return all permutations");

  is(scalar loop_randperm(16,undef,10), 10, "randperm(16) can return multiple permutations");
  is(scalar loop_randperm(4,8,6), 6, "randperm(4,8) can return multiple permutations");
  is(scalar loop_randperm(42,1,6), 6, "randperm(42,1) can return multiple permutations");
  is(scalar loop_randperm(1024,2,6), 6, "randperm(1024,2) can return multiple permutations");
  is(scalar loop_randperm(75,6,16), 16, "randperm(75,6) can return multiple permutations");
  is(scalar loop_randperm(30,12,2), 2, "randperm(30,12) can return multiple permutations");
  is(scalar loop_randperm(54321,10,2), 2, "randperm(54321,10) can return multiple permutations");
  is(scalar loop_randperm(123456789,37,2), 2, "randperm(123456789,37) can return multiple permutations");
};


subtest 'shuffle', sub {
  is_deeply([shuffle()], [], "shuffle() = ()");
  is_deeply([shuffle(-277)], [-277], "shuffle(x) = (x)");

  my @d128 = (1..128);
  my @s128 = shuffle(@d128);
  my @t128 = sort { $a<=>$b } @s128;
  is(scalar @s128, scalar @d128, "shuffle n items returns n items");
  isnt("@s128","@d128", "shuffled 128-element array isn't identical");
  is("@t128","@d128", "outputs are the same elements as input");

  {
    my @L = (3,1,2);
    my %seen;
    my $tries = 0;
    # Average tries is 14.7, but could be 100+.
    for (1 .. 1000) {
      $tries++;
      my @S = shuffle(@L);
      undef $seen{"@S"};
      last if 6 <= scalar keys %seen;
    }
    my $fs = join " / ", sort keys %seen;
    is($fs, "1 2 3 / 1 3 2 / 2 1 3 / 2 3 1 / 3 1 2 / 3 2 1",
       "shuffle(a,b,c) selected each permutation at least once ($tries tries)");
  }
};


subtest 'vecsample', sub {
  is_deeply([map {[vecsample($_)]} 0,1,1000], [[],[],[]], "vecsample(k) = ()");
  is_deeply([map {[vecsample($_,())]} 0,1,999],[[],[],[]], "vecsample(k,()) = ()");
  is_deeply([map {[vecsample($_,[])]} 0,1,999],[[],[],[]], "vecsample(k,[]) = ()");
  is_deeply([map {[vecsample(1,($_))]} 0,1,999],[[0],[1],[999]], "vecsample(1,(n)) = (n)");
  is_deeply([map {[vecsample(1,[$_])]} 0,1,999],[[0],[1],[999]], "vecsample(1,(n)) = (n)");
  is(scalar @{[vecsample(2,[80..100])]}, 2, "returns k items with a large list");
  is(scalar @{[vecsample(45,[8,9,10])]}, 3, "returns all items with large k");
  is(scalar @{[vecsample(4,[8..11])]}, 4, "returns all items with exact k");
  is_deeply([sort {$a<=>$b} vecsample(5,[177,888,15,4,-2])],
            [-2,4,15,177,888], "returns all items");
  ok(is_one_of(vecsample(1,500..600), 500..600), "vecsample(1,L) returns something from L");
  {
    my $L = [87,388,657,890];
    my %seen;
    my $tries = 0;
    # Average tries is 3.8, but could be 20+.
    for (1 .. 300) {
      $tries++;
      my @S = vecsample(2,$L);
      undef @seen{@S};
      last if 4 <= scalar keys %seen;
    }
    my @S = sort { $a<=>$b } keys %seen;
    is_deeply(\@S, $L, "vecsample(2,a,b,c,d) selected each value at least once ($tries tries)");
  }
  {
    my @L = (1588..1620);
    vecsample(5,@L);
    is_deeply(\@L, [1588..1620], "Input list is not modified");
    vecsample(5,\@L);
    is_deeply(\@L, [1588..1620], "Input aref is not modified");
  }
};

subtest 'using csrand', sub {
  my @x;
  for (1..4) {
    csrand(14);
    push @x, join " ",shuffle(1..128);
  }
  ok($x[0] eq $x[1] && $x[0] eq $x[2] && $x[0] eq $x[3],
     "shuffles are repeatable with csrand");
  # Entropy seed the RNG for further tests
  csrand();
};


# Destruction test, from List::MoreUtils
{
  my $dcount;
  sub TObjSample::DESTROY { $dcount++; }
  my @ret = vecsample(3, map { bless [], "TObjSample" } 1..10);
  is($dcount, 7, "vecsample unselected items destroyed");
  @ret=();
  is($dcount, 10, "vecsample all items destroyed");
}



sub loop_randperm {
  my($n, $k, $min) = @_;
  my $explen = !defined $k ? $n : $k > $n ? $n : $k;

  my %seen;
  my $tries = 0;
  for (1..10000) {
    $tries++;
    my @S = defined $k ? randperm($n,$k) : randperm($n);
    return () if scalar @S != $explen;
    #diag "??? expected $explen but got ",scalar(@S),"" if scalar @S != $explen;
    undef $seen{"@S"};
    last if $min <= scalar keys %seen;
  }
  diag "loop_randperm($n,$k,$min) took $tries tries" if $tries > 60;
  my @S = sort keys %seen;
  @S;
}

# veccontains should do this, except for the undef part
sub is_one_of {
  my($n, @list) = @_;
  if (!defined $n) {
    for (@list) { return 1 if !defined $_; }
    return 0;
  }
  for (@list) { return 1 if defined $_ && $n eq $_; }
  0;
}
