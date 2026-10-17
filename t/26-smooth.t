#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_smooth is_rough smooth_count rough_count
                         has_factors_above has_factors_below
                         has_no_factors_above has_no_factors_below
                         factor vecnone/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

plan tests => 4          # small is_smooth / is_rough
            + 4+3        # special case is_smooth / is_rough
            + 4+4        # is_smooth large
            + 4+4        # smooth_count / rough_count
            + 3;         # factor predicates

###### is_smooth / is_rough

{
  my(@sexp, @sgot, @sngot,  @rexp, @rgot, @rngot);
  for my $k (0..9,11,13,17,29) {
    for my $n (0..12, 143, 187, 253, 319, 341, 851, 1073, 1147) {
      push @sexp, fac_is_smooth($n, $k);
      push @rexp, fac_is_rough($n, $k);
      push @sgot, is_smooth($n, $k);
      push @rgot, is_rough($n, $k);
      push @sngot, is_smooth(-$n, $k);
      push @rngot, is_rough(-$n, $k);
    }
  }
  is_deeply( \@sgot, \@sexp, "is_smooth(n,k) for small inputs" );
  is_deeply( \@sngot, \@sexp, "is_smooth(-n,k) for small inputs" );
  is_deeply( \@rgot, \@rexp, "is_rough(n,k) for small inputs" );
  is_deeply( \@rngot, \@rexp, "is_rough(-n,k) for small inputs" );
}

is(is_smooth(1000000,10000),1,"1000000 is 10000-smooth");
is(is_smooth(1000127,10000),0,"1000127 is not 10000-smooth");
is(is_rough(1000127,3000),0,"1000127 is not 3000-rough");
is(is_rough(1000157,3000),0,"1000157 is not 3000-rough");

is(is_rough("137438953481",3000),1,"137438953481 is 3000-rough");
is(is_rough("137438953493",3000),0,"137438953493 is not 3000-rough");
is(is_rough("137438953529",3000),1,"137438953529 is 3000-rough");

{
  my $n = "1377276413364943226363244108454842276965894752197358387200000"; # 97
  is( is_smooth($n, 23), 0, "large 97-smooth number" );
  is( is_smooth($n, 96), 0, "large 97-smooth number" );
  is( is_smooth($n, 97), 1, "large 97-smooth number" );
  is( is_smooth($n, 98), 1, "large 97-smooth number" );
}
{
  my $n = "172864518041328651521584134678230948270774322090771071422829"; # 2081
  is( is_smooth($n, 4073), 1, "large 4073-smooth, 2081-rough number" );
  is( is_rough($n, 2080), 1, "large 4073-smooth, 2081-rough number" );
  is( is_rough($n, 2081), 1, "large 4073-smooth, 2081-rough number" );
  is( is_rough($n, 2082), 0, "large 4073-smooth, 2081-rough number" );
}

subtest 'factor predicates and boundary conventions' => sub {
  # n, k, has_above, has_below, has_no_above, has_no_below
  my @cases = (
    [  0, 0, 0, 0, 1, 1 ],
    [  0, 1, 0, 1, 1, 0 ],
    [  0, 8, 0, 1, 1, 0 ],
    [  1, 0, 0, 0, 1, 1 ],
    [  1, 1, 0, 0, 1, 1 ],
    [  1, 8, 0, 0, 1, 1 ],
    [  2, 0, 1, 0, 0, 1 ],
    [  2, 1, 1, 0, 0, 1 ],
    [  2, 2, 0, 0, 1, 1 ],
    [  2, 3, 0, 1, 1, 0 ],
    [ 49, 7, 0, 0, 1, 1 ],
    [ 77, 7, 1, 0, 0, 1 ],
    [ 77, 8, 1, 1, 0, 0 ],
    [121, 8, 1, 0, 0, 1 ],
    [210, 8, 0, 1, 1, 0 ],
  );
  for my $case (@cases) {
    my($n, $k, @exp) = @$case;
    for my $signed_n ($n ? ($n, -$n) : (0)) {
      is_deeply([has_factors_above($signed_n, $k),
                 has_factors_below($signed_n, $k),
                 has_no_factors_above($signed_n, $k),
                 has_no_factors_below($signed_n, $k)],
                \@exp, "factor predicates($signed_n,$k)");
    }
  }
};

subtest 'factor predicates with large integers' => sub {
  require Math::BigInt;
  my $smooth = '1377276413364943226363244108454842276965894752197358387200000';
  my $rough = '172864518041328651521584134678230948270774322090771071422829';
  my $bigk = '18446744073709551616';
  my @cases = (
    [$smooth,   96, 1, 1, 0, 0],
    [$smooth,   97, 0, 1, 1, 0],
    [$rough,  2081, 1, 0, 0, 1],
    [$rough,  2082, 1, 1, 0, 0],
    [     0, $bigk, 0, 1, 1, 0],
    [     1, $bigk, 0, 0, 1, 1],
    [    12, $bigk, 0, 1, 1, 0],
  );
  for my $case (@cases) {
    my($n, $k, @exp) = @$case;
    for my $signed_n ($n ? ($n, "-$n") : (0)) {
      for my $args ([$signed_n, $k],
                    [Math::BigInt->new($signed_n), Math::BigInt->new($k)]) {
        is_deeply([has_factors_above($args->[0], $args->[1]),
                   has_factors_below($args->[0], $args->[1]),
                   has_no_factors_above($args->[0], $args->[1]),
                   has_no_factors_below($args->[0], $args->[1])],
                  \@exp, 'factor predicates with large strings or objects');
      }
    }
  }
};

subtest 'factor predicate argument validation' => sub {
  for my $name (qw/has_factors_above has_factors_below
                   has_no_factors_above has_no_factors_below/) {
    my $func = Math::Prime::Util->can($name);
    for my $args ([undef, 2], ['1.5', 2], ['invalid', 2],
                  [2, undef], [2, -1], [2, '1.5']) {
      my $ok = eval { $func->(@$args); 1 };
      ok(!$ok && $@, "$name rejects invalid arguments");
    }
  }
};

###### smooth_count
{
  # mpu 'for $n (0..5){for $k (0..5){push @v,vecsum(map{is_smooth($_,$k)}1..$n)}} say join ",",@v;'
  my @exp = (0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,1,1,2,3,3,3,1,1,3,4,4,4,1,1,3,4,4,5);
  my @got;
  for my $n (0..5) { for my $k (0..5) { push @got, smooth_count($n,$k); } }
  is_deeply( \@got, \@exp, "smooth_count(0..5, 0..5)" );
}
is(smooth_count(100,17), 67, "smooth_count(100,17)");
is(smooth_count(1980627498,9), 5832, "smooth_count(1980627498,9)");
SKIP: {
  skip "skipping slow smooth count test with PP", 1 unless $usexs || $extra;
  is(smooth_count(10000000,400), 1132424, "smooth_count(10000000,400)");
}

###### rough_count
{
  # mpu 'for $n (0..5){for $k (0..5){push @v,vecsum(map{is_rough($_,$k)}1..$n)}} say join ",",@v;'
  my @exp = (0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,1,1,1,3,3,3,2,1,1,4,4,4,2,1,1,5,5,5,3,2,2);
  my @got;
  for my $n (0..5) { for my $k (0..5) { push @got, rough_count($n,$k); } }
  is_deeply( \@got, \@exp, "rough_count(0..5, 0..5)" );
}
is(rough_count(3700621409,15), 709809501, "rough_count(3700621409,15)");
is("".rough_count("10000000000000000000000000000",13), "2077922077922077922077922077",
   "rough_count(10^28,13) bigint result");
is("".rough_count("10000000000000000000000000000",19), "1805253569959452312393488861",
   "rough_count(10^28,19) bigint result");


###### ---- helper functions ----

sub fac_is_smooth {
  my($n, $k) = @_;
  # True if no prime factors of n are larger than k
  return 0+(vecnone { $_ > $k } factor($n));
}

sub fac_is_rough {
  my($n, $k) = @_;
  # True if no prime factors of n are smaller than k
  return 0+(vecnone { $_ < $k } factor($n));
}
