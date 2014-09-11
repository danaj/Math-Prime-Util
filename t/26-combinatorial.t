#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/binomial factorial forcomb forperm/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

use Math::BigInt try => "GMP,Pari";

my %perms = (
 0 => [[]],
 1 => [[0]],
 2 => [[0,1],[1,0]],
 3 => [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]],
 4 => [[0,1,2,3],[0,1,3,2],[0,2,1,3],[0,2,3,1],[0,3,1,2],[0,3,2,1],[1,0,2,3],[1,0,3,2],[1,2,0,3],[1,2,3,0],[1,3,0,2],[1,3,2,0],[2,0,1,3],[2,0,3,1],[2,1,0,3],[2,1,3,0],[2,3,0,1],[2,3,1,0],[3,0,1,2],[3,0,2,1],[3,1,0,2],[3,1,2,0],[3,2,0,1],[3,2,1,0]],
);

# TODO: Add a bunch of combs here:  "5,3" => [[..],[..],[..]],

plan tests => 1                        # Factorial
            + 6 + 4                    # Combinations
            + scalar(keys(%perms)) + 1 # Permutations
            ;

sub fact { my $n = Math::BigInt->new("$_[0]"); $n->bfac; }
{
  my @result = map { factorial($_) } 0 .. 100;
  my @expect = map { fact($_) } 0 .. 100;
  is_deeply( \@result, \@expect, "Factorials 0 to 100" );
}


{ my @p = (); forcomb { push @p, [@_] } 0;
  is_deeply( [@p], [[]], "forcomb 0" ); }
{ my @p = (); forcomb { push @p, [@_] } 1;
  is_deeply( [@p], [[0]], "forcomb 1" ); }
{ my @p = (); forcomb { push @p, [@_] } 0,0;
  is_deeply( [@p], [[]], "forcomb 0,0" ); }
{ my @p = (); forcomb { push @p, [@_] } 5,0;
  is_deeply( [@p], [[]], "forcomb 5,0" ); }
{ my @p = (); forcomb { push @p, [@_] } 5,6;
  is_deeply( [@p], [], "forcomb 5,6" ); }
{ my @p = (); forcomb { push @p, [@_] } 5,5;
  is_deeply( [@p], [[0,1,2,3,4]], "forcomb 5,5" ); }

{ my @data = (qw/apple bread curry/);
  my @p = (); forcomb { push @p, [@data[@_]] } @data,2;
  my @e = ([qw/apple bread/],[qw/apple curry/],[qw/bread curry/]);
  is_deeply( \@p,\@e, "forcomb 3,2" ); }
{ my @data = (qw/ant bee cat dog/);
  my @p = (); forcomb { push @p, [@data[@_]] } @data,3;
  my @e = ([qw/ant bee cat/],[qw/ant bee dog/],[qw/ant cat dog/],[qw/bee cat dog/]);
  is_deeply( \@p,\@e, "forcomb 4,3" ); }

{ my $b = binomial(20,15);
  my $s = 0; forcomb { $s++ } 20,15;
  is($b, 15504, "binomial(20,15) is 15504");
  is($s, $b, "forcomb 20,15 yields binomial(20,15) combinations"); }


while (my($n, $expect) = each (%perms)) {
  my @p = (); forperm { push @p, [@_] } $n;
  is_deeply( \@p, $expect, "forperm $n" );
}

{ my $s = 0; forperm { $s++ } 7;
  is($s, factorial(7), "forperm 7 yields factorial(7) permutations"); }
