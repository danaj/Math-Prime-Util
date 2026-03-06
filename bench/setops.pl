#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/vecequal toset set_is_subset setcontains
                         setintersect setunion setminus setdelta/;
use Math::Prime::Util::PP;   # Have it available for comparison
use Benchmark qw/:all/;


my $N = shift || 100000;

sub _hash32 {
  use integer;
  my $x = shift;
  $x = (($x >> 16) ^ $x) * 0x45d9f3b;
  $x = (($x >> 16) ^ $x) * 0x45d9f3b;
  $x = ($x >> 16) ^ $x;
  return $x & 0xFFFFFFFF;
}

my @set4 = map {_hash32($_)} 0..$N-1;       # random1
my @set5 = map {_hash32(10*$N+$_)} 0..$N-1; # random2
my @iset4 = @{toset(@set4)};   print "iset4 entries:  ",$#iset4+1,"\n";
my @iset5 = @{toset(@set5)};   print "iset5 entries:  ",$#iset5+1,"\n";

use Set::Tiny;
my $st4 = Set::Tiny->new(@iset4);
my $st5 = Set::Tiny->new(@iset5);

my $R;

my $ts=0;
cmpthese(-1, {
  "is_subset" => sub { $ts += set_is_subset(\@iset4, \@iset5); },
  "contains"  => sub { $ts += setcontains(\@iset4, \@iset5); },
}) if (0);

cmpthese(-1, {
  "intersect list" => sub { $R=setintersect(\@set4,\@set5); },
  "intersect iset" => sub { $R=setintersect(\@iset4,\@iset5); },
  "intersect Set::Tiny" => sub { $R=$st4->intersection($st5);},
  #"MPUPP iset"=>sub { $R=Math::Prime::Util::PP::setintersect(\@iset4,\@iset5); },
}) if (1);

cmpthese(-1, {
  "union list"      => sub { $R=setunion(\@set4,\@set5); },
  "union iset"      => sub { $R=setunion(\@iset4,\@iset5); },
  "union Set::Tiny" => sub { $R=$st4->union($st5);},
}) if (1);

cmpthese(-1, {
  "minus list"       => sub { $R=setminus(\@set4,\@set5); },
  "minus iset"       => sub { $R=setminus(\@iset4,\@iset5); },
  "minus Set::Tiny"  => sub { $R=$st4->difference($st5);},
}) if (1);

cmpthese(-1, {
  "delta list"       => sub { $R=setdelta(\@set4,\@set5); },
  "delta iset"       => sub { $R=setdelta(\@iset4,\@iset5); },
  "delta Set::Tiny"  => sub { $R=$st4->symmetric_difference($st5);},
  #"MPUPP iset"=>sub { $R=Math::Prime::Util::PP::setdelta(\@iset4,\@iset5); },
}) if (1);
