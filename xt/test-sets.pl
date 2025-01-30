#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/vecequal toset setintersect setunion setminus setdelta/;
use Benchmark qw/:all/;


my $N = 100000;

my @set1 = map { 2*$_+1 } 0..$N-1;          # Odds
my @set2 = map { 2*$_   } 0..$N-1;          # Evens
my @set3 = reverse @set1;                   # Odds descending
my @set4 = map { _hash32($_) } 0..$N-1;     # random1
my @set5 = map {_hash32(10*$N+$_)} 0..$N-1; # random2


sub _hash32 {
  use integer;
  my $x = shift;
  $x = (($x >> 16) ^ $x) * 0x45d9f3b;
  $x = (($x >> 16) ^ $x) * 0x45d9f3b;
  $x = ($x >> 16) ^ $x;
  return $x & 0xFFFFFFFF;
}

# Verify toset

{
  my @s = toset(\@set3);
  die "toset on descending odds didn't work right" unless vecequal(\@set1,\@s);
}

my @iset4 = toset(\@set4);   print "iset4 entries:  ",$#iset4+1,"\n";
my @iset5 = toset(\@set5);   print "iset5 entries:  ",$#iset5+1,"\n";

# Verify functions

if (1) {
  my @s1 = setintersect(\@set1,\@set2);
  die "intersect of odds and evens should be null" unless @s1 = 0;
  my @s2 = setintersect(\@set1,\@set3);
  die "intersect of odds and odds should be odds" unless vecequal(\@set1,\@s2);
  my @s3 = setintersect(\@set2,\@set3);
  die "intersect of evens and odds should be null" unless @s3 = 0;
}
if (1) {
  my @s1 = setunion(\@set1,\@set2);
  die "union of odds and evens should be all" unless vecequal(\@s1,[0..2*$N-1]);
  my @s2 = setunion(\@set1,\@set3);
  die "union of odds and odds should be odds" unless vecequal(\@set1,\@s2);
  my @s3 = setunion(\@set2,\@set3);
  die "union of evens and odds should be all" unless vecequal(\@s3,[0..2*$N-1]);
}
if (1) {
  my @s1 = setminus(\@set1,\@set2);
  die "diff of odds and evens should be odds" unless vecequal(\@s1,\@set1);
  my @s2 = setminus(\@set1,\@set3);
  die "diff of odds and odds should be null" unless @s2 = 0;
  my @s3 = setminus(\@set2,\@set3);
  die "diff of evens and odds should be evens" unless vecequal(\@s3,\@set2);
}
if (1) {
  my @s1 = setdelta(\@set1,\@set2);
  die "delta of odds and evens should be all" unless vecequal(\@s1,[0..2*$N-1]);
  my @s2 = setdelta(\@set1,\@set3);
  die "delta of odds and odds should be null" unless @s2 = 0;
  my @s3 = setdelta(\@set2,\@set3);
  die "delta of evens and odds should be all" unless vecequal(\@s3,[0..2*$N-1]);
}

my @R;

use Array::Set qw/set_intersect set_union set_diff set_symdiff/;

use Set::SortedArray;
my $saset4 = Set::SortedArray->new(@iset4);
my $saset5 = Set::SortedArray->new(@iset5);
my $R;

use Set::IntSpan::Fast;
my $sisf4 = Set::IntSpan::Fast->new(@iset4);
my $sisf5 = Set::IntSpan::Fast->new(@iset5);

use Set::Functional;
my @sf4 = Set::Functional::setify(@iset4);
my @sf5 = Set::Functional::setify(@iset5);


cmpthese(-1, {
  "intersect odds/odds"   => sub { @R=setintersect(\@set1,\@set1); },
  "intersect odds/evens"  => sub { @R=setintersect(\@set1,\@set2); },
  "intersect rodds/evens" => sub { @R=setintersect(\@set2,\@set3); },
  "intersect hashes"      => sub { @R=setintersect(\@set4,\@set5); },
  "intersect iset hashes" => sub { @R=setintersect(\@iset4,\@iset5); },
  #"Set::SortedArray iset" => sub { $R=$saset4->intersection($saset5); },
  #"Array::Set iset"       => sub { $R=set_intersect(\@iset4,\@iset5); },
  #"Set::IntSpan::Fast iset"=>sub { $R=$sisf4->intersection($sisf5); },
  "Set::Functional"       => sub {@R=Set::Functional::intersection(\@sf4,\@sf5);},
  #"inter1 iset hashes" => sub { @R=inter1(\@iset4,\@iset5); },
  #"inter2 iset hashes" => sub { @R=inter2(\@iset4,\@iset5); },
  #"inter3 iset hashes" => sub { @R=inter3(\@iset4,\@iset5); },
});
exit(0);

cmpthese(-1, {
  "union odds/odds"   => sub { @R=setunion(\@set1,\@set1); },
  "union odds/evens"  => sub { @R=setunion(\@set1,\@set2); },
  "union rodds/evens" => sub { @R=setunion(\@set2,\@set3); },
  "union hashes"      => sub { @R=setunion(\@set4,\@set5); },
  "union iset hashes" => sub { @R=setunion(\@iset4,\@iset5); },
  #"Set::SortedArray iset" => sub { $R=$saset4->union($saset5); },
  #"Array::Set iset"       => sub { $R=set_union(\@iset4,\@iset5); },
  #"Set::IntSpan::Fast iset"=>sub { $R=$sisf4->union($sisf5); },
});

cmpthese(-1, {
  "minus odds/odds"   => sub { @R=setminus(\@set1,\@set1); },
  "minus odds/evens"  => sub { @R=setminus(\@set1,\@set2); },
  "minus rodds/evens" => sub { @R=setminus(\@set2,\@set3); },
  "minus hashes"      => sub { @R=setminus(\@set4,\@set5); },
  "minus iset hashes" => sub { @R=setminus(\@iset4,\@iset5); },
  #"Set::SortedArray iset" => sub { $R=$saset4->difference($saset5); },
  #"Array::Set iset"       => sub { $R=set_diff(\@iset4,\@iset5); },
  #"Set::IntSpan::Fast iset"=>sub { $R=$sisf4->diff($sisf5); },
});

cmpthese(-1, {
  "delta odds/odds"   => sub { @R=setdelta(\@set1,\@set1); },
  "delta odds/evens"  => sub { @R=setdelta(\@set1,\@set2); },
  "delta rodds/evens" => sub { @R=setdelta(\@set2,\@set3); },
  "delta hashes"      => sub { @R=setdelta(\@set4,\@set5); },
  "delta iset hashes" => sub { @R=setdelta(\@iset4,\@iset5); },
  "Set::SortedArray iset" => sub { $R=$saset4->symmetric_difference($saset5);},
  "Array::Set iset"       => sub { $R=set_symdiff(\@iset4,\@iset5); },
  "Set::IntSpan::Fast iset"=>sub { $R=$sisf4->diff($sisf5); },
});


sub inter1 {
  my($a1,$a2) = @_;
  my %counts;
  ++$counts{$_} for @$a1;
  my @c = grep { --$counts{$_} >= 0 } @$a2;
  return @c;
}
sub inter2 {
  my($a1,$a2) = @_;
  my %count= ();
  foreach my $element (@$a1, @$a2) {
    $count{$element}++;
  }
  my @c;
  foreach my $element (keys %count) {
    push @c, $element if $count{$element} > 1;
  }
  return @c;
}
sub inter3 {
  my($ra,$rb) = @_;
  #croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY';
  #croak 'Not an array reference' unless (ref($rb) || '') eq 'ARRAY';
  ($ra,$rb) = ($rb,$ra) if scalar(@$ra) > scalar(@$rb);  # Performance
  my(%ina,%seen,$k);
  $ina{$_}=undef for @$ra;
  my @set = grep { exists $ina{$_} && not $seen{$k=$_}++ } @$rb;
  #for (@set) { return vecsort(@set) if !ref($_) && ($_ >= INTMAX || $_ <= INTMIN); }
  @set = sort { $a<=>$b } @set;
  return @set;
}
