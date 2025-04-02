#!/usr/bin/perl
$|=1;
use warnings;
use strict;
use Math::Prime::Util qw/:all/;
use Math::Prime::Util::PP;

csrand(4);

#for my $bits (3..65,70) {
for my $bits (3..63) {
  print "$bits ";
  for my $n (1..800) {
    my(@x,@y);
    my @a = map { urandomb($bits); } 1..$n;
    my @b = map { urandomb($bits); } 1..$n;
    checksetops(\@a,\@b,$bits,$n);

    my $mp = powint(2,$bits-1);
    @a = map { subint($_,$mp) } @a;
    @b = map { subint($_,$mp) } @a;
    checksetops(\@a,\@b,$bits,$n);
  }
}
print "\n";

sub checksetops {
  my($a,$b,$bits,$n) = @_;
  my(@x,@y);

  @x = setunion($a,$b);
  @y = Math::Prime::Util::PP::setunion($a,$b);
  die "wrong for setunion $bits $n [@$a] [@$b]  [@x]  [@y]" unless vecequal(\@x,\@y);
  @x = setintersect($a,$b);
  @y = Math::Prime::Util::PP::setintersect($a,$b);
  die "wrong for setintersect $bits $n" unless vecequal(\@x,\@y);
  @x = setminus($a,$b);
  @y = Math::Prime::Util::PP::setminus($a,$b);
  die "wrong for setminus $bits $n [@$a] [@$b]  [@x] [@y]" unless vecequal(\@x,\@y);
  @x = setdelta($a,$b);
  @y = Math::Prime::Util::PP::setdelta($a,$b);
  die "wrong for setdelta $bits $n [@$a] [@$b]  [@x]  [@y]" unless vecequal(\@x,\@y);

  # First set must be in set form for insert, contains, remove, invert
  $a = [toset($a)];
  # Second list must not have duplicates for set_is_subset
  $b = [vecuniq(@$b)];
  my $s1 = Math::Prime::Util::set_is_subset($a,$b);
  my $s2 = Math::Prime::Util::setcontains($a,$b);
  my $s3 = Math::Prime::Util::PP::set_is_subset($a,$b);
  my $s4 = Math::Prime::Util::PP::setcontains($a,$b);
  die "wrong for contains $bits $n [@$a] [@$b]" unless $s2 == $s1;
  die "wrong for PP subset $bits $n [@$a] [@$b]" unless $s3 == $s1;
  die "wrong for PP contains $bits $n" unless $s4 == $s1;
}
