#!/usr/bin/perl
use warnings;
use strict;
use v5.16;
use ntheory ":all";
use Math::Prime::Util::PP;
#*rootmod = \&Math::Prime::Util::rootmod;

my $lpr = primes(1e10,1e10+1000);

test("primes", 2, $lpr);
test("primes", 3, $lpr);
test("primes", 5, $lpr);
test("primes", 4, $lpr);
test("primes", 6, $lpr);
test("primes", 9, $lpr);
test("primes",27, $lpr);
test("primes",81, $lpr);

say "";

my $lp2 = [map { $_*$_ } 2,3,5,101,103,107,109];
test("prime squares", 2, $lp2);
test("prime squares", 3, $lp2);
test("prime squares", 5, $lp2);
my $lp3 = [map { powint($_,3) } 2,3,5,101,103,107,109];
test("prime cubes", 2, $lp3);
test("prime cubes", 3, $lp3);
test("prime cubes", 5, $lp3);
my $lpp = [8,16,32,243,625,14641,130321,371293,24137569];
test("prime powers", 2, $lpp);
test("prime powers", 3, $lpp);
test("prime powers", 5, $lpp);

say "";

my $lsf = [grep { is_square_free($_) } 1e10..1e10+100];
test("square free", 2, $lsf);
test("square free", 3, $lsf);
test("square free", 5, $lsf);
test("square free", 4, $lsf);
test("square free", 6, $lsf);
test("square free", 9, $lsf);
test("square free",27, $lsf);
test("square free",81, $lsf);

say "";

my $lcp = [grep { !is_prime($_) && !is_square_free($_) } 1e10 .. 1e10+200];
test("composites", 2, $lcp);
test("composites", 3, $lcp);
test("composites", 5, $lcp);
test("composites", 4, $lcp);
test("composites", 6, $lcp);
test("composites", 9, $lcp);
test("composites",27, $lcp);
test("composites",81, $lcp);

say "";



sub test {
  my($name, $k, $list) = @_;

  my($t,$bad,$und) = (0,0,0);

  for my $p (@$list) {
    my %h;
    my $lim = $p > 1000 ? 1000 : $p-1;
    $h{ powmod($p-$_, $k, $p) }++ for 1..$lim;
    for my $a (keys %h) {
      #my $r = Math::Prime::Util::PP::rootmod($a, $k, $p);
      my $r = rootmod($a, $k, $p);
      if (!defined $r) { $und++; }
      elsif (powmod($r,$k,$p) != $a) { $bad++; }
      $t++;
    }
  }

  printf "Total: %8u  undef: %7u  bad: %5u  for x^%2u = a mod $name\n",$t,$und,$bad,$k;
  #say "Total: $t  undef: $und  bad: $bad  for x^$k = a mod $name";
}
