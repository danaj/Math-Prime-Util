#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use v5.16;
$|=1;

my @cmap = (0..9, 'a'..'z', 'A'..'Z','!','@','#','$','%');
while (1) {
  my $a = urandomm(1e13);
  my $l = 1 + urandomm(50000);
  my $k = 1 + urandomm(64);
  my $b = $a + $l - 1;

  #$b=~0;  $a=$b-$l;

  my(@a1, @a2);

  print "$cmap[$k]";
  #@a1 = @{Math::Prime::Util::almost_prime_sieve($k, $a, $b)};
  foralmostprimes { push @a1,$_; } $k, $a, $b;
  for ($a .. $b) { push @a2, $_ if prime_bigomega($_)==$k; }
  #die "k $k  beg $a  end $b" unless areq(\@a1,\@a2);
  #die "\nforalmostprimes { say } $k, $a, $b;\n" unless areq(\@a1,\@a2);
  die "\nforalmostprimes { say } $k, $a, $b;\n" unless vecequal(\@a1,\@a2);
}


# We have vecequal now
sub areq {
  my($a,$b) = @_;
  return 1 if !defined $a && !defined $b;
  return 0 unless defined $a && defined $b;
  return 0 unless scalar(@$a) == scalar(@$b);
  for my $i (0 .. $#$a) {
    next if !defined $a->[$i] && !defined $b->[$i];
    return 0 if !defined $a->[$i] || !defined $b->[$i];
    return 0 if $a->[$i] != $b->[$i];
  }
  1;
}

# mmpu 'for my $k (1..12) { for my $l (0 .. 100) { for my $s (0 .. 1000) { say "---- $k $s $l"; foralmostprimes { say } $k,$s,$s+$l; } } }' | shasum
# mmpu 'for my $k (1..61) { say "$k:"; foralmostprimes { say } $k,nth_almost_prime($k,10),nth_almost_prime($k,11); }' | shasum
