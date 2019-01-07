#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use v5.16;

# formultiperm is done by doing a require of the Perl code.  Since this is
# inside the forprimes, it does the require inside a call_sv.  This can lead
# to a stack overflow if the required module does too many static variables.
# This is probably a bug in Perl, but regardless this is something we should
# work around.

sub t1 {  my $n = shift;
  forprimes {
    my $ok = 1;
    formultiperm {
      if (!is_prime(fromdigits([@_]))) { $ok = 0; lastfor; }
    } [todigits($_)];
    say if $ok;
  } $n;
}

sub t2 {  my $n = shift;
  forprimes {
    if (!/[^1379]/) {
      my $ok = 1;
      formultiperm {
        if (!is_prime(fromdigits([@_]))) { $ok = 0; lastfor; }
      } [todigits($_)];
      say if $ok;
    }
  } $n;
}

# M. F. Hasler method.
# This runs about 2x faster than Pari/GP
sub t3 { my $digits = shift;
  for my $n (1 .. $digits) {
    my @S;
    my $r = divint(powint(10,$n),9);
    for my $a (1 .. (($n<=1) ? 1 : 9)) {
      for my $b ( (($n>2) ? 1-$a : 0) .. 9-$a ) {
        my $v = mulint($a,$r);
        if ($b == 0) {
          push @S, $v if is_prime($v);
        } else {
          next unless vecall { is_prime(addint($v,mulint($b,powint(10,$_)))) } 0 .. $n-1;
          push @S, map { addint($v,mulint($b,powint(10,$_))) } 0 .. $n-1;
        }
      }
    }
    say for sort { $a<=>$b } @S;
  }
}

t2(10000);
#t3(1031);
