#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use threads::shared;
use Math::Prime::Util qw/is_prime is_strong_pseudoprime forcomposites/;
my $nthreads = 4;

# Single base.

my @composites;
forcomposites { push @composites, $_ if $_ % 2; } 1_000_000;

# Serial:
# my $base = 2;
# my $maxn = 2;
# while (1) {
#   for my $n (@composites) {
#     if (is_strong_pseudoprime($n,$base)) {
#       if ($n > $maxn) {
#         print "base $base good up to $n\n";
#         $maxn = $n;
#       }
#       last;
#     }
#   }
#   $base++;
# }

# Parallel:
my $maxn :shared;
my $start = int(2**60+2**41);  # People have mined below 2^55
$maxn = 2047;
my @threads;
push @threads, threads->create('search_bases', $start, $_) for 1..$nthreads;
# We should sit here doing cond_waits on a results array.
$_->join() for (@threads);

sub search_bases {
  my($start, $t) = @_;
  for (my $base = $start + $t - 1; 1; $base += $t) {
    next if is_strong_pseudoprime(4, $base) || is_strong_pseudoprime(6, $base);
    for my $n (@composites) {
      if (is_strong_pseudoprime($n,$base)) {
        if ($n > $maxn) {
          lock($maxn);
          print "base $base good up to $n\n" if $n > $maxn;
          $maxn = $n;
        }
        last;
      }
    }
  }
}

__END__

base 2 good up to 2047
base 3273 good up to 2209
base 4414 good up to 2443
base 5222 good up to 2611
base 8286 good up to 4033
base 10822 good up to 5411
base 13011 good up to 6505
base 67910 good up to 9073
base 82967 good up to 10371
base 254923 good up to 18299
base 2974927 good up to 18721
base 4095086 good up to 38323
base 70903283 good up to 38503

(best results known, not found with this program)
2011-02-12  base 814494960528 good up to 132239
2012-07-02  base 64390572806844 good up to 161701
2012-10-15  base 1769236083487960 good up to 192001
2012-10-17  base 1948244569546278 good up to 212321
2013-01-14  base 34933608779780163 good up to 218245
2013-03-03  base 9345883071009581737 good up to 341531
