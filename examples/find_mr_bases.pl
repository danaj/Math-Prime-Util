#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use threads::shared;
use Math::Prime::Util qw/is_prime is_strong_pseudoprime/;
my $nthreads = 12;

# Single base.

my @composites;
for my $n (3 .. 1000000) {
  push @composites, $n if $n % 2 && !is_prime($n);
}

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
my $start = int(2**59+2**41);  # People have mined below 2^55
$maxn = 2047;
my $nextn = 2049;
my @threads;
push @threads, threads->create('search_bases', $start, $_) for (0..$nthreads-1);
# We should sit here doing cond_waits on a results array.
$_->join() for (@threads);

sub search_bases {
  my($start, $t) = @_;
  my $base = $start + $t;
  while (1) {
    do { $base += $t; next; } if is_strong_pseudoprime($nextn, $base);
    for my $n (@composites) {
      if (is_strong_pseudoprime($n,$base)) {
        if ($n > $maxn) {
          lock($maxn);
          print "base $base good up to $n\n";
          $maxn = $n;
          $nextn = $n+2;  $nextn++ while is_prime($nextn);
        }
        last;
      }
    }
    $base += $t;
  }
}

__END__

base 2 good up to 2047
base 1320 good up to 4097
base 4712 good up to 4711
base 5628 good up to 5627
base 7252 good up to 7251
base 7852 good up to 7851
base 14787 good up to 9409
base 17340 good up to 10261
base 61380 good up to 11359
base 78750 good up to 13747
base 254923 good up to 18299
base 486605 good up to 25761
base 1804842 good up to 32761
base 4095086 good up to 38323
base 12772344 good up to 40501
base 42162995 good up to 97921

(best results known, not found with this program)
2011-02-12  base 814494960528 good up to 132239
2012-07-02  base 64390572806844 good up to 161701
2012-10-15  base 1769236083487960 good up to 192001
2012-10-17  base 1948244569546278 good up to 212321
2013-01-14  base 34933608779780163 good up to 218245
