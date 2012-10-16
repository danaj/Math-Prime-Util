#!/usr/bin/env perl
use strict;
use warnings;
use Math::BigInt;
use Getopt::Long;
use Math::Prime::Util qw/primes next_prime is_prime/;
$| = 1;

my $show_safe = 0;
my $show_sophie = 0;
my $show_twin = 0;
my $show_usage = 0;
my $segment_size = 30 * 128_000;   # 128kB

GetOptions( "safe"   => \$show_safe,
            "sophie" => \$show_sophie,
            "twin"   => \$show_twin,
            "help"   => \$show_usage,
          ) || die_usage();
die_usage() if $show_usage;

# Get the start and end values.  Verify they're positive integers.
die_usage() unless @ARGV == 2;
my ($start, $end) = @ARGV;
die "$start isn't a positive integer" if $start =~ tr/0123456789//c;
die "$end isn't a positive integer" if $end =~ tr/0123456789//c;

# Turn start and end into bigints if they're very large
if ( ($start >= ~0 && $start ne ~0) || ($end >= ~0 && $end ne ~0) ) {
  $start = Math::BigInt->new($start);
  $end = Math::BigInt->new($end);
}

while ($start <= $end) {

    # small segment size if we're doing bigints
    $segment_size = 10000 if $start > ~0;

    my $seg_start = $start;
    my $seg_end = $start + $segment_size;
    $seg_end = $end if $end < $seg_end;
    $start = $seg_end+1;

    # Get a list of all primes in the segment.
    my @p = @{primes($seg_start, $seg_end)};
    #warn "-- ", scalar @p, " primes between $start and $segment_end\n";
    next unless scalar @p > 0;


    # Restrict to twin primes if requested.
    if ($show_twin) {
      # Add a last element so we can look at it.  Simple: push next_prime.
      push @p, is_prime($p[-1]+2) ? $p[-1]+2 : 0;
      # Trivial but slow:
      #@p = grep { is_prime( $_+2 ); } @p;
      # Loop over @p looking for values with prime == next+2
      my @twin;
      my $prime = shift @p;
      foreach my $next (@p) {
        push @twin, $prime if $prime+2 == $next;
        $prime = $next;
      }
      @p = @twin;
      #warn "   reduced to ", scalar @p, " twin primes\n";
    }


    # Restrict to safe primes if requested.
    if ($show_safe) {
      @p = grep { is_prime( ($_-1) >> 1 ); } @p;
      #warn "   reduced to ", scalar @p, " safe primes\n";
    }


    # Restrict to Sophie Germain primes if requested.
    if ($show_sophie) {
      @p = grep { is_prime( 2*$_+1 ); } @p;
      #warn "   reduced to ", scalar @p, " SG primes\n";
    }

    # print this segment
    print join("\n", @p), "\n"  if scalar @p > 0;
}

sub die_usage {
  die <<EOU;
Usage: $0 [options]  START  END

Displays all primes between the positive integers START and END, inclusive.

Options:
  --help      displays this help message
  --safe      displays only safe primes, where (p-1)/2 is also prime
  --twin      displays only twin primes, where p+2 is also prime
  --sophie    displays only Sophie Germain primes, where 2p+1 is also prime

Note that options can be combined, e.g. display all safe twin primes.
EOU
}
