#!/usr/bin/env perl
use strict;
use warnings;
use Math::BigInt;
use Getopt::Long;
use Math::Prime::Util qw/primes next_prime is_prime/;
$| = 1;

# For many more types, see:
#   http://en.wikipedia.org/wiki/List_of_prime_numbers
#   http://mathworld.wolfram.com/IntegerSequencePrimes.html

my $show_safe = 0;
my $show_sophie = 0;
my $show_twin = 0;
my $show_lucas = 0;
my $show_fibonacci = 0;
my $show_palindromic = 0;
my $show_usage = 0;
my $segment_size = 30 * 128_000;   # 128kB

GetOptions( "safe"      => \$show_safe,
            "sophie|sg" => \$show_sophie,
            "twin"      => \$show_twin,
            "lucas"     => \$show_lucas,
            "fibonacci" => \$show_fibonacci,
            "palindromic|palindrome|palendrome" => \$show_palindromic,
            "help"      => \$show_usage,
          ) || die_usage();
die_usage() if $show_usage;

# Get the start and end values.  Verify they're positive integers.
die_usage() unless @ARGV == 2;
my ($start, $end) = @ARGV;
# Allow "10**100" as arguments
$start =~ s/^(\d+)\*\*(\d+)$/Math::BigInt->new($1)->bpow($2)/e;
$end   =~ s/^(\d+)\*\*(\d+)$/Math::BigInt->new($1)->bpow($2)/e;
die "$start isn't a positive integer" if $start =~ tr/0123456789//c;
die "$end isn't a positive integer" if $end =~ tr/0123456789//c;

# Turn start and end into bigints if they're very large
if ( ($start >= ~0 && $start ne ~0) || ($end >= ~0 && $end ne ~0) ) {
  $start = Math::BigInt->new($start);
  $end = Math::BigInt->new($end);
}

# Fibonacci numbers
{
  my @fibs = (Math::BigInt->new(0), Math::BigInt->new(1));
  sub fib {
    my $n = shift;
    return $n if $n < 2;
    if (!defined $fibs[$n]) {
      my ($nm2, $nm1) = ($fibs[-2],$fibs[-1]);
      for (scalar @fibs .. $n) {
        ($nm2, $nm1) = ($nm1, $nm2 + $nm1);
        push @fibs, $nm1;
      }
    }
    return $fibs[$n];
  }
}

# Return all Lucas primes between start and end, using identity:
#     L_n = F_n-1 + F_n+1
sub lucas_primes {
  my ($start, $end) = @_;
  my @lprimes;
  my $prime = 2;
  my $k = 0;
  while ($prime < $start) {
    $k++;
    $prime = fib($k) + fib($k+2);
  }
  while ($prime <= $end) {
    push @lprimes, $prime if is_prime($prime);
    $k++;
    $prime = fib($k) + fib($k+2);
  }
  @lprimes;
}

sub fibonacci_primes {
  my ($start, $end) = @_;
  my @fprimes;
  my $prime = 2;
  my $k = 3;
  while ($prime < $start) {
    $k++;
    $prime = fib($k);
  }
  while ($prime <= $end) {
    push @fprimes, $prime if is_prime($prime);
    # For all but k=4, F_k is prime only when k is prime.
    $k = ($k <= 4)  ?  $k+1  :  next_prime($k);
    $prime = fib($k);
  }
  @fprimes;
}


my @p;
while ($start <= $end) {

    # This will need to be made more generic if we add more options like this.
    # FIXME: It also doesn't work for our twin prime selection.
    if ($show_lucas && $show_fibonacci) {
      my %l;
      undef @l{ lucas_primes($start, $end) };
      @p = grep { exists $l{$_} } fibonacci_primes($start, $end);
      $start = $end+1;

    } elsif ($show_lucas) {
      @p = lucas_primes($start, $end);
      $start = $end+1;

    } elsif ($show_fibonacci) {
      @p = fibonacci_primes($start, $end);
      $start = $end+1;

    } else {
      # small segment size if we're doing bigints
      $segment_size = 10000 if $start > ~0;

      if ( $show_palindromic && $start >= 100 && (length($start) % 2) == 0 ) {
        # anything with an even number of digits is divisible by 11
        $start = 10 ** (length($start)) + 1;
      }

      my $seg_start = $start;
      my $seg_end = $start + $segment_size;
      $seg_end = $end if $end < $seg_end;
      $start = $seg_end+1;
      # Skip ranges where there are no palindromic primes.
      # - anything starting with 24568 won't be a prime when reversed
      if ( $show_palindromic && $seg_start >= 100 &&
           length($seg_start) == length($seg_end)) {
        next if $seg_start =~ /^2/ && $seg_end =~ /^2/;
        next if $seg_start =~ /^8/ && $seg_end =~ /^8/;
        next if $seg_start =~ /^[456]/ && $seg_end =~ /^[456]/;
      }
      @p = @{primes($seg_start, $seg_end)};
      #warn "-- ", scalar @p, " primes between $start and $seg_end\n";
    }

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
      @p = grep { is_prime( ($_-1) >> 1 ); }
           grep { ($_ <= 7) || ($_ % 12) == 11; }      # Optimization
           @p;
      #warn "   reduced to ", scalar @p, " safe primes\n";
    }

    # Restrict to Sophie Germain primes if requested.
    if ($show_sophie) {
      @p = grep { is_prime( 2*$_+1 ); }
           grep { my $m30 = $_ % 30;
                  $_ <= 5 || $m30 == 11 || $m30 == 23 || $m30 == 29 ; }
           @p;
      #warn "   reduced to ", scalar @p, " SG primes\n";
    }

    # Restrict to Palendromic primes if requested.
    if ($show_palindromic) {
      @p = grep { $_ eq reverse $_; } @p;
      #warn "   reduced to ", scalar @p, " Palindromic primes\n";
    }

    # print this segment
    print join("\n", @p), "\n"  if scalar @p > 0;
}




sub die_usage {
  die <<EOU;
Usage: $0 [options]  START  END

Displays all primes between the positive integers START and END, inclusive.
The START and END values must be integers, however the shortcut "x**y" may
be used, which allows very large values (e.g. '10**500' or '2**64')

Options:
  --help      displays this help message
  --twin      displays only twin primes, where p+2 is also prime
  --safe      displays only safe primes, where (p-1)/2 is also prime
  --sophie    displays only Sophie Germain primes, where 2p+1 is also prime
  --lucas     displays only Lucas primes, where L_n is prime
  --fibonacci displays only Fibonacci primes, where F_n is prime
  --palindr   displays only Palindromic primes, where p eq reverse p

Note that options can be combined, e.g. display all safe twin primes.
EOU
}
