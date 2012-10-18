#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Math::BigInt try => 'GMP';
use Math::Prime::Util qw/primes is_prime next_prime prev_prime
                         prime_count primorial pn_primorial/;
$| = 1;

# For many more types, see:
#   http://en.wikipedia.org/wiki/List_of_prime_numbers
#   http://mathworld.wolfram.com/IntegerSequencePrimes.html

my $segment_size = 30 * 128_000;   # 128kB
my %opts;
GetOptions(\%opts,
           'safe|A005385',
           'sophie|sg|A005384',
           'twin|A001359',
           'lucas|A005479',
           'fibonacci|A005478',
           'triplet|A007529',
           'quadruplet|A007530',
           'cousin|A023200',
           'sexy|A023201',
           'mersenne|A000668',
           'palindromic|palindrome|palendrome|A002385',
           'pillai|A063980',
           'good|A028388',
           'cuban1|A002407',
           'cuban2|A002648',
           'pnp1|A005234',
           'pnm1|A006794',
           'euclid|A018239',
           'help',
          ) || die_usage();
die_usage() if exists $opts{'help'};


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

sub mersenne_primes {
  my ($start, $end) = @_;
  my @mprimes;
  my $p = 1;
  while (1) {
    $p = next_prime($p);  # Mp is not prime if p is not prime
    next if $p > 3 && ($p % 4) == 3 && is_prime(2*$p+1);
    my $Mp = Math::BigInt->bone->blsft($p)->bsub(1);
    last if $Mp > $end;
    # Lucas-Lehmer test would be faster
    push @mprimes, $Mp if $Mp >= $start && is_prime($Mp);
  }
  @mprimes;
}

sub euclid_primes {
  my ($start, $end, $add) = @_;
  my @eprimes;
  my $k = 0;
  while (1) {
    my $primorial = pn_primorial(Math::BigInt->new($k)) + $add;
    last if $primorial > $end;
    push @eprimes, $primorial if $primorial >= $start && is_prime($primorial);
    $k++;
  }
  @eprimes;
}

# This is not a general palindromic digit function!
sub ndig_palindromes {
  my $digits = shift;
  return (2,3,5,7,9) if $digits == 1;
  return (11) if $digits == 2;
  return () if ($digits % 2) == 0;

  my @prefixes = (1,3,7,9);
  my $inner_digits = ($digits-1) / 2 - 1;
  foreach my $d (1 .. $inner_digits) {
    @prefixes = map { ($_.'0', $_.'1', $_.'2', $_.'3', $_.'4',
                       $_.'5', $_.'6', $_.'7', $_.'8', $_.'9',); } @prefixes;
  }
  return map { my $r = reverse($_);
               ($_.'0'.$r, $_.'1'.$r, $_.'2'.$r, $_.'3'.$r,
                $_.'4'.$r, $_.'5'.$r, $_.'6'.$r, $_.'7'.$r,
                $_.'8'.$r, $_.'9'.$r,);
             } @prefixes;
}

# See: http://en.wikipedia.org/wiki/Pillai_prime
# This is quite slow.
sub is_pillai {
  my $p = shift;
  return 0 if $p < 2;
  # return 0 unless is_prime($p);
  my $n_factorial_mod_p = Math::BigInt->bone();
  for my $n (2 .. $p-1) {
    $n_factorial_mod_p->bmul($n)->bmod($p);
    next if $p % $n == 1;
    return 1 if $n_factorial_mod_p == ($p-1);
  }
  0;
}

# Also slow.
sub is_good_prime {
  my $p = shift;
  return 0 if $p <= 2; # 2 isn't a good prime
  my $lower = $p;
  my $upper = $p;
  while ($lower > 2) {
    $lower = prev_prime($lower);
    $upper = next_prime($upper);
    return 0 if ($p*$p) <= ($upper * $lower);
  }
  1;
}

sub merge_primes {
  my ($genref, $pref, $name, @primes) = @_;
  if (!defined $$genref) {
    @$pref = @primes;
    $$genref = $name;
  } else {
    my %f;
    undef @f{ @primes };
    @$pref = grep { exists $f{$_} } @$pref;
  }
}

# This is used for things that can generate a filtered list faster than
# searching through all primes in the range.
sub gen_and_filter {
  my ($start, $end) = @_;
  my $gen;
  my @p;

  if (exists $opts{'lucas'}) {
    merge_primes(\$gen, \@p, 'lucas', lucas_primes($start, $end));
  }
  if (exists $opts{'fibonacci'}) {
    merge_primes(\$gen, \@p, 'fibonacci', fibonacci_primes($start, $end));
  }
  if (exists $opts{'mersenne'}) {
    merge_primes(\$gen, \@p, 'mersenne', mersenne_primes($start, $end));
  }
  if (exists $opts{'euclid'}) {
    merge_primes(\$gen, \@p, 'euclid', euclid_primes($start, $end, 1));
  }
  if (exists $opts{'palindromic'}) {
    if (!defined $gen) {
      foreach my $d (length($start) .. length($end)) {
        push @p, grep { $_ >= $start && $_ <= $end && is_prime($_) }
                 ndig_palindromes($d);
      }
      $gen = 'palindromic';
    } else {
      @p = grep { $_ eq reverse $_; } @p;
    }
  }

  if (!defined $gen) {
    @p = @{primes($start, $end)};
    $gen = 'primes';
  }

  if (exists $opts{'twin'}) {
    if ($gen ne 'primes') {
      @p = grep { is_prime( $_+2 ); } @p;
    } elsif (scalar @p > 0) {
      # All primes in the range are here, so just look in the array.
      push @p, is_prime($p[-1]+2) ? $p[-1]+2 : 0;
      my @twin;
      my $prime = shift @p;
      foreach my $next (@p) {
        push @twin, $prime if $prime+2 == $next;
        $prime = $next;
      }
      @p = @twin;
    }
  }

  if (exists $opts{'triplet'}) {   # could be optimized like twin
    @p = grep { is_prime($_+6) && (is_prime($_+2) || is_prime($_+4)); } @p;
  }

  if (exists $opts{'quadruplet'}) {   # could be optimized like twin
    @p = grep { is_prime($_+2) && is_prime($_+6) && is_prime($_+8); }
         grep { $_ <= 5 || ($_ % 30) == 11; }
         @p;
  }

  if (exists $opts{'cousin'}) {   # could be optimized like twin
    @p = grep { is_prime($_+4); }
         grep { ($_ <= 3) || ($_ % 6) == 1; }
         @p;
  }

  if (exists $opts{'sexy'}) {   # could be optimized like twin
    @p = grep { is_prime($_+6); } @p;
  }

  if (exists $opts{'safe'}) {
    @p = grep { is_prime( ($_-1) >> 1 ); }
         grep { ($_ <= 7) || ($_ % 12) == 11; }
         @p;
  }
  if (exists $opts{'sophie'}) {
    @p = grep { is_prime( 2*$_+1 ); }
         grep { my $m30 = $_ % 30;
                $_ <= 5 || $m30 == 11 || $m30 == 23 || $m30 == 29 ; }
         @p;
  }
  if (exists $opts{'cuban1'}) {
    @p = grep { my $n = sqrt((4*$_-1)/3);  $n == int($n); } @p;
  }
  if (exists $opts{'cuban2'}) {
    @p = grep { my $n = sqrt(($_-1)/3);  $n == int($n); } @p;
  }
  if (exists $opts{'pnm1'}) {
    @p = grep { is_prime( primorial(Math::BigInt->new($_))-1 ) } @p;
  }
  if (exists $opts{'pnp1'}) {
    @p = grep { is_prime( primorial(Math::BigInt->new($_))+1 ) } @p;
  }
  if (exists $opts{'pillai'}) {
    @p = grep { is_pillai($_); } @p;
  }
  if (exists $opts{'good'}) {
    @p = grep { is_good_prime($_); } @p;
  }
  @p;
}

if (exists $opts{'lucas'} ||
    exists $opts{'fibonacci'} ||
    exists $opts{'palindromic'} ||
    exists $opts{'euclid'} ||
    exists $opts{'mersenne'}) {
  my @p = gen_and_filter($start, $end);
  print join("\n", @p), "\n"  if scalar @p > 0;
} else {
  my @p;
  while ($start <= $end) {

    # Adjust segment sizes for some cases
    $segment_size = 10000 if $start > ~0;   # small if doing bigints
    if (exists $opts{'pillai'}) {
      $segment_size = ($start < 10000) ? 100 : 1000;  # very small for Pillai
    }

    my $seg_start = $start;
    my $seg_end = $start + $segment_size;
    $seg_end = $end if $end < $seg_end;
    $start = $seg_end+1;

    @p = gen_and_filter($seg_start, $seg_end);

    # print this segment
    print join("\n", @p), "\n"  if scalar @p > 0;
  }
}




sub die_usage {
  die <<EOU;
Usage: $0 [options]  START  END

Displays all primes between the positive integers START and END, inclusive.
The START and END values must be integers, however the shortcut "x**y" may
be used, which allows very large values (e.g. '10**500' or '2**64')

General options:

  --help       displays this help message

Filter options, which will cause the list of primes to be further filtered
to only those primes additionally meeting these conditions:

  --twin       Twin             p+2 is prime
  --triplet    Triplet          p+6 and (p+2 or p+4) are prime
  --quadruplet Quadruplet       p+2, p+6, and p+8 are prime
  --cousin     Cousin           p+4 is prime
  --sexy       Sexy             p+6 is prime
  --safe       Safe             (p-1)/2 is also prime
  --sophie     Sophie Germain   2p+1 is also prime
  --lucas      Lucas            L_p is prime
  --fibonacci  Fibonacci        F_p is prime
  --mersenne   Mersenne         M_p = 2^p-1 is prime
  --palindr    Palindromic      p is equal to p with its base-10 digits reversed
  --pillai     Pillai           n! % p = p-1 and p % n != 1 for some n
  --good       Good             p_n^2 > p_{n-i}*p_{n+i} for all i in (1..n-1)
  --cuban1     Cuban (y+1)      p = (x^3 - y^3)/(x-y), x=y+1
  --cuban2     Cuban (y+2)      p = (x^3 - y^3)/(x-y), x=y+2
  --pnp1       Primorial+1      p#+1 is prime
  --pnm1       Primorial-1      p#-1 is prime
  --euclid     Euclid           pn#+1 is prime

Note that options can be combined, e.g. display only safe twin primes.
In all cases involving multiples (twin, triplet, etc.), the value returned
is p -- the least value of the set.

EOU
}
