package Math::Prime::Util;
use strict;
use warnings;
use Carp qw/croak confess/;

BEGIN {
  $Math::Prime::Util::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::VERSION = '0.01';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw(
                     prime_precalc prime_free
                     is_prime
                     primes
                     next_prime  prev_prime
                     prime_count prime_count_lower prime_count_upper prime_count_approx
                     nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
                   );

BEGIN {
  eval {
    require XSLoader;
    XSLoader::load(__PACKAGE__, $Math::Prime::Util::VERSION);
    prime_precalc(0);
    1;
  } or do {
    # We could insert a Pure Perl implementation here.
    croak "XS Code not available: $@";
  }
}


sub primes {
  my $optref = {};  $optref = shift if ref $_[0] eq 'HASH';
  croak "no parameters to primes" unless scalar @_ > 0;
  croak "too many parameters to primes" unless scalar @_ <= 2;
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;
  my $sref = [];

  # Validate parameters
  if ( (!defined $low) || (!defined $high) ||
       #($low < 0) || ($high < 0) ||
       ($low =~ tr/0123456789//c) || ($high =~ tr/0123456789//c)
     ) {
    croak "Parameters must be positive integers";
  }
  return $sref if ($low > $high) || ($high < 2);

  my $method = $optref->{'method'};
  $method = 'Dynamic' unless defined $method;

  if ($method =~ /^(Dyn\w*|Default|Generate)$/i) {
    # Dynamic -- we should try to do something smart.

    # Tiny range?
    if (($low+1) >= $high) {
      $method = 'Trial';

    # Fast for cached sieve?
    } elsif (($high <= (65536*30)) || ($high <= _get_prime_cache_size)) {
      $method = 'Sieve';

    # More memory than we should reasonably use for base sieve?
    } elsif ($high > (32*1024*1024*30)) {
      $method = 'Segment';

    # Only want half or less of the range low-high ?
    } elsif ( int($high / ($high-$low)) >= 2 ) {
      $method = 'Segment';

    } else {
      $method = 'Sieve';
    }
  }

  if    ($method =~ /^Trial$/i)     { $sref = trial_primes($low, $high); }
  elsif ($method =~ /^Erat\w*$/i)   { $sref = erat_primes($low, $high); }
  elsif ($method =~ /^Simple\w*$/i) { $sref = erat_simple_primes($low, $high); }
  elsif ($method =~ /^Seg\w*$/i)    { $sref = segment_primes($low, $high); }
  elsif ($method =~ /^Sieve$/i)     { $sref = sieve_primes($low, $high); }
  else { croak "Unknown prime method: $method"; }

  #return (wantarray) ? @{$sref} : $sref;
  return $sref;
}


1;

__END__


# ABSTRACT: Utilities related to prime numbers, including fast generators / sievers

=pod

=encoding utf8


=head1 NAME

Math::Prime::Util - Utilities related to prime numbers, including fast generators / sievers


=head1 VERSION

Version 0.01


=head1 SYNOPSIS

  use Math::Prime::Util qw/primes/;

  # Get a big array reference of many primes
  my $aref = primes( 100_000_000 );

  # All the primes between 5k and 10k inclusive
  my $aref = primes( 5_000, 10_000 );



=head1 DESCRIPTION

A set of utilities related to prime numbers.  These include multiple sieving
methods, is_prime, prime_count, nth_prime, approximations and bounds for
the prime_count and nth prime, next_prime and prev_prime, factoring utilities,
and more.

The default sieving and factoring are intended to be the fastest on CPAN,
including Math::Prime::XS, Math::Prime::FastSieve, and Math::Factor::XS.



=head1 FUNCTIONS

=head2 prime_precalc

  prime_precalc( 1_000_000_000 );

Let the module prepare for fast operation up to a specific number.  It is not
necessary to call this, but it gives you more control over when memory is
allocated and gives faster results for multiple calls in some cases.  In the
current implementation this will calculate a sieve for all numbers up to the
specified number.


=head2 prime_free

  prime_free;

Frees any extra memory the module may have allocated.  Like with
C<prime_precalc>, it is not necessary to call this, but if you're done
making calls, or want things cleanup up, you can use this.


=head2 is_prime

Returns true if the number is prime, false if not.

  print "$n is prime" if is_prime($n);


=head2 primes

Returns all the primes between the lower and upper limits (inclusive), with
a lower limit of C<2> if none is given.

An array reference is returned (with large lists this is much faster and uses
less memory than returning an array directly).

  my $aref1 = primes( 1_000_000 );
  my $aref2 = primes( 1_000_000_000_000, 1_000_000_001_000 );

  my @primes = @{ primes( 500 ) };

  print "$_\n" for (@{primes( 20, 100 )});

Sieving will be done if required.  The algorithm used will depend on the range
and whether a sieve result already exists.  Possibilities include trial
division (for ranges with only one expected prime), a Sieve of Eratosthenes
using wheel factorization, or a segmented sieve.


=head2 next_prime

  $n = next_prime($n);

Returns the next prime greater than the input number.


=head2 prev_prime

  $n = prev_prime($n);

Returns the prime smaller than the input number.  0 is returned if the
input is C<2> or lower.


=head2 prime_count

  my $number_of_primes = prime_count( 1_000 );

Returns the Prime Count function C<Pi(n)>.  The current implementation relies
on sieving to find the primes within the interval, so will take some time and
memory.  There are slightly faster ways to handle the sieving (e.g. maintain
a list of counts from C<2 - j> for various C<j>, then do a segmented sieve
between C<j> and C<n>), and for very large numbers the methods of Meissel,
Lehmer, or Lagarias-Miller-Odlyzko-Deleglise-Rivat may be appropriate.


=head2 prime_count_upper

=head2 prime_count_lower

  my $lower_limit = prime_count_lower($n);
  die unless prime_count($n) >= $lower_limit;

  my $upper_limit = prime_count_upper($n);
  die unless prime_count($n) <= $upper_limit;

Returns an upper or lower bound on the number of primes below the input number.
These are analytical routines, so will take a fixed amount of time and no
memory.  The actual C<prime_count> will always be on or between these numbers.

A common place these would be used is sizing an array to hold the first C<$n>
primes.  It may be desirable to use a bit more memory than is necessary, to
avoid calling C<prime_count>.

These routines use hand-verified tight limits below a range at least C<2^32>,
and fall back to the proven Dusart bounds of

    x/logx * (1 + 1/logx + 1.80/log^2x) <= Pi(x)

    x/logx * (1 + 1/logx + 2.51/log^2x) >= Pi(x)

above that range.


=head2 prime_count_approx

  print "there are about ",
        prime_count_approx( 10 ** 18 ),
        " primes below one quintillion.\n";

Returns an approximation to the C<prime_count> function, without having to
generate any primes.  The results are very close for small numbers, but less
so with large ranges.  The current implementation is 0.00033% too small
for the example, but takes under a microsecond and no memory to get the
result.


=head2 nth_prime

  say "The ten thousandth prime is ", nth_prime(10_000);

Returns the prime that lies in index C<n> in the array of prime numbers.  Put
another way, this returns the smallest C<p> such that C<Pi(p) E<gt>= n>.

This relies on generating primes, so can require a lot of time and space for
large inputs.


=head2 nth_prime_upper

=head2 nth_prime_lower

  my $lower_limit = nth_prime_lower($n);
  die unless nth_prime($n) >= $lower_limit;

  my $upper_limit = nth_prime_upper($n);
  die unless nth_prime($n) <= $upper_limit;

Returns an analytical upper or lower bound on the Nth prime.  This will be
very fast.  The lower limit uses the Dusart 1999 bounds for all C<n>, while
the upper bound uses one of the two Dusart 1999 bounds for C<n E<gt>= 27076>,
the Robin 1983 bound for C<n E<gt>= 7022>, and the simple bound of
C<n * (logn + loglogn)> for C<n E<lt> 7022>.


=head2 nth_prime_approx

  say "The one trillionth prime is ~ ", nth_prime_approx(10**12);

Returns an approximation to the C<nth_prime> function, without having to
generate any primes.  Uses the Cipolla 1902 approximation with two
polynomials, plus a correction term for small values to reduce the error.


=head1 LIMITATIONS

Some of the functions have not yet transitioned to using a segmented sieve,
so will use too much memory to be practical when called with very large
numbers (C<10^11> and up).

I have not completed testing all the functions near the word size limit
(e.g. C<2^32> for 32-bit machines).  Please report any problems you find.


=head1 PERFORMANCE

Counting the primes to C<10^10> (10 billion), with time in seconds.
Pi(10^10) = 455,052,511.

       1.9  primesieve 3.6 forced to use only a single thread
       5.6  Tomás Oliveira e Silva's segmented sieve v2 (Sep 2010)
       6.6  primegen (optimized Sieve of Atkin)
      11.2  Tomás Oliveira e Silva's segmented sieve v1 (May 2003)
  
      15.6  My Sieve of Eratosthenes using a mod-30 wheel
      17.2  A slightly modified verion of Terje Mathisen's mod-30 sieve
      35.5  Basic Sieve of Eratosthenes on odd numbers
      33.4  Sieve of Atkin, from Praxis (not correct)
      72.8  Sieve of Atkin, 10-minute fixup of basic algorithm
      91.6  Sieve of Atkin, Wikipedia-like

Perl modules, counting the primes to C<800_000_000> (800 million), in seconds:

       0.9  Math::Prime::Util           0.01
       2.9  Math::Prime::FastSieve      0.12
      11.7  Math::Prime::XS             0.29
      15.0  Bit::Vector                 7.2
   (hours)  Math::Primality             0.04
    

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 ACKNOWLEDGEMENTS

Eratosthenes of Cyrene provided the elegant and simple algorithm for finding
the primes.

Terje Mathisen, A.R. Quesada, and B. Van Pelt all had useful ideas which I
used in my wheel sieve.

Tomás Oliveira e Silva has released the source for a very fast segmented sieve.
The current implementation does not use these ideas, but future versions likely
will.


=head1 COPYRIGHT

Copyright 2011-2012 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
