package Math::Prime::Util;
use strict;
use warnings;
use Carp qw/croak confess/;

BEGIN {
  $Math::Prime::Util::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::VERSION = '0.03';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw(
                     prime_precalc prime_free
                     is_prime is_prob_prime miller_rabin
                     primes
                     next_prime  prev_prime
                     prime_count prime_count_lower prime_count_upper prime_count_approx
                     nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
                     factor
                   );
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);

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


my $_maxparam = (_maxbits == 32) ? 4294967295 : 18446744073709551615;

sub primes {
  my $optref = {};  $optref = shift if ref $_[0] eq 'HASH';
  croak "no parameters to primes" unless scalar @_ > 0;
  croak "too many parameters to primes" unless scalar @_ <= 2;
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;
  my $sref = [];

  # Validate parameters
  if ( (!defined $low) || (!defined $high) ||
       ($low =~ tr/0123456789//c) || ($high =~ tr/0123456789//c)
     ) {
    croak "Parameters [ $low $high ] must be positive integers";
  }

  # Verify the parameters are in range.
  croak "Parameters [ $low $high ] not in range 0-$_maxparam" unless $low <= $_maxparam && $high <= $_maxparam;

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

# We use this object to let them control when memory is freed.
package Math::Prime::Util::MemFree;
use Carp qw/croak confess/;
my $memfree_instances = 0;
sub new {
  my $self = bless {}, shift;
  $memfree_instances++;
  return $self;
}
sub DESTROY {
  my $self = shift;
  confess "instances count mismatch" unless $memfree_instances > 0;
  Math::Prime::Util::prime_free if --$memfree_instances == 0;
}
package Math::Prime::Util;



1;

__END__


# ABSTRACT: Utilities related to prime numbers, including fast generators / sievers

=pod

=encoding utf8


=head1 NAME

Math::Prime::Util - Utilities related to prime numbers, including fast sieves and factoring


=head1 VERSION

Version 0.03


=head1 SYNOPSIS

  # Normally you would just import the functions you are using.
  # Nothing is exported by default.
  use Math::Prime::Util ':all';


  # Get a big array reference of many primes
  my $aref = primes( 100_000_000 );

  # All the primes between 5k and 10k inclusive
  my $aref = primes( 5_000, 10_000 );

  # If you want them in an array instead
  my @primes = @{primes( 500 )};


  # is_prime returns 0 for composite, 1 for prime
  say "$n is prime"  if is_prime($n);


  # step to the next prime
  $n = next_prime($n);

  # step back (returns 0 if given input less than 2)
  $n = prev_prime($n);


  # Return Pi(n) -- the number of primes E<lt>= n.
  my $number_of_primes = prime_count( 1_000_000 );

  # Quickly return an approximation to Pi(n)
  my $approx_number_of_primes = prime_count_approx( 10**17 );

  # Lower and upper bounds.  lower <= Pi(n) <= upper for all n
  die unless prime_count_lower($n) <= prime_count($n);
  die unless prime_count_upper($n) >= prime_count($n);


  # Return p_n, the nth prime
  say "The ten thousandth prime is ", nth_prime(10_000);

  # Return a quick approximation to the nth prime
  say "The one trillionth prime is ~ ", nth_prime_approx(10**12);

  # Lower and upper bounds.   lower <= nth_prime(n) <= upper for all n
  die unless nth_prime_lower($n) <= nth_prime($n);
  die unless nth_prime_upper($n) >= nth_prime($n);


  # Get the prime factors of a number
  @prime_factors = factor( $n );


  # Precalculate a sieve, possibly speeding up later work.
  prime_precalc( 1_000_000_000 );

  # Free any memory used by the module.
  prime_free;

  # Alternate way to free.  When this leaves scope, memory is freed.
  my $mf = Math::Prime::Util::MemFree->new;


=head1 DESCRIPTION

A set of utilities related to prime numbers.  These include multiple sieving
methods, is_prime, prime_count, nth_prime, approximations and bounds for
the prime_count and nth prime, next_prime and prev_prime, factoring utilities,
and more.

The default sieving and factoring are intended to be the fastest on CPAN,
including Math::Prime::XS, Math::Prime::FastSieve, and Math::Factor::XS.



=head1 FUNCTIONS

=head2 is_prime

  print "$n is prime" if is_prime($n);

Returns 2 if the number is prime, 0 if not.  Also note there are
probabilistic prime testing functions available.


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
memory.  It uses a segmented sieve if the primary sieve is not large enough,
so is very memory efficient.  However, time growth is slightly more than
linear with C<n>, so it can take a very long time.  A hybrid table approach
is the usual means taken to speed this up, which a later version may do.  For
very large inputs the methods of Meissel, Lehmer, or
Lagarias-Miller-Odlyzko-Deleglise-Rivat may be appropriate.


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
large inputs.  A segmented sieve is used for large inputs, so it is memory
efficient.  On my machine it will return the 203,280,221st prime (the largest
that fits in 32-bits) in 2.5 seconds.  The 10^9th prime takes 15 seconds to
find, while the 10^10th prime takes nearly four minutes.


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


=head2 miller_rabin

  my $maybe_prime = miller_rabin($n, 2);
  my $probably_prime = miller_rabin($n, 2, 3, 5, 7, 11, 13, 17);

Takes a positive number as input and one or more bases.  The bases must be
between C<2> and C<n - 2>.  Returns 2 is C<n> is definitely prime, 1 if C<n>
is probably prime, and 0 if C<n> is definitely composite.  Since this is
just the Miller-Rabin test, a value of 2 is only returned for inputs of
2 and 3, which are shortcut.  If 0 is returned, then the number really is a
composite.  If 1 is returned, we aren't sure.

This is usually used in combination with other tests to make either stronger
tests or deterministic results for numbers less than some verified limit.

=head2 is_prob_prime

  my $prob_prime = is_prob_prime($n);
  # Returns 0 (composite), 2 (prime), or 1 (probably prime)

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).

This is done with a tuned set of Miller-Rabin tests such that the result
should be deterministic for 64-bit input.  Either 2, 3, 4, 5, or 7 Miller-Rabin
tests are performed (no more than 3 for 32-bit input), and the result will
then always be 0 (composite) or 2 (prime).  A later implementation may switch
to a BPSW test, depending on speed.


=head1 UTILITY FUNCTIONS

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
making calls, or want things cleanup up, you can use this.  The object method
might be a better choice for complicated uses.

=head2 Math::Prime::Util::MemFree->new

  my $mf = Math::Prime::Util::MemFree->new;
  # perform operations.  When $mf goes out of scope, memory will be recovered.

This is a more robust way of making sure any cached memory is freed, as it
will be handled by the last C<MemFree> object leaving scope.  This means if
your routines were inside an eval that died, things will still get cleaned up.
If you call another function that uses a MemFree object, the cache will stay
in place because you still have an object.



=head1 FACTORING FUNCTIONS

=head2 factor

  my @factors = factor(3_369_738_766_071_892_021);

Produces the prime factors of a positive number input.  They will typically
but necessarily be in numerical order.  The special cases of C<n = 0> and
C<n = 1> will return C<n>, which guarantees multiplying the factors together
will always result in the input value, though those are the only cases where
the returned factors are not prime.

The current algorithm is to use trial division for 32-bit numbers, while for
larger numbers a small set of trial divisions is performed, followed by a
single run of SQUFOF, then trial division of the results.  This results in
faster factorization of most large numbers.  More sophisticated methods could
be used.

=head2 fermat_factor

  my @factors = fermat_factor($n);

Produces factors, not necessarily prime, of the positive number input.  The
particular algorithm is Knuth's algorithm C.  For small inputs this will be
very fast, but it slows down quite rapidly as the number of digits increases.
If there is a factor close to the midpoint (e.g. a semiprime p*q where p and
q are the same number of digits), then this will be very fast.

=head2 squfof_factor

  my @factors = squfof_factor($n);

Produces factors, not necessarily prime, of the positive number input.  It
is possible the function will be unable to find a factor, in which case a
single element, the input, is returned.

=head2 prho_factor

=head2 pbrent_factor

=head2 pminus1_factor

Attempts to find a factor using one of the probabilistic algorigthms of
Pollard Rho, Brent's modification of Pollard Rho, or Pollard's C<p - 1>.
These are more specialized algorithms usually used for pre-factoring very
large inputs, or checking very large inputs for naive mistakes.  If given
a prime input, or if just unlucky, these will take a long time to return
back the single input value.  There are cases where these can result in
finding a factor or very large inputs in remarkably short time, similar to
how Fermat's method works very well for factors near C<sqrt(n)>.  They are
also amenable to massively parallel searching.

For 64-bit input, these are unlikely to be of much use.  An optimized SQUFOF
implementation takes under 20 milliseconds to find a factor for any 62-bit
input on modern desktop computers.  Lightweight quadratic sieves are
typically much faster for inputs in the 19+ digit range.


=head1 LIMITATIONS

I have not completed testing all the functions near the word size limit
(e.g. C<2^32> for 32-bit machines).  Please report any problems you find.

The extra factoring algorithms are mildly interesting but really limited by
not being big number aware.

Perl versions earlier than 5.8.0 have issues with 64-bit.  The test suite will
try to determine if your Perl is broken.  This will show up in factoring tests.


=head1 PERFORMANCE

Counting the primes to C<10^10> (10 billion), with time in seconds.
Pi(10^10) = 455,052,511.

       1.9  primesieve 3.6 forced to use only a single thread
       5.6  Tomás Oliveira e Silva's segmented sieve v2 (Sep 2010)
       6.6  primegen (optimized Sieve of Atkin)
      11.2  Tomás Oliveira e Silva's segmented sieve v1 (May 2003)

       5.9  My SoE called in segments
      15.6  My Sieve of Eratosthenes using a mod-30 wheel
      17.2  A slightly modified verion of Terje Mathisen's mod-30 sieve
      35.5  Basic Sieve of Eratosthenes on odd numbers
      33.4  Sieve of Atkin, from Praxis (not correct)
      72.8  Sieve of Atkin, 10-minute fixup of basic algorithm
      91.6  Sieve of Atkin, Wikipedia-like

Perl modules, counting the primes to C<800_000_000> (800 million), in seconds:

  Time (s)  Module                      Version
  --------  --------------------------  ------
       0.4  Math::Prime::Util           0.02
       0.9  Math::Prime::Util           0.01
       2.9  Math::Prime::FastSieve      0.12
      11.7  Math::Prime::XS             0.29
      15.0  Bit::Vector                 7.2
   [hours]  Math::Primality             0.04


I have not done extensive timing on factoring.  The difference in speed for
most inputs between Math::Factor::XS and Math::Prime::Util is small.  For
some ranges M::F::XS is faster, and some ranges it is slower.  Factoring
random semiprimes is about 1.1x to 6x faster with M::P::U, and some inputs
are even faster (e.g. C<204518747 * 16476429743>, which is ~50x faster).


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

The SQUFOF implementation being used is my modifications to Ben Buhrow's
modifications to Bob Silverman's code.  I may experiment with some other
implementations (Ben Buhrows and Jason Papadopoulos both have published
excellent versions in the public domain).



=head1 COPYRIGHT

Copyright 2011-2012 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
