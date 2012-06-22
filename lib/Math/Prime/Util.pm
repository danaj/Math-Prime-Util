package Math::Prime::Util;
use strict;
use warnings;
use Carp qw/croak confess carp/;

BEGIN {
  $Math::Prime::Util::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::VERSION = '0.07';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw(
                     prime_precalc prime_memfree
                     is_prime is_prob_prime miller_rabin
                     primes
                     next_prime  prev_prime
                     prime_count prime_count_lower prime_count_upper prime_count_approx
                     nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
                     random_prime random_ndigit_prime
                     factor all_factors
                     ExponentialIntegral LogarithmicIntegral RiemannR
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
END {
  _prime_memfreeall;
}


my $_maxparam  = (_maxbits == 32) ? 4294967295 : 18446744073709551615;
my $_maxdigits = (_maxbits == 32) ? 10 : 20;
my $_maxprime  = (_maxbits == 32) ? 4294967291 : 18446744073709551557;
my $_maxprimeidx=(_maxbits == 32) ?  203280221 :   425656284035217743;

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
    $low  = 'undef' unless defined $low;
    $high = 'undef' unless defined $high;
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

  if ($method =~ /^Simple\w*$/i) {
    carp "Method 'Simple' is deprecated.";
    $method = 'Erat';
  }

  if    ($method =~ /^Trial$/i)     { $sref = trial_primes($low, $high); }
  elsif ($method =~ /^Erat\w*$/i)   { $sref = erat_primes($low, $high); }
  elsif ($method =~ /^Seg\w*$/i)    { $sref = segment_primes($low, $high); }
  elsif ($method =~ /^Sieve$/i)     { $sref = sieve_primes($low, $high); }
  else { croak "Unknown prime method: $method"; }

  # Using this line:
  #   return (wantarray) ? @{$sref} : $sref;
  # would allow us to return an array ref in scalar context, and an array
  # in array context.  Handy for people who might write:
  #   @primes = primes(100);
  # but I think the dual interface could bite us later.
  return $sref;
}

sub random_prime {
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;
  if ( (!defined $low) || (!defined $high) ||
       ($low =~ tr/0123456789//c) || ($high =~ tr/0123456789//c)
     ) {
    $low  = 'undef' unless defined $low;
    $high = 'undef' unless defined $high;
    croak "Parameters [ $low $high ] must be positive integers";
  }
  croak "Parameters [ $low $high ] not in range 0-$_maxparam" unless $low <= $_maxparam && $high <= $_maxparam;
  $low = 2 if $low < 2;

  # Make sure we have a valid range.
  # TODO: this is is killing performance with large numbers
  $low = next_prime($low-1);
  $high = ($high < ~0) ? prev_prime($high+1) : prev_prime($high);
  return $low if ($low == $high) && is_prime($low);
  return if $low >= $high;

  # At this point low and high are both primes, and low < high.
  my $range = $high - $low + 1;
  my $prime;

  # Note:  I was using rand($range), but Math::Random::MT ignores the argument
  #        instead of following its documentation.
  my $irandf = (exists &::rand) ? sub { return int(::rand()*shift); }
                                : sub { return int(rand()*shift); };

  if ($high < 30000) {
    # nice deterministic solution, but gets very costly with large values.
    my $li = ($low == 2) ? 1 : prime_count($low);
    my $hi = prime_count($high);
    my $irange = $hi - $li + 1;
    my $rand = $irandf->($irange);
    $prime = nth_prime($li + $rand);
  } else {
    # random loop
    my $loop_limit = 2000 * 1000;  # To protect against broken rand
    if ($range <= 4294967295) {
      do {
        $prime = $low + $irandf->($range);
        croak "Random function broken?" if $loop_limit-- < 0;
      }  while ( !($prime%2) || !($prime%3) || !is_prime($prime) );
    } else {
      do {
        my $rand = ( ($irandf->(4294967295) << 32) + $irandf->(4294967295) ) % $range;
        $prime = $low + $rand;
        croak "Random function broken?" if $loop_limit-- < 0;
      }  while ( !($prime%2) || !($prime%3) || !is_prime($prime) );
    }
  }
  return $prime;
}

sub random_ndigit_prime {
  my $digits = shift;
  if ((!defined $digits) || ($digits > $_maxdigits) || ($digits < 1)) {
    croak "Digits must be between 1 and $_maxdigits";
  }
  my $low = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  return random_prime($low, $max);
}

# Perhaps a random_nbit_prime ?   Definition?

sub all_factors {
  my $n = shift;
  my @factors = factor($n);
  my %all_factors;
  foreach my $f1 (@factors) {
    next if $f1 >= $n;
    my @all = keys %all_factors;;
    foreach my $f2 (@all) {
      $all_factors{$f1*$f2} = 1 if ($f1*$f2) < $n;
    }
    $all_factors{$f1} = 1;
  }
  @factors = sort {$a<=>$b} keys %all_factors;
  return @factors;
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
  Math::Prime::Util::prime_memfree if --$memfree_instances == 0;
  return;
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

Version 0.07


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


  # is_prime returns 0 for composite, 2 for prime
  say "$n is prime"  if is_prime($n);

  # is_prob_prime returns 0 for composite, 2 for prime, and 1 for maybe prime
  say "$n is ", qw(composite maybe_prime? prime)[is_prob_prime($n)];


  # step to the next prime (returns 0 if the next one is more than ~0)
  $n = next_prime($n);

  # step back (returns 0 if given input less than 2)
  $n = prev_prime($n);


  # Return Pi(n) -- the number of primes E<lt>= n.
  $primepi = prime_count( 1_000_000 );
  $primepi = prime_count( 10**14, 10**14+1000 );  # also does ranges

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
  prime_memfree;

  # Alternate way to free.  When this leaves scope, memory is freed.
  my $mf = Math::Prime::Util::MemFree->new;


  # Random primes
  my $small_prime = random_prime(1000);      # random prime <= limit
  my $rand_prime = random_prime(100, 10000); # random prime within a range
  my $rand_prime = random_ndigit_prime(6);   # random 6-digit prime


=head1 DESCRIPTION

A set of utilities related to prime numbers.  These include multiple sieving
methods, is_prime, prime_count, nth_prime, approximations and bounds for
the prime_count and nth prime, next_prime and prev_prime, factoring utilities,
and more.

All routines currently work in native integers (32-bit or 64-bit).  Bignum
support may be added later.  If you need bignum support for these types of
functions inside Perl now, I recommend L<Math::Pari>.

The default sieving and factoring are intended to be (and currently are)
the fastest on CPAN, including L<Math::Prime::XS>, L<Math::Prime::FastSieve>,
and L<Math::Factor::XS>.  It seems to be faster than L<Math::Pari> for
everything except factoring certain 16-20 digit numbers.

The module is thread-safe and allows concurrency between Perl threads while
still sharing a prime cache.  It is not itself multithreaded.


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

Returns the next prime greater than the input number.  0 is returned if the
next prime is larger than a native integer type (the last representable
primes being C<4,294,967,291> in 32-bit Perl and
C<18,446,744,073,709,551,557> in 64-bit).


=head2 prev_prime

  $n = prev_prime($n);

Returns the prime smaller than the input number.  0 is returned if the
input is C<2> or lower.


=head2 prime_count

  my $primepi = prime_count( 1_000 );
  my $pirange = prime_count( 1_000, 10_000 );

Returns the Prime Count function C<Pi(n)>, also called C<primepi> in some
math packages.  When given two arguments, it returns the inclusive
count of primes between the ranges (e.g. C<(13,17)> returns 2, C<14,17>
and C<13,16> return 1, and C<14,16> returns 0).

The current implementation relies on sieving to find the primes within the
interval, so will take some time and memory.  It uses a segmented sieve so
is very memory efficient, and also allows fast results even with large
base values.  The complexity for C<prime_count(a, b)> is approximately
C<O(sqrt(a) + (b-a))>, where the first term is typically negligible below
C<~ 10^11>.  Memory use is proportional only to C<sqrt(a)>, with total
memory use under 1MB for any base under C<10^14>.

A later implementation may work on improving performance for values, both
in reducing memory use (the current maximum is 140MB at C<2^64>) and improving
speed.  Possibilities include a hybrid table approach, using an explicit
formula with C<li(x)> or C<R(x)>, or one of the Meissel, Lehmer,
or Lagarias-Miller-Odlyzko-Deleglise-Rivat methods.


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

These routines use hand-verified tight limits below a range at least C<2^35>,
and fall back to the Dusart bounds of

    x/logx * (1 + 1/logx + 1.80/log^2x) <= Pi(x)

    x/logx * (1 + 1/logx + 2.51/log^2x) >= Pi(x)

above that range.


=head2 prime_count_approx

  print "there are about ",
        prime_count_approx( 10 ** 18 ),
        " primes below one quintillion.\n";

Returns an approximation to the C<prime_count> function, without having to
generate any primes.  The current implementation uses the Riemann R function
which is quite accurate: an error of less than C<0.0005%> is typical for
input values over C<2^32>.  A slightly faster (0.1ms vs. 1ms), but much less
accurate, answer can be obtained by averaging the upper and lower bounds.


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
composite.  If 1 is returned, there is a good chance the number is prime
(depending on the input and the bases), but we cannot be sure.

This is usually used in combination with other tests to make either stronger
tests (e.g. the strong BPSW test) or deterministic results for numbers less
than some verified limit (such as the C<is_prob_prime> function in this module).


=head2 is_prob_prime

  my $prob_prime = is_prob_prime($n);
  # Returns 0 (composite), 2 (prime), or 1 (probably prime)

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).

This is done with a tuned set of Miller-Rabin tests such that the result
will be deterministic for 64-bit input.  Either 2, 3, 4, 5, or 7 Miller-Rabin
tests are performed (no more than 3 for 32-bit input), and the result will
then always be 0 (composite) or 2 (prime).  A later implementation may switch
to a BPSW test, depending on speed.


=head2 random_prime

  my $small_prime = random_prime(1000);      # random prime <= limit
  my $rand_prime = random_prime(100, 10000); # random prime within a range

Returns a psuedo-randomly selected prime that will be greater than or equal
to the lower limit and less than or equal to the upper limit.  If no lower
limit is given, 2 is implied.  Returns undef if no primes exist within the
range.  The L<rand> function is called one or more times for selection.

This will return a uniform distribution of the primes in the range, meaning
for each prime in the range, the chances are equally likely that it will be
seen.

The current algorithm does a random index selection for small numbers, which
is deterministic.  For larger numbers, this can be very slow, so the
obvious Monte Carlo method is used, where random numbers in the range are
selected until one is prime.  That also gets slow as the number of digits
increases, but not something that impacts us in 64-bit.

If you want cryptographically secure primes, I suggest looking at
L<Crypt::Primes> or something similar.  The current L<Math::Prime::Util>
module does not use strong randomness, and its primes are ridiculously small
by cryptographic standards.

Perl's L<rand> function is normally called, but if the sub C<main::rand>
exists, it will be used instead.  When called with no arguments it should
return a float value between 0 and 1-epsilon, with 32 bits of randomness.
Examples:

  # Use Mersenne Twister
  use Math::Random::MT::Auto qw/rand/;

  # Use my custom random function
  sub rand { ... }

=head2 random_ndigit_prime

  say "My 4-digit prime number is: ", random_ndigit_prime(4);

Selects a random n-digit prime, where the input is an integer number of
digits between 1 and the maximum native type (10 for 32-bit, 20 for 64-bit).
One of the primes within that range (e.g. 1000 - 9999 for 4-digits) will be
uniformly selected using the L<rand> function.



=head1 UTILITY FUNCTIONS

=head2 prime_precalc

  prime_precalc( 1_000_000_000 );

Let the module prepare for fast operation up to a specific number.  It is not
necessary to call this, but it gives you more control over when memory is
allocated and gives faster results for multiple calls in some cases.  In the
current implementation this will calculate a sieve for all numbers up to the
specified number.


=head2 prime_memfree

  prime_memfree;

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

Produces the prime factors of a positive number input.  They may not be in
numerical order.  The special cases of C<n = 0> and C<n = 1> will
return C<n>, which guarantees multiplying the factors together will
always result in the input value, though those are the only cases where
the returned factors are not prime.

The current algorithm is to use trial division for small numbers, while large
numbers go through a sequence of small trials, SQUFOF, Pollard's Rho, Hart's
one line factorization, and finally trial division for any survivors.  This
process is repeated for each non-prime factor.


=head2 all_factors

  my @divisors = all_factors(30);   # returns (2, 3, 5, 6, 10, 15)

Produces all the divisors of a positive number input.  1 and the input number
are excluded (which implies that an empty list is returned for any prime
number input).  The divisors are a power set of multiplications of the prime
factors, returned as a uniqued sorted list.


=head2 trial_factor

  my @factors = trial_factor($n);

Produces the prime factors of a positive number input.  The factors will be
in numerical order.  The special cases of C<n = 0> and C<n = 1> will return
C<n>, while with all other inputs the factors are guaranteed to be prime.
For large inputs this will be very slow.

=head2 fermat_factor

  my @factors = fermat_factor($n);

Produces factors, not necessarily prime, of the positive number input.  The
particular algorithm is Knuth's algorithm C.  For small inputs this will be
very fast, but it slows down quite rapidly as the number of digits increases.
It is very fast for inputs with a factor close to the midpoint
(e.g. a semiprime p*q where p and q are the same number of digits).

=head2 holf_factor

  my @factors = holf_factor($n);

Produces factors, not necessarily prime, of the positive number input.  An
optional number of rounds can be given as a second parameter.  It is possible
the function will be unable to find a factor, in which case a single element,
the input, is returned.  This uses Hart's One Line Factorization with no
premultiplier.  It is an interesting alternative to Fermat's algorithm,
and there are some inputs it can rapidly factor.  In the long run it has the
same advantages and disadvantages as Fermat's method.

=head2 squfof_factor

  my @factors = squfof_factor($n);

Produces factors, not necessarily prime, of the positive number input.  An
optional number of rounds can be given as a second parameter.  It is possible
the function will be unable to find a factor, in which case a single element,
the input, is returned.  This function typically runs very fast.

=head2 prho_factor

=head2 pbrent_factor

=head2 pminus1_factor

  my @factors = prho_factor($n);

  # Use a very small number of rounds
  my @factors = prho_factor($n, 1000);

Produces factors, not necessarily prime, of the positive number input.  An
optional number of rounds can be given as a second parameter.  These attempt
to find a single factor using one of the probabilistic algorigthms of
Pollard Rho, Brent's modification of Pollard Rho, or Pollard's C<p - 1>.
These are more specialized algorithms usually used for pre-factoring very
large inputs, or checking very large inputs for naive mistakes.  If the
input is prime or they run out of rounds, they will return the single
input value.  On some inputs they will take a very long time, while on
others they succeed in a remarkably short time.



=head1 MATHEMATICAL FUNCTIONS

=head2 ExponentialIntegral

  my $Ei = ExponentialIntegral($x);

Given a non-zero floating point input C<x>, this returns the real-valued
exponential integral of C<x>, defined as the integral of C<e^t/t dt>
from C<-infinity> to C<x>.
Depending on the input, the integral is calculated using
continued fractions (C<x E<lt> -1>),
rational Chebyshev approximation (C< -1 E<lt> x E<lt> 0>),
a convergent series (small positive C<x>),
or an asymptotic divergent series (large positive C<x>).

Accuracy should be at least 14 digits.


=head2 LogarithmicIntegral

  my $li = LogarithmicIntegral($x)

Given a positive floating point input, returns the floating point logarithmic
integral of C<x>, defined as the integral of C<dt/ln t> from C<0> to C<x>.
If given a negative input, the function will croak.  The function returns
0 at C<x = 0>, and C<-infinity> at C<x = 1>.

This is often known as C<li(x)>.  A related function is the offset logarithmic
integral, sometimes known as C<Li(x)> which avoids the singularity at 1.  It
may be defined as C<Li(x) = li(x) - li(2)>.

This function is implemented as C<li(x) = Ei(ln x)> after handling special
values.

Accuracy should be at least 14 digits.


=head2 RiemannR

  my $r = RiemannR($x);

Given a positive non-zero floating point input, returns the floating
point value of Riemann's R function.  Riemann's R function gives a very close
approximation to the prime counting function.

Accuracy should be at least 14 digits.



=head1 LIMITATIONS

I have not completed testing all the functions near the word size limit
(e.g. C<2^32> for 32-bit machines).  Please report any problems you find.

The extra factoring algorithms are mildly interesting but really limited by
not being big number aware.  Assuming a desktop PC, every 32-bit number
should be factored by the main routine in a few microseconds, and 64-bit
numbers should be a few milliseconds at worst.

Perl versions earlier than 5.8.0 have issues with 64-bit.  The test suite will
try to determine if your Perl is broken.  This will show up in factoring tests.
Perl 5.6.2 32-bit works fine, as do later versions with 32-bit and 64-bit.

The module is thread-safe and should allow good concurrency.


=head1 PERFORMANCE

Counting the primes to C<10^10> (10 billion), with time in seconds.
Pi(10^10) = 455,052,511.

   External C programs in C / C++:

       1.9  primesieve 3.6 forced to use only a single thread
       2.2  yafu 1.31
       5.6  Tomás Oliveira e Silva's segmented sieve v2 (Sep 2010)
       6.6  primegen (optimized Sieve of Atkin)
      11.2  Tomás Oliveira e Silva's segmented sieve v1 (May 2003)
      17.0  Pari 2.3.5 (primepi)

   Small portable functions suitable for plugging into XS:

       5.3  My segmented SoE used in this module
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



C<is_prime>: my impressions:

   Module                    Small inputs   Large inputs (10-20dig)
   -----------------------   -------------  ----------------------
   Math::Prime::Util         Very fast      Pretty fast
   Math::Prime::XS           Very fast      Very, very slow if no small factors
   Math::Pari                Slow           OK
   Math::Prime::FastSieve    Very fast      N/A (too much memory)
   Math::Primality           Very slow      Very slow

The differences are in the implementations:
   - L<Math::Prime::FastSieve> only works in a sieved range, which is really
     fast if you can do it (M::P::U will do the same if you call
     C<prime_precalc>).  Larger inputs just need too much time and memory
     for the sieve.
   - L<Math::Primality> uses a GMP BPSW test which is overkill for our 64-bit
     range.  It's generally an order of magnitude or two slower than any
     of the others.  
   - L<Math::Pari> has some very effective code, but it has some overhead to get
     to it from Perl.  That means for small numbers it is relatively slow: an
     order of magnitude slower than M::P::XS and M::P::Util (though arguably
     this is only important for benchmarking since "slow" is ~2 microseconds).
     Large numbers transition over to smarter tests so don't slow down much.
   - L<Math::Prime::XS> does trial divisions, which is wonderful if the input
     has a small factor (or is small itself).  But it can take 1000x longer
     if given a large prime.
   - L<Math::Prime::Util> looks in the sieve for a fast bit lookup if that
     exists (default up to 30,000 but it can be expanded, e.g.
     C<prime_precalc>), uses trial division for numbers higher than this but
     not too large (0.1M on 64-bit machines, 100M on 32-bit machines), and a
     deterministic set of Miller-Rabin tests for large numbers.



Factoring performance depends on the input, and the algorithm choices used
are still being tuned.  Compared to Math::Factor::XS, it is a tiny bit faster
for most input under 10M or so, and rapidly gets faster.  For numbers
larger than 32 bits it's 10-100x faster (depending on the number -- a power
of two will be identical, while a semiprime with large factors will be on
the extreme end).  Pari's underlying algorithms and code are very
sophisticated, and will always be more so than this module, and of course
supports bignums which is a huge advantage.  Small numbers factor much, much
faster with Math::Prime::Util.  Pari passes M::P::U in speed somewhere in the
16 digit range and rapidly increases its lead.

The presentation here:
 L<http://math.boisestate.edu/~liljanab/BOISECRYPTFall09/Jacobsen.pdf>
has a lot of data on 64-bit and GMP factoring performance I collected in 2009.
Assuming you do not know anything about the inputs, trial division and
optimized Fermat work very well for small numbers (<= 10 digits), while
native SQUFOF is typically the method of choice for 11-18 digits (I've
seen claims that a lightweight QS can be faster for 15+ digits).  Some form
of Quadratic Sieve is usually used for inputs in the 19-100 digit range, and
beyond that is the Generalized Number Field Sieve.  For serious factoring,
I recommend looking info C<yafu>, C<msieve>, C<Pari>, and C<GGNFS>.



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
