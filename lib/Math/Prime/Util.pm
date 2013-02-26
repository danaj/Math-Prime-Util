package Math::Prime::Util;
use strict;
use warnings;
use Carp qw/croak confess carp/;
use Bytes::Random::Secure;

BEGIN {
  $Math::Prime::Util::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::VERSION = '0.21';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw(
                     prime_get_config prime_set_config
                     prime_precalc prime_memfree
                     is_prime is_prob_prime is_provable_prime
                     is_strong_pseudoprime is_strong_lucas_pseudoprime
                     is_aks_prime
                     miller_rabin
                     primes
                     next_prime  prev_prime
                     prime_count prime_count_lower prime_count_upper prime_count_approx
                     nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
                     random_prime random_ndigit_prime random_nbit_prime
                     random_strong_prime random_maurer_prime
                     primorial pn_primorial
                     factor all_factors
                     moebius mertens euler_phi jordan_totient exp_mangoldt
                     divisor_sum
                     ExponentialIntegral LogarithmicIntegral RiemannZeta RiemannR
                   );
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);

my %_Config;

# Similar to how boolean handles its option
sub import {
    my @options = grep $_ ne '-nobigint', @_;
    $_[0]->_import_nobigint if @options != @_;
    @_ = @options;
    goto &Exporter::import;
}

sub _import_nobigint {
  $_Config{'nobigint'} = 1;
  return unless $_Config{'xs'};
  undef *factor;        *factor          = \&_XS_factor;
  undef *is_prime;      *is_prime        = \&_XS_is_prime;
  undef *is_prob_prime; *is_prob_prime   = \&_XS_is_prob_prime;
  undef *next_prime;    *next_prime      = \&_XS_next_prime;
  undef *prev_prime;    *prev_prime      = \&_XS_prev_prime;
 #undef *prime_count;   *prime_count     = \&_XS_prime_count;
  undef *nth_prime;     *nth_prime       = \&_XS_nth_prime;
  undef *is_strong_pseudoprime;  *is_strong_pseudoprime = \&_XS_miller_rabin;
  undef *miller_rabin;  *miller_rabin    = \&_XS_miller_rabin;
  undef *moebius;       *moebius         = \&_XS_moebius;
  undef *mertens;       *mertens         = \&_XS_mertens;
  undef *euler_phi;     *euler_phi       = \&_XS_totient;
}

BEGIN {

  # Load PP code.  Nothing exported.
  require Math::Prime::Util::PP;  Math::Prime::Util::PP->import();

  eval {
    return 0 if defined $ENV{MPU_NO_XS} && $ENV{MPU_NO_XS} == 1;
    require XSLoader;
    XSLoader::load(__PACKAGE__, $Math::Prime::Util::VERSION);
    prime_precalc(0);
    $_Config{'xs'} = 1;
    $_Config{'maxbits'} = _XS_prime_maxbits();
    1;
  } or do {
    #carp "Using Pure Perl implementation: $@";

    $_Config{'xs'} = 0;
    $_Config{'maxbits'} = Math::Prime::Util::PP::_PP_prime_maxbits();

    *_prime_memfreeall = \&Math::Prime::Util::PP::_prime_memfreeall;
    *prime_memfree  = \&Math::Prime::Util::PP::prime_memfree;
    *prime_precalc  = \&Math::Prime::Util::PP::prime_precalc;

    # These probably shouldn't even be exported
    *trial_factor   = \&Math::Prime::Util::PP::trial_factor;
    *fermat_factor  = \&Math::Prime::Util::PP::fermat_factor;
    *holf_factor    = \&Math::Prime::Util::PP::holf_factor;
    *squfof_factor  = \&Math::Prime::Util::PP::squfof_factor;
    *rsqufof_factor = \&Math::Prime::Util::PP::squfof_factor;
    *pbrent_factor  = \&Math::Prime::Util::PP::pbrent_factor;
    *prho_factor    = \&Math::Prime::Util::PP::prho_factor;
    *pminus1_factor = \&Math::Prime::Util::PP::pminus1_factor;
  };

  $_Config{'nobigint'} = 0;
  $_Config{'gmp'} = 0;
  # See if they have the GMP module and haven't requested it not to be used.
  if (!defined $ENV{MPU_NO_GMP} || $ENV{MPU_NO_GMP} != 1) {
    $_Config{'gmp'} = 1 if eval { require Math::Prime::Util::GMP;
                                  Math::Prime::Util::GMP->import();
                                  1; };
  }

  # Try to figure out a system rand configuration that works for us.
  # Using something other than the craptastic system rand would be best.
  use Config;
  $_Config{'system_randbits'} = $Config{'randbits'};
  # Keep things in integer range.
  $_Config{'system_randbits'} = $_Config{'maxbits'}-1 if $_Config{'system_randbits'} >= $_Config{'maxbits'};
  # drand48 has an alternating last bit on almost every system.
  $_Config{'system_randbits'}-- if $_Config{'system_randbits'} == 48;
  no Config;

}
END {
  _prime_memfreeall;
}

if ($_Config{'maxbits'} == 32) {
  $_Config{'maxparam'}    = 4294967295;
  $_Config{'maxdigits'}   = 10;
  $_Config{'maxprime'}    = 4294967291;
  $_Config{'maxprimeidx'} = 203280221;
} else {
  $_Config{'maxparam'}    = 18446744073709551615;
  $_Config{'maxdigits'}   = 20;
  $_Config{'maxprime'}    = 18446744073709551557;
  $_Config{'maxprimeidx'} = 425656284035217743;
}
$_Config{'assume_rh'} = 0;
$_Config{'verbose'} = 0;
$_Config{'irand'} = undef;

# used for code like:
#    return _XS_foo($n)  if $n <= $_XS_MAXVAL
# which builds into one scalar whether XS is available and if we can call it.
my $_XS_MAXVAL = $_Config{'xs'}  ?  $_Config{'maxparam'}  :  -1;
my $_HAVE_GMP = $_Config{'gmp'};

# Infinity in Perl is rather O/S specific.
our $_Infinity = 0+'inf';
$_Infinity = 20**20**20 if 65535 > $_Infinity;   # E.g. Windows
our $_Neg_Infinity = -$_Infinity;

# Notes on how we're dealing with big integers:
#
#  1) if (ref($n) eq 'Math::BigInt')
#     $n is a bigint, so do bigint stuff
#
#  2) if (defined $bigint::VERSION && $n > ~0)
#     make $n into a bigint.  This is debatable, but they *did* hand us a
#     string with a big integer in it.  The big gotcha here is that
#     is_strong_lucas_pseudoprime does bigint computations, so it will load
#     up bigint and there is no way to unload it.
#
#  3) if (ref($n) =~ /^Math::Big/)
#     $n is a big int, float, or rat.  We probably want this as an int.
#
#  $n = $n->numify if $n < ~0 && ref($n) =~ /^Math::Big/;
#     get us out of big math if we can
#
# Sadly, non-modern versions of bignum (5.12.4 and earlier) completely make a
# mess of things like BigInt::numify and int(BigFloat).  Using int($x->bstr)
# seems to work.
# E.g.:
#    $n = 33662485846146713;  $n->numify;   $n is now 3.36624858461467e+16


sub prime_get_config {
  my %config = %_Config;

  $config{'precalc_to'} = ($_Config{'xs'})
                        ? _get_prime_cache_size()
                        : Math::Prime::Util::PP::_get_prime_cache_size();

  return \%config;
}

# Note: You can cause yourself pain if you turn on xs or gmp when they're not
# loaded.  Your calls will probably die horribly.
sub prime_set_config {
  my %params = (@_);  # no defaults
  while (my($param, $value) = each %params) {
    $param = lc $param;
    # dispatch table should go here.
    if      ($param eq 'xs') {
      $_Config{'xs'} = ($value) ? 1 : 0;
      $_XS_MAXVAL = $_Config{'xs'}  ?  $_Config{'maxparam'}  :  -1;
    } elsif ($param eq 'gmp') {
      $_Config{'gmp'} = ($value) ? 1 : 0;
      $_HAVE_GMP = $_Config{'gmp'};
    } elsif ($param eq 'nobigint') {
      $_Config{'nobigint'} = ($value) ? 1 : 0;
    } elsif ($param eq 'irand') {
      croak "irand must supply a sub" unless ref($value) eq 'CODE';
      $_Config{'irand'} = $value;
    } elsif ($param =~ /^(assume[_ ]?)?[ge]?rh$/ || $param =~ /riemann\s*h/) {
      $_Config{'assume_rh'} = ($value) ? 1 : 0;
    } elsif ($param eq 'verbose') {
      if    ($value =~ /^\d+$/) { }
      elsif ($value =~ /^[ty]/i) { $value = 1; }
      elsif ($value =~ /^[fn]/i) { $value = 0; }
      else { croak("Invalid setting for verbose.  0, 1, 2, etc."); }
      $_Config{'verbose'} = $value;
      _XS_set_verbose($value) if $_Config{'xs'};
      Math::Prime::Util::GMP::_GMP_set_verbose($value) if $_Config{'gmp'};
    } else {
      croak "Unknown or invalid configuration setting: $param\n";
    }
  }
  1;
}

my $_bigint_small;
sub _validate_positive_integer {
  my($n, $min, $max) = @_;
  croak "Parameter must be defined" if !defined $n;
  if (ref($n) eq 'Math::BigInt') {
    croak "Parameter '$n' must be a positive integer" unless $n->sign() eq '+';
  } else {
    croak "Parameter '$n' must be a positive integer"
          if $n eq '' || $n =~ tr/0123456789//c;
  }
  croak "Parameter '$n' must be >= $min" if defined $min && $n < $min;
  croak "Parameter '$n' must be <= $max" if defined $max && $n > $max;

  if (ref($_[0])) {
    $_[0] = Math::BigInt->new("$_[0]") unless ref($_[0]) eq 'Math::BigInt';
    $_bigint_small = Math::BigInt->new("$_Config{'maxparam'}")
                   unless defined $_bigint_small;
    if ($_[0]->bacmp($_bigint_small) <= 0) {
      $_[0] = int($_[0]->bstr);
    } else {
      $_[0]->upgrade(undef) if $_[0]->upgrade();  # Stop BigFloat upgrade
    }
  } else {
    # The second term is used instead of '<=' to fix strings like ~0+delta.
    if ( ! ($n < $_Config{'maxparam'} || int($n) eq $_Config{'maxparam'}) ) {
      # We were handed a string representing a big number.
      croak "Parameter '$n' outside of integer range" if !defined $bigint::VERSION;
      $_[0] = Math::BigInt->new("$n"); # Make $n a proper bigint object
      $_[0]->upgrade(undef) if $_[0]->upgrade();  # Stop BigFloat upgrade
    }
  }
  # One of these will be true:
  #     1) $n <= max and $n is not a bigint
  #     2) $n  > max and $n is a bigint
  1;
}

# If you use bigint then call one of the approx/bounds/math functions, you'll
# end up with full bignum turned on.  This seems non-optimal.  However, if I
# don't do this, then you'll get wrong results and end up with it turned on
# _anyway_.  As soon as anyone does something like log($n) where $n is a
# Math::BigInt, it auto-upgrade and loads up Math::BigFloat.
#
# Ideally we'd notice we were causing this, and turn off Math::BigFloat after
# we were done.
sub _upgrade_to_float {
  my($n) = @_;
  return $n unless defined $Math::BigInt::VERSION || defined $Math::BigFloat::VERSION;
  do { require Math::BigFloat; Math::BigFloat->import() }
     if defined $Math::BigInt::VERSION && !defined $Math::BigFloat::VERSION;
  return Math::BigFloat->new($n);   # $n is a Math::BigInt
}

my @_primes_small = (
   0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499);
my @_prime_count_small = (
   0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9,10,10,
   11,11,11,11,11,11,12,12,12,12,13,13,14,14,14,14,15,15,15,15,15,15,
   16,16,16,16,16,16,17,17,18,18,18,18,18,18,19);
#my @_prime_next_small = (
#   2,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,
#   29,29,29,29,29,29,31,31,37,37,37,37,37,37,41,41,41,41,43,43,47,
#   47,47,47,53,53,53,53,53,53,59,59,59,59,59,59,61,61,67,67,67,67,67,67,71);





#############################################################################

sub primes {
  my $optref = (ref $_[0] eq 'HASH')  ?  shift  :  {};
  croak "no parameters to primes" unless scalar @_ > 0;
  croak "too many parameters to primes" unless scalar @_ <= 2;
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;

  _validate_positive_integer($low);
  _validate_positive_integer($high);

  my $sref = [];
  return $sref if ($low > $high) || ($high < 2);

  if ($high > $_XS_MAXVAL) {
    if ($_HAVE_GMP) {
      $sref = Math::Prime::Util::GMP::primes($low,$high);
      if ($high > ~0) {
        # Convert the returned strings into BigInts
        croak "Internal error: large value without bigint loaded."
              unless defined $Math::BigInt::VERSION;
        @$sref = map { Math::BigInt->new("$_") } @$sref;
      } else {
        @$sref = map { int($_) } @$sref;
      }
      return $sref;
    }
    return Math::Prime::Util::PP::primes($low,$high);
  }

  my $method = $optref->{'method'};
  $method = 'Dynamic' unless defined $method;

  if ($method =~ /^(Dyn\w*|Default|Generate)$/i) {
    # Dynamic -- we should try to do something smart.

    # Tiny range?
    if (($low+1) >= $high) {
      $method = 'Trial';

    # Fast for cached sieve?
    } elsif (($high <= (65536*30)) || ($high <= _get_prime_cache_size())) {
      $method = 'Sieve';

    # At some point the segmented sieve is faster than the base sieve, not
    # to mention using much less memory.
    } elsif ($high > (1024*1024*30)) {
      $method = 'Segment';
      # The segment sieve doesn't itself use a segmented sieve for the base,
      # so it will slow down for very large endpoints (larger than 10^16).
      # Make a crude predictor of segment and trial and decide.
      if ($high > 10**14) {
        my $est_trial = ($high-$low) / 1_000_000;  # trial estimate 1s per 1M
        # segment is exponential on high, plus very fast scan.
        my $est_segment = 0.2 * 3.3**(log($high / 10**15) / log(10))
                          + ($high-$low) / 1_000_000_000_000;
        $method = 'Trial' if $est_trial <= $est_segment;
      }

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


# For random primes, there are two good papers that should be examined:
#
#  "Fast Generation of Prime Numbers and Secure Public-Key Cryptographic Parameters"
#  by Ueli M. Maurer, 1995
#  http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.26.2151
#  related discussions:
#      http://www.daimi.au.dk/~ivan/provableprimesproject.pdf
#      Handbook of Applied Cryptography by Menezes, et al.
#
#  "Close to Uniform Prime Number Generation With Fewer Random Bits"
#   by Pierre-Alain Fouque and Mehdi Tibouchi, 2011
#   http://eprint.iacr.org/2011/481
#
#
#  Some things to note:
#
#    1) Joye and Paillier have patents on their methods.  Never use them.
#
#    2) The easy-peasy method of next_prime(random number) is fast but gives
#       a terribly distribution, and not only in the obvious positive bias.
#       The probability for a prime is proportional to its gap, which is
#       really a bad distribution.
#
# In this code, for ranges within randbits (typically 48 on UNIX system rand,
# 31 for user-provided rand, and 16 for most Win32 systems), the results
# are completely uniform.  For larger ranges it is close.
#
# The random_maurer_prime function uses Maurer's FastPrime algorithm.
#
# These functions are quite fast for native size inputs, and reasonably fast
# for bigints.  Some factors that make a significant difference:
#   - Is Math::Prime::Util::GMP installed?
#   - Using Math::BigInt::GMP or Math::BigInt::Pari?  Very important.
#   - Which platform?  Typically x86_64 is best optimized.
#   - If using system rand, is RANDBITS large?
#   - What RNG?
#
# Timings using Math::BigInt::GMP, x86_64, system rand with 32+ randbits.
#
#                   random_nbit_prime         random_maurer_prime
#    n-bits       no GMP   w/ MPU::GMP        no GMP   w/ MPU::GMP
#    ----------  --------  -----------       --------  -----------
#       24-bit       25uS      same             same       same
#       64-bit       87uS      same             same       same
#      128-bit     0.032s      0.0049s         0.098s      0.056s
#      256-bit     0.062s      0.0097s         0.25s       0.15s
#      512-bit     0.13s       0.019s          0.65s       0.30s
#     1024-bit     0.28s       0.058s          1.3s        0.94s
#     2048-bit     0.91s       0.4s            3.2s        3.1s
#     4096-bit     6.6s        4.0s           23s         12.0s
#
# Writing these entirely in GMP has a problem, which is that we want to use
# a user-supplied rand function, which means a lot of callbacks.  One
# possibility is to, if they do not supply a rand function, use the GMP MT
# function with an appropriate seed.
#
# It will generate primes with more bits, but it slows down a lot.  The
# time variation becomes quite extreme once bit sizes get over 6000 or so.
#
# Random timings for 1M calls:
#   0.032   system rand
#   0.20    Math::Random::MT::Auto
#   0.96    Math::Random::Secure  irand  (with Math::Random::ISAAC::XS)
#   1.71    Math::Random::Secure         (with Math::Random::ISAAC::XS)
#   1.90   *Bytes::Random::Secure irand  (with Math::Random::ISAAC::XS)
#   2.40    Math::Random::Secure  irand
#   3.37    Math::Random::Secure
#   2.21   *Bytes::Random::Secure        (with Math::Random::ISAAC::XS)
#   3.32   *Bytes::Random::Secure irand
#   3.62   *Bytes::Random::Secure
# 123.8    *Crypt::Random         irand
# BRS: sub irand { return unpack("L", random_bytes(4)); }
#      sub rand  {return ($_[0]||1)*(unpack("L",random_bytes(4))/4294967296.0)}
# CR:  makerandom(Size=>32, Uniform=>1, Strength=>0)
#      haveged daemon running to stop /dev/random blocking
# Both BRS and CR have more features that this isn't measuring.
#
# To verify distribution:
#   perl -Iblib/lib -Iblib/arch -MMath::Prime::Util=:all -E 'my %freq; $n=1000000; $freq{random_nbit_prime(6)}++ for (1..$n); printf("%4d %6.3f%%\n", $_, 100.0*$freq{$_}/$n) for sort {$a<=>$b} keys %freq;'
#   perl -Iblib/lib -Iblib/arch -MMath::Prime::Util=:all -E 'my %freq; $n=1000000; $freq{random_prime(1260437,1260733)}++ for (1..$n); printf("%4d %6.3f%%\n", $_, 100.0*$freq{$_}/$n) for sort {$a<=>$b} keys %freq;'

{
  # These are much faster than straightforward trial division when n is big.
  # You'll want to first do a test up to and including 23.
  my @_big_gcd;
  my $_big_gcd_top = 20046;
  my $_big_gcd_use = -1;
  sub _make_big_gcds {
    croak "Internal error: make_big_gcds needs Math::BigInt!" unless defined $Math::BigInt::VERSION;
    my $p0 = primorial(Math::BigInt->new( 520));
    my $p1 = primorial(Math::BigInt->new(2052));
    my $p2 = primorial(Math::BigInt->new(6028));
    my $p3 = primorial(Math::BigInt->new($_big_gcd_top));
    $_big_gcd[0] = $p0->bdiv(223092870)->bfloor->as_int;
    $_big_gcd[1] = $p1->bdiv($p0)->bfloor->as_int;
    $_big_gcd[2] = $p2->bdiv($p1)->bfloor->as_int;
    $_big_gcd[3] = $p3->bdiv($p2)->bfloor->as_int;
  }

  # Returns a function that will get a uniform random number
  # between 0 and $max inclusive.  $max can be a bigint.
  my $_BRS;
  sub _get_rand_func {
    # First define a function $irandf that returns a 32-bit integer.  This
    # corresponds to the irand function of many CPAN modules:
    #    Math::Random::MT
    #    Math::Random::ISAAC
    #    Math::Random::Xorshift
    #    Math::Random::Secure
    # (but not Math::Random::MT::Auto which will return 64-bits)
    #
    # See if they passed one in via prime_set_config(irand=> \&irand).
    # If not, make a Bytes::Random::Secure object with non-blocking seed, and
    # use its irand method.
    #
    # This gives us a good starting point to make arbitrary size random
    # numbers.  Bytes::Random::Secure will get us excellent quality 32-bit
    # numbers on any platform, which means we can avoid possible nightmares
    # with bad system rand functions.
    my $irandf = $_Config{'irand'};
    if (!defined $irandf) {
      $_BRS = Bytes::Random::Secure->new(NonBlocking=>1) unless defined $_BRS;
      $irandf = sub { return $_BRS->irand(); };
    }
    # OK, now we have a function irandf.  Use it.
    my $randf = sub {
      my($max) = @_;
      return 0 if $max <= 0;
      my $range = $max+1;
      my $U;
      if (ref($range) eq 'Math::BigInt') {
        my $zero = $range->copy->bzero;
        my $rbits = length($range->as_bin) - 2;   # bits in range
        my $rwords = int($rbits/32) + (($rbits % 32) ? 1 : 0);
        # Generate more bits so we rarely have to loop.
        my $rmax = (($zero+2) ** ($rwords*32)) - 1;
        my $remainder = $rmax % $range;
        do {
          $U = $zero;
          $U = ($U << 32) + $irandf->()  for 1 .. $rwords;
        } while $U >= ($rmax - $remainder);
      } elsif ($range <= 4294967295) {
        my $remainder = 4294967295 % $range;
        do {
          $U = $irandf->();
        } while $U >= (4294967295 - $remainder);
      } else {
        croak "randf given max out of bounds: $max" if $range > ~0;
        my $remainder = 18446744073709551615 % $range;
        do {
          $U = ($irandf->() << 32) + $irandf->();
        } while $U >= (18446744073709551615 - $remainder);
      }
      return $U % $range;
    };
    return $randf;
  }

  # Sub to call with low and high already primes and verified range.
  my $_random_prime = sub {
    my($low,$high) = @_;
    my $prime;

    # irandf->($n) gives numbers in the range [0, $n].
    my $irandf = _get_rand_func();

    #{ my $bsize = 100; my @bins; my $counts = 10000000;
    #  for my $c (1..$counts) { $bins[ $irandf->($bsize-1) ]++; }
    #  for my $b (0..$bsize) {printf("%4d %8.5f%%\n", $b, $bins[$b]/$counts);} }

    # low and high are both primes, and low < high.

    # This is fast for small values, low memory, perfectly uniform, and consumes
    # the absolute minimum amount of randomness needed.  But it isn't feasible
    # with large values.
    if ($high <= 131072 && $high <= $_XS_MAXVAL) {
      my $li     = _XS_prime_count(2, $low);
      my $irange = _XS_prime_count($low, $high);
      my $rand = $irandf->($irange-1);
      return _XS_nth_prime($li + $rand);
    }

    $low-- if $low == 2;  # Low of 2 becomes 1 for our program.
    confess "Invalid _random_prime parameters: $low, $high" if ($low % 2) == 0 || ($high % 2) == 0;

    # We're going to look at the odd numbers only.
    #my $range = $high - $low + 1;
    my $oddrange = (($high - $low) >> 1) + 1;

    # If $low is large (e.g. >10 digits) and $range is small (say ~10k), it
    # would be fastest to call primes in the range and randomly pick one.  I'm
    # not implementing it now because it seems like a rare case.

    # If the range is reasonably small, generate using simple Monte Carlo
    # method (aka the 'trivial' method).  Completely uniform.
    if ($oddrange < $_Config{'maxparam'}) {
      my $loop_limit = 2000 * 1000;  # To protect against broken rand
      if ($low > 11) {
        while ($loop_limit-- > 0) {
          $prime = $low + 2 * $irandf->($oddrange-1);
          next if !($prime % 3) || !($prime % 5) || !($prime % 7) || !($prime % 11);
          return $prime if is_prob_prime($prime);
        }
      } else {
        while ($loop_limit-- > 0) {
          $prime = $low + 2 * $irandf->($oddrange-1);
          next if $prime > 11 && (!($prime % 3) || !($prime % 5) || !($prime % 7) || !($prime % 11));
          return 2 if $prime == 1;  # Remember the special case for 2.
          return $prime if is_prob_prime($prime);
        }
      }
      croak "Random function broken?";
    }

    # We have an ocean of range, and a teaspoon to hold randomness.

    # Since we have an arbitrary range and not a power of two, I don't see how
    # Fouque's algorithm A1 could be used (where we generate lower bits and
    # generate random sets of upper).  Similarly trying to simply generate
    # upper bits is full of ways to trip up and get non-uniform results.
    #
    # What I'm doing here is:
    #
    #   1) divide the range into semi-evenly sized partitions, where each part
    #      is as close to $rand_max_val as we can.
    #   2) randomly select one of the partitions.
    #   3) iterate choosing random values within the partition.
    #
    # The downside is that we're skewing a _lot_ farther from uniformity than
    # we'd like.  Imagine we started at 0 with 1e18 partitions of size 100k each.
    # Probability of '5' being returned =
    #   1.04e-22 = 1e-18 (chose first partition) * 1/9592 (chose '5')
    # Probability of '100003' being returned =
    #   1.19e-22 = 1e-18 (chose second partition) * 1/8392 (chose '100003')
    # Probability of '99999999999999999999977' being returned =
    #   5.20e-22 = 1e-18 (chose last partition)  *  1/1922 (chose '99...77')
    # So the primes in the last partition will show up 5x more often.
    # The partitions are selected uniformly, and the primes within are selected
    # uniformly, but the number of primes in each bucket is _not_ uniform.
    # Their individual probability of being selected is the probability of the
    # partition (uniform) times the probability of being selected inside the
    # partition (uniform with respect to all other primes in the same
    # partition, but each partition is different and skewed).
    #
    # When selecting n-bit or n-digit primes, this effect is _much_ smaller, as
    # the skew becomes approx lg(2^n) / lg(2^(n-1)) which is pretty close to 1.
    # Note that we really want big partitions to even out any local skews, which
    # worries me on systems with randbits of 16.  In fact, we should probably
    # just get two numbers on those systems.
    #
    # Another idea I'd like to try sometime is:
    #  pclo = prime_count_lower(low);
    #  pchi = prime_count_upper(high);
    #  do {
    #    $nth = random selection between pclo and pchi
    #    $prguess = nth_prime_approx($nth);
    #  } while ($prguess >= low) && ($prguess <= high);
    #  monte carlo select a prime in $prguess-2**24 to $prguess+2**24
    # which accounts for the prime distribution.

    my($binsize, $nparts);
    my $rand_part_size = 1 << (($_Config{'maxbits'} > 32) ? 32 : 31);
    if (ref($oddrange) eq 'Math::BigInt') {
      # Go to some trouble here because some systems are wonky, such as
      # giving us +a/+b = -r.  Also note the quotes for the bigint argument.
      # Without that, Math::BigInt::GMP on 32-bit Win32 will return garbage.
      my($nbins, $rem);
      ($nbins, $rem) = $oddrange->copy->bdiv( "$rand_part_size" );
      $nbins++ if $rem > 0;
      $nbins = $nbins->as_int();
      ($binsize,$rem) = $oddrange->copy->bdiv($nbins);
      $binsize++ if $rem > 0;
      $binsize = $binsize->as_int();
      $nparts  = $oddrange->copy->bdiv($binsize)->as_int();
      $low = $high->copy->bzero->badd($low) if ref($low) ne 'Math::BigInt';
    } else {
      my $nbins = int($oddrange / $rand_part_size);
      $nbins++ if $nbins * $rand_part_size != $oddrange;
      $binsize = int($oddrange / $nbins);
      $binsize++ if $binsize * $nbins != $oddrange;
      $nparts = int($oddrange/$binsize);
    }
    $nparts-- if ($nparts * $binsize) == $oddrange;

    my $rpart = $irandf->($nparts);

    my $primelow = $low + 2 * $binsize * $rpart;
    my $partsize = ($rpart < $nparts) ? $binsize
                                      : $oddrange - ($nparts * $binsize);
    $partsize = int($partsize->bstr) if ref($partsize) eq 'Math::BigInt';
    #warn "range $oddrange  = $nparts * $binsize + ", $oddrange - ($nparts * $binsize), "\n";
    #warn "  chose part $rpart size $partsize\n";
    #warn "  primelow is $low + 2 * $binsize * $rpart = $primelow\n";
    #die "Result could be too large" if ($primelow + 2*($partsize-1)) > $high;

    # Generate random numbers in the interval until one is prime.
    my $loop_limit = 2000 * 1000;  # To protect against broken rand

    # Simply things for non-bigints.
    if (ref($low) ne 'Math::BigInt') {
      while ($loop_limit-- > 0) {
        my $rand = $irandf->($partsize-1);
        $prime = $primelow + $rand + $rand;
        croak "random prime failure, $prime > $high" if $prime > $high;
        if ($prime <= 23) {
          $prime = 2 if $prime == 1;  # special case for low = 2
          next unless (0,0,1,1,0,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,1)[$prime];
          return $prime;
        }
        next if !($prime % 3) || !($prime % 5) || !($prime % 7) || !($prime % 11);
        # It looks promising.  Check it.
        next unless is_prob_prime($prime);
        return $prime;
      }
      croak "Random function broken?";
    }

    # By checking a wheel 30 mod, we can skip anything that would be a multiple
    # of 2, 3, or 5, without even having to create the bigint prime.
    my @w30 = (1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0);
    my $primelow30 = $primelow % 30;
    $primelow30 = int($primelow30->bstr) if ref($primelow30) eq 'Math::BigInt';

    # Big GCD's are hugely fast with GMP or Pari, but super slow with Calc.
    if ($_big_gcd_use < 0) {
      $_big_gcd_use = 0;
      my $lib = Math::BigInt->config()->{lib};
      $_big_gcd_use = 1 if $lib =~ /^Math::BigInt::(GMP|Pari)/;
      _make_big_gcds() if $_big_gcd_use;
    }

    while ($loop_limit-- > 0) {
      my $rand = $irandf->($partsize-1);
      # Check wheel-30 mod
      my $rand30 = $rand % 30;
      next if $w30[($primelow30 + 2*$rand30) % 30]
              && ($rand > 3 || $primelow > 5);
      # Construct prime
      $prime = $primelow + $rand + $rand;
      croak "random prime failure, $prime > $high" if $prime > $high;
      if ($prime <= 23) {
        $prime = 2 if $prime == 1;  # special case for low = 2
        next unless (0,0,1,1,0,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,1)[$prime];
        return $prime;
      }
      # Perform quick trial division
      next unless Math::BigInt::bgcd($prime, 7436429) == 1;
      if ($_big_gcd_use && $prime > $_big_gcd_top) {
        next unless Math::BigInt::bgcd($prime, $_big_gcd[0]) == 1;
        next unless Math::BigInt::bgcd($prime, $_big_gcd[1]) == 1;
        next unless Math::BigInt::bgcd($prime, $_big_gcd[2]) == 1;
        next unless Math::BigInt::bgcd($prime, $_big_gcd[3]) == 1;
      }
      # It looks promising.  Check it.
      next unless is_prob_prime($prime);
      return $prime;
    }
    croak "Random function broken?";
  };

  # Cache of tight bounds for each digit.  Helps performance a lot.
  my @_random_ndigit_ranges = (undef, [2,7], [11,97] );
  my @_random_nbit_ranges   = (undef, undef, [2,3],[5,7] );

  sub random_prime {
    my $low = (@_ == 2)  ?  shift  :  2;
    my $high = shift;
    _validate_positive_integer($low);
    _validate_positive_integer($high);

    # Tighten the range to the nearest prime.
    $low = ($low <= 2)  ?  2  :  next_prime($low-1);
    $high = ($high < ~0)  ?  prev_prime($high + 1)  :  prev_prime($high);
    return $low if ($low == $high) && is_prob_prime($low);
    return if $low >= $high;

    # At this point low and high are both primes, and low < high.
    return $_random_prime->($low, $high);
  }

  sub random_ndigit_prime {
    my($digits) = @_;
    _validate_positive_integer($digits, 1);

    my $bigdigits = $digits >= $_Config{'maxdigits'};
    croak "Large random primes not supported on old Perl" if $] < 5.008 && $_Config{'maxbits'} > 32 && !$bigdigits && $digits > 15;
    if ($bigdigits && $_Config{'nobigint'}) {
      _validate_positive_integer($digits, 1, $_Config{'maxdigits'});
      # Special case for nobigint and threshold digits
      if (!defined $_random_ndigit_ranges[$digits]) {
        my $low  = int(10 ** ($digits-1));
        my $high = ~0;
        $_random_ndigit_ranges[$digits] = [next_prime($low),prev_prime($high)];
      }
    }

    if (!defined $_random_ndigit_ranges[$digits]) {
      if ($bigdigits) {
        if (!defined $Math::BigInt::VERSION) {
          eval { require Math::BigInt; Math::BigInt->import(try=>'GMP,Pari'); 1; }
          or do { croak "Cannot load Math::BigInt"; };
        }
        my $low  = Math::BigInt->new('10')->bpow($digits-1);
        my $high = Math::BigInt->new('10')->bpow($digits);
        # Just pull the range in to the nearest odd.
        $_random_ndigit_ranges[$digits] = [$low+1, $high-1];
      } else {
        my $low  = int(10 ** ($digits-1));
        my $high = int(10 ** $digits);
        # Note: Perl 5.6.2 cannot represent 10**15 as an integer, so things
        # will crash all over the place if you try.  We can stringify it, but
        # will just fail tests later.
        $_random_ndigit_ranges[$digits] = [next_prime($low),prev_prime($high)];
      }
    }
    my ($low, $high) = @{$_random_ndigit_ranges[$digits]};
    return $_random_prime->($low, $high);
  }

  sub random_nbit_prime {
    my($bits) = @_;
    _validate_positive_integer($bits, 2);

    if (!defined $_random_nbit_ranges[$bits]) {
      my $bigbits = $bits > $_Config{'maxbits'};
      croak "Large random primes not supported on old Perl" if $] < 5.008 && $_Config{'maxbits'} > 32 && !$bigbits && $bits > 49;
      if ($bigbits) {
        if (!defined $Math::BigInt::VERSION) {
          eval { require Math::BigInt; Math::BigInt->import(try=>'GMP,Pari'); 1; }
          or do { croak "Cannot load Math::BigInt"; };
        }
        my $low  = Math::BigInt->new('2')->bpow($bits-1);
        my $high = Math::BigInt->new('2')->bpow($bits);
        # Don't pull the range in to primes, just odds
        $_random_nbit_ranges[$bits] = [$low+1, $high-1];
      } else {
        my $low  = 1 << ($bits-1);
        my $high = ($bits == $_Config{'maxbits'})
                   ? ~0-1
                   : ~0 >> ($_Config{'maxbits'} - $bits);
        $_random_nbit_ranges[$bits] = [next_prime($low-1),prev_prime($high+1)];
        # Example: bits = 7.
        #    low = 1<<6 = 64.            next_prime(64-1)  = 67
        #    high = ~0 >> (64-7) = 127.  prev_prime(127+1) = 127
      }
    }
    my ($low, $high) = @{$_random_nbit_ranges[$bits]};
    return $_random_prime->($low, $high);
  }

  sub random_maurer_prime {
    my($k) = @_;
    _validate_positive_integer($k, 2);
    if ($] < 5.008 && $_Config{'maxbits'} > 32) {
      return random_nbit_prime($k) if $k <= 49;
      croak "Random Maurer not supported on old Perl";
    }

    # Results for random_nbit_prime are proven for all native bit sizes.  We
    # could go even higher if we used is_provable_prime or looked for is_prime
    # returning 2.  This should be reasonably fast to ~128 bits with MPU::GMP.
    my $p0 = $_Config{'maxbits'};

    return random_nbit_prime($k) if $k <= $p0;

    if (!defined $Math::BigInt::VERSION) {
      eval { require Math::BigInt; Math::BigInt->import(try=>'GMP,Pari'); 1; }
      or do { croak "Cannot load Math::BigInt"; };
    }
    if (!defined $Math::BigFloat::VERSION) {
      eval { require Math::BigFloat; Math::BigFloat->import(); 1; }
      or do { croak "Cannot load Math::BigFloat"; };
    }

    # Set verbose to 3 to get pretty output like Crypt::Primes
    my $verbose = $_Config{'verbose'};
    local $| = 1 if $verbose > 2;

    my $c = Math::BigFloat->new("0.065"); # higher = more trial divisions
    my $r = Math::BigFloat->new("0.5");   # relative size of the prime q
    my $m = 20;                           # makes sure R is big enough
    my $B = ($c * $k * $k)->bfloor;
    my $irandf = _get_rand_func();

    # Generate a random prime q of size $r*$k, where $r >= 0.5.  Try to
    # cleverly select r to match the size of a typical random factor.
    if ($k > 2*$m) {
      do {
        my $s = Math::BigFloat->new($irandf->(2147483647))->bdiv(2147483648);
        $r = Math::BigFloat->new(2)->bpow($s-1);
      } while ($k*$r >= $k-$m);
    }

    # I've seen +0, +1, and +2 here.  Maurer uses +0.  Menezes uses +1.
    my $q = random_maurer_prime( ($r * $k)->bfloor + 1 );
    $q = Math::BigInt->new("$q") unless ref($q) eq 'Math::BigInt';
    my $I = Math::BigInt->new(2)->bpow($k-2)->bdiv($q)->bfloor->as_int();
    print "B = $B  r = $r  k = $k  q = $q  I = $I\n" if $verbose && $verbose != 3;

    # Big GCD's are hugely fast with GMP or Pari, but super slow with Calc.
    if ($_big_gcd_use < 0) {
      $_big_gcd_use = 0;
      my $lib = Math::BigInt->config()->{lib};
      $_big_gcd_use = 1 if $lib =~ /^Math::BigInt::(GMP|Pari)/;
      _make_big_gcds() if $_big_gcd_use;
    }

    my $loop_limit = 1_000_000 + $k * 1_000;
    while ($loop_limit-- > 0) {
      # R is a random number between $I+1 and 2*$I
      my $R = $I + 1 + $irandf->( $I - 1 );
      #my $n = 2 * $R * $q + 1;
      my $n = Math::BigInt->new(2)->bmul($R)->bmul($q)->badd(1);
      # We constructed a promising looking $n.  Now test it.
      print "." if $verbose > 2;

      # Trial divisions, trying to quickly weed out non-primes.
      next unless Math::BigInt::bgcd($n, 111546435) == 1;
      if ($_big_gcd_use && $n > $_big_gcd_top) {
        next unless Math::BigInt::bgcd($n, $_big_gcd[0]) == 1;
        next unless Math::BigInt::bgcd($n, $_big_gcd[1]) == 1;
        next unless Math::BigInt::bgcd($n, $_big_gcd[2]) == 1;
        next unless Math::BigInt::bgcd($n, $_big_gcd[3]) == 1;
      }
      print "+" if $verbose > 2;
      if ($_HAVE_GMP) {
        next unless Math::Prime::Util::GMP::is_strong_pseudoprime($n, 2);
      }
      print "*" if $verbose > 2;

      # Now we do Lemma 1 -- a special case of the Pocklington test.
      # Let F = q where q is prime, and n = 2RF+1.
      # If F > sqrt(n) or F odd and F > R, and a^((n-1)/F)-1 mod n = 1, n prime.

      # Based on our construction, this should always be true.  Check anyway.
      next unless $q > $R;

      # Select random 'a' values.  If n is prime, then almost any 'a' value
      # will work, so just try two small ones instead of generating a giant
      # random 'a' between 2 and n-2.  This makes the powmods run faster.
      foreach my $try_a (2, 7) {
        # my $a = 2 + $irandf->( $n - 4 );
        my $a = Math::BigInt->new($try_a);
        my $b = $a->copy->bmodpow($n-1, $n);
        next unless $b == 1;

        # Now do the one gcd check we need to do.
        $b = $a->copy->bmodpow(2*$R, $n);
        next unless Math::BigInt::bgcd($b-1, $n) == 1;
        if ($verbose) {
          print "", ($verbose == 3) ? "($k)" : "$n passed final gcd\n";
        }

        # Instead of the previous gcd, we could check q >= n**1/3 and also do
        # some tests on x & y from 2R = xq+y (see Lemma 2 from Maurer's paper).
        # Crypt::Primes does the q test but doesn't do the x/y tests.
        #   next if ($q <= $n->copy->broot(3));
        #   my $x = (2*$R)->bdiv($q)->bfloor;
        #   my $y = 2*$R - $x*$q;
        #   my $z = $y*$y - 4*$x;
        #   next if $z == 0;
        #   next if $z->bsqrt->bfloor->bpow(2) == $z;  # perfect square
        # Menezes seems to imply only the q test needs to be done, but this
        # doesn't follow from Lemma 2.  Also note the entire POINT of going to
        # Lemma 2 is that we now allow r to be 0.3334, making q smaller.  If we
        # run this without changing r, then x will typically be 0 and this fails.

        # Verify with a BPSW test on the result.  This could:
        #  1) save us from accidently outputing a non-prime due to some mistake
        #  2) make history by finding the first known BPSW pseudo-prime
        croak "Maurer prime $n=2*$R*$q+1 failed BPSW" unless is_prob_prime($n);

        return $n;
      }
      # Didn't pass the selected a values.  Try another R.
    }
    croak "Failure in random_maurer_prime, could not find a prime\n";
  } # End of random_maurer_prime

  # Gordon's algorithm for generating a strong prime.
  sub random_strong_prime {
    my($t) = @_;
    _validate_positive_integer($t, 128);
    croak "Random strong primes must be >= 173 bits on old Perl" if $] < 5.008 && $_Config{'maxbits'} > 32 && $t < 173;

    if (!defined $Math::BigInt::VERSION) {
      eval { require Math::BigInt; Math::BigInt->import(try=>'GMP,Pari'); 1; }
      or do { croak "Cannot load Math::BigInt"; };
    }
    my $irandf = _get_rand_func();

    my $l   = (($t+1) >> 1) - 2;
    my $lp  = int($t/2) - 20;
    my $lpp = $l - 20;
    while (1) {
      my $qp  = random_nbit_prime($lp);
      my $qpp = random_nbit_prime($lpp);
      $qp  = Math::BigInt->new("$qp")  unless ref($qp)  eq 'Math::BigInt';
      $qpp = Math::BigInt->new("$qpp") unless ref($qpp) eq 'Math::BigInt';
      my ($il, $rem) = Math::BigInt->new(2)->bpow($l-1)->bsub(1)->bdiv(2*$qpp);
      $il++ if $rem > 0;
      $il = $il->as_int();
      my $iu = Math::BigInt->new(2)->bpow($l)->bsub(2)->bdiv(2*$qpp)->as_int();
      my $istart = $il + $irandf->($iu - $il);
      for (my $i = $istart; $i <= $iu; $i++) {  # Search for q
        my $q = 2 * $i * $qpp + 1;
        next unless is_prob_prime($q);
        my $pp = $qp->copy->bmodpow($q-2, $q)->bmul(2)->bmul($qp)->bsub(1);
        my ($jl, $rem) = Math::BigInt->new(2)->bpow($t-1)->bsub($pp)->bdiv(2*$q*$qp);
        $jl++ if $rem > 0;
        $jl = $jl->as_int();
        my $ju = Math::BigInt->new(2)->bpow($t)->bsub(1)->bsub($pp)->bdiv(2*$q*$qp)->as_int();
        my $jstart = $jl + $irandf->($ju - $jl);
        for (my $j = $jstart; $j <= $ju; $j++) {  # Search for p
          my $p = $pp + 2 * $j * $q * $qp;
          return $p if is_prob_prime($p);
        }
      }
    }
  }

} # end of the random prime section

sub primorial {
  my($n) = @_;
  _validate_positive_integer($n);

  my $pn = 1;
  if ($n >= (($_Config{'maxbits'} == 32) ? 29 : 53)) {
    if (!defined $Math::BigInt::VERSION) {
      eval { require Math::BigInt; Math::BigInt->import(try=>'GMP,Pari'); 1; }
      or do { croak "Cannot load Math::BigInt"; };
    }
    $pn = Math::BigInt->bone();
  }
  # Make sure we use their type if they passed one in.
  $pn = $_[0]->copy->bone() if ref($_[0]) eq 'Math::BigInt';

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::primorial) {
    if (ref($pn) eq 'Math::BigInt') {
      $pn->bzero->badd( Math::Prime::Util::GMP::primorial($n) );
    } else {
      $pn = int( Math::Prime::Util::GMP::primorial($n) );
    }
  } else {
    foreach my $p ( @{ primes($n) } ) {
      $pn *= $p;
    }
  }
  return $pn;
}

sub pn_primorial {
  my($n) = @_;
  return primorial(nth_prime($n));
}


sub all_factors {
  my $n = shift;
  my $use_bigint = defined $bigint::VERSION
             || !($n < $_Config{'maxparam'} || int($n) eq $_Config{'maxparam'});
  my @factors = factor($n);
  my %all_factors;
  if ($use_bigint) {
    foreach my $f1 (@factors) {
      next if $f1 >= $n;
      my $big_f1 = Math::BigInt->new("$f1");
      my @to_add = map { ($_ <= ~0) ? int($_->bstr) : $_ }
                   grep { $_ < $n }
                   map { $big_f1 * $_ }
                   keys %all_factors;
      undef @all_factors{ $f1, @to_add };
    }
  } else {
    foreach my $f1 (@factors) {
      next if $f1 >= $n;
      my @to_add = grep { $_ < $n }
                   map { int($f1 * $_) }
                   keys %all_factors;
      undef @all_factors{ $f1, @to_add };
    }
  }
  @factors = sort {$a<=>$b} keys %all_factors;
  return @factors;
}


# A008683 Moebius function mu(n)
# A030059, A013929, A030229, A002321, A005117, A013929 all relate.
sub moebius {
  my($n, $nend) = @_;
  _validate_positive_integer($n);

  if (defined $nend) {
    _validate_positive_integer($nend);
    return if $nend < $n;
  } else {
    $nend = $n;
  }
  return _XS_moebius($n, $nend) if $nend <= $_XS_MAXVAL;

  # Moebius over a range.
  if ($nend != $n) {
    my ($lo,$hi) = ($n,$nend);
    my @mu = map { 1 } $lo .. $hi;
    $mu[0] = 0 if $lo == 0;
    my $sqrtn = int(sqrt($hi)+0.5);
    foreach my $p ( @{ primes($sqrtn) } ) {
      my $i = $p * $p;
      $i = $i * int($lo/$i) + (($lo % $i)  ? $i : 0)  if $i < $lo;
      while ($i <= $hi) {
        $mu[$i-$lo] = 0;
        $i += $p * $p;
      }
      $i = $p;
      $i = $i * int($lo/$i) + (($lo % $i)  ? $i : 0)  if $i < $lo;
      while ($i <= $hi) {
        $mu[$i-$lo] *= -$p;
        $i += $p;
      }
    }
    foreach my $i ($lo .. $hi) {
      my $m = $mu[$i-$lo];
      $m *= -1 if abs($m) != $i;
      $mu[$i-$lo] = ($m>0) - ($m<0);
    }
    return @mu;
  }

  return $n if $n <= 1;
  # Quick check for small replicated factors
  return 0 if ($n >= 25) && (!($n % 4) || !($n % 9) || !($n % 25));
  my @factors = factor($n);
  my $lastf = 0;
  foreach my $factor (@factors) {
    return 0 if $factor == $lastf;
    $lastf = $factor;
  }
  return (((scalar @factors) % 2) == 0) ? 1 : -1;
}

# A002321 Mertens' function.  mertens(n) = sum(moebius(1,n))
sub mertens {
  my($n) = @_;
  _validate_positive_integer($n);
  return _XS_mertens($n) if $n <= $_XS_MAXVAL;
  # This is the most basic Deléglise and Rivat algorithm.  u = n^1/2
  # and no segmenting is done.  Their algorithm uses u = n^1/3, breaks
  # the summation into two parts, and calculates those in segments.  Their
  # computation time growth is half of this code.
  return $n if $n <= 1;
  my $u = int(sqrt($n));
  my @mu = (0, moebius(1, $u));          # Hold values of mu for 0-u
  my $musum = 0;
  my @M = map { $musum += $_; } @mu;     # Hold values of M for 0-u
  my $sum = $M[$u];
  foreach my $m (1 .. $u) {
    next if $mu[$m] == 0;
    my $inner_sum = 0;
    my $lower = int($u/$m) + 1;
    my $last_nmk = int($n/($m*$lower));
    my ($this_k, $next_k) = (0, int($n/($m*1)));
    for my $nmk (1 .. $last_nmk) {
      $this_k = $next_k;
      $next_k = int($n/($m*($nmk+1)));
      $inner_sum += $M[$nmk] * ($this_k - $next_k);
    }
    $sum -= $mu[$m] * $inner_sum;
  }
  return $sum;
}


# A000010 Euler Phi, aka Euler Totient
sub euler_phi {
  my($n, $nend) = @_;
  # SAGE defines this to be 0 for all n <= 0.  Others choose differently.
  # I am following SAGE's decision for n <= 0.
  return 0 if defined $n && $n < 0;
  _validate_positive_integer($n);
  if (defined $nend) {
    _validate_positive_integer($nend);
    return if $nend < $n;
  } else {
    $nend = $n;
  }
  return _XS_totient($n, $nend) if $nend <= $_XS_MAXVAL;

  # Totient over a range.  Could be improved, as this can use huge memory.
  if ($nend != $n) {
    return () if $nend < $n;
    my @totients = (0 .. $nend);
    foreach my $i (2 .. $nend) {
      next unless $totients[$i] == $i;
      $totients[$i] = $i-1;
      foreach my $j (2 .. int($nend / $i)) {
        $totients[$i*$j] = ($totients[$i*$j]*($i-1))/$i;
      }
    }
    splice(@totients, 0, $n) if $n > 0;
    return @totients;
  }

  return $n if $n <= 1;
  my %factor_mult;
  my @factors = grep { !$factor_mult{$_}++ } factor($n);

  if (ref($n) ne 'Math::BigInt') {
    my $totient = 1;
    foreach my $factor (@factors) {
      $totient *= ($factor - 1);
      $totient *= $factor for (2 .. $factor_mult{$factor});
    }
    return $totient;
  }

  # Some real wackiness to solve issues with Math::BigInt::GMP (not seen with
  # Pari or Calc).  Results of the multiply will go negative if we don't do
  # this.  To see if you hit the standalone bug:
  #      perl -E 'my $a = 2931542417; use bigint lib=>'GMP'; my $n = 49754396241690624; my $x = $n*$a; say $x;'
  # This may be related to RT 71548 of Math::BigInt::GMP.
  my $totient = $n->copy->bone;
  foreach my $factor (@factors) {
    my $f = $n->copy->bzero->badd("$factor");
    $totient->bmul($f->copy->bsub(1));
    $totient->bmul($f)  for (2 .. $factor_mult{$factor});
  }
  return $totient;
}

# Jordan's totient -- a generalization of Euler's totient.
sub jordan_totient {
  my($k, $n) = @_;
  _validate_positive_integer($k, 1);
  return euler_phi($n) if $k == 1;

  return 0 if defined $n && $n <= 0;  # Following SAGE's logic here.
  _validate_positive_integer($n);
  return 1 if $n <= 1;

  my %factor_mult;
  my @factors = grep { !$factor_mult{$_}++ }
                ($n <= $_XS_MAXVAL) ? _XS_factor($n) : factor($n);

  my $totient = $n - $n + 1;

  if (ref($n) ne 'Math::BigInt') {
    foreach my $factor (@factors) {
      my $fmult = int($factor ** $k);
      $totient *= ($fmult - 1);
      $totient *= $fmult for (2 .. $factor_mult{$factor});
    }
  } else {
    my $zero = $n->copy->bzero;
    foreach my $factor (@factors) {
      my $fmult = $zero->copy->badd("$factor")->bpow($k);
      $totient->bmul($fmult->copy->bsub(1));
      $totient->bmul($fmult) for (2 .. $factor_mult{$factor});
    }
  }
  return $totient;
}

# Mathematica and Pari both have functions like this.
sub divisor_sum {
  my($n, $sub) = @_;
  return 0 if defined $n && $n < 1;
  _validate_positive_integer($n);

  if (!defined $sub) {
    return $n if $n == 1;
    my $sum = 1 + $n;
    foreach my $f ( all_factors($n) ) {
      $sum += $f;
    }
    return $sum;
  }

  croak "Second argument must be a code ref" unless ref($sub) eq 'CODE';
  my $sum = $sub->(1);
  return $sum if $n == 1;
  foreach my $f (all_factors($n), $n ) {
    $sum += $sub->($f);
  }
  return $sum;
}

# Omega function A001221.  Just an example.
sub _omega {
  my($n) = @_;
  return 0 if defined $n && $n <= 1;
  _validate_positive_integer($n);
  my %factor_mult;
  my @factors = grep { !$factor_mult{$_}++ } factor($n);
  return scalar @factors;
}

# Exponential of Mangoldt function (A014963).
# Return p if n = p^m [p prime, m >= 1], 1 otherwise.
sub exp_mangoldt {
  my($n) = @_;
  return 1 if defined $n && $n <= 1;
  _validate_positive_integer($n);

  #my $is_prime = ($n<=$_XS_MAXVAL) ? _XS_is_prob_prime($n) : is_prob_prime($n);
  #return $n if $is_prime;

  my %factor_mult;
  my @factors = grep { !$factor_mult{$_}++ }
                ($n <= $_XS_MAXVAL) ? _XS_factor($n) : factor($n);

  return 1 unless scalar @factors == 1;
  return $factors[0];
}


#############################################################################
# Front ends to functions.
#
# These will do input validation, then call the appropriate internal function
# based on the input (XS, GMP, PP).
#############################################################################

# Doing a sub here like:
#
#   sub foo {  my($n) = @_;  _validate_positive_integer($n);
#              return _XS_... if $_Config{'xs'} && $n <= $_Config{'maxparam'}; }
#
# takes about 0.7uS on my machine.  Operations like is_prime and factor run
# on small input (under 100_000) typically take a lot less time than this.  So
# the overhead for these is significantly more than just the XS call itself.
#
# The plan for some of these functions will be to invert the operation.  That
# is, the XS functions will look at the input and make a call here if the input
# is large.

sub is_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_positive_integer($n);

  return _XS_is_prime($n) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_prime($n) if $_HAVE_GMP;
  return is_prob_prime($n);
}

sub is_aks_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_positive_integer($n);

  return _XS_is_aks_prime($n) if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_aks_prime($n) if $_HAVE_GMP
                       && defined &Math::Prime::Util::GMP::is_aks_prime;
  return Math::Prime::Util::PP::is_aks_prime($n);
}


sub next_prime {
  my($n) = @_;
  _validate_positive_integer($n);

  # If we have XS and n is either small or bigint is unknown, then use XS.
  return _XS_next_prime($n) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL
             && (!defined $bigint::VERSION || $n < $_Config{'maxprime'} );

  # Try to stick to the plan with respect to maximum return values.
  return 0 if ref($_[0]) ne 'Math::BigInt' && $n >= $_Config{'maxprime'};

  if ($_HAVE_GMP) {
    # If $n is a bigint object, try to make the return value the same
    return (ref($_[0]) eq 'Math::BigInt')
        ?  $_[0]->copy->bzero->badd(Math::Prime::Util::GMP::next_prime($n))
        :  Math::Prime::Util::GMP::next_prime($n);
  }
  return Math::Prime::Util::PP::next_prime($n);
}

sub prev_prime {
  my($n) = @_;
  _validate_positive_integer($n);

  return _XS_prev_prime($n) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL;
  if ($_HAVE_GMP) {
    # If $n is a bigint object, try to make the return value the same
    return (ref($n) eq 'Math::BigInt')
        ?  $n->copy->bzero->badd(Math::Prime::Util::GMP::prev_prime($n))
        :  Math::Prime::Util::GMP::prev_prime($n);
  }
  return Math::Prime::Util::PP::prev_prime($n);
}

sub prime_count {
  my($low,$high) = @_;
  if (defined $high) {
    _validate_positive_integer($low);
    _validate_positive_integer($high);
  } else {
    ($low,$high) = (2, $low);
    _validate_positive_integer($high);
  }
  return 0 if $high < 2  ||  $low > $high;

  if ($high <= $_XS_MAXVAL) {
    if ($high > 4_000_000) {
      # These estimates need a lot of work.
      #my $est_segment = 10.0 * 1.5**(log($high / 10**16) / log(10))
      #                  + (($high-$low)/10**12);
      #my $est_lehmer = 0.0000000057 * $high**0.72
      #                 + 0.0000000057 * $low**0.72;
      #if ($est_lehmer < $est_segment) {
      if ( ($high / ($high-$low+1)) < 100 ) {
        $low = 2 if $low < 2;
        return _XS_lehmer_pi($high) - _XS_lehmer_pi($low-1);
      }
    }
    return _XS_prime_count($low,$high);
  }
  # We can relax these constraints if MPU::GMP gets a Lehmer implementation.
  return Math::Prime::Util::GMP::prime_count($low,$high) if $_HAVE_GMP
                       && defined &Math::Prime::Util::GMP::prime_count
                       && (   (ref($high) eq 'Math::BigInt')
                           || (($high-$low) < int($low/1_000_000))
                          );
  return Math::Prime::Util::PP::prime_count($low,$high);
}

sub nth_prime {
  my($n) = @_;
  _validate_positive_integer($n);

  return _XS_nth_prime($n) if $_Config{'xs'} && $n <= $_Config{'maxprimeidx'};
  return Math::Prime::Util::PP::nth_prime($n);
}

sub factor {
  my($n) = @_;
  _validate_positive_integer($n);

  return _XS_factor($n) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL;

  if ($_HAVE_GMP) {
    my @factors = Math::Prime::Util::GMP::factor($n);
    if (ref($n) eq 'Math::BigInt') {
      @factors = map { ($_ > ~0) ? $n->copy->bzero->badd($_) : $_ } @factors;
    }
    return @factors;
  }

  return Math::Prime::Util::PP::factor($n);
}

sub is_strong_pseudoprime {
  my($n) = shift;
  _validate_positive_integer($n);
  # validate bases?
  return _XS_miller_rabin($n, @_) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_strong_pseudoprime($n, @_) if $_HAVE_GMP;
  return Math::Prime::Util::PP::miller_rabin($n, @_);
}

sub is_strong_lucas_pseudoprime {
  my($n) = shift;
  _validate_positive_integer($n);
  return Math::Prime::Util::GMP::is_strong_lucas_pseudoprime("$n") if $_HAVE_GMP;
  return Math::Prime::Util::PP::is_strong_lucas_pseudoprime($n);
}

sub miller_rabin {
  #warn "miller_rabin() is deprecated. Use is_strong_pseudoprime instead.";
  return is_strong_pseudoprime(@_);
}

#############################################################################

  # Oct 2012 note:  these numbers are old.
  #
  # Timings for various combinations, given the current possibilities of:
  #    1) XS MR optimized (either x86-64, 32-bit on 64-bit mach, or half-word)
  #    2) XS MR non-optimized (big input not on 64-bit machine)
  #    3) PP MR with small input (non-bigint Perl)
  #    4) PP MR with large input (using functions for mulmod)
  #    5) PP MR with full bigints
  #    6) PP Lucas with small input
  #    7) PP Lucas with large input
  #    8) PP Lucas with full bigints
  #
  # Time for one test:
  #       0.5uS  XS MR with small input
  #       0.8uS  XS MR with large input
  #       7uS    PP MR with small input
  #     400uS    PP MR with large input
  #    5000uS    PP MR with bigint
  #    2700uS    PP LP with small input
  #    6100uS    PP LP with large input
  #    7400uS    PP LP with bigint

sub is_prob_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_positive_integer($n);

  return _XS_is_prob_prime($n) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_prob_prime($n) if $_HAVE_GMP;

  return 2 if $n == 2 || $n == 3 || $n == 5 || $n == 7;
  return 0 if $n < 11;
  return 0 if !($n % 2) || !($n % 3) || !($n % 5) || !($n % 7);
  foreach my $i (qw/11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71/) {
    return 2 if $i*$i > $n;   return 0 if !($n % $i);
  }

  if ($n < 105936894253) {   # BPSW seems to be faster after this
    # Deterministic set of Miller-Rabin tests.  If the MR routines can handle
    # bases greater than n, then this can be simplified.
    my @bases;
    if    ($n <          9080191) { @bases = (31,       73); }
    elsif ($n <         19471033) { @bases = ( 2,   299417); }
    elsif ($n <         38010307) { @bases = ( 2,  9332593); }
    elsif ($n <        316349281) { @bases = ( 11000544, 31481107); }
    elsif ($n <       4759123141) { @bases = ( 2, 7, 61); }
    elsif ($n <     105936894253) { @bases = ( 2, 1005905886, 1340600841); }
    elsif ($n <   31858317218647) { @bases = ( 2, 642735, 553174392, 3046413974); }
    elsif ($n < 3071837692357849) { @bases = ( 2, 75088, 642735, 203659041, 3613982119); }
    else                          { @bases = ( 2, 325, 9375, 28178, 450775, 9780504, 1795265022); }
    return Math::Prime::Util::PP::miller_rabin($n, @bases)  ?  2  :  0;
  }

  # BPSW probable prime.  No composites are known to have passed this test
  # since it was published in 1980, though we know infinitely many exist.
  # It has also been verified that no 64-bit composite will return true.
  # Slow since it's all in PP, but it's the Right Thing To Do.

  return 0 unless Math::Prime::Util::PP::miller_rabin($n, 2);
  return 0 unless Math::Prime::Util::PP::is_strong_lucas_pseudoprime($n);
  return ($n <= 18446744073709551615)  ?  2  :  1;
}

sub is_provable_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_positive_integer($n);

  # Shortcut some of the calls.
  return _XS_is_prime($n) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_provable_prime($n) if $_HAVE_GMP
                       && defined &Math::Prime::Util::GMP::is_provable_prime;

  my $is_prob_prime = is_prob_prime($n);
  return $is_prob_prime unless $is_prob_prime == 1;

  # At this point we know it is almost certainly a prime, but we need to
  # prove it.  We should do ECPP or APR-CL now, or failing that, do the
  # Brillhart-Lehmer-Selfridge test, or Pocklington-Lehmer.  Until those
  # are written here, we'll do a Lucas test, which is super simple but may
  # be very slow.  We have AKS code, but it's insanely slow.
  # See http://primes.utm.edu/prove/merged.html or other sources.

  # It shouldn't be possible to get here without BigInt already loaded.
  if (!defined $Math::BigInt::VERSION) {
    eval { require Math::BigInt;   Math::BigInt->import(try=>'GMP,Pari'); 1; }
    or do { croak "Cannot load Math::BigInt"; };
  }
  my $nm1 = $n-1;
  my @factors = factor($nm1);
  # Remember that we have to prove the primality of every factor.
  if ( (scalar grep { is_provable_prime($_) != 2 } @factors) > 0) {
    carp "could not prove primality of $n.\n";
    return 1;
  }

  for (my $a = 2; $a < $nm1; $a++) {
    my $ap = Math::BigInt->new("$a");
    # 1. a^(n-1) = 1 mod n.
    next if $ap->copy->bmodpow($nm1, $n) != 1;
    # 2. a^((n-1)/f) != 1 mod n for all f.
    next if (scalar grep { $_ == 1 }
             map { $ap->copy->bmodpow(int($nm1/$_),$n); }
             @factors) > 0;
    return 2;
  }
  carp "proved $n is not prime\n";
  return 0;
}


#############################################################################

sub prime_count_approx {
  my($x) = @_;
  _validate_positive_integer($x);

  return $_prime_count_small[$x] if $x <= $#_prime_count_small;

  # Below 2^58th or so, all differences between the high precision and C double
  # precision result are less than 0.5.
  if ($x <= $_XS_MAXVAL && $x <= 144115188075855872) {
    return int(_XS_RiemannR($x) + 0.5);
  }

  # Turn on high precision FP if they gave us a big number.
  $x = _upgrade_to_float($x) if ref($x) eq 'Math::BigInt';

  #    Method             10^10 %error  10^19 %error
  #    -----------------  ------------  ------------
  #    average bounds      .01%          .0002%
  #    li(n)               .0007%        .00000004%
  #    li(n)-li(n^.5)/2    .0004%        .00000001%
  #    R(n)                .0004%        .00000001%
  #
  # Also consider: http://trac.sagemath.org/sage_trac/ticket/8135

  # my $result = int( (prime_count_upper($x) + prime_count_lower($x)) / 2);

  # my $result = int( LogarithmicIntegral($x) );

  # my $result = int(LogarithmicIntegral($x) - LogarithmicIntegral(sqrt($x))/2);

  if (ref($x) eq 'Math::BigFloat') {
    # Make sure we get enough accuracy, and also not too much more than needed
    $x->accuracy(length($x->bfloor->bstr())+2);
  }
  my $result = RiemannR($x) + 0.5;

  return Math::BigInt->new($result->bfloor->bstr()) if ref($result) eq 'Math::BigFloat';
  return int($result);
}

sub prime_count_lower {
  my($x) = @_;
  _validate_positive_integer($x);

  return $_prime_count_small[$x] if $x <= $#_prime_count_small;

  $x = _upgrade_to_float($x) if ref($x) eq 'Math::BigInt';

  my $flogx = log($x);

  # Chebyshev:            1*x/logx       x >= 17
  # Rosser & Schoenfeld:  x/(logx-1/2)   x >= 67
  # Dusart 1999:          x/logx*(1+1/logx+1.8/logxlogx)  x >= 32299
  # Dusart 2010:          x/logx*(1+1/logx+2.0/logxlogx)  x >= 88783
  # The Dusart (1999 or 2010) bounds are far, far better than the others.

  my $result;
  if ($x > 1000_000_000_000 && $_Config{'assume_rh'}) {
    my $lix = LogarithmicIntegral($x);
    my $sqx = sqrt($x);
    # Schoenfeld bound:    (constant is 8 * Pi)
    $result = $lix - (($sqx*$flogx) / 25.13274122871834590770114707);
  } elsif ($x < 599) {
    $result = $x / ($flogx - 0.7);   # For smaller numbers this works out well.
  } else {
    my $a;
    # Hand tuned for small numbers (< 60_000M)
    if    ($x <       2700) { $a = 0.30; }
    elsif ($x <       5500) { $a = 0.90; }
    elsif ($x <      19400) { $a = 1.30; }
    elsif ($x <      32299) { $a = 1.60; }
    elsif ($x <     176000) { $a = 1.80; }
    elsif ($x <     315000) { $a = 2.10; }
    elsif ($x <    1100000) { $a = 2.20; }
    elsif ($x <    4500000) { $a = 2.31; }
    elsif ($x <  233000000) { $a = 2.36; }
    elsif ($x < 5433800000) { $a = 2.32; }
    elsif ($x <60000000000) { $a = 2.15; }
    else                    { $a = 2.00; } # Dusart 2010, page 2
    $result = ($x/$flogx) * (1.0 + 1.0/$flogx + $a/($flogx*$flogx));
  }

  return Math::BigInt->new($result->bfloor->bstr()) if ref($result) eq 'Math::BigFloat';
  return int($result);
}

sub prime_count_upper {
  my($x) = @_;
  _validate_positive_integer($x);

  return $_prime_count_small[$x] if $x <= $#_prime_count_small;

  $x = _upgrade_to_float($x) if ref($x) eq 'Math::BigInt';

  # Chebyshev:            1.25506*x/logx       x >= 17
  # Rosser & Schoenfeld:  x/(logx-3/2)         x >= 67
  # Dusart 1999:          x/logx*(1+1/logx+2.51/logxlogx)   x >= 355991
  # Dusart 2010:          x/logx*(1+1/logx+2.334/logxlogx)  x >= 2_953_652_287

  # As with the lower bounds, Dusart bounds are best by far.

  # Another possibility here for numbers under 3000M is to use Li(x)
  # minus a correction.

  my $flogx = log($x);

  my $result;
  if ($x > 10000_000_000_000 && $_Config{'assume_rh'}) {
    my $lix = LogarithmicIntegral($x);
    my $sqx = sqrt($x);
    # Schoenfeld bound:    (constant is 8 * Pi)
    $result = $lix + (($sqx*$flogx) / 25.13274122871834590770114707);
  } elsif ($x <  1621) { $result = ($x / ($flogx - 1.048)) + 1.0; }
    elsif ($x <  5000) { $result = ($x / ($flogx - 1.071)) + 1.0; }
    elsif ($x < 15900) { $result = ($x / ($flogx - 1.098)) + 1.0; }
    else {
    my $a;
    # Hand tuned for small numbers (< 60_000M)
    if    ($x <      24000) { $a = 2.30; }
    elsif ($x <      59000) { $a = 2.48; }
    elsif ($x <     350000) { $a = 2.52; }
    elsif ($x <     355991) { $a = 2.54; }
    elsif ($x <     356000) { $a = 2.51; }
    elsif ($x <    3550000) { $a = 2.50; }
    elsif ($x <    3560000) { $a = 2.49; }
    elsif ($x <    5000000) { $a = 2.48; }
    elsif ($x <    8000000) { $a = 2.47; }
    elsif ($x <   13000000) { $a = 2.46; }
    elsif ($x <   18000000) { $a = 2.45; }
    elsif ($x <   31000000) { $a = 2.44; }
    elsif ($x <   41000000) { $a = 2.43; }
    elsif ($x <   48000000) { $a = 2.42; }
    elsif ($x <  119000000) { $a = 2.41; }
    elsif ($x <  182000000) { $a = 2.40; }
    elsif ($x <  192000000) { $a = 2.395; }
    elsif ($x <  213000000) { $a = 2.390; }
    elsif ($x <  271000000) { $a = 2.385; }
    elsif ($x <  322000000) { $a = 2.380; }
    elsif ($x <  400000000) { $a = 2.375; }
    elsif ($x <  510000000) { $a = 2.370; }
    elsif ($x <  682000000) { $a = 2.367; }
    elsif ($x < 2953652287) { $a = 2.362; }
    else                    { $a = 2.334; } # Dusart 2010, page 2
    #elsif ($x <60000000000) { $a = 2.362; }
    #else                    { $a = 2.51;  } # Dusart 1999, page 14

    # Old versions of Math::BigFloat will do the Wrong Thing with this.
    #return int( ($x/$flogx) * (1.0 + 1.0/$flogx + $a/($flogx*$flogx)) + 1.0 );
    $result = ($x/$flogx) * (1.0 + 1.0/$flogx + $a/($flogx*$flogx)) + 1.0;
  }

  return Math::BigInt->new($result->bfloor->bstr()) if ref($result) eq 'Math::BigFloat';
  return int($result);
}

#############################################################################

sub nth_prime_approx {
  my($n) = @_;
  _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if ref($n) eq 'Math::BigInt';

  my $flogn  = log($n);
  my $flog2n = log($flogn);

  # Cipolla 1902:
  #    m=0   fn * ( flogn + flog2n - 1 );
  #    m=1   + ((flog2n - 2)/flogn) );
  #    m=2   - (((flog2n*flog2n) - 6*flog2n + 11) / (2*flogn*flogn))
  #    + O((flog2n/flogn)^3)
  #
  # Shown in Dusart 1999 page 12, as well as other sources such as:
  #   http://www.emis.de/journals/JIPAM/images/153_02_JIPAM/153_02.pdf
  # where the main issue you run into is that you're doing polynomial
  # interpolation, so it oscillates like crazy with many high-order terms.
  # Hence I'm leaving it at m=2.
  #

  my $approx = $n * ( $flogn + $flog2n - 1
                      + (($flog2n - 2)/$flogn)
                      - ((($flog2n*$flog2n) - 6*$flog2n + 11) / (2*$flogn*$flogn))
                    );

  # Apply a correction to help keep values close.
  my $order = $flog2n/$flogn;
  $order = $order*$order*$order * $n;

  if    ($n <        259) { $approx += 10.4 * $order; }
  elsif ($n <        775) { $approx +=  7.52* $order; }
  elsif ($n <       1271) { $approx +=  5.6 * $order; }
  elsif ($n <       2000) { $approx +=  5.2 * $order; }
  elsif ($n <       4000) { $approx +=  4.3 * $order; }
  elsif ($n <      12000) { $approx +=  3.0 * $order; }
  elsif ($n <     150000) { $approx +=  2.1 * $order; }
  elsif ($n <  200000000) { $approx +=  0.0 * $order; }
  else                    { $approx += -0.010 * $order; }
  # $approx = -0.025 is better for the last, but it gives problems with some
  # other code that always wants the asymptotic approximation to be >= actual.

  if ( ($approx >= ~0) && (ref($approx) ne 'Math::BigFloat') ) {
    return $_Config{'maxprime'} if $n <= $_Config{'maxprimeidx'};
    croak "nth_prime_approx($n) overflow";
  }

  return int($approx + 0.5);
}

# The nth prime will be greater than or equal to this number
sub nth_prime_lower {
  my($n) = @_;
  _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if ref($n) eq 'Math::BigInt';

  my $flogn  = log($n);
  my $flog2n = log($flogn);  # Note distinction between log_2(n) and log^2(n)

  # Dusart 1999 page 14, for all n >= 2
  #my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.25)/$flogn));
  # Dusart 2010 page 2, for all n >= 3
  my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.10)/$flogn));

  if ( ($lower >= ~0) && (ref($lower) ne 'Math::BigFloat') ) {
    return $_Config{'maxprime'} if $n <= $_Config{'maxprimeidx'};
    croak "nth_prime_lower($n) overflow";
  }

  return int($lower);
}

# The nth prime will be less or equal to this number
sub nth_prime_upper {
  my($n) = @_;
  _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if ref($n) eq 'Math::BigInt';

  my $flogn  = log($n);
  my $flog2n = log($flogn);  # Note distinction between log_2(n) and log^2(n)

  my $upper;
  if      ($n >= 688383) {   # Dusart 2010 page 2
    $upper = $n * ( $flogn  +  $flog2n - 1.0 + (($flog2n-2.00)/$flogn) );
  } elsif ($n >= 178974) {   # Dusart 2010 page 7
    $upper = $n * ( $flogn  +  $flog2n - 1.0 + (($flog2n-1.95)/$flogn) );
  } elsif ($n >=  39017) {   # Dusart 1999 page 14
    $upper = $n * ( $flogn  +  $flog2n - 0.9484 );
  } elsif ($n >=      6) {   # Modified Robin 1983, for 6-39016 only
    $upper = $n * ( $flogn  +  0.6000 * $flog2n );
  } else {
    $upper = $n * ( $flogn  +  $flog2n );
  }

  if ( ($upper >= ~0) && (ref($upper) ne 'Math::BigFloat') ) {
    return $_Config{'maxprime'} if $n <= $_Config{'maxprimeidx'};
    croak "nth_prime_upper($n) overflow";
  }

  return int($upper + 1.0);
}


#############################################################################


#############################################################################

sub RiemannZeta {
  my($n) = @_;
  croak("Invalid input to ReimannZeta:  x must be > 0") if $n <= 0;

  return _XS_RiemannZeta($n)
    if !defined $bignum::VERSION && ref($n) ne 'Math::BigFloat' && $n <= $_XS_MAXVAL;
  return Math::Prime::Util::PP::RiemannZeta($n);
}

sub RiemannR {
  my($n) = @_;
  croak("Invalid input to ReimannR:  x must be > 0") if $n <= 0;

  return _XS_RiemannR($n)
    if !defined $bignum::VERSION && ref($n) ne 'Math::BigFloat' && $n <= $_XS_MAXVAL;
  return Math::Prime::Util::PP::RiemannR($n);
}

sub ExponentialIntegral {
  my($n) = @_;
  return $_Neg_Infinity if $n == 0;
  return 0              if $n == $_Neg_Infinity;
  return $_Infinity     if $n == $_Infinity;

  return _XS_ExponentialIntegral($n)
   if !defined $bignum::VERSION && ref($n) ne 'Math::BigFloat' && $_Config{'xs'};

  return Math::Prime::Util::PP::ExponentialIntegral($n);
}

sub LogarithmicIntegral {
  my($n) = @_;
  return 0              if $n == 0;
  return $_Neg_Infinity if $n == 1;
  return $_Infinity     if $n == $_Infinity;

  croak("Invalid input to LogarithmicIntegral:  x must be >= 0") if $n <= 0;

  if (!defined $bignum::VERSION && ref($n) ne 'Math::BigFloat' && $_Config{'xs'}) {
    return 1.045163780117492784844588889194613136522615578151 if $n == 2;
    return _XS_LogarithmicIntegral($n);
  }

  return Math::Prime::Util::PP::LogarithmicIntegral($n);
}

#############################################################################

use Math::Prime::Util::MemFree;

1;

__END__


# ABSTRACT: Utilities related to prime numbers, including fast generators / sievers

=pod

=encoding utf8


=head1 NAME

Math::Prime::Util - Utilities related to prime numbers, including fast sieves and factoring


=head1 VERSION

Version 0.21


=head1 SYNOPSIS

  # Normally you would just import the functions you are using.
  # Nothing is exported by default.  List the functions, or use :all.
  use Math::Prime::Util ':all';


  # Get a big array reference of many primes
  my $aref = primes( 100_000_000 );

  # All the primes between 5k and 10k inclusive
  my $aref = primes( 5_000, 10_000 );

  # If you want them in an array instead
  my @primes = @{primes( 500 )};


  # For non-bigints, is_prime and is_prob_prime will always be 0 or 2.
  # They return return 0 (composite), 2 (prime), or 1 (probably prime)
  say "$n is prime"  if is_prime($n);
  say "$n is ", (qw(composite maybe_prime? prime))[is_prob_prime($n)];

  # Strong pseudoprime test with multiple bases, using Miller-Rabin
  say "$n is a prime or 2/7/61-psp" if is_strong_pseudoprime($n, 2, 7, 61);

  # Strong Lucas-Selfridge test
  say "$n is a prime or slpsp" if is_strong_lucas_pseudoprime($n);

  # step to the next prime (returns 0 if not using bigints and we'd overflow)
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

  # Get all factors
  @divisors = all_factors( $n );

  # Euler phi (Euler's totient) on a large number
  use bigint;  say euler_phi( 801294088771394680000412 );
  say jordan_totient(5, 1234);  # Jordan's totient

  # Moebius function used to calculate Mertens
  $sum += moebius($_) for (1..200); say "Mertens(200) = $sum";
  # Mertens function directly (more efficient for large values)
  say mertens(10_000_000);

  # Exponential of Mangoldt function
  say "lamba(49) = ", log(exp_mangoldt(49));

  # divisor sum
  $sigma  = divisor_sum( $n );
  $sigma2 = divisor_sum( $n, sub { $_[0]*$_[0] } );

  # The primorial n# (product of all primes <= n)
  say "15# (2*3*5*7*11*13) is ", primorial(15);
  # The primorial p(n)# (product of first n primes)
  say "P(9)# (2*3*5*7*11*13*17*19*23) is ", pn_primorial(9);

  # Ei, li, and Riemann R functions
  my $ei = ExponentialIntegral($x);    # $x a real: $x != 0
  my $li = LogarithmicIntegral($x);    # $x a real: $x >= 0
  my $R  = RiemannR($x)                # $x a real: $x > 0


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
  my $rand_prime = random_nbit_prime(128);   # random 128-bit prime
  my $rand_prime = random_strong_prime(256); # random 256-bit strong prime
  my $rand_prime = random_maurer_prime(256); # random 256-bit provable prime


=head1 DESCRIPTION

A set of utilities related to prime numbers.  These include multiple sieving
methods, is_prime, prime_count, nth_prime, approximations and bounds for
the prime_count and nth prime, next_prime and prev_prime, factoring utilities,
and more.

The default sieving and factoring are intended to be (and currently are)
the fastest on CPAN, including L<Math::Prime::XS>, L<Math::Prime::FastSieve>,
L<Math::Factor::XS>, L<Math::Prime::TiedArray>, L<Math::Big::Factors>,
L<Math::Factoring>, and L<Math::Primality> (when the GMP module is available).
For numbers in the 10-20 digit range, it is often orders of magnitude faster.
Typically it is faster than L<Math::Pari> for 64-bit operations.

All operations support both Perl UV's (32-bit or 64-bit) and bignums.  It
requires no external software for big number support, as there are Perl
implementations included that solely use Math::BigInt and Math::BigFloat.
However, performance will be improved for most big number functions by
installing L<Math::Prime::Util::GMP>, and is definitely recommended if you
do many bignum operations.  Also look into L<Math::Pari> as an alternative.

The module is thread-safe and allows concurrency between Perl threads while
still sharing a prime cache.  It is not itself multithreaded.  See the
L<Limitations|/"LIMITATIONS"> section if you are using Win32 and threads in
your program.


=head1 BIGNUM SUPPORT

By default all functions support bignums.  With a few exceptions, the module
will not turn on bignum support for you -- you will need to C<use bigint>,
C<use bignum>, or pass in a L<Math::BigInt> or L<Math::BigFloat> object as
your input.  The functions take some care to perform all bignum operations
using the same class as was passed in, allowing the module to work properly
with Calc, FastCalc, GMP, Pari, etc.  You should try to install
L<Math::Prime::Util::GMP> if you plan to use bigints with this module, as
it will make it run much faster.


Some of the functions, notably:

  factor
  is_prime
  is_prob_prime
  is_strong_pseudoprime
  next_prime
  prev_prime
  nth_prime

work very fast (under 1 microsecond) on small inputs, but the wrappers for
input validation and bigint support take more time than the function itself.
Using the flag '-bigint', e.g.:

  use Math::Prime::Util qw(-bigint);

will turn off bigint support for those functions.  Those functions will then
go directly to the XS versions, which will speed up very small inputs a B<lot>.
This is useful if you're using the functions in a loop, but since the difference
is less than a millisecond, it's really not important in general (also, a
future implementation may find a way to speed this up without the option).


If you are using bigints, here are some performance suggestions:

=over 4

=item *

Install L<Math::Prime::Util::GMP>, as that will vastly increase the speed
of many of the functions.  This does require the L<GMP|gttp://gmplib.org>
library be installed on your system, but this increasingly comes
pre-installed or easily available using the OS vendor package installation tool.

=item *

Install and use L<Math::BigInt::GMP> or L<Math::BigInt::Pari>, then use
C<use bigint try =E<gt> 'GMP,Pari'> in your script, or on the command line
C<-Mbigint=lib,GMP>.  Large modular exponentiation is much faster using the
GMP or Pari backends, as are the math and approximation functions when
called with very large inputs.

=item *

Install L<Math::MPFR> if you use the Ei, li, Zeta, or R functions.  If that
module can be loaded, these functions will run much faster on bignum inputs,
and are able to provide higher accuracy.

=item *

Having run these functions on many versions of Perl, if you're using anything
older than Perl 5.14, I would recommend you upgrade if you are using bignums
a lot.  There are some brittle behaviors on 5.12.4 and earlier with bignums.

=back


=head1 FUNCTIONS

=head2 is_prime

  print "$n is prime" if is_prime($n);

Returns 2 if the number is prime, 0 if not.  For numbers larger than C<2^64>
it will return 0 for composite and 1 for probably prime, using a strong BPSW
test.  If L<Math::Prime::Util::GMP> is installed, some quick primality proofs
are run on larger numbers, so will return 2 for many of those also.

Also see the L</"is_prob_prime"> function, which will never do additional
tests, and the L</"is_provable_prime"> function which will try very hard to
return only 0 or 2 for any input.

For native precision numbers (anything smaller than C<2^64>, all three
functions are identical and use a deterministic set of Miller-Rabin tests.
While L</"is_prob_prime"> and L</"is_prime"> return probable prime results
for larger numbers, they use the strong Baillie-PSW test, which has had
no counterexample found since it was published in 1980 (though certainly they
exist).


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

Returns the next prime greater than the input number.  If the input is not a
bigint, then 0 is returned if the next prime is larger than a native integer
type (the last representable primes being C<4,294,967,291> in 32-bit Perl and
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

The current implementation decides based on the ranges whether to use a
segmented sieve with a fast bit count, or Lehmer's algorithm.  The former
is preferred for small sizes as well as small ranges.  The latter is much
faster for large ranges.

The segmented sieve is very memory efficient and is quite fast even with
large base values.  Its complexity is approximately C<O(sqrt(a) + (b-a))>,
where the first term is typically negligible below C<~ 10^11>.  Memory use
is proportional only to C<sqrt(a)>, with total memory use under 1MB for any
base under C<10^14>.

Lehmer's method has complexity approximately C<O(b^0.7) + O(a^0.7)>.  It
does use more memory however.  A calculation of C<Pi(10^14)> completes in
under 1 minute, C<Pi(10^15)> in under 5 minutes, and C<Pi(10^16)> in under
30 minutes, however using nearly 1400MB of peak memory for the last.
In contrast, even primesieve using 12 cores would take over a week on this
same computer to determine C<Pi(10^16)>.

Also see the function L</"prime_count_approx"> which gives a very good
approximation to the prime count, and L</"prime_count_lower"> and
L</"prime_count_upper"> which give tight bounds to the actual prime count.
These functions return quickly for any input, including bigints.


=head2 prime_count_upper

=head2 prime_count_lower

  my $lower_limit = prime_count_lower($n);
  my $upper_limit = prime_count_upper($n);
  #   $lower_limit  <=  prime_count(n)  <=  $upper_limit

Returns an upper or lower bound on the number of primes below the input number.
These are analytical routines, so will take a fixed amount of time and no
memory.  The actual C<prime_count> will always be equal to or between these
numbers.

A common place these would be used is sizing an array to hold the first C<$n>
primes.  It may be desirable to use a bit more memory than is necessary, to
avoid calling C<prime_count>.

These routines use verified tight limits below a range at least C<2^35>, and
use the Dusart (2010) bounds of

    x/logx * (1 + 1/logx + 2.000/log^2x) <= Pi(x)

    x/logx * (1 + 1/logx + 2.334/log^2x) >= Pi(x)

above that range.  These bounds do not assume the Riemann Hypothesis.  If the
configuration option C<assume_rh> has been set (it is off by default), then
the Schoenfeld (1976) bounds are used for large values.


=head2 prime_count_approx

  print "there are about ",
        prime_count_approx( 10 ** 18 ),
        " primes below one quintillion.\n";

Returns an approximation to the C<prime_count> function, without having to
generate any primes.  The current implementation uses the Riemann R function
which is quite accurate: an error of less than C<0.0005%> is typical for
input values over C<2^32>.  A slightly faster (0.1ms vs. 1ms) but much less
accurate answer can be obtained by averaging the upper and lower bounds.


=head2 nth_prime

  say "The ten thousandth prime is ", nth_prime(10_000);

Returns the prime that lies in index C<n> in the array of prime numbers.  Put
another way, this returns the smallest C<p> such that C<Pi(p) E<gt>= n>.

For relatively small inputs (below 2 million or so), this does a sieve over
a range containing the nth prime, then counts up to the number.  This is fairly
efficient in time and memory.  For larger values, the Dusart 2010 bounds are
calculated, Lehmer's fast prime counting method is used to calculate the
count up to that point, then sieving is done in the range between the bounds.

While this method is hundreds of times faster than generating primes, and
doesn't involve big tables of precomputed values, it still can take a fair
amount of time and space for large inputs.  Calculating the C<10^11th> prime
takes a bit over 2 seconds, the C<10^12th> prime takes 20 seconds, and the
C<10^13th> prime (323780508946331) takes 4 minutes.  Think about whether
a bound or approximation would be acceptable, as they can be computed
analytically.

If the bigint or bignum module is not in use, this will generate an overflow
exception if the number requested would result in a prime that cannot fit in
a native type.  If bigints are in use, then the calculation will proceed,
though it will be exceedingly slow.  A later version of
L<Math::Prime::Util::GMP> may include this functionality which would help for
32-bit machines.


=head2 nth_prime_upper

=head2 nth_prime_lower

  my $lower_limit = nth_prime_lower($n);
  my $upper_limit = nth_prime_upper($n);
  #   $lower_limit  <=  nth_prime(n)  <=  $upper_limit

Returns an analytical upper or lower bound on the Nth prime.  These are very
fast as they do not need to sieve or search through primes or tables.  An
exact answer is returned for tiny values of C<n>.  The lower limit uses the
Dusart 2010 bound for all C<n>, while the upper bound uses one of the two
Dusart 2010 bounds for C<n E<gt>= 178974>, a Dusart 1999 bound for
C<n E<gt>= 39017>, and a simple bound of C<n * (logn + 0.6 * loglogn)>
for small C<n>.


=head2 nth_prime_approx

  say "The one trillionth prime is ~ ", nth_prime_approx(10**12);

Returns an approximation to the C<nth_prime> function, without having to
generate any primes.  Uses the Cipolla 1902 approximation with two
polynomials, plus a correction for small values to reduce the error.


=head2 is_strong_pseudoprime

  my $maybe_prime = is_strong_pseudoprime($n, 2);
  my $probably_prime = is_strong_pseudoprime($n, 2, 3, 5, 7, 11, 13, 17);

Takes a positive number as input and one or more bases.  The bases must be
greater than C<1>.  Returns 1 if the input is a prime or a strong
pseudoprime to all of the bases, and 0 if not.

If 0 is returned, then the number really is a composite.  If 1 is returned,
then it is either a prime or a strong pseudoprime to all the given bases.
Given enough distinct bases, the chances become very, very strong that the
number is actually prime.

This is usually used in combination with other tests to make either stronger
tests (e.g. the strong BPSW test) or deterministic results for numbers less
than some verified limit (e.g. it has long been known that no more than three
selected bases are required to give correct primality test results for any
32-bit number).  Given the small chances of passing multiple bases, there
are some math packages that just use multiple MR tests for primality testing.

Even numbers other than 2 will always return 0 (composite).  While the
algorithm does run with even input, most sources define it only on odd input.
Returning composite for all non-2 even input makes the function match most
other implementations including L<Math::Primality>'s C<is_strong_pseudoprime>
function.

=head2 miller_rabin

An alias for C<is_strong_pseudoprime>.  This name is being deprecated.


=head2 is_strong_lucas_pseudoprime

Takes a positive number as input, and returns 1 if the input is a strong
Lucas pseudoprime using the Selfridge method of choosing D, P, and Q (some
sources call this a strong Lucas-Selfridge pseudoprime).  This is one half
of the BPSW primality test (the Miller-Rabin strong pseudoprime test with
base 2 being the other half).


=head2 is_prob_prime

  my $prob_prime = is_prob_prime($n);
  # Returns 0 (composite), 2 (prime), or 1 (probably prime)

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).

For 64-bit input (native or bignum), this uses a tuned set of Miller-Rabin
tests such that the result will be deterministic.  Either 2, 3, 4, 5, or 7
Miller-Rabin tests are performed (no more than 3 for 32-bit input), and the
result will then always be 0 (composite) or 2 (prime).  A later implementation
may change the internals, but the results will be identical.

For inputs larger than C<2^64>, a strong Baillie-PSW primality test is
performed (aka BPSW or BSW).  This is a probabilistic test, so only
0 (composite) and 1 (probably prime) are returned.  There is a possibility that
composites may be returned marked prime, but since the test was published in
1980, not a single BPSW pseudoprime has been found, so it is extremely likely
to be prime.
While we believe (Pomerance 1984) that an infinite number of counterexamples
exist, there is a weak conjecture (Martin) that none exist under 10000 digits.


=head2 is_provable_prime

  say "$n is definitely prime" if is_provable_prime($n) == 2;

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).  This gives it the same return
values as C<is_prime> and C<is_prob_prime>.

The current implementation uses a Lucas test requiring a complete factorization
of C<n-1>, which may not be possible in a reasonable amount of time.  The GMP
version uses the BLS (Brillhart-Lehmer-Selfridge) method, requiring C<n-1> to
be factored to the cube root of C<n>, which is more likely to succeed and will
usually take less time, but can still fail.  Hence you should always test that
the result is C<2> to ensure the prime is proven.


=head2 is_aks_prime

  say "$n is definitely prime" if is_aks_prime($n);

Takes a positive number as input, and returns 1 if the input passes the
Agrawal-Kayal-Saxena (AKS) primality test.  This is a deterministic
unconditional primality test which runs in polynomial time for general input.

This function is only included for completeness and as an example.  While the
implementation is fast compared to the only other Perl implementation available
(in L<Math::Primality>), it is slow compared to others.  However, even
optimized AKS implementations are far slower than ECPP or other modern
primality tests.


=head2 moebius

  say "$n is square free" if moebius($n) != 0;
  $sum += moebius($_) for (1..200); say "Mertens(200) = $sum";

Returns the Möbius function (also called the Moebius, Mobius, or MoebiusMu
function) for a non-negative integer input.  This function is 1 if
C<n = 1>, 0 if C<n> is not square free (i.e. C<n> has a repeated factor),
and C<-1^t> if C<n> is a product of C<t> distinct primes.  This is an
important function in prime number theory.  Like SAGE, we define
C<moebius(0) = 0> for convenience.

If called with two arguments, they define a range C<low> to C<high>, and the
function returns an array with the value of the Möbius function for every n
from low to high inclusive.  Large values of high will result in a lot of
memory use.  The algorithm used is Deléglise and Rivat (1996) algorithm 4.1,
which is a segmented version of Lioen and van de Lune (1994) algorithm 3.2.


=head2 mertens

  say "Mertens(10M) = ", mertens(10_000_000);   # = 1037

Returns the Mertens function of the positive non-zero integer input.  This
function is defined as C<sum(moebius(1..n))>.  This is a much more efficient
solution for larger inputs.  For example, computing Mertens(100M) takes:

   time    approx mem
     0.36s     0.1MB   mertens(100_000_000)
    74.8s   7000MB     List::Util::sum(moebius(1,100_000_000))
   325.7s      0MB     $sum += moebius($_) for 1..100_000_000

The summation of individual terms via factoring is quite expensive in time,
though uses O(1) space.  Computing the summation via a moebius sieve to C<n>
is much faster, though a segmented sieve must be used for large C<n> to
control the memory taken.  Benito and Varona (2008) show a simple C<n/3>
summation which is much faster and uses less memory.  Better yet is a
simple C<n^1/2> version of Deléglise and Rivat (1996), which is what the
current implementation uses.  Deléglise and Rivat's full segmented C<n^1/3>
algorithm is faster.  Kuznetsov (2011) gives an alternate method that he
indicates is even faster.  In theory using one of the advanced prime count
algorithms can lead to a faster solution.


=head2 euler_phi

  say "The Euler totient of $n is ", euler_phi($n);

Returns the Euler totient function (also called Euler's phi or phi function)
for an integer value.  This is an arithmetic function that counts the number
of positive integers less than or equal to C<n> that are relatively prime to
C<n>.  Given the definition used, C<euler_phi> will return 0 for all
C<n E<lt> 1>.  This follows the logic used by SAGE.  Mathematic/WolframAlpha
also returns 0 for input 0, but returns C<euler_phi(-n)> for C<n E<lt> 0>.

If called with two arguments, they define a range C<low> to C<high>, and the
function returns an array with the totient of every n from low to high
inclusive.  Large values of high will result in a lot of memory use.


=head2 jordan_totient

  say "Jordan's totient J_$k($n) is ", jordan_totient($k, $n);

Returns Jordan's totient function for a given integer value.  Jordan's totient
is a generalization of Euler's totient, where
  C<jordan_totient(1,$n) == euler_totient($n)>
This counts the number of k-tuples less than or equal to n that form a coprime
tuple with n.  As with C<euler_phi>, 0 is returned for all C<n E<lt> 1>.
This function can be used to generate some other useful functions, such as
the Dedikind psi function, where C<psi(n) = J(2,n) / J(1,n)>.


=head2 exp_mangoldt

  say "exp(lambda($_)) = ", exp_mangoldt($_) for 1 .. 100;

The Mangoldt function Λ(n) (also known as von Mangoldt's function) is equal
to log p if n is prime or a power of a prime, and 0 otherwise.  We return
the exponential so all results are integers.  Hence the return value
for C<exp_mangoldt> is:
   C<p> if C<n = p^m> for some prime C<p> and integer C<m E<gt>= 1>
   1 otherwise.


=head2 divisor_sum

  say "Sum of divisors of $n:", divisor_sum( $n );

This function takes a positive integer as input and returns the sum of all
the divisors of the input, including 1 and itself.  This is known as the
sigma function (see Hardy and Wright section 16.7, or OEIS A000203).

The more general form takes a code reference as a second parameter, which
is applied to each divisor before the summation.  This allows computation
of numerous functions such as OEIS A000005 [d(n), sigma_0(n), tau(n)]:

  divisor_sum( $n, sub { 1 } );

OEIS A001157 [sigma_2(n)]:

  divisor_sum( $n, sub { $_[0]*$_[0] } )

the general sigma_k (OEIS A00005, A000203, A001157, A001158, etc.):

  divisor_sum( $n, sub { $_[0] ** $k } );

the 5th Jordan totient (OEIS A059378):

  divisor_sum( $n, sub { my $d=shift; $d**5 * moebius($n/$d); } );

though in the last case we have a function L</"jordan_totient"> to compute
it more efficiently.

This function is useful for calculating things like aliquot sums, abundant
numbers, perfect numbers, etc.

The summation is done as a bigint if the input was a bigint object.  You may
need to ensure the result of the subroutine does not overflow a native int.


=head2 primorial

  $prim = primorial(11); #        11# = 2*3*5*7*11 = 2310

Returns the primorial C<n#> of the positive integer input, defined as the
product of the prime numbers less than or equal to C<n>.  This is the
L<OEIS series A034386|http://oeis.org/A034386>: primorial numbers second
definition.

  primorial(0)  == 1
  primorial($n) == pn_primorial( prime_count($n) )

The result will be a L<Math::BigInt> object if it is larger than the native
bit size.

Be careful about which version (C<primorial> or C<pn_primorial>) matches the
definition you want to use.  Not all sources agree on the terminology, though
they should give a clear definition of which of the two versions they mean.
OEIS, Wikipedia, and Mathworld are all consistent, and these functions should
match that terminology.


=head2 pn_primorial

  $prim = pn_primorial(5); #      p_5# = 2*3*5*7*11 = 2310

Returns the primorial number C<p_n#> of the positive integer input, defined as
the product of the first C<n> prime numbers (compare to the factorial, which
is the product of the first C<n> natural numbers).  This is the
L<OEIS series A002110|http://oeis.org/A002110>: primorial numbers first
definition.

  pn_primorial(0)  == 1
  pn_primorial($n) == primorial( nth_prime($n) )

The result will be a L<Math::BigInt> object if it is larger than the native
bit size.


=head2 random_prime

  my $small_prime = random_prime(1000);      # random prime <= limit
  my $rand_prime = random_prime(100, 10000); # random prime within a range

Returns a psuedo-randomly selected prime that will be greater than or equal
to the lower limit and less than or equal to the upper limit.  If no lower
limit is given, 2 is implied.  Returns undef if no primes exist within the
range.

The goal is to return a uniform distribution of the primes in the range,
meaning for each prime in the range, the chances are equally likely that it
will be seen.  This is removes from consideration such algorithms as
C<PRIMEINC>, which although efficient, gives very non-random output.

For small numbers, a random index selection is done, which gives ideal
uniformity and is very efficient with small inputs.  For ranges larger than
this ~16-bit threshold but within the native bit size, a Monte Carlo method
is used (multiple calls to C<rand> may be made if necessary).  This also
gives ideal uniformity and can be very fast for reasonably sized ranges.
For even larger numbers, we partition the range, choose a random partition,
then select a random prime from the partition.  This gives some loss of
uniformity but results in many fewer bits of randomness being consumed as
well as being much faster.

If an C<irand> function has been set via L</"prime_set_config">, it will be
used to construct any ranged random numbers needed.  The function should
return a uniformly random 32-bit integer, which is how the irand functions
exported by L<Math::Random::Secure>, L<Math::Random::MT>,
L<Math::Random::ISAAC> and most other modules behave.

If no C<irand> function was set, then L<Bytes::Random::Secure> is used with
a non-blocking seed.  This will create good quality random numbers, so there
should be little reason to change unless one is generating long-term keys,
where using the blocking random source may be preferred.

Examples of irand functions:

  # Math::Random::Secure.  Uses ISAAC and strong seed methods.
  use Math::Random::Secure;
  prime_set_config(irand => \&Math::Random::Secure::irand);

  # Bytes::Random::Secure.  Also uses ISAAC and strong seed methods.
  use Bytes::Random::Secure qw/random_bytes/;
  prime_set_config(irand => sub { return unpack("L", random_bytes(4)); });

  # Crypt::Random.  Uses Pari and /dev/random.  Very slow.
  use Crypt::Random qw/makerandom/;
  prime_set_config(irand => sub { makerandom(Size=>32, Uniform=>1); });

  # Mersenne Twister.  Very fast, decent RNG, auto seeding.
  use Math::Random::MT::Auto;
  prime_set_config(irand=>sub {Math::Random::MT::Auto::irand() & 0xFFFFFFFF});


=head2 random_ndigit_prime

  say "My 4-digit prime number is: ", random_ndigit_prime(4);

Selects a random n-digit prime, where the input is an integer number of
digits between 1 and the maximum native type (10 for 32-bit, 20 for 64-bit,
10000 if bigint is active).  One of the primes within that range
(e.g. 1000 - 9999 for 4-digits) will be uniformly selected using the
C<irand> function as described above.

If the number of digits is greater than or equal to the maximum native type,
then the result will be returned as a BigInt.  However, if the '-nobigint'
tag was used, then numbers larger than the threshold will be flagged as an
error, and numbers on the threshold will be restricted to native numbers.


=head2 random_nbit_prime

  my $bigprime = random_nbit_prime(512);

Selects a random n-bit prime, where the input is an integer number of bits
between 2 and the maximum representable bits (32, 64, or 100000 for native
32-bit, native 64-bit, and bigint respectively).  A prime with the nth bit
set will be uniformly selected, with randomness supplied via calls to the
C<irand> function as described above.

Since this uses the random_prime function, all uniformity properties of that
function apply to this.  The n-bit range is partitioned into nearly equal
segments less than C<2^32>, a segment is randomly selected, then the trivial
Monte Carlo algorithm is used to select a prime from within the segment.
This gives a reasonably uniform distribution, doesn't use excessive random
source, and can be very fast.

The result will be a BigInt if the number of bits is greater than the native
bit size.  For better performance with very large bit sizes, install
L<Math::BigInt::GMP>.


=head2 random_strong_prime

  my $bigprime = random_strong_prime(512);

Constructs an n-bit strong prime using Gordon's algorithm.  We consider a
strong prime I<p> to be one where

=over

=item * I<p> is large.   This function requires at least 128 bits.

=item * I<p-1> has a large prime factor I<r>.

=item * I<p+1> has a large prime factor I<s>

=item * I<r-1> has a large prime factor I<t>

=back

Using a strong prime in cryptography guards against easy factoring with
algorithms like Pollard's Rho.  Rivest and Silverman (1999) present a case
that using strong primes is unnecessary, and most modern cryptographic systems
agree.  First, the smoothness does not affect more modern factoring methods
such as ECM.  Second, modern factoring methods like GNFS are far faster than
either method so make the point moot.  Third, due to key size growth and
advances in factoring and attacks, for practical purposes, using large random
primes offer security equivalent to using strong primes.


=head2 random_maurer_prime

  my $bigprime = random_maurer_prime(512);

Construct an n-bit provable prime, using the FastPrime algorithm of
Ueli Maurer (1995).  This is the same algorithm used by L<Crypt::Primes>.
Similar to L</"random_nbit_prime">, the result will be a BigInt if the
number of bits is greater than the native bit size.

The differences between this function and that in L<Crypt::Primes> include

=over

=item *

Version 0.50 of Crypt::Primes can return composites.

=item *

Version 0.50 of Crypt::Primes uses the C<PRIMEINC> algorithm for the base
generator, which gives a very non-uniform distribution.  This differs
from Maurer's algorithm which uses the Monte Carlo algorithm (which is what
this module uses).

=item *

No external libraries are needed for this module, while C::P requires
L<Math::Pari>.  See the next item however.

=item *

Crypt::Primes is quite fast for all sizes since it uses Pari for all heavy
lifting.  M::P::U is really fast for native bit sizes.  It is similar speed
to Crypt::Primes if the BigInt package in use is GMP or Pari, e.g.

   use Math::BigInt lib=>'GMP';

but a lot slower without.  Having the L<Math::Prime::Util::GMP> module
installed helps in any case.

=item *

Crypt::Primes has some useful options for cryptography.

=item *

Crypt::Primes is hardcoded to use L<Crypt::Random>, while M::P::U uses
L<Bytes::Random::Secure>, and also allows plugging in a random function.
This is more flexible, faster, has fewer dependencies, and uses a CSPRNG
for security.

=back

Any feedback on this function would be greatly appreciated.





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


=head2 prime_get_config

  my $cached_up_to = prime_get_config->{'precalc_to'};

Returns a reference to a hash of the current settings.  The hash is copy of
the configuration, so changing it has no effect.  The settings include:

  precalc_to      primes up to this number are calculated
  maxbits         the maximum number of bits for native operations
  xs              0 or 1, indicating the XS code is available
  gmp             0 or 1, indicating GMP code is available
  maxparam        the largest value for most functions, without bigint
  maxdigits       the max digits in a number, without bigint
  maxprime        the largest representable prime, without bigint
  maxprimeidx     the index of maxprime, without bigint
  assume_rh       whether to assume the Riemann hypothesis (default 0)

=head2 prime_set_config

  prime_set_config( assume_rh => 1 );

Allows setting of some parameters.  Currently the only parameters are:

  xs              Allows turning off the XS code, forcing the Pure Perl code
                  to be used.  Set to 0 to disable XS, set to 1 to re-enable.
                  You probably will never want to do this.

  gmp             Allows turning off the use of L<Math::Prime::Util::GMP>,
                  which means using Pure Perl code for big numbers.  Set to
                  0 to disable GMP, set to 1 to re-enable.
                  You probably will never want to do this.

  assume_rh       Allows functions to assume the Riemann hypothesis is true
                  if set to 1.  This defaults to 0.  Currently this setting
                  only impacts prime count lower and upper bounds, but could
                  later be applied to other areas such as primality testing.
                  A later version may also have a way to indicate whether
                  no RH, RH, GRH, or ERH is to be assumed.

  irand           Takes a code ref to an irand function returning a uniform
                  number between 0 and 2**32-1.  This will be used for all
                  random number generation in the module.


=head1 FACTORING FUNCTIONS

=head2 factor

  my @factors = factor(3_369_738_766_071_892_021);
  # returns (204518747,16476429743)

Produces the prime factors of a positive number input, in numerical order.
The special cases of C<n = 0> and C<n = 1> will return C<n>, which
guarantees multiplying the factors together will always result in the
input value, though those are the only cases where the returned factors
are not prime.

The current algorithm for non-bigints is a sequence of small trial division,
a few rounds of Pollard's Rho, SQUFOF, Pollard's p-1, Hart's OLF, a long
run of Pollard's Rho, and finally trial division if anything survives.  This
process is repeated for each non-prime factor.  In practice, it is very rare
to require more than the first Rho + SQUFOF to find a factor.

Factoring bigints works with pure Perl, and can be very handy on 32-bit
machines for numbers just over the 32-bit limit, but it can be B<very> slow
for "hard" numbers.  Installing the L<Math::Prime::Util::GMP> module will speed
up bigint factoring a B<lot>, and all future effort on large number factoring
will be in that module.  If you do not have that module for some reason, use
the GMP or Pari version of bigint if possible
(e.g. C<use bigint try =E<gt> 'GMP,Pari'>), which will run 2-3x faster (though
still 100x slower than the real GMP code).


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

=head2 rsqufof_factor

  my @factors = squfof_factor($n);
  my @factors = rsqufof_factor($n);  # racing multiplier version

Produces factors, not necessarily prime, of the positive number input.  An
optional number of rounds can be given as a second parameter.  It is possible
the function will be unable to find a factor, in which case a single element,
the input, is returned.  This function typically runs very fast.

=head2 prho_factor

=head2 pbrent_factor

  my @factors = prho_factor($n);
  my @factors = pbrent_factor($n);

  # Use a very small number of rounds
  my @factors = prho_factor($n, 1000);

Produces factors, not necessarily prime, of the positive number input.  An
optional number of rounds can be given as a second parameter.  These attempt
to find a single factor using Pollard's Rho algorithm, either the original
version or Brent's modified version.  These are more specialized algorithms
usually used for pre-factoring very large inputs, as they are very fast at
finding small factors.


=head2 pminus1_factor

  my @factors = pminus1_factor($n);
  my @factors = pminus1_factor($n, 1_000);          # set B1 smoothness
  my @factors = pminus1_factor($n, 1_000, 50_000);  # set B1 and B2

Produces factors, not necessarily prime, of the positive number input.  This
is Pollard's C<p-1> method, using two stages with default smoothness
settings of 1_000_000 for B1, and C<10 * B1> for B2.  This method can rapidly
find a factor C<p> of C<n> where C<p-1> is smooth (it has no large factors).



=head1 MATHEMATICAL FUNCTIONS

=head2 ExponentialIntegral

  my $Ei = ExponentialIntegral($x);

Given a non-zero floating point input C<x>, this returns the real-valued
exponential integral of C<x>, defined as the integral of C<e^t/t dt>
from C<-infinity> to C<x>.

If the bignum module has been loaded, all inputs will be treated as if they
were Math::BigFloat objects.

For non-BigInt/BigFloat objects, the result should be accurate to at least 14
digits.

For BigInt / BigFloat objects, we first check to see if the Math::MPFR module
is installed.  If so, then it is used, as it will return results much faster
and can be more accurate.  Accuracy when using MPFR will be equal to the
C<accuracy()> value of the input (or the default BigFloat accuracy, which
is 40 by default).

MPFR is used for positive inputs only.  If Math::MPFR is not installed or the
input is negative, then other methods are used:
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
may be defined as C<Li(x) = li(x) - li(2)>.  Crandall and Pomerance use the
term C<li0> for this function, and define C<li(x) = Li0(x) - li0(2)>.  Due to
this terminology confusion, it is important to check which exact definition is
being used.

If the bignum module has been loaded, all inputs will be treated as if they
were Math::BigFloat objects.

For non-BigInt/BigFloat objects, the result should be accurate to at least 14
digits.

For BigInt / BigFloat objects, we first check to see if the Math::MPFR module
is installed.  If so, then it is used, as it will return results much faster
and can be more accurate.  Accuracy when using MPFR will be equal to the
C<accuracy()> value of the input (or the default BigFloat accuracy, which
is 40 by default).

MPFR is used for inputs greater than 1 only.  If Math::MPFR is not installed or
the input is less than 1, results will be calculated as C<Ei(ln x)>.


=head2 RiemannZeta

  my $z = RiemannZeta($s);

Given a floating point input C<s> where C<s E<gt> 0>, returns the floating
point value of ζ(s)-1, where ζ(s) is the Riemann zeta function.  One is
subtracted to ensure maximum precision for large values of C<s>.  The zeta
function is the sum from k=1 to infinity of C<1 / k^s>.  This function only
uses real arguments, so is basically the Euler Zeta function.

If the bignum module has been loaded, all inputs will be treated as if they
were Math::BigFloat objects.

For non-BigInt/BigFloat objects, the result should be accurate to at least 14
digits.  The XS code uses a rational Chebyshev approximation between 0.5 and 5,
and a series for larger values.

For BigInt / BigFloat objects, we first check to see if the Math::MPFR module
is installed.  If so, then it is used, as it will return results much faster
and can be more accurate.  Accuracy when using MPFR will be equal to the
C<accuracy()> value of the input (or the default BigFloat accuracy, which
is 40 by default).

If Math::MPFR is not installed, then results are calculated using either
Borwein (1991) algorithm 2, or the basic series.  Full input accuracy is
attempted, but there are defects in Math::BigFloat with high accuracy
computations that make this difficult.  It is also very slow.  I highly
recommend installing Math::MPFR for BigFloat computations.


=head2 RiemannR

  my $r = RiemannR($x);

Given a positive non-zero floating point input, returns the floating
point value of Riemann's R function.  Riemann's R function gives a very close
approximation to the prime counting function.

If the bignum module has been loaded, all inputs will be treated as if they
were Math::BigFloat objects.

For non-BigInt/BigFloat objects, the result should be accurate to at least 14
digits.

For BigInt / BigFloat objects, we first check to see if the Math::MPFR module
is installed.  If so, then it is used, as it will return results much faster
and can be more accurate.  Accuracy when using MPFR will be equal to the
C<accuracy()> value of the input (or the default BigFloat accuracy, which
is 40 by default).  Accuracy without MPFR should be 35 digits.



=head1 EXAMPLES

Print pseudoprimes base 17:

    perl -MMath::Prime::Util=:all -E 'my $n=$base|1; while(1) { print "$n " if is_strong_pseudoprime($n,$base) && !is_prime($n); $n+=2; } BEGIN {$|=1; $base=17}'

Print some primes above 64-bit range:

    perl -MMath::Prime::Util=:all -Mbigint -E 'my $start=100000000000000000000; say join "\n", @{primes($start,$start+1000)}'
    # Similar code using Pari:
    # perl -MMath::Pari=:int,PARI,nextprime -E 'my $start = PARI "100000000000000000000"; my $end = $start+1000; my $p=nextprime($start); while ($p <= $end) { say $p; $p = nextprime($p+1); }'


Project Euler, problem 3 (Largest prime factor):

  use Math::Prime::Util qw/factor/;
  use bigint;  # Only necessary for 32-bit machines.
  say "", (factor(600851475143))[-1]

Project Euler, problem 7 (10001st prime):

  use Math::Prime::Util qw/nth_prime/;
  say nth_prime(10_001);

Project Euler, problem 10 (summation of primes):

  use Math::Prime::Util qw/primes/;
  my $sum = 0;
  $sum += $_ for @{primes(2_000_000)};
  say $sum;

Project Euler, problem 21 (Amicable numbers):

  use Math::Prime::Util qw/divisor_sum/;
  sub dsum { my $n = shift; divisor_sum($n) - $n; }
  my $sum = 0;
  foreach my $a (1..10000) {
    my $b = dsum($a);
    $sum += $a + $b if $b > $a && dsum($b) == $a;
  }
  say $sum;

Project Euler, problem 41 (Pandigital prime), brute force command line:

  perl -MMath::Prime::Util=:all -E 'my @p = grep { /1/&&/2/&&/3/&&/4/&&/5/&&/6/&&/7/} @{primes(1000000,9999999)}; say $p[-1];

Project Euler, problem 47 (Distinct primes factors):

  use Math::Prime::Util qw/pn_primorial factor/;
  use List::MoreUtils qw/distinct/;
  sub nfactors { scalar distinct factor(shift); }
  my $n = pn_primorial(4);  # Start with the first 4-factor number
  $n++ while (nfactors($n) != 4 || nfactors($n+1) != 4 || nfactors($n+2) != 4 || nfactors($n+3) != 4);
  say $n;

Project Euler, problem 69, stupid brute force solution (about 5 seconds):

  use Math::Prime::Util qw/euler_phi/;
  my ($n, $max) = (0,0);
  do {
    my $ndivphi = $_ / euler_phi($_);
    ($n, $max) = ($_, $ndivphi) if $ndivphi > $max;
  } for 1..1000000;
  say "$n  $max";

Here's the right way to do PE problem 69 (under 0.03s):

  use Math::Prime::Util qw/pn_primorial/;
  my $n = 0;
  $n++ while pn_primorial($n+1) < 1000000;
  say pn_primorial($n);'

Project Euler, problem 187, stupid brute force solution, ~3 minutes:

  use Math::Prime::Util qw/factor -nobigint/;
  my $nsemis = 0;
  do { my @f = factor($_); $nsemis++ if scalar @f == 2; }
     for 1 .. int(10**8)-1;
  say $nsemis;

Here's the best way for PE187.  Under 30 milliseconds from the command line:

  use Math::Prime::Util qw/primes prime_count -nobigint/;
  use List::Util qw/sum/;
  my $limit = shift || int(10**8);
  my @primes = @{primes(int(sqrt($limit)))};
  say sum( map { prime_count(int(($limit-1)/$primes[$_-1])) - $_ + 1 }
               1 .. scalar @primes );


=head1 LIMITATIONS

I have not completed testing all the functions near the word size limit
(e.g. C<2^32> for 32-bit machines).  Please report any problems you find.

Perl versions earlier than 5.8.0 have issues with 64-bit that show up in the
factoring tests.  The test suite will try to determine if your Perl is broken.
If you use later versions of Perl, or Perl 5.6.2 32-bit, or Perl 5.6.2 64-bit
and keep numbers below C<~ 2^52>, then everything works.  The best solution is
to update to a more recent Perl.

The module is thread-safe and should allow good concurrency on all platforms
that support Perl threads except Win32 (Cygwin works).  With Win32, either
don't use threads or make sure C<prime_precalc> is called before using
C<primes>, C<prime_count>, or C<nth_prime> with large inputs.  This is B<only>
an issue if you use non-Cygwin Win32 and call these routines from within
Perl threads.



=head1 PERFORMANCE

Counting the primes to C<10^10> (10 billion), with time in seconds.
Pi(10^10) = 455,052,511.
The numbers below are for sieving.  Calculating C<Pi(10^10)> takes 0.064
seconds using the Lehmer algorithm in version 0.12.

   External C programs in C / C++:

       1.9  primesieve 3.6 forced to use only a single thread
       2.2  yafu 1.31
       3.8  primegen (optimized Sieve of Atkin, conf-word 8192)
       5.6  Tomás Oliveira e Silva's unoptimized segmented sieve v2 (Sep 2010)
       6.7  Achim Flammenkamp's prime_sieve (32k segments)
       9.3  http://tverniquet.com/prime/ (mod 2310, single thread)
      11.2  Tomás Oliveira e Silva's unoptimized segmented sieve v1 (May 2003)
      17.0  Pari 2.3.5 (primepi)

   Small portable functions suitable for plugging into XS:

       4.1  My segmented SoE used in this module (with unrolled inner loop)
      15.6  My Sieve of Eratosthenes using a mod-30 wheel
      17.2  A slightly modified verion of Terje Mathisen's mod-30 sieve
      35.5  Basic Sieve of Eratosthenes on odd numbers
      33.4  Sieve of Atkin, from Praxis (not correct)
      72.8  Sieve of Atkin, 10-minute fixup of basic algorithm
      91.6  Sieve of Atkin, Wikipedia-like

Perl modules, counting the primes to C<800_000_000> (800 million):

  Time (s)   Module                      Version  Notes
  ---------  --------------------------  -------  -----------
       0.03  Math::Prime::Util           0.12     using Lehmer's method
       0.28  Math::Prime::Util           0.17     segmented mod-30 sieve
       0.47  Math::Prime::Util::PP       0.14     Perl (Lehmer's method)
       0.9   Math::Prime::Util           0.01     mod-30 sieve
       2.9   Math::Prime::FastSieve      0.12     decent odd-number sieve
      11.7   Math::Prime::XS             0.29     "" but needs a count API
      15.0   Bit::Vector                 7.2
      48.9   Math::Prime::Util::PP       0.14     Perl (fastest I know of)
     170.0   Faster Perl sieve (net)     2012-01  array of odds
     548.1   RosettaCode sieve (net)     2012-06  simplistic Perl
    3048.1   Math::Primality             0.08     Perl + Math::GMPz
  ~11000     Math::Primality             0.04     Perl + Math::GMPz
  >20000     Math::Big                   1.12     Perl, > 26GB RAM used

Python can do this in 2.8s using an RWH numpy function, 14.3s using an RWH
pure Python function.  However the standard modules are far slower.  mpmath
v0.17 primepi takes 169.5s and 25+ GB of RAM.  sympi 0.7.1 primepi takes
292.2s.


C<is_prime>: my impressions for various sized inputs:

   Module                   1-10 digits  10-20 digits  BigInts
   -----------------------  -----------  ------------  --------------
   Math::Prime::Util        Very fast    Pretty fast   Slow to Fast (3)
   Math::Prime::XS          Very fast    Very slow (1) --
   Math::Prime::FastSieve   Very fast    N/A (2)       --
   Math::Primality          Very slow    Very slow     Fast
   Math::Pari               Slow         OK            Fast

   (1) trial division only.  Very fast if every factor is tiny.
   (2) Too much memory to hold the sieve (11dig = 6GB, 12dig = ~50GB)
   (3) If L<Math::Prime::Util::GMP> is installed, then all three of the
       BigInt capable modules run at reasonably similar speeds, capable of
       performing the BPSW test on a 3000 digit input in ~ 1 second.  Without
       that module all computations are done in Perl, so this module using
       GMP bigints runs 2-3x slower, using Pari bigints about 10x slower,
       and using the default bigints (Calc) it can run much slower.

The differences are in the implementations:

=over 4

=item L<Math::Prime::Util>

looks in the sieve for a fast bit lookup if that exists (default up to 30,000
but it can be expanded, e.g.  C<prime_precalc>), uses trial division for
numbers higher than this but not too large (0.1M on 64-bit machines,
100M on 32-bit machines), a deterministic set of Miller-Rabin tests for
64-bit and smaller numbers, and a BPSW test for bigints.

=item L<Math::Prime::XS>

does trial divisions, which is wonderful if the input has a small factor
(or is small itself).  But if given a large prime it can take orders of
magnitude longer.  It does not support bigints.

=item L<Math::Prime::FastSieve>

only works in a sieved range, which is really fast if you can do it
(M::P::U will do the same if you call C<prime_precalc>).  Larger inputs
just need too much time and memory for the sieve.

=item L<Math::Primality>

uses GMP for all work.  Under ~32-bits it uses 2 or 3 MR tests, while
above 4759123141 it performs a BPSW test.  This is is fantastic for
bigints over 2^64, but it is significantly slower than native precision
tests.  With 64-bit numbers it is generally an order of magnitude or more
slower than any of the others.  Once bigints are being used, its
performance is quite good.  It is faster than this module unless
L<Math::Prime::Util::GMP> has been installed, in which case this module
is just a little bit faster.

=item L<Math::Pari>

has some very effective code, but it has some overhead to get to it from
Perl.  That means for small numbers it is relatively slow: an order of
magnitude slower than M::P::XS and M::P::Util (though arguably this is
only important for benchmarking since "slow" is ~2 microseconds).  Large
numbers transition over to smarter tests so don't slow down much.  The
C<ispseudoprime(n,0)> function will perform the BPSW test and is fast
even for large inputs.

=back


Factoring performance depends on the input, and the algorithm choices used
are still being tuned.  L<Math::Factor::XS> is very fast when given input with
only small factors, but it slows down rapidly as the smallest factor increases
in size.  For numbers larger than 32 bits, L<Math::Prime::Util> can be 100x or
more faster (a number with only very small factors will be nearly identical,
while a semiprime with large factors will be the extreme end).  L<Math::Pari>'s
underlying algorithms and code are much more mature than this module, and
for 21+ digit numbers will be a better choice.  Small numbers factor much
faster with Math::Prime::Util.  For 30+ digit numbers, L<Math::Pari> is much
faster.  Without the L<Math::Prime::Util::GMP> module, almost all actions on
numbers greater than native scalars will be much faster in Pari.

The presentation here:
 L<http://math.boisestate.edu/~liljanab/BOISECRYPTFall09/Jacobsen.pdf>
has a lot of data on 64-bit and GMP factoring performance I collected in 2009.
Assuming you do not know anything about the inputs, trial division and
optimized Fermat or Lehman work very well for small numbers (<= 10 digits),
while native SQUFOF is typically the method of choice for 11-18 digits (I've
seen claims that a lightweight QS can be faster for 15+ digits).  Some form
of Quadratic Sieve is usually used for inputs in the 19-100 digit range, and
beyond that is the General Number Field Sieve.  For serious factoring,
I recommend looking at
L<yafu|http://sourceforge.net/projects/yafu/>,
L<msieve|http://sourceforge.net/projects/msieve/>,
L<gmp-ecm|http://ecm.gforge.inria.fr/>,
L<GGNFS|http://sourceforge.net/projects/ggnfs/>,
and L<Pari|http://pari.math.u-bordeaux.fr/>.


The primality proving algorithms leave much to be desired.  If you have
numbers larger than C<2^128>, I recommend Pari's C<isprime(n, 2)> which
will run a fast APRCL test, or
L<GMP-ECPP|http://http://sourceforge.net/projects/gmp-ecpp/>.  Either one
will be much faster than the Lucas or BLS algorithms used in MPU for large
inputs.


=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 ACKNOWLEDGEMENTS

Eratosthenes of Cyrene provided the elegant and simple algorithm for finding
the primes.

Terje Mathisen, A.R. Quesada, and B. Van Pelt all had useful ideas which I
used in my wheel sieve.

Tomás Oliveira e Silva has released the source for a very fast segmented sieve.
The current implementation does not use these ideas.  Future versions might.

The SQUFOF implementation being used is a slight modification to the public
domain racing version written by Ben Buhrow.  Enhancements with ideas from
Ben's later code as well as Jason Papadopoulos's public domain implementations
are planned for a later version.  The old SQUFOF implementation, still included
in the code, is my modifications to Ben Buhrow's modifications to Bob
Silverman's code.


=head1 REFERENCES

=over 4

=item *

Pierre Dusart, "Estimates of Some Functions Over Primes without R.H.", preprint, 2010.  L<http://arxiv.org/abs/1002.0442/>

=item *

Pierre Dusart, "Autour de la fonction qui compte le nombre de nombres premiers", PhD thesis, 1998.  In French, but the mathematics is readable and highly recommended reading if you're interesting in prime number bounds.  L<http://www.unilim.fr/laco/theses/1998/T1998_01.html>

=item *

Gabriel Mincu, "An Asymptotic Expansion", <Journal of Inequalities in Pure and Applied Mathematics>, v4, n2, 2003.  A very readable account of Cipolla's 1902 nth prime approximation.  L<http://www.emis.de/journals/JIPAM/images/153_02_JIPAM/153_02.pdf>

=item *

David M. Smith, "Multiple-Precision Exponential Integral and Related Functions", I<ACM Transactions on Mathematical Software>, v37, n4, 2011.  L<http://myweb.lmu.edu/dmsmith/toms2011.pdf>

=item *

Vincent Pegoraro and Philipp Slusallek, "On the Evaluation of the Complex-Valued Exponential Integral", I<Journal of Graphics, GPU, and Game Tools>, v15, n3, pp 183-198, 2011.  L<http://www.cs.utah.edu/~vpegorar/research/2011_JGT/paper.pdf>

=item *

William H. Press et al., "Numerical Recipes", 3rd edition.

=item *

W. J. Cody and Henry C. Thacher, Jr., "Chebyshev approximations for the exponential integral Ei(x)", I<Mathematics of Computation>, v23, pp 289-303, 1969.  L<http://www.ams.org/journals/mcom/1969-23-106/S0025-5718-1969-0242349-2/>

=item *

W. J. Cody and Henry C. Thacher, Jr., "Rational Chebyshev Approximations for the Exponential Integral E_1(x)", I<Mathematics of Computation>, v22, pp 641-649, 1968.

=item *

W. J. Cody, K. E. Hillstrom, and Henry C. Thacher Jr., "Chebyshev Approximations for the Riemann Zeta Function", L<Mathematics of Computation>, v25, n115, pp 537-547, July 1971.

=item *

Ueli M. Maurer, "Fast Generation of Prime Numbers and Secure Public-Key Cryptographic Parameters", 1995.  L<http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.26.2151>

=item *

Pierre-Alain Fouque and Mehdi Tibouchi, "Close to Uniform Prime Number Generation With Fewer Random Bits", pre-print, 2011.  L<http://eprint.iacr.org/2011/481>

=item *

Douglas A. Stoll and Patrick Demichel , "The impact of ζ(s) complex zeros on π(x) for x E<lt> 10^{10^{13}}", L<Mathematics of Computation>, v80, n276, pp 2381-2394, October 2011.  L<http://www.ams.org/journals/mcom/2011-80-276/S0025-5718-2011-02477-4/home.html>

=item *

L<OEIS: Primorial|http://oeis.org/wiki/Primorial>.

=item *

Walter M. Lioen and Jan van de Lune, "Systematic Computations on Mertens' Conjecture and Dirichlet's Divisor Problem by Vectorized Sieving", in I<From Universal Morphisms to Megabytes>, Centrum voor Wiskunde en Informatica, pp. 421-432, 1994.  Describes a nice way to compute a range of Möbius values.  L<http://walter.lioen.com/papers/LL94.pdf>

=item *

Marc Deléglise and Joöl Rivat, "Computing the summation of the Möbius function", I<Experimental Mathematics>, v5, n4, pp 291-295, 1996.  Enhances the Möbius computation in Lioen/van de Lune, and gives a very efficient way to compute the Mertens function.  L<http://projecteuclid.org/euclid.em/1047565447>

=item *

Manuel Benito and Juan L. Varona, "Recursive formulas related to the summation of the Möbius function", I<The Open Mathematics Journal>, v1, pp 25-34, 2007.  Presents a nice algorithm to compute the Mertens functions with only n/3 Möbius values.  L<http://www.unirioja.es/cu/jvarona/downloads/Benito-Varona-TOMATJ-Mertens.pdf>

=back


=head1 COPYRIGHT

Copyright 2011-2012 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
