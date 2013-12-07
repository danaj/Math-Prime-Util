package Math::Prime::Util;
use strict;
use warnings;
use Carp qw/croak confess carp/;
use Math::BigInt try=>"GMP,Pari";

BEGIN {
  $Math::Prime::Util::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::VERSION = '0.35';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK =
  qw( prime_get_config prime_set_config
      prime_precalc prime_memfree
      is_prime is_prob_prime is_provable_prime is_provable_prime_with_cert
      prime_certificate verify_prime
      is_pseudoprime is_strong_pseudoprime
      is_lucas_pseudoprime
      is_strong_lucas_pseudoprime
      is_extra_strong_lucas_pseudoprime
      is_almost_extra_strong_lucas_pseudoprime
      is_frobenius_underwood_pseudoprime
      is_aks_prime
      miller_rabin miller_rabin_random
      lucas_sequence
      primes
      forprimes prime_iterator prime_iterator_object
      next_prime  prev_prime
      prime_count
      prime_count_lower prime_count_upper prime_count_approx
      nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
      random_prime random_ndigit_prime random_nbit_prime random_strong_prime
      random_proven_prime random_proven_prime_with_cert
      random_maurer_prime random_maurer_prime_with_cert
      primorial pn_primorial consecutive_integer_lcm
      factor factor_exp all_factors
      moebius mertens euler_phi jordan_totient exp_mangoldt liouville
      partitions
      chebyshev_theta chebyshev_psi
      divisor_sum
      carmichael_lambda znorder
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
  undef *factor;          *factor            = \&_XS_factor;
  undef *factor_exp;      *factor_exp        = \&_XS_factor_exp;
 #undef *prime_count;     *prime_count       = \&_XS_prime_count;
  undef *nth_prime;       *nth_prime         = \&_XS_nth_prime;
  undef *is_pseudoprime;  *is_pseudoprime    = \&_XS_is_pseudoprime;
  undef *is_strong_pseudoprime;  *is_strong_pseudoprime = \&_XS_miller_rabin;
  undef *moebius;         *moebius           = \&_XS_moebius;
  undef *mertens;         *mertens           = \&_XS_mertens;
  undef *euler_phi;       *euler_phi         = \&_XS_totient;
  undef *chebyshev_theta; *chebyshev_theta   = \&_XS_chebyshev_theta;
  undef *chebyshev_psi;   *chebyshev_psi     = \&_XS_chebyshev_psi;
  # These should be fast anyway, but this skips validation.
  undef *is_prime;        *is_prime          = \&_XS_is_prime;
  undef *is_prob_prime;   *is_prob_prime     = \&_XS_is_prob_prime;
  undef *next_prime;      *next_prime        = \&_XS_next_prime;
  undef *prev_prime;      *prev_prime        = \&_XS_prev_prime;
  1;
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

    *_validate_num = \&Math::Prime::Util::PP::_validate_num;
    *is_prob_prime = \&Math::Prime::Util::_generic_is_prob_prime;
    *is_prime      = \&Math::Prime::Util::_generic_is_prime;
    *next_prime    = \&Math::Prime::Util::_generic_next_prime;
    *prev_prime    = \&Math::Prime::Util::_generic_prev_prime;
    *exp_mangoldt  = \&Math::Prime::Util::_generic_exp_mangoldt;
    *forprimes     = sub (&$;$) { _generic_forprimes(@_); }; ## no critic qw(ProhibitSubroutinePrototypes)

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
    *pplus1_factor  = \&Math::Prime::Util::PP::pminus1_factor;   # TODO: implement PP p+1.
  };

  $_Config{'nobigint'} = 0;
  $_Config{'gmp'} = 0;
  # See if they have the GMP module and haven't requested it not to be used.
  if (!defined $ENV{MPU_NO_GMP} || $ENV{MPU_NO_GMP} != 1) {
    $_Config{'gmp'} = 1 if eval { require Math::Prime::Util::GMP;
                                  Math::Prime::Util::GMP->import();
                                  1; };
  }

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
_XS_set_callgmp($_HAVE_GMP) if $_Config{'xs'};

# Infinity in Perl is rather O/S specific.
our $_Infinity = 0+'inf';
$_Infinity = 20**20**20 if 65535 > $_Infinity;   # E.g. Windows
our $_Neg_Infinity = -$_Infinity;

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
      _XS_set_callgmp($_HAVE_GMP) if $_Config{'xs'};
    } elsif ($param eq 'nobigint') {
      $_Config{'nobigint'} = ($value) ? 1 : 0;
    } elsif ($param eq 'irand') {
      croak "irand must supply a sub" unless (!defined $value) || (ref($value) eq 'CODE');
      $_Config{'irand'} = $value;
      _clear_randf();  # Force a new randf to be generated
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

sub _validate_positive_integer {
  my($n, $min, $max) = @_;
  # We've gone through _validate_num already, so we just need to handle bigints
  croak "Parameter '$n' must be a positive integer"
     if ref($n) eq 'Math::BigInt' && $n->sign() ne '+';
  croak "Parameter '$n' must be >= $min" if defined $min && $n < $min;
  croak "Parameter '$n' must be <= $max" if defined $max && $n > $max;

  $_[0] = Math::BigInt->new("$_[0]") unless ref($_[0]) eq 'Math::BigInt';
  if ($_[0]->bacmp(''.~0) <= 0) {
    $_[0] = int($_[0]->bstr);
  } else {
    $_[0]->upgrade(undef) if $_[0]->upgrade();  # Stop BigFloat upgrade
  }
  # One of these will be true:
  #     1) $n <= ~0 and $n is not a bigint
  #     2) $n  > ~0 and $n is a bigint
  1;
}

sub _upgrade_to_float {
  do { require Math::BigFloat; Math::BigFloat->import(); }
    if !defined $Math::BigFloat::VERSION;
  return Math::BigFloat->new($_[0]);
}

my @_primes_small = (
   0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509);
sub _tiny_prime_count {
  my($n) = @_;
  return if $n >= $_primes_small[-1];
  my $j = $#_primes_small;
  my $i = 1 + ($n >> 4);
  while ($i < $j) {
    my $mid = ($i+$j)>>1;
    if ($_primes_small[$mid] <= $n) { $i = $mid+1; }
    else                            { $j = $mid;   }
  }
  return $i-1;
}





#############################################################################

sub primes {
  my $optref = (ref $_[0] eq 'HASH')  ?  shift  :  {};
  croak "no parameters to primes" unless scalar @_ > 0;
  croak "too many parameters to primes" unless scalar @_ <= 2;
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;

  _validate_num($low) || _validate_positive_integer($low);
  _validate_num($high) || _validate_positive_integer($high);

  my $sref = [];
  return $sref if ($low > $high) || ($high < 2);

  if ($high > $_XS_MAXVAL) {
    if ($_HAVE_GMP) {
      $sref = Math::Prime::Util::GMP::primes($low,$high);
      if ($high > ~0) {
        # Convert the returned strings into BigInts
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
#  "Fast Generation of Prime Numbers and Secure Public-Key
#   Cryptographic Parameters" by Ueli M. Maurer, 1995
#  http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.26.2151
#  related discussions:
#      http://www.daimi.au.dk/~ivan/provableprimesproject.pdf
#      Handbook of Applied Cryptography by Menezes, et al.
#
#  "Close to Uniform Prime Number Generation With Fewer Random Bits"
#   by Pierre-Alain Fouque and Mehdi Tibouchi, 2011
#   http://eprint.iacr.org/2011/481
#
#  Some things to note:
#
#    1) Joye and Paillier have patents on their methods.  Never use them.
#
#    2) The easy method of next_prime(random number), known as PRIMEINC, is
#       fast but gives a terrible distribution.  It has a positive bias and
#       most importantly the probability for a prime is proportional to its
#       gap, which makes a terrible distribution (some numbers in the range
#       will be thousands of times more likely than others).
#
# We use:
#   TRIVIAL range within native integer size (2^32 or 2^64)
#   FTA1    random_nbit_prime with 65+ bits
#   INVA1   other ranges with 65+ bit range
# where
#   TRIVIAL = monte-carlo method or equivalent, perfect uniformity.
#   FTA1    = Fouque/Tibouchi A1, very close to uniform
#   INVA1   = inverted FTA1, less uniform but works with arbitrary ranges
#
# The random_maurer_prime function uses Maurer's FastPrime algorithm.
#
# If Math::Prime::Util::GMP is installed, these functions will be many times
# faster than other methods (e.g. Math::Pari monte-carlo or Crypt::Primes).
#
# Timings on x86_64, with Math::BigInt::GMP and Math::Random::ISAAC::XS.
#
#                   random_nbit_prime         random_maurer_prime
#    n-bits       no GMP   w/ MPU::GMP        no GMP   w/ MPU::GMP
#    ----------  --------  -----------       --------  -----------
#       24-bit       22uS      same             same       same
#       64-bit       94uS      same             same       same
#      128-bit     0.017s      0.0020s         0.098s      0.056s
#      256-bit     0.033s      0.0033s         0.25s       0.15s
#      512-bit     0.066s      0.0093s         0.65s       0.30s
#     1024-bit     0.16s       0.060s          1.3s        0.94s
#     2048-bit     0.83s       0.5s            3.2s        3.1s
#     4096-bit     6.6s        4.0s           23s         12.0s
#
# Writing these entirely in GMP has a problem, which is that we want to use
# a user-supplied rand function, which means a lot of callbacks.  One
# possibility is to, if they do not supply a rand function, use the GMP MT
# function with an appropriate seed.
#
# Random timings for 10M calls:
#    1.92    system rand
#    2.62    Math::Random::MT::Auto
#   12.0     Math::Random::Secure           w/ISAAC::XS
#   12.6     Bytes::Random::Secure OO       w/ISAAC::XS     <==== our
#   31.1     Bytes::Random::Secure OO                       <==== default
#   44.5     Bytes::Random::Secure function w/ISAAC::XS
#   44.8     Math::Random::Secure
#   71.5     Bytes::Random::Secure function
# 1840       Crypt::Random
#
# time perl -E 'sub irand {int(rand(4294967296));} irand() for 1..10000000;'
# time perl -E 'use Math::Random::MT::Auto qw/irand/; irand() for 1..10000000;'
# time perl -E 'use Math::Random::Secure qw/irand/; irand() for 1..10000000;'
# time perl -E 'use Bytes::Random::Secure qw/random_bytes/; sub irand {return unpack("L",random_bytes(4));} irand() for 1..10000000;'
# time perl -E 'use Bytes::Random::Secure; my $rng = Bytes::Random::Secure->new(); sub irand {return $rng->irand;} irand() for 1..10000000;'
# time perl -E 'use Crypt::Random qw/makerandom/; sub irand {makerandom(Size=>32, Uniform=>1, Strength=>0)} irand() for 1..100_000;'
# > haveged daemon running to stop /dev/random blocking
# > Both BRS and CR have more features that this isn't measuring.
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
    return if $_big_gcd_use >= 0;
    if ($_HAVE_GMP) {
      $_big_gcd_use = 0;
      return;
    }
    if (Math::BigInt->config()->{lib} !~ /^Math::BigInt::(GMP|Pari)/) {
      $_big_gcd_use = 0;
      return;
    }
    $_big_gcd_use = 1;
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
  my $_RANDF;
  my $_RANDF_NBIT;
  sub _set_randf {
    # First define a function $irandf that returns a 32-bit integer.  This
    # corresponds to the irand function of many CPAN modules:
    #    Math::Random::MT
    #    Math::Random::ISAAC
    #    Math::Random::Xorshift
    #    Math::Random::Secure
    # (but not Math::Random::MT::Auto which will return 64-bits)
    my $irandf = $_Config{'irand'};
    if (!defined $irandf) {   # Default irand: BRS nonblocking
      require Bytes::Random::Secure;
      $_BRS = Bytes::Random::Secure->new(NonBlocking=>1) unless defined $_BRS;
      $_RANDF_NBIT = sub {
        my($bits) = @_;
        return 0 if $bits <= 0;
        return ($_BRS->irand() >> (32-$bits))
          if $bits <= 32;
        return ((($_BRS->irand() << 32) + $_BRS->irand()) >> (64-$bits))
          if $bits <= 64 && ~0 > 4294967295;
        my $bytes = int(($bits+7)/8);
        my $n = Math::BigInt->from_hex('0x' . $_BRS->bytes_hex($bytes));
        $n->brsft( 8*$bytes - $bits ) unless ($bits % 8) == 0;
        return $n;
      };
      $_RANDF = sub {
        my($max) = @_;
        my $range = $max+1;
        my $U;
        if (ref($range) eq 'Math::BigInt') {
          my $bits = length($range->as_bin) - 2;   # bits in range
          my $bytes = 1 + int(($bits+7)/8);  # extra byte to reduce ave. loops
          my $rmax = Math::BigInt->bone->blsft($bytes*8)->bdec();
          my $overflow = $rmax - ($rmax % $range);
          do {
            $U = Math::BigInt->from_hex('0x' . $_BRS->bytes_hex($bytes));
          } while $U >= $overflow;
        } elsif ($range <= 4294967295) {
          my $overflow = 4294967295 - (4294967295 % $range);
          do {
            $U = $_BRS->irand();
          } while $U >= $overflow;
        } else {
          croak "randf given max out of bounds: $max" if $range > ~0;
          my $overflow = 18446744073709551615 - (18446744073709551615 % $range);
          do {
            $U = ($_BRS->irand() << 32) + $_BRS->irand();
          } while $U >= $overflow;
        }
        return $U % $range;
      };
    } else { # Custom irand
      $_RANDF_NBIT = sub {
        my($bits) = @_;
        return 0 if $bits <= 0;
        return ($irandf->() >> (32-$bits))
          if $bits <= 32;
        return ((($irandf->() << 32) + $irandf->()) >> (64-$bits))
          if $bits <= 64 && ~0 > 4294967295;
        my $words = int(($bits+31)/32);
        my $n = Math::BigInt->from_hex
          ("0x" . join '', map { sprintf("%08X", $irandf->()) } 1 .. $words );
        $n->brsft( 32*$words - $bits ) unless ($bits % 32) == 0;
        return $n;
      };
      $_RANDF = sub {
        my($max) = @_;
        return 0 if $max <= 0;
        my $range = $max+1;
        my $U;
        if (ref($range) eq 'Math::BigInt') {
          my $zero = $range->copy->bzero;
          my $rbits = length($range->as_bin) - 2;   # bits in range
          my $rwords = int($rbits/32) + (($rbits % 32) ? 1 : 0);
          my $rmax = Math::BigInt->bone->blsft($rwords*32)->bdec();
          my $overflow = $rmax - ($rmax % $range);
          do {
            $U = $range->copy->from_hex
              ("0x" . join '', map { sprintf("%08X", $irandf->()) } 1 .. $rwords);
          } while $U >= $overflow;
        } elsif ($range <= 4294967295) {
          my $overflow = 4294967295 - (4294967295 % $range);
          do {
            $U = $irandf->();
          } while $U >= $overflow;
        } else {
          croak "randf given max out of bounds: $max" if $range > ~0;
          my $overflow = 18446744073709551615 - (18446744073709551615 % $range);
          do {
            $U = ($irandf->() << 32) + $irandf->();
          } while $U >= $overflow;
        }
        return $U % $range;
      };
    }
  }
  sub _clear_randf {
    undef $_RANDF;
    undef $_RANDF_NBIT;
    undef $_BRS;
  }
  sub _get_randf {
    _set_randf() unless defined $_RANDF;
    return $_RANDF;
  }
  sub _get_randf_nbit {
    _set_randf() unless defined $_RANDF_NBIT;
    return $_RANDF_NBIT;
  }

  # Sub to call with low and high already primes and verified range.
  my $_random_prime = sub {
    my($low,$high) = @_;
    my $prime;

    _set_randf() unless defined $_RANDF;

    #{ my $bsize = 100; my @bins; my $counts = 10000000;
    #  for my $c (1..$counts) { $bins[ $irandf->($bsize-1) ]++; }
    #  for my $b (0..$bsize) {printf("%4d %8.5f%%\n", $b, $bins[$b]/$counts);} }

    # low and high are both odds, and low < high.

    # This is fast for small values, low memory, perfectly uniform, and
    # consumes the minimum amount of randomness needed.  But it isn't feasible
    # with large values.  Also note that low must be a prime.
    if ($high <= 262144 && $high <= $_XS_MAXVAL) {
      my $li     = _XS_prime_count(2, $low);
      my $irange = _XS_prime_count($low, $high);
      my $rand = $_RANDF->($irange-1);
      return _XS_nth_prime($li + $rand);
    }

    $low-- if $low == 2;  # Low of 2 becomes 1 for our program.
    confess "Invalid _random_prime parameters: $low, $high" if ($low % 2) == 0 || ($high % 2) == 0;

    # We're going to look at the odd numbers only.
    #my $range = $high - $low + 1;
    my $oddrange = (($high - $low) >> 1) + 1;

    croak "Large random primes not supported on old Perl" if $] < 5.008 && $_Config{'maxbits'} > 32 && $oddrange > 4294967295;

    # If $low is large (e.g. >10 digits) and $range is small (say ~10k), it
    # would be fastest to call primes in the range and randomly pick one.  I'm
    # not implementing it now because it seems like a rare case.

    # If the range is reasonably small, generate using simple Monte Carlo
    # method (aka the 'trivial' method).  Completely uniform.
    if ($oddrange < $_Config{'maxparam'}) {
      my $loop_limit = 2000 * 1000;  # To protect against broken rand
      if ($low > 11) {
        while ($loop_limit-- > 0) {
          $prime = $low + 2 * $_RANDF->($oddrange-1);
          next if !($prime % 3) || !($prime % 5) || !($prime % 7) || !($prime % 11);
          return $prime if is_prob_prime($prime);
        }
      } else {
        while ($loop_limit-- > 0) {
          $prime = $low + 2 * $_RANDF->($oddrange-1);
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
    # we'd like.  Imagine we started at 0 with 1e18 partitions of size 100k
    # each.
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
    # Partitions are typically much larger than 100k, but with a huge range
    # we still see this (e.g. ~3x from 0-10^30, ~10x from 0-10^100).
    #
    # When selecting n-bit or n-digit primes, this effect is MUCH smaller, as
    # the skew becomes approx lg(2^n) / lg(2^(n-1)) which is pretty close to 1.
    #
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

    my $rpart = $_RANDF->($nparts);

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
        my $rand = $_RANDF->($partsize-1);
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
    _make_big_gcds() if $_big_gcd_use < 0;

    while ($loop_limit-- > 0) {
      my $rand = $_RANDF->($partsize-1);
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
      # With GMP, the fastest thing to do is check primality.
      if ($_HAVE_GMP) {
        next unless Math::Prime::Util::GMP::is_prime($prime);
        return $prime;
      }
      # No MPU:GMP, so primality checking is slow.  Skip some composites here.
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
  my %_random_cache_small;

  # For fixed small ranges with XS, e.g. 6-digit, 18-bit
  sub _random_xscount_prime {
    my($low,$high) = @_;
    my($istart, $irange);
    my $cachearef = $_random_cache_small{$low,$high};
    if (defined $cachearef) {
      ($istart, $irange) = @$cachearef;
    } else {
      my $beg = ($low <= 2)  ?  2  :  next_prime($low-1);
      my $end = ($high < ~0)  ?  prev_prime($high + 1)  :  prev_prime($high);
      ($istart, $irange) = ( _XS_prime_count(2, $beg), _XS_prime_count($beg, $end) );
      $_random_cache_small{$low,$high} = [$istart, $irange];
    }
    _set_randf() unless defined $_RANDF;
    my $rand = $_RANDF->($irange-1);
    return _XS_nth_prime($istart + $rand);
  }

  sub random_prime {
    my $low = (@_ == 2)  ?  shift  :  2;
    my $high = shift;
    _validate_num($low) || _validate_positive_integer($low);
    _validate_num($high) || _validate_positive_integer($high);

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
    _validate_num($digits, 1) || _validate_positive_integer($digits, 1);

    return _random_xscount_prime( int(10 ** ($digits-1)), int(10 ** $digits) )
      if $digits <= 6 && int(10**$digits) <= $_XS_MAXVAL;

    my $bigdigits = $digits >= $_Config{'maxdigits'};
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

  my @_random_nbit_m;
  my @_random_nbit_lambda;
  my @_random_nbit_arange;

  sub random_nbit_prime {
    my($bits) = @_;
    _validate_num($bits, 2) || _validate_positive_integer($bits, 2);

    _set_randf() unless defined $_RANDF_NBIT;

    # Very small size, use the nth-prime method
    if ($bits <= 18 && int(2**$bits) <= $_XS_MAXVAL) {
      if ($bits <= 4) {
        return (2,3)[$_RANDF_NBIT->(1)] if $bits == 2;
        return (5,7)[$_RANDF_NBIT->(1)] if $bits == 3;
        return (11,13)[$_RANDF_NBIT->(1)] if $bits == 4;
      }
      return _random_xscount_prime( 1 << ($bits-1), 1 << $bits );
    }

    croak "Mid-size random primes not supported on broken old Perl"
      if $] < 5.008 && $bits > 49
      && $_Config{'maxbits'} > 32 && $bits <= $_Config{'maxbits'};

    # Fouque and Tibouchi (2011) Algorithm 1 (basic)
    # Modified to make sure the nth bit is always set.
    #
    # Example for random_nbit_prime(512) on 64-bit Perl:
    # p:  1aaaaaaaabbbbbbbbbbbbbbbbbbbb1
    #     ^^       ^                   ^--- Trailing 1 so p is odd
    #     ||       +--- 512-63-2 = 447 lower bits selected before loop
    #     |+--- 63 upper bits selected in loop, repeated until p is prime
    #     +--- upper bit is 1 so we generate an n-bit prime
    # total: 1 + 63 + 447 + 1 = 512 bits
    #
    # Algorithm 2 is implemented in a previous commit on github.  The problem
    # is that it doesn't set the nth bit, and making that change requires a
    # modification of the algorithm.  It was not a lot faster than this A1
    # with the native int trial division.  If the irandf function was very
    # slow, then A2 would look more promising.
    #
    if (1 && $bits > 64) {
      my $l = ($_Config{'maxbits'} > 32 && $bits > 79)  ?  63  :  31;
      $l = 49 if $l == 63 && $] < 5.008;  # Fix for broken Perl 5.6
      $l = $bits-2 if $bits-2 < $l;

      my $brand = $_RANDF_NBIT->($bits-$l-2);
      $brand = Math::BigInt->new("$brand") unless ref($brand) eq 'Math::BigInt';
      my $b = $brand->blsft(1)->binc();

      # Precalculate some modulii so we can do trial division on native int
      # 9699690 = 2*3*5*7*11*13*17*19, so later operations can be native ints
      my @premod;
      my $bpremod = int($b->copy->bmod(9699690)->bstr);
      my $twopremod = int(Math::BigInt->new(2)->bmodpow($bits-$l-1, 9699690)->bstr);
      foreach my $zi (0 .. 19-1) {
        foreach my $pm (3, 5, 7, 11, 13, 17, 19) {
          next if $zi >= $pm || defined $premod[$pm];
          $premod[$pm] = $zi if ( ($twopremod*$zi+$bpremod) % $pm) == 0;
        }
      }
      _make_big_gcds() if $_big_gcd_use < 0;
      my $loop_limit = 1_000_000;
      while ($loop_limit-- > 0) {
        my $a = (1 << $l) + $_RANDF_NBIT->($l);
        # $a % s == $premod[s]  =>  $p % s == 0  =>  p will be composite
        next if $a %  3 == $premod[ 3] || $a %  5 == $premod[ 5]
             || $a %  7 == $premod[ 7] || $a % 11 == $premod[11]
             || $a % 13 == $premod[13] || $a % 17 == $premod[17]
             || $a % 19 == $premod[19];
        my $p = Math::BigInt->new("$a")->blsft($bits-$l-1)->badd($b);
        #die " $a $b $p" if $a % 11 == $premod[11] && $p % 11 != 0;
        #die "!$a $b $p" if $a % 11 != $premod[11] && $p % 11 == 0;
        if ($_HAVE_GMP) {
          next unless Math::Prime::Util::GMP::is_prime($p);
        } else {
          next unless Math::BigInt::bgcd($p, 1348781387) == 1; # 23-43
          if ($_big_gcd_use && $p > $_big_gcd_top) {
            next unless Math::BigInt::bgcd($p, $_big_gcd[0]) == 1;
            next unless Math::BigInt::bgcd($p, $_big_gcd[1]) == 1;
            next unless Math::BigInt::bgcd($p, $_big_gcd[2]) == 1;
            next unless Math::BigInt::bgcd($p, $_big_gcd[3]) == 1;
          }
          # We know we don't have GMP and are > 2^64, so skip all the middle.
          #next unless is_prob_prime($p);
          next unless Math::Prime::Util::PP::miller_rabin($p, 2);
          next unless Math::Prime::Util::PP::is_extra_strong_lucas_pseudoprime($p);
        }
        return $p;
      }
      croak "Random function broken?";
    }

    # The Trivial method.  Great uniformity, and fine for small sizes.  It
    # gets very slow as the bit size increases, but that is why we have the
    # method above for bigints.
    if (1) {

      my $loop_limit = 2_000_000;
      if ($bits > $_Config{'maxbits'}) {
        my $p = Math::BigInt->bone->blsft($bits-1)->binc();
        while ($loop_limit-- > 0) {
          my $n = Math::BigInt->new(''.$_RANDF_NBIT->($bits-2))->blsft(1)->badd($p);
          return $n if is_prob_prime($n);
        }
      } else {
        my $p = (1 << ($bits-1)) + 1;
        while ($loop_limit-- > 0) {
          my $n = $p + ($_RANDF_NBIT->($bits-2) << 1);
          return $n if is_prob_prime($n);
        }
      }
      croak "Random function broken?";

    } else {

      # Send through the generic random_prime function.  Decently fast, but
      # quite a bit slower than the F&T A1 method above.
      if (!defined $_random_nbit_ranges[$bits]) {
        if ($bits > $_Config{'maxbits'}) {
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
  }

  sub random_maurer_prime {
    my $k = shift;
    _validate_num($k, 2) || _validate_positive_integer($k, 2);
    if ($k <= $_Config{'maxbits'} && $] >= 5.008) {
      return random_nbit_prime($k);
    }
    my ($n, $cert) = random_maurer_prime_with_cert($k);
    croak "maurer prime $n failed certificate verification!"
          unless verify_prime($cert);
    return $n;
  }

  sub random_maurer_prime_with_cert {
    my($k) = @_;
    _validate_num($k, 2) || _validate_positive_integer($k, 2);

    # Results for random_nbit_prime are proven for all native bit sizes.
    my $p0 = $_Config{'maxbits'};
    $p0 = 49 if $] < 5.008 && $_Config{'maxbits'} > 49;

    if ($k <= $p0) {
      my $n = random_nbit_prime($k);
      my ($isp, $cert) = is_provable_prime_with_cert($n);
      croak "small nbit prime could not be proven" if $isp != 2;
      return ($n, $cert);
    }

    # Set verbose to 3 to get pretty output like Crypt::Primes
    my $verbose = $_Config{'verbose'};
    local $| = 1 if $verbose > 2;

    do { require Math::BigFloat; Math::BigFloat->import(); }
      if !defined $Math::BigFloat::VERSION;

    # Ignore Maurer's g and c that controls how much trial division is done.
    my $r = Math::BigFloat->new("0.5");   # relative size of the prime q
    my $m = 20;                           # makes sure R is big enough
    _set_randf() unless defined $_RANDF;

    # Generate a random prime q of size $r*$k, where $r >= 0.5.  Try to
    # cleverly select r to match the size of a typical random factor.
    if ($k > 2*$m) {
      do {
        my $s = Math::BigFloat->new($_RANDF->(2147483647))->bdiv(2147483648);
        $r = Math::BigFloat->new(2)->bpow($s-1);
      } while ($k*$r >= $k-$m);
    }

    # I've seen +0, +1, and +2 here.  Maurer uses +0.  Menezes uses +1.
    my ($q, $qcert) = random_maurer_prime_with_cert( ($r * $k)->bfloor + 1 );
    $q = Math::BigInt->new("$q") unless ref($q) eq 'Math::BigInt';
    my $I = Math::BigInt->new(2)->bpow($k-2)->bdiv($q)->bfloor->as_int();
    print "r = $r  k = $k  q = $q  I = $I\n" if $verbose && $verbose != 3;
    $qcert = ($q < Math::BigInt->new("18446744073709551615"))
             ? "" : _strip_proof_header($qcert);

    # Big GCD's are hugely fast with GMP or Pari, but super slow with Calc.
    _make_big_gcds() if $_big_gcd_use < 0;

    my $loop_limit = 1_000_000 + $k * 1_000;
    while ($loop_limit-- > 0) {
      # R is a random number between $I+1 and 2*$I
      my $R = $I + 1 + $_RANDF->( $I - 1 );
      #my $n = 2 * $R * $q + 1;
      my $n = Math::BigInt->new(2)->bmul($R)->bmul($q)->binc();
      # We constructed a promising looking $n.  Now test it.
      print "." if $verbose > 2;
      if ($_HAVE_GMP) {
        # MPU::GMP::is_prob_prime has fast tests built in.
        next unless Math::Prime::Util::GMP::is_prob_prime($n);
      } else {
        # No GMP, so first do trial divisions, then a SPSP test.
        next unless Math::BigInt::bgcd($n, 111546435) == 1;
        if ($_big_gcd_use && $n > $_big_gcd_top) {
          next unless Math::BigInt::bgcd($n, $_big_gcd[0]) == 1;
          next unless Math::BigInt::bgcd($n, $_big_gcd[1]) == 1;
          next unless Math::BigInt::bgcd($n, $_big_gcd[2]) == 1;
          next unless Math::BigInt::bgcd($n, $_big_gcd[3]) == 1;
        }
        print "+" if $verbose > 2;
        next unless Math::Prime::Util::is_strong_pseudoprime($n, 3);
      }
      print "*" if $verbose > 2;

      # We could set r = 0.3333 earlier, then use BLS75 theorem 5 here.
      # The chain would be shorter, requiring less overall work for
      # large inputs.  Maurer's paper discusses the idea.

      # Use BLS75 theorem 3.  This is easier and possibly faster than
      # BLS75 theorem 4 (Pocklington) used by Maurer's paper.

      # Check conditions -- these should be redundant.
      my $m = 2 * $R;
      if (! ($q->is_odd && $q > 2 && $m > 0 &&
             $m * $q + 1 == $n && 2*$q+1 > $n->copy->bsqrt()) ) {
        carp "Maurer prime failed BLS75 theorem 3 conditions.  Retry.";
        next;
      }
      # Find a suitable a.  Move on if one isn't found quickly.
      foreach my $trya (2, 3, 5, 7, 11, 13) {
        my $a = Math::BigInt->new($trya);
        # m/2 = R    (n-1)/2 = (2*R*q)/2 = R*q
        next unless $a->copy->bmodpow($R, $n) != $n-1;
        next unless $a->copy->bmodpow($R*$q, $n) == $n-1;
        print "($k)" if $verbose > 2;
        croak "Maurer prime $n=2*$R*$q+1 failed BPSW" unless is_prob_prime($n);
        my $cert = "[MPU - Primality Certificate]\nVersion 1.0\n\n" .
                   "Proof for:\nN $n\n\n" .
                   "Type BLS3\nN $n\nQ $q\nA $a\n" .
                   $qcert;
        return ($n, $cert);
      }
      # Didn't pass the selected a values.  Try another R.
    }
    croak "Failure in random_maurer_prime, could not find a prime\n";
  } # End of random_maurer_prime

  # Gordon's algorithm for generating a strong prime.
  sub random_strong_prime {
    my($t) = @_;
    _validate_num($t, 128) || _validate_positive_integer($t, 128);
    croak "Random strong primes must be >= 173 bits on old Perl" if $] < 5.008 && $_Config{'maxbits'} > 32 && $t < 173;

    _set_randf() unless defined $_RANDF;

    my $l   = (($t+1) >> 1) - 2;
    my $lp  = int($t/2) - 20;
    my $lpp = $l - 20;
    while (1) {
      my $qp  = random_nbit_prime($lp);
      my $qpp = random_nbit_prime($lpp);
      $qp  = Math::BigInt->new("$qp")  unless ref($qp)  eq 'Math::BigInt';
      $qpp = Math::BigInt->new("$qpp") unless ref($qpp) eq 'Math::BigInt';
      my ($il, $rem) = Math::BigInt->new(2)->bpow($l-1)->bdec()->bdiv(2*$qpp);
      $il++ if $rem > 0;
      $il = $il->as_int();
      my $iu = Math::BigInt->new(2)->bpow($l)->bsub(2)->bdiv(2*$qpp)->as_int();
      my $istart = $il + $_RANDF->($iu - $il);
      for (my $i = $istart; $i <= $iu; $i++) {  # Search for q
        my $q = 2 * $i * $qpp + 1;
        next unless is_prob_prime($q);
        my $pp = $qp->copy->bmodpow($q-2, $q)->bmul(2)->bmul($qp)->bdec();
        my ($jl, $rem) = Math::BigInt->new(2)->bpow($t-1)->bsub($pp)->bdiv(2*$q*$qp);
        $jl++ if $rem > 0;
        $jl = $jl->as_int();
        my $ju = Math::BigInt->new(2)->bpow($t)->bdec()->bsub($pp)->bdiv(2*$q*$qp)->as_int();
        my $jstart = $jl + $_RANDF->($ju - $jl);
        for (my $j = $jstart; $j <= $ju; $j++) {  # Search for p
          my $p = $pp + 2 * $j * $q * $qp;
          return $p if is_prob_prime($p);
        }
      }
    }
  }

  sub random_proven_prime {
    my $k = shift;
    my ($n, $cert) = random_proven_prime_with_cert($k);
    croak "maurer prime $n failed certificate verification!"
          unless verify_prime($cert);
    return $n;
  }

  sub random_proven_prime_with_cert {
    my $k = shift;
    _validate_num($k, 2) || _validate_positive_integer($k, 2);

    if ($_HAVE_GMP && $k <= 450) {
      my $n = random_nbit_prime($k);
      my ($isp, $cert) = is_provable_prime_with_cert($n);
      croak "small nbit prime could not be proven" if $isp != 2;
      return ($n, $cert);
    }
    return random_maurer_prime_with_cert($k);
  }

} # end of the random prime section

sub primorial {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return Math::BigInt->new(''.Math::Prime::Util::GMP::primorial($n))
    if $_HAVE_GMP && defined &Math::Prime::Util::GMP::primorial;

  my $max = ($_Config{'maxbits'} == 32) ? 29 : 53;
  my $pn = (ref($_[0]) eq 'Math::BigInt') ? $_[0]->copy->bone()
         : ($n >= $max) ? Math::BigInt->bone()
         : 1;
  forprimes { $pn *= $_ } $n;
  return $pn;
}

sub pn_primorial {
  my($n) = @_;

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::pn_primorial) {
    _validate_num($n) || _validate_positive_integer($n);
    return Math::BigInt->new(''.Math::Prime::Util::GMP::pn_primorial($n))
  }

  return primorial(nth_prime($n));
}

sub consecutive_integer_lcm {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  return 0 if $n < 1;

  my $max = ($_Config{'maxbits'} == 32) ? 22 : ($] < 5.008) ? 43 : 46;

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::consecutive_integer_lcm) {
    my $clcm = Math::Prime::Util::GMP::consecutive_integer_lcm($n);
    return ($n < $max) ? int($clcm) : Math::BigInt->new("$clcm");
  }

  my $pn = (ref($_[0]) eq 'Math::BigInt') ? $_[0]->copy->bone()
         : ($n >= $max) ? Math::BigInt->bone()
         : 1;
  forprimes {
    my($p_power, $pmin) = ($_, int($n/$_));
    $p_power *= $_ while $p_power <= $pmin;
    $pn *= $p_power;
  } $n;

  return $pn;
}


sub all_factors {
  my $n = shift;

  # In scalar context, returns sigma_0(n).  Very fast.
  return divisor_sum($n,0) unless wantarray;

  my $use_bigint = !($n < $_Config{'maxparam'} || int($n) eq $_Config{'maxparam'});
  my @factors = factor($n);
  if ($n <= 0) { @factors = (); return @factors; }
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
  # Add 1 and n
  undef $all_factors{1};
  undef $all_factors{$n};
  @factors = sort {$a<=>$b} keys %all_factors;
  return @factors;
}


# A008683 Moebius function mu(n)
# A030059, A013929, A030229, A002321, A005117, A013929 all relate.
sub moebius {
  my($n, $nend) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  if (defined $nend) {
    _validate_num($nend) || _validate_positive_integer($nend);
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
    forprimes {
      my $p = $_;
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
    } $sqrtn;
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
  foreach my $i (1 .. $#factors) {
    return 0 if $factors[$i] == $factors[$i-1];
  }
  return (((scalar @factors) % 2) == 0) ? 1 : -1;
}

# A002321 Mertens' function.  mertens(n) = sum(moebius(1,n))
sub mertens {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);
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
    my ($denom, $this_k, $next_k) = ($m, 0, int($n/($m*1)));
    for my $nmk (1 .. $last_nmk) {
      $denom += $m;
      $this_k = int($n/$denom);
      next if $this_k == $next_k;
      ($this_k, $next_k) = ($next_k, $this_k);
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
  _validate_num($n) || _validate_positive_integer($n);
  if (defined $nend) {
    _validate_num($nend) || _validate_positive_integer($nend);
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
        $totients[$i*$j] -= $totients[$i*$j]/$i;
      }
    }
    splice(@totients, 0, $n) if $n > 0;
    return @totients;
  }

  return $n if $n <= 1;

  my @pe = factor_exp($n);
  my $totient = $n - $n + 1;

  if (ref($n) ne 'Math::BigInt') {
    foreach my $f (@pe) {
      my ($p, $e) = @$f;
      $totient *= $p - 1;
      $totient *= $p for 2 .. $e;
    }
  } else {
    my $zero = $n->copy->bzero;
    foreach my $f (@pe) {
      my ($p, $e) = @$f;
      $p = $zero->copy->badd("$p");
      $totient->bmul($p->copy->bdec());
      $totient->bmul($p) for 2 .. $e;
    }
  }
  return $totient;
}

# Jordan's totient -- a generalization of Euler's totient.
sub jordan_totient {
  my($k, $n) = @_;
  _validate_num($k, 1) || _validate_positive_integer($k, 1);
  return euler_phi($n) if $k == 1;

  return 0 if defined $n && $n <= 0;  # Following SAGE's logic here.
  _validate_num($n) || _validate_positive_integer($n);
  return 1 if $n <= 1;

  my @pe = ($n <= $_XS_MAXVAL) ? _XS_factor_exp($n) : factor_exp($n);

  my $totient = $n - $n + 1;

  if (ref($n) ne 'Math::BigInt') {
    foreach my $f (@pe) {
      my ($p, $e) = @$f;
      my $fmult = int($p ** $k);
      $totient *= ($fmult - 1);
      $totient *= $fmult for (2 .. $e);
    }
  } else {
    my $zero = $n->copy->bzero;
    foreach my $f (@pe) {
      my ($p, $e) = @$f;
      my $fmult = $zero->copy->badd("$p")->bpow($k);
      $totient->bmul($fmult->copy->bdec());
      $totient->bmul($fmult) for (2 .. $e);
    }
  }
  return $totient;
}

my @_ds_overflow =  # We'll use BigInt math if the input is larger than this.
  (~0 > 4294967295)
   ? (124, 3000000000000000000, 3000000000, 2487240, 64260, 7026)
   : ( 50,           845404560,      52560,    1548,   252,   84);
sub divisor_sum {
  my($n, $k) = @_;
  return (0,1)[$n] if defined $n && $n <= 1;
  _validate_num($n) || _validate_positive_integer($n);

  # Call the XS routine for k=0 and k=1 immediately if possible.
  if ($n <= $_XS_MAXVAL) {
    if (defined $k) {
      return _XS_divisor_sum($n, 0) if $k eq '0';
      return _XS_divisor_sum($n, 1) if $k eq '1' && $n < $_ds_overflow[1];
    } else {
      return _XS_divisor_sum($n, 1) if $n < $_ds_overflow[1];
    }
  }

  if (defined $k && ref($k) eq 'CODE') {
    my $sum = $n-$n;
    foreach my $f (all_factors($n)) {
      $sum += $k->($f);
    }
    return $sum;
  }

  croak "Second argument must be a code ref or number"
    unless !defined $k || _validate_num($k) || _validate_positive_integer($k);
  $k = 1 if !defined $k;

  my $will_overflow = ($k == 0) ? (length($n) >= $_ds_overflow[0])
                    : ($k <= 5) ? ($n >= $_ds_overflow[$k])
                    : 1;

  return _XS_divisor_sum($n, $k) if $n <= $_XS_MAXVAL && !$will_overflow;

  # The standard way is:
  #    my $pk = $f ** $k;  $product *= ($pk ** ($e+1) - 1) / ($pk - 1);
  # But we get less overflow using:
  #    my $pk = $f ** $k;  $product *= $pk**E for E in 0 .. e
  # Also separate BigInt and do fiddly bits for better performance.

  my $product = 1;
  my @pe = ($n <= $_XS_MAXVAL) ? _XS_factor_exp($n) : factor_exp($n);

  if (!$will_overflow) {
    foreach my $f (@pe) {
      my ($p, $e) = @$f;
      if ($k == 0) {
        $product *= ($e+1);
      } else {
        my $pk = $p ** $k;
        my $fmult = $pk + 1;
        foreach my $E (2 .. $e) { $fmult += $pk**$E }
        $product *= $fmult;
      }
    }
  } else {
    $product = Math::BigInt->bone;
    my $bik = Math::BigInt->new("$k");
    foreach my $f (@pe) {
      my ($p, $e) = @$f;
      my $pk = Math::BigInt->new("$p")->bpow($bik);
      if    ($e == 1) { $pk->binc(); $product->bmul($pk); }
      elsif ($e == 2) { $pk->badd($pk*$pk)->binc(); $product->bmul($pk); }
      else {
        my $fmult = $pk;
        foreach my $E (2 .. $e) { $fmult += $pk->copy->bpow($E) }
        $fmult->binc();
        $product *= $fmult;
      }
    }
  }
  return $product;
}

                                   # Need proto for the block
sub _generic_forprimes (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, $beg, $end) = @_;
  if (!defined $end) { $end = $beg; $beg = 2; }
  _validate_num($beg) || _validate_positive_integer($beg);
  _validate_num($end) || _validate_positive_integer($end);
  $beg = 2 if $beg < 2;
  {
    my $pp;
    local *_ = \$pp;
    for (my $p = next_prime($beg-1);  $p <= $end;  $p = next_prime($p)) {
      $pp = $p;
      $sub->();
    }
  }
}

sub prime_iterator {
  my($start) = @_;
  $start = 0 unless defined $start;
  _validate_num($start) || _validate_positive_integer($start);
  my $p = ($start > 0) ? $start-1 : 0;
  if (ref($p) ne 'Math::BigInt' && $p <= $_XS_MAXVAL) {
    return sub { $p = _XS_next_prime($p); return $p; };
  } elsif ($_HAVE_GMP) {
    return sub { $p = $p-$p+Math::Prime::Util::GMP::next_prime($p); return $p;};
  } else {
    return sub { $p = Math::Prime::Util::PP::next_prime($p); return $p; }
  }
}

sub prime_iterator_object {
  my($start) = @_;
  eval { require Math::Prime::Util::PrimeIterator; 1; }
  or do { croak "Cannot load Math::Prime::Util::PrimeIterator"; };
  return Math::Prime::Util::PrimeIterator->new($start);
}

# Exponential of Mangoldt function (A014963).
# Return p if n = p^m [p prime, m >= 1], 1 otherwise.
sub _generic_exp_mangoldt {
  my($n) = @_;
  return 1 if defined $n && $n <= 1;  # n <= 1
  _validate_num($n) || _validate_positive_integer($n);

  return 2 if ($n & ($n-1)) == 0;     # n power of 2
  return 1 unless $n & 1;             # even n can't be p^m

  my @pe = ($n <= $_XS_MAXVAL) ? _XS_factor_exp($n) : factor_exp($n);
  return 1 if scalar @pe > 1;
  return $pe[0]->[0];
}

sub liouville {
  my($n) = @_;
  # Note the special behavior for n = 1
  return 1 if defined $n && $n <= 1;
  _validate_num($n) || _validate_positive_integer($n);
  my $l = ($n <= $_XS_MAXVAL)  ?  (-1) ** scalar _XS_factor($n)
                               :  (-1) ** scalar factor($n);
  return $l;
}

# See 2011+ FLINT and Fredrik Johansson's work for state of the art.
# Very crude timing estimates (ignores growth rates).
#   Perl-comb   partitions(10^5)  ~ 300 seconds  ~200,000x slower than Pari
#   GMP-comb    partitions(10^6)  ~ 120 seconds    ~1,000x slower than Pari
#   Pari        partitions(10^8)  ~ 100 seconds
#   Bober 0.6   partitions(10^9)  ~  20 seconds       ~50x faster than Pari
#   Johansson   partitions(10^12) ~  10 seconds     >1000x faster than Pari
sub partitions {
  my($n) = @_;
  return 1 if defined $n && $n <= 0;
  _validate_num($n) || _validate_positive_integer($n);

  return Math::BigInt->new(''.Math::Prime::Util::GMP::partitions($n))
    if $_HAVE_GMP && defined &Math::Prime::Util::GMP::partitions;

  my $d = int(sqrt($n+1));
  my @pent = (1, map { (($_*(3*$_+1))>>1, (($_+1)*(3*$_+2))>>1) } 1 .. $d);
  my @part = (Math::BigInt->bone);
  foreach my $j (scalar @part .. $n) {
    my ($psum1, $psum2, $k) = (Math::BigInt->bzero, Math::BigInt->bzero, 1);
    foreach my $p (@pent) {
      last if $p > $j;
      if ((++$k) & 2) { $psum1->badd( $part[ $j - $p ] ); }
      else            { $psum2->badd( $part[ $j - $p ] ); }
    }
    $part[$j] = $psum1 - $psum2;
  }
  return $part[$n];
}

sub chebyshev_theta {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  return _XS_chebyshev_theta($n) if $n <= $_XS_MAXVAL;
  my $sum = 0.0;
  forprimes { $sum += log($_); } $n;
  return $sum;
}
sub chebyshev_psi {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  return 0 if $n <= 1;
  return _XS_chebyshev_psi($n) if $n <= $_XS_MAXVAL;
  my ($sum, $logn, $sqrtn) = (0.0, log($n), int(sqrt($n)));
  forprimes {
    my $logp = log($_);
    $sum += $logp * int($logn/$logp+1e-15);
  } $sqrtn;
  forprimes {
    $sum += log($_);
  } $sqrtn+1, $n;
  return $sum;
}

sub carmichael_lambda {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  # lambda(n) = phi(n) for n < 8
  return euler_phi($n) if $n < 8;
  # lambda(n) = phi(n)/2 for powers of two greater than 4
  return euler_phi($n)/2 if ($n & ($n-1)) == 0;

  my @pe = ($n <= $_XS_MAXVAL) ? _XS_factor_exp($n) : factor_exp($n);
  $pe[0]->[1]-- if $pe[0]->[0] == 2 && $pe[0]->[1] > 2;

  my $lcm = Math::BigInt::blcm(
    map { $_->[0]->copy->bpow($_->[1]->copy->bdec)->bmul($_->[0]->copy->bdec) }
    map { [ map { Math::BigInt->new("$_") } @$_ ] }
    @pe
  );
  $lcm = int($lcm->bstr) if $lcm->bacmp(''.~0) <= 0;
  return $lcm;
}

sub znorder {
  my($a, $n) = @_;
  _validate_num($a) || _validate_positive_integer($a);
  _validate_num($n) || _validate_positive_integer($n);
  return if $n <= 0;
  return (undef,1)[$a] if $a <= 1;
  return 1 if $n == 1;

  # Sadly, Calc/FastCalc are horrendously slow for this function.
  return if Math::BigInt::bgcd($a, $n) > 1;
  # Method 1:  check all a^k 1 .. $n-1.
  #            Naive and terrible slow.
  # Method 2:  check all k in (all_factors(euler_phi($n), $n).
  # Method 3:  check all k in (all_factors(carmichael_lambda($n), $n).
  #            Good for most cases, but slow when factor quantity is large.
  # Method 4:  Das algorithm 1.7, just enough multiples of p
  #            Fastest.
  #
  # Most of the time is spent factoring $n and $phi.  We could do the phi
  # construction here and build its factors to save a little more time.

  $a = Math::BigInt->new("$a") unless ref($a) eq 'Math::BigInt';
  my $phi = Math::BigInt->new('' . euler_phi($n));
  my @pe = ($phi <= $_XS_MAXVAL) ? _XS_factor_exp($phi) : factor_exp($phi);
  my $k = Math::BigInt->bone;
  foreach my $f (@pe) {
    my($pi, $ei, $enum) = (Math::BigInt->new("$f->[0]"), $f->[1], 0);
    my $phidiv = $phi / ($pi**$ei);
    my $b = $a->copy->bmodpow($phidiv, $n);
    while ($b != 1) {
      return if $enum++ >= $ei;
      $b = $b->copy->bmodpow($pi, $n);
      $k *= $pi;
    }
  }
  return $k;

  # Method 3:
  # my $cl = carmichael_lambda($n);
  # foreach my $k (all_factors($cl), $cl) {
  #   my $t = $a->copy->bmodpow("$k", $n);
  #   return $k if $t->is_one;
  # }
  # return;
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
# on small inputs typically take a lot less time than this.  So the overhead
# for these is significantly more than just the XS call itself.  For some
# functions we have inverted the operation, so the XS function gets called,
# does validation, and calls the _generic_* sub here if it doesn't know what
# to do.  This removes all the overhead, making functions like is_prime run
# about 5x faster for very small inputs.

sub _generic_is_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_is_prime($n) if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_prime($n) if $_HAVE_GMP;

  if ($n < 7) { return ($n == 2) || ($n == 3) || ($n == 5) ? 2 : 0; }
  return 0 if !($n % 2) || !($n % 3) || !($n % 5);
  return Math::Prime::Util::PP::_is_prime7($n);
}

sub _generic_is_prob_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_is_prob_prime($n) if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_prob_prime($n) if $_HAVE_GMP;

  if ($n < 7) { return ($n == 2) || ($n == 3) || ($n == 5) ? 2 : 0; }
  return 0 if !($n % 2) || !($n % 3) || !($n % 5);
  return Math::Prime::Util::PP::_is_prime7($n);
}

sub _generic_next_prime {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_next_prime($n)
    if $n <= $_XS_MAXVAL && $n < $_Config{'maxprime'};

  if ($_HAVE_GMP) {
    my $r = Math::Prime::Util::GMP::next_prime($n);
    return (ref($n) eq 'Math::BigInt' || $n >= $_Config{'maxprime'})
           ?  Math::BigInt->new("$r")  :  int($r);
  }

  return Math::Prime::Util::PP::next_prime($_[0]);
}

sub _generic_prev_prime {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_prev_prime($n)
    if $n <= $_XS_MAXVAL;

  if ($_HAVE_GMP) {
    my $r = Math::Prime::Util::GMP::prev_prime($n);
    return (ref($n) eq 'Math::BigInt' && $r > $_Config{'maxprime'})
           ?  Math::BigInt->new("$r")  :  int($r);
  }

  return Math::Prime::Util::PP::prev_prime($_[0]);
}

sub prime_count {
  my($low,$high) = @_;
  if (defined $high) {
    _validate_num($low) || _validate_positive_integer($low);
    _validate_num($high) || _validate_positive_integer($high);
  } else {
    ($low,$high) = (2, $low);
    _validate_num($high) || _validate_positive_integer($high);
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
        my $count;
        $count  =  _XS_LMO_pi($high);
        $count -=  _XS_LMO_pi($low-1) if $low > 2;
        return $count;
      }
    }
    return _XS_prime_count($low,$high);
  }
  # We can relax these constraints if MPU::GMP gets a fast implementation.
  return Math::Prime::Util::GMP::prime_count($low,$high) if $_HAVE_GMP
                       && defined &Math::Prime::Util::GMP::prime_count
                       && (   (ref($high) eq 'Math::BigInt')
                           || (($high-$low) < int($low/1_000_000))
                          );
  return Math::Prime::Util::PP::prime_count($low,$high);
}

sub nth_prime {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_nth_prime($n)
    if $n <= $_XS_MAXVAL && $n < $_Config{'maxprimeidx'};

  return Math::Prime::Util::PP::nth_prime($n);
}

sub factor {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_factor($n) if $n <= $_XS_MAXVAL;

  if ($_HAVE_GMP) {
    my @factors = Math::Prime::Util::GMP::factor($n);
    if (ref($_[0]) eq 'Math::BigInt') {
      @factors = map { ($_ > ~0) ? Math::BigInt->new(''.$_) : $_ } @factors;
    }
    return @factors;
  }

  return Math::Prime::Util::PP::factor($n);
}

sub factor_exp {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_factor_exp($n) if ref($n) ne 'Math::BigInt' && $n <= $_XS_MAXVAL;

  my %exponents;
  my @factors = grep { !$exponents{$_}++ } factor($n);
  return scalar @factors unless wantarray;
  return (map { [$_, $exponents{$_}] } @factors);
}

#sub trial_factor {
#  my($n, $limit) = @_;
#  _validate_num($n) || _validate_positive_integer($n);
#  $limit = 2147483647 unless defined $limit;
#  _validate_num($limit) || _validate_positive_integer($limit);
#  return _XS_trial_factor($n, $limit) if $n <= $_XS_MAXVAL;
#  if ($_HAVE_GMP) {
#    my @factors;
#    while (1) {
#      last if $n <= 1 || is_prob_prime($n);
#      my @f = sort { $a <=> $b }
#              Math::Prime::Util::GMP::trial_factor($n, $limit);
#      pop(@f);  # Remove remainder
#      last unless scalar @f > 0;
#      foreach my $f (@f) {
#        push @factors, $f;
#        $n /= $f;
#      }
#    }
#    push @factors, $n if $n > 1;
#    return @factors;
#  }
#  return Math::Prime::Util::PP::trial_factor($n, $limit);
#}

sub is_pseudoprime {
  my($n, $a) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  _validate_num($a) || _validate_positive_integer($a);
  return _XS_is_pseudoprime($n, $a)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_pseudoprime($n, $a)
    if $_HAVE_GMP && defined &Math::Prime::Util::GMP::is_pseudoprime;
  return Math::Prime::Util::PP::is_pseudoprime($n, $a);
}

sub is_strong_pseudoprime {
  my($n) = shift;
  _validate_num($n) || _validate_positive_integer($n);
  # validate bases?
  return _XS_miller_rabin($n, @_)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_strong_pseudoprime($n, @_) if $_HAVE_GMP;
  return Math::Prime::Util::PP::miller_rabin($n, @_);
}

sub is_lucas_pseudoprime {
  my($n) = shift;
  _validate_num($n) || _validate_positive_integer($n);
  return _XS_is_lucas_pseudoprime($n)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_lucas_pseudoprime("$n")
    if $_HAVE_GMP && defined &Math::Prime::Util::GMP::is_lucas_pseudoprime;
  return Math::Prime::Util::PP::is_lucas_pseudoprime($n);
}

sub is_strong_lucas_pseudoprime {
  my($n) = shift;
  _validate_num($n) || _validate_positive_integer($n);
  return _XS_is_strong_lucas_pseudoprime($n)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_strong_lucas_pseudoprime("$n")
    if $_HAVE_GMP;
  return Math::Prime::Util::PP::is_strong_lucas_pseudoprime($n);
}

sub is_extra_strong_lucas_pseudoprime {
  my($n) = shift;
  _validate_num($n) || _validate_positive_integer($n);
  return _XS_is_extra_strong_lucas_pseudoprime($n)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_extra_strong_lucas_pseudoprime("$n")
    if $_HAVE_GMP
    && defined &Math::Prime::Util::GMP::is_extra_strong_lucas_pseudoprime;
  return Math::Prime::Util::PP::is_extra_strong_lucas_pseudoprime($n);
}

sub is_almost_extra_strong_lucas_pseudoprime {
  my($n, $inc) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  if (!defined $inc) {
    $inc = 1;
  } elsif ($inc ne '1' && $inc ne '2') {
    (_validate_num($inc) && $inc > 0 && $inc < 257)
    || _validate_positive_integer($inc, 1, 256);
  }
  return _XS_is_almost_extra_strong_lucas_pseudoprime($n, $inc)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_almost_extra_strong_lucas_pseudoprime("$n", $inc)
    if $_HAVE_GMP
    && defined &Math::Prime::Util::GMP::is_almost_extra_strong_lucas_pseudoprime;
  return Math::Prime::Util::PP::is_almost_extra_strong_lucas_pseudoprime($n, $inc);
}

sub is_frobenius_underwood_pseudoprime {
  my($n) = shift;
  _validate_num($n) || _validate_positive_integer($n);
  return _XS_is_frobenius_underwood_pseudoprime($n)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_frobenius_underwood_pseudoprime("$n")
    if $_HAVE_GMP
    && defined &Math::Prime::Util::GMP::is_frobenius_underwood_pseudoprime;
  return Math::Prime::Util::PP::is_frobenius_underwood_pseudoprime($n);
}

sub miller_rabin {
  #warn "miller_rabin() is deprecated. Use is_strong_pseudoprime instead.";
  return is_strong_pseudoprime(@_);
}

sub lucas_sequence {
  my($n, $P, $Q, $k) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  _validate_num($k) || _validate_positive_integer($k);
  { my $testP = (!defined $P || $P >= 0) ? $P : -$P;
    _validate_num($testP) || _validate_positive_integer($testP); }
  { my $testQ = (!defined $Q || $Q >= 0) ? $Q : -$Q;
    _validate_num($testQ) || _validate_positive_integer($testQ); }
  return _XS_lucas_sequence($n, $P, $Q, $k)
    if ref($_[0]) ne 'Math::BigInt' && $n <= $_XS_MAXVAL
    && ref($_[3]) ne 'Math::BigInt' && $k <= $_XS_MAXVAL;
  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::lucas_sequence) {
    return map { ($_ > ~0) ? Math::BigInt->new(''.$_) : $_ }
           Math::Prime::Util::GMP::lucas_sequence($n, $P, $Q, $k);
  }
  return Math::Prime::Util::PP::lucas_sequence($n, $P, $Q, $k);
}

sub miller_rabin_random {
  my($n, $k, $seed) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  _validate_num($k) || _validate_positive_integer($k);

  return 1 if $k <= 0;
  return (is_prob_prime($n) > 0) if $n < 100;
  return 0 unless $n & 1;
  return Math::Prime::Util::GMP::miller_rabin_random($n, $k)
    if $_HAVE_GMP
    && defined &Math::Prime::Util::GMP::miller_rabin_random;

  # Testing this many bases is silly, but let's pretend they have some
  # good reason.  A composite n > 9 must have at least n/4 witnesses,
  # hence we need to check only floor(3/4)+1 at most.  We could improve
  # this is $_Config{'assume_rh'} is true, to 1 .. 2(logn)^2.
  if ($k >= int(3*$n/4)) {
    return is_strong_pseudoprime($n, 2 .. int(3*$n/4)+1+2 );
  }

  my $brange = $n-3;
  my $irandf = _get_randf();
  # Do one first before doing batches
  return 0 unless is_strong_pseudoprime($n, $irandf->($brange)+2 );
  $k--;
  while ($k > 0) {
    my $nbases = ($k >= 20) ? 20 : $k;
    my @bases = map { $irandf->($brange)+2 } 1..$nbases;
    if ($n <= $_XS_MAXVAL) {
      return 0 unless _XS_miller_rabin($n, @bases);
    } else {
      return 0 unless is_strong_pseudoprime($n, @bases);
    }
    $k -= $nbases;
  }
  1;
}


#############################################################################

  # Simple timings, with 0.32 code (Montgomery Reduction for 64-bit)
  # my $n=2**20+1;  do { is_strong_pseudoprime($n,2); $n+=2 } for 1..1000000;
  # my $n=2**47+1;  do { is_strong_pseudoprime($n,2); $n+=2 } for 1..1000000;
  #
  #      1.0 uS  XS 32-bit input, is_strong_pseudoprime
  #      1.8 uS  XS 64-bit input, is_strong_pseudoprime
  #      1.3 uS  XS 32-bit input, is_strong_lucas_pseudoprime
  #      1.9 uS  XS 64-bit input, is_strong_lucas_pseudoprime
  #      1.2 uS  XS 32-bit input, is_almost_extra_strong_lucas_pseudoprime
  #      1.8 uS  XS 64-bit input, is_almost_extra_strong_lucas_pseudoprime
  #
  #     13 uS  Perl 32-bit input, is_strong_pseudoprime
  #    945 uS  Perl 64-bit input, is_strong_pseudoprime
  #   1700 uS  Perl 32-bit input, is_strong_lucas_pseudoprime
  #   3000 uS  Perl 64-bit input, is_strong_lucas_pseudoprime
  #   1500 uS  Perl 32-bit input, is_almost_extra_strong_lucas_pseudoprime
  #   3400 uS  Perl 64-bit input, is_almost_extra_strong_lucas_pseudoprime

sub is_aks_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_is_aks_prime($n) if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_aks_prime($n) if $_HAVE_GMP
                       && defined &Math::Prime::Util::GMP::is_aks_prime;
  return Math::Prime::Util::PP::is_aks_prime($n);
}

# For stripping off the header on certificates so they can be combined.
sub _strip_proof_header {
  my $proof = shift;
  $proof =~ s/^\[MPU - Primality Certificate\]\nVersion \S+\n+Proof for:\nN (\d+)\n+//ms;
  return $proof;
}

# Return just the non-cert portion.
sub is_provable_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_num($n) || _validate_positive_integer($n);

  return _XS_is_prime($n)
    if $n <= $_XS_MAXVAL;
  return Math::Prime::Util::GMP::is_provable_prime($n)
         if $_HAVE_GMP && defined &Math::Prime::Util::GMP::is_provable_prime;

  my ($is_prime, $cert) = is_provable_prime_with_cert($n);
  return $is_prime;
}

# Return just the cert portion.
sub prime_certificate {
  my($n) = @_;
  my ($is_prime, $cert) = is_provable_prime_with_cert($n);
  return $cert;
}


sub is_provable_prime_with_cert {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_num($n) || _validate_positive_integer($n);
  my $header = "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN $n\n\n";

  if ($n <= $_XS_MAXVAL) {
    my $isp = _XS_is_prime("$n");
    return ($isp, '') unless $isp == 2;
    return (2, "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN $n\n\nType Small\nN $n\n");
  }

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::is_provable_prime_with_cert) {
    my ($isp, $cert) = Math::Prime::Util::GMP::is_provable_prime_with_cert($n);
    # New version that returns string format.
    #return ($isp, $cert) if ref($cert) ne 'ARRAY';
    if (ref($cert) ne 'ARRAY') {
      # Fix silly 0.13 mistake (TODO: deprecate this)
      $cert =~ s/^Type Small\n(\d+)/Type Small\nN $1/smg;
      return ($isp, $cert);
    }
    # Old version.  Convert.
    if (!defined $Math::Prime::Util::PrimalityProving::VERSION) {
      eval { require Math::Prime::Util::PrimalityProving; 1; }
      or do { croak "Cannot load Math::Prime::Util::PrimalityProving"; };
    }
    return ($isp, Math::Prime::Util::PrimalityProving::convert_array_cert_to_string($cert));
  }

  {
    my $isp = is_prob_prime($n);
    return ($isp, '') if $isp == 0;
    return (2, "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN $n\n\nType Small\nN $n\n") if $isp == 2;
  }

  # Choice of methods for proof:
  #   ECPP         needs a fair bit of programming work
  #   APRCL        needs a lot of programming work
  #   BLS75 combo  Corollary 11 of BLS75.  Trial factor n-1 and n+1 to B, find
  #                factors F1 of n-1 and F2 of n+1.  Quit when:
  #                B > (N/(F1*F1*(F2/2)))^1/3 or B > (N/((F1/2)*F2*F2))^1/3
  #   BLS75 n+1    Requires factoring n+1 to (n/2)^1/3 (theorem 19)
  #   BLS75 n-1    Requires factoring n-1 to (n/2)^1/3 (theorem 5 or 7)
  #   Pocklington  Requires factoring n-1 to n^1/2 (BLS75 theorem 4)
  #   Lucas        Easy, requires factoring of n-1 (BLS75 theorem 1)
  #   AKS          horribly slow
  # See http://primes.utm.edu/prove/merged.html or other sources.

  if (!defined $Math::Prime::Util::PrimalityProving::VERSION) {
    eval { require Math::Prime::Util::PrimalityProving; 1; }
    or do { croak "Cannot load Math::Prime::Util::PrimalityProving"; };
  }

  #my ($isp, $pref) = Math::Prime::Util::PrimalityProving::primality_proof_lucas($n);
  my ($isp, $pref) = Math::Prime::Util::PrimalityProving::primality_proof_bls75($n);
  carp "proved $n is not prime\n" if !$isp;
  return ($isp, $pref);
}


sub verify_prime {
  my @cdata = @_;

  if (!defined $Math::Prime::Util::PrimalityProving::VERSION) {
    eval { require Math::Prime::Util::PrimalityProving; 1; }
    or do { croak "Cannot load Math::Prime::Util::PrimalityProving"; };
  }

  my $cert = '';
  if (scalar @cdata == 1 && ref($cdata[0]) eq '') {
    $cert = $cdata[0];
  } else {
    # We've been given an old array cert
    $cert = Math::Prime::Util::PrimalityProving::convert_array_cert_to_string(@cdata);
    if ($cert eq '') {
      print "primality fail: error converting old certificate" if $_Config{'verbose'};
      return 0;
    }
  }
  return 0 if $cert eq '';
  return Math::Prime::Util::PrimalityProving::verify_cert($cert);
}


#############################################################################

sub prime_count_approx {
  my($x) = @_;
  _validate_num($x) || _validate_positive_integer($x);

  # With XS using 30k tables, this is super fast.
  return _XS_prime_count(2,$x) if $x < $_XS_MAXVAL && $x < 3_000_000;
  # Give an exact answer for what we have in our little table.
  return _tiny_prime_count($x) if $x < $_primes_small[-1];

  # Below 2^58th or so, all differences between the high precision and C double
  # precision result are less than 0.5.
  if ($x <= $_XS_MAXVAL && $x <= 144115188075855872) {
    return int(_XS_RiemannR($x) + 0.5);
  }

  # Turn on high precision FP if they gave us a big number.
  $x = _upgrade_to_float($x) if ref($_[0]) eq 'Math::BigInt';

  #    Method             10^10 %error  10^19 %error
  #    -----------------  ------------  ------------
  #    n/(log(n)-1)        .22%          .06%
  #    average bounds      .01%          .0002%
  #    li(n)               .0007%        .00000004%
  #    li(n)-li(n^.5)/2    .0004%        .00000001%
  #    R(n)                .0004%        .00000001%
  #
  # Also consider: http://trac.sagemath.org/sage_trac/ticket/8135

  # my $result = int( (prime_count_upper($x) + prime_count_lower($x)) / 2);
  # my $result = int( LogarithmicIntegral($x) );
  # my $result = int(LogarithmicIntegral($x) - LogarithmicIntegral(sqrt($x))/2);
  # my $result = RiemannR($x) + 0.5;

  # Sadly my Perl RiemannR function is really slow for big values.  If MPFR
  # is available, then use it -- it rocks.  Otherwise, switch to LiCorr for
  # very big values.  This is hacky and shouldn't be necessary.
  my $result;
  if ( $x < 1e36 || Math::Prime::Util::PP::_MPFR_available() ) {
    if (ref($x) eq 'Math::BigFloat') {
      # Make sure we get enough accuracy, and also not too much more than needed
      $x->accuracy(length($x->bfloor->bstr())+2);
    }
    $result = RiemannR($x) + 0.5;
  } else {
    $result = int(LogarithmicIntegral($x) - LogarithmicIntegral(sqrt($x))/2);
  }

  return Math::BigInt->new($result->bfloor->bstr()) if ref($result) eq 'Math::BigFloat';
  return int($result);
}

sub prime_count_lower {
  my($x) = @_;
  _validate_num($x) || _validate_positive_integer($x);

  # With XS using 30k tables, this is super fast.
  return _XS_prime_count(2,$x) if $x < $_XS_MAXVAL && $x < 3_000_000;
  # Give an exact answer for what we have in our little table.
  return _tiny_prime_count($x) if $x < $_primes_small[-1];

  $x = _upgrade_to_float($x)
    if ref($x) eq 'Math::BigInt' || ref($_[0]) eq 'Math::BigInt';

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
  _validate_num($x) || _validate_positive_integer($x);

  # With XS using 30k tables, this is super fast.
  return _XS_prime_count(2,$x) if $x < $_XS_MAXVAL && $x < 3_000_000;
  # Give an exact answer for what we have in our little table.
  return _tiny_prime_count($x) if $x < $_primes_small[-1];

  $x = _upgrade_to_float($x)
    if ref($x) eq 'Math::BigInt' || ref($_[0]) eq 'Math::BigInt';

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
    $result = ($x/$flogx) * (1.0 + 1.0/$flogx + $a/($flogx*$flogx)) + 1.0;
  }

  return Math::BigInt->new($result->bfloor->bstr()) if ref($result) eq 'Math::BigFloat';
  return int($result);
}

#############################################################################

sub nth_prime_approx {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n)
    if ref($n) eq 'Math::BigInt' || $n >= $_Config{'maxprimeidx'};

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

  return int($approx + 0.5);
}

# The nth prime will be greater than or equal to this number
sub nth_prime_lower {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if $n > $_Config{'maxprimeidx'} || $n > 2**45;

  my $flogn  = log($n);
  my $flog2n = log($flogn);  # Note distinction between log_2(n) and log^2(n)

  # Dusart 1999 page 14, for all n >= 2
  #my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.25)/$flogn));
  # Dusart 2010 page 2, for all n >= 3
  my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.10)/$flogn));

  return int($lower);
}

# The nth prime will be less or equal to this number
sub nth_prime_upper {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if $n > $_Config{'maxprimeidx'} || $n > 2**45;

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

=for stopwords forprimes Möbius Deléglise totient moebius mertens liouville znorder irand primesieve uniqued k-tuples von SoE pari yafu fonction qui compte le nombre nombres voor PhD superset sqrt(N) gcd(A^M k-th (10001st primegen libtommath


=head1 NAME

Math::Prime::Util - Utilities related to prime numbers, including fast sieves and factoring


=head1 VERSION

Version 0.35


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

  # You can do something for every prime in a range.  Twin primes to 10k:
  forprimes { say if is_prime($_+2) } 10000;

  # For non-bigints, is_prime and is_prob_prime will always be 0 or 2.
  # They return 0 (composite), 2 (prime), or 1 (probably prime)
  say "$n is prime"  if is_prime($n);
  say "$n is ", (qw(composite maybe_prime? prime))[is_prob_prime($n)];

  # Strong pseudoprime test with multiple bases, using Miller-Rabin
  say "$n is a prime or 2/7/61-psp" if is_strong_pseudoprime($n, 2, 7, 61);

  # Standard and strong Lucas-Selfridge, and extra strong Lucas tests
  say "$n is a prime or lpsp"   if is_lucas_pseudoprime($n);
  say "$n is a prime or slpsp"  if is_strong_lucas_pseudoprime($n);
  say "$n is a prime or eslpsp" if is_extra_strong_lucas_pseudoprime($n);

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

  # Return ([p1,e1],[p2,e2], ...) for $n = p1^e1 * p2*e2 * ...
  @pe = factor_exp( $n );

  # Get all divisors other than 1 and n
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
  # Some more number theoretical functions
  say liouville(4292384);
  say chebyshev_psi(234984);
  say chebyshev_theta(92384234);
  say partitions(1000);

  # divisor sum
  $sigma  = divisor_sum( $n );       # sum of divisors
  $sigma0 = divisor_sum( $n, 0 );    # count of divisors
  $sigmak = divisor_sum( $n, $k );
  $sigmaf = divisor_sum( $n, sub { log($_[0]) } ); # arbitrary func

  # primorial n#, primorial p(n)#, and lcm
  say "The product of primes below 47 is ",     primorial(47);
  say "The product of the first 47 primes is ", pn_primorial(47);
  say "lcm(1..1000) is ", consecutive_integer_lcm(1000);

  # Ei, li, and Riemann R functions
  my $ei   = ExponentialIntegral($x);   # $x a real: $x != 0
  my $li   = LogarithmicIntegral($x);   # $x a real: $x >= 0
  my $R    = RiemannR($x)               # $x a real: $x > 0
  my $Zeta = RiemannZeta($x)            # $x a real: $x >= 0


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

All operations support both Perl UV's (32-bit or 64-bit) and bignums.  If
you want high performance with big numbers (larger than Perl's native 32-bit
or 64-bit size), you should install L<Math::Prime::Util::GMP> and
L<Math::BigInt::GMP>.  This will be a recurring theme throughout this
documentation -- while all bignum operations are supported in pure Perl,
most methods will be much slower than the C+GMP alternative.

The module is thread-safe and allows concurrency between Perl threads while
still sharing a prime cache.  It is not itself multi-threaded.  See the
L<Limitations|/"LIMITATIONS"> section if you are using Win32 and threads in
your program.

Two scripts are also included and installed by default:

=over 4

=item *

primes.pl displays primes between start and end values or expressions,
with many options for filtering (e.g. twin, safe, circular, good, lucky,
etc.).  Use C<--help> to see all the options.

=item *

factor.pl operates similar to the GNU C<factor> program.  It supports
bigint and expression inputs.

=back


=head1 BIGNUM SUPPORT

By default all functions support bignums.  For performance, you should
install and use L<Math::BigInt::GMP> or L<Math::BigInt::Pari>, and
L<Math::Prime::Util::GMP>.

Some of the functions, including:

  factor
  factor_exp
  is_pseudoprime
  is_strong_pseudoprime
  nth_prime
  moebius
  mertens
  euler_phi
  chebyshev_theta
  chebyshev_psi
  is_prime
  is_prob_prime
  next_prime
  prev_prime

work very fast (under 1 microsecond) on small inputs, but the wrappers for
input validation and bigint support take more time than the function itself.
Using the flag '-bigint', e.g.:

  use Math::Prime::Util qw(-bigint);

will turn off bigint support for those functions.  Those functions will then
go directly to the XS versions, which will speed up very small inputs a B<lot>.
This is useful if you're using the functions in a loop, but since the difference
is less than a millisecond, it's really not important in general.  The last
five functions have shortcuts by default so will only skip validation.


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

I have run these functions on many versions of Perl, and my experience is that
if you're using anything older than Perl 5.14, I would recommend you upgrade
if you are using bignums a lot.  There are some brittle behaviors on 5.12.4
and earlier with bignums.

=back


=head1 PRIMALITY TESTING

This module provides three functions for general primality testing, as
well as numerous specialized functions.  The three main functions are:
L</is_prob_prime> and L</is_prime> for general use, and L</is_provable_prime>
for proofs.  For inputs below C<2^64> the functions are identical and
fast deterministic testing is performed.  That is, the results will always
be correct and should take at most a few microseconds for any input.  This
is hundreds to thousands of times faster than other CPAN modules.  For
inputs larger than C<2^64>, an extra-strong
L<BPSW test|http://en.wikipedia.org/wiki/Baillie-PSW_primality_test>
is used.  See the L</PRIMALITY TESTING NOTES> section for more
discussion.


=head1 FUNCTIONS

=head2 is_prime

  print "$n is prime" if is_prime($n);

Returns 2 if the number is prime, 0 if not.  For numbers larger than C<2^64>
it will return 0 for composite and 1 for probably prime, using an
extra-strong BPSW test.  If L<Math::Prime::Util::GMP> is installed, some
additional primality tests are also performed on large inputs, and a
quick attempt is made to perform a primality proof, so it will
return 2 for many other inputs.

Also see the L</is_prob_prime> function, which will never do additional
tests, and the L</is_provable_prime> function which will construct a proof
that the input is number prime and returns 2 for almost all primes (at the
expense of speed).

For native precision numbers (anything smaller than C<2^64>, all three
functions are identical and use a deterministic set of tests (selected
Miller-Rabin bases or BPSW).  For larger inputs both L</is_prob_prime> and
L</is_prime> return probable prime results using the extra-strong
Baillie-PSW test, which has had no counterexample found since it was
published in 1980.

For cryptographic key generation, you may want even more testing for probable
primes (NIST recommends some additional M-R tests).  This can be done using
additional random bases with L</is_strong_pseudoprime>, or a different test
such as L</is_frobenius_underwood_pseudoprime>.  Even better, make sure
L<Math::Prime::Util::GMP> is installed and use L</is_provable_prime> which
should be reasonably fast for sizes under 2048 bits.
Another possibility is to use L<Math::Prime::Util/random_maurer_prime> which
constructs a random provable prime.


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

Returns the prime preceding the input number (i.e. the largest prime that is
strictly less than the input).  0 is returned if the input is C<2> or lower.


=head2 forprimes

  forprimes { say } 100,200;                  # print primes from 100 to 200

  $sum=0;  forprimes { $sum += $_ } 100000;   # sum primes to 100k

  forprimes { say if is_prime($_+2) } 10000;  # print twin primes to 10k

Given a block and either an end count or a start and end pair, calls the
block for each prime in the range.  Compared to getting a big array of primes
and iterating through it, this is more memory efficient and perhaps more
convenient.  This will almost always be the fastest way to loop over a range
of primes.  Nesting and using in threads are allowed.

Math::BigInt objects may be used for the range.

For some uses an iterator (L</prime_iterator>, L</prime_iterator_object>)
or a tied array (L<Math::Prime::Util::PrimeArray>) may be more convenient.
Objects can be passed to functions, and allow early loop exits without
exceptions.  Here is a clumsy L</forprimes> exception example:

  use bigint;
  eval { forprimes { die "$_\n" if $_ % 123 == 1 } 2**100, 2**101 };
  my $n = 0+$@;


=head2 prime_iterator

  my $it = prime_iterator;
  $sum += $it->() for 1..100000;

Returns a closure-style iterator.  The start value defaults to the first
prime (2) but an initial value may be given as an argument, which will result
in the first value returned being the next prime greater than or equal to the
argument.  For example, this:

  my $it = prime_iterator(200);  say $it->();  say $it->();

will return 211 followed by 223, as those are the next primes E<gt>= 200.
On each call, the iterator returns the current value and increments to
the next prime.

Other options include L</forprimes> (more efficiency, less flexibility),
L<Math::Prime::Util::PrimeIterator> (an iterator with more functionality),
or L<Math::Prime::Util::PrimeArray> (a tied array).


=head2 prime_iterator_object

  my $it = prime_iterator_object;
  while ($it->value < 100) { say $it->value; $it->next; }
  $sum += $it->iterate for 1..100000;

Returns a L<Math::Prime::Util::PrimeIterator> object.  A shortcut that loads
the package if needed, calls new, and returns the object.  See the
documentation for that package for details.  This object has more features
than the simple one above (e.g. the iterator is bi-directional), and also
handles iterating across bigints.


=head2 prime_count

  my $primepi = prime_count( 1_000 );
  my $pirange = prime_count( 1_000, 10_000 );

Returns the Prime Count function C<Pi(n)>, also called C<primepi> in some
math packages.  When given two arguments, it returns the inclusive
count of primes between the ranges.  E.g. C<(13,17)> returns 2, C<(14,17)>
and C<(13,16)> return 1, C<(14,16)> returns 0.

The current implementation decides based on the ranges whether to use a
segmented sieve with a fast bit count, or the extended LMO algorithm.
The former is preferred for small sizes as well as small ranges.
The latter is much faster for large ranges.

The segmented sieve is very memory efficient and is quite fast even with
large base values.  Its complexity is approximately C<O(sqrt(a) + (b-a))>,
where the first term is typically negligible below C<~ 10^11>.  Memory use
is proportional only to C<sqrt(a)>, with total memory use under 1MB for any
base under C<10^14>.

The extended LMO method has complexity approximately
C<O(b^(2/3)) + O(a^(2/3))>, and also uses low memory.
A calculation of C<Pi(10^14)> completes in a few seconds, C<Pi(10^15)>
in well under a minute, and C<Pi(10^16)> in about one minute.  In
contrast, even parallel primesieve would take over a week on a
similar machine to determine C<Pi(10^16)>.

Also see the function L</prime_count_approx> which gives a very good
approximation to the prime count, and L</prime_count_lower> and
L</prime_count_upper> which give tight bounds to the actual prime count.
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
generate any primes.  For values under C<10^36> this uses the Riemann R
function, which is quite accurate: an error of less than C<0.0005%> is typical
for input values over C<2^32>, and decreases as the input gets larger.  If
L<Math::MPFR> is installed, the Riemann R function is used for all values, and
will be very fast.  If not, then values of C<10^36> and larger will use the
approximation C<li(x) - li(sqrt(x))/2>.  While not as accurate as the Riemann
R function, it still should have error less than C<0.00000000000000001%>.

A slightly faster but much less accurate answer can be obtained by averaging
the upper and lower bounds.


=head2 nth_prime

  say "The ten thousandth prime is ", nth_prime(10_000);

Returns the prime that lies in index C<n> in the array of prime numbers.  Put
another way, this returns the smallest C<p> such that C<Pi(p) E<gt>= n>.

For relatively small inputs (below 1 million or so), this does a sieve over
a range containing the nth prime, then counts up to the number.  This is fairly
efficient in time and memory.  For larger values, create a low-biased estimate
using the inverse logarithmic integral, use a fast prime count, then sieve in
the small difference.

While this method is thousands of times faster than generating primes, and
doesn't involve big tables of precomputed values, it still can take a fair
amount of time for large inputs.  Calculating the C<10^12th> prime takes
about 1 second, the C<10^13th> prime takes under 10 seconds, and the
C<10^14th> prime (3475385758524527) takes under one minute.  Think about
whether a bound or approximation would be acceptable, as they can be
computed analytically.

If the result is larger than a native integer size (32-bit or 64-bit), the
result will take a very long time.  A later version of
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


=head2 is_pseudoprime

Takes a positive number C<n> and a base C<a> as input, and returns 1 if
C<n> is a probable prime to base C<a>.  This is the simple Fermat primality
test.  Removing primes, given base 2 this produces the sequence
L<OEIS A001567|http://oeis.org/A001567>.

=head2 is_strong_pseudoprime

  my $maybe_prime = is_strong_pseudoprime($n, 2);
  my $probably_prime = is_strong_pseudoprime($n, 2, 3, 5, 7, 11, 13, 17);

Takes a positive number as input and one or more bases.  The bases must be
greater than C<1>.  Returns 1 if the input is a strong probable prime
to all of the bases, and 0 if not.

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

Even inputs other than 2 will always return 0 (composite).  While the
algorithm does run with even input, most sources define it only on odd input.
Returning composite for all non-2 even input makes the function match most
other implementations including L<Math::Primality>'s C<is_strong_pseudoprime>
function.

=head2 miller_rabin

An alias for C<is_strong_pseudoprime>.  This name is deprecated.

=head2 is_lucas_pseudoprime

Takes a positive number as input, and returns 1 if the input is a standard
Lucas probable prime using the Selfridge method of choosing D, P, and Q (some
sources call this a Lucas-Selfridge pseudoprime).
Removing primes, this produces the sequence
L<OEIS A217120|http://oeis.org/A217120>.

=head2 is_strong_lucas_pseudoprime

Takes a positive number as input, and returns 1 if the input is a strong
Lucas probable prime using the Selfridge method of choosing D, P, and Q (some
sources call this a strong Lucas-Selfridge pseudoprime).  This is one half
of the BPSW primality test (the Miller-Rabin strong pseudoprime test with
base 2 being the other half).  Removing primes, this produces the sequence
L<OEIS A217255|http://oeis.org/A217255>.

=head2 is_extra_strong_lucas_pseudoprime

Takes a positive number as input, and returns 1 if the input passes the extra
strong Lucas test (as defined in
L<Grantham 2000|http://www.ams.org/mathscinet-getitem?mr=1680879>).  This
test has more stringent conditions than the strong Lucas test, and produces
about 60% fewer pseudoprimes.  Performance is typically 20-30% I<faster>
than the strong Lucas test.

The parameters are selected using the
L<Baillie-OEIS method|http://oeis.org/A217719>
method: increment C<P> from C<3> until C<jacobi(D,n) = -1>.
Removing primes, this produces the sequence
L<OEIS A217719|http://oeis.org/A217719>.

=head2 is_almost_extra_strong_lucas_pseudoprime

This is similar to the L</is_extra_strong_lucas_pseudoprime> function, but
does not calculate C<U>, so is a little faster, but also weaker.
With the current implementations, there is little reason to prefer this unless
trying to reproduce specific results.  The extra-strong implementation has been
optimized to use similar features, removing most of the performance advantage.

An optional second argument (an integer between 1 and 256) indicates the
increment amount for C<P> parameter selection.  The default value of 1 yields
the parameter selection described in L</is_extra_strong_lucas_pseudoprime>,
creating a pseudoprime sequence which is a superset of the latter's
pseudoprime sequence L<OEIS A217719|http://oeis.org/A217719>.
A value of 2 yields the method used by
L<Pari|http://pari.math.u-bordeaux.fr/faq.html#primetest>.

Because the C<U = 0> condition is ignored, this produces about 5% more
pseudoprimes than the extra-strong Lucas test.  However this is still only
66% of the number produced by the strong Lucas-Selfridge test.  No BPSW
counterexamples have been found with any of the Lucas tests described.

=head2 is_frobenius_underwood_pseudoprime

Takes a positive number as input, and returns 1 if the input passes the minimal
lambda+2 test (see Underwood 2012 "Quadratic Compositeness Tests"), where
C<(L+2)^(n-1) = 5 + 2x mod (n, L^2 - Lx + 1)>.  The computational cost for this
is between the cost of 2 and 3 strong pseudoprime tests.  There are no known
counterexamples, but this is not a well studied test.

=head2 miller_rabin_random

Takes a positive number (C<n>) as input and a positive number (C<k>) of bases
to use.  Performs C<k> Miller-Rabin tests using uniform random bases
between 2 and C<n-2>.

This should not be used in place of L</is_prob_prime>, L</is_prime>,
or L</is_provable_prime>.  Those functions will be faster and provide
better results than running C<k> Miller-Rabin tests.  This function can
be used if one wants more assurances for non-proven primes, such as for
cryptographic uses where the size is large enough that proven primes are
not desired.



=head2 is_prob_prime

  my $prob_prime = is_prob_prime($n);
  # Returns 0 (composite), 2 (prime), or 1 (probably prime)

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).

For 64-bit input (native or bignum), this uses either a deterministic set of
Miller-Rabin tests (1, 2, or 3 tests) or a strong BPSW test consisting of a
single base-2 strong probable prime test followed by a strong Lucas test.
This has been verified with Jan Feitsma's 2-PSP database to produce no false
results for 64-bit inputs.  Hence the result will always be 0 (composite) or
2 (prime).

For inputs larger than C<2^64>, an extra-strong Baillie-PSW primality test is
performed (also called BPSW or BSW).  This is a probabilistic test, so only
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
values as L</is_prime> and L</is_prob_prime>.  Note that numbers below 2^64
are considered proven by the deterministic set of Miller-Rabin bases or the
BPSW test.  Both of these have been tested for all small (64-bit) composites
and do not return false positives.

Using the L<Math::Prime::Util::GMP> module is B<highly recommended> for doing
primality proofs, as it is much, much faster.  The pure Perl code is just not
fast for this type of operation, nor does it have the best algorithms.
It should suffice for proofs of up to 40 digit primes, while the latest
MPU::GMP works for primes of hundreds of digits (thousands with an optional
larger polynomial set).

The pure Perl implementation uses theorem 5 of BLS75 (Brillhart, Lehmer, and
Selfridge's 1975 paper), an improvement on the Pocklington-Lehmer test.
This requires C<n-1> to be factored to C<(n/2)^(1/3))>.  This is often fast,
but as C<n> gets larger, it takes exponentially longer to find factors.

L<Math::Prime::Util::GMP> implements both the BLS75 theorem 5 test as well
as ECPP (elliptic curve primality proving).  It will typically try a quick
C<n-1> proof before using ECPP.  Certificates are available with either method.
This results in proofs of 200-digit primes in under 1 second on average, and
many hundreds of digits are possible.  This makes it significantly faster
than Pari 2.1.7's C<is_prime(n,1)> which is the default for L<Math::Pari>.


=head2 prime_certificate

  my $cert = prime_certificate($n);
  say verify_prime($cert) ? "proven prime" : "not prime";

Given a positive integer C<n> as input, returns a primality certificate
as a multi-line string.  If we could not prove C<n> prime, an empty
string is returned (C<n> may or may not be composite).
This may be examined or given to L</verify_prime> for verification.  The latter
function contains the description of the format.


=head2 is_provable_prime_with_cert

Given a positive integer as input, returns a two element array containing
the result of L</is_provable_prime>:
  0  definitely composite
  1  probably prime
  2  definitely prime
and a primality certificate like L</prime_certificate>.
The certificate will be an empty string if the first element is not 2.


=head2 verify_prime

  my $cert = prime_certificate($n);
  say verify_prime($cert) ? "proven prime" : "not prime";

Given a primality certificate, returns either 0 (not verified)
or 1 (verified).  Most computations are done using pure Perl with
Math::BigInt, so you probably want to install and use Math::BigInt::GMP,
and ECPP certificates will be faster with Math::Prime::Util::GMP for
its elliptic curve computations.

If the certificate is malformed, the routine will carp a warning in addition
to returning 0.  If the C<verbose> option is set (see L</prime_set_config>)
then if the validation fails, the reason for the failure is printed in
addition to returning 0.  If the C<verbose> option is set to 2 or higher, then
a message indicating success and the certificate type is also printed.

A certificate may have arbitrary text before the beginning (the primality
routines from this module will not have any extra text, but this way
verbose output from the prover can be safely stored in a certificate).
The certificate begins with the line:

  [MPU - Primality Certificate]

All lines in the certificate beginning with C<#> are treated as comments
and ignored, as are blank lines.  A version number may follow, such as:

  Version 1.0

For all inputs, base 10 is the default, but at any point this may be
changed with a line like:

  Base 16

where allowed bases are 10, 16, and 62.  This module will only use base 10,
so its routines will not output Base commands.

Next, we look for (using "100003" as an example):

  Proof for:
  N 100003

where the text C<Proof for:> indicates we will read an C<N> value.  Skipping
comments and blank lines, the next line should be "N " followed by the number.

After this, we read one or more blocks.  Each block is a proof of the form:

  If Q is prime, then N is prime.

Some of the blocks have more than one Q value associated with them, but most
only have one.  Each block has its own set of conditions which must be
verified, and this can be done completely self-contained.  That is, each
block is independent of the other blocks and may be processed in any order.
To be a complete proof, each block must successfully verify.  The block
types and their conditions are shown below.

Finally, when all blocks have been read and verified, we must ensure we
can construct a proof tree from the set of blocks.  The root of the tree
is the initial C<N>, and for each node (block), all C<Q> values must
either have a block using that value as its C<N> or C<Q> must be less
than C<2^64> and pass BPSW.

Some other certificate formats (e.g. Primo) use an ordered chain, where
the first block must be for the initial C<N>, a single C<Q> is given which
is the implied C<N> for the next block, and so on.  This simplifies
validation implementation somewhat, and removes some redundant
information from the certificate, but has no obvious way to add proof
types such as Lucas or the various BLS75 theorems that use multiple
factors.  I decided that the most general solution was to have the
certificate contain the set in any order, and let the verifier do the
work of constructing the tree.

The blocks begin with the text "Type ..." where ... is the type.  One or
more values follow.  The defined types are:

=over 4

=item C<Small>

  Type Small
  N 5791

N must be less than 2^64 and be prime (use BPSW or deterministic M-R).

=item C<BLS3>

  Type BLS3
  N  2297612322987260054928384863
  Q  16501461106821092981
  A  5

A simple n-1 style proof using BLS75 theorem 3.  This block verifies if:
  a  Q is odd
  b  Q > 2
  c  Q divides N-1
  .  Let M = (N-1)/Q
  d  MQ+1 = N
  e  M > 0
  f  2Q+1 > sqrt(N)
  g  A^((N-1)/2) mod N = N-1
  h  A^(M/2) mod N != N-1

=item C<Pocklington>

  Type Pocklington
  N  2297612322987260054928384863
  Q  16501461106821092981
  A  5

A simple n-1 style proof using generalized Pocklington.  This is more
restrictive than BLS3 and much more than BLS5.  This is Primo's type 1,
and this module does not currently generate these blocks.
This block verifies if:
  a  Q divides N-1
  .  Let M = (N-1)/Q
  b  M > 0
  c  M < Q
  d  MQ+1 = N
  e  A > 1
  f  A^(N-1) mod N = 1
  g  gcd(A^M - 1, N) = 1

=item C<BLS15>

  Type BLS15
  N  8087094497428743437627091507362881
  Q  175806402118016161687545467551367
  LP 1
  LQ 22

A simple n+1 style proof using BLS75 theorem 15.  This block verifies if:
  a  Q is odd
  b  Q > 2
  c  Q divides N+1
  .  Let M = (N+1)/Q
  d  MQ-1 = N
  e  M > 0
  f  2Q-1 > sqrt(N)
  .  Let D = LP*LP - 4*LQ
  g  D != 0
  h  Jacobi(D,N) = -1
  .  Note: V_{k} indicates the Lucas V sequence with LP,LQ
  i  V_{m/2} mod N != 0
  j  V_{(N+1)/2} mod N == 0

=item C<BLS5>

  Type BLS5
  N  8087094497428743437627091507362881
  Q[1]  98277749
  Q[2]  3631
  A[0]  11
  ----

A more sophisticated n-1 proof using BLS theorem 5.  This requires N-1 to
be factored only to C<(N/2)^(1/3)>.  While this looks much more complicated,
it really isn't much more work.  The biggest drawback is just that we have
multiple Q values to chain rather than a single one.  This block verifies if:

  a  N > 2
  b  N is odd
  .  Note: the block terminates on the first line starting with a C<->.
  .  Let Q[0] = 2
  .  Let A[i] = 2 if Q[i] exists and A[i] does not
  c  For each i (0 .. maxi):
  c1   Q[i] > 1
  c2   Q[i] < N-1
  c3   A[i] > 1
  c4   A[i] < N
  c5   Q[i] divides N-1
  . Let F = N-1 divided by each Q[i] as many times as evenly possible
  . Let R = (N-1)/F
  d  F is even
  e  gcd(F, R) = 1
  . Let s = integer    part of R / 2F
  . Let f = fractional part of R / 2F
  . Let P = (F+1) * (2*F*F + (r-1)*F + 1)
  f  n < P
  g  s = 0  OR  r^2-8s is not a perfect square
  h  For each i (0 .. maxi):
  h1   A[i]^(N-1) mod N = 1
  h2   gcd(A[i]^((N-1)/Q[i])-1, N) = 1

=item C<ECPP>

  Type ECPP
  N  175806402118016161687545467551367
  A  96642115784172626892568853507766
  B  111378324928567743759166231879523
  M  175806402118016177622955224562171
  Q  2297612322987260054928384863
  X  3273750212
  Y  82061726986387565872737368000504

An elliptic curve primality block, typically generated with an Atkin/Morain
ECPP implementation, but this should be adequate for anything using the
Atkin-Goldwasser-Kilian-Morain style certificates.
Some basic elliptic curve math is needed for these.
This block verifies if:

  .  Note: A and B are allowed to be negative, with -1 not uncommon.
  .  Let A = A % N
  .  Let B = B % N
  a  N > 0
  b  gcd(N, 6) = 1
  c  gcd(4*A^3 + 27*B^2, N) = 1
  d  Y^2 mod N = X^3 + A*X + B mod N
  e  M >= N - 2*sqrt(N) + 1
  f  M <= N + 2*sqrt(N) + 1
  g  Q > (N^(1/4)+1)^2
  h  Q < N
  i  M != Q
  j  Q divides M
  .  Note: EC(A,B,N,X,Y) is the point (X,Y) on Y^2 = X^3 + A*X + B, mod N
  .        All values work in affine coordinates, but in theory other
  .        representations work just as well.
  .  Let POINT1 = (M/Q) * EC(A,B,N,X,Y)
  .  Let POINT2 = M * EC(A,B,N,X,Y)  [ = Q * POINT1 ]
  k  POINT1 is not the identity
  l  POINT2 is the identity

=back

=head2 is_aks_prime

  say "$n is definitely prime" if is_aks_prime($n);

Takes a positive number as input, and returns 1 if the input passes the
Agrawal-Kayal-Saxena (AKS) primality test.  This is a deterministic
unconditional primality test which runs in polynomial time for general input.

This function is only included for completeness and as an example.  The Perl
implementation is fast compared to the only other Perl implementation
available (in L<Math::Primality>), and the implementation in
L<Math::Prime::Util::GMP> compares favorably to others in the literature.
However AKS in general is far too slow to be of practical use.  R.P. Brent,
2010: "AKS is not a practical algorithm.  ECPP is much faster."


=head2 lucas_sequence

  my($U, $V, $Qk) = lucas_sequence($n, $P, $Q, $k)

Computes C<U_k>, C<V_k>, and C<Q_k> for the Lucas sequence defined by
C<P>,C<Q>, modulo C<n>.  The modular Lucas sequence is used in a
number of primality tests and proofs.

The following conditions must hold:
  - C<< D = P*P - 4*Q  !=  0 >>
  - C<< P > 0 >>
  - C<< P < n >>
  - C<< Q < n >>
  - C<< k >= 0 >>
  - C<< n >= 2 >>


=head2 moebius

  say "$n is square free" if moebius($n) != 0;
  $sum += moebius($_) for (1..200); say "Mertens(200) = $sum";

Returns μ(n), the Möbius function (also called the Moebius, Mobius, or
MoebiusMu function) for a non-negative integer input.  This function is 1 if
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

Returns M(n), the Mertens function for a non-negative integer input.  This
function is defined as C<sum(moebius(1..n))>, but calculated more efficiently
for large inputs.  For example, computing Mertens(100M) takes:

   time    approx mem
     0.4s      0.1MB   mertens(100_000_000)
    74.8s   7000MB     List::Util::sum(moebius(1,100_000_000))
    88.5s      0MB     $sum += moebius($_) for 1..100_000_000   [-nobigint]
   181.8s      0MB     $sum += moebius($_) for 1..100_000_000

The summation of individual terms via factoring is quite expensive in time,
though uses O(1) space.  This function will generate the equivalent output
via a sieving method, which will use some more memory, but be much faster.
The current method is a simple C<n^1/2> version of Deléglise and Rivat (1996),
which involves calculating all moebius values to C<n^1/2>, which in turn will
require prime sieving to C<n^1/4>.

Various algorithms exist for this, using differing quantities of μ(n).  The
simplest way is to efficiently sum all C<n> values.  Benito and Varona (2008)
show a clever and simple method that only requires C<n/3> values.  Deléglise
and Rivat (1996) describe a segmented method using only C<n^1/3> values.  The
current implementation does a simple non-segmented C<n^1/2> version of their
method.  Kuznetsov (2011) gives an alternate method that he indicates is even
faster.  Lastly, one of the advanced prime count algorithms could be
theoretically used to create a faster solution.


=head2 euler_phi

  say "The Euler totient of $n is ", euler_phi($n);

Returns φ(n), the Euler totient function (also called Euler's phi or phi
function) for an integer value.  This is an arithmetic function that counts
the number of positive integers less than or equal to C<n> that are relatively
prime to C<n>.  Given the definition used, C<euler_phi> will return 0 for all
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

Returns EXP(Λ(n)), the exponential of the Mangoldt function (also known
as von Mangoldt's function) for an integer value.
It is equal to log p if n is prime or a power of a prime,
and 0 otherwise.  We return the exponential so all results are integers.
Hence the return value for C<exp_mangoldt> is:

   p   if n = p^m for some prime p and integer m >= 1
   1   otherwise.


=head2 liouville

Returns λ(n), the Liouville function for a non-negative integer input.
This is -1 raised to Ω(n) (the total number of prime factors).


=head2 chebyshev_theta

  say chebyshev_theta(10000);

Returns θ(n), the first Chebyshev function for a non-negative integer input.
This is the sum of the logarithm of each prime where C<p E<lt>= n>.  An
alternate computation is as the logarithm of n primorial.
Hence these functions:

  use List::Util qw/sum/;  use Math::BigFloat;

  sub c1a { 0+sum( map { log($_) } @{primes(shift)} ) }
  sub c1b { Math::BigFloat->new(primorial(shift))->blog }

yield similar results, albeit slower and using more memory.


=head2 chebyshev_psi

  say chebyshev_psi(10000);

Returns ψ(n), the second Chebyshev function for a non-negative integer input.
This is the sum of the logarithm of each prime where C<p^k E<lt>= n> for an
integer k.  An alternate computation is as the summatory Mangoldt function.
Another alternate computation is as the logarithm of LCM(1,2,...,n).
Hence these functions:

  use List::Util qw/sum/;  use Math::BigFloat;

  sub c2a { 0+sum( map { log(exp_mangoldt($_)) } 1 .. shift ) }
  sub c2b { Math::BigFloat->new(consecutive_integer_lcm(shift))->blog }

yield similar results, albeit slower and using more memory.


=head2 divisor_sum

  say "Sum of divisors of $n:", divisor_sum( $n );

This function takes a positive integer as input and returns the sum of the
k-th powers of the divisors of the input, including 1 and itself.  If the
second argument (C<k>) is omitted it is assumed to be 1.  This is known as
the sigma function (see Hardy and Wright section 16.7, or OEIS A000203).
The API is identical to Pari/GP's C<sigma> function.

The second argument can be a code reference, which is called for each divisor
and the results are summed.  This allows computation of other functions,
but will be less efficient than using the numeric second argument.

An example of the 5th Jordan totient (OEIS A059378):

  divisor_sum( $n, sub { my $d=shift; $d**5 * moebius($n/$d); } );

though we have a function L</jordan_totient> which is more efficient.

This function is useful for calculating things like aliquot sums, abundant
numbers, perfect numbers, etc.

For numeric second arguments (sigma computations), the result will be a bigint
if necessary.  For the code reference case, the user must take care to return
bigints if overflow will be a concern.


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
match that terminology.  This function should return the same result as the
C<mpz_primorial_ui> function added in GMP 5.1.


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


=head2 consecutive_integer_lcm

  $lcm = consecutive_integer_lcm($n);

Given an unsigned integer argument, returns the least common multiple of all
integers from 1 to C<n>.  This can be done by manipulation of the primes up
to C<n>, resulting in much faster and memory-friendly results than using
a factorial.


=head2 partitions

Calculates the partition function p(n) for a non-negative integer input.
This is the number of ways of writing the integer n as a sum of positive
integers, without restrictions.  This corresponds to Pari's C<numbpart>
function and Mathematica's C<PartitionsP> function.  The values produced
in order are L<OEIS series A000041|http://oeis.org/A000041>.

This uses a combinatorial calculation, which means it will not be very
fast compared to Pari, Mathematica, or FLINT which use the Rademacher
formula using multi-precision floating point.  In 10 seconds, the pure
Perl version can produce C<partitions(10000)> while with
L<Math::Prime::Util::GMP> it can do C<partitions(200000)>.  In contrast,
in about 10 seconds Pari can solve C<numbpart(22000000)>.

If you want the enumerated partitions, see L<Integer::Partition>.  It is
very fast and uses an extremely memory efficient iterator.  It is not,
however, practical for producing the partition I<number> for values
over 100 or so.


=head2 carmichael_lambda

Returns the Carmichael function (also called the reduced totient function,
or Carmichael λ(n)) of a positive integer argument.  It is the smallest
positive integer m such that a^m = 1 mod n for every integer a coprime to n.
This is L<OEIS series A002322|http://oeis.org/A002322>.


=head2 znorder

  $order = znorder(2, next_prime(10**19)-6);

Given two positive integers C<a> and C<n>, returns the multiplicative order
of C<a> modulo <n>.  This is the smallest positive integer C<k> such that
C<a^k ≡ 1 mod n>.  Returns 1 if C<a = 1>.  Return undef if C<a = 0> or if
C<a> and C<n> are not coprime, since no value will result in 1 mod n.


=head2 random_prime

  my $small_prime = random_prime(1000);      # random prime <= limit
  my $rand_prime = random_prime(100, 10000); # random prime within a range

Returns a pseudo-randomly selected prime that will be greater than or equal
to the lower limit and less than or equal to the upper limit.  If no lower
limit is given, 2 is implied.  Returns undef if no primes exist within the
range.

The goal is to return a uniform distribution of the primes in the range,
meaning for each prime in the range, the chances are equally likely that it
will be seen.  This is removes from consideration such algorithms as
C<PRIMEINC>, which although efficient, gives very non-random output.  This
also implies that the numbers will not be evenly distributed, since the
primes are not evenly distributed.  Stated again, the random prime functions
return a uniformly selected prime from the set of primes within the range.
Hence given C<random_prime(1000)>, the numbers 2, 3, 487, 631, and 997 all
have the same probability of being returned.

For small numbers, a random index selection is done, which gives ideal
uniformity and is very efficient with small inputs.  For ranges larger than
this ~16-bit threshold but within the native bit size, a Monte Carlo method
is used (multiple calls to C<irand> will be made if necessary).  This also
gives ideal uniformity and can be very fast for reasonably sized ranges.
For even larger numbers, we partition the range, choose a random partition,
then select a random prime from the partition.  This gives some loss of
uniformity but results in many fewer bits of randomness being consumed as
well as being much faster.

If an C<irand> function has been set via L</prime_set_config>, it will be
used to construct any ranged random numbers needed.  The function should
return a uniformly random 32-bit integer, which is how the irand functions
exported by L<Math::Random::Secure>, L<Math::Random::MT>,
L<Math::Random::ISAAC>, and most other modules behave.

If no C<irand> function was set, then L<Bytes::Random::Secure> is used with
a non-blocking seed.  This will create good quality random numbers, so there
should be little reason to change unless one is generating long-term keys,
where using the blocking random source may be preferred.

Examples of various ways to set your own irand function:

  # Math::Random::Secure.  Uses ISAAC and strong seed methods.
  use Math::Random::Secure;
  prime_set_config(irand => \&Math::Random::Secure::irand);

  # Bytes::Random::Secure (OO interface with full control of options):
  use Bytes::Random::Secure ();
  BEGIN {
    my $rng = Bytes::Random::Secure->new( Bits => 512 );
    sub irand { return $rng->irand; }
  }
  prime_set_config(irand => \&irand);

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
then the result will be returned as a BigInt.  However, if the C<-nobigint>
tag was used, then numbers larger than the threshold will be flagged as an
error, and numbers on the threshold will be restricted to native numbers.
For better performance with large bit sizes, install L<Math::Prime::Util::GMP>.


=head2 random_nbit_prime

  my $bigprime = random_nbit_prime(512);

Selects a random n-bit prime, where the input is an integer number of bits
between 2 and the maximum representable bits (32, 64, or 100000 for native
32-bit, native 64-bit, and bigint respectively).  A prime with the nth bit
set will be uniformly selected, with randomness supplied via calls to the
C<irand> function as described above.

For bit sizes of 64 and lower, L</random_prime> is used, which gives completely
uniform results in this range.  For sizes larger than 64, Algorithm 1 of
Fouque and Tibouchi (2011) is used, wherein we select a random odd number
for the lower bits, then loop selecting random upper bits until the result
is prime.  This allows a more uniform distribution than the general
L</random_prime> case while running slightly faster (in contrast, for large
bit sizes L</random_prime> selects a random upper partition then loops
on the values within the partition, which very slightly skews the results
towards smaller numbers).

The C<irand> function is used for randomness, so all the discussion in
L</random_prime> about that applies here.
The result will be a BigInt if the number of bits is greater than the native
bit size.  For better performance with large bit sizes, install
L<Math::Prime::Util::GMP>.


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
primes offer security equivalent to strong primes.

Similar to L</random_nbit_prime>, the result will be a BigInt if the
number of bits is greater than the native bit size.  For better performance
with large bit sizes, install L<Math::Prime::Util::GMP>.


=head2 random_proven_prime

  my $bigprime = random_proven_prime(512);

Constructs an n-bit random proven prime.  Internally this may use
L</is_provable_prime>(L</random_nbit_prime>) or
L</random_maurer_prime> depending on the platform and bit size.


=head2 random_proven_prime_with_cert

  my($n, $cert) = random_proven_prime_with_cert(512)

Similar to L</random_proven_prime>, but returns a two-element array containing
the n-bit provable prime along with a primality certificate.  The certificate
is the same as produced by L</prime_certificate> or
L</is_provable_prime_with_cert>, and can be parsed by L</verify_prime> or
any other software that understands MPU primality certificates.


=head2 random_maurer_prime

  my $bigprime = random_maurer_prime(512);

Construct an n-bit provable prime, using the FastPrime algorithm of
Ueli Maurer (1995).  This is the same algorithm used by L<Crypt::Primes>.
Similar to L</random_nbit_prime>, the result will be a BigInt if the
number of bits is greater than the native bit size.  For better performance
with large bit sizes, install L<Math::Prime::Util::GMP>.

The differences between this function and that in L<Crypt::Primes> are
described in the L</"SEE ALSO"> section.

Internally this additionally runs the BPSW probable prime test on every
partial result, and constructs a primality certificate for the final
result, which is verified.  These provide additional checks that the resulting
value has been properly constructed.

An alternative to this function is to run L</is_provable_prime> on the
result of L</random_nbit_prime>, which will provide more diversity and
will be faster up to 512 or so bits.  Maurer's method should be much
faster for large bit sizes (larger than 2048).  If you don't need absolutely
proven results, then using L</random_nbit_prime> followed by additional
tests (L</is_strong_pseudoprime> and/or L</is_frobenius_underwood_pseudoprime>)
should be much faster.


=head2 random_maurer_prime_with_cert

  my($n, $cert) = random_maurer_prime_with_cert(512)

As with L</random_maurer_prime>, but returns a two-element array containing
the n-bit provable prime along with a primality certificate.  The certificate
is the same as produced by L</prime_certificate> or
L</is_provable_prime_with_cert>, and can be parsed by L</verify_prime> or
any other software that understands MPU primality certificates.
The proof construction consists of a single chain of C<BLS3> types.



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

  xs              Allows turning off the XS code, forcing the Pure Perl
                  code to be used.  Set to 0 to disable XS, set to 1 to
                  re-enable.  You probably will never want to do this.

  gmp             Allows turning off the use of L<Math::Prime::Util::GMP>,
                  which means using Pure Perl code for big numbers.  Set
                  to 0 to disable GMP, set to 1 to re-enable.
                  You probably will never want to do this.

  assume_rh       Allows functions to assume the Riemann hypothesis is
                  true if set to 1.  This defaults to 0.  Currently this
                  setting only impacts prime count lower and upper
                  bounds, but could later be applied to other areas such
                  as primality testing.  A later version may also have a
                  way to indicate whether no RH, RH, GRH, or ERH is to
                  be assumed.

  irand           Takes a code ref to an irand function returning a
                  uniform number between 0 and 2**32-1.  This will be
                  used for all random number generation in the module.


=head1 FACTORING FUNCTIONS

=head2 factor

  my @factors = factor(3_369_738_766_071_892_021);
  # returns (204518747,16476429743)

Produces the prime factors of a positive number input, in numerical order.
The special cases of C<n = 0> and C<n = 1> will return C<n>, which
guarantees multiplying the factors together will always result in the
input value, though those are the only cases where the returned factors
are not prime.

In scalar context, returns Ω(n), the total number of prime factors
(L<OEIS A001222|http://oeis.org/A001222>).
This corresponds to Pari's C<bigomega(n)> function and Mathematica's
C<PrimeOmega[n]> function.
This is same result that we would get if we evaluated the resulting
array in scalar context.
Do note that the inputs of C<0> and C<1> will return C<1>, contrary
to the standard definition of Omega.

The current algorithm for non-bigints is a sequence of small trial division,
a few rounds of Pollard's Rho, SQUFOF, Pollard's p-1, Hart's OLF, a long
run of Pollard's Rho, and finally trial division if anything survives.  This
process is repeated for each non-prime factor.  In practice, it is very rare
to require more than the first Rho + SQUFOF to find a factor, and I have not
seen anything go to the last step.

Factoring bigints works with pure Perl, and can be very handy on 32-bit
machines for numbers just over the 32-bit limit, but it can be B<very> slow
for "hard" numbers.  Installing the L<Math::Prime::Util::GMP> module will speed
up bigint factoring a B<lot>, and all future effort on large number factoring
will be in that module.  If you do not have that module for some reason, use
the GMP or Pari version of bigint if possible
(e.g. C<use bigint try =E<gt> 'GMP,Pari'>), which will run 2-3x faster (though
still 100x slower than the real GMP code).


=head2 factor_exp

  my @factor_exponent_pairs = factor_exp(29513484000);
  # returns ([2,5], [3,4], [5,3], [7,2], [11,1], [13,2])
  # factor(29513484000)
  # returns (2,2,2,2,2,3,3,3,3,5,5,5,7,7,11,13,13)

Produces pairs of prime factors and exponents in numerical factor order.
This may be more convenient for some algorithms.  This is the same form
that Mathematica's C<FactorInteger[n]> function returns.

In scalar context, returns ω(n), the number of unique prime factors
(L<OEIS A001221|http://oeis.org/A001221>).
This corresponds to Pari's C<omega(n)> function and Mathematica's
C<PrimeNu[n]> function.
This is same result that we would get if we evaluated the resulting
array in scalar context.
Do note that the inputs of C<0> and C<1> will return C<1>, contrary
to the standard definition of omega.

The internals are identical to L</factor>, so all comments there apply.
Just the way the factors are arranged is different.


=head2 all_factors

  my @divisors = all_factors(30);   # returns (1, 2, 3, 5, 6, 10, 15, 30)

Produces all the divisors of a positive number input, including 1 and
the input number.  The divisors are a power set of multiplications of
the prime factors, returned as a uniqued sorted list.  The result is
identical to that of Pari's C<divisors> function.

In scalar context this returns the sigma0 function,
the sigma function (see Hardy and Wright section 16.7, or OEIS A000203).


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

=head2 pplus1_factor

  my @factors = pplus1_factor($n);
  my @factors = pplus1_factor($n, 1_000);          # set B1 smoothness

Produces factors, not necessarily prime, of the positive number input.  This
is Williams' C<p+1> method, using one stage and two predefined initial points.



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

For BigInt / BigFloat objects, we first check to see if L<Math::MPFR> is
available.  If so, then it is used since it is very fast and has high accuracy.
Accuracy when using MPFR will be equal to the C<accuracy()> value of the
input (or the default BigFloat accuracy, which is 40 by default).

MPFR is used for positive inputs only.  If L<Math::MPFR> is not available
or the input is negative, then other methods are used:
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

For BigInt / BigFloat objects, we first check to see if L<Math::MPFR> is
available.  If so, then it is used, as it will return results much faster
and can be more accurate.  Accuracy when using MPFR will be equal to the
C<accuracy()> value of the input (or the default BigFloat accuracy, which
is 40 by default).

MPFR is used for inputs greater than 1 only.  If L<Math::MPFR> is not
installed or the input is less than 1, results will be calculated as
C<Ei(ln x)>.


=head2 RiemannZeta

  my $z = RiemannZeta($s);

Given a floating point input C<s> where C<s E<gt>= 0>, returns the floating
point value of ζ(s)-1, where ζ(s) is the Riemann zeta function.  One is
subtracted to ensure maximum precision for large values of C<s>.  The zeta
function is the sum from k=1 to infinity of C<1 / k^s>.  This function only
uses real arguments, so is basically the Euler Zeta function.

If the bignum module has been loaded, all inputs will be treated as if they
were Math::BigFloat objects.

For non-BigInt/BigFloat objects, the result should be accurate to at least 14
digits.  The XS code uses a rational Chebyshev approximation between 0.5 and 5,
and a series for other values.  The PP code uses an identical series for all
values.

For BigInt / BigFloat objects, we first check to see if the Math::MPFR module
is installed.  If so, then it is used, as it will return results much faster
and can be more accurate.  Accuracy when using MPFR will be equal to the
C<accuracy()> value of the input (or the default BigFloat accuracy, which
is 40 by default).

If Math::MPFR is not installed, then results are calculated using either
Borwein (1991) algorithm 2, or the basic series.  Full input accuracy is
attempted, but Math::BigFloat
L<RT 43692|https://rt.cpan.org/Ticket/Display.html?id=43692>
produces incorrect high-accuracy computations without the fix.
It is also very slow.  I highly
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

Print strong pseudoprimes to base 17 up to 10M:

    # Similar to A001262's isStrongPsp function, but over 4x faster
    perl -MMath::Prime::Util=:all -E 'my $n=3; while($n <= 10000000) { print "$n " if is_strong_pseudoprime($n,$base) && !is_prime($n); $n+=2; } BEGIN {$|=1; $base=17}'

or, slightly faster, use forprimes and loop over the odds between primes:

   perl -MMath::Prime::Util=:all -E '$|=1; $base=17; my $prev = 1; forprimes { $prev += 2; while ($prev < $_) { print "$prev " if is_strong_pseudoprime($prev,$base); $prev += 2; } } 3,10000000'

Print some primes above 64-bit range:

    perl -MMath::Prime::Util=:all -Mbigint -E 'my $start=100000000000000000000; say join "\n", @{primes($start,$start+1000)}'
    # Similar code using Pari:
    # perl -MMath::Pari=:int,PARI,nextprime -E 'my $start = PARI "100000000000000000000"; my $end = $start+1000; my $p=nextprime($start); while ($p <= $end) { say $p; $p = nextprime($p+1); }'

Examining the η3(x) function of Planat and Solé (2011):

  sub nu3 {
    my $n = shift;
    my $phix = chebyshev_psi($n);
    my $nu3 = 0;
    foreach my $nu (1..3) {
      $nu3 += (moebius($nu)/$nu)*LogarithmicIntegral($phix**(1/$nu));
    }
    return $nu3;
  }
  say prime_count(1000000);
  say prime_count_approx(1000000);
  say nu3(1000000);

Construct and use a Sophie-Germain prime iterator:

  sub make_sophie_germain_iterator {
    my $p = shift || 2;
    my $it = prime_iterator($p);
    return sub {
      do { $p = $it->() } while !is_prime(2*$p+1);
      $p;
    };
  }
  my $sgit = make_sophie_germain_iterator();
  for (1 .. 10000) {
    print $sgit->(), "\n";
  }

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

  perl -MMath::Prime::Util=:all -E 'my @p = grep { /1/&&/2/&&/3/&&/4/&&/5/&&/6/&&/7/} @{primes(1000000,9999999)}; say $p[-1];'

Project Euler, problem 47 (Distinct primes factors):

  use Math::Prime::Util qw/pn_primorial factor/;
  use List::MoreUtils qw/distinct/;
  sub nfactors { scalar distinct factor(shift); }
  my $n = pn_primorial(4);  # Start with the first 4-factor number
  $n++ while (nfactors($n) != 4 || nfactors($n+1) != 4 || nfactors($n+2) != 4 || nfactors($n+3) != 4);
  say $n;

Project Euler, problem 69, stupid brute force solution (about 2 seconds):

  use Math::Prime::Util qw/euler_phi/;
  my ($n, $max) = (0,0);
  do {
    my $ndivphi = $_ / euler_phi($_);
    ($n, $max) = ($_, $ndivphi) if $ndivphi > $max;
  } for 1..1000000;
  say "$n  $max";

Here is the right way to do PE problem 69 (under 0.03s):

  use Math::Prime::Util qw/pn_primorial/;
  my $n = 0;
  $n++ while pn_primorial($n+1) < 1000000;
  say pn_primorial($n);

Project Euler, problem 187, stupid brute force solution, ~4 minutes:

  use Math::Prime::Util qw/factor -nobigint/;
  my $nsemis = 0;
  do { $nsemis++ if scalar factor($_) == 2; }
     for 1 .. int(10**8)-1;
  say $nsemis;

Here is the best way for PE187.  Under 30 milliseconds from the command line:

  use Math::Prime::Util qw/primes prime_count/;
  use List::Util qw/sum/;
  my $limit = shift || int(10**8);
  my @primes = @{primes(int(sqrt($limit)))};
  say sum( map { prime_count(int(($limit-1)/$primes[$_-1])) - $_ + 1 }
               1 .. scalar @primes );


=head1 PRIMALITY TESTING NOTES

Above C<2^64>, L</is_prob_prime> performs an extra-strong
L<BPSW test|http://en.wikipedia.org/wiki/Baillie-PSW_primality_test>
which is fast (a little less than the time to perform 3 Miller-Rabin
tests) and has no known counterexamples.  If you believe the primality
testing done by Pari, Maple, SAGE, FLINT, etc., then this function
should be appropriate for you.  L</is_prime> will do the same BPSW
test as well as some additional testing, making it slightly more time
consuming but less likely to produce a false result.  This is a little
more stringent than Mathematica.  L</is_provable_prime> constructs a
primality proof.  If a certificate is requested, then either BLS75
theorem 5 or ECPP is performed.  Without a certificate, the method
is implementation specific (currently it is identical, but later
releases may use APRCL).  With L<Math::Prime::Util::GMP> installed,
this is quite fast through 300 or so digits.

Math systems 30 years ago used to use Miller-Rabin tests with C<k>
bases (typically fixed bases, sometimes random) for primality
testing, but these have generally been replaced by some form of BPSW
as used in this module.  See Pinch's 1993 paper for examples of why
using C<k> M-R tests leads to poor results.  The three exceptions in
common contemporary use I am aware of are:

=over 4

=item libtommath

Uses the first C<k> prime bases.  This is unacceptable for
cryptographic use, as there are known methods (e.g. Arnault 1994) for
constructing counterexamples.  The number of bases required to avoid
false results is unreasonably high, hence performance is slow even
if one ignores counterexamples.  Unfortunately this is the
multi-precision math library used for Perl 6 and at least one CPAN
Crypto module.

=item GMP/MPIR

Uses a set of C<k> static-random bases.  The bases are randomly chosen
using a PRNG that is seeded identically each call (the seed changes
with each release).  This offers a very slight advantage over using
the first C<k> prime bases, but not much.  See, for example, Nicely's
L<mpz_probab_prime_p pseudoprimes|http://www.trnicely.net/misc/mpzspsp.html>
page.

=item L<Math::Pari>

Pari 2.1.7 is the default version installed with the L<Math::Pari>
module.  It uses 10 random M-R bases (the PRNG uses a fixed seed
set at compile time).  Pari 2.3.0 was released in May 2006 and it,
like all later releases through at least 2.6.1, use BPSW / APRCL,
after complaints of false results from using M-R tests.

=back

Basically the problem is that it is just too easy to get counterexamples
from running C<k> M-R tests, forcing one to use a very large number of
tests (at least 20) to avoid frequent false results.  Using the BPSW test
results in no known counterexamples after 30+ years and runs much faster.
It can be enhanced with one or more random bases if one desires, and
will I<still> be much faster.

Using C<k> fixed bases has another problem, which is that in any
adversarial situation we can assume the inputs will be selected such
that they are one of our counterexamples.  Now we need absurdly large
numbers of tests.  This is like playing "pick my number" but the
number is fixed forever at the start, the guesser gets to know
everyone else's guesses and results who has every played, and can
keep playing as long as they like.  It's only valid if the players
are completely oblivious to what is happening.


=head1 LIMITATIONS

Perl versions earlier than 5.8.0 have a rather broken 64-bit implementation,
in that the values are accessed as doubles.  Hence any value larger
than C<~ 2^49> will start losing bottom bits.  This causes numerous functions
to not work properly.  The test suite will try to determine if your Perl is
broken (this only applies to really old versions of Perl compiled for 64-bit
when using numbers larger than C<~ 2^49>).  The best solution is updating to
a more recent Perl.

The module is thread-safe and should allow good concurrency on all platforms
that support Perl threads except Win32.  With Win32, either don't use threads
or make sure C<prime_precalc> is called before using C<primes>,
C<prime_count>, or C<nth_prime> with large inputs.  This is B<only>
an issue if you use non-Cygwin Win32 B<and> call these routines from within
Perl threads.


=head1 SEE ALSO

This section describes other CPAN modules available that have some feature
overlap with this one.  Also see the L</REFERENCES> section.  Please let me
know if any of this information is inaccurate.  Also note that just because
a module doesn't match what I believe are the best set of features, doesn't
mean it isn't perfect for someone else.

I will use SoE to indicate the Sieve of Eratosthenes, and MPU to denote this
module (L<Math::Prime::Util>).  Some quick alternatives I can recommend if
you don't want to use MPU:

=over 4

=item * L<Math::Prime::FastSieve> is the alternative module I use for basic
functionality with small integers.  It's fast and simple, and has a good
set of features.

=item * L<Math::Primality> is the alternative module I use for primality
testing on bigints.  The downside is that it can be slow, and the functions
other than primality tests are I<very> slow.

=item * L<Math::Pari> if you want the kitchen sink and can install it and
handle using it.  There are still some functions it doesn't do well
(e.g. prime count and nth_prime).

=back


L<Math::Prime::XS> has C<is_prime> and C<primes> functionality.  There is
no bigint support.  The C<is_prime> function uses well-written trial
division, meaning it is very fast for small numbers, but terribly slow for
large 64-bit numbers.  Because MPU does input validation and bigint
conversion, there is about 20 microseconds of additional overhead making
MPXS a little faster for tiny inputs, but once over 700k MPU is faster
for all values.  MPXS's prime sieve is an unoptimized non-segmented SoE
which returns an array.  Sieve bases larger than C<10^7> start taking
inordinately long and using a lot of memory (gigabytes beyond C<10^10>).
E.g. C<primes(10**9, 10**9+1000)> takes 36 seconds with MPXS, but only
0.00015 seconds with MPU.

L<Math::Prime::FastSieve> supports C<primes>, C<is_prime>, C<next_prime>,
C<prev_prime>, C<prime_count>, and C<nth_prime>.  The caveat is that all
functions only work within the sieved range, so are limited to about C<10^10>.
It uses a fast SoE to generate the main sieve.  The sieve is 2-3x slower than
the base sieve for MPU, and is non-segmented so cannot be used for
larger values.  Since the functions work with the sieve, they are very fast.
The fast bit-vector-lookup functionality can be replicated in MPU using
C<prime_precalc> but is not required.

L<Bit::Vector> supports the C<primes> and C<prime_count> functionality in a
somewhat similar way to L<Math::Prime::FastSieve>.  It is the slowest of all
the XS sieves, and has the most memory use.  It is faster than pure Perl code.

L<Crypt::Primes> supports C<random_maurer_prime> functionality.  MPU has
more options for random primes (n-digit, n-bit, ranged, and strong) in
addition to Maurer's algorithm.  MPU does not have the critical bug
L<RT81858|https://rt.cpan.org/Ticket/Display.html?id=81858>.  MPU should have
a more uniform distribution as well as return a larger subset of primes
(L<RT81871|https://rt.cpan.org/Ticket/Display.html?id=81871>).
MPU does not depend on L<Math::Pari> though can run slow for bigints unless
the L<Math::BigInt::GMP> or L<Math::BigInt::Pari> modules are installed.
Having L<Math::Prime::Util::GMP> installed also helps performance for MPU.
Crypt::Primes is hardcoded to use L<Crypt::Random>, while MPU uses
L<Bytes::Random::Secure>, and also allows plugging in a random function.
This is more flexible, faster, has fewer dependencies, and uses a CSPRNG
for security.  MPU can return a primality certificate.
What Crypt::Primes has that MPU does not is support for returning a generator.

L<Math::Factor::XS> calculates prime factors and factors, which correspond to
the L</factor> and L</all_factors> functions of MPU.  These functions do
not support bigints.  Both are implemented with trial division, meaning they
are very fast for really small values, but quickly become unusably slow
(factoring 19 digit semiprimes is over 700 times slower).  It has additional
functions C<count_prime_factors> (possible in MPU using C<scalar factor($n)>)
and C<matches> which has no equivalent.

L<Math::Big> version 1.12 includes C<primes> functionality.  The current
code is only usable for very tiny inputs as it is incredibly slow and uses
lots of memory.  L<RT81986|https://rt.cpan.org/Ticket/Display.html?id=81986>
has a patch to make it run much faster and use much less memory.  Since it is
in pure Perl it will still run quite slow compared to MPU.

L<Math::Big::Factors> supports factorization using wheel factorization (smart
trial division).  It supports bigints.  Unfortunately it is extremely slow on
any input that isn't the product of just small factors.  Even 7 digit inputs
can take hundreds or thousands of times longer to factor than MPU or
L<Math::Factor::XS>.  19-digit semiprimes will take I<hours> versus MPU's
single milliseconds.

L<Math::Factoring> is a placeholder module for bigint factoring.  Version 0.02
only supports trial division (the Pollard-Rho method does not work).

L<Math::Prime::TiedArray> allows random access to a tied primes array, almost
identically to what MPU provides in L<Math::Prime::Util::PrimeArray>.  MPU
has attempted to fix Math::Prime::TiedArray's shift bug
(L<RT58151|https://rt.cpan.org/Ticket/Display.html?id=58151>).  MPU is
typically much faster and will use less memory, but there are some cases where
MP:TA is faster (MP:TA stores all entries up to the largest request, while
MPU:PA stores only a window around the last request).

L<Math::Primality> supports C<is_prime>, C<is_pseudoprime>,
C<is_strong_pseudoprime>, C<is_strong_lucas_pseudoprime>, C<next_prime>,
C<prev_prime>, C<prime_count>, and C<is_aks_prime> functionality.
This is a great little module that implements
primality functionality.  It was the first module to support the BPSW and
AKS tests.  All inputs are processed using GMP, so it of course supports
bigints.  In fact, Math::Primality was made originally with bigints in mind,
while MPU was originally targeted to native integers, but both have added
better support for the other.  The main differences are extra functionality
(MPU has more functions) and performance.  With native integer inputs, MPU
is generally much faster, especially with L</prime_count>.  For bigints,
MPU is slower unless the L<Math::Prime::Util::GMP> module is installed, in
which case MPU is ~2x faster.  L<Math::Primality> also installs
a C<primes.pl> program, but it has much less functionality than the one
included with MPU.

L<Math::NumSeq> is more a related module rather than one with direct
functionality.  It does however offer a way to get similar results such as
primes, twin primes, Sophie-Germain primes, lucky primes, moebius, divisor
count, factor count, Euler totient, primorials, etc.  Math::NumSeq is mainly
set up for accessing these values in order, rather than for arbitrary values,
though some sequences support that.  The primary advantage I see is the
uniform access mechanism for a I<lot> of sequences.  For those methods that
overlap, MPU is usually much faster.  Importantly, most of the sequences in
Math::NumSeq are limited to 32-bit indices.

L<Math::Pari> supports a lot of features, with a great deal of overlap.  In
general, MPU will be faster for native 64-bit integers, while it's differs
for bigints (Pari will always be faster if L<Math::Prime::Util::GMP> is not
installed; with it, it varies by function).
Some of the highlights:

=over 4

=item C<isprime>

Similar to MPU's L</is_prob_prime> or L</is_prime> functions.
MPU is deterministic for native integers, and uses a strong
BPSW test for bigints (with a quick primality proof tried as well).  The
default version of Pari used by L<Math::Pari> (2.1.7) uses 10 random M-R
bases, which is a probable prime test usually considered much weaker than the
BPSW test used by MPU and newer versions of Pari (though better than a fixed
set of bases).  Calling as C<isprime($n,1)> performs a Pocklington-Lehmer
C<n-1> proof.  This is comparable in performance to MPU:GMP's C<n-1> proof
implementation, and is reasonably fast for about 70 digits, but much slower
than ECPP.

If L<Math::Pari> is compiled with version 2.3.5 of Pari (this is not easy to
do on many platforms), then the algorithms are completely different.  The
C<isprime> function now acts like L</is_provable_prime> -- an APRCL proof
is performed, which is quite efficient though requires using a larger stack
for numbers of 300+ digits.  It is somewhat comparable in speed to MPU:GMP's
ECPP proof method, but without a certificate.  Using the C<ispseudoprime>
function will perform a BPSW test similar to L</is_prob_prime>.

=item C<primepi>

Only available with version 2.3 of Pari.  Similar to MPU's L</prime_count>
function in API, but uses a naive counting algorithm with its precalculated
primes, so is not of practical use.  Incidently, Pari 2.6 (not usable from
Perl) has fixed the pre-calculation requirement so it is more useful, but is
still hundreds of times slower than MPU.

=item C<primes>

Doesn't support ranges, requires bumping up the precalculated
primes for larger numbers, which means knowing in advance the upper limit
for primes.  Support for numbers larger than 400M requires using Pari
version 2.3.5.  If that is used, sieving is about 2x faster than MPU, but
doesn't support segmenting.

=item C<factorint>

Similar to MPU's L</factor> though with a different return.  MPU offers
L</factor> for a linear array of prime factors where
   n = p1 * p2 * p3 * ...   as (p1,p2,p3,...)
and L</factor_exp> for an array of factor/exponent pairs where:
   n = p1^e1 * p2^e2 * ...  as ([p1,e1],[p2,e2],...)
while Pari returns a 2D Pari matrix which may be interpreted as:
   n = p1^e1 * p2^e2 * ...  as ([p1,p2,...],[e1,e2,...])
Slower than MPU for all 64-bit inputs on an x86_64 platform, it may be
faster for large values on other platforms.  With the newer
L<Math::Prime::Util::GMP> releases, bigint factoring is slightly
faster on average in MPU.

=item C<divisors>

Similar to MPU's L</all_factors>.

=item C<eulerphi>

Similar to MPU's L</euler_phi>.  MPU is 2-20x faster for native integers.
There is also support for a range, which can be much more efficient.
Without L<Math::Prime::Util::GMP> installed, MPU is very slow with bigints.
With it installed, it is about 2x slower than Math::Pari.

=item C<moebius>

Similar to MPU's L</moebius>.  Comparisons are similar to C<eulerphi>.

=item C<sigma>

Similar to MPU's L</divisor_sum>.  MPU is ~10x faster for native integers
and about 2x slower for bigints.

=item C<numbpart>

Similar to MPU's L</partitions>.  This function is not in Pari 2.1, which
is the default version used by Math::Pari.  With Pari 2.3 or newer, the
functions produce identical results, but Pari is much, much faster.

=item C<eint1>

Similar to MPU's L</ExponentialIntegral>.

=item C<zeta>

A more feature-rich version MPU's L</RiemannZeta> function (supports negative
and complex inputs).

=back

Overall, L<Math::Pari> supports a huge variety of functionality and has a
sophisticated and mature code base behind it (noting that the default version
of Pari used is about 10 years old now).
For native integers often using Math::Pari will be slower, but bigints are
often superior and it rarely has any performance surprises.  Some of the
unique features MPU offers include super fast prime counts, nth_prime,
ECPP primality proofs with certificates, approximations and limits for both,
random primes, fast Mertens calculations, Chebyshev theta and psi functions,
and the logarithmic integral and Riemann R functions.  All with fairly
minimal installation requirements.


=head1 PERFORMANCE

First, for those looking for the state of the art non-Perl solutions:

=over 4

=item Primality testing

L<PFGW|http://sourceforge.net/projects/openpfgw/> is the fastest primality
testing software I'm aware of once past 2000 or so digits, has fast trial
division, and is especially fast on many special forms.  It does not have
a BPSW test however, and there are quite a few counterexamples for a given
base of its PRP test, so for primality testing it is most useful for fast
filtering of very large candidates.  A test such as the BPSW test in this
module is then recommended.

=item Primality proofs

L<Primo|http://www.ellipsa.eu/> is the best method for open source primality
proving for inputs over 1000 digits.  Primo also does well below that size,
but other good alternatives are
L<WraithX APRCL|http://sourceforge.net/projects/mpzaprcl/>,
the APRCL from the modern L<Pari|http://pari.math.u-bordeaux.fr/> package,
or the standalone ECPP from this module with large polynomial set.

=item Factoring

L<yafu|http://sourceforge.net/projects/yafu/>,
L<msieve|http://sourceforge.net/projects/msieve/>, and
L<gmp-ecm|http://ecm.gforge.inria.fr/> are all good choices for large
inputs.  The factoring code in this module (and all other CPAN modules) is
very limited compared to those.

=item Primes

L<primesieve|http://code.google.com/p/primesieve/> is the fastest publically
available code I am aware of.  It is much faster than any of the alternatives,
and even more so when run multi-threaded.  Tomás Oliveira e Silva's private
code may be faster for very large values, but it isn't available for testing.

Note that the Sieve of Atkin is I<not> faster than the Sieve of Eratosthenes
when both are well implemented.  The only Sieve of Atkin that is even
competitive is Bernstein's super optimized I<primegen>, which runs about
10% faster than the simple SoE in this module, slower than Pari and yafu's
SoE implementations, and 2x slower than primesieve.

=item Prime Counts and Nth Prime

Up to a limit, extensive use of tables plus a good segmented sieve will
produce the fastest results, but the number of tables needed to maintain
good performance grows exponentially.  The code in this module approaches
the best publically available results, with the notable exception of
Christian Bau's L-M-O implementation.  The author of primesieve is also
working on an L-M-O implementation.  None of these are state of the art
compared to private research methods.

=back


=head2 PRIME COUNTS

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
       0.002 Math::Prime::Util           0.35     using extended LMO
       0.007 Math::Prime::Util           0.12     using Lehmer's method
       0.27  Math::Prime::Util           0.17     segmented mod-30 sieve
       0.39  Math::Prime::Util::PP       0.24     Perl (Lehmer's method)
       0.9   Math::Prime::Util           0.01     mod-30 sieve
       2.9   Math::Prime::FastSieve      0.12     decent odd-number sieve
      11.7   Math::Prime::XS             0.26     needs some optimization
      15.0   Bit::Vector                 7.2
      48.9   Math::Prime::Util::PP       0.14     Perl (fastest I know of)
     170.0   Faster Perl sieve (net)     2012-01  array of odds
     548.1   RosettaCode sieve (net)     2012-06  simplistic Perl
    3048.1   Math::Primality             0.08     Perl + Math::GMPz
  ~11000     Math::Primality             0.04     Perl + Math::GMPz
  >20000     Math::Big                   1.12     Perl, > 26GB RAM used

Python's standard modules are very slow: MPMATH v0.17 C<primepi> takes 169.5s
and 25+ GB of RAM.  SymPy 0.7.1 C<primepi> takes 292.2s.  However there are
very fast solutions written by Robert William Hanks (included in the xt/
directory of this distribution): pure Python in 12.1s and NUMPY in 2.8s.

=head2 PRIMALITY TESTING

C<is_prime>: my impressions for various sized inputs:

   Module                   1-10 digits  10-20 digits  BigInts
   -----------------------  -----------  ------------  --------------
   Math::Prime::Util        Very fast    Very fast     Slow / Very Fast (1)
   Math::Prime::XS          Very fast    Very slow (2) --
   Math::Prime::FastSieve   Very fast    N/A (3)       --
   Math::Primality          Very slow    Very slow     Fast
   Math::Pari               Slow         OK            OK / Fast (4)

   (1) Without / With L<Math::Prime::Util::GMP> installed.
   (2) Trial division only.  Very fast if every factor is tiny.
   (3) Too much memory to hold the sieve (11dig = 6GB, 12dig = ~50GB)
   (4) By default L<Math::Pari> installs Pari 2.1.7, which uses 10 M-R tests
       for isprime and is not fast.  See notes below for 2.3.5.


The differences are in the implementations:

=over 4

=item L<Math::Prime::Util>

first does simple divisibility tests to quickly recognize composites, then
looks in the sieve for a fast bit lookup if possible (default up to 30,000
but can be expanded via C<prime_precalc>).  Next, for relatively small inputs,
a deterministic set of Miller-Rabin tests are used, while for larger inputs
a strong BPSW test is performed.  For native integers, this is faster than
any of the other modules.  With Bigints, you need the L<Math::Prime::Util::GMP>
module installed to get good performance.  With that installed, it is about
2x faster than Math::Primality and 10x faster than Math::Pari (default 2.1.7).

=item L<Math::Prime::XS>

does trial divisions, which is wonderful if the input has a small factor
(or is small itself).  If given a large prime it can be tens of thousands of
times slower than MPU.  It does not support bigints.

=item L<Math::Prime::FastSieve>

only works in a sieved range, which is really fast if you can do it
(M::P::U will do the same if you call C<prime_precalc>).  Larger inputs
just need too much time and memory for the sieve.

=item L<Math::Primality>

uses GMP (in Perl) for all work.  Under ~32-bits it uses 2 or 3 MR tests,
while above 4759123141 it performs a BPSW test.  This is great for
bigints over 2^64, but it is significantly slower than native precision
tests.  With 64-bit numbers it is generally an order of magnitude or more
slower than any of the others.  Once bigints are being used, its
performance is quite good.  It is faster than this module unless
L<Math::Prime::Util::GMP> has been installed, in which case Math::Prime::Util
is faster.

=item L<Math::Pari>

has some very effective code, but it has some overhead to get to it from
Perl.  That means for small numbers it is relatively slow: an order of
magnitude slower than M::P::XS and M::P::Util (though arguably this is
only important for benchmarking since "slow" is ~2 microseconds).  Large
numbers transition over to smarter tests so don't slow down much.  With
the default Pari version, C<isprime> will do M-R tests for 10 randomly
chosen bases, but can perform a Pocklington-Lehmer proof if requested using
C<isprime(x,1)>.  Both could fail to identify a composite.  If pari 2.3.5
is used instead (this requires hand-building the Math::Pari module) then
the options are quite different.  C<ispseudoprime(x,0)> performs a strong
BPSW test, while C<isprime> now performs a primality proof using a fast
implementation of the APRCL method.  While the APRCL method is very fast
it is still much, much slower than a BPSW probable prime test for large inputs.

=back

=head2 FACTORING

Factoring performance depends on the input, and the algorithm choices used
are still being tuned.  L<Math::Factor::XS> is very fast when given input with
only small factors, but it slows down rapidly as the smallest factor increases
in size.  For numbers larger than 32 bits, L<Math::Prime::Util> can be 100x or
more faster (a number with only very small factors will be nearly identical,
while a semiprime with large factors will be the extreme end).  L<Math::Pari>
is much slower with native sized inputs, probably due to calling
overhead.  For bigints, the L<Math::Prime::Util::GMP> module is needed or
performance will be far worse than Math::Pari.  With the GMP module,
performance is pretty similar from 20 through 70 digits, which the caveat
that the current MPU factoring uses more memory for 60+ digit numbers.


L<This slide presentation|http://math.boisestate.edu/~liljanab/BOISECRYPTFall09/Jacobsen.pdf>
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
and L<Pari|http://pari.math.u-bordeaux.fr/>.  The latest yafu should cover most
uses, with GGNFS likely only providing a benefit for numbers large enough to
warrant distributed processing.

=head2 PRIMALITY PROVING

The C<n-1> proving algorithm in L<Math::Prime::Util::GMP> compares well to
the version including in Pari.  Both are pretty fast to about 60 digits, and
work reasonably well to 80 or so before starting to take over many minutes per
number on a fast computer.  Version 0.09 and newer of MPU::GMP contain an
ECPP implementation that, while not state of the art compared to closed source
solutions, works quite well.
It averages less than a second for proving 200-digit primes
including creating a certificate.  Times below 200 digits are faster than
Pari 2.3.5's APR-CL proof.  For larger inputs the bottleneck is a limited set
of discriminants, and time becomes more variable.  There is a larger set of
discriminants on github that help, with 300-digit primes taking ~5 seconds on
average and typically under a minute for 500-digits.  For primality proving
with very large numbers, I recommend L<Primo|http://www.ellipsa.eu/>.

=head2 RANDOM PRIME GENERATION

Seconds per prime for random prime generation on a circa-2009 workstation,
with L<Math::BigInt::GMP>, L<Math::Prime::Util::GMP>, and
L<Math::Random::ISAAC::XS> installed.

  bits    random   +testing  rand_prov   Maurer   CPMaurer
  -----  --------  --------  ---------  --------  --------
     64    0.0001  +0.000008   0.0002     0.0001    0.022
    128    0.0020  +0.00023    0.011      0.063     0.057
    256    0.0034  +0.0004     0.058      0.13      0.16
    512    0.0097  +0.0012     0.28       0.28      0.41
   1024    0.060   +0.0060     0.65       0.65      2.19
   2048    0.57    +0.039      4.8        4.8      10.99
   4096    6.24    +0.25      31.9       31.9      79.71
   8192   58.6     +1.61     234.0      234.0     947.3

  random    = random_nbit_prime  (results pass BPSW)
  random+   = additional time for 3 M-R and a frobenius test
  rand_prov = random_proven_prime
  maurer    = random_maurer_prime
  CPMaurer  = Crypt::Primes::maurer

L</random_nbit_prime> is reasonably fast, and for most purposes should
suffice.  For cryptographic purposes, one may want additional tests or a
proven prime.  Additional tests are quite cheap, as shown by the time for
three extra M-R and a Frobenius test.  At these bit sizes, the chances a
composite number passes BPSW, three more M-R tests, and a Frobenius test
is I<extraordinarily> small.

L</random_proven_prime> provides a randomly selected prime with an optional
certificate, without specifying the particular method.  Below 512 bits,
using L</is_provable_prime>(L</random_nbit_prime>) is typically faster
than Maurer's algorithm, but becomes quite slow as the bit size increases.
This leaves the decision of the exact method of proving the result to the
implementation.

L</random_maurer_prime> constructs a provable prime.  A primality test is
run on each intermediate, and it also constructs a complete primality
certificate which is verified at the end (and can be returned).  While the
result is uniformly distributed, only about 10% of the primes in the range
are selected for output.  This is a result of the FastPrime algorithm and
is usually unimportant.

L<Crypt::Primes/maurer> times are included for comparison.  It is pretty
fast for small sizes but gets slow as the size increases.  It does not
perform any primality checks on the intermediate results or the final
result (I highly recommended you run a primality test on the output).
Additionally important for servers, L<Crypt::Primes/maurer> uses excessive
system entropy and can grind to a halt if C</dev/random> is exhausted
(it can take B<days> to return).  The times above are on a machine running
L<HAVEGED|http://www.issihosts.com/haveged/>
so never waits for entropy.  Without this, the times would be much higher.


=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 ACKNOWLEDGEMENTS

Eratosthenes of Cyrene provided the elegant and simple algorithm for finding
primes.

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

Pierre Dusart, "Estimates of Some Functions Over Primes without R.H.", preprint, 2010.  Updates to the best non-RH bounds for prime count and nth prime.  L<http://arxiv.org/abs/1002.0442/>

=item *

Pierre Dusart, "Autour de la fonction qui compte le nombre de nombres premiers", PhD thesis, 1998.  In French.  The mathematics is readable and highly recommended reading if you're interesting in prime number bounds.  L<http://www.unilim.fr/laco/theses/1998/T1998_01.html>

=item *

Gabriel Mincu, "An Asymptotic Expansion", I<Journal of Inequalities in Pure and Applied Mathematics>, v4, n2, 2003.  A very readable account of Cipolla's 1902 nth prime approximation.  L<http://www.emis.de/journals/JIPAM/images/153_02_JIPAM/153_02.pdf>

=item *

Christian Bau, "The Extended Meissel-Lehmer Algorithm", 2003, preprint with example C++ implementation.  Very detailed implementation-specific paper which was used for the implementation here.  Highly recommended for implementing a sieve-based LMO.  L<http://cs.swan.ac.uk/~csoliver/ok-sat-library/OKplatform/ExternalSources/sources/NumberTheory/ChristianBau/>

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

Ueli M. Maurer, "Fast Generation of Prime Numbers and Secure Public-Key Cryptographic Parameters", 1995.  Generating random provable primes by building up the prime.  L<http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.26.2151>

=item *

Pierre-Alain Fouque and Mehdi Tibouchi, "Close to Uniform Prime Number Generation With Fewer Random Bits", pre-print, 2011.  Describes random prime distributions, their algorithm for creating random primes using few random bits, and comparisons to other methods.  Definitely worth reading for the discussions of uniformity.  L<http://eprint.iacr.org/2011/481>

=item *

Douglas A. Stoll and Patrick Demichel , "The impact of ζ(s) complex zeros on π(x) for x E<lt> 10^{10^{13}}", L<Mathematics of Computation>, v80, n276, pp 2381-2394, October 2011.  L<http://www.ams.org/journals/mcom/2011-80-276/S0025-5718-2011-02477-4/home.html>

=item *

L<OEIS: Primorial|http://oeis.org/wiki/Primorial>

=item *

Walter M. Lioen and Jan van de Lune, "Systematic Computations on Mertens' Conjecture and Dirichlet's Divisor Problem by Vectorized Sieving", in I<From Universal Morphisms to Megabytes>, Centrum voor Wiskunde en Informatica, pp. 421-432, 1994.  Describes a nice way to compute a range of Möbius values.  L<http://walter.lioen.com/papers/LL94.pdf>

=item *

Marc Deléglise and Joöl Rivat, "Computing the summation of the Möbius function", I<Experimental Mathematics>, v5, n4, pp 291-295, 1996.  Enhances the Möbius computation in Lioen/van de Lune, and gives a very efficient way to compute the Mertens function.  L<http://projecteuclid.org/euclid.em/1047565447>

=item *

Manuel Benito and Juan L. Varona, "Recursive formulas related to the summation of the Möbius function", I<The Open Mathematics Journal>, v1, pp 25-34, 2007.  Among many other things, shows a simple formula for computing the Mertens functions with only n/3 Möbius values (not as fast as Deléglise and Rivat, but really simple).  L<http://www.unirioja.es/cu/jvarona/downloads/Benito-Varona-TOMATJ-Mertens.pdf>

=item *

John Brillhart, D. H. Lehmer, and J. L. Selfridge, "New Primality Criteria and Factorizations of 2^m +/- 1", Mathematics of Computation, v29, n130, Apr 1975, pp 620-647.  L<http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf>

=back


=head1 COPYRIGHT

Copyright 2011-2013 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
