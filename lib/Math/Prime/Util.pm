package Math::Prime::Util;
use strict;
use warnings;
use Carp qw/croak confess carp/;

BEGIN {
  $Math::Prime::Util::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::VERSION = '0.41';
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
      is_power
      miller_rabin miller_rabin_random
      lucas_sequence
      primes
      forprimes forcomposites fordivisors forpart
      prime_iterator prime_iterator_object
      next_prime  prev_prime
      prime_count
      prime_count_lower prime_count_upper prime_count_approx
      nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
      twin_prime_count twin_prime_count_approx
      nth_twin_prime nth_twin_prime_approx
      random_prime random_ndigit_prime random_nbit_prime random_strong_prime
      random_proven_prime random_proven_prime_with_cert
      random_maurer_prime random_maurer_prime_with_cert
      random_shawe_taylor_prime random_shawe_taylor_prime_with_cert
      primorial pn_primorial consecutive_integer_lcm
      gcd lcm factor factor_exp all_factors divisors valuation invmod vecsum
      moebius mertens euler_phi jordan_totient exp_mangoldt liouville
      partitions
      chebyshev_theta chebyshev_psi
      divisor_sum carmichael_lambda
      kronecker binomial znorder znprimroot znlog legendre_phi
      ExponentialIntegral LogarithmicIntegral RiemannZeta RiemannR
  );
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);

# These are only exported if specifically asked for
push @EXPORT_OK, (qw/trial_factor fermat_factor holf_factor squfof_factor prho_factor pbrent_factor pminus1_factor pplus1_factor/);

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
  1;
}

BEGIN {

  # Separate lines to keep compatible with default from 5.6.2.
  # We could alternately use Config's $Config{uvsize} for MAXBITS
  use constant OLD_PERL_VERSION=> $] < 5.008;
  use constant MPU_MAXBITS     => (~0 == 4294967295) ? 32 : 64;
  use constant MPU_32BIT       => MPU_MAXBITS == 32;
  use constant MPU_MAXPARAM    => MPU_32BIT ? 4294967295 : 18446744073709551615;
  use constant MPU_MAXDIGITS   => MPU_32BIT ?         10 : 20;
  use constant MPU_MAXPRIME    => MPU_32BIT ? 4294967291 : 18446744073709551557;
  use constant MPU_MAXPRIMEIDX => MPU_32BIT ?  203280221 :   425656284035217743;
  use constant UVPACKLET       => MPU_32BIT ?        'L' : 'Q';

  eval {
    return 0 if defined $ENV{MPU_NO_XS} && $ENV{MPU_NO_XS} == 1;
    require XSLoader;
    XSLoader::load(__PACKAGE__, $Math::Prime::Util::VERSION);
    prime_precalc(0);
    $_Config{'maxbits'} = _XS_prime_maxbits();
    $_Config{'xs'} = 1;
    1;
  } or do {
    carp "Using Pure Perl implementation: $@"
      unless defined $ENV{MPU_NO_XS} && $ENV{MPU_NO_XS} == 1;

    $_Config{'xs'} = 0;
    $_Config{'maxbits'} = MPU_MAXBITS;

    # Load PP front end code
    require Math::Prime::Util::PPFE;

    *next_prime    = \&Math::Prime::Util::_generic_next_prime;
    *prev_prime    = \&Math::Prime::Util::_generic_prev_prime;
    *prime_count   = \&Math::Prime::Util::_generic_prime_count;
    *factor        = \&Math::Prime::Util::_generic_factor;
    *factor_exp    = \&Math::Prime::Util::_generic_factor_exp;
  };

  # aliases for deprecated names.  Will eventually be removed.
  *all_factors = \&divisors;
  *miller_rabin = \&is_strong_pseudoprime;

  $_Config{'nobigint'} = 0;
  $_Config{'gmp'} = 0;
  # See if they have the GMP module and haven't requested it not to be used.
  if (!defined $ENV{MPU_NO_GMP} || $ENV{MPU_NO_GMP} != 1) {
    $_Config{'gmp'} = 1 if eval { require Math::Prime::Util::GMP;
                                  Math::Prime::Util::GMP->import();
                                  1; };
  }
}

croak "Perl and XS don't agree on bit size"
      if $_Config{'xs'} && MPU_MAXBITS != _XS_prime_maxbits();

$_Config{'maxparam'}    = MPU_MAXPARAM;
$_Config{'maxdigits'}   = MPU_MAXDIGITS;
$_Config{'maxprime'}    = MPU_MAXPRIME;
$_Config{'maxprimeidx'} = MPU_MAXPRIMEIDX;
$_Config{'assume_rh'}   = 0;
$_Config{'verbose'}     = 0;
$_Config{'irand'}       = undef;

# used for code like:
#    return _XS_foo($n)  if $n <= $_XS_MAXVAL
# which builds into one scalar whether XS is available and if we can call it.
my $_XS_MAXVAL = $_Config{'xs'}  ?  MPU_MAXPARAM  :  -1;
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
  foreach my $param (keys %params) {
    my $value = $params{$param};
    $param = lc $param;
    # dispatch table should go here.
    if      ($param eq 'xs') {
      $_Config{'xs'} = ($value) ? 1 : 0;
      $_XS_MAXVAL = $_Config{'xs'}  ?  MPU_MAXPARAM  :  -1;
    } elsif ($param eq 'gmp') {
      $_Config{'gmp'} = ($value) ? 1 : 0;
      $_HAVE_GMP = $_Config{'gmp'};
      _XS_set_callgmp($_HAVE_GMP) if $_Config{'xs'};
    } elsif ($param eq 'nobigint') {
      $_Config{'nobigint'} = ($value) ? 1 : 0;
    } elsif ($param eq 'irand') {
      croak "irand must supply a sub" unless (!defined $value) || (ref($value) eq 'CODE');
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

sub _bigint_to_int {
  return (OLD_PERL_VERSION) ? unpack(UVPACKLET,pack(UVPACKLET,$_[0]->bstr))
                            : int($_[0]->bstr);
}
sub _to_bigint {
  do { require Math::BigInt;  Math::BigInt->import(try=>"GMP,Pari"); }
    unless defined $Math::BigInt::VERSION;
  return Math::BigInt->new("$_[0]");
}
sub _reftyped {
  my $ref0 = ref($_[0]);
  if ($ref0) {
    return  ($ref0 eq ref($_[1])) ?  $_[1]  :  $ref0->new("$_[1]");
  }
  my $strn = "$_[1]";
  return $_[1] if $strn <= ~0;
  do { require Math::BigInt;  Math::BigInt->import(try=>"GMP,Pari"); }
    unless defined $Math::BigInt::VERSION;
  return Math::BigInt->new($strn);
}


#*_validate_positive_integer = \&Math::Prime::Util::PP::_validate_positive_integer;
sub _validate_positive_integer {
  my($n, $min, $max) = @_;
  croak "Parameter must be defined" if !defined $n;
  if (ref($n) eq 'CODE') {
    $_[0] = $_[0]->();
    $n = $_[0];
  }
  if (ref($n) eq 'Math::BigInt') {
    croak "Parameter '$n' must be a positive integer"
      if $n->sign() ne '+' || !$n->is_int();
    $_[0] = _bigint_to_int($_[0])
      if $n <= (OLD_PERL_VERSION ? 562949953421312 : ''.~0);
  } else {
    my $strn = "$n";
    croak "Parameter '$strn' must be a positive integer"
      if $strn =~ tr/0123456789//c && $strn !~ /^\+?\d+$/;
    if ($n <= (OLD_PERL_VERSION ? 562949953421312 : ''.~0)) {
      $_[0] = $strn if ref($n);
    } else {
      #$_[0] = Math::BigInt->new($strn)
      $_[0] = _to_bigint($strn);
    }
  }
  $_[0]->upgrade(undef) if ref($_[0]) && $_[0]->upgrade();
  croak "Parameter '$_[0]' must be >= $min" if defined $min && $_[0] < $min;
  croak "Parameter '$_[0]' must be <= $max" if defined $max && $_[0] > $max;
  1;
}


#############################################################################

sub primes {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    _validate_num($low) || _validate_positive_integer($low);
    _validate_num($high) || _validate_positive_integer($high);
  } else {
    ($low,$high) = (2, $low);
    _validate_num($high) || _validate_positive_integer($high);
  }

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
    require Math::Prime::Util::PP;
    return Math::Prime::Util::PP::primes($low,$high);
  }

  # Decide the method to use.  We have four to choose from:
  #  1. Trial     No memory, no overhead, but more time per prime.
  #  2. Sieve     Monolithic cached sieve.
  #  3. Erat      Monolithic uncached sieve.
  #  4. Segment   Segment sieve.  Never a bad decision.

  if (($low+1) >= $high ||                      # Tiny range, or
      $high > 10**14 && ($high-$low) < 50000) { # Small relative range

      $sref = trial_primes($low, $high);

  } elsif ($high <= (65536*30) ||                # Very small, or
           $high <= _get_prime_cache_size()) {   # already in the main cache.

      $sref = sieve_primes($low, $high);

  } else {

      $sref = segment_primes($low, $high);

  }

  # We could return an array ref in scalar context, array in list context with:
  #   return (wantarray) ? @{$sref} : $sref;
  # but I think the dual interface could be confusing, albeit often handy.
  return $sref;
}

#############################################################################
# Random primes.  These are front end functions that do input validation,
# load the RandomPrimes module, and call its function.

sub random_prime {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    _validate_num($low) || _validate_positive_integer($low);
    _validate_num($high) || _validate_positive_integer($high);
  } else {
    ($low,$high) = (2, $low);
    _validate_num($high) || _validate_positive_integer($high);
  }
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_prime($low,$high);
}

sub random_ndigit_prime {
  my($digits) = @_;
  _validate_num($digits, 1) || _validate_positive_integer($digits, 1);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_ndigit_prime($digits);
}

sub random_nbit_prime {
  my($bits) = @_;
  _validate_num($bits, 2) || _validate_positive_integer($bits, 2);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_nbit_prime($bits);
}

sub random_maurer_prime {
  my($bits) = @_;
  _validate_num($bits, 2) || _validate_positive_integer($bits, 2);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_maurer_prime($bits);
}

sub random_maurer_prime_with_cert {
  my($bits) = @_;
  _validate_num($bits, 2) || _validate_positive_integer($bits, 2);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_maurer_prime_with_cert($bits);
}

sub random_shawe_taylor_prime {
  my($bits) = @_;
  _validate_num($bits, 2) || _validate_positive_integer($bits, 2);
  require Math::Prime::Util::RandomPrimes;
  my ($n, $cert) = Math::Prime::Util::RandomPrimes::random_shawe_taylor_prime_with_cert($bits);
  return $n;
}

sub random_shawe_taylor_prime_with_cert {
  my($bits) = @_;
  _validate_num($bits, 2) || _validate_positive_integer($bits, 2);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_shawe_taylor_prime_with_cert($bits);
}

sub random_strong_prime {
  my($bits) = @_;
  _validate_num($bits, 128) || _validate_positive_integer($bits, 128);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_strong_prime($bits);
}

sub random_proven_prime {
  my($bits) = @_;
  _validate_num($bits, 2) || _validate_positive_integer($bits, 2);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_proven_prime($bits);
}

sub random_proven_prime_with_cert {
  my($bits) = @_;
  _validate_num($bits, 2) || _validate_positive_integer($bits, 2);
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_proven_prime_with_cert($bits);
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

  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::miller_rabin_random($n, $k, $seed);
}


#############################################################################
# These functions almost always return bigints, so there is no XS
# implementation.  Try to run the GMP version, and if it isn't available,
# load PP and call it.

sub primorial {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::primorial) {
    return _reftyped($_[0], Math::Prime::Util::GMP::primorial($n));
  }
  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::primorial($n);
}

sub pn_primorial {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::pn_primorial) {
    return _reftyped($_[0], Math::Prime::Util::GMP::pn_primorial($n));
  }

  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::primorial(nth_prime($n));
}

sub consecutive_integer_lcm {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  return 0 if $n < 1;

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::consecutive_integer_lcm) {
    return _reftyped($_[0],Math::Prime::Util::GMP::consecutive_integer_lcm($n));
  }
  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::consecutive_integer_lcm($n);
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

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::partitions) {
    return _reftyped($_[0],Math::Prime::Util::GMP::partitions($n));
  }

  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::partitions($n);
}


#############################################################################
# forprimes, forcomposites, fordivisors.
# These are used when the XS code can't handle it.

sub _generic_forprimes {
  my($sub, $beg, $end) = @_;
  if (!defined $end) { $end = $beg; $beg = 2; }
  _validate_positive_integer($beg);
  _validate_positive_integer($end);
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

sub _generic_forcomposites {
  my($sub, $beg, $end) = @_;
  if (!defined $end) { $end = $beg; $beg = 4; }
  _validate_positive_integer($beg);
  _validate_positive_integer($end);
  $beg = 4 if $beg < 4;
  $end = Math::BigInt->new(''.~0) if ref($end) ne 'Math::BigInt' && $end == ~0;
  {
    my $pp;
    local *_ = \$pp;
    for ( ; $beg <= $end ; $beg++ ) {
      if (!is_prime($beg)) {
        $pp = $beg;
        $sub->();
      }
    }
  }
}

sub _generic_fordivisors {
  my($sub, $n) = @_;
  _validate_positive_integer($n);
  my @divisors = divisors($n);
  {
    my $pp;
    local *_ = \$pp;
    foreach my $d (@divisors) {
      $pp = $d;
      $sub->();
    }
  }
}

#############################################################################
# Iterators

sub prime_iterator {
  my($start) = @_;
  $start = 0 unless defined $start;
  _validate_num($start) || _validate_positive_integer($start);
  my $p = ($start > 0) ? $start-1 : 0;
  # This works fine:
  #   return sub { $p = next_prime($p); return $p; };
  # but we can optimize a little
  if (ref($p) ne 'Math::BigInt' && $p <= $_XS_MAXVAL) {
    return sub { $p = next_prime($p); return $p; };
  } elsif ($_HAVE_GMP) {
    return sub { $p = $p-$p+Math::Prime::Util::GMP::next_prime($p); return $p;};
  } else {
    require Math::Prime::Util::PP;
    return sub { $p = Math::Prime::Util::PP::next_prime($p); return $p; }
  }
}

sub prime_iterator_object {
  my($start) = @_;
  require Math::Prime::Util::PrimeIterator;
  return Math::Prime::Util::PrimeIterator->new($start);
}

#############################################################################
# Front ends to functions.
#
# These will do input validation, then call the appropriate internal function
# based on the input (XS, GMP, PP).
#############################################################################

sub _generic_next_prime {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  if ($_HAVE_GMP) {
    return _reftyped($_[0], Math::Prime::Util::GMP::next_prime($n));
  }

  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::next_prime($_[0]);
}

sub _generic_prev_prime {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  if ($_HAVE_GMP) {
    return _reftyped($_[0], Math::Prime::Util::GMP::prev_prime($n));
  }

  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::prev_prime($_[0]);
}

sub _generic_prime_count {
  my($low,$high) = @_;
  if (defined $high) {
    _validate_num($low) || _validate_positive_integer($low);
    _validate_num($high) || _validate_positive_integer($high);
  } else {
    ($low,$high) = (2, $low);
    _validate_num($high) || _validate_positive_integer($high);
  }
  return 0 if $high < 2  ||  $low > $high;

  # We can relax these constraints if MPU::GMP gets a fast implementation.
  return Math::Prime::Util::GMP::prime_count($low,$high) if $_HAVE_GMP
                       && defined &Math::Prime::Util::GMP::prime_count
                       && (   (ref($high) eq 'Math::BigInt')
                           || (($high-$low) < int($low/1_000_000))
                          );
  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::prime_count($low,$high);
}

sub _generic_factor {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  if ($_HAVE_GMP) {
    my @factors;
    if ($n != 1) {
      @factors = Math::Prime::Util::GMP::factor($n);
      if (ref($_[0]) eq 'Math::BigInt') {
        @factors = map { ($_ > ~0) ? Math::BigInt->new(''.$_) : $_ } @factors;
      }
    }
    return @factors;
  }

  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::factor($n);
}

sub _generic_factor_exp {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  my %exponents;
  my @factors = grep { !$exponents{$_}++ } factor($n);
  return scalar @factors unless wantarray;
  return (map { [$_, $exponents{$_}] } @factors);
}


sub lucas_sequence {
  my($n, $P, $Q, $k) = @_;
  _validate_num($n) || _validate_positive_integer($n);
  croak("Invalid input to lucas_sequence:  modulus n must be > 1") if $n <= 1;
  _validate_num($k) || _validate_positive_integer($k);
  { my $testP = (!defined $P || $P >= 0) ? $P : -$P;
    _validate_num($testP) || _validate_positive_integer($testP); }
  { my $testQ = (!defined $Q || $Q >= 0) ? $Q : -$Q;
    _validate_num($testQ) || _validate_positive_integer($testQ); }

  return _XS_lucas_sequence($n, $P, $Q, $k)
    if ref($_[0]) ne 'Math::BigInt' && $n <= $_XS_MAXVAL
    && ref($_[3]) ne 'Math::BigInt' && $k <= $_XS_MAXVAL;

  if ($_HAVE_GMP && defined &Math::Prime::Util::GMP::lucas_sequence) {
    return map { ($_ > ''.~0) ? Math::BigInt->new(''.$_) : $_ }
           Math::Prime::Util::GMP::lucas_sequence($n, $P, $Q, $k);
  }
  require Math::Prime::Util::PP;
  return map { ($_ <= ''.~0) ? _bigint_to_int($_) : $_ }
         Math::Prime::Util::PP::lucas_sequence($n, $P, $Q, $k);
}


#############################################################################

# Return just the non-cert portion.
sub is_provable_prime {
  my($n) = @_;
  return 0 if defined $n && $n < 2;
  _validate_num($n) || _validate_positive_integer($n);

  return is_prime($n) if $n <= $_XS_MAXVAL;
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
    my $isp = is_prime($n);
    return ($isp, '') unless $isp == 2;
    return (2, $header . "Type Small\nN $n\n");
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
    require Math::Prime::Util::PrimalityProving;
    return ($isp, Math::Prime::Util::PrimalityProving::convert_array_cert_to_string($cert));
  }

  {
    my $isp = is_prob_prime($n);
    return ($isp, '') if $isp == 0;
    return (2, $header . "Type Small\nN $n\n") if $isp == 2;
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

  require Math::Prime::Util::PrimalityProving;
  #my ($isp, $pref) = Math::Prime::Util::PrimalityProving::primality_proof_lucas($n);
  my ($isp, $pref) = Math::Prime::Util::PrimalityProving::primality_proof_bls75($n);
  carp "proved $n is not prime\n" if !$isp;
  return ($isp, $pref);
}


sub verify_prime {
  my @cdata = @_;

  require Math::Prime::Util::PrimalityProving;
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

#############################################################################

sub RiemannZeta {
  my($n) = @_;
  croak("Invalid input to ReimannZeta:  x must be > 0") if $n <= 0;

  return _XS_RiemannZeta($n)
    if !defined $bignum::VERSION && ref($n) ne 'Math::BigFloat' && $n <= $_XS_MAXVAL;
  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::RiemannZeta($n);
}

sub RiemannR {
  my($n) = @_;
  croak("Invalid input to ReimannR:  x must be > 0") if $n <= 0;

  return _XS_RiemannR($n)
    if !defined $bignum::VERSION && ref($n) ne 'Math::BigFloat' && $n <= $_XS_MAXVAL;
  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::RiemannR($n);
}

sub ExponentialIntegral {
  my($n) = @_;
  return $_Neg_Infinity if $n == 0;
  return 0              if $n == $_Neg_Infinity;
  return $_Infinity     if $n == $_Infinity;

  return _XS_ExponentialIntegral($n)
   if !defined $bignum::VERSION && ref($n) ne 'Math::BigFloat' && $_Config{'xs'};

  require Math::Prime::Util::PP;
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

  require Math::Prime::Util::PP;
  return Math::Prime::Util::PP::LogarithmicIntegral($n);
}

#############################################################################

use Math::Prime::Util::MemFree;

1;

__END__


# ABSTRACT: Utilities related to prime numbers, including fast generators / sievers

=pod

=encoding utf8

=for stopwords forprimes forcomposites fordivisors forpart Möbius Deléglise totient moebius mertens liouville znorder irand primesieve uniqued k-tuples von SoE pari yafu fonction qui compte le nombre nombres voor PhD superset sqrt(N) gcd(A^M k-th (10001st primegen libtommath kronecker znprimroot znlog gcd lcm invmod untruncated vecsum


=head1 NAME

Math::Prime::Util - Utilities related to prime numbers, including fast sieves and factoring


=head1 VERSION

Version 0.41


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
  # Or for the composites in a range
  forcomposites { say if is_strong_pseudoprime($_,2) } 10000, 10**6;

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
  @divisors = divisors( $n );
  # Or just apply a block for each one
  fordivisors  { $sum += $_ + $_*$_ }  $n;

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
  # Show all prime partitions of 25
  forpart { say "@_" unless scalar grep { !is_prime($_) } @_ } 25;

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
  my $rand_prime = random_shawe_taylor_prime(256);  # as above


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
your program.  Also note that L<Math::Pari> is not thread-safe (and will
crash as soon as it is loaded in threads), so if you use
L<Math::BigInt::Pari> rather than L<Math::BigInt::GMP> or the default backend,
things will go pear-shaped.

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
and earlier with bignums.  For example, the default BigInt backend in older
versions of Perl will sometimes convert small results to doubles, resulting
in corrupted output.

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

Returns 0 is the number is composite, 1 if it is probably prime, and 2 if
it is definitely prime.  For numbers smaller than C<2^64> it will only
return 0 (composite) or 2 (definitely prime), as this range has been
exhaustively tested and has no counterexamples.  For larger numbers,
an extra-strong BPSW test is used.
If L<Math::Prime::Util::GMP> is installed, some additional primality tests
are also performed, and a quick attempt is made to perform a primality
proof, so it will return 2 for many other inputs.

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
a different test (e.g. L</is_frobenius_underwood_pseudoprime>) or using
additional M-R tests with random bases with L</miller_rabin_random>.
Even better, make sure L<Math::Prime::Util::GMP> is installed and use
L</is_provable_prime> which should be reasonably fast for sizes under
2048 bits.  Another possibility is to use
L<Math::Prime::Util/random_maurer_prime> or
L<Math::Prime::Util/random_shawe_taylor_prime> which construct random
provable primes.


=head2 primes

Returns all the primes between the lower and upper limits (inclusive), with
a lower limit of C<2> if none is given.

An array reference is returned (with large lists this is much faster and uses
less memory than returning an array directly).

  my $aref1 = primes( 1_000_000 );
  my $aref2 = primes( 1_000_000_000_000, 1_000_000_001_000 );

  my @primes = @{ primes( 500 ) };

  print "$_\n" for @{primes(20,100)};

Sieving will be done if required.  The algorithm used will depend on the range
and whether a sieve result already exists.  Possibilities include primality
testing (for very small ranges), a Sieve of Eratosthenes using wheel
factorization, or a segmented sieve.


=head2 next_prime

  $n = next_prime($n);

Returns the next prime greater than the input number.  The result will be a
bigint if it can not be exactly represented in the native int type
(larger than C<4,294,967,291> in 32-bit Perl;
larger than C<18,446,744,073,709,551,557> in 64-bit).


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
of primes.  Nesting and use in threads are allowed.

Math::BigInt objects may be used for the range.

For some uses an iterator (L</prime_iterator>, L</prime_iterator_object>)
or a tied array (L<Math::Prime::Util::PrimeArray>) may be more convenient.
Objects can be passed to functions, and allow early loop exits.


=head2 forcomposites

  forcomposites { say } 1000;
  forcomposites { say } 2000,2020;

Given a block and either an end number or a start and end pair, calls the
block for each composite in the inclusive range.  The composites,
L<OEIS A002808|http://oeis.org/A002808>, are the numbers greater than 1
which are not prime:  C<4, 6, 8, 9, 10, 12, 14, 15, ...>


=head2 fordivisors

  fordivisors { $prod *= $_ } $n;

Given a block and a non-negative number C<n>, the block is called with
C<$_> set to each divisor in sorted order.  Also see L</divisor_sum>.


=head2 forpart

  forpart { say "@_" } 25;           # unrestricted partitions
  forpart { say "@_" } 25,{n=>5}     # ... with exactly 5 values
  forpart { say "@_" } 25,{nmax=>5}  # ... with <=5 values

Given a non-negative number C<n>, the block is called with C<@_> set to
the array of additive integer partitions.  The operation is very similar
to the C<forpart> function in Pari/GP 2.6.x, though the ordering is
different.  The algorithm is ZS1 from Zoghbi and Stojmenović (1998), hence
the ordering is identical to that of L<Integer::Partition>.
Use L</partitions> to get just the count of unrestricted partitions.


An optional hash reference may be given to produce restricted partitions.
Each value must be a non-negative integer.  The allowable keys are:

  n       restrict to exactly this many values
  amin    all elements must be at least this value
  amax    all elements must be at most this value
  nmin    the array must have at least this many values
  nmax    the array must have at most this many values

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


=head2 twin_prime_count

Similar to prime count, but returns the count of twin primes (primes C<p>
where C<p+2> is also prime).  Takes either a single number indicating a count
from 2 to the argument, or two numbers indicating a range.

The primes being counted are the first value, so a range of C<(3,5)> will
return a count of two, because both C<3> and C<5> are counted as twin primes.
A range of C<(12,13)> will return a count of zero, because neither C<12+2>
nor C<13+2> are prime.  In contrast, C<primesieve> requires all elements of
a constellation to be within the range to be counted, so would return one for
the first example (C<5> is not counted because its pair C<7> is not in the
range).

There is no useful formula known for this, unlike prime counts.  We sieve
for the answer, using some small table acceleration.

=head2 twin_prime_count_approx

Returns an approximation to the twin prime count of C<n>.  This returns
quickly and has a very small error for large values.  The method used is
conjecture B of Hardy and Littlewood 1922, as stated in
Sebah and Gourdon 2002.  For inputs under 10M, a correction factor is
additionally applied to reduce the mean squared error.


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
C<10^14th> prime (3475385758524527) takes under 30 seconds.  Think about
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
generate any primes.  For values where the nth prime is smaller than
C<2^64>, an inverse Riemann R function is used.  For larger values, uses the
Cipolla 1902 approximation with up to 2nd order terms, plus a third order
correction.


=head2 nth_twin_prime

Returns the Nth twin prime.  This is done via sieving and counting, so
is not very fast for large values.

=head2 nth_twin_prime_approx

Returns an approximation to the Nth twin prime.  A curve fit is used for
small inputs (under 1200), while for larger inputs a binary search is done
on the approximate twin prime count.


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


=head2 is_bpsw_prime

Given a positive number input, returns 0 (composite), 2 (definitely prime),
or 1 (probably prime), using the BPSW primality test (extra-strong variant).
Normally one of the L<Math::Prime::Util/is_prime> or
L<Math::Prime::Util/is_prob_prime> functions will suffice, but those
functions do pre-tests to find easy composites.  If you know this is not
necessary, then calling L</is_bpsw_prime> may save a small amount of time.


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

While this is an important theoretical algorithm, and makes an interesting
example, it is hard to overstate just how impractically slow it is in
practice.  It is not used for any purpose in non-theoretical work, as it is
literally B<millions> of times slower than other algorithms.  From R.P.
Brent, 2010:  "AKS is not a practical algorithm.  ECPP is much faster."
We have ECPP, and indeed it is much faster.

This implementation includes the v6 improvements from Lenstra as well as
further improvements from Bernstein and Voloch.  It runs substantially
faster than the original or v6 versions.  The GMP implementation uses
a binary segmentation method for modular polynomial multiplication
(see Bernstein's 2007 Quartic paper), which reduces to a single scalar
multiplication, at which GMP excels.  Because of this, the GMP
implementation is likely to be faster once the input is larger than C<2^32>.


=head2 is_power

  say "$n is a perfect square" if is_power($n, 2);
  say "$n is a perfect cube" if is_power($n, 3);
  say "$n is a ", is_power($n), "-th power";

Given a single positive integer input C<n>, returns k if C<n = p^k> for
some integer C<p E<gt> 1, k E<gt> 1>, and 0 otherwise.  The k returned is
the largest possible.  This can be used in a boolean statement to
determine if C<n> is a perfect power.

If given two arguments C<n> and C<k>, returns 1 if C<n> is a C<k-th> power,
and 0 otherwise.  For example, if C<k=2> then this detects perfect squares.

This corresponds to Pari/GP's C<ispower> function, with the limitations of
only integer arguments and no third argument may be given to return the root.


=head2 lucas_sequence

  my($U, $V, $Qk) = lucas_sequence($n, $P, $Q, $k)

Computes C<U_k>, C<V_k>, and C<Q_k> for the Lucas sequence defined by
C<P>,C<Q>, modulo C<n>.  The modular Lucas sequence is used in a
number of primality tests and proofs.
The following conditions must hold:
C< |P| E<lt> n>  ;
C< |Q| E<lt> n>  ;
C< k E<gt>= 0>  ;
C< n E<gt>= 2>.


=head2 gcd

Given a list of integers, returns the greatest common divisor.  This is
often used to test for L<coprimality|https://oeis.org/wiki/Coprimality>.

=head2 lcm

Given a list of integers, returns the least common multiple.  Note that we
follow the semantics of Mathematica, Pari, and Perl 6, re:

  lcm(0, n) = 0              Any zero in list results in zero return
  lcm(n,-m) = lcm(n, m)      We use the absolute values


=head2 vecsum

  say "Totient sum 500,000: ", vecsum(euler_phi(0,500_000));

Returns the sum of all arguments, each of which must be an integer.  This
is similar to List::Util's L<List::Util/sum0> function, but has a very
important difference.  List::Util turns all inputs into doubles and returns
a double, which will mean incorrect results with large integers.  C<vecsum>
sums (signed) integers and returns the untruncated result.  Processing is
done on native integers while possible.


=head2 invmod

  say "The inverse of 42 mod 2017 = ", invmod(42,2017);

Given two integers C<a> and C<n>, return the inverse of C<a> modulo C<n>.
If not defined, undef is returned.  If defined, then the return value
multiplied by C<a> equals C<1> modulo C<n>.

This results correspond to the Pari result of C<lift(Mod(1/a,n))>.  The
semantics with respect to negative arguments match Pari.  Notably, a
negative C<n> is negated, which is different from Math::BigInt, but in both
cases the return value is still congruent to C<1> modulo C<n> as expected.

=head2 valuation

  say "$n is divisible by 2 ", valuation($n,2), " times.";

Given integers C<n> and C<k>, returns the numbers of times C<n> is divisible
by C<k>.  This is a very limited version of the algebraic valuation meaning,
just applied to integers.
This corresponds to Pari's C<valuation> function.
C<0> is returned if C<n> or C<k> is one of the values C<-1>, C<0>, or C<1>.

=head2 moebius

  say "$n is square free" if moebius($n) != 0;
  $sum += moebius($_) for (1..200); say "Mertens(200) = $sum";
  say "Mertens(2000) = ", vecsum(moebius(0,2000));

Returns μ(n), the Möbius function (also known as the Moebius, Mobius, or
MoebiusMu function) for an integer input.  This function is 1 if
C<n = 1>, 0 if C<n> is not square free (i.e. C<n> has a repeated factor),
and C<-1^t> if C<n> is a product of C<t> distinct primes.  This is an
important function in prime number theory.  Like SAGE, we define
C<moebius(0) = 0> for convenience.

If called with two arguments, they define a range C<low> to C<high>, and the
function returns an array with the value of the Möbius function for every n
from low to high inclusive.  Large values of high will result in a lot of
memory use.  The algorithm used for ranges is Deléglise and Rivat (1996)
algorithm 4.1, which is a segmented version of Lioen and van de Lune (1994)
algorithm 3.2.

The return values are read-only constants.  This should almost never come up,
but it means trying to modify aliased return values will cause an
exception (modifying the returned scalar or array is fine).


=head2 mertens

  say "Mertens(10M) = ", mertens(10_000_000);   # = 1037

Returns M(n), the Mertens function for a non-negative integer input.  This
function is defined as C<sum(moebius(1..n))>, but calculated more efficiently
for large inputs.  For example, computing Mertens(100M) takes:

   time    approx mem
     0.4s      0.1MB   mertens(100_000_000)
     5.6s    880MB     vecsum(moebius(1,100_000_000))
   102s        0MB     $sum += moebius($_) for 1..100_000_000

The summation of individual terms via factoring is quite expensive in time,
though uses O(1) space.  Using the range version of moebius is much faster,
but returns a 100M element array which, even though they are shared constants,
is not good for memory at this size.
In comparison, this function will generate the equivalent output
via a sieving method that is relatively memory frugal and very fast.
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
function) for an integer value.  This is an arithmetic function which counts
the number of positive integers less than or equal to C<n> that are relatively
prime to C<n>.  Given the definition used, C<euler_phi> will return 0 for all
C<n E<lt> 1>.  This follows the logic used by SAGE.  Mathematica and Pari
return C<euler_phi(-n)> for C<n E<lt> 0>.  Mathematica returns 0 for C<n = 0>
while Pari raises an exception.

If called with two arguments, they define a range C<low> to C<high>, and the
function returns an array with the totient of every n from low to high
inclusive.


=head2 jordan_totient

  say "Jordan's totient J_$k($n) is ", jordan_totient($k, $n);

Returns Jordan's totient function for a given integer value.  Jordan's totient
is a generalization of Euler's totient, where
  C<jordan_totient(1,$n) == euler_totient($n)>
This counts the number of k-tuples less than or equal to n that form a coprime
tuple with n.  As with C<euler_phi>, 0 is returned for all C<n E<lt> 1>.
This function can be used to generate some other useful functions, such as
the Dedekind psi function, where C<psi(n) = J(2,n) / J(1,n)>.


=head2 exp_mangoldt

  say "exp(lambda($_)) = ", exp_mangoldt($_) for 1 .. 100;

Returns EXP(Λ(n)), the exponential of the Mangoldt function (also known
as von Mangoldt's function) for an integer value.
The Mangoldt function is equal to log p if n is prime or a power of a prime,
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
This is the sum of the logarithm of each prime where C<p E<lt>= n>.  This
is effectively:

  my $s = 0;  forprimes { $s += log($_) } $n;  return $s;

=head2 chebyshev_psi

  say chebyshev_psi(10000);

Returns ψ(n), the second Chebyshev function for a non-negative integer input.
This is the sum of the logarithm of each prime power where C<p^k E<lt>= n>
for an integer k.  An alternate but slower computation is as the summatory
Mangoldt function, such as:

  my $s = 0;  for (1..$n) { $s += log(exp_mangoldt($_)) }  return $s;


=head2 divisor_sum

  say "Sum of divisors of $n:", divisor_sum( $n );
  say "sigma_2($n) = ", divisor_sum($n, 2);
  say "Number of divisors: sigma_0($n) = ", divisor_sum($n, 0);

This function takes a positive integer as input and returns the sum of
its divisors, including 1 and itself.  An optional second argument C<k>
may be given, which will result in the sum of the C<k-th> powers of the
divisors to be returned.

This is known as the sigma function (see Hardy and Wright section 16.7,
or OEIS A000203).  The API is identical to Pari/GP's C<sigma> function.
This function is useful for calculating things like aliquot sums, abundant
numbers, perfect numbers, etc.

The second argument may also be a code reference, which is called for each
divisor and the results are summed.  This allows computation of other
functions, but will be less efficient than using the numeric second argument.
This corresponds to Pari/GP's C<sumdiv> function.

An example of the 5th Jordan totient (OEIS A059378):

  divisor_sum( $n, sub { my $d=shift; $d**5 * moebius($n/$d); } );

though we have a function L</jordan_totient> which is more efficient.

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
formula using multi-precision floating point.  In 10 seconds:

           70    Integer::Partition
           90    MPU forpart { $n++ }
       10_000    MPU pure Perl partitions
      250_000    MPU GMP partitions
   35_000_000    Pari's numbpart
  500_000_000    Jonathan Bober's partitions_c.cc v0.6

If you want the enumerated partitions, see L</forpart>.


=head2 carmichael_lambda

Returns the Carmichael function (also called the reduced totient function,
or Carmichael λ(n)) of a positive integer argument.  It is the smallest
positive integer C<m> such that C<a^m = 1 mod n> for every integer C<a>
coprime to C<n>.  This is L<OEIS series A002322|http://oeis.org/A002322>.

=head2 kronecker

Returns the Kronecker symbol C<(a|n)> for two integers.  The possible
return values with their meanings for odd prime C<n> are:

   0   a = 0 mod n
   1   a is a quadratic residue mod n       (a = x^2 mod n for some x)
  -1   a is a quadratic non-residue mod n   (no a where a = x^2 mod n)

The Kronecker symbol is an extension of the Jacobi symbol to all integer
values of C<n> from the latter's domain of positive odd values of C<n>.
The Jacobi symbol is itself an extension of the Legendre symbol, which is
only defined for odd prime values of C<n>.  This corresponds to Pari's
C<kronecker(a,n)> function, Mathematica's C<KroneckerSymbol[n,m]>
function, and GMP's C<mpz_kronecker(a,n)>, C<mpz_jacobi(a,n)>, and
C<mpz_legendre(a,n)> functions.


=head2 binomial

Given integer arguments C<n> and C<k>, returns the binomial coefficient
C<n*(n-1)*...*(n-k+1)/k!>, also known as the choose function.  Negative
arguments use the L<Kronenburg extensions|http://arxiv.org/abs/1105.3689/>.
This corresponds to Pari's C<binomial(n,k)> function, Mathematica's
C<Binomial[n,k]> function, and GMP's C<mpz_bin_ui> function.

For negative arguments, this matches Mathematica.  Pari does not implement
the C<n E<lt> 0, k E<lt>= n> extension and instead returns C<0> for this
case.  GMP's API does not allow negative C<k> but otherwise matches.
L<Math::BigInt> does not implement any extensions and the results for
C<n E<lt> 0, k > 0> are undefined.


=head2 znorder

  $order = znorder(2, next_prime(10**16)-6);

Given two positive integers C<a> and C<n>, returns the multiplicative order
of C<a> modulo C<n>.  This is the smallest positive integer C<k> such that
C<a^k ≡ 1 mod n>.  Returns 1 if C<a = 1>.  Returns undef if C<a = 0> or if
C<a> and C<n> are not coprime, since no value will result in 1 mod n.
This corresponds to Pari's C<znorder(Mod(a,n))> function and Mathematica's
C<MultiplicativeOrder[a,n]> function.

=head2 znprimroot

Given a positive integer C<n>, returns the smallest primitive root
of C<(Z/nZ)^*>, or C<undef> if no root exists.  A root exists when
C<euler_phi($n) == carmichael_lambda($n)>, which will be true for
all prime C<n> and some composites.

L<OEIS A033948|http://oeis.org/A033948> is a sequence of integers where
the primitive root exists, while L<OEIS A046145|http://oeis.org/A046145>
is a list of the smallest primitive roots, which is what this function
produces.

=head2 znlog

  $k = znlog($a, $g, $p)

Returns the integer C<k> that solves the equation C<a = g^k mod p>, or
undef if no solution is found.  This is the discrete logarithm problem.
The implementation in this version is not very useful, but may be improved.


=head2 legendre_phi

  $phi = legendre_phi(1000000000, 41);

Given a non-negative integer C<n> and a non-negative prime number C<a>,
returns the Legendre phi function (also called Legendre's sum).  This is
the count of positive integers E<lt>= C<n> which are not divisible by any
of the first C<a> primes.


=head1 RANDOM PRIMES

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
primes are not evenly distributed.  Stated differently, the random prime
functions return a uniformly selected prime from the set of primes within
the range.  Hence given C<random_prime(1000)>, the numbers 2, 3, 487, 631,
and 997 all have the same probability of being returned.

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

  # System rand.  You probably don't want to do this.
  prime_set_config(irand => sub { int(rand(4294967296)) });

  # Math::Random::MTwist.  Fastest RNG by quite a bit.
  use Math::Random::MTwist;
  prime_set_config(irand => \&Math::Random::MTwist::_irand32);

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

  # Crypt::Random.  Uses Pari and /dev/random.  *VERY* slow.
  use Crypt::Random qw/makerandom/;
  prime_set_config(irand => sub { makerandom(Size=>32, Uniform=>1); });

  # Net::Random.  You probably don't want to use this, but if you do:
  use Net::Random;
  { my $rng = Net::Random->new(src=>"fourmilab.ch",max=>0xFFFFFFFF);
    sub nr_irand { return $rng->get(1); } }
  prime_set_config(irand => \&nr_irand);

  # Go back to MPU's default configuration
  prime_set_config(irand => undef);


=head2 random_ndigit_prime

  say "My 4-digit prime number is: ", random_ndigit_prime(4);

Selects a random n-digit prime, where the input is an integer number of
digits.  One of the primes within that range (e.g. 1000 - 9999 for
4-digits) will be uniformly selected using the C<irand> function as
described above.

If the number of digits is greater than or equal to the maximum native type,
then the result will be returned as a BigInt.  However, if the C<nobigint>
configuration option is on, then output will be restricted to native size
numbers, and requests for more digits than natively supported will result
in an error.
For better performance with large bit sizes, install L<Math::Prime::Util::GMP>.


=head2 random_nbit_prime

  my $bigprime = random_nbit_prime(512);

Selects a random n-bit prime, where the input is an integer number of bits.
A prime with the nth bit set will be uniformly selected, with randomness
supplied via calls to the C<irand> function as described above.

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
with large bit sizes, install L<Math::Prime::Util::GMP>.  Also see
L</random_shawe_taylor_prime>.

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


=head2 random_shawe_taylor_prime

  my $bigprime = random_shawe_taylor_prime(8192);

Construct an n-bit provable prime, using the Shawe-Taylor algorithm in
section C.6 of FIPS 186-4.  This uses 512 bits of randomness and SHA-256
as the hash.  This is a slightly simpler and older (1986) method than
Maurer's 1999 construction.  It is a bit faster than Maurer's method, and
uses less system entropy for large sizes.  The primary reason to use this
rather than Maurer's method is to use the FIPS 186-4 algorithm.

Similar to L</random_nbit_prime>, the result will be a BigInt if the
number of bits is greater than the native bit size.  For better performance
with large bit sizes, install L<Math::Prime::Util::GMP>.  Also see
L</random_maurer_prime> and L</random_proven_prime>.

Internally this additionally runs the BPSW probable prime test on every
partial result, and constructs a primality certificate for the final
result, which is verified.  These provide additional checks that the resulting
value has been properly constructed.


=head2 random_shawe_taylor_prime_with_cert

  my($n, $cert) = random_shawe_taylor_prime_with_cert(4096)

As with L</random_shawe_taylor_prime>, but returns a two-element array
containing the n-bit provable prime along with a primality certificate.
The certificate is the same as produced by L</prime_certificate> or
L</is_provable_prime_with_cert>, and can be parsed by L</verify_prime> or
any other software that understands MPU primality certificates.
The proof construction consists of a single chain of C<Pocklington> types.



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

  verbose         verbose level.  1 or more will result in extra output.
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

  verbose      The default setting of 0 will generate no extra output.
               Setting to 1 or higher results in extra output.  For
               example, at setting 1 the AKS algorithm will indicate
               the chosen r and s values.  At setting 2 it will output
               a sequence of dots indicating progress.  Similarly, for
               random_maurer_prime, setting 3 shows real time progress.
               Factoring large numbers is another place where verbose
               settings can give progress indications.

  xs           Allows turning off the XS code, forcing the Pure Perl
               code to be used.  Set to 0 to disable XS, set to 1 to
               re-enable.  You probably will never want to do this.

  gmp          Allows turning off the use of L<Math::Prime::Util::GMP>,
               which means using Pure Perl code for big numbers.  Set
               to 0 to disable GMP, set to 1 to re-enable.
               You probably will never want to do this.

  assume_rh    Allows functions to assume the Riemann hypothesis is
               true if set to 1.  This defaults to 0.  Currently this
               setting only impacts prime count lower and upper
               bounds, but could later be applied to other areas such
               as primality testing.  A later version may also have a
               way to indicate whether no RH, RH, GRH, or ERH is to
               be assumed.

  irand        Takes a code ref to an irand function returning a
               uniform number between 0 and 2**32-1.  This will be
               used for all random number generation in the module.


=head1 FACTORING FUNCTIONS

=head2 factor

  my @factors = factor(3_369_738_766_071_892_021);
  # returns (204518747,16476429743)

Produces the prime factors of a positive number input, in numerical order.
The product of the returned factors will be equal to the input.  C<n = 1>
will return an empty list, and C<n = 0> will return 0.  This matches Pari.

In scalar context, returns Ω(n), the total number of prime factors
(L<OEIS A001222|http://oeis.org/A001222>).
This corresponds to Pari's C<bigomega(n)> function and Mathematica's
C<PrimeOmega[n]> function.
This is same result that we would get if we evaluated the resulting
array in scalar context.

The current algorithm does a little trial division, a check for perfect
powers, followed by combinations of Pollard's Rho, SQUFOF, and Pollard's
p-1.  The combination is applied to each non-prime factor found.

Factoring bigints works with pure Perl, and can be very handy on 32-bit
machines for numbers just over the 32-bit limit, but it can be B<very> slow
for "hard" numbers.  Installing the L<Math::Prime::Util::GMP> module will
speed up bigint factoring a B<lot>, and all future effort on large number
factoring will be in that module.  If you do not have that module for
some reason, use the GMP or Pari version of bigint if possible
(e.g. C<use bigint try =E<gt> 'GMP,Pari'>), which will run 2-3x faster
(though still 100x slower than the real GMP code).


=head2 factor_exp

  my @factor_exponent_pairs = factor_exp(29513484000);
  # returns ([2,5], [3,4], [5,3], [7,2], [11,1], [13,2])
  # factor(29513484000)
  # returns (2,2,2,2,2,3,3,3,3,5,5,5,7,7,11,13,13)

Produces pairs of prime factors and exponents in numerical factor order.
This is more convenient for some algorithms.  This is the same form that
Mathematica's C<FactorInteger[n]> and Pari/GP's C<factorint> functions
return.  Note that L<Math::Pari> transposes the Pari result matrix.

In scalar context, returns ω(n), the number of unique prime factors
(L<OEIS A001221|http://oeis.org/A001221>).
This corresponds to Pari's C<omega(n)> function and Mathematica's
C<PrimeNu[n]> function.
This is same result that we would get if we evaluated the resulting
array in scalar context.

The internals are identical to L</factor>, so all comments there apply.
Just the way the factors are arranged is different.


=head2 divisors

  my @divisors = divisors(30);   # returns (1, 2, 3, 5, 6, 10, 15, 30)

Produces all the divisors of a positive number input, including 1 and
the input number.  The divisors are a power set of multiplications of
the prime factors, returned as a uniqued sorted list.  The result is
identical to that of Pari's C<divisors> and Mathematica's C<Divisors[n]>
functions.

In scalar context this returns the sigma0 function,
the sigma function (see Hardy and Wright section 16.7, or OEIS A000203).
This is the same result as evaluating the array in scalar context.

Also see the L</for_divisors> functions for looping over the divisors.

C<all_factors> is the deprecated name for this function.


=head2 trial_factor

  my @factors = trial_factor($n);

Produces the prime factors of a positive number input.
The factors will be in numerical order.
For large inputs this will be very slow.
Like all the specific-algorithm C<*_factor> routines, this is not exported
unless explicitly requested.

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
and there are some inputs it can rapidly factor.  Overall it has the
same advantages and disadvantages as Fermat's method.

=head2 squfof_factor

  my @factors = squfof_factor($n);

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

    # Similar to A001262's isStrongPsp function, but much faster
    perl -MMath::Prime::Util=:all -E 'forcomposites { say if is_strong_pseudoprime($_,17) } 10000000;'

Print some primes above 64-bit range:

    perl -MMath::Prime::Util=:all -Mbigint -E 'my $start=100000000000000000000; say join "\n", @{primes($start,$start+1000)}'

    # Another way
    perl -MMath::Prime::Util=:all -E 'forprimes { say } "100000000000000000039", "100000000000000000993"'

    # Similar using Math::Pari:
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
  print $sgit->(), "\n"  for 1 .. 10000;

Project Euler, problem 3 (Largest prime factor):

  use Math::Prime::Util qw/factor/;
  use bigint;  # Only necessary for 32-bit machines.
  say 0+(factor(600851475143))[-1]

Project Euler, problem 7 (10001st prime):

  use Math::Prime::Util qw/nth_prime/;
  say nth_prime(10_001);

Project Euler, problem 10 (summation of primes):

  use Math::Prime::Util qw/forprimes/;
  my $sum = 0;
  forprimes { $sum += $_ } 2_000_000;
  say $sum;
  # Or (fine for this size, not good for very large values)
  use Math::Prime::Util qw/vecsum primes/;
  say vecsum( @{primes(2_000_000)} );

Project Euler, problem 21 (Amicable numbers):

  use Math::Prime::Util qw/divisor_sum/;
  my $sum = 0;
  foreach my $x (1..10000) {
    my $y = divisor_sum($x)-$x;
    $sum += $x + $y if $y > $x && $x == divisor_sum($y)-$y;
  }
  say $sum;
  # Or using a pipeline:
  use Math::Prime::Util qw/vecsum divisor_sum/;
  say vecsum( map { divisor_sum($_) }
              grep { my $y = divisor_sum($_)-$_;
                     $y > $_ && $_==(divisor_sum($y)-$y) }
              1 .. 10000 );

Project Euler, problem 41 (Pandigital prime), brute force command line:

  perl -MMath::Prime::Util=primes -MList::Util=first -E 'say first { /1/&&/2/&&/3/&&/4/&&/5/&&/6/&&/7/} reverse @{primes(1000000,9999999)};'

Project Euler, problem 47 (Distinct primes factors):

  use Math::Prime::Util qw/pn_primorial factor_exp/;
  my $n = pn_primorial(4);  # Start with the first 4-factor number
  # factor_exp in scalar context returns the number of distinct prime factors
  $n++ while (factor_exp($n) != 4 || factor_exp($n+1) != 4 || factor_exp($n+2) != 4 || factor_exp($n+3) != 4);
  say $n;

Project Euler, problem 69, stupid brute force solution (about 1 second):

  use Math::Prime::Util qw/euler_phi/;
  my ($maxn, $maxratio) = (0,0);
  foreach my $n (1..1000000) {
    my $ndivphi = $n / euler_phi($n);
    ($maxn, $maxratio) = ($n, $ndivphi) if $ndivphi > $maxratio;
  }
  say "$maxn  $maxratio";

Here is the right way to do PE problem 69 (under 0.03s):

  use Math::Prime::Util qw/pn_primorial/;
  my $n = 0;
  $n++ while pn_primorial($n+1) < 1000000;
  say pn_primorial($n);

Project Euler, problem 187, stupid brute force solution, 2 to 3 minutes:

  use Math::Prime::Util qw/forcomposites factor/;
  my $nsemis = 0;
  forcomposites { $nsemis++ if scalar factor($_) == 2; } int(10**8)-1;
  say $nsemis;

Here is one of the best ways for PE187:  under 20 milliseconds from the
command line.  Much faster than Pari, and competitive with Mathematica.

  use Math::Prime::Util qw/forprimes prime_count/;
  my $limit = shift || int(10**8);
  $limit--;
  my ($sum, $pc) = (0, 1);
  forprimes {
    $sum += prime_count(int($limit/$_)) + 1 - $pc++;
  } int(sqrt($limit));
  say $sum;

To get the result of L<Math::Factor::XS/matches>:

  use Math::Prime::Util qw/divisors/;
  sub matches {
    my @d = divisors(shift);
    return map { [$d[$_],$d[$#d-$_]] } 1..(@d-1)>>1;
  }
  my $n = 139650;
  say "$n = ", join(" = ", map { "$_->[0]·$_->[1]" } matches($n));

or its C<matches> function with the C<skip_multiples> option:

  sub matches {
    my @d = divisors(shift);
    return map { [$d[$_],$d[$#d-$_]] }
           grep { my $div=$d[$_]; !scalar(grep {!($div % $d[$_])} 1..$_-1) }
           1..(@d-1)>>1; }
  }

Compute L<OEIS A054903|http://oeis.org/A054903> just like CRG4's Pari example:

  use Math::Prime::Util qw/forcomposite divisor_sum/;
  forcomposites {
    say if divisor_sum($_)+6 == divisor_sum($_+6)
  } 9,1e7;

Construct the table shown in L<OEIS A046147|http://oeis.org/A046147>:

  use Math::Prime::Util qw/znorder euler_phi gcd/;
  foreach my $n (1..100) {
    if (!znprimroot($n)) {
      say "$n -";
    } else {
      my $phi = euler_phi($n);
      my @r = grep { gcd($_,$n) == 1 && znorder($_,$n) == $phi } 1..$n-1;
      say "$n ", join(" ", @r);
    }
  }

=head1 PRIMALITY TESTING NOTES

Above C<2^64>, L</is_prob_prime> performs an extra-strong
L<BPSW test|http://en.wikipedia.org/wiki/Baillie-PSW_primality_test>
which is fast (a little less than the time to perform 3 Miller-Rabin
tests) and has no known counterexamples.  If you trust the primality
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

Math systems 30 years ago typically used Miller-Rabin tests with C<k>
bases (usually fixed bases, sometimes random) for primality
testing, but these have generally been replaced by some form of BPSW
as used in this module.  See Pinch's 1993 paper for examples of why
using C<k> M-R tests leads to poor results.  The three exceptions in
common contemporary use I am aware of are:

=over 4

=item libtommath

Uses the first C<k> prime bases.  This is problematic for
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

=item L<Math::Pari> (not recent Pari/GP)

Pari 2.1.7 is the default version installed with the L<Math::Pari>
module.  It uses 10 random M-R bases (the PRNG uses a fixed seed
set at compile time).  Pari 2.3.0 was released in May 2006 and it,
like all later releases through at least 2.6.1, use BPSW / APRCL,
after complaints of false results from using M-R tests.  For example,
it will indicate 9 is prime about 1 out of every 276k calls.

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
everyone else's guesses and results, and can keep playing as long as
they like.  It's only valid if the players are completely oblivious to
what is happening.


=head1 LIMITATIONS

Perl versions earlier than 5.8.0 have problems doing exact integer math.
Some operations will flip signs, and many operations will convert intermediate
or output results to doubles, which loses precision on 64-bit systems.
This causes numerous functions to not work properly.  The test suite will
try to determine if your Perl is broken (this only applies to really old
versions of Perl compiled for 64-bit when using numbers larger than
C<~ 2^49>).  The best solution is updating to a more recent Perl.

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
a module doesn't match what I believe are the best set of features doesn't
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
large 64-bit numbers.  MPU is similarly fast with small numbers, but becomes
faster as the size increases.
MPXS's prime sieve is an unoptimized non-segmented SoE
which returns an array.  Sieve bases larger than C<10^7> start taking
inordinately long and using a lot of memory (gigabytes beyond C<10^10>).
E.g. C<primes(10**9, 10**9+1000)> takes 36 seconds with MPXS, but only
0.0001 seconds with MPU.

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
L<RT81858|https://rt.cpan.org/Ticket/Display.html?id=81858>.  MPU has
a more uniform distribution as well as return a larger subset of primes
(L<RT81871|https://rt.cpan.org/Ticket/Display.html?id=81871>).
MPU does not depend on L<Math::Pari> though can run slow for bigints unless
the L<Math::BigInt::GMP> or L<Math::BigInt::Pari> modules are installed.
Having L<Math::Prime::Util::GMP> installed also helps performance for MPU.
Crypt::Primes is hardcoded to use L<Crypt::Random>, while MPU uses
L<Bytes::Random::Secure>, and also allows plugging in a random function.
This is more flexible, faster, has fewer dependencies, and uses a CSPRNG
for security.  MPU can return a primality certificate.
What Crypt::Primes has that MPU does not is the ability to return a generator.

L<Math::Factor::XS> calculates prime factors and factors, which correspond to
the L</factor> and L</divisors> functions of MPU.  These functions do
not support bigints.  Both are implemented with trial division, meaning they
are very fast for really small values, but become very slow as the input
gets larger (factoring 19 digit semiprimes is over 1000 times slower).  The
function C<count_prime_factors> can be done in MPU using C<scalar factor($n)>.
See the L</"EXAMPLES"> section for a 2-line function replicating C<matches>.

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
primality functionality.  It was the first CPAN module to support the BPSW
test.  All inputs are processed using GMP, so it of course supports
bigints.  In fact, Math::Primality was made originally with bigints in mind,
while MPU was originally targeted to native integers, but both have added
better support for the other.  The main differences are extra functionality
(MPU has more functions) and performance.  With native integer inputs, MPU
is generally much faster, especially with L</prime_count>.  For bigints,
MPU is slower unless the L<Math::Prime::Util::GMP> module is installed, in
which case MPU is ~2x faster.  L<Math::Primality> also installs
a C<primes.pl> program, but it has much less functionality than the one
included with MPU.

L<Math::NumSeq> does not have a one-to-one mapping between functions in MPU,
but it does offer a way to get many similar results such as
primes, twin primes, Sophie-Germain primes, lucky primes, moebius, divisor
count, factor count, Euler totient, primorials, etc.  Math::NumSeq is
set up for accessing these values in order rather than for arbitrary values,
though a few sequences support random access.  The primary advantage I see
is the uniform access mechanism for a I<lot> of sequences.  For those methods
that overlap, MPU is usually much faster.  Importantly, most of the sequences
in Math::NumSeq are limited to 32-bit indices.

L<Math::Pari> supports a lot of features, with a great deal of overlap.  In
general, MPU will be faster for native 64-bit integers, while it's differs
for bigints (Pari will always be faster if L<Math::Prime::Util::GMP> is not
installed; with it, it varies by function).  Note that Pari extends many of
these functions to other spaces (Gaussian integers, complex numbers, vectors,
matrices, polynomials, etc.) which are beyond the realm of this module.
Some of the highlights:

=over 4

=item C<isprime>

The default L<Math::Pari> is built with Pari 2.1.7.  This uses 10 M-R
tests with randomly chosen bases (fixed seed, but doesn't reset each
invocation like GMP's C<is_probab_prime>).  This has a greater chance
of false positives compared to the BPSW test -- some composites such as
C<9>, C<88831>, C<38503>, etc.
(L<OEIS A141768|http://oeis.org/A141768>)
have a surprisingly high chance of being indicated prime.
Using C<isprime($n,1)> will perform an C<n-1> proof,
but this becomes unreasonably slow past 70 or so digits.

If L<Math::Pari> is built using Pari 2.3.5 (this requires manual
configuration) then the primality tests are completely different.  Using
C<ispseudoprime> will perform a BPSW test and is quite a bit faster than
the older test.  C<isprime> now does an APR-CL proof (fast, but no
certificate).

L<Math::Primality> uses a strong BPSW test, which is the standard BPSW
test based on the 1980 paper.  It has no known counterexamples (though
like all these tests, we know some exist).  Pari 2.3.5 (and through at
least 2.6.2) uses an almost-extra-strong BPSW test for its
C<ispseudoprime> function.  This is deterministic for native integers,
and should be excellent for bigints, with a slightly lower chance of
counterexamples than the traditional strong test.
L<Math::Prime::Util> uses the
full extra-strong BPSW test, which has an even lower chance of
counterexample.
With L<Math::Prime::Util::GMP>, C<is_prime> adds 1 to 5 extra M-R tests
using random bases, which further reduces the probability of a composite
being allowed to pass.

=item C<primepi>

Only available with version 2.3 of Pari.  Similar to MPU's L</prime_count>
function in API, but uses a naive counting algorithm with its precalculated
primes, so is not of practical use.  Incidently, Pari 2.6 (not usable from
Perl) has fixed the pre-calculation requirement so it is more useful, but is
still thousands of times slower than MPU.

=item C<primes>

Doesn't support ranges, requires bumping up the precalculated
primes for larger numbers, which means knowing in advance the upper limit
for primes.  Support for numbers larger than 400M requires using Pari
version 2.3.5.  If that is used, sieving is about 2x faster than MPU, but
doesn't support segmenting.

=item C<factorint>

Similar to MPU's L</factor_exp> though with a slightly different return.
MPU offers L</factor> for a linear array of prime factors where
   n = p1 * p2 * p3 * ...   as (p1,p2,p3,...)
and L</factor_exp> for an array of factor/exponent pairs where:
   n = p1^e1 * p2^e2 * ...  as ([p1,e1],[p2,e2],...)
Pari/GP returns an array similar to the latter.  L<Math::Pari> returns
a transposed matrix like:
   n = p1^e1 * p2^e2 * ...  as ([p1,p2,...],[e1,e2,...])
Slower than MPU for all 64-bit inputs on an x86_64 platform, it may be
faster for large values on other platforms.  With the newer
L<Math::Prime::Util::GMP> releases, bigint factoring is slightly
faster on average in MPU.

=item C<divisors>

Similar to MPU's L</divisors>.

=item C<forprime>, C<forcomposite>, C<fordiv>, C<sumdiv>

Similar to MPU's L</forprimes>, L</forcomposites>, L</fordivisors>, and
L</divisor_sum>.

=item C<eulerphi>, C<moebius>

Similar to MPU's L</euler_phi> and L</moebius>.  MPU is 2-20x faster for
native integers.  MPU also supported range inputs, which can be much
more efficient.  Without L<Math::Prime::Util::GMP> installed, MPU is
very slow with bigints.  With it installed, it is about 2x slower than
Math::Pari.

=item C<gcd>, C<lcm>, C<kronecker>, C<znorder>, C<znprimroot>, C<znlog>

Similar to MPU's L</gcd>, L</lcm>, L</kronecker>, L</znorder>,
L</znprimroot>, and L</znlog>.  Pari's C<znprimroot> only returns the
smallest root for prime powers.  The behavior is undefined when the group is
not cyclic (sometimes it throws an exception, sometimes it returns
an incorrect answer, sometimes it hangs).  MPU's L</znprimroot> will always
return the smallest root if it exists, and C<undef> otherwise.

=item C<sigma>

Similar to MPU's L</divisor_sum>.  MPU is ~10x faster for native integers
and about 2x slower for bigints.

=item C<numbpart>, C<forpart>

Similar to MPU's L</partitions> and L</forpart>.  These functions were
introduced in Pari 2.3 and 2.6, hence are not in Math::Pari.  C<numbpart>
produce identical results to C<partitions>, but Pari is I<much> faster.
L<forpart> is very similar to Pari's function, but produces a different
ordering (MPU is the standard anti-lexicographical, Pari uses a size sort).
Currently Pari is somewhat faster due to Perl function call overhead.  When
using restrictions, Pari has much better optimizations.

=item C<eint1>

Similar to MPU's L</ExponentialIntegral>.

=item C<zeta>

MPU has L</RiemannZeta> which takes non-negative real inputs, while Pari's
function supports negative and complex inputs.

=back

Overall, L<Math::Pari> supports a huge variety of functionality and has a
sophisticated and mature code base behind it (noting that the Pari library
used is about 10 years old now).
For native integers, typically Math::Pari will be slower than MPU.  For
bigints, Math::Pari may be superior and it rarely has any performance
surprises.  Some of the
unique features MPU offers include super fast prime counts, nth_prime,
ECPP primality proofs with certificates, approximations and limits for both,
random primes, fast Mertens calculations, Chebyshev theta and psi functions,
and the logarithmic integral and Riemann R functions.  All with fairly
minimal installation requirements.


=head1 PERFORMANCE

First, for those looking for the state of the art non-Perl solutions:

=over 4

=item Primality testing

For general numbers smaller than 2000 or so digits, MPU is the fastest
solution I am aware of (it is faster than Pari 2.7, PFGW, and FLINT).
For very large inputs,
L<PFGW|http://sourceforge.net/projects/openpfgw/> is the fastest primality
testing software I'm aware of.  It has fast trial division, and is especially
fast on many special forms.  It does not have a BPSW test however, and there
are quite a few counterexamples for a given base of its PRP test, so it is
commonly used for fast filtering of large candidates.
A test such as the BPSW test in this module is then recommended.

=item Primality proofs

L<Primo|http://www.ellipsa.eu/> is the best method for open source primality
proving for inputs over 1000 digits.  Primo also does well below that size,
but other good alternatives are
David Cleaver's L<mpzaprcl|http://sourceforge.net/projects/mpzaprcl/>,
the APRCL from the modern L<Pari|http://pari.math.u-bordeaux.fr/> package,
or the standalone ECPP from this module with large polynomial set.

=item Factoring

L<yafu|http://sourceforge.net/projects/yafu/>,
L<msieve|http://sourceforge.net/projects/msieve/>, and
L<gmp-ecm|http://ecm.gforge.inria.fr/> are all good choices for large
inputs.  The factoring code in this module (and all other CPAN modules) is
very limited compared to those.

=item Primes

L<primesieve|http://code.google.com/p/primesieve/> and
L<yafu|http://sourceforge.net/projects/yafu/>
are the fastest publically available code I am aware of.  Primesieve
will additionally take advantage of multiple cores with excellent
efficiency.
Tomás Oliveira e Silva's private code may be faster for very large
values, but isn't available for testing.

Note that the Sieve of Atkin is I<not> faster than the Sieve of Eratosthenes
when both are well implemented.  The only Sieve of Atkin that is even
competitive is Bernstein's super optimized I<primegen>, which runs on par
with the SoE in this module.  The SoE's in Pari, yafu, and primesieve
are all faster.

=item Prime Counts and Nth Prime

Outside of private research implementations doing prime counts for
C<n E<gt> 2^64>, this module should be close to state of the art in
performance, and supports results up to C<2^64>.  Further performance
improvements are planned, as well as expansion to larger values.

The fastest solution for small inputs is a hybrid table/sieve method.
This module does this for values below 60M.  As the inputs get larger,
either the tables have to grow exponentially or speed must be
sacrificed.  Hence this is not a good general solution for most uses.

=back


=head2 PRIME COUNTS

Counting the primes to C<800_000_000> (800 million):

  Time (s)   Module                      Version  Notes
  ---------  --------------------------  -------  -----------
       0.001 Math::Prime::Util           0.37     using extended LMO
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
  >20000     Math::Big                   1.12     Perl, > 26GB RAM used

Python's standard modules are very slow: MPMATH v0.17 C<primepi> takes 169.5s
and 25+ GB of RAM.  SymPy 0.7.1 C<primepi> takes 292.2s.  However there are
very fast solutions written by Robert William Hanks (included in the xt/
directory of this distribution): pure Python in 12.1s and NUMPY in 2.8s.

=head2 PRIMALITY TESTING

=over 4

=item Small inputs:  is_prime from 1 to 20M

    2.6s  Math::Prime::Util      (sieve lookup if prime_precalc used)
    3.4s  Math::Prime::FastSieve (sieve lookup)
    4.4s  Math::Prime::Util      (trial + deterministic M-R)
   10.9s  Math::Prime::XS        (trial)
   36.5s  Math::Pari w/2.3.5     (BPSW)
   78.2s  Math::Pari             (10 random M-R)
  501.3s  Math::Primality        (deterministic M-R)

=item Large native inputs:  is_prime from 10^16 to 10^16 + 20M

    7.0s  Math::Prime::Util      (BPSW)
   42.6s  Math::Pari w/2.3.5     (BPSW)
  144.3s  Math::Pari             (10 random M-R)
  664.0s  Math::Primality        (BPSW)
  30 HRS  Math::Prime::XS        (trial)

  These inputs are too large for Math::Prime::FastSieve.

=item bigints:  is_prime from 10^100 to 10^100 + 0.2M

    2.5s  Math::Prime::Util          (BPSW + 1 random M-R)
    3.0s  Math::Pari w/2.3.5         (BPSW)
   12.9s  Math::Primality            (BPSW)
   35.3s  Math::Pari                 (10 random M-R)
   53.5s  Math::Prime::Util w/o GMP  (BPSW)
   94.4s  Math::Prime::Util          (n-1 or ECPP proof)
  102.7s  Math::Pari w/2.3.5         (APR-CL proof)

=back

=over 4

=item *

MPU is consistently the fastest solution, and performs the most
stringent probable prime tests on bigints.

=item *

Math::Primality has a lot of overhead that makes it quite slow for
native size integers.  With bigints we finally see it work well.

=item *

Math::Pari built with 2.3.5 not only has a better primality test versus
the default 2.1.7, but runs faster.  It still has quite a bit of overhead
with native size integers.  Pari/GP 2.5.0 takes 11.3s, 16.9s, and 2.9s
respectively for the tests above.  MPU is still faster, but clearly the
time for native integers is dominated by the calling overhead.

=back

=head2 FACTORING

Factoring performance depends on the input, and the algorithm choices used
are still being tuned.  L<Math::Factor::XS> is very fast when given input with
only small factors, but it slows down rapidly as the smallest factor increases
in size.  For numbers larger than 32 bits, L<Math::Prime::Util> can be 100x or
more faster (a number with only very small factors will be nearly identical,
while a semiprime may be 3000x faster).  L<Math::Pari>
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
the version included in Pari.  Both are pretty fast to about 60 digits, and
work reasonably well to 80 or so before starting to take many minutes per
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

  bits    random   +testing  rand_prov   Maurer   Shw-Tylr  CPMaurer
  -----  --------  --------  ---------  --------  --------  --------
     64    0.0001  +0.000008   0.0002     0.0001    0.010     0.022
    128    0.0020  +0.00023    0.011      0.063     0.028     0.057
    256    0.0034  +0.0004     0.058      0.13      0.042     0.16
    512    0.0097  +0.0012     0.28       0.28      0.085     0.41
   1024    0.060   +0.0060     0.65       0.65      0.24      2.19
   2048    0.57    +0.039      4.8        4.8       1.0      10.99
   4096    6.24    +0.25      31.9       31.9       8.2      79.71
   8192   58.6     +1.61     234.0      234.0     112.9     947.3

  random    = random_nbit_prime  (results pass BPSW)
  random+   = additional time for 3 M-R and a Frobenius test
  rand_prov = random_proven_prime
  maurer    = random_maurer_prime
  Shw-Tylr  = random_shawe_taylor_prime
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

L</random_shawe_taylor_prime> similarly constructs a provable prime.  It
uses a simpler construction method.  The implementation uses a single large
random seed followed by SHA-256 as specified by FIPS 186-4.  As seen, it
is a bit faster than the Maurer implementation.

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
are planned for a later version.

The LMO implementation is based on the 2003 preprint from Christian Bau,
as well as the 2006 paper from Tomás Oliveira e Silva.  I also want to
thank Kim Walisch for the many discussions about prime counting.


=head1 REFERENCES

=over 4

=item *

Christian Bau, "The Extended Meissel-Lehmer Algorithm", 2003, preprint with example C++ implementation.  Very detailed implementation-specific paper which was used for the implementation here.  Highly recommended for implementing a sieve-based LMO.  L<http://cs.swan.ac.uk/~csoliver/ok-sat-library/OKplatform/ExternalSources/sources/NumberTheory/ChristianBau/>

=item *

Manuel Benito and Juan L. Varona, "Recursive formulas related to the summation of the Möbius function", I<The Open Mathematics Journal>, v1, pp 25-34, 2007.  Among many other things, shows a simple formula for computing the Mertens functions with only n/3 Möbius values (not as fast as Deléglise and Rivat, but really simple).  L<http://www.unirioja.es/cu/jvarona/downloads/Benito-Varona-TOMATJ-Mertens.pdf>

=item *

John Brillhart, D. H. Lehmer, and J. L. Selfridge, "New Primality Criteria and Factorizations of 2^m +/- 1", Mathematics of Computation, v29, n130, Apr 1975, pp 620-647.  L<http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf>

=item *

W. J. Cody and Henry C. Thacher, Jr., "Rational Chebyshev Approximations for the Exponential Integral E_1(x)", I<Mathematics of Computation>, v22, pp 641-649, 1968.

=item *

W. J. Cody and Henry C. Thacher, Jr., "Chebyshev approximations for the exponential integral Ei(x)", I<Mathematics of Computation>, v23, pp 289-303, 1969.  L<http://www.ams.org/journals/mcom/1969-23-106/S0025-5718-1969-0242349-2/>

=item *

W. J. Cody, K. E. Hillstrom, and Henry C. Thacher Jr., "Chebyshev Approximations for the Riemann Zeta Function", L<Mathematics of Computation>, v25, n115, pp 537-547, July 1971.

=item *

Henri Cohen, "A Course in Computational Algebraic Number Theory", Springer, 1996.  Practical computational number theory from the team lead of L<Pari|http://pari.math.u-bordeaux.fr/>.  Lots of explicit algorithms.

=item *

Marc Deléglise and Joöl Rivat, "Computing the summation of the Möbius function", I<Experimental Mathematics>, v5, n4, pp 291-295, 1996.  Enhances the Möbius computation in Lioen/van de Lune, and gives a very efficient way to compute the Mertens function.  L<http://projecteuclid.org/euclid.em/1047565447>

=item *

Pierre Dusart, "Autour de la fonction qui compte le nombre de nombres premiers", PhD thesis, 1998.  In French.  The mathematics is readable and highly recommended reading if you're interesting in prime number bounds.  L<http://www.unilim.fr/laco/theses/1998/T1998_01.html>

=item *

Pierre Dusart, "Estimates of Some Functions Over Primes without R.H.", preprint, 2010.  Updates to the best non-RH bounds for prime count and nth prime.  L<http://arxiv.org/abs/1002.0442/>

=item *

Pierre-Alain Fouque and Mehdi Tibouchi, "Close to Uniform Prime Number Generation With Fewer Random Bits", pre-print, 2011.  Describes random prime distributions, their algorithm for creating random primes using few random bits, and comparisons to other methods.  Definitely worth reading for the discussions of uniformity.  L<http://eprint.iacr.org/2011/481>

=item *

Walter M. Lioen and Jan van de Lune, "Systematic Computations on Mertens' Conjecture and Dirichlet's Divisor Problem by Vectorized Sieving", in I<From Universal Morphisms to Megabytes>, Centrum voor Wiskunde en Informatica, pp. 421-432, 1994.  Describes a nice way to compute a range of Möbius values.  L<http://walter.lioen.com/papers/LL94.pdf>

=item *

Ueli M. Maurer, "Fast Generation of Prime Numbers and Secure Public-Key Cryptographic Parameters", 1995.  Generating random provable primes by building up the prime.  L<http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.26.2151>

=item *

Gabriel Mincu, "An Asymptotic Expansion", I<Journal of Inequalities in Pure and Applied Mathematics>, v4, n2, 2003.  A very readable account of Cipolla's 1902 nth prime approximation.  L<http://www.emis.de/journals/JIPAM/images/153_02_JIPAM/153_02.pdf>

=item *

L<OEIS: Primorial|http://oeis.org/wiki/Primorial>

=item *

Vincent Pegoraro and Philipp Slusallek, "On the Evaluation of the Complex-Valued Exponential Integral", I<Journal of Graphics, GPU, and Game Tools>, v15, n3, pp 183-198, 2011.  L<http://www.cs.utah.edu/~vpegorar/research/2011_JGT/paper.pdf>

=item *

William H. Press et al., "Numerical Recipes", 3rd edition.

=item *

Hans Riesel, "Prime Numbers and Computer Methods for Factorization", Birkh?user, 2nd edition, 1994.  Lots of information, some code, easy to follow.

=item *

David M. Smith, "Multiple-Precision Exponential Integral and Related Functions", I<ACM Transactions on Mathematical Software>, v37, n4, 2011.  L<http://myweb.lmu.edu/dmsmith/toms2011.pdf>

=item *

Douglas A. Stoll and Patrick Demichel , "The impact of ζ(s) complex zeros on π(x) for x E<lt> 10^{10^{13}}", L<Mathematics of Computation>, v80, n276, pp 2381-2394, October 2011.  L<http://www.ams.org/journals/mcom/2011-80-276/S0025-5718-2011-02477-4/home.html>

=back


=head1 COPYRIGHT

Copyright 2011-2014 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
