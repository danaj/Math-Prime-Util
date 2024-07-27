package Math::Prime::Util::PPFE;
use strict;
use warnings;
use Math::Prime::Util::PP;
use Math::Prime::Util::Entropy;

# The PP front end, only loaded if XS is not used.
# It is intended to load directly into the MPU namespace.

package Math::Prime::Util;
use Carp qw/carp croak confess/;

*_validate_integer = \&Math::Prime::Util::PP::_validate_integer;
*_validate_integer_nonneg = \&Math::Prime::Util::PP::_validate_integer_nonneg;
*_validate_integer_positive = \&Math::Prime::Util::PP::_validate_integer_positive;
*_validate_integer_abs = \&Math::Prime::Util::PP::_validate_integer_abs;

*_prime_memfreeall = \&Math::Prime::Util::PP::_prime_memfreeall;
*prime_memfree  = \&Math::Prime::Util::PP::prime_memfree;
*prime_precalc  = \&Math::Prime::Util::PP::prime_precalc;

use Math::Prime::Util::ChaCha;
*_is_csprng_well_seeded = \&Math::Prime::Util::ChaCha::_is_csprng_well_seeded;
*_csrand = \&Math::Prime::Util::ChaCha::csrand;
*_srand = \&Math::Prime::Util::ChaCha::srand;
*random_bytes = \&Math::Prime::Util::ChaCha::random_bytes;
*irand = \&Math::Prime::Util::ChaCha::irand;
*irand64 = \&Math::Prime::Util::ChaCha::irand64;

sub srand {
  my($seed) = @_;
  croak "secure option set, manual seeding disabled" if prime_get_config()->{'secure'};
  if (!defined $seed) {
    my $nbytes = (~0 == 4294967295) ? 4 : 8;
    $seed = entropy_bytes( $nbytes );
    $seed = unpack(($nbytes==4) ? "L" : "Q", $seed);
  }
  Math::Prime::Util::GMP::seed_csprng(8,pack("LL",$seed))
    if $Math::Prime::Util::_GMPfunc{"seed_csprng"};
  Math::Prime::Util::_srand($seed);
}
sub csrand {
  my($seed) = @_;
  croak "secure option set, manual seeding disabled" if defined $seed && prime_get_config()->{'secure'};
  $seed = entropy_bytes( 64 ) unless defined $seed;
  Math::Prime::Util::GMP::seed_csprng(length($seed),$seed)
    if $Math::Prime::Util::_GMPfunc{"seed_csprng"};
  Math::Prime::Util::_csrand($seed);
  1; # Don't return the seed
}
sub entropy_bytes {
  my($bytes) = @_;
  croak "entropy_bytes: input must be integer bytes between 1 and 4294967295"
    if !defined($bytes) || $bytes < 1 || $bytes > 4294967295 || $bytes != int($bytes);
  my $data = Math::Prime::Util::Entropy::entropy_bytes($bytes);
  if (!defined $data) {
    # We can't find any entropy source!  Highly unusual.
    Math::Prime::Util::_srand();
    $data = random_bytes($bytes);
  }
  croak "entropy_bytes internal got wrong amount!" unless length($data) == $bytes;
  $data;
}

# Fill all the mantissa bits for our NV, regardless of 32-bit or 64-bit Perl.
{
  use Config;
  my $nvbits = (defined $Config{nvmantbits})  ? $Config{nvmantbits}
             : (defined $Config{usequadmath}) ? 112
             : 53;
  my $uvbits = (~0 > 4294967295) ? 64 : 32;
  my $rsub;
  my $_tonv_32  = 1.0;        $_tonv_32 /= 2.0 for 1..32;
  my $_tonv_64  = $_tonv_32;  $_tonv_64 /= 2.0 for 1..32;
  my $_tonv_96  = $_tonv_64;  $_tonv_96 /= 2.0 for 1..32;
  my $_tonv_128 = $_tonv_96;  $_tonv_128/= 2.0 for 1..32;
  if ($uvbits == 64) {
    if ($nvbits <= 32) {
      *drand = sub { my $d = irand() * $_tonv_32;  $d *= $_[0] if $_[0];  $d; };
    } elsif ($nvbits <= 64) {
      *drand = sub { my $d = irand64() * $_tonv_64;  $d *= $_[0] if $_[0];  $d; };
    } else {
      *drand = sub { my $d = irand64() * $_tonv_64 + irand64() * $_tonv_128;  $d *= $_[0] if $_[0];  $d; };
    }
  } else {
    if ($nvbits <= 32) {
      *drand = sub { my $d = irand() * $_tonv_32;  $d *= $_[0] if $_[0];  $d; };
    } elsif ($nvbits <= 64) {
      *drand = sub { my $d = ((irand() >> 5) * 67108864.0 + (irand() >> 6)) / 9007199254740992.0;  $d *= $_[0] if $_[0];  $d; };
    } else {
      *drand = sub { my $d = irand() * $_tonv_32 + irand() * $_tonv_64 + irand() * $_tonv_96 + irand() * $_tonv_128;  $d *= $_[0] if $_[0];  $d; };
    }
  }
  *rand = \&drand;
}


# These functions all do input validation within the PP code.
# Therefore we can send user input straight to them.

# The advantage is simplicity and speed for a single user call.
#
# The disadvantage is that we're doing very expensive PP validation
# for each function call within the PP code itself.

# Rules of thumb:
#   if a function is expensive, no harm in validation
#   if a function is cheap and often called, consider validation here.

# TODO: revisit decision for all of these

*urandomb = \&Math::Prime::Util::PP::urandomb;
*urandomm = \&Math::Prime::Util::PP::urandomm;

*sumdigits = \&Math::Prime::Util::PP::sumdigits;
*todigits = \&Math::Prime::Util::PP::todigits;
*todigitstring = \&Math::Prime::Util::PP::todigitstring;
*fromdigits = \&Math::Prime::Util::PP::fromdigits;
*inverse_li = \&Math::Prime::Util::PP::inverse_li;
*sieve_prime_cluster = \&Math::Prime::Util::PP::sieve_prime_cluster;
*sieve_range = \&Math::Prime::Util::PP::sieve_range;
*lucky_numbers = \&Math::Prime::Util::PP::lucky_numbers;
*powerful_numbers = \&Math::Prime::Util::PP::powerful_numbers;

*prime_power_count = \&Math::Prime::Util::PP::prime_power_count;
*twin_prime_count = \&Math::Prime::Util::PP::twin_prime_count;
*semiprime_count = \&Math::Prime::Util::PP::semiprime_count;
*almost_prime_count = \&Math::Prime::Util::PP::almost_prime_count;
*omega_prime_count = \&Math::Prime::Util::PP::omega_prime_count;
*ramanujan_prime_count = \&Math::Prime::Util::PP::ramanujan_prime_count;
*lucky_count = \&Math::Prime::Util::PP::lucky_count;

*sum_primes = \&Math::Prime::Util::PP::sum_primes;
*print_primes = \&Math::Prime::Util::PP::print_primes;

*is_prime          = \&Math::Prime::Util::PP::is_prime;
*is_prob_prime     = \&Math::Prime::Util::PP::is_prob_prime;
*is_provable_prime = \&Math::Prime::Util::PP::is_provable_prime;
*is_bpsw_prime     = \&Math::Prime::Util::PP::is_bpsw_prime;
*is_pseudoprime    = \&Math::Prime::Util::PP::is_pseudoprime;
*is_euler_pseudoprime = \&Math::Prime::Util::PP::is_euler_pseudoprime;
*is_strong_pseudoprime = \&Math::Prime::Util::PP::is_strong_pseudoprime;
*is_euler_plumb_pseudoprime = \&Math::Prime::Util::PP::is_euler_plumb_pseudoprime;

*is_carmichael = \&Math::Prime::Util::PP::is_carmichael;
*is_quasi_carmichael = \&Math::Prime::Util::PP::is_quasi_carmichael;
*is_pillai = \&Math::Prime::Util::PP::is_pillai;
*is_fundamental = \&Math::Prime::Util::PP::is_fundamental;
*is_semiprime = \&Math::Prime::Util::PP::is_semiprime;
*is_almost_prime = \&Math::Prime::Util::PP::is_almost_prime;
*is_chen_prime = \&Math::Prime::Util::PP::is_chen_prime;
*is_omega_prime = \&Math::Prime::Util::PP::is_omega_prime;
*is_totient = \&Math::Prime::Util::PP::is_totient;
*is_square = \&Math::Prime::Util::PP::is_square;
*is_lucky = \&Math::Prime::Util::PP::is_lucky;
*is_gaussian_prime = \&Math::Prime::Util::PP::is_gaussian_prime;
*is_smooth = \&Math::Prime::Util::PP::is_smooth;
*is_rough = \&Math::Prime::Util::PP::is_rough;
*is_perfect_power = \&Math::Prime::Util::PP::is_perfect_power;
*is_powerful = \&Math::Prime::Util::PP::is_powerful;
*is_odd = \&Math::Prime::Util::PP::is_odd;
*is_even = \&Math::Prime::Util::PP::is_even;
*is_divisible = \&Math::Prime::Util::PP::is_divisible;
*is_congruent = \&Math::Prime::Util::PP::is_congruent;
*is_congruent_number = \&Math::Prime::Util::PP::is_congruent_number;
*is_perfect_number = \&Math::Prime::Util::PP::is_perfect_number;
*is_delicate_prime = \&Math::Prime::Util::PP::is_delicate_prime;
*powerful_count = \&Math::Prime::Util::PP::powerful_count;
*nth_powerful = \&Math::Prime::Util::PP::nth_powerful;
*sumpowerful = \&Math::Prime::Util::PP::sumpowerful;
*perfect_power_count = \&Math::Prime::Util::PP::perfect_power_count;
*is_square_free = \&Math::Prime::Util::PP::is_square_free;
*is_powerfree = \&Math::Prime::Util::PP::is_powerfree;
*powerfree_count = \&Math::Prime::Util::PP::powerfree_count;
*nth_powerfree = \&Math::Prime::Util::PP::nth_powerfree;
*powerfree_sum = \&Math::Prime::Util::PP::powerfree_sum;
*powerfree_part = \&Math::Prime::Util::PP::powerfree_part;
*powerfree_part_sum = \&Math::Prime::Util::PP::powerfree_part_sum;
*squarefree_kernel = \&Math::Prime::Util::PP::squarefree_kernel;
# TODO: Should this do validation here?
*powersum = \&Math::Prime::Util::PP::powersum;
*next_chen_prime = \&Math::Prime::Util::PP::next_chen_prime;

*random_prime = \&Math::Prime::Util::PP::random_prime;
*random_ndigit_prime = \&Math::Prime::Util::PP::random_ndigit_prime;
*random_nbit_prime = \&Math::Prime::Util::PP::random_nbit_prime;
*random_proven_prime = \&Math::Prime::Util::PP::random_maurer_prime; # redir
*random_safe_prime = \&Math::Prime::Util::PP::random_safe_prime;
*random_strong_prime = \&Math::Prime::Util::PP::random_strong_prime;
*random_maurer_prime = \&Math::Prime::Util::PP::random_maurer_prime;
*random_shawe_taylor_prime =\&Math::Prime::Util::PP::random_shawe_taylor_prime;
*miller_rabin_random = \&Math::Prime::Util::PP::miller_rabin_random;
*random_semiprime = \&Math::Prime::Util::PP::random_semiprime;
*random_unrestricted_semiprime = \&Math::Prime::Util::PP::random_unrestricted_semiprime;
*random_factored_integer = \&Math::Prime::Util::PP::random_factored_integer;

*next_prime = \&Math::Prime::Util::PP::next_prime;
*prev_prime = \&Math::Prime::Util::PP::prev_prime;
*next_prime_power = \&Math::Prime::Util::PP::next_prime_power;
*prev_prime_power = \&Math::Prime::Util::PP::prev_prime_power;
*next_perfect_power = \&Math::Prime::Util::PP::next_perfect_power;
*prev_perfect_power = \&Math::Prime::Util::PP::prev_perfect_power;

*numtoperm = \&Math::Prime::Util::PP::numtoperm;
*permtonum = \&Math::Prime::Util::PP::permtonum;
*randperm = \&Math::Prime::Util::PP::randperm;
*shuffle = \&Math::Prime::Util::PP::shuffle;

*prime_bigomega = \&Math::Prime::Util::PP::prime_bigomega;
*prime_omega = \&Math::Prime::Util::PP::prime_omega;
*moebius = \&Math::Prime::Util::PP::moebius;
*euler_phi = \&Math::Prime::Util::PP::euler_phi;
*inverse_totient = \&Math::Prime::Util::PP::inverse_totient;
*sumtotient = \&Math::Prime::Util::PP::sumtotient;
*valuation = \&Math::Prime::Util::PP::valuation;
*chinese = \&Math::Prime::Util::PP::chinese;
*chinese2 = \&Math::Prime::Util::PP::chinese2;
*cornacchia = \&Math::Prime::Util::PP::cornacchia;
*primorial = \&Math::Prime::Util::PP::primorial;
*pn_primorial = \&Math::Prime::Util::PP::pn_primorial;
*divisors = \&Math::Prime::Util::PP::divisors;
*partitions = \&Math::Prime::Util::PP::partitions;
*consecutive_integer_lcm = \&Math::Prime::Util::PP::consecutive_integer_lcm;
*frobenius_number = \&Math::Prime::Util::PP::frobenius_number;
*subfactorial = \&Math::Prime::Util::PP::subfactorial;
*fubini = \&Math::Prime::Util::PP::fubini;
*falling_factorial = \&Math::Prime::Util::PP::falling_factorial;
*rising_factorial = \&Math::Prime::Util::PP::rising_factorial;

*divint = \&Math::Prime::Util::PP::divint;
*modint = \&Math::Prime::Util::PP::modint;
*cdivint = \&Math::Prime::Util::PP::cdivint;
*divrem = \&Math::Prime::Util::PP::divrem;
*tdivrem = \&Math::Prime::Util::PP::tdivrem;
*fdivrem = \&Math::Prime::Util::PP::fdivrem;
*cdivrem = \&Math::Prime::Util::PP::cdivrem;
*absint = \&Math::Prime::Util::PP::absint;
*negint = \&Math::Prime::Util::PP::negint;
*signint = \&Math::Prime::Util::PP::signint;
*cmpint = \&Math::Prime::Util::PP::cmpint;
*add1int = \&Math::Prime::Util::PP::add1int;
*sub1int = \&Math::Prime::Util::PP::sub1int;

*negmod = \&Math::Prime::Util::PP::negmod;
*sqrtmod = \&Math::Prime::Util::PP::sqrtmod;
*allsqrtmod = \&Math::Prime::Util::PP::allsqrtmod;
*rootmod = \&Math::Prime::Util::PP::rootmod;
*allrootmod = \&Math::Prime::Util::PP::allrootmod;
*factorialmod = \&Math::Prime::Util::PP::factorialmod;
*binomialmod = \&Math::Prime::Util::PP::binomialmod;
*lucasumod = \&Math::Prime::Util::PP::lucasumod;
*lucasvmod = \&Math::Prime::Util::PP::lucasvmod;
*lucasuvmod = \&Math::Prime::Util::PP::lucasuvmod;
*pisano_period = \&Math::Prime::Util::PP::pisano_period;
*znlog = \&Math::Prime::Util::PP::znlog;
*znorder = \&Math::Prime::Util::PP::znorder;
*znprimroot = \&Math::Prime::Util::PP::znprimroot;
*is_primitive_root = \&Math::Prime::Util::PP::is_primitive_root;
*qnr = \&Math::Prime::Util::PP::qnr;
*is_qr = \&Math::Prime::Util::PP::is_qr;


*vecequal = \&Math::Prime::Util::PP::vecequal;
*tozeckendorf = \&Math::Prime::Util::PP::tozeckendorf;
*fromzeckendorf = \&Math::Prime::Util::PP::fromzeckendorf;

*prime_count_approx = \&Math::Prime::Util::PP::prime_count_approx;
*prime_count_lower = \&Math::Prime::Util::PP::prime_count_lower;
*prime_count_upper = \&Math::Prime::Util::PP::prime_count_upper;
*prime_power_count_approx = \&Math::Prime::Util::PP::prime_power_count_approx;
*prime_power_count_lower = \&Math::Prime::Util::PP::prime_power_count_lower;
*prime_power_count_upper = \&Math::Prime::Util::PP::prime_power_count_upper;
*perfect_power_count_approx = \&Math::Prime::Util::PP::perfect_power_count_approx;
*perfect_power_count_lower = \&Math::Prime::Util::PP::perfect_power_count_lower;
*perfect_power_count_upper = \&Math::Prime::Util::PP::perfect_power_count_upper;
*lucky_count_approx = \&Math::Prime::Util::PP::lucky_count_approx;
*lucky_count_lower = \&Math::Prime::Util::PP::lucky_count_lower;
*lucky_count_upper = \&Math::Prime::Util::PP::lucky_count_upper;
*nth_prime_power = \&Math::Prime::Util::PP::nth_prime_power;
*nth_prime_power_approx = \&Math::Prime::Util::PP::nth_prime_power_approx;
*nth_prime_power_lower = \&Math::Prime::Util::PP::nth_prime_power_lower;
*nth_prime_power_upper = \&Math::Prime::Util::PP::nth_prime_power_upper;
*nth_perfect_power = \&Math::Prime::Util::PP::nth_perfect_power;
*nth_perfect_power_approx = \&Math::Prime::Util::PP::nth_perfect_power_approx;
*nth_perfect_power_lower = \&Math::Prime::Util::PP::nth_perfect_power_lower;
*nth_perfect_power_upper = \&Math::Prime::Util::PP::nth_perfect_power_upper;
#*nth_lucky = \&Math::Prime::Util::PP::nth_lucky;
*nth_lucky_approx = \&Math::Prime::Util::PP::nth_lucky_approx;
*nth_lucky_lower = \&Math::Prime::Util::PP::nth_lucky_lower;
*nth_lucky_upper = \&Math::Prime::Util::PP::nth_lucky_upper;

*semiprime_count_approx = \&Math::Prime::Util::PP::semiprime_count_approx;
*nth_semiprime_approx = \&Math::Prime::Util::PP::nth_semiprime_approx;
*twin_prime_count_approx = \&Math::Prime::Util::PP::twin_prime_count_approx;
*nth_twin_prime_approx = \&Math::Prime::Util::PP::nth_twin_prime_approx;
*nth_semiprime = \&Math::Prime::Util::PP::nth_semiprime;

*almost_prime_count_approx = \&Math::Prime::Util::PP::almost_prime_count_approx;

*prho_factor = \&Math::Prime::Util::PP::prho_factor;
*pbrent_factor = \&Math::Prime::Util::PP::pbrent_factor;
*ecm_factor = \&Math::Prime::Util::PP::ecm_factor;
*trial_factor = \&Math::Prime::Util::PP::trial_factor;
*fermat_factor = \&Math::Prime::Util::PP::fermat_factor;
*holf_factor = \&Math::Prime::Util::PP::holf_factor;
*squfof_factor = \&Math::Prime::Util::PP::squfof_factor;
*lehman_factor = \&Math::Prime::Util::PP::lehman_factor;
*pminus1_factor = \&Math::Prime::Util::PP::pminus1_factor;
*pplus1_factor = \&Math::Prime::Util::PP::pplus1_factor;
*cheb_factor = \&Math::Prime::Util::PP::cheb_factor;

*primes = \&Math::Prime::Util::PP::primes;
*prime_powers = \&Math::Prime::Util::PP::prime_powers;
*twin_primes = \&Math::Prime::Util::PP::twin_primes;
*semi_primes = \&Math::Prime::Util::PP::semi_primes;
*ramanujan_primes = \&Math::Prime::Util::PP::ramanujan_primes;
*almost_primes = \&Math::Prime::Util::PP::almost_primes;
*omega_primes = \&Math::Prime::Util::PP::omega_primes;


# We are doing the validation here so the PP code doesn't have to do it.

sub jordan_totient {
  my($k, $n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::jordan_totient($k, $n);
}
sub ramanujan_sum {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::ramanujan_sum($k, $n);
}
sub carmichael_lambda {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::carmichael_lambda($n);
}
sub mertens {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::mertens($n);
}
sub liouville {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::liouville($n);
}
sub sumliouville {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::sumliouville($n);
}
sub exp_mangoldt {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return 1 if $n <= 1;
  return Math::Prime::Util::PP::exp_mangoldt($n);
}
sub hclassno {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::hclassno($n);
}


sub nth_prime {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_prime($n);
}
sub nth_prime_lower {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_prime_lower($n);
}
sub nth_prime_upper {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_prime_upper($n);
}
sub nth_prime_approx {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_prime_approx($n);
}
sub almost_prime_count_lower {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::almost_prime_count_lower($k,$n);
}
sub almost_prime_count_upper {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::almost_prime_count_upper($k,$n);
}
sub ramanujan_prime_count_lower {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::ramanujan_prime_count_lower($n);
}
sub ramanujan_prime_count_upper {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::ramanujan_prime_count_upper($n);
}
sub ramanujan_prime_count_approx {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::ramanujan_prime_count_approx($n);
}
sub nth_twin_prime {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_twin_prime($n);
}
sub nth_almost_prime {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_almost_prime($k,$n);
}
sub nth_almost_prime_approx {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_almost_prime_approx($k,$n);
}
sub nth_almost_prime_lower {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_almost_prime_lower($k,$n);
}
sub nth_almost_prime_upper {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_almost_prime_upper($k,$n);
}
sub nth_omega_prime {
  my($k,$n) = @_;
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_omega_prime($k,$n);
}
sub nth_ramanujan_prime {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_ramanujan_prime($n);
}
sub nth_ramanujan_prime_lower {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_ramanujan_prime_lower($n);
}
sub nth_ramanujan_prime_upper {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_ramanujan_prime_upper($n);
}
sub nth_ramanujan_prime_approx {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_ramanujan_prime_approx($n);
}
sub nth_lucky {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::nth_lucky($n);
}
sub smooth_count {
  my($n, $k) = @_;
  _validate_integer_nonneg($n);
  _validate_integer_nonneg($k);
  return Math::Prime::Util::PP::smooth_count($n, $k);
}
sub rough_count {
  my($n, $k) = @_;
  _validate_integer_nonneg($n);
  _validate_integer_nonneg($k);
  return Math::Prime::Util::PP::rough_count($n, $k);
}


sub is_lucas_pseudoprime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_lucas_pseudoprime($n);
}
sub is_strong_lucas_pseudoprime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_strong_lucas_pseudoprime($n);
}
sub is_extra_strong_lucas_pseudoprime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_extra_strong_lucas_pseudoprime($n);
}
sub is_almost_extra_strong_lucas_pseudoprime {
  my($n, $increment) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  if (defined $increment) { _validate_integer_nonneg($increment); }
  else                    { $increment = 1; }
  croak "aes lucas pseudoprime parameter must be 1-256"
    if $increment < 1 || $increment > 256;
  return Math::Prime::Util::PP::is_almost_extra_strong_lucas_pseudoprime($n, $increment);
}
sub is_perrin_pseudoprime {
  my($n,$restrict) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  if (defined $restrict) { _validate_integer_nonneg($restrict); }
  else                   { $restrict = 0; }
  return Math::Prime::Util::PP::is_perrin_pseudoprime($n, $restrict);
}
sub is_catalan_pseudoprime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_catalan_pseudoprime($n);
}
sub is_frobenius_pseudoprime {
  my($n, $P, $Q) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  # TODO: validate P & Q
  return Math::Prime::Util::PP::is_frobenius_pseudoprime($n, $P, $Q);
}
sub is_frobenius_underwood_pseudoprime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_frobenius_underwood_pseudoprime($n);
}
sub is_frobenius_khashin_pseudoprime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_frobenius_khashin_pseudoprime($n);
}
sub is_aks_prime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_aks_prime($n);
}
sub is_ramanujan_prime {
  my($n) = @_;
  _validate_integer($n);
  return 0 if $n < 0;
  return Math::Prime::Util::PP::is_ramanujan_prime($n);
}
sub is_mersenne_prime {
  my($p) = @_;
  _validate_integer_nonneg($p);
  return Math::Prime::Util::PP::is_mersenne_prime($p);
}
sub lucas_sequence {
  my($n, $P, $Q, $k) = @_;
  my ($vp, $vq) = ($P, $Q);
  _validate_integer_positive($n);
  _validate_integer($vp);
  _validate_integer($vq);
  _validate_integer_nonneg($k);
  $vp = -$vp if $vp < 0;
  $vq = -$vq if $vq < 0;
  return Math::Prime::Util::PP::lucas_sequence(@_);
}
sub lucasu {
  my($P, $Q, $k) = @_;
  my ($vp, $vq) = ($P, $Q);
  _validate_integer($vp);
  _validate_integer($vq);
  _validate_integer_nonneg($k);
  $vp = -$vp if $vp < 0;
  $vq = -$vq if $vq < 0;
  return Math::Prime::Util::PP::lucasu($P,$Q,$k);
}
sub lucasv {
  my($P, $Q, $k) = @_;
  my ($vp, $vq) = ($P, $Q);
  _validate_integer($vp);
  _validate_integer($vq);
  _validate_integer_nonneg($k);
  $vp = -$vp if $vp < 0;
  $vq = -$vq if $vq < 0;
  return Math::Prime::Util::PP::lucasv($P,$Q,$k);
}
sub lucasuv {
  my($P, $Q, $k) = @_;
  _validate_integer($P);
  _validate_integer($Q);
  _validate_integer_nonneg($k);
  return Math::Prime::Util::PP::lucasuv($P,$Q,$k);
}

sub kronecker {
  my($a, $b) = @_;
  my ($va, $vb) = ($a, $b);
  _validate_integer($va);
  _validate_integer($vb);
  $vb = -$va if $va < 0;
  $vb = -$vb if $vb < 0;
  return Math::Prime::Util::PP::kronecker(@_);
}

sub factorial {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::factorial($n);
}

sub binomial {
  my($n, $k) = @_;
  _validate_integer($n);
  _validate_integer($k);
  return Math::Prime::Util::PP::binomial($n, $k);
}

sub stirling {
  my($n, $k, $type) = @_;
  _validate_integer_nonneg($n);
  _validate_integer_nonneg($k);
  _validate_integer_nonneg($type) if defined $type;
  return Math::Prime::Util::PP::stirling($n, $k, $type);
}

sub divisor_sum {
  my($n, $k) = @_;
  _validate_integer_nonneg($n);
  _validate_integer_nonneg($k) if defined $k && ref($k) ne 'CODE';
  return Math::Prime::Util::PP::divisor_sum($n, $k);
}

sub gcd {
  my(@v) = @_;
  _validate_integer($_) for @v;
  return Math::Prime::Util::PP::gcd(@v);
}
sub lcm {
  my(@v) = @_;
  _validate_integer($_) for @v;
  return Math::Prime::Util::PP::lcm(@v);
}
sub gcdext {
  my($a,$b) = @_;
  _validate_integer($a);
  _validate_integer($b);
  return Math::Prime::Util::PP::gcdext($a,$b);
}
sub vecsum {
  my(@v) = @_;
  _validate_integer($_) for @v;
  return Math::Prime::Util::PP::vecsum(@v);
}
sub vecprod {
  my(@v) = @_;
  _validate_integer($_) for @v;
  return Math::Prime::Util::PP::vecprod(@v);
}
sub vecmin {
  my(@v) = @_;
  _validate_integer($_) for @v;
  return Math::Prime::Util::PP::vecmin(@v);
}
sub vecmax {
  my(@v) = @_;
  _validate_integer($_) for @v;
  return Math::Prime::Util::PP::vecmax(@v);
}
sub vecmex {
  my(@v) = @_;
  _validate_integer_nonneg($_) for @v;
  return Math::Prime::Util::PP::vecmex(@v);
}
sub vecpmex {
  my(@v) = @_;
  for (@v) {
    _validate_integer_nonneg($_);
    croak "parameter must be a positive integer (x > 0)" if $_ <= 0;
  }
  return Math::Prime::Util::PP::vecpmex(@v);
}
sub invmod {
  my ($a, $n) = @_;
  _validate_integer($a);
  _validate_integer($n);
  return Math::Prime::Util::PP::invmod($a,$n);
}
sub addmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b); _validate_integer($n);
  return Math::Prime::Util::PP::addmod($a,$b, $n);
}
sub submod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b); _validate_integer($n);
  return Math::Prime::Util::PP::submod($a,$b, $n);
}
sub mulmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b); _validate_integer($n);
  return Math::Prime::Util::PP::mulmod($a,$b, $n);
}
sub divmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b); _validate_integer($n);
  return Math::Prime::Util::PP::divmod($a,$b, $n);
}
sub powmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b); _validate_integer($n);
  return Math::Prime::Util::PP::powmod($a,$b, $n);
}
sub muladdmod {
  my ($a, $b, $c, $n) = @_;
  _validate_integer($a); _validate_integer($b); _validate_integer($c); _validate_integer($n);
  return Math::Prime::Util::PP::muladdmod($a,$b,$c, $n);
}
sub mulsubmod {
  my ($a, $b, $c, $n) = @_;
  _validate_integer($a); _validate_integer($b); _validate_integer($c); _validate_integer($n);
  return Math::Prime::Util::PP::mulsubmod($a,$b,$c, $n);
}
sub sqrtint {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::sqrtint($n);
}
sub rootint {
  my($n, $k, $refp) = @_;
  _validate_integer_nonneg($n);
  _validate_integer_nonneg($k);
  return Math::Prime::Util::PP::rootint($n, $k, $refp);
}
sub logint {
  my($n, $b, $refp) = @_;
  _validate_integer_positive($n);
  _validate_integer_nonneg($b);
  return Math::Prime::Util::PP::logint($n, $b, $refp);
}
sub powint {
  my($a, $b) = @_;
  _validate_integer($a);
  _validate_integer($b);
  return Math::Prime::Util::PP::powint($a, $b);
}
sub mulint {
  my($a, $b) = @_;
  _validate_integer($a);
  _validate_integer($b);
  return Math::Prime::Util::PP::mulint($a, $b);
}
sub addint {
  my($a, $b) = @_;
  _validate_integer($a);
  _validate_integer($b);
  return Math::Prime::Util::PP::addint($a, $b);
}
sub subint {
  my($a, $b) = @_;
  _validate_integer($a);
  _validate_integer($b);
  return Math::Prime::Util::PP::subint($a, $b);
}
sub lshiftint {
  my($n, $k) = @_;
  _validate_integer($n);
  if (!defined $k) { $k = 1; } else { _validate_integer_nonneg($k); }
  return Math::Prime::Util::PP::lshiftint($n, $k);
}
sub rshiftint {
  my($n, $k) = @_;
  _validate_integer($n);
  if (!defined $k) { $k = 1; } else { _validate_integer_nonneg($k); }
  return Math::Prime::Util::PP::rshiftint($n, $k);
}
sub rashiftint {
  my($n, $k) = @_;
  _validate_integer($n);
  if (!defined $k) { $k = 1; } else { _validate_integer_nonneg($k); }
  return Math::Prime::Util::PP::rashiftint($n, $k);
}
sub legendre_phi {
  my($x, $a) = @_;
  _validate_integer_nonneg($x);
  _validate_integer_nonneg($a);
  return Math::Prime::Util::PP::legendre_phi($x, $a);
}

sub chebyshev_theta {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::chebyshev_theta($n);
}
sub chebyshev_psi {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::chebyshev_psi($n);
}
sub ramanujan_tau {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::ramanujan_tau($n);
}

sub is_power {
  my($n, $a, $refp) = @_;
  _validate_integer($n);
  _validate_integer_nonneg($a) if defined $a;
  return Math::Prime::Util::PP::is_power($n, $a, $refp);
}
sub is_prime_power {
  my($n, $refp) = @_;
  _validate_integer($n);
  return Math::Prime::Util::PP::is_prime_power($n, $refp);
}
sub is_polygonal {
  my($n, $s, $refp) = @_;
  _validate_integer($n);
  _validate_integer_nonneg($s);
  return Math::Prime::Util::PP::is_polygonal($n, $s, $refp);
}
sub is_sum_of_squares {
  my($n, $k) = @_;
  _validate_integer($n);
  $n = -$n if $n < 0;
  if (!defined $k) { $k = 2; } else { _validate_integer_nonneg($k); }
  return Math::Prime::Util::PP::is_sum_of_squares($n, $k);
}
sub is_practical {
  my($n) = @_;
  _validate_integer_nonneg($n);
  return Math::Prime::Util::PP::is_practical($n);
}
sub hammingweight {
  my($n) = @_;
  _validate_integer($n);
  return Math::Prime::Util::PP::hammingweight($n);
}

sub Pi {
  my($digits) = @_;
  _validate_integer_nonneg($digits) if defined $digits;
  return Math::Prime::Util::PP::Pi($digits);
}

#############################################################################

my $_exitloop = 0;
sub lastfor { $_exitloop = 1; }
sub _get_forexit { $_exitloop; }
sub _start_for_loop { my $old = $_exitloop; $_exitloop = 0; $old; }
sub _end_for_loop { $_exitloop = shift; }

sub forprimes (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, $beg, $end) = @_;
  if (scalar @_ > 2) {
    _validate_integer_nonneg($beg);
  } else {
    ($beg,$end) = (2,$beg);
  }
  _validate_integer_nonneg($end);
  $beg = 2 if $beg < 2;
  my $oldforexit = _start_for_loop();
  {
    my $pp;
    local *_ = \$pp;
    for (my $p = next_prime($beg-1);  $p <= $end;  $p = next_prime($p)) {
      $pp = $p;
      $sub->();
      last if $_exitloop;
    }
  }
  _end_for_loop($oldforexit);
}

sub forcomposites(&$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::_generic_forcomp_sub('composites', @_);
}
sub foroddcomposites(&$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::_generic_forcomp_sub('oddcomposites', @_);
}
sub forsemiprimes(&$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::_generic_forcomp_sub('semiprimes', @_);
}
sub foralmostprimes(&$$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::PP::foralmostprimes(@_);
}

sub forfactored(&$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::_generic_forfac(0, @_);
}
sub forsquarefree(&$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::_generic_forfac(1, @_);
}

sub fordivisors (&$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, $n) = @_;
  _validate_integer_nonneg($n);
  my @divisors = divisors($n);
  my $oldforexit = _start_for_loop();
  {
    my $pp;
    local *_ = \$pp;
    foreach my $d (@divisors) {
      $pp = $d;
      $sub->();
      last if $_exitloop;
    }
  }
  _end_for_loop($oldforexit);
}

sub forpart (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::PP::forpart(@_);
}
sub forcomp (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::PP::forcomp(@_);
}
sub forcomb (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::PP::forcomb(@_);
}
sub forperm (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::PP::forperm(@_);
}
sub forderange (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  Math::Prime::Util::PP::forderange(@_);
}

sub forsetproduct (&@) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, @v) = @_;
  croak 'Not a subroutine reference' unless (ref($sub) || '') eq 'CODE';
  croak 'Not an array reference' if grep {(ref($_) || '') ne 'ARRAY'} @v;
  # Exit if no arrays or any are empty.
  return if scalar(@v) == 0 || grep { !@$_ } @v;

  my @outv = map { $v[$_]->[0] } 0 .. $#v;
  my @cnt = (0) x @v;

  my $oldforexit = _start_for_loop();
  my $i = 0;
  while ($i >= 0) {
    $sub->(@outv);
    last if $_exitloop;
    for ($i = $#v; $i >= 0; $i--) {
      if ($cnt[$i] >= $#{$v[$i]}) { $cnt[$i] = 0; $outv[$i] = $v[$i]->[0]; }
      else { $outv[$i] = $v[$i]->[++$cnt[$i]]; last; }
    }
  }
  _end_for_loop($oldforexit);
}

sub vecreduce (&@) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, @v) = @_;

  # Mastering Perl page 162, works with old Perl
  my $caller = caller();
  no strict 'refs'; ## no critic(strict)
  local(*{$caller.'::a'}) = \my $a;
  local(*{$caller.'::b'}) = \my $b;
  $a = shift @v;
  for my $v (@v) {
    $b = $v;
    $a = $sub->();
  }
  $a;
}

sub vecany (&@) {       ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  $sub->() and return 1 foreach @_;
  0;
}
sub vecall (&@) {       ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  $sub->() or return 0 foreach @_;
  1;
}
sub vecnone (&@) {      ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  $sub->() and return 0 foreach @_;
  1;
}
sub vecnotall (&@) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  $sub->() or return 1 foreach @_;
  undef;
}

sub vecfirst (&@) {     ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  $sub->() and return $_ foreach @_;
  undef;
}

sub vecfirstidx (&@) {     ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  my $i = 0;
  ++$i and $sub->() and return $i-1 foreach @_;
  -1;
}

sub vecextract {
  my($aref, $mask) = @_;
  croak "vecextract first argument must be an array reference"
    unless ref($aref) eq 'ARRAY';
  return Math::Prime::Util::PP::vecextract(@_);
}

1;

__END__

=pod

=head1 NAME

Math::Prime::Util::PPFE - PP front end for Math::Prime::Util

=head1 SYNOPSIS

This loads the PP code and adds input validation front ends.  It is only
meant to be used when XS is not used.

=head1 DESCRIPTION

Loads PP module and implements PP front-end functions for all XS code.
This is used only if the XS code is not loaded.

=head1 SEE ALSO

L<Math::Prime::Util>

L<Math::Prime::Util::PP>

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2014-2024 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
