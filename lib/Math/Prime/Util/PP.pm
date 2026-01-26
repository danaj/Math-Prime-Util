package Math::Prime::Util::PP;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::PP::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PP::VERSION = '0.73';
}

BEGIN {
  do { require Math::BigInt;  Math::BigInt->import(try=>"GMPz,GMP,LTM,Pari"); }
    unless defined $Math::BigInt::VERSION;
}

# The Pure Perl versions of all the Math::Prime::Util routines.
#
# Some of these will be relatively similar in performance, some will be
# very slow in comparison.
#
# Most of these are pretty simple.  Also, you really should look at the C
# code for more detailed comments, including references to papers.

BEGIN {
  use constant OLD_PERL_VERSION=> $] < 5.008;
  use constant MPU_MAXBITS     => (~0 == 4294967295) ? 32 : 64;
  use constant MPU_64BIT       => MPU_MAXBITS == 64;
  use constant MPU_32BIT       => MPU_MAXBITS == 32;
 #use constant MPU_MAXPARAM    => MPU_32BIT ? 4294967295 : 18446744073709551615;
 #use constant MPU_MAXDIGITS   => MPU_32BIT ? 10 : 20;
  use constant MPU_MAXPRIME    => MPU_32BIT ? 4294967291 : 18446744073709551557;
  use constant MPU_MAXPRIMEIDX => MPU_32BIT ?  203280221 :  425656284035217743;
  use constant MPU_HALFWORD    => MPU_32BIT ? 65536 : OLD_PERL_VERSION ? 33554432 : 4294967296;
  use constant UVPACKLET       => MPU_32BIT ? 'L' : 'Q';
  use constant MPU_INFINITY    => (65535 > 0+'inf') ? 20**20**20 : 0+'inf';
  use constant INTMAX          => (!OLD_PERL_VERSION || MPU_32BIT) ? ~0 : 562949953421312;
  use constant INTMIN          => (MPU_32BIT ? -2147483648 : !OLD_PERL_VERSION ? -9223372036854775808 : -562949953421312);
  use constant SINTMAX         => (INTMAX >> 1);
  use constant B_PRIM235       => Math::BigInt->new("30");
  use constant PI_TIMES_8      => 25.13274122871834590770114707;
}

# TODO: Change this whole file to use this / tobigint
our $_BIGINT;
*_BIGINT = \$Math::Prime::Util::_BIGINT;


# By using these aliases, we call into the main code instead of
# to the PP function.
#
# If we have turned off XS, then this will call the PPFE or direct function.
# This might be the same, but if the PPFE does input validation it will
# be slower (albeit every call will be validated).
#
# Otherwise, we'll go to the XS function, which will either handle it
# directly (e.g. we've broken down the input into smaller values which
# the XS code can handle), or call the GMP backend, otherwise call here.
#
# For the usual case where we have XS, this is significantly faster.  The
# aliases make the code here much easier to read.  An alternate
# implementation would be to make the perl subs here use a pp_{...} prefix.


*validate_integer = \&Math::Prime::Util::_validate_integer;
*validate_integer_nonneg = \&Math::Prime::Util::_validate_integer_nonneg;
*validate_integer_positive = \&Math::Prime::Util::_validate_integer_positive;
*validate_integer_abs = \&Math::Prime::Util::_validate_integer_abs;
*_bigint_to_int = \&Math::Prime::Util::_bigint_to_int;
*reftyped = \&Math::Prime::Util::_reftyped;
#*load_bigint = \&Math::Prime::Util::_load_bigint;
*tobigint = \&Math::Prime::Util::_to_bigint;
*maybetobigint = \&Math::Prime::Util::_to_bigint_if_needed;
*maybetobigintall = \&Math::Prime::Util::_maybe_bigint_allargs;
*getconfig = \&Math::Prime::Util::prime_get_config;

*Maddint = \&Math::Prime::Util::addint;
*Msubint = \&Math::Prime::Util::subint;
*Mmulint = \&Math::Prime::Util::mulint;
*Mdivint = \&Math::Prime::Util::divint;
*Mpowint = \&Math::Prime::Util::powint;
*Mmodint = \&Math::Prime::Util::modint;
*Mcdivint = \&Math::Prime::Util::cdivint;
*Mabsint = \&Math::Prime::Util::absint;
*Msqrtint = \&Math::Prime::Util::sqrtint;
*Mrootint = \&Math::Prime::Util::rootint;
*Mlogint = \&Math::Prime::Util::logint;
*Mnegint = \&Math::Prime::Util::negint;
*Mcmpint = \&Math::Prime::Util::cmpint;
*Mlshiftint = \&Math::Prime::Util::lshiftint;
*Mrshiftint = \&Math::Prime::Util::rshiftint;
*Mdivrem  = \&Math::Prime::Util::divrem;
*Mtdivrem = \&Math::Prime::Util::tdivrem;

*Maddmod = \&Math::Prime::Util::addmod;
*Msubmod = \&Math::Prime::Util::submod;
*Mmulmod = \&Math::Prime::Util::mulmod;
*Mdivmod = \&Math::Prime::Util::divmod;
*Mpowmod = \&Math::Prime::Util::powmod;
*Minvmod = \&Math::Prime::Util::invmod;
*Mrootmod = \&Math::Prime::Util::rootmod;
*Mmuladdmod = \&Math::Prime::Util::muladdmod;
*Mmulsubmod = \&Math::Prime::Util::mulsubmod;

*Mgcd = \&Math::Prime::Util::gcd;
*Mlcm = \&Math::Prime::Util::lcm;
*Mgcdext = \&Math::Prime::Util::gcdext;
*Mfactor = \&Math::Prime::Util::factor;
*Mfactor_exp = \&Math::Prime::Util::factor_exp;
*Mtrial_factor = \&Math::Prime::Util::trial_factor;
*Mdivisors = \&Math::Prime::Util::divisors;
*Mdivisor_sum = \&Math::Prime::Util::divisor_sum;
*Mis_prime = \&Math::Prime::Util::is_prime;
*Mis_semiprime = \&Math::Prime::Util::is_semiprime;
*Mis_prime_power = \&Math::Prime::Util::is_prime_power;
*Mis_power = \&Math::Prime::Util::is_power;
*Mis_square_free = \&Math::Prime::Util::is_square_free;
*Mis_odd = \&Math::Prime::Util::is_odd;
*Mis_even = \&Math::Prime::Util::is_even;
*Mis_congruent = \&Math::Prime::Util::is_congruent;
*Mis_divisible = \&Math::Prime::Util::is_divisible;
*Mchinese = \&Math::Prime::Util::chinese;
*Mvaluation = \&Math::Prime::Util::valuation;
*Mkronecker = \&Math::Prime::Util::kronecker;
*Mmoebius = \&Math::Prime::Util::moebius;
*Mtotient = \&Math::Prime::Util::euler_phi;
*Mfactorial = \&Math::Prime::Util::factorial;
*Mfalling_factorial = \&Math::Prime::Util::falling_factorial;
*Mprimorial = \&Math::Prime::Util::primorial;
*Mpn_primorial = \&Math::Prime::Util::pn_primorial;
*Mbinomial = \&Math::Prime::Util::binomial;
*Mstirling = \&Math::Prime::Util::stirling;
*Mpowersum = \&Math::Prime::Util::powersum;
*Murandomm = \&Math::Prime::Util::urandomm;
*Murandomb = \&Math::Prime::Util::urandomb;
*Mnext_prime = \&Math::Prime::Util::next_prime;
*Mprev_prime = \&Math::Prime::Util::prev_prime;
*Mprime_count = \&Math::Prime::Util::prime_count;
*Mlucasumod = \&Math::Prime::Util::lucasumod;
*Mznorder = \&Math::Prime::Util::znorder;
*Mhclassno = \&Math::Prime::Util::hclassno;

*Mvecall = \&Math::Prime::Util::vecall;
*Mvecany = \&Math::Prime::Util::vecany;
*Mvecnone = \&Math::Prime::Util::vecnone;
*Mvecsum = \&Math::Prime::Util::vecsum;
*Mvecprod = \&Math::Prime::Util::vecprod;
*Mvecmin = \&Math::Prime::Util::vecmin;
*Mvecmax = \&Math::Prime::Util::vecmax;
*Mvecfirst = \&Math::Prime::Util::vecfirst;
*Mvecsort = \&Math::Prime::Util::vecsort;
*Mvecsorti = \&Math::Prime::Util::vecsorti;
*Mtoset = \&Math::Prime::Util::toset;
*Msetinsert = \&Math::Prime::Util::setinsert;
*Msetcontains = \&Math::Prime::Util::setcontains;
*Msetunion = \&Math::Prime::Util::setunion;
*Msetintersect = \&Math::Prime::Util::setintersect;

*Mfromdigits = \&Math::Prime::Util::fromdigits;
*Mtodigits = \&Math::Prime::Util::todigits;
*Mtodigitstring = \&Math::Prime::Util::todigitstring;

*Mprimes = \&Math::Prime::Util::primes;
*Mfordivisors = \&Math::Prime::Util::fordivisors;
*Mforprimes = \&Math::Prime::Util::forprimes;
*MLi = \&Math::Prime::Util::LogarithmicIntegral;
*Mprime_omega = \&Math::Prime::Util::prime_omega;
*Mnth_prime_upper = \&Math::Prime::Util::nth_prime_upper;

if (defined $Math::Prime::Util::GMP::VERSION && $Math::Prime::Util::GMP::VERSION >= 0.53) {
  *Saddint = \&Math::Prime::Util::GMP::addint;
  *Ssubint = \&Math::Prime::Util::GMP::subint;
  *Smulint = \&Math::Prime::Util::GMP::mulint;
  *Sdivint = \&Math::Prime::Util::GMP::divint;
  *Spowint = \&Math::Prime::Util::GMP::powint;
} else {
  *Saddint = \&Math::Prime::Util::addint;
  *Ssubint = \&Math::Prime::Util::subint;
  *Smulint = \&Math::Prime::Util::mulint;
  *Sdivint = \&Math::Prime::Util::divint;
  *Spowint = \&Math::Prime::Util::powint;
}

# We don't have this function yet.  Use a simple version for now.
*Mtoint = \&_toint_simple;

my $_precalc_size = 0;
sub prime_precalc {
  my($n) = @_;
  croak "Parameter '$n' must be a non-negative integer" unless _is_nonneg_int($n);
  $_precalc_size = $n if $n > $_precalc_size;
}
sub prime_memfree {
  $_precalc_size = 0;
  eval { Math::Prime::Util::GMP::_GMP_memfree(); }
    if defined $Math::Prime::Util::GMP::VERSION && $Math::Prime::Util::GMP::VERSION >= 0.49;
}
sub _get_prime_cache_size { $_precalc_size }
sub _prime_memfreeall { prime_memfree; }


sub _is_nonneg_int {
  ((defined $_[0]) && $_[0] ne '' && ($_[0] !~ tr/0123456789//c));
}

sub _upgrade_to_float {
  do { require Math::BigFloat; Math::BigFloat->import(); }
    if !defined $Math::BigFloat::VERSION;
  Math::BigFloat->new(@_);
}

# Get the accuracy of variable x, or the max default from BigInt/BigFloat
# One might think to use ref($x)->accuracy() but numbers get upgraded and
# downgraded willy-nilly, and it will do the wrong thing from the user's
# perspective.
sub _find_big_acc {
  my($x) = @_;
  my $b;

  $b = $x->accuracy() if ref($x) =~ /^Math::Big/;
  return $b if defined $b;

  my ($i,$f) = (Math::BigInt->accuracy(), Math::BigFloat->accuracy());
  return (($i > $f) ? $i : $f)   if defined $i && defined $f;
  return $i if defined $i;
  return $f if defined $f;

  ($i,$f) = (Math::BigInt->div_scale(), Math::BigFloat->div_scale());
  return (($i > $f) ? $i : $f)   if defined $i && defined $f;
  return $i if defined $i;
  return $f if defined $f;
  return 18;
}

sub _bfdigits {
  my($wantbf, $xdigits) = (0, 17);
  if (defined $bignum::VERSION || ref($_[0]) =~ /^Math::Big/) {
    do { require Math::BigFloat; Math::BigFloat->import(); }
      if !defined $Math::BigFloat::VERSION;
    if (ref($_[0]) eq 'Math::BigInt') {
      my $xacc = ($_[0])->accuracy();
      $_[0] = Math::BigFloat->new($_[0]);
      ($_[0])->accuracy($xacc) if $xacc;
    }
    $_[0] = Math::BigFloat->new("$_[0]") if ref($_[0]) ne 'Math::BigFloat';
    $wantbf = _find_big_acc($_[0]);
    $xdigits = $wantbf;
  }
  ($wantbf, $xdigits);
}


sub _validate_integer {
  my($n) = @_;
  croak "Parameter must be defined" if !defined $n;

  my $refn = ref($n);

  if (!$refn) {   # Typical case, an integer or string
    croak "Parameter '$n' must be an integer"
      if $n eq '' || ($n =~ tr/0123456789//c && $n !~ /^([+-]?)\d+\z/);
    substr($_[0],0,1,'') if $1 && (substr($n,0,1) eq '+' || $n eq '-0');
    $_[0] = maybetobigint($n) if $n >= INTMAX || $n <= INTMIN;
  } elsif ($refn eq 'CODE') {
    $_[0] = $_[0]->();
    return _validate_integer($_[0]);
  } elsif ($refn !~ /^Math::/) {
     $_[0] = "$_[0]";
     return _validate_integer($_[0]);
  } else {
    if ($refn =~ /^Math::Big(Int|Float)$/) {
      croak "Parameter '$n' must be an integer" unless $n->is_int();
      my $bits = length($n->as_bin) - 2;
      $_[0] = _bigint_to_int($_[0])
        if $bits <= MPU_MAXBITS || ($bits == MPU_MAXBITS+1 && $n == INTMIN);
    } else {
      $_[0] = _bigint_to_int($_[0]) if $n <= INTMAX && $n >= INTMIN;
    }
  }
  $_[0]->upgrade(undef) if ref($_[0]) eq 'Math::BigInt' && $_[0]->upgrade();
  1;
}
sub _validate_integer_nonneg {
  my($n) = @_;
  croak "Parameter must be defined" if !defined $n;

  my $refn = ref($n);

  if (!$refn) {   # Typical case, an integer or string
    croak "Parameter '$n' must be a non-negative integer"
      if $n eq '' || ($n =~ tr/0123456789//c && $n !~ /^(\+?)\d+\z/) || $n < 0;
    substr($_[0],0,1,'') if $1 && substr($n,0,1) eq '+';
    # If probably a bigint, do the upgrade, then verify for edge cases.
    $_[0] = maybetobigint($n) if $n >= INTMAX;
  } elsif ($refn eq 'CODE') {
    $_[0] = $_[0]->();
    return _validate_integer_nonneg($_[0]);
  } elsif ($refn !~ /^Math::/) {
     $_[0] = "$_[0]";
     return _validate_integer_nonneg($_[0]);
  } else {
    croak "Parameter '$n' must be a non-negative integer"
      if ($refn =~ /^Math::Big(Int|Float)$/ && !$n->is_int()) || $n < 0;
    $_[0] = _bigint_to_int($_[0]) if $n <= INTMAX;
  }
  $_[0]->upgrade(undef) if ref($_[0]) eq 'Math::BigInt' && $_[0]->upgrade();
  1;
}
sub _validate_integer_positive {
  _validate_integer($_[0]);
  croak "Parameter '$_[0]' must be a positive integer"
    if $_[0] < 1;
  1;
}
sub _validate_integer_abs {
  _validate_integer($_[0]);
  $_[0] = -$_[0] if $_[0] < 0;
  1;
}

# If we try to call the function in any normal way, just loading this module
# will auto-vivify an empty sub.  So we do a string eval to keep it hidden.
sub _gmpcall {
  my($fname, @args) = @_;
  my $call = "Math::Prime::Util::GMP::$fname(".join(",",map {"\"$_\""} @args).");";
  return eval $call ## no critic qw(ProhibitStringyEval)
}

sub _binary_search {
  my($n, $lo, $hi, $sub, $exitsub) = @_;
  while ($lo < $hi) {
    my $mid = $lo + int(($hi-$lo) >> 1);
    return $mid if defined $exitsub && $exitsub->($n,$lo,$hi);
    if ($sub->($mid) < $n) { $lo = $mid+1; }
    else                   { $hi = $mid;   }
  }
  return $lo-1;
}

################################################################################

# TODO: this is in progress.
#   It's TBD what should be done on failures (undef? croak?)
#   Handling of trivial floats is terrible.
#   A single native int should be as fast as possible
sub _toint {
  my @v = @_;  # copy them all
  my @out;
  for my $v (@v) {
    if (!defined $v) { push @out, 0; next; }
    if (ref($v)) {
      $v = $v->as_int() if ref($v) eq 'Math::BigFloat';
    } elsif ($v =~ /^[+-]?\d+\z/) {
      # Good as-is
    } elsif ($v =~ /e/i || $v =~ /\./) {
      $v = _upgrade_to_float($v)->as_int();
    } else {
      $v = int($v);
    }
    if ($v =~ /^nan\z/i) { push @out, undef; next; }

    validate_integer($v);
    push @out, $v;
  }
  @out;
}

sub _toint_simple {
  my $n = shift;
  if ($n >= 0) {
    my $max = MPU_32BIT ? 4294967295 : 70368744177664;  # 2^46
    if ($n =~ /^[+]?\d+\z/) {
      return int("$n") if $n < $max;
    } elsif ($n < $max) {
      return int("$n");
    } else {
      $n = _upgrade_to_float("$n")->as_int;
    }
  } else {
    my $min = MPU_32BIT ? -2147483648 : -35184372088832;  # -2^45
    if ($n =~ /^[-]\d+\z/) {
      return int($n) if $n > $min;
    } elsif ($n > $min) {
      return int($n);
    } else {
      $n = _upgrade_to_float($n)->as_int;
    }
  }
  validate_integer($n);
  $n = tobigint($n) if ref($n) && defined $_BIGINT && ref($n) ne $_BIGINT;
  $n;
}

sub _frombinary {
  my($bstr) = @_;
  $bstr =~ s/^0//;
  return oct('0b' . $bstr) if length($bstr) <= 32;
  # Avoid the useless portable warning that can't be silenced.
  if (MPU_MAXBITS >= 64 && length($bstr) <= 64) {  # 64-bit Perl, 33-64 bit str
    my $low = substr($bstr,-32,32,'');
    return oct('0b'.$bstr) << 32 + oct('0b'.$low);
  }
  # Length is bigger than word size, so must be a bigint
  if (!defined $_BIGINT) {
    return Math::BigInt->new("0b$bstr");
  } elsif ($_BIGINT =~ /^Math::(BigInt|GMPz|GMP)$/) {
    return $_BIGINT->new("0b$bstr");
  } else {
    return tobigint( Math::BigInt->new("0b$bstr") );
  }
}

################################################################################
################################################################################

my @_primes_small = (0,2);
{
  my($n, $s, $sieveref) = (7-2, 3, _sieve_erat_string(5003));
  push @_primes_small, 2*pos($$sieveref)-1 while $$sieveref =~ m/0/g;
}
my @_prime_next_small = (
   2,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,
   29,29,29,29,29,29,31,31,37,37,37,37,37,37,41,41,41,41,43,43,47,
   47,47,47,53,53,53,53,53,53,59,59,59,59,59,59,61,61,67,67,67,67,67,67,71);

# For wheel-30
my @_prime_indices = (1, 7, 11, 13, 17, 19, 23, 29);
my @_nextwheel30 = (1,7,7,7,7,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,29,29,29,29,29,29,1);
my @_prevwheel30 = (29,29,1,1,1,1,1,1,7,7,7,7,11,11,13,13,13,13,17,17,19,19,19,19,23,23,23,23,23,23);
my @_wheeladvance30 = (1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2);
my @_wheelretreat30 = (1,2,1,2,3,4,5,6,1,2,3,4,1,2,1,2,3,4,1,2,1,2,3,4,1,2,3,4,5,6);

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

sub _is_prime7 {  # n must not be divisible by 2, 3, or 5
  my($n) = @_;

  $n = _bigint_to_int($n) if ref($n) && $n <= INTMAX;

  if (ref($n)) {
    # Check div by 7,11,13,17,19,23,29;  then by 31,37,...,109,113
    return 0 unless Mgcd($n,215656441) == 1;
    return 0 unless Mgcd($n,'4885866070719029716366506343847722513') == 1;
    return 0 unless _miller_rabin_2($n);
    if (Mcmpint($n,"18446744073709551615") <= 0) {
      return is_almost_extra_strong_lucas_pseudoprime($n) ? 2 : 0;
    }
    return is_extra_strong_lucas_pseudoprime($n) ? 1 : 0;
  }

  if ($n < 61*61) {
    foreach my $i (qw/7 11 13 17 19 23 29 31 37 41 43 47 53 59/) {
      return 2 if $i*$i > $n;
      return 0 if !($n % $i);
    }
    return 2;
  }

  return 0 if !($n %  7) || !($n % 11) || !($n % 13) || !($n % 17) ||
              !($n % 19) || !($n % 23) || !($n % 29) || !($n % 31) ||
              !($n % 37) || !($n % 41) || !($n % 43) || !($n % 47) ||
              !($n % 53) || !($n % 59);

  # We could do:
  #   return is_strong_pseudoprime($n, (2,299417)) if $n < 19471033;
  # or:
  #   foreach my $p (@_primes_small[18..168]) {
  #     last if $p > $limit;
  #     return 0 unless $n % $p;
  #   }
  #   return 2;

  if ($n <= 1_500_000) {
    my $limit = int(sqrt($n));
    my $i = 61;
    while (($i+30) <= $limit) {
      return 0 unless ($n% $i    ) && ($n%($i+ 6)) &&
                      ($n%($i+10)) && ($n%($i+12)) &&
                      ($n%($i+16)) && ($n%($i+18)) &&
                      ($n%($i+22)) && ($n%($i+28));
      $i += 30;
    }
    for my $inc (6,4,2,4,2,4,6,2) {
      last if $i > $limit;
      return 0 if !($n % $i);
      $i += $inc;
    }
    return 2;
  }

  if ($n < 47636622961201) {   # BPSW seems to be faster after this
    # Deterministic set of Miller-Rabin tests.  If the MR routines can handle
    # bases greater than n, then this can be simplified.
    my @bases;
    # n > 1_000_000 because of the previous block.
    if    ($n <         19471033) { @bases = ( 2,  299417); }
    elsif ($n <         38010307) { @bases = ( 2,  9332593); }
    elsif ($n <        316349281) { @bases = ( 11000544, 31481107); }
    elsif ($n <       4759123141) { @bases = ( 2, 7, 61); }
    elsif ($n <     154639673381) { @bases = ( 15, 176006322, 4221622697); }
    elsif ($n <   47636622961201) { @bases = ( 2, 2570940, 211991001, 3749873356); }
    elsif ($n < 3770579582154547) { @bases = ( 2, 2570940, 880937, 610386380, 4130785767); }
    else                          { @bases = ( 2, 325, 9375, 28178, 450775, 9780504, 1795265022); }
    return is_strong_pseudoprime($n, @bases)  ?  2  :  0;
  }

  # Inlined BPSW
  return 0 unless _miller_rabin_2($n);
  return is_almost_extra_strong_lucas_pseudoprime($n) ? 2 : 0;
}

sub is_prime {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n < 2;

  if (ref($n) eq 'Math::BigInt') {
    return 0 unless Math::BigInt::bgcd($n, B_PRIM235)->is_one;
  } else {
    if ($n < 7) { return ($n == 2) || ($n == 3) || ($n == 5) ? 2 : 0; }
    return 0 if !($n % 2) || !($n % 3) || !($n % 5);
  }
  return _is_prime7($n);
}

# is_prob_prime is the same thing for us.
*is_prob_prime = \&is_prime;

# BPSW probable prime.  No composites are known to have passed this test
# since it was published in 1980, though we know infinitely many exist.
# It has also been verified that no 64-bit composite will return true.
# Slow since it's all in PP and uses bigints.
sub is_bpsw_prime {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n < 2;
  return 0 unless _miller_rabin_2($n);
  if ($n <= 18446744073709551615) {
    return is_almost_extra_strong_lucas_pseudoprime($n) ? 2 : 0;
  }
  return is_extra_strong_lucas_pseudoprime($n) ? 1 : 0;
}

sub is_provable_prime {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n < 2;
  if ($n <= 18446744073709551615) {
    return 0 unless _miller_rabin_2($n);
    return 0 unless is_almost_extra_strong_lucas_pseudoprime($n);
    return 2;
  }
  my($is_prime, $cert) = Math::Prime::Util::is_provable_prime_with_cert($n);
  $is_prime;
}

# Possible sieve storage:
#   1) vec with mod-30 wheel:   8 bits  / 30
#   2) vec with mod-2 wheel :  15 bits  / 30
#   3) str with mod-30 wheel:   8 bytes / 30
#   4) str with mod-2 wheel :  15 bytes / 30
#
# It looks like using vecs is about 2x slower than strs, and the strings also
# let us do some fast operations on the results.  E.g.
#   Count all primes:
#      $count += $$sieveref =~ tr/0//;
#   Loop over primes:
#      foreach my $s (split("0", $$sieveref, -1)) {
#        $n += 2 + 2 * length($s);
#        .. do something with the prime $n
#      }
#
# We're using method 4, though sadly it is memory intensive relative to the
# other methods.  I will point out that it is 30-60x less memory than sieves
# using an array, and the performance of this function is over 10x that
# of naive sieves.

sub _sieve_erat_string {
  my($end) = @_;
  $end-- if ($end & 1) == 0;
  my $s_end = $end >> 1;

  my $whole = int( $s_end / 15);   # Prefill with 3 and 5 already marked.
  croak "Sieve too large" if $whole > 1_145_324_612;  # ~32 GB string
  my $sieve = '100010010010110' . '011010010010110' x $whole;
  substr($sieve, $s_end+1) = '';   # Ensure we don't make too many entries
  my ($n, $limit) = ( 7, int(sqrt($end)) );
  while ( $n <= $limit ) {
    for (my $s = ($n*$n) >> 1; $s <= $s_end; $s += $n) {
      substr($sieve, $s, 1) = '1';
    }
    do { $n += 2 } while substr($sieve, $n>>1, 1);
  }
  return \$sieve;
}

# TODO: this should be plugged into precalc, memfree, etc. just like the C code
{
  my $primary_size_limit = 15000;
  my $primary_sieve_size = 0;
  my $primary_sieve_ref;
  sub _sieve_erat {
    my($end) = @_;

    return _sieve_erat_string($end) if $end > $primary_size_limit;

    if ($primary_sieve_size == 0) {
      $primary_sieve_size = $primary_size_limit;
      $primary_sieve_ref = _sieve_erat_string($primary_sieve_size);
    }
    my $sieve = substr($$primary_sieve_ref, 0, ($end+1)>>1);
    return \$sieve;
  }
}


sub _sieve_segment {
  my($beg,$end,$limit) = @_;
  ($beg, $end) = map { _bigint_to_int($_) } ($beg, $end)
    if ref($end) && $end <= INTMAX;
  croak "Internal error: segment beg is even" if ($beg % 2) == 0;
  croak "Internal error: segment end is even" if ($end % 2) == 0;
  croak "Internal error: segment end < beg" if $end < $beg;
  croak "Internal error: segment beg should be >= 3" if $beg < 3;
  my $range = int( ($end - $beg) / 2 ) + 1;

  # Prefill with 3 and 5 already marked, and offset to the segment start.
  my $whole = int( ($range+14) / 15);
  my $startp = ($beg % 30) >> 1;
  my $sieve = substr('011010010010110', $startp) . '011010010010110' x $whole;
  # Set 3 and 5 to prime if we're sieving them.
  substr($sieve,0,2) = '00' if $beg == 3;
  substr($sieve,0,1) = '0'  if $beg == 5;
  # Get rid of any extra we added.
  substr($sieve, $range) = '';

  # If the end value is below 7^2, then the pre-sieve is all we needed.
  return \$sieve if $end < 49;

  my $sqlimit = Msqrtint($end);
  $limit = $sqlimit if !defined $limit || $sqlimit < $limit;
  # For large value of end, it's a huge win to just walk primes.

  my($p, $s, $primesieveref) = (7-2, 3, _sieve_erat($limit));
  my $sieve_end = ($end - $beg) >> 1;
  while ( (my $nexts = 1 + index($$primesieveref, '0', $s)) > 0 ) {
    $p += 2 * ($nexts - $s);
    $s = $nexts;
    my $p2 = $p*$p;

    if ($p2 < $beg) {  # Make p2 the next odd multiple of p >= beg
      if ($beg < 2**49) {
        my $f = int(($beg+$p-1)/$p);
        $p2 = $p * ($f + (1-($f&1)));
      } else {
        my $f = Mcdivint($beg,$p);
        $p2 = Mmulint($p, $f + (1-($f&1)));
      }
    }

    # Large bases and small segments often don't hit the segment at all.
    next if $p2 > $end;

    # Inner loop marking multiples of p, divide by 2 to keep loop simpler.
    for ($p2 = ($p2 - $beg) >> 1;  $p2 <= $sieve_end;  $p2 += $p) {
      substr($sieve, $p2, 1) = '1';
    }
  }
  \$sieve;
}

sub trial_primes {
  my($low,$high) = @_;
  if (!defined $high) {
    $high = $low;
    $low = 2;
  }
  validate_integer_nonneg($low);
  validate_integer_nonneg($high);
  return if $low > $high;
  my @primes;

  # For a tiny range, just use next_prime calls
  if (($high-$low) < 1000) {
    $low-- if $low >= 2;
    my $curprime = Mnext_prime($low);
    while ($curprime <= $high) {
      push @primes, $curprime;
      $curprime = Mnext_prime($curprime);
    }
    return \@primes;
  }

  # Sieve to 10k then BPSW test
  push @primes, 2  if ($low <= 2) && ($high >= 2);
  push @primes, 3  if ($low <= 3) && ($high >= 3);
  push @primes, 5  if ($low <= 5) && ($high >= 5);
  $low = 7 if $low < 7;
  $low++ if ($low % 2) == 0;
  $high-- if ($high % 2) == 0;
  my $sieveref = _sieve_segment($low, $high, 10000);
  my $n = $low-2;
  while ($$sieveref =~ m/0/g) {
    my $p = $n+2*pos($$sieveref);
    push @primes, $p if _miller_rabin_2($p) && is_extra_strong_lucas_pseudoprime($p);
  }
  return \@primes;
}

sub primes {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    validate_integer_nonneg($low);
    $low = 2 if $low < 2;
  } else {
    ($low,$high) = (2, $low);
  }
  validate_integer_nonneg($high);
  my $sref = [];
  return $sref if ($low > $high) || ($high < 2);
  return [grep { $_ >= $low && $_ <= $high } @_primes_small]
    if $high <= $_primes_small[-1];

  if ($Math::Prime::Util::_GMPfunc{"sieve_primes"} && $Math::Prime::Util::GMP::VERSION >= 0.34) {
    my @pr = Math::Prime::Util::GMP::sieve_primes($low, $high, 0);
    return ref($high) ? [maybetobigintall(@pr)] : \@pr;
  }

  # At some point even the pretty-fast pure perl sieve is going to be a
  # dog, and we should move to trials.  This is typical with a small range
  # on a large base.  More thought on the switchover should be done.
  return trial_primes($low, $high) if ref($low)  eq 'Math::BigInt'
                                   || ref($high) eq 'Math::BigInt'
                                   || ($low > 1_000_000_000_000 && ($high-$low) < int($low/1_000_000));

  push @$sref, 2  if ($low <= 2) && ($high >= 2);
  push @$sref, 3  if ($low <= 3) && ($high >= 3);
  push @$sref, 5  if ($low <= 5) && ($high >= 5);
  $low = 7 if $low < 7;
  $low++ if ($low % 2) == 0;
  $high-- if ($high % 2) == 0;
  return $sref if $low > $high;

  my($n,$sieveref);
  if ($low == 7) {
    $n = 0;
    $sieveref = _sieve_erat($high);
    substr($$sieveref,0,3,'111');
  } else {
    $n = $low-1;
    $sieveref = _sieve_segment($low,$high);
  }
  push @$sref, $n+2*pos($$sieveref)-1 while $$sieveref =~ m/0/g;
  $sref;
}

sub sieve_range {
  my($n, $width, $depth) = @_;
  validate_integer_nonneg($n);
  validate_integer_nonneg($width);
  validate_integer_nonneg($depth);

  my @candidates;
  my $start = $n;

  if ($n < 5) {
    push @candidates, (2-$n) if $n <= 2 && $n+$width-1 >= 2;
    push @candidates, (3-$n) if $n <= 3 && $n+$width-1 >= 3;
    push @candidates, (4-$n) if $n <= 4 && $n+$width-1 >= 4 && $depth < 2;
    $start = 5;
    $width -= ($start - $n);
  }

  return @candidates, map {$start+$_-$n } 0 .. $width-1 if $depth < 2;
  return @candidates, map { $_ - $n }
                      grep { ($_ & 1) && ($depth < 3 || ($_ % 3)) }
                      map { $start+$_ }
                      0 .. $width-1                     if $depth < 5;

  if (!($start & 1)) { $start++; $width--; }
  $width-- if !($width&1);
  return @candidates if $width < 1;

  my $sieveref = _sieve_segment($start, $start+$width-1, $depth);
  my $offset = $start - $n - 2;
  while ($$sieveref =~ m/0/g) {
    push @candidates, $offset + (pos($$sieveref) << 1);
  }
  return @candidates;
}

sub sieve_prime_cluster {
  my($lo,$hi,@cl) = @_;
  my $_verbose = getconfig()->{'verbose'};
  validate_integer_nonneg($lo);
  validate_integer_nonneg($hi);

  if ($Math::Prime::Util::_GMPfunc{"sieve_prime_cluster"}) {
    return maybetobigintall(
             Math::Prime::Util::GMP::sieve_prime_cluster($lo,$hi,@cl)
           );
  }

  return @{Mprimes($lo,$hi)} if scalar(@cl) == 0;

  unshift @cl, 0;
  for my $i (1 .. $#cl) {
    validate_integer_nonneg($cl[$i]);
    croak "sieve_prime_cluster: values must be even" if $cl[$i] & 1;
    croak "sieve_prime_cluster: values must be increasing" if $cl[$i] <= $cl[$i-1];
  }
  my($p,$sievelim,@p) = (17, 3000);
  if (defined $_BIGINT && (ref($lo) || ref($hi))) {
    ($lo,$hi) = map {tobigint($_)} ($lo,$hi) if ref($lo) ne $_BIGINT || ref($hi) ne $_BIGINT;
  }
  $p = 13 if ($hi-$lo) < 50_000_000;
  $p = 11 if ($hi-$lo) <  1_000_000;
  $p =  7 if ($hi-$lo) <     20_000 && $lo < INTMAX;

  # Add any cases under our sieving point.
  if ($lo <= $sievelim) {
    $sievelim = $hi if $sievelim > $hi;
    for my $n (@{Mprimes($lo,$sievelim)}) {
      my $ac = 1;
      for my $ci (1 .. $#cl) {
        if (!Mis_prime($n+$cl[$ci])) { $ac = 0; last; }
      }
      push @p, $n if $ac;
    }
    $lo = Mnext_prime($sievelim);
  }
  return @p if $lo > $hi;

  # Compute acceptable residues.
  my $pr = Mprimorial($p);
  my $startpr = _bigint_to_int($lo % $pr);

  my @acc = grep { ($_ & 1) && $_%3 }  ($startpr .. $startpr + $pr - 1);
  for my $c (@cl) {
    if ($p >= 7) {
      @acc = grep { (($_+$c)%3) && (($_+$c)%5) && (($_+$c)%7) } @acc;
    } else {
      @acc = grep { (($_+$c)%3)  && (($_+$c)%5) } @acc;
    }
  }
  for my $c (@cl) {
    @acc = grep { Mgcd($_+$c,$pr) == 1 } @acc;
  }
  @acc = map { $_-$startpr } @acc;

  print "cluster sieve using ",scalar(@acc)," residues mod $pr\n" if $_verbose;
  return @p if scalar(@acc) == 0;

  # Prepare table for more sieving.
  my @mprimes = @{Mprimes( $p+1, $sievelim)};
  my(@lorem,@vprem);
  for my $pidx (0..$#mprimes) {
    my $p = $mprimes[$pidx];
    $lorem[$pidx] = _bigint_to_int($lo % $p);
    for my $c (@cl) {
      $vprem[$pidx]->[ ($p-($c%$p)) % $p ] = 1;
    }
  }

  # Walk the range in primorial chunks, doing primality tests.
  my($nummr, $numlucas) = (0,0);
  while ($lo <= $hi) {

    my @racc = @acc;

    # Make sure we don't do anything past the limit
    if (($lo+$acc[-1]) > $hi) {
      my $max = _bigint_to_int($hi-$lo);
      @racc = grep { $_ <= $max } @racc;
    }

    # Sieve more values using native math
    for my $pidx (0 .. $#mprimes) {
      my $p = $mprimes[$pidx];
      my $rem = $lorem[$pidx];
      @racc = grep { !$vprem[$pidx]->[ ($rem+$_) % $p ] } @racc;
      last unless scalar(@racc);
    }

    # Do final primality tests.
    if ($lo < 1e13) {
      for my $r (@racc) {
        my($good, $p) = (1, $lo + $r);
        for my $c (@cl) {
          $nummr++;
          if (!Mis_prime($p+$c)) { $good = 0; last; }
        }
        push @p, $p if $good;
      }
    } else {
      for my $r (@racc) {
        my($good, $p) = (1, $lo + $r);
        for my $c (@cl) {
          $nummr++;
          if (!_miller_rabin_2($p+$c)) { $good = 0; last; }
        }
        next unless $good;
        for my $c (@cl) {
          $numlucas++;
          if (!Math::Prime::Util::is_extra_strong_lucas_pseudoprime($p+$c)) { $good = 0; last; }
        }
        push @p, $p if $good;
      }
    }

    $lo += $pr;
    if ($lo <= $hi) { # update native remainders
      $lorem[$_] = ($lorem[$_] + $pr) % $mprimes[$_] for 0..$#mprimes;
    }
  }
  print "cluster sieve ran $nummr MR and $numlucas Lucas tests\n" if $_verbose;
  @p;
}

sub prime_powers {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    validate_integer_nonneg($low);
    $low = 2 if $low < 2;
  } else {
    ($low,$high) = (2, $low);
  }
  validate_integer_nonneg($high);

  if ($high > 1e18 || ($high-$low) < 10) {
    my $sref = [];
    while ($low <= $high) {
      push @$sref, $low if Mis_prime_power($low);
      $low = Maddint($low, 1);
    }
    return $sref;
  } else {
    my @powers;
    for my $k (2 .. Mlogint($high,2)) {
      my $P = Mpowint(2,$k);
      push @powers, $P if $P >= $low;
    }
    for my $k (2 .. Mlogint($high,3)) {
      my $P = Mpowint(3,$k);
      push @powers, $P if $P >= $low;
    }
    for my $k (2 .. Mlogint($high,5)) {
      my $rootn = Mrootint($high, $k);
      Mforprimes( sub {
        my $P = Mpowint($_,$k);
        push @powers, $P if $P >= $low;
      }, 5, $rootn);
    }
    push @powers, @{Mprimes($low,$high)};
    return Mvecsorti(\@powers);
  }
}

sub twin_primes {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    validate_integer_nonneg($low);
    $low = 2 if $low < 2;
  } else {
    ($low,$high) = (2, $low);
  }
  validate_integer_nonneg($high);
  my @tp;
  if ($Math::Prime::Util::_GMPfunc{"twin_twin_primes"}) {
    @tp = Math::Prime::Util::GMP::sieve_twin_primes($low, $high);
  } else {
    @tp = sieve_prime_cluster($low, $high, 2);
  }
  return ref($high) ? [maybetobigintall(@tp)] : \@tp;
}

sub semi_primes {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    validate_integer_nonneg($low);
    $low = 4 if $low < 4;
  } else {
    ($low,$high) = (4, $low);
  }
  validate_integer_nonneg($high);
  my @sp;
  Math::Prime::Util::forsemiprimes(sub { push @sp,$_; }, $low, $high);
  \@sp;
}

# TODO: Port n_range_ramanujan_primes to replace this.
#       export it as a function
#
# For now, let's ignore it, this is only used for the PP.

sub _n_ramanujan_primes {
  my($n) = @_;
  return [] if $n <= 0;
  my $max = Mnth_prime_upper(int(48/19*$n)+1);
  my @L = (2, (0) x $n-1);
  my $s = 1;
  for (my $k = 7; $k <= $max; $k += 2) {
    $s++ if Mis_prime($k);
    $L[$s] = $k+1 if $s < $n;
    $s-- if ($k&3) == 1 && Mis_prime(($k+1)>>1);
    $L[$s] = $k+2 if $s < $n;
  }
  \@L;
}

sub ramanujan_primes {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    validate_integer_nonneg($low);
    $low = 2 if $low < 2;
  } else {
    ($low,$high) = (2, $low);
  }
  validate_integer_nonneg($high);
  return [] if ($low > $high) || ($high < 2);
  my $nn = Math::Prime::Util::prime_count_upper($high) >> 1;
  my $L = _n_ramanujan_primes($nn);
  shift @$L while @$L && $L->[0] < $low;
  pop @$L while @$L && $L->[-1] > $high;
  $L;
}

sub is_ramanujan_prime {
  my($n) = @_;
  return 1 if $n == 2;
  return 0 if $n < 11;
  my $L = Math::Prime::Util::ramanujan_primes($n,$n);
  return (scalar(@$L) > 0) ? 1 : 0;
}

sub nth_ramanujan_prime {
  my($n) = @_;
  validate_integer_nonneg($n);
  return undef if $n <= 0;  ## no critic qw(ProhibitExplicitReturnUndef)
  my $L = _n_ramanujan_primes($n);
  return $L->[$n-1];
}

sub next_prime {
  my($n) = @_;
  validate_integer_nonneg($n);
  return $_prime_next_small[$n] if $n <= $#_prime_next_small;
  # This turns out not to be faster.
  # return $_primes_small[1+_tiny_prime_count($n)] if $n < $_primes_small[-1];

  return tobigint(MPU_32BIT ? "4294967311" : "18446744073709551629") if !ref($n) && $n >= MPU_MAXPRIME;
  # n is now either 1) not bigint and < maxprime, or (2) bigint and >= uvmax

  if ($n > 4294967295 && getconfig()->{'gmp'}) {
    return reftyped($_[0], Math::Prime::Util::GMP::next_prime($n));
  }

  do {
    $n += $_wheeladvance30[$n%30];
  } while !($n%7) || !_is_prime7($n);

  $n;
}

sub prev_prime {
  my($n) = @_;
  validate_integer_nonneg($n);
  return (undef,undef,undef,2,3,3,5,5,7,7,7,7)[$n] if $n <= 11;
  if ($n > 4294967295 && getconfig()->{'gmp'}) {
    return reftyped($_[0], Math::Prime::Util::GMP::prev_prime($n));
  }

  do {
    $n -= $_wheelretreat30[$n%30];
  } while !($n%7) || !_is_prime7($n);

  $n = _bigint_to_int($n) if ref($n) && $n <= INTMAX;
  $n;
}

sub next_prime_power {
  my($n) = @_;
  validate_integer_nonneg($n);
  return (2,2,3,4,5,7,7,8,9)[$n] if $n <= 8;
  while (1) {
    $n = Maddint($n, 1);
    return $n if Mis_prime_power($n);
  }
}
sub prev_prime_power {
  my($n) = @_;
  validate_integer_nonneg($n);
  return (undef,undef,undef,2,3,4,5,5,7)[$n] if $n <= 8;
  while (1) {
    $n = Msubint($n, 1);
    return $n if Mis_prime_power($n);
  }
}

sub partitions {
  my $n = shift;
  validate_integer_nonneg($n);

  my $d = Msqrtint(Maddint($n,1));
  my @pent = (1, map { (($_*(3*$_+1))>>1, (($_+1)*(3*$_+2))>>1) } 1 .. $d);
  my $bigpn = (~0 > 4294967295) ? 400 : 270;
  my($ZERO,$ONE) = map { $n >= $bigpn ? tobigint($_) : $_ } (0,1);
  my @part = ($ONE);
  foreach my $j (scalar @part .. $n) {
    my ($psum1, $psum2) = ($ZERO, $ZERO);
    my $k = 1;
    foreach my $p (@pent) {
      last if $p > $j;
      if ((++$k) & 2) { $psum1 += $part[ $j - $p ] }
      else            { $psum2 += $part[ $j - $p ] }
    }
    $part[$j] = $psum1 - $psum2;
  }
  return $part[$n];
}

my @_lf63 = (0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,1,1,1,0,0,1,0,0,1,1,1,1,0,0,1,0,0,1,0,0,1,1,1,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0);
my @_small_lucky = (undef,1,3,7,9,13,15,21,25,31,33,37,43,49,51,63,67,69,73,75,79,87,93,99,105,111,115,127,129,133,135,141,151,159,163,169,171,189,193,195);

sub lucky_numbers {
  my($lo,$hi) = @_;
  if (defined $hi) { validate_integer_nonneg($lo); }
  else             { ($lo,$hi) = (1, $lo);         }
  validate_integer_nonneg($hi);
  return [] if $hi < $lo || $hi == 0;

  my @lucky;
  # This wheel handles the evens and every 3rd by a mod 6 wheel,
  # then uses the mask to skip every 7th and 9th remaining value.
  for (my $k = 1;  $k <= $hi;  $k += 6) {
    my $m63 = $k % 63;
    push @lucky, $k unless $_lf63[$m63];
    push @lucky, $k+2 unless $_lf63[$m63+2];
  }
  delete $lucky[-1] if $lucky[-1] > $hi;

  # Do the standard lucky sieve.
  for (my $k = 4; $k <= $#lucky && $lucky[$k]-1 <= $#lucky; $k++) {
    for (my $skip = my $index = $lucky[$k]-1;  $index <= $#lucky;  $index += $skip) {
      splice(@lucky, $index, 1);
    }
  }

  if ($lo > 1) { @lucky = grep { $_ >= $lo } @lucky; }

  \@lucky;
}

sub lucky_count {
  my($lo,$hi) = @_;
  if (defined $hi) { validate_integer_nonneg($lo); }
  else             { ($lo,$hi) = (1, $lo);         }
  validate_integer_nonneg($hi);
  return 0 if $hi < $lo || $hi == 0;

  # Return from our static data if very small.
  return scalar(grep { defined $_ && $_ >= $lo && $_ <= $hi } @_small_lucky) if $hi <= $_small_lucky[-1];

  # Trivial but slow way:
  # return scalar(@{Math::Prime::Util::lucky_numbers($lo, $hi)});

  $lo-- if $lo & 1;
  $hi++ if $hi & 1;
  my $lsize = 1 + lucky_count_upper($hi);
  my ($locount, $hicount) = ($lo >> 1, $hi >> 1);
  my $ln = Math::Prime::Util::lucky_numbers($lsize);
  shift @$ln;
  if ($lo <= 1) {
    $hicount -= int($hicount / $_) for @$ln;
  } else {
    for my $l (@$ln) {
      last if $l > $hicount;
      $locount -= int($locount / $l) if $l <= $lo;
      $hicount -= int($hicount / $l);
    }
  }
  return $hicount - $locount;
}
sub _simple_lucky_count_approx {
  my $n = shift;
  $n = "$n" if ref($n);
  return 0 + ($n > 0) + ($n > 2) if $n < 7;
  return 0.9957 * $n/log($n) if $n <= 1000000;
  return (1.03670 - log($n)/299) * $n/log($n);
}
sub _simple_lucky_count_upper {
  my $n = shift;
  $n = "$n" if ref($n);
  return 0 + ($n > 0) + ($n > 2) if $n < 7;
  return int(5 + 1.039 * $n/log($n)) if $n <= 7000;
  my $a = ($n < 10017000) ?   0.58003 - 3.00e-9 * ($n-7000)   : 0.55;
  return int($n/(1.065*log($n) - $a - 3.1/log($n) - 2.85/(log($n)*log($n))));
}
sub _simple_lucky_count_lower {
  my $n = shift;
  my $approx = _simple_lucky_count_approx($n);
  my $est = $approx * (($n <= 10000) ? 0.9 : 0.99);
  int($est);
}
sub lucky_count_approx {
  my $n = shift;
  validate_integer_nonneg($n);
  return scalar(grep { defined $_ && $_ <= $n } @_small_lucky) if $n <= $_small_lucky[-1];
  my($lo,$hi) = (_simple_lucky_count_lower($n), _simple_lucky_count_upper($n));
  _binary_search($n, $lo, $hi,
                 sub{Math::Prime::Util::nth_lucky_approx(shift)});
}
sub lucky_count_upper {
  my $n = shift;
  validate_integer_nonneg($n);
  return scalar(grep { defined $_ && $_ <= $n } @_small_lucky) if $n <= $_small_lucky[-1];
  my($lo,$hi) = (_simple_lucky_count_lower($n), _simple_lucky_count_upper($n));
  1+_binary_search($n, $lo, $hi,
                   sub{Math::Prime::Util::nth_lucky_lower(shift)});
}
sub lucky_count_lower {
  my $n = shift;
  validate_integer_nonneg($n);
  return scalar(grep { defined $_ && $_ <= $n } @_small_lucky) if $n <= $_small_lucky[-1];
  my($lo,$hi) = (_simple_lucky_count_lower($n), _simple_lucky_count_upper($n));
  _binary_search($n, $lo, $hi,
                 sub{Math::Prime::Util::nth_lucky_upper(shift)});
}

sub nth_lucky {
  my $n = shift;
  validate_integer_nonneg($n);
  return $_small_lucky[$n] if $n <= $#_small_lucky;
  my $k = $n-1;
  my $ln = lucky_numbers($n);
  shift @$ln;
  $k += int($k / ($_-1)) for reverse @$ln;
  2*$k+1;
}
sub nth_lucky_approx {
  my $n = shift;
  validate_integer_nonneg($n);
  return $_small_lucky[$n] if $n <= $#_small_lucky;
  $n = "$n" if ref($n);
  my($logn, $loglogn, $mult) = (log($n), log(log($n)), 1);
  if ($n <= 80000) {
    my $c = ($n <= 10000) ? 0.2502 : 0.2581;
    $mult = $logn + 0.5 * $loglogn + $c * $loglogn * $loglogn;
  } else {
    my $c = ($n <=    10000) ? -0.0173 : ($n <=   100000) ? -0.0318 :
            ($n <=  1000000) ? -0.0384 : ($n <= 10000000) ? -0.0422 : -0.0440;
    $mult = $logn + (0.5 + $c) * $loglogn *$loglogn;
  }
  return int( $n * $mult + 0.5);
}
sub nth_lucky_upper {
  my $n = shift;
  validate_integer_nonneg($n);
  return $_small_lucky[$n] if $n <= $#_small_lucky;
  my $c = ($n <= 100) ? 1.05 : ($n <= 300) ? 1.03 : ($n <= 800) ? 1.01 : 1.0033;
  return 1 + int( $c * nth_lucky_approx($n) + 0.5 );
}
sub nth_lucky_lower {
  my $n = shift;
  validate_integer_nonneg($n);
  return $_small_lucky[$n] if $n <= $#_small_lucky;
  my $c = ($n <= 130) ? 0.985 : ($n <= 1000) ? 0.992 : 0.996;
  return int( $c * nth_lucky_approx($n) );
}

sub is_lucky {
  my $n = shift;

  # Pretests
  return 0 if $n <= 0 || !($n % 2) || ($n % 6) == 5 || $_lf63[$n % 63];
  return 1 if $n < 45;

  # Really simple but slow:
  # return lucky_numbers($n)->[-1] == $n;

  my $upper = int(200 + 0.994 * $n / log($n));
  my $lucky = lucky_numbers($upper);
  my $pos = ($n+1) >> 1;
  my $i = 1;
  while (1) {
    my $l = ($i <= $#$lucky) ? $lucky->[$i++] : nth_lucky($i++);
    return 1 if $pos < $l;
    my $quo = int($pos / $l);
    return 0 if $pos == $quo * $l;
    $pos -= $quo;
  }
}

sub minimal_goldbach_pair {
  my $n = shift;
  validate_integer_nonneg($n);
  return undef if $n < 4;
  return Mis_prime($n-2) ? 2 : undef  if $n & 1 || $n == 4;
  my($p,$H)=(3,$n >> 1);
  while ($p <= $H) {
    return $p if Mis_prime($n-$p);
    $p = next_prime($p);
  }
  undef;
}
sub goldbach_pair_count {
  my $n = shift;
  validate_integer_nonneg($n);
  return 0 if $n < 4;
  return Mis_prime($n-2) ? 1 : 0  if $n & 1 || $n == 4;
  my $s = 0;
  Mforprimes( sub {
    $s++ if Mis_prime($n-$_);
  }, Mrshiftint($n,1), $n-3);
  $s;
}
sub goldbach_pairs {
  my $n = shift;
  return goldbach_pair_count($n) unless wantarray;
  validate_integer_nonneg($n);
  return () if $n < 4;
  return Mis_prime($n-2) ? (2) : ()  if $n & 1 || $n == 4;
  my @L;
  Mforprimes( sub {
    push @L,$n-$_ if Mis_prime($n-$_);
  }, Mrshiftint($n,1), $n-3);
  reverse @L;
}


sub primorial {
  my $n = shift;

  my @plist = @{Mprimes($n)};
  my $max = (MPU_32BIT) ? 29 : (OLD_PERL_VERSION) ? 43 : 53;

  # If small enough, multiply the small primes.
  if ($n < $max) {
    my $pn = 1;
    $pn *= $_ for @plist;
    return $pn;
  }

  # Otherwise, combine them as UVs, then combine using product tree.
  my $i = 0;
  while ($i < $#plist) {
    my $m = $plist[$i] * $plist[$i+1];
    if ($m <= INTMAX) { splice(@plist, $i, 2, $m); }
    else              { $i++;                      }
  }
  Mvecprod(@plist);
}

sub pn_primorial {
  my $n = shift;
  return (1,2,6,30,210,2310,30030,510510,9699690,223092870)[$n] if $n < 10;
  Mprimorial(nth_prime($n));
}

sub consecutive_integer_lcm {
  my $n = shift;
  validate_integer_nonneg($n);

  return (1,1,2)[$n] if $n <= 2;
  my @powers;
  for (my $p = 2; $p <= $n; $p = Mnext_prime($p)) {
    my($p_power, $pmin) = ($p, int($n/$p));
    $p_power = Mmulint($p_power,$p) while $p_power <= $pmin;
    push @powers, $p_power;
  }
  my $pn = Mvecprod(@powers);
  $pn = _bigint_to_int($pn) if $pn <= INTMAX;
  return $pn;
}

sub frobenius_number {
  my(@A) = @_;
  return undef if scalar(@A) == 0;
  validate_integer_positive($_) for @A;
  Mvecsorti(\@A);
  return -1 if $A[0] == 1;
  return undef if $A[0] <= 1 || scalar(@A) <= 1;
  croak "Frobenius number set must be coprime" unless Mgcd(@A) == 1;

  return Msubint(Msubint(Mmulint($A[0],$A[1]),$A[0]),$A[1]) if scalar(@A) == 2;

  # Basic Round Robin algorithm from Böcker and Lipták
  # https://bio.informatik.uni-jena.de/wp/wp-content/uploads/2024/01/BoeckerLiptak_FastSimpleAlgorithm_reprint_2007.pdf

  my $nlen = $A[0];
  my @N = (0, (undef) x ($nlen-1));
  for my $i (1 .. $#A) {
    { # Optimization 3, skip redundant bases
      my $ai = $A[$i];
      my $np = $N[Mmodint($ai,$nlen)];
      next if defined $np && $np <= $ai;
    }
    my $d = Mgcd($A[0], $A[$i]);
    my $nlend = Mdivint($nlen,$d);
    for my $r (0 .. $d-1) {
      my $n = ($r == 0) ? 0
            : Mvecmin(grep {defined} @N[map { $r+$_*$d } 0..$nlend]);
      if (defined $n) {
        if (Maddint($n,Mmulint($A[$i],$nlend-1)) <= INTMAX) {
          for (1 .. $nlend-1) {
            $n += $A[$i];
            my $p = $n % $nlen;
            if (!defined $N[$p] || $N[$p] >= $n) {$N[$p]=$n;} else {$n=$N[$p];}
          }
        } else {
          for (1 .. $nlend-1) {
            $n = Maddint($n,$A[$i]);
            my $p = Mmodint($n,$nlen);
            if (!defined $N[$p] || $N[$p] >= $n) {$N[$p]=$n;} else {$n=$N[$p];}
          }
        }
      }
    }
  }
  my $max = Mvecmax(grep { defined } @N);
  $max -= $nlen if defined $max;
  $max;
}

sub jordan_totient {
  my($k, $n) = @_;
  return ($n == 1) ? 1 : 0  if $k == 0;
  return Mtotient($n)       if $k == 1;
  return ($n == 1) ? 1 : 0  if $n <= 1;

  return reftyped($_[0], Math::Prime::Util::GMP::jordan_totient($k, $n))
    if $Math::Prime::Util::_GMPfunc{"jordan_totient"};

  my $totient = 1;
  foreach my $f (Mfactor_exp($n)) {
    my ($p, $e) = @$f;
    $p = Mpowint($p,$k);
    $totient = Mmulint($totient, $p-1);
    $totient = Mmulint($totient, $p) for 2 .. $e;
  }
  $totient;
}

sub euler_phi {
  return euler_phi_range(@_) if scalar @_ > 1;
  my($n) = @_;
  return 0 if defined $n && $n < 0;

  return reftyped($_[0],Math::Prime::Util::GMP::totient($n))
    if $Math::Prime::Util::_GMPfunc{"totient"};

  validate_integer_nonneg($n);
  return $n if $n <= 1;

  my ($t2, $tot) = (1,1);

  if ($n % 2 == 0) {
    my $totk = 0;
    while (($n % 4) == 0) { $n >>= 1;  $totk++; }
    $n >>= 1;
    $t2 = $totk < 32 ? 1 << $totk : Mlshiftint(1,$totk) if $totk > 0;
  }

  if ($n < INTMAX) {
    foreach my $f (Mfactor_exp($n)) {
      my ($p, $e) = @$f;
      $tot *= $p-1;
      $tot *= $p for 2 .. $e;
    }
   } else {
    foreach my $f (Mfactor_exp($n)) {
      my ($p, $e) = @$f;
      $tot = Mmulint($tot, $p-1);
      $tot = Mmulint($tot, $p) for 2 .. $e;
    }
   }
   Mmulint($t2, $tot);
}

sub inverse_totient {
  my($n) = @_;
  validate_integer_nonneg($n);

  return wantarray ? (1,2) : 2 if $n == 1;
  return wantarray ? () : 0 if $n < 1 || ($n & 1);

  if (Mis_prime($n >> 1)) {   # Coleman Remark 3.3 (Thm 3.1) and Prop 6.2
    my $np1 = Maddint($n,1);
    return wantarray ? () : 0                      if !Mis_prime($np1);
    return wantarray ? ($np1, Mmulint($np1,2)) : 2 if $n >= 10;
  }

  if (!wantarray) {
    my %r = ( 1 => 1 );
    Mfordivisors(sub { my $d = $_;
      my $p = Maddint($d,1);
      if (Mis_prime($p)) {
        my($dp,@sumi,@sumv) = ($d);
        for my $v (0 .. Mvaluation($n, $p)) {
          Mfordivisors(sub { my $d2 = $_;
            if (defined $r{$d2}) { push @sumi, Mmulint($d2,$dp); push @sumv, $r{$d2}; }
          }, Mdivint($n,$dp));
          $dp = Mmulint($dp,$p);
        }
        $r{ $sumi[$_] } += $sumv[$_]  for 0 .. $#sumi;
      }
    }, $n);
    return (defined $r{$n}) ? $r{$n} : 0;

  } else {

    # To save memory, we split this into two steps.

    my $_verbose = getconfig()->{'verbose'};
    my %r = ( 1 => [1] );
    my %needed = ( $n => 0 );
    my @DIVINFO;

    # 1. For each divisor from 1 .. n, track which values are needed.
    for my $d (divisors($n)) {
      my $p = Maddint($d,1);
      next unless Mis_prime($p);
      my @L;
      for my $v (0 .. Mvaluation($n, $p)) {
        my $pv = Mpowint($p, $v);
        my($dp,$pp) = map { Mmulint($_,$pv) } ($d,$p);
        Mfordivisors(sub { my $d2 = $_;
          my $F = Mmulint($d2,$dp);
          # In phase 2, we will look at the list in d2 to add to list in F.
          # If F isn't needed later then we ignore it completely.
          if (defined $needed{$F} && $needed{$F} < $d) {
            $needed{$d2} = $d unless defined $needed{$d2};
            push @L, [$d2,$pp,$F];
          }
        }, Mdivint($n, $dp));
      }
      push @DIVINFO, [$d, @L];
    }

    print "   ... inverse_totient phase 1 complete ...\n" if $_verbose;

    # 2. Process the divisors in reverse order.
    for my $dinfo (reverse @DIVINFO) {
      my($d,@L) = @$dinfo;
      my %todelete;
      my @T;
      # Multiply through by $pp
      for my $dset (@L) {
        if (defined $r{$dset->[0]}) {
          my($d2,$pp,$F) = @$dset;
          push @T, [$F, [map { Mmulint($pp,$_) } @{$r{$d2}}]];
          $todelete{$d2} = 1 if $needed{$d2} >= $d;
        }
      }
      # Delete intermediate data that isn't needed any more
      delete $r{$_} for keys %todelete;
      # Append the multiplied lists.
      push @{$r{$_->[0]}}, @{$_->[1]} for @T;
    }
    undef %needed;
    print "   ... inverse_totient phase 2 complete ...\n" if $_verbose;

    return (defined $r{$n}) ? @{Mvecsorti($r{$n})} : ();
  }
}

sub euler_phi_range {
  my($lo, $hi) = @_;
  validate_integer($lo);
  validate_integer($hi);

  my @totients;
  while ($lo < 0 && $lo <= $hi) {
    push @totients, 0;
    $lo++;
  }
  return @totients if $hi < $lo;

  if ($hi > 2**30 || $hi-$lo < 100) {
    ($lo,$hi) = (tobigint($lo),tobigint($hi)) if $hi > 2**49;
    push @totients, euler_phi($lo++)  while $lo <= $hi;
  } else {
    my @tot = (0 .. $hi);
    foreach my $i (2 .. $hi) {
      next unless $tot[$i] == $i;
      $tot[$i] = $i-1;
      foreach my $j (2 .. int($hi / $i)) {
        $tot[$i*$j] -= $tot[$i*$j]/$i;
      }
    }
    splice(@tot, 0, $lo) if $lo > 0;
    push @totients, @tot;
  }
  @totients;
}

sub _sumtot {
  my($n, $cdata, $ecache) = @_;
  return $cdata->[$n] if $n <= $#$cdata;
  return $ecache->{$n} if defined $ecache->{$n};

  my $sum = Mmulint($n, $n+1) >> 1;
  my $s = sqrtint($n);
  my $lim = Mdivint($n, $s+1);

  my($x, $nextx) = ($n, Mdivint($n,2));
  $sum -= Mmulint($x - $nextx, $cdata->[1]);
  for my $k (2 .. $lim) {
    ($x,$nextx) = ($nextx, Mdivint($n,$k+1));
    $sum -= ($x <= $#$cdata) ? $cdata->[$x] : _sumtot($x, $cdata, $ecache);
    $sum -= Mmulint($x - $nextx,
                    ($k <= $#$cdata) ? $cdata->[$k] : _sumtot($k, $cdata, $ecache));
  }
  if ($s > $lim) {
    ($x,$nextx) = ($nextx, Mdivint($n,$s+1));
    $sum -= Mmulint($x - $nextx,
                    ($s <= $#$cdata) ? $cdata->[$s] : _sumtot($s, $cdata, $ecache));
  }
  $ecache->{$n} = $sum;
  $sum;
}

sub sumtotient {
  my($n) = @_;
  validate_integer_nonneg($n);
  return $n if $n <= 2;

  if ($n < 900) {     # Simple linear sum for small values.
    my $sum = 0;
    $sum += $_ for Mtotient(1,$n);
    return $sum;
  }

  my $cbrt = Mrootint($n,3);
  my $csize = Mvecprod(4, $cbrt, $cbrt);
  $csize = 50_000_000 if $csize > 50_000_000;  # Limit memory use to ~2.5GB
  my @sumcache = Mtotient(0,$csize);
  $sumcache[$_] += $sumcache[$_-1] for 2 .. $csize;
  _sumtot($n, \@sumcache, {});
}


sub prime_bigomega {
  my($n) = @_;
  validate_integer_abs($n);
  return scalar(Mfactor($n));
}
sub prime_omega {
  my($n) = @_;
  validate_integer_abs($n);
  return scalar(Mfactor_exp($n));
}

sub moebius {
  return moebius_range(@_) if scalar @_ > 1;
  my($n) = @_;
  validate_integer_abs($n);
  return ($n == 1) ? 1 : 0  if $n <= 1;
  return 0 if ($n >= 49) && (!($n % 4) || !($n % 9) || !($n % 25) || !($n%49) );
  my @factors = Mfactor($n);
  foreach my $i (1 .. $#factors) {
    return 0 if $factors[$i] == $factors[$i-1];
  }
  return ((scalar @factors) % 2) ? -1 : 1;
}
sub is_square_free {
  return (Mmoebius($_[0]) != 0) ? 1 : 0;
}

sub is_odd {
  my($n) = @_;
  validate_integer($n);
  # Note:  If n is a Math::BigInt (Calc), then performance:
  #    0.25x  $n->is_odd()                           # method is fastest
  #    0.34x  (substr($n,-1,1) =~ tr/13579/13579/)   # Perl look at string
  #    0.46x  is_odd($n)                             # XS looks at the string
  #    1.0x   $n % 2 ? 1 : 0
  #    1.6x   $n & 1 : 1 : 0
  # Using LTM backend:
  #    0.21   $n->is_odd
  #    0.41   (substr($n,-1,1) =~ tr/13579/13579/)
  #    0.64   is_odd($n)
  #    0.9    $n & 1 ? 1 : 0
  #    1.0    $n % 2 ? 1 : 0
  #
  # Math::GMPz (30x faster baseline)
  #    0.23   Math::GMPz::Rmpz_odd_p($n)
  #    0.73   (substr($n,-1,1) =~ tr/13579/13579/)
  #    0.95   $n & 1 ? 1 : 0
  #    1.0    $n % 2 ? 1 : 0
  #    1.5    is_odd($n)
  my $R = ref($n);
  return $n->is_odd() if defined $R && $R eq 'Math::BigInt';
  return $n % 2 ? 1 : 0;
}
sub is_even {
  my($n) = @_;
  validate_integer($n);
  return $n % 2 ? 0 : 1;
}

sub is_divisible {
  my($n,@d) = @_;
  validate_integer_abs($n);
  for my $d (@d) {
    validate_integer_abs($d);
    if ($d == 0) { return 1 if $n == 0; }
    else         { return 1 if $n % $d == 0; }
  }
  0;
}
sub is_congruent {
  my($n,$c,$d) = @_;
  validate_integer($n);
  validate_integer($c);
  validate_integer_abs($d);
  if ($d != 0) {
    $n = Mmodint($n,$d) if $n < 0 || $n >= $d;
    $c = Mmodint($c,$d) if $c < 0 || $c >= $d;
  }
  return 0+($n == $c);
}

sub is_smooth {
  my($n, $k) = @_;
  validate_integer_abs($n);
  validate_integer_nonneg($k);

  return 1 if $n <= 1;
  return 0 if $k <= 1;
  return 1 if $n <= $k;

  return Math::Prime::Util::GMP::is_smooth($n,$k)
    if $Math::Prime::Util::_GMPfunc{"is_smooth"};

  if ($k <= 10000000) {
    my @f;
    while (1) {
      @f = Mtrial_factor($n, $k);
      last if scalar(@f) <= 1;
      return 0 if $f[-2] > $k;
      $n = $f[-1];
    }
    return 0 + ($f[0] <= $k);
  }

  return (Mvecnone(sub { $_ > $k }, Mfactor($n))) ? 1 : 0;
}
sub is_rough {
  my($n, $k) = @_;
  validate_integer_abs($n);
  validate_integer_nonneg($k);

  return 0+($k == 0) if $n == 0;
  return 1 if $n == 1 || $k <= 1;
  return 0 if $k > $n;
  return 0+($n >= 1) if $k == 2;

  return Math::Prime::Util::GMP::is_rough($n,$k)
    if $Math::Prime::Util::_GMPfunc{"is_rough"};

  if ($k < 50000) {
    my @f = Mtrial_factor($n, $k-1);
    return 0 + ($f[0] >= $k);
  }

  return (Mvecnone(sub { $_ < $k }, Mfactor($n))) ? 1 : 0;
}
sub is_powerful {
  my($n, $k) = @_;
  validate_integer($n);
  if (defined $k) { validate_integer_nonneg($k); } else { $k = 2; }

  return 0 if $n < 1;
  return 1 if $n == 1 || $k <= 1;

  return _gmpcall("is_powerful",$n,$k)
    if $Math::Prime::Util::_GMPfunc{"is_powerful"};

  # First quick checks for inadmissibility.
  if ($k == 2) {
    return 0 if ($n%3)  == 0 && ($n%9) != 0;
    return 0 if ($n%5)  == 0 && ($n%25) != 0;
    return 0 if ($n%7)  == 0 && ($n%49) != 0;
    return 0 if ($n%11) == 0 && ($n%121) != 0;
  } else {
    return 0 if ($n%3)  == 0 && ($n%27) != 0;
    return 0 if ($n%5)  == 0 && ($n%125) != 0;
    return 0 if ($n%7)  == 0 && ($n%343) != 0;
    return 0 if ($n%11) == 0 && ($n%1331) != 0;
  }

  # Next, check and remove all primes under 149 with three 64-bit gcds.
  for my $GCD ("614889782588491410","3749562977351496827","4343678784233766587") {
    my $g = Mgcd($n, $GCD);
    if ($g != 1) {
      # Check anything that divides n also divides k times (and remove).
      my $gk = Mpowint($g, $k);
      return 0 if ($n % $gk) != 0;
      $n = Mdivint($n, $gk);
      # Now remove any possible further amounts of these divisors.
      $g = Mgcd($n, $g);
      while ($n > 1 && $g > 1) {
        $n = Mdivint($n, $g);
        $g = Mgcd($n, $g);
      }
      return 1 if $n == 1;
    }
  }

  # For small primes, check for perfect powers and thereby limit the search
  # to divisibiilty conditions on primes less than n^(1/(2k)).  This is
  # usually faster than full factoring.
  #
  # But ... it's possible this will take far too long (e.g. n=2^256+1).  So
  # limit to something reasonable.

  return 1 if $n == 1 || Mis_power($n) >= $k;
  return 0 if $n < Mpowint(149, 2*$k);

  my $lim_actual = Mrootint($n, 2*$k);
  my $lim_effect = ($lim_actual > 10000) ? 10000 : $lim_actual;

  if ($Math::Prime::Util::_GMPfunc{"trial_factor"}) {
    while (1) {
      my @fac = Math::Prime::Util::GMP::trial_factor($n, $lim_effect);
      last if scalar(@fac) <= 1;
      my $f = $fac[0];
      my $fk = ($k==2) ? $f*$f : Mpowint($f,$k);
      return 0 if ($n % $fk) != 0;
      $n = Mdivint($n, $fk);
      $n = Mdivint($n, $f) while !($n % $f);
      return 1 if $n == 1 || Mis_power($n) >= $k;
      return 0 if $n < $fk*$fk;
    }
  } else {
    Mforprimes( sub {
      my $pk = ($k==2) ? $_*$_ : Mpowint($_,$k);
      Math::Prime::Util::lastfor(),return if $n < $pk*$pk;
      if (($n%$_) == 0) {
        Math::Prime::Util::lastfor(),return if ($n % $pk) != 0;
        $n = Mdivint($n, $pk);
        $n = Mdivint($n, $_) while ($n % $_) == 0;
        Math::Prime::Util::lastfor(),return if $n == 1 || Mis_power($n) >= $k;
      }
    }, 149, $lim_effect);
  }
  return 1 if $n == 1 || Mis_power($n) >= $k;
  return 0 if $n <= Mpowint($lim_effect, 2*$k);

  # Taking too long.  Factor what is left.
  return (Mvecall(sub { $_->[1] >= $k }, Mfactor_exp($n))) ? 1 : 0;
}

sub _powerful_count_recurse {
  my($n, $k, $m, $r) = @_;
  my $lim = Mrootint(Mdivint($n, $m), $r);

  return $lim if $r <= $k;

  my $sum = 0;
  for my $i (1 .. $lim) {
    if (Mgcd($m,$i) == 1 && Mis_square_free($i)) {
      $sum += _powerful_count_recurse($n, $k, Mmulint($m,Mpowint($i,$r)), $r-1);
    }
  }
  $sum;
}

sub powerful_count {
  my($n, $k) = @_;
  validate_integer($n);   $n = 0 if $n < 0;
  if (defined $k) { validate_integer_nonneg($k); } else { $k = 2; }

  return $n if $k <= 1 || $n <= 1;

  if ($k == 2) {
    my $sum = 0;
    for my $i (1 .. Mrootint($n,3)) {
      $sum += Msqrtint(Mdivint($n,Mpowint($i,3)))  if Mis_square_free($i);
    }
    return $sum;
  }

  _powerful_count_recurse($n, $k, 1, 2*$k-1);
}

sub nth_powerful {
  my($n, $k) = @_;
  validate_integer_nonneg($n);
  if (defined $k) { validate_integer_nonneg($k); } else { $k = 2; }

  return undef if $n == 0;
  return $n if $k <= 1 || $n <= 1;
  return Mpowint(2,$k) if $n == 2;
  return Mpowint(2,$k+1) if $n == 3;

  # For small n, we can generate k-powerful numbers rapidly.  But without
  # a reasonable upper limit, it's not clear how to effectively do it.
  # E.g. nth_powerful(100,60) = 11972515182562019788602740026717047105681

  my $lo = Mpowint(2, $k+1);
  my $hi = ~0;
  if ($k == 2) {
    $lo = int( $n*$n/4.72303430688484 + 0.3 * $n**(5/3) );
    $hi = int( $n*$n/4.72303430688484 + 0.5 * $n**(5/3) );  # for n >= 170
    $hi = ~0 if $hi > ~0;
    $lo = $hi >> 1 if $lo > $hi;
  }
  # We should use some power estimate here.

  # hi could be too low.
  while (Math::Prime::Util::powerful_count($hi,$k) < $n) {
    $lo = Maddint($hi,1);
    $hi = Mmulint($k, $hi);
  }

  # Simple binary search
  while ($lo < $hi) {
    my $mid = $lo + (($hi-$lo) >> 1);
    if (Math::Prime::Util::powerful_count($mid,$k) < $n) { $lo = $mid+1; }
    else                                                 { $hi = $mid; }
  }
  $hi;
}

# Not currently used
sub _genpowerful {
  my($m, $r, $n, $k, $arr) = @_;
  if ($r < $k) { push @$arr, $m; return; }
  my $rootdiv = Mrootint(Mdivint($n, $m), $r);
  if ($r == $k) {
    push @$arr, Mmulint($m, Mpowint($_,$k))  for 1 .. $rootdiv;
  } else {
    for my $i (1 .. $rootdiv) {
      if (Mgcd($m,$i) == 1 && Mis_square_free($i)) {
        _genpowerful(Mmulint($m, Mpowint($i,$r)), $r-1, $n, $k, $arr);
      }
    }
  }
}

sub _sumpowerful {
  my($m, $r, $n, $k) = @_;
  return $m if $r < $k;

  my $rootdiv = Mrootint(Mdivint($n, $m), $r);

  return Mmulint($m, Mpowersum($rootdiv, $k))  if $r == $k;

  # Generating then summing the list turns out to be MUCH faster than
  # summing as we go.
  my @v;
  for my $i (1 .. $rootdiv) {
    if (Mgcd($m,$i) == 1 && Mis_square_free($i)) {
      push @v, _sumpowerful(Mmulint($m, Mpowint($i,$r)), $r-1, $n, $k);
    }
  }
  Mvecsum(@v);
}

sub sumpowerful {
  my($n, $k) = @_;
  validate_integer($n);   $n = 0 if $n < 0;
  if (defined $k) { validate_integer_nonneg($k); } else { $k = 2; }

  return $n if $n <= 1;
  return Mrshiftint(Mmulint($n,Maddint($n,1)),1) if $k <= 1;

  # Alternate method for testing.
  # my @a;  _genpowerful(1, 2*$k-1, $n, $k, \@a);  return Mvecsum(@a);

  return _sumpowerful(1, 2*$k-1, $n, $k);
}


# Generate k-powerful numbers.  See Trizen, Feb 2020 and Feb 2024

sub _pcg {
  my($lo, $hi, $k, $m, $r, $pn) = @_;
  my($beg,$end) = (1, Mrootint(Mdivint($hi,$m), $r));

  if ($r <= $k) {
    if ($lo > $m) {
      my $lom = Mcdivint($lo,$m);
      if ( ($lom >> $r) == 0) {
        $beg = 2;
      } else {
        $beg = Mrootint($lom,$r);
        $beg++ if Mpowint($beg,$r) != $lom;
      }
    }
    push @$pn, $m * Mpowint($_,$r) for ($beg .. $end);
    return;
  }

  for my $v ($beg .. $end) {
    _pcg($lo, $hi, $k, $m * Mpowint($v,$r), $r-1, $pn)
       if Mgcd($m,$v) == 1 && Mis_square_free($v);
  }
}
sub powerful_numbers {
  my($lo, $hi, $k) = @_;
  if (defined $k) { validate_integer_nonneg($k); } else { $k = 2; }
  if (defined $hi) {
    validate_integer_nonneg($lo);
  } else {
    ($lo, $hi) = (1, $lo);
  }
  validate_integer_nonneg($hi);
  return [] if $hi < $lo;
  return [$lo .. $hi] if $k <= 1;

  my $npn = Math::Prime::Util::powerful_count($hi,$k);
  $npn -= Math::Prime::Util::powerful_count($lo-1, $k) if $lo > 1;

  my $pn = [];
  _pcg($lo, $hi, $k, 1, 2*$k-1, $pn);
  Mvecsorti($pn);
}

sub is_powerfree {
  my($n, $k) = @_;
  validate_integer_abs($n);
  if (defined $k) { validate_integer_nonneg($k); }
  else            { $k = 2; }

  return (($n == 1) ? 1 : 0)  if $k < 2 || $n <= 1;
  #return 1 if $n < Mpowint(2,$k);
  return 1 if $n < 4;

  if ($k == 2) {
    return 0 if !($n % 4) || !($n % 9) || !($n % 25);
    return 1 if $n < 49;   # 7^2
  } elsif ($k == 3) {
    return 0 if !($n % 8) || !($n % 27) || !($n % 125);
    return 1 if $n < 343;  # 7^3
  }

  # return (Mvecall(sub { $_->[1] < $k }, Mfactor_exp($n))) ? 1 : 0;
  for my $pe (Mfactor_exp($n)) {
    return 0 if $pe->[1] >= $k;
  }
  1;
}

sub powerfree_count {
  my($n, $k) = @_;
  validate_integer_abs($n);
  if (defined $k) { validate_integer_nonneg($k); }
  else            { $k = 2; }

  return (($n >= 1) ? 1 : 0)  if $k < 2 || $n <= 1;

  my $count = 0;
  my $nk = Mrootint($n, $k);

  # If we can do everything native, do that.
  if ($n < SINTMAX && $nk < 20000) {
    use integer;
    my @mu = Mmoebius(0, $nk);
    foreach my $i (2 .. $nk) {
      $count += $mu[$i] * $n/($i**$k) if $mu[$i];
    }
    return Maddint($count,$n);
  } elsif ($n < SINTMAX && $nk < 1e8) {
    # Split out the trailing n/i^k = 1, saves memory and time if large enough.
    use integer;
    my $L1 = Mrootint($n/2,$k);
    my @mu = Mmoebius(0, $L1);
    foreach my $i (2 .. $L1) {
      $count += $mu[$i] * $n/($i**$k) if $mu[$i];
    }
    #@mu = Mmoebius($L1+1, $nk);   my $c1 = 0;   $c1 += $_ for @mu;
    my $c1 = Math::Prime::Util::mertens($nk) - Math::Prime::Util::mertens($L1);
    return Mvecsum($count,$c1,$n);
  }

  # Simple way.  All the bigint math kills performance.
  # Math::Prime::Util::forsquarefree(
  #   sub {
  #     my $t = Mdivint($n, Mpowint($_, $k));
  #     $count = (scalar(@_) & 1) ? Msubint($count,$t) : Maddint($count,$t);
  #   },
  #   2, $nk
  # );

  # Optimization 1:  pull out all the ranges at the end with small constant
  #                  multiplications.
  # Optimization 2:  Use GMP basic arithmetic functions if possible, saving
  #                  all the bigint object overhead.  Can be 10x faster.

  my $A = Msqrtint($nk);
  my @L = (0, $nk, map { Mrootint(Mdivint($n,$_),$k) } 2..$A);
  my @C;

  Math::Prime::Util::forsquarefree(
    sub {
      $count = (scalar(@_) & 1)
             ? Ssubint($count, Sdivint($n, Spowint($_, $k)))
             : Saddint($count, Sdivint($n, Spowint($_, $k)));
    },
    2, $L[$A]
  );
  for my $i (2 .. $A) {
    my($c, $lo, $hi) = (0, $L[$i], $L[$i-1]);
    if ($i < 15) {
      $c = Math::Prime::Util::mertens($hi) - Math::Prime::Util::mertens($lo);
    } else {
      $c += $_ for Mmoebius( Maddint($lo,1), $hi );
    }
    push @C, $c * ($i-1);
    @C = (Mvecsum(@C)) if scalar(@C) > 100000;  # Save/restrict memory.
  }
  my $ctot = Mvecsum(@C); # Can typically be done in native math.
  Mvecsum($count, $n, $ctot);
}

sub nth_powerfree {
  my($n, $k) = @_;
  validate_integer_nonneg($n);
  if (defined $k) { validate_integer_nonneg($k); }
  else            { $k = 2; }

  return undef if $n == 0 || $k < 2;
  return $n if $n < 4;

  # 1. zm is the zeta multiplier (float), qk is the expected value (integer).
  my($zm, $qk);
  if ($n <= 2**52) {
    $zm = ($k == 2) ? 1.644934066848226 : 1.0 + RiemannZeta($k);
  } else {
    do { require Math::BigFloat; Math::BigFloat->import(); }
      if !defined $Math::BigFloat::VERSION;
    require Math::Prime::Util::ZetaBigFloat;
    my $acc = length("$n")+10;
    my $bk = Math::BigFloat->new($k);  $bk->accuracy($acc);
    $zm = Math::Prime::Util::ZetaBigFloat::RiemannZeta($bk)->badd(1)->numify;
  }
  my $verbose = getconfig()->{'verbose'};

  $qk = Mtoint($zm * "$n");
  print "nth_powerfree: zm $zm  qk $qk\n" if $verbose;

  my($count, $diff);
  # In practice this converges very rapidly, usually needing only one iteration.
  for (1 .. 10) {
    # 2. Get the actual count at qk and the difference from our goal.
    $count = Math::Prime::Util::powerfree_count($qk,$k);
    $diff = ($count >= $n)  ?  $count-$n  :  $n-$count;
    print "nth_powerfree: iter $_, count $count diff $diff\n" if $verbose;
    last if $diff <= 300;   # Threshold could be improved.

    # 3. If not close, update the estimate using the expected density zm.
    my $delta = Mtoint($zm * "$diff");
    $qk = $count > $n  ?  Msubint($qk,$delta)  :  Maddint($qk,$delta);
  }
  print "nth_powerfree: $qk, moving down to a powerfree number\n" if $verbose;

  # 4. Make sure we're on a powerfree number.
  $qk-- while !Math::Prime::Util::is_powerfree($qk,$k);
  print "nth_powerfree: $qk, need to move ",abs($n-$count)," steps\n" if $verbose;

  # 5. Walk forward or backward to next/prev powerfree number.
  my $adder = ($count < $n) ? 1 : -1;
  while ($count != $n) {
    do { $qk += $adder; } while !Math::Prime::Util::is_powerfree($qk,$k);
    $count += $adder;
  }
  $qk;
}

sub powerfree_sum {
  my($n, $k) = @_;
  validate_integer_nonneg($n);
  if (defined $k) { validate_integer_nonneg($k); }
  else            { $k = 2; }

  return (($n >= 1) ? 1 : 0)  if $k < 2 || $n <= 1;

  my $sum = 0;
  my($ik, $nik, $T);
  Math::Prime::Util::forsquarefree(
    sub {
      $ik = Mpowint($_, $k);
      $nik = Mdivint($n, $ik);
      $T = Mrshiftint(Mmulint($nik, Maddint($nik,1)), 1);
      $sum = (scalar(@_) & 1) ? Msubint($sum, Mmulint($ik,$T)) :
                                Maddint($sum, Mmulint($ik,$T));
    },
    Mrootint($n, $k)
  );
  $sum;
}

sub powerfree_part {
  my($n, $k) = @_;
  my $negmul = ($n < 0) ? -1 : 1;
  validate_integer_abs($n);
  if (defined $k) { validate_integer_nonneg($k); }
  else            { $k = 2; }

  return $negmul if $n == 1;
  return 0 if $k < 2 || $n == 0;

  #return Mvecprod(map { Mpowint($_->[0], $_->[1] % $k) } Mfactor_exp($n));

  # Rather than build with k-free section, we will remove excess powers
  my $P = $n;
  for my $pe (Mfactor_exp($n)) {
    $P = Mdivint($P, Mpowint($pe->[0], $pe->[1] - ($pe->[1] % $k)))
      if $pe->[1] >= $k;
  }
  $P = Mnegint($P) unless $negmul == 1;
  $P;
}

sub _T {
  my($n)=shift;
  Mdivint(Mmulint($n, Maddint($n, 1)), 2);
}
sub _fprod {
  my($n,$k)=@_;
  Mvecprod(map { 1 - Mpowint($_->[0], $k) } Mfactor_exp($n));
}

sub powerfree_part_sum {
  my($n, $k) = @_;
  validate_integer_abs($n);
  if (defined $k) { validate_integer_nonneg($k); }
  else            { $k = 2; }

  return (($n >= 1) ? 1 : 0)  if $k < 2 || $n <= 1;

  my $sum = _T($n);
  for (2 .. Mrootint($n,$k)) {
    $sum = Maddint($sum, Mmulint(_fprod($_,$k), _T(Mdivint($n, Mpowint($_, $k)))));
  }
  $sum;
}

sub squarefree_kernel {
  my($n) = @_;
  validate_integer($n);
  return -1 * Mlcm(Mfactor(-$n)) if $n < 0;
  Mlcm(Mfactor($n));
}

sub is_perfect_power {
  my($n) = @_;
  validate_integer($n);
  if ($n < 0) {
    my $res = Mis_power(Mnegint($n));
    return ($n == -1 || ($res > 2 && (($res & ($res-1)) != 0)))  ?  1  :  0;
  }
  return (1,1,0,0,1,0,0,0,1,1)[$n] if $n <= 9;
  return (Mis_power($n) > 1) ? 1 : 0;
}

sub _perfect_power_count {
  my($n) = @_;
  return 0+($n>=1)+($n>=4) if $n < 8;
  #return reftyped($_[0], Math::Prime::Util::GMP::perfect_power_count($n))
  #  if $Math::Prime::Util::_GMPfunc{"perfect_power_count"};
  my @T = (1);

  my $log2n = Mlogint($n,2);
  for my $k (2 .. $log2n) {
    my $m = Mmoebius($k);
    next if $m == 0;
    push @T, Mmulint(-$m, Msubint(Mrootint($n,$k),1));
  }
  Mvecsum(@T);
}
sub perfect_power_count {
  my($lo,$hi) = @_;
  if (defined $hi) { validate_integer_nonneg($lo); }
  else             { ($lo,$hi) = (1, $lo);         }
  validate_integer_nonneg($hi);
  return 0 if $hi < $lo || $hi == 0;
  return _perfect_power_count($hi) - (($lo <= 1) ? 0 : _perfect_power_count($lo-1));
}

sub perfect_power_count_approx {
  my($n) = @_;
  validate_integer_nonneg($n);
  _perfect_power_count($n);
}
sub perfect_power_count_lower {
  my($n) = @_;
  validate_integer_nonneg($n);
  _perfect_power_count($n);
}
sub perfect_power_count_upper {
  my($n) = @_;
  validate_integer_nonneg($n);
  _perfect_power_count($n);
}

sub _next_perfect_power {
  my($n, $only_oddpowers) = @_;
  croak "_npp must have positive n" if $n < 0;

  return 1 if $n == 0;
  return ($only_oddpowers ? 8 : 4) if $n == 1;

  my $log2n = Mlogint($n,2);
  my $kinit = $only_oddpowers ? 3 : 2;
  my $kinc  = $only_oddpowers ? 2 : 1;

  my $best = Mpowint(Maddint(Mrootint($n,$kinit),1),$kinit);
  for (my $k = $kinit+$kinc; $k <= 1+$log2n; $k += $kinc) {
    my $r = Mrootint($n,$k);
    my $c = Mpowint(Maddint($r,1),$k);
    $best = $c if $c < $best && $c > $n;
  }
  $best;
}
sub _prev_perfect_power {
  my($n, $only_oddpowers) = @_;
  croak "_ppp must have positive n" if $n < 0;

  return 0 + ($n>1) - ($n==0) if $n <= 4;
  return $only_oddpowers ? 1 : 4 if $n <= 8;

  my $log2n = Mlogint($n,2);
  my $kinit = $only_oddpowers ? 3 : 2;
  my $kinc  = $only_oddpowers ? 2 : 1;

  my $best = 8;
  for (my $k = $kinit; $k <= $log2n; $k += $kinc) {
    my $r = Mrootint($n,$k);
    if ($r > 1) {
      my $c = Mpowint($r,$k);
      $c = Mpowint(Msubint($r,1),$k) if $c >= $n;
      $best = $c if $c > $best && $c < $n;
    }
  }
  $best;
}

sub next_perfect_power {
  my($n) = @_;
  validate_integer($n);

  return 0 + ($n>=0) - ($n<-1)  if $n < 1 && $n >= -4;

  return Mnegint( _prev_perfect_power( Mnegint($n), 1 ) )  if $n < 0;
  _next_perfect_power($n, 0);
}

sub prev_perfect_power {
  my($n) = @_;
  validate_integer($n);

  return 0 + ($n>1) - ($n==0)  if $n <= 4 && $n >= 0;

  return Mnegint( _next_perfect_power( Mnegint($n), 1 ) )  if $n < 0;
  _prev_perfect_power($n, 0);
}


sub nth_perfect_power_approx {
  my($n) = @_;
  validate_integer_nonneg($n);
  return (undef,1,4,8,9,16,25,27)[$n] if $n < 8;

  # See  https://www.emis.de/journals/JIS/VOL15/Jakimczuk/jak29.pdf
  # See  https://www.researchgate.net/publication/268998744_Sums_of_perfect_powers

  # This is more accurate and about 200x faster than using BigFloat.
  if ($n > 2**32 && $Math::Prime::Util::_GMPfunc{"powreal"}) {
    *Gaddreal = \&Math::Prime::Util::GMP::addreal;
    *Gmulreal = \&Math::Prime::Util::GMP::mulreal;
    *Gpowreal = \&Math::Prime::Util::GMP::powreal;
    my $d = 2 * length($n) + 2;
    my $pp = Gmulreal($n,$n,$d);
    $pp = Gaddreal($pp, Gmulreal(13/3 ,Gpowreal($n, 4/3 ,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal(32/15,Gpowreal($n,16/15,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal(-2   ,Gpowreal($n, 5/3 ,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal(-2   ,Gpowreal($n, 7/5 ,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal(-2   ,Gpowreal($n, 9/7 ,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal( 2   ,Gpowreal($n,12/10,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal(-2   ,Gpowreal($n,13/11,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal(-2   ,Gpowreal($n,15/13,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal( 2   ,Gpowreal($n,16/14,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal( 2   ,Gpowreal($n,17/15,$d),$d),$d);
    $pp = Gaddreal($pp, Gmulreal(-0.48,Gpowreal($n,19/17,$d),$d),$d);
    $pp = Gaddreal($pp, -1.5,$d);
    $pp =~ s/\..*//;
    return Mtoint("$pp");
  }

  # Without this upgrade, it will return non-integers.
  $n = _upgrade_to_float($n) if $n > 2**32;

  if (!ref($n)) {
    my $pp = $n*$n  +  (13/3)*$n**(4/3)  +  (32/15)*$n**(16/15);
    $pp += -2*$n**( 5/ 3) + -2*$n**( 7/ 5);
    $pp += -2*$n**( 9/ 7) +  2*$n**(12/10);
    $pp += -2*$n**(13/11) + -2*$n**(15/13);
    $pp +=  2*$n**(16/14) +  2*$n**(17/15);
    $pp -= 0.48*$n**(19/17);
    return Mtoint($pp - 1.5);
  }

  # Taking roots is very expensive with Math::BigFloat, so minimize.
  my $n143 = $n->copy->broot(143);
  my $n105 = $n->copy->broot(105);

  my $n15 = $n105->copy->bpow(7);
  my $n13  = $n143->copy->bpow(11);
  my $n11  = $n143->copy->bpow(13);
  my $n7  = $n105->copy->bpow(15);
  my $n5  = $n105->copy->bpow(21);
  my $n3  = $n105->copy->bpow(35);

  my $pp = $n*$n  +  (13/3)*$n*$n3 +  (32/15)*$n*$n15;
  $pp += -2*$n*$n3**2   + -2*$n*$n5**2;
  $pp += -2*$n*$n7**2   +  2*$n*$n5;
  $pp += -2*$n*$n11**2  + -2*$n*$n13**2;
  $pp +=  2*$n*$n7      +  2*$n*$n15**2;
  $pp -= 0.48*$n*$n143**16.82352941176470588;  # close to 2/17
  $pp -= 1.5;
  $pp = $pp->as_int();
  Mtoint($pp);
}

sub nth_perfect_power_lower {
  my($n) = @_;
  validate_integer_nonneg($n);
  return (undef,1,4,8,9,16,25,27)[$n] if $n < 8;
  $n = _upgrade_to_float($n) if ref($n) || $n > 2**32;

  my $pp = $n*$n  +  (13/3)*$n**(4/3)  +  (32/15)*$n**(16/15);
  $pp += -2*$n**( 5/ 3) + -2*$n**( 7/ 5);
  $pp += -2*$n**( 9/ 7) +  2*$n**(12/10);
  $pp += -2*$n**(13/11) + -2*$n**(15/13);
  $pp += 1.5;
  Mtoint($pp);
}
sub nth_perfect_power_upper {
  my($n) = @_;
  validate_integer_nonneg($n);
  return (undef,1,4,8,9,16,25,27)[$n] if $n < 8;
  $n = _upgrade_to_float($n) if ref($n) || $n > 2**32;

  my $pp = $n*$n  +  (13/3)*$n**(4/3)  +  (32/15)*$n**(16/15);
  $pp += -2*$n**( 5/ 3) + -2*$n**( 7/ 5);
  $pp += -2*$n**( 9/ 7) +  2*$n**(12/10);
  $pp +=  2*$n**(16/14);
  $pp -= 3.5;
  Mtoint($pp);
}

sub nth_perfect_power {
  my($n) = @_;
  validate_integer_nonneg($n);
  return (undef,1,4,8,9,16,25,27)[$n] if $n < 8;
  my($g,$c,$apn,$gn);

  $gn = 1;
  $g = $apn = nth_perfect_power_approx($n);
  $c = _perfect_power_count($g);
  while ($n != $c && abs($n-$c) > 1000) {
    $g += $apn - nth_perfect_power_approx($c);
    $c = _perfect_power_count($g);
    last if $gn++ >= 20;
  }
  if ($c >= $n) {
    for ($g = Math::Prime::Util::prev_perfect_power($g+1);  $c > $n;  $c--) {
      $g = Math::Prime::Util::prev_perfect_power($g);
    }
  } else {
    for ( ; $c < $n; $c++) {
      $g = Math::Prime::Util::next_perfect_power($g);
     }
  }
  $g;
}

sub _prime_power_count {
  my($n) = @_;
  return (0,0,1,2,3,4)[$n] if $n <= 5;
  Mvecsum(
    map { Mprime_count( Mrootint($n, $_)) }  1 .. Mlogint($n,2)
  );
}
sub prime_power_count {
  my($lo,$hi) = @_;
  if (defined $hi) { validate_integer_nonneg($lo); }
  else             { ($lo,$hi) = (2, $lo);         }
  validate_integer_nonneg($hi);
  return 0 if $hi < $lo || $hi == 0;
  return _prime_power_count($hi) - (($lo <= 2) ? 0 : _prime_power_count($lo-1));
}
sub prime_power_count_lower {
  my $n = shift;
  validate_integer_nonneg($n);
  return (0,0,1,2,3,4)[$n] if $n <= 5;
  Mvecsum(
    map { Math::Prime::Util::prime_count_lower( Mrootint($n, $_)) }  1 .. Mlogint($n,2)
  );
}
sub prime_power_count_upper {
  my $n = shift;
  validate_integer_nonneg($n);
  return (0,0,1,2,3,4)[$n] if $n <= 5;
  Mvecsum(
    map { Math::Prime::Util::prime_count_upper( Mrootint($n, $_)) }  1 .. Mlogint($n,2)
  );
}
sub prime_power_count_approx {
  my $n = shift;
  validate_integer_nonneg($n);
  return (0,0,1,2,3,4)[$n] if $n <= 5;
  Mvecsum(
    map { Math::Prime::Util::prime_count_approx( Mrootint($n, $_)) }  1 .. Mlogint($n,2)
  );
}

sub _simple_nth_prime_power_upper {
  my $n = shift;
  Mnth_prime_upper($n);
}
sub _simple_nth_prime_power_lower {
  my $n = shift;
  return nth_prime_lower(int(0.65*$n)) if $n < 90;
  int( 0.98 * Math::Prime::Util::nth_prime_lower($n) - 400 );
}
sub nth_prime_power_lower {
  my $n = shift;
  validate_integer_nonneg($n);
  return (undef,2,3,4,5,7,8,9)[$n] if $n < 8;
  my($lo,$hi) = (_simple_nth_prime_power_lower($n), _simple_nth_prime_power_upper($n));
  _binary_search($n, $lo, $hi,
                 sub{Math::Prime::Util::prime_power_count_upper(shift)});
}
sub nth_prime_power_upper {
  my $n = shift;
  validate_integer_nonneg($n);
  return (undef,2,3,4,5,7,8,9)[$n] if $n < 8;
  my($lo,$hi) = (_simple_nth_prime_power_lower($n), _simple_nth_prime_power_upper($n));
  1+_binary_search($n, $lo, $hi,
                   sub{Math::Prime::Util::prime_power_count_lower(shift)});
}
sub nth_prime_power_approx {
  my $n = shift;
  validate_integer_nonneg($n);
  return (undef,2,3,4,5,7,8,9)[$n] if $n < 8;
  my($lo,$hi) = (_simple_nth_prime_power_lower($n), _simple_nth_prime_power_upper($n));
  _binary_search($n, $lo, $hi,
                 sub{Math::Prime::Util::prime_power_count_approx(shift)});
}
sub nth_prime_power {
  my $n = shift;
  validate_integer_nonneg($n);
  return (undef,2,3,4,5,7,8,9)[$n] if $n < 8;
  # TODO: This is a good candidte for the approx interpolation method
  my($lo,$hi) = (_simple_nth_prime_power_lower($n), _simple_nth_prime_power_upper($n));
  1+_binary_search($n, $lo, $hi,
                   sub{Math::Prime::Util::prime_power_count(shift)});
}


sub smooth_count {
  my($n, $k) = @_;
  return 0 if $n < 1;
  return 1 if $k <= 1;
  return $n if $k >= $n;

  my $sum = 1 + Mlogint($n,2);
  if ($k >= 3) {
    my $n3 = Mdivint($n, 3);
    while ($n3 > 3) {
      $sum += 1 + Mlogint($n3,2);
      $n3 = Mdivint($n3, 3);
    }
    $sum += $n3;
  }
  if ($k >= 5) {
    my $n5 = Mdivint($n, 5);
    while ($n5 > 5) {
      $sum += 1 + Mlogint($n5,2);
      my $n3 = Mdivint($n5, 3);
      while ($n3 > 3) {
        $sum += 1 + Mlogint($n3,2);
        $n3 = Mdivint($n3, 3);
      }
      $sum += $n3;
      $n5 = Mdivint($n5, 5);
    }
    $sum += $n5;
  }
  my $p = 7;
  while ($p <= $k) {
    my $np = Mdivint($n, $p);
    $sum += ($p >= $np) ? $np : Math::Prime::Util::smooth_count($np, $p);
    $p = Mnext_prime($p);
  }
  $sum;
}

sub rough_count {
  my($n, $k) = @_;
  return $n if $k <= 2;
  return $n-($n>>1) if $k <= 3;
  Math::Prime::Util::legendre_phi($n, Mprime_count($k-1));
}


# Recursive almost primes from Trizen.
sub _genkap {
  my($A, $B, $k, $m, $p, $cb) = @_;
  if ($k == 1) {
    Mforprimes( sub {
      $cb->(Mmulint($m, $_));
    }, Mvecmax($p, Mcdivint($A, $m)), Mdivint($B, $m));
  } else {
    my $s = Mrootint(Mdivint($B, $m), $k);
    while ($p <= $s) {
      my $t = mulint($m, $p);
      _genkap($A, $B, $k-1, $t, $p, $cb)
        if Mcdivint($A, $t) <= Mdivint($B, $t);  # Faster for tight ranges
      $p = next_prime($p);
    }
  }
}

sub _generate_almost_primes {
  my($A, $B, $k, $cb) = @_;
  $A = Mvecmax($A, Mpowint(2, $k));
  _genkap($A, $B, $k, 1, 2, $cb)  if $A <= $B;
}


sub almost_primes {
  my($k, $low, $high) = @_;
  validate_integer_nonneg($k);
  if (defined $high) { validate_integer_nonneg($low); }
  else               { ($low,$high) = (1, $low);      }
  validate_integer_nonneg($high);

  if ($k == 0) { return ($low <= 1 && $high >= 1) ? [1] : [] }
  if ($k == 1) { return Mprimes($low,$high); }
  # Don't call this, we could end up back here
  #if ($k == 2) { return Math::Prime::Util::semi_primes($low,$high); }

  my $minlow = Mpowint(2,$k);
  $low = $minlow if $low < $minlow;
  return [] if $low > $high;
  my @ap;

  if ($low > 1e9) {
    #while ($low <= $high) {
    #  push @ap, $low if is_almost_prime($k, $low);
    #  $low = add1int($low);
    #}
    Math::Prime::Util::forfactored(sub { push @ap,$_ if scalar(@_) == $k }, $low, $high);
    return \@ap;
  }

  _generate_almost_primes($low, $high, $k, sub { push @ap,$_[0]; });
  Mvecsorti(\@ap);
}

sub _rec_omega_primes {
  my($k, $lo, $hi, $m, $p, $opl) = @_;
  my $s = Mrootint(Mdivint($hi, $m), $k);
  foreach my $q (@{Mprimes($p, $s)}) {
    next if $m % $q == 0;
    for (my $v = Mmulint($m, $q); $v <= $hi ; $v = Mmulint($v, $q)) {
      if ($k == 1) {
        push @$opl, $v  if $v >= $lo;
      } else {
        _rec_omega_primes($k-1,$lo,$hi,$v,$q,$opl)  if Mmulint($v,$q) <= $hi;
      }
    }
  }
}

sub omega_primes {
  my($k, $low, $high) = @_;
  validate_integer_nonneg($k);
  if (defined $high) { validate_integer_nonneg($low); }
  else               { ($low,$high) = (1, $low);      }
  validate_integer_nonneg($high);

  if ($k == 0) { return ($low <= 1 && $high >= 1) ? [1] : [] }
  if ($k == 1) { return Math::Prime::Util::prime_powers($low,$high); }

  $low = Mvecmax($low, Mpn_primorial($k));
  return [] if $low > $high;

  my $opl = [];

  if ($high-$low > 1000000000 || ($k >= 10 && $high-$low > 10000000)) {
    # Recursive method from trizen
    _rec_omega_primes($k, $low, $high, 1, 2, $opl);
    Mvecsorti($opl);
  } else {
    # Simple iteration
    while ($low <= $high) {
      push @$opl, $low if Mprime_omega($low) == $k;
      $low++;
    }
  }
  $opl;
}

sub is_semiprime {
  my($n) = @_;
  validate_integer($n);
  return 0+($n == 4) if $n < 6;
  if ($n > 15) {
    return 0 if ($n %  4) == 0 || ($n %  6) == 0 || ($n %  9) == 0
             || ($n % 10) == 0 || ($n % 14) == 0 || ($n % 15) == 0;
  }
  return (Math::Prime::Util::is_prob_prime($n>>1) ? 1 : 0) if ($n % 2) == 0;
  return (Math::Prime::Util::is_prob_prime($n/3)  ? 1 : 0) if ($n % 3) == 0;
  return (Math::Prime::Util::is_prob_prime($n/5)  ? 1 : 0) if ($n % 5) == 0;

if (0) {  # TODO:  This is all REALLY slow without GMP
  # TODO: Something with GMP.  If nothing else, just factor.
  {
    my @f = trial_factor($n, 4999);
    return 0 if @f > 2;
    return (_is_prime7($f[1]) ? 1 : 0) if @f == 2;
  }
  return 0 if _is_prime7($n);
  {
    my @f = pminus1_factor ($n, 250_000);
    return 0 if @f > 2;
    return (_is_prime7($f[1]) ? 1 : 0) if @f == 2;
  }
  {
    my @f = pbrent_factor ($n, 128*1024, 3, 1);
    return 0 if @f > 2;
    return (_is_prime7($f[1]) ? 1 : 0) if @f == 2;
  }
}
  return (scalar(Mfactor($n)) == 2) ? 1 : 0;
}

sub is_almost_prime {
  my($k, $n) = @_;
  validate_integer_nonneg($k);
  validate_integer($n);
  return 0 if $n <= 0;

  return 0+($n==1) if $k == 0;
  return (Mis_prime($n) ? 1 : 0) if $k == 1;
  return Mis_semiprime($n) if $k == 2;
  return 0 if ($n >> $k) == 0;

  # TODO: Optimization here
  if (0) {  # This seems to just be slower
    while ($k > 0 && !($n % 2)) { $k--;  $n >>= 1; }
    while ($k > 0 && !($n % 3)) { $k--;  $n /= 3; }
    while ($k > 0 && !($n % 5)) { $k--;  $n /= 5; }
    while ($k > 0 && !($n % 7)) { $k--;  $n /= 7; }
    return 0+($n == 1) if $k == 0;
    return (Mis_prime($n) ? 1 : 0) if $k == 1;
    return Mis_semiprime($n) if $k == 2;
    return 0 if $n < Mpowint(11,$k);
  }

  return (scalar(Mfactor($n)) == $k) ? 1 : 0;
}
sub is_chen_prime {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n < 2;
  my $n2 = Maddint($n,2);
  return (Mis_prime($n) && (Mis_prime($n2) || Mis_semiprime($n2)));
}
sub next_chen_prime {
  my($n) = @_;
  validate_integer_nonneg($n);
  $n = Mnext_prime($n);
  while (1) {
    my $n2 = Maddint($n,2);
    return $n if Mis_prime($n2) || Mis_semiprime($n2);
    $n = Mnext_prime($n2);
  }
}

sub is_omega_prime {
  my($k, $n) = @_;
  validate_integer_nonneg($k);
  validate_integer($n);
  return 0 if $n <= 0;

  return 0+($n==1) if $k == 0;

  return (Mprime_omega($n) == $k) ? 1 : 0;
}

sub is_practical {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n <= 0;

  return $n==1?1:0 if $n % 2;
  return 1 if ($n & ($n-1)) == 0;
  return 0 if ($n % 6) && ($n % 20) && ($n % 28) && ($n % 88) && ($n % 104) && ($n % 16);

  my $prod = 1;
  my @pe = Mfactor_exp($n);
  for my $i (1 .. $#pe) {
    my($f,$e) = @{$pe[$i-1]};
    my $fmult = $f + 1;
    if ($e >= 2) {
      my $pke = $f;
      for (2 .. $e) {
        $pke = Mmulint($pke, $f);
        $fmult = Maddint($fmult, $pke);
      }
    }
    $prod = Mmulint($prod, $fmult);
    return 0 if $pe[$i]->[0] > (1 + $prod);
  }
  1;
}

sub is_delicate_prime {
  my($n, $b) = @_;
  validate_integer_nonneg($n);
  if (defined $b) { validate_integer_nonneg($b); } else { $b = 10; }
  croak "is_delicate_prime base must be >= 2" if $b < 2;

  return 0 if $b == 10 && $n < 100;   # Easy shown.
  return 1 if $b ==  3 && $n == 2;
  return 0 unless Mis_prime($n);

  if ($b == 10) {
    # String replacement method.  Faster in Perl.
    my $ndigits = length($n);
    for my $d (0 .. $ndigits-1) {
      my $N = "$n";
      my $dold = substr($N,$d,1);
      for my $dnew (0 .. 9) {
        next if $dnew == $dold;
        if ($d == 0 && $dnew == 0) {  # Leading zeros
          (my $T = substr($N,1)) =~ s/^0*//;
          return 0 if Mis_prime($T);
        } else {
          substr($N,$d,1) = $dnew;
          return 0 if Mis_prime($N);
        }
      }
    }
  } else {
    # Using todigitstring is slightly faster for bases < 10, but this is
    # decent and works for all 32-bit bases.
    # This is faster than Stamm's algorithm (in Perl, for possible bigints).
    my $D = [Mtodigits($n, $b)];
    for my $d (0 .. $#$D) {
      my $dold = $D->[$d];
      for my $dnew (0 .. $b-1) {
        next if $dnew == $dold;
        $D->[$d] = $dnew;
        return 0 if Mis_prime(Mfromdigits($D,$b));
      }
      $D->[$d] = $dold;
    }
  }
  1;
}

sub _totpred {
  my($n, $maxd) = @_;
  return 0 if $maxd <= 1 || Mis_odd($n);
  return 1 if ($n & ($n-1)) == 0;
  $n >>= 1;
  return 1 if $n == 1 || ($n < $maxd && Mis_prime(2*$n+1));
  for my $d (Mdivisors($n)) {
    last if $d >= $maxd;
    my $p = ($d < (INTMAX >> 2))  ?  ($d << 1) + 1 :
            Maddint(Mlshiftint($d,1),1);
    next unless Mis_prime($p);
    my $r = Mdivint($n,$d);
    while (1) {
      return 1 if $r == $p || _totpred($r, $d);
      my($Q,$R) = divrem($r,$p);
      last if $R != 0;
      $r = $Q;
    }
  }
  0;
}
sub is_totient {
  my($n) = @_;
  validate_integer($n);
  return 0+($n==1) if $n <= 1;
  return _totpred($n,$n);
}


sub moebius_range {
  my($lo, $hi) = @_;
  validate_integer($lo);
  validate_integer($hi);
  return () if $hi < $lo;
  return moebius($lo) if $lo == $hi;
  if ($lo < 0) {
    if ($hi < 0) {
      return reverse(moebius_range(-$hi,-$lo));
    } else {
      return (reverse(moebius_range(1,-$lo)), moebius_range(0,$hi));
    }
  }
  my @mu;
  if ($hi > 2**32) {
    ($lo,$hi) = (tobigint($lo),tobigint($hi)) if $hi > 2**49;
    push @mu, moebius($lo++)  while $lo <= $hi;
    return @mu;
  }
  for (my $i = $lo; $i <= $hi; $i++) { push @mu, 1; }
  $mu[0] = 0 if $lo == 0;
  my($p, $sqrtn) = (2, Msqrtint($hi));
  while ($p <= $sqrtn) {
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
    $p = Mnext_prime($p);
  }
  for (my $i = $lo; $i <= $hi; $i++) {
    my $m = $mu[$i-$lo];
    $m *= -1 if abs($m) != $i;
    $mu[$i-$lo] = ($m>0) - ($m<0);
  }
  return @mu;
}

sub _omertens {
  my($n) = @_;
  # This is the most basic Deléglise and Rivat algorithm.  u = n^1/2
  # and no segmenting is done.  Their algorithm uses u = n^1/3, breaks
  # the summation into two parts, and calculates those in segments.  Their
  # computation time growth is half of this code.
  return $n if $n <= 1;
  my $u = int(sqrt($n));
  my @mu = (0, Mmoebius(1, $u)); # Hold values of mu for 0-u
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

sub _rmertens {
  my($n, $Mref, $Href, $size) = @_;
  return $Mref->[$n] if $n <= $size;
  return $Href->{$n} if exists $Href->{$n};
  my $s = Msqrtint($n);
  my $ns = int($n/($s+1));

  my ($nk, $nk1) = ($n, Mrshiftint($n,1));
  my $SUM = Msubint(1,Msubint($nk,$nk1));
  my @S;
  foreach my $k (2 .. $ns) {
    ($nk, $nk1) = ($nk1, Mdivint($n,$k+1));
    push @S, ($nk <= $size) ? $Mref->[$nk]
                            : _rmertens($nk, $Mref, $Href, $size);
    push @S, $Mref->[$k] * ($nk - $nk1);
  }
  push @S, Mmulint($Mref->[$s], Mdivint($n,$s) - $ns) if $s > $ns;
  $SUM = Msubint($SUM, Mvecsum(@S));

  $Href->{$n} = $SUM;
  $SUM;
}

sub mertens {
  my($n) = @_;

  return _omertens($n) if $n < 20000;

  # Larger size would be faster, but more memory.
  my $size = (Mrootint($n, 3)**2) >> 2;
  $size = Msqrtint($n) if $size < Msqrtint($n);

  my @M = (0);
  push @M, $M[-1] + $_ for Mmoebius(1, $size);

  my %seen;
  return _rmertens($n, \@M, \%seen, $size);
}


sub ramanujan_sum {
  my($k,$n) = @_;
  return 0 if $k < 1 || $n <  1;
  my $g = $k / Mgcd($k,$n);
  my $m = Mmoebius($g);
  return $m if $m == 0 || $k == $g;
  $m * (Mtotient($k) / Mtotient($g));
}

sub liouville {
  my($n) = @_;
  my $l = (-1) ** scalar Mfactor($n);
  return $l < 0 ? -1 : $l > 0 ? 1 : 0;
}

sub sumliouville {
  my($n) = @_;
  return (0,1,0,-1,0,-1,0,-1,-2,-1,0,-1,-2,-3,-2,-1)[$n] if $n < 16;

  # Build the Mertens lookup info once.
  my $sqrtn = Msqrtint($n);
  my $size = (Mrootint($n, 3)**2) >> 2;
  $size = $sqrtn if $size < $sqrtn;
  my %seen;
  my @M = (0);
  push @M, $M[-1] + $_ for Mmoebius(1, $size);

  # L(n) = sum[k=1..sqrt(n)](Mertens(n/(k^2)))
  my $L = 0;
  for my $k (1 .. $sqrtn) {
    #my $nk = Mdivint($n, Mmulint($k,$k));
    my $nk = int($n/($k*$k));
    return Mvecsum($L, $sqrtn, -$k, 1) if $nk == 1;
    $L = Maddint($L,($nk <= $size) ? $M[$nk] : _rmertens($nk,\@M,\%seen,$size));
  }
  return $L;
}

# Exponential of Mangoldt function (A014963).
# Return p if n = p^m [p prime, m >= 1], 1 otherwise.
sub exp_mangoldt {
  my($n) = @_;
  my $p;
  return 1 unless Mis_prime_power($n,\$p);
  $p;
}

sub carmichael_lambda {
  my($n) = @_;
  return Mtotient($n) if $n < 8;           # = phi(n) for n < 8
  return $n >> 2 if ($n & ($n-1)) == 0;    # = phi(n)/2 = n/4 for 2^k, k>2

  my @pe = Mfactor_exp($n);
  $pe[0]->[1]-- if $pe[0]->[0] == 2 && $pe[0]->[1] > 2;

  my $lcm;
  if (!ref($n)) {
    $lcm = Mlcm(
      map { ($_->[0] ** ($_->[1]-1)) * ($_->[0]-1) } @pe
    );
  } else {
    $lcm = Mlcm(
      map { Mmulint(Mpowint($_->[0],$_->[1]-1),$_->[0]-1) } @pe
    );
  }
  $lcm;
}

sub is_cyclic {
  my($n) = @_;
  validate_integer($n);

  return 0+($n > 0) if $n < 4;
  return 0 if ($n % 2) == 0;
  return 0 if (!($n % 9) || !($n%25) || !($n%49) || !($n%121));
  return 0 if (!($n %21) || !($n%39) || !($n%55) || !($n%57) || !($n%93));

  return 1 if Mgcd($n,Mtotient($n)) == 1;
  0;
}

sub is_carmichael {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n < 561 || ($n % 2) == 0;

  return reftyped($_[0], Math::Prime::Util::GMP::is_carmichael($n))
    if $Math::Prime::Util::_GMPfunc{"is_carmichael"};

  # This works fine, but very slow
  # return !is_prime($n) && ($n % carmichael_lambda($n)) == 1;

  return 0 if (!($n % 9) || !($n % 25) || !($n%49) || !($n%121));

  # Check Korselt's criterion for small divisors
  my $fn = $n;
  for my $a (5,7,11,13,17,19,23,29,31,37,41,43) {
    if (($fn % $a) == 0) {
      return 0 if (($n-1) % ($a-1)) != 0;   # Korselt
      $fn /= $a;
      return 0 unless $fn % $a;             # not square free
    }
  }
  return 0 if Mpowmod(2, $n-1, $n) != 1;

  # After pre-tests, it's reasonably likely $n is a Carmichael number or prime

  # Use probabilistic test if too large to reasonably factor.
  if (length($fn) > 50) {
    return 0 if Mis_prime($n);
    for my $t (13 .. 150) {
      my $a = $_primes_small[$t];
      my $gcd = Mgcd($a, $fn);
      if ($gcd == 1) {
        return 0 if Mpowmod($a, $n-1, $n) != 1;
      } else {
        return 0 if $gcd != $a;              # Not square free
        return 0 if (($n-1) % ($a-1)) != 0;  # factor doesn't divide
        $fn /= $a;
      }
    }
    return 1;
  }

  # Verify with factoring.
  my @pe = Mfactor_exp($n);
  return 0 if scalar(@pe) < 3;
  for my $pe (@pe) {
    return 0 if $pe->[1] > 1 || (($n-1) % ($pe->[0]-1)) != 0;
  }
  1;
}

sub is_quasi_carmichael {
  my($n) = @_;
  validate_integer_nonneg($n);

  return 0 if $n < 35;
  return 0 if (!($n % 4) || !($n % 9) || !($n % 25) || !($n%49) || !($n%121));

  my @pe = Mfactor_exp($n);
  # Not quasi-Carmichael if prime
  return 0 if scalar(@pe) < 2;
  # Not quasi-Carmichael if not square free
  for my $pe (@pe) {
    return 0 if $pe->[1] > 1;
  }
  my @f = map { $_->[0] } @pe;
  my $nbases = 0;
  if ($n < 2000) {
    # In theory for performance, but mainly keeping to show direct method.
    my $lim = $f[-1];
    $lim = (($n-$lim*$lim) + $lim - 1) / $lim;
    for my $b (1 .. $f[0]-1) {
      my $nb = $n - $b;
      $nbases++ if Mvecall(sub { $nb % ($_-$b) == 0 }, @f);
    }
    if (scalar(@f) > 2) {
      for my $b (1 .. $lim-1) {
        my $nb = $n + $b;
        $nbases++ if Mvecall(sub { $nb % ($_+$b) == 0 }, @f);
      }
    }
  } else {
    my($spf,$lpf) = ($f[0], $f[-1]);
    if (scalar(@f) == 2) {
      foreach my $d (Mdivisors($n/$spf - 1)) {
        my $k = $spf - $d;
        my $p = $n - $k;
        last if $d >= $spf;
        $nbases++ if Mvecall(sub { my $j = $_-$k;  $j && ($p % $j) == 0 }, @f);
      }
    } else {
      foreach my $d (Mdivisors($lpf * ($n/$lpf - 1))) {
        my $k = $lpf - $d;
        my $p = $n - $k;
        next if $k == 0 || $k >= $spf;
        $nbases++ if Mvecall(sub { my $j = $_-$k;  $j && ($p % $j) == 0 }, @f);
      }
    }
  }
  $nbases;
}

sub is_pillai {
  my($p) = @_;
  validate_integer($p);
  return 0 if $p < 23;
  return 0 unless $p % 2 && $p % 3 && $p % 5 && $p % 7;

  my $pm1 = $p-1;
  my $nfac = 5040 % $p;
  for (my $n = 8; $n < $p; $n++) {
    $nfac = Mmulmod($nfac, $n, $p);
    return $n if $nfac == $pm1 && ($p % $n) != 1;
  }
  0;
}

sub is_fundamental {
  my($n) = @_;
  validate_integer($n);
  my $neg = ($n < 0);
  $n = -$n if $neg;
  my $r = $n & 15;
  if ($r) {
    my $r4 = $r & 3;
    if (!$neg) {
      return (($r ==  4) ? 0 : Mis_square_free($n >> 2)) if $r4 == 0;
      return Mis_square_free($n) if $r4 == 1;
    } else {
      return (($r == 12) ? 0 : Mis_square_free($n >> 2)) if $r4 == 0;
      return Mis_square_free($n) if $r4 == 3;
    }
  }
  0;
}

my @_ds_overflow =  # We'll use BigInt math if the input is larger than this.
  (~0 > 4294967295)
   ? (124, 3000000000000000000, 3000000000, 2487240, 64260, 7026)
   : ( 50,           845404560,      52560,    1548,   252,   84);
sub divisor_sum {
  my($n, $k) = @_;
  validate_integer_nonneg($n);
  return 0 if $n == 0;

  if (defined $k && ref($k) eq 'CODE') {
    my $sum = $n-$n;
    my $refn = ref($n);
    foreach my $d (Mdivisors($n)) {
      $sum += $k->( $refn ? $refn->new("$d") : $d );
    }
    return $sum;
  }
  return 1 if $n == 1;

  croak "Second argument must be a code ref or number"
    unless !defined $k || validate_integer_nonneg($k);
  $k = 1 if !defined $k;

  return reftyped($_[0], Math::Prime::Util::GMP::sigma($n, $k))
    if $Math::Prime::Util::_GMPfunc{"sigma"};

  my @factors = Mfactor_exp($n);

  return Mvecprod(map { $_->[1]+1 } @factors) if $k == 0;

  my @prod;

  if ($k == 1) {
    foreach my $f (@factors) {
      my ($p, $e) = @$f;
      if ($e == 1) {
        push @prod, $p+1;
      } elsif ($e == 2 && $p < 65536) {
        push @prod, ($p+1) + $p * $p;
      } else {
        push @prod, Mvecsum($p+1, map { Mpowint($p,$_) } 2..$e);
      }
    }
  } else {
    foreach my $f (@factors) {
      my ($p, $e) = @$f;
      my $pk = Mpowint($p,$k);
      if ($e == 1) {
        push @prod, Maddint($pk,1);
      } else {
        push @prod, Mvecsum(Maddint($pk,1), map { Mpowint($pk,$_) } 2..$e);
      }
    }
  }
  return $prod[0] if @prod == 1;
  if ($k == 1 && $n < 845404560) {   # divisor_sum(n) < 2^32
    my $r = 1;
    $r *= $_ for @prod;
    return $r;
  }
  return Mmulint($prod[0],$prod[1]) if @prod == 2;
  Mvecprod(@prod);
}

#############################################################################
#                       Lehmer prime count
#
#my @_s0 = (0);
#my @_s1 = (0,1);
#my @_s2 = (0,1,1,1,1,2);
#my @_s3 = (0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7,8);
#my @_s4 = (0,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7,7,8,8,8,8,8,8,9,9,9,9,10,10,11,11,11,11,12,12,12,12,12,12,13,13,13,13,13,13,14,14,15,15,15,15,15,15,16,16,16,16,17,17,18,18,18,18,18,18,19,19,19,19,20,20,20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,22,23,23,24,24,24,24,25,25,26,26,26,26,27,27,27,27,27,27,27,27,28,28,28,28,28,28,29,29,29,29,30,30,30,30,30,30,31,31,32,32,32,32,33,33,33,33,33,33,34,34,35,35,35,35,35,35,36,36,36,36,36,36,37,37,37,37,38,38,39,39,39,39,40,40,40,40,40,40,41,41,42,42,42,42,42,42,43,43,43,43,44,44,45,45,45,45,46,46,47,47,47,47,47,47,47,47,47,47,48);
my(@_s3,@_s4);
my @_pred5 = (1,0,1,2,3,4,5,0,1,2,3,0,1,0,1,2,3,0,1,0,1,2,3,0,1,2,3,4,5,0);

sub _tablephi {
  my($x, $a) = @_;
  if ($a == 0) { return $x; }
  elsif ($a == 1) { return $x-int($x/2); }
  elsif ($a == 2) { return $x-int($x/2) - int($x/3) + int($x/6); }
  elsif ($a == 3) { return  8 * int($x /  30) + $_s3[$x %  30]; }
  elsif ($a == 4) { return 48 * int($x / 210) + $_s4[$x % 210]; }
  elsif ($a == 5) { my $xp = int($x/11);
                    return ( (48 * int($x   / 210) + $_s4[$x   % 210]) -
                             (48 * int($xp  / 210) + $_s4[$xp  % 210]) ); }
  else            { my ($xp,$x2) = (int($x/11),int($x/13));
                    my $x2p = int($x2/11);
                    return ( (48 * int($x   / 210) + $_s4[$x   % 210]) -
                             (48 * int($xp  / 210) + $_s4[$xp  % 210]) -
                             (48 * int($x2  / 210) + $_s4[$x2  % 210]) +
                             (48 * int($x2p / 210) + $_s4[$x2p % 210]) ); }
}

sub legendre_phi {
  my ($x, $a, $primes) = @_;
  validate_integer_nonneg($x);
  validate_integer_nonneg($a);
  if ($#_s3 == -1) {
    @_s3 = (0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7,8);
    @_s4 = (0,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7,7,8,8,8,8,8,8,9,9,9,9,10,10,11,11,11,11,12,12,12,12,12,12,13,13,13,13,13,13,14,14,15,15,15,15,15,15,16,16,16,16,17,17,18,18,18,18,18,18,19,19,19,19,20,20,20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,22,23,23,24,24,24,24,25,25,26,26,26,26,27,27,27,27,27,27,27,27,28,28,28,28,28,28,29,29,29,29,30,30,30,30,30,30,31,31,32,32,32,32,33,33,33,33,33,33,34,34,35,35,35,35,35,35,36,36,36,36,36,36,37,37,37,37,38,38,39,39,39,39,40,40,40,40,40,40,41,41,42,42,42,42,42,42,43,43,43,43,44,44,45,45,45,45,46,46,47,47,47,47,47,47,47,47,47,47,48);
  }
  return _tablephi($x,$a) if $a <= 6;
  $primes = Mprimes(Mnth_prime_upper($a+1)) unless defined $primes;
  return ($x > 0 ? 1 : 0) if $x < $primes->[$a];

  my $sum = 0;
  my %vals = ( $x => 1 );
  while ($a > 6) {
    my $primea = $primes->[$a-1];
    my %newvals;
    while (my($v,$c) = each %vals) {
      my $sval = int($v / $primea);
      $sval -= $_pred5[$sval % 30];   # Reduce sval to one with same phi.
      if ($sval < $primea) {
        $sum -= $c;
      } else {
        $newvals{$sval} -= $c;
      }
    }
    # merge newvals into vals
    while (my($v,$c) = each %newvals) {
      $vals{$v} += $c;
      delete $vals{$v} if $vals{$v} == 0;
    }
    $a--;
  }
  while (my($v,$c) = each %vals) {
    $sum += $c * _tablephi($v, $a);
  }
  return $sum;
}

sub _sieve_prime_count {
  my $high = shift;
  return (0,0,1,2,2,3,3)[$high] if $high < 7;
  $high-- unless ($high % 2);
  return 1 + ${_sieve_erat($high)} =~ tr/0//;
}

sub _count_with_sieve {
  my ($sref, $low, $high) = @_;
  ($low, $high) = (2, $low) if !defined $high;
  my $count = 0;
  if   ($low < 3) { $low = 3; $count++; }
  else            { $low |= 1; }
  $high-- unless ($high % 2);
  return $count if $low > $high;
  my $sbeg = $low >> 1;
  my $send = $high >> 1;

  if ( !defined $sref || $send >= length($$sref) ) {
    # outside our range, so call the segment siever.
    my $seg_ref = _sieve_segment($low, $high);
    return $count + $$seg_ref =~ tr/0//;
  }
  return $count + substr($$sref, $sbeg, $send-$sbeg+1) =~ tr/0//;
}

sub _lehmer_pi {
  my $x = shift;
  return _sieve_prime_count($x) if $x < 1_000;

  my $z = Msqrtint($x);
  my $a = _lehmer_pi(Msqrtint($z));
  my $b = _lehmer_pi($z);
  my $c = _lehmer_pi(Mrootint($x,3));

  # Generate at least b primes.
  my $bth_prime_upper = ($b <= 10) ? 29 : int("$b"*(log("$b")+log(log("$b")))) + 1;
  my $primes = Mprimes( $bth_prime_upper );

  my $sum = Mmulint(Mvecsum($b,$a,-2),Mvecsum($b,-$a,1)) >> 1;
  $sum += legendre_phi($x, $a, $primes);

  # Get a big sieve for our primecounts.  The C code compromises with either
  # b*10 or x^3/5, as that cuts out all the inner loop sieves and about half
  # of the big outer loop counts.
  # Our sieve count isn't nearly as optimized here, so error on the side of
  # more primes.  This uses a lot more memory but saves a lot of time.
  my $sref = _sieve_erat( Mdivint(Mdivint($x,$primes->[$a]),5) );

  my ($lastw, $lastwpc) = (0,0);
  foreach my $i (reverse $a+1 .. $b) {
    my $w = int($x / $primes->[$i-1]);
    $lastwpc += _count_with_sieve($sref,$lastw+1, $w);
    $lastw = $w;
    $sum -= $lastwpc;
    #$sum -= _count_with_sieve($sref,$w);
    if ($i <= $c) {
      my $bi = _count_with_sieve($sref,int(sqrt($w)+0.5));
      foreach my $j ($i .. $bi) {
        $sum = $sum - _count_with_sieve($sref,int($w / $primes->[$j-1])) + $j - 1;
      }
    }
  }
  $sum;
}
#############################################################################


sub prime_count {
  my($low,$high) = @_;
  if (defined $high) { validate_integer_nonneg($low); }
  else               { ($low,$high) = (2, $low);      }
  validate_integer_nonneg($high);

  return 0 if $high < 2 || $low > $high;

  return reftyped($high, Math::Prime::Util::GMP::prime_count($low,$high))
    if $Math::Prime::Util::_GMPfunc{"prime_count"}
    && (ref($high) eq 'Math::BigInt' || ($high-$low) < int($low/1_000_000));

  my $count = 0;
  $count++ if ($low <= 2) && ($high >= 2);   # Count 2
  $low = 3 if $low < 3;
  $low++ if ($low % 2) == 0;   # Make low go to odd number.
  $high-- if ($high % 2) == 0; # Make high go to odd number.
  return $count if $low > $high;

  if (   ref($low) eq 'Math::BigInt' || ref($high) eq 'Math::BigInt'
      || ($high-$low) < 10
      || ($high-$low) < int($low/100_000_000_000) ) {
    # Trial primes seems best.  Needs some tuning.
    my $curprime = Mnext_prime($low-1);
    while ($curprime <= $high) {
      $count++;
      $curprime = Mnext_prime($curprime);
    }
    return $count;
  }

  # TODO: Needs tuning
  if ($high > 50_000) {
    if ( ($high / ($high-$low+1)) < 100 ) {
      $count += _lehmer_pi($high);
      $count -= ($low == 3) ? 1 : _lehmer_pi($low-1);
      return $count;
    }
  }

  return (_sieve_prime_count($high) - 1 + $count) if $low == 3;

  my $sieveref = _sieve_segment($low,$high);
  $count += $$sieveref =~ tr/0//;
  return $count;
}


sub nth_prime {
  my($n) = @_;
  validate_integer_nonneg($n);

  return undef if $n <= 0;  ## no critic qw(ProhibitExplicitReturnUndef)
  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if ref($n) || $n > MPU_MAXPRIMEIDX || $n > 2**45;

  my $prime = 0;
  my $count = 1;
  my $start = 3;

  my $logn = log($n);
  my $loglogn = log($logn);
  my $nth_prime_upper = ($n <= 10) ? 29 : int($n*($logn + $loglogn)) + 1;
  if ($nth_prime_upper > 100000) {
    # Use fast Lehmer prime count combined with lower bound to get close.
    my $nth_prime_lower = int($n * ($logn + $loglogn - 1.0 + (($loglogn-2.10)/$logn)));
    $nth_prime_lower-- unless $nth_prime_lower % 2;
    $count = _lehmer_pi($nth_prime_lower);
    $start = $nth_prime_lower + 2;
  }

  {
    # Make sure incr is an even number.
    my $incr = ($n < 1000) ? 1000 : ($n < 10000) ? 10000 : 100000;
    my $sieveref;
    while (1) {
      $sieveref = _sieve_segment($start, $start+$incr);
      my $segcount = $$sieveref =~ tr/0//;
      last if ($count + $segcount) >= $n;
      $count += $segcount;
      $start += $incr+2;
    }
    # Our count is somewhere in this segment.  Need to look for it.
    $prime = $start - 2;
    while ($count < $n) {
      $prime += 2;
      $count++ if !substr($$sieveref, ($prime-$start)>>1, 1);
    }
  }
  $prime;
}

# The nth prime will be less or equal to this number
sub nth_prime_upper {
  my($n) = @_;
  validate_integer_nonneg($n);

  return undef if $n <= 0;  ## no critic qw(ProhibitExplicitReturnUndef)
  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if ref($n) || $n > MPU_MAXPRIMEIDX || $n > 2**45;

  my $flogn  = log($n);
  my $flog2n = log($flogn);  # Note distinction between log_2(n) and log^2(n)

  my $upper;
  if      ($n >= 46254381) {  # Axler 2017 Corollary 1.2
    $upper = $n * ( $flogn  +  $flog2n-1.0  +  (($flog2n-2.00)/$flogn)  -  (($flog2n*$flog2n - 6*$flog2n + 10.667)/(2*$flogn*$flogn)) );
  } elsif ($n >=  8009824) {  # Axler 2013 page viii Korollar G
    $upper = $n * ( $flogn  +  $flog2n-1.0  +  (($flog2n-2.00)/$flogn)  -  (($flog2n*$flog2n - 6*$flog2n + 10.273)/(2*$flogn*$flogn)) );
  } elsif ($n >=  688383) {   # Dusart 2010 page 2
    $upper = $n * ( $flogn  +  $flog2n - 1.0 + (($flog2n-2.00)/$flogn) );
  } elsif ($n >=  178974) {   # Dusart 2010 page 7
    $upper = $n * ( $flogn  +  $flog2n - 1.0 + (($flog2n-1.95)/$flogn) );
  } elsif ($n >=   39017) {   # Dusart 1999 page 14
    $upper = $n * ( $flogn  +  $flog2n - 0.9484 );
  } elsif ($n >=       6) {   # Modified Robin 1983, for 6-39016 only
    $upper = $n * ( $flogn  +  0.6000 * $flog2n );
  } else {
    $upper = $n * ( $flogn  +  $flog2n );
  }

  Mtoint($upper + 1.0);
}

# The nth prime will be greater than or equal to this number
sub nth_prime_lower {
  my($n) = @_;
  validate_integer_nonneg($n);

  return undef if $n <= 0;  ## no critic qw(ProhibitExplicitReturnUndef)
  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if ref($n) || $n > MPU_MAXPRIMEIDX || $n > 2**45;

  my $flogn  = log($n);
  my $flog2n = log($flogn);  # Note distinction between log_2(n) and log^2(n)

  # Dusart 1999 page 14, for all n >= 2
  #my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.25)/$flogn));
  # Dusart 2010 page 2, for all n >= 3
  #my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.10)/$flogn));
  # Axler 2013 page viii Korollar I, for all n >= 2
  #my $lower = $n * ($flogn + $flog2n-1.0 + (($flog2n-2.00)/$flogn) - (($flog2n*$flog2n - 6*$flog2n + 11.847)/(2*$flogn*$flogn)) );
  # Axler 2017 Corollary 1.4
  my $lower = $n * ($flogn + $flog2n-1.0 + (($flog2n-2.00)/$flogn) - (($flog2n*$flog2n - 6*$flog2n + 11.508)/(2*$flogn*$flogn)) );

  my $plower = Mtoint($lower + 0.999999999);
  # We clamp to the max UV representable.
  if (MPU_32BIT) {
    $plower = 4294967291 if $n >= 203280221 && $plower < 4294967291;
  } else {
    $plower = 18446744073709551557 if $n >= 425656284035217743 && $plower < 18446744073709551557;
  }
  $plower;
}

sub inverse_li_nv {
  my($n) = @_;
  $n = 0.0 + "$n";
  my $t = $n * log($n);

  # Iterate Halley's method until error term grows
  my $old_term = MPU_INFINITY;
  for my $iter (1 .. 10000) {
    my $dn = MLi($t) - $n;
    $dn = 0.0 + "$dn" if ref($dn);
    my $term = $dn * log($t) / (1.0 + $dn/(2*$t));
    last if abs($term) >= abs($old_term);
    $old_term = $term;
    $t -= $term;
    last if abs($term) < 1e-6;
  }
  $t;
}

sub inverse_li {
  my($n) = @_;
  validate_integer_nonneg($n);

  return (0,2,3,5,6,8)[$n] if $n <= 5;
  my $t = Math::Prime::Util::inverse_li_nv(0.0 + "$n");

  $t = Mtoint($t + 0.5);

  # Make it an exact answer
  my $inc = ($n > 4e16) ? 2048 : 128;
  if (int(MLi($t-1)) >= $n) {
    $t -= $inc while int(MLi($t-$inc)) >= $n;
    for ($inc = $inc >> 1;  $inc > 0;  $inc >>= 1) {
      $t -= $inc if int(MLi($t-$inc)) >= $n;
    }
  } else {
    $t += $inc while int(MLi($t+$inc-1)) < $n;
    for ($inc = $inc >> 1;  $inc > 0;  $inc >>= 1) {
      $t += $inc if int(MLi($t+$inc-1)) < $n;
    }
  }

  $t;
}
# Not currently used
sub _inverse_R {
  my($n) = @_;
  validate_integer_nonneg($n);

  return (0,2,3,5,6,8)[$n] if $n <= 5;
  $n = _upgrade_to_float($n) if ref($n) || $n > MPU_MAXPRIMEIDX || $n > 2**45;
  my $t = $n * log($n);

  # Iterate Halley's method until error term grows
  my $old_term = MPU_INFINITY;
  for my $iter (1 .. 10000) {
    my $dn = Math::Prime::Util::RiemannR($t) - $n;
    my $term = $dn * log($t) / (1.0 + $dn/(2*$t));
    last if abs($term) >= abs($old_term);
    $old_term = $term;
    $t -= $term;
    last if abs($term) < 1e-6;
  }
  Mtoint( ref($t) ? $t->bceil->bstr : $t+0.99999 );
}

sub nth_prime_approx {
  my($n) = @_;
  validate_integer_nonneg($n);

  return undef if $n <= 0;  ## no critic qw(ProhibitExplicitReturnUndef)
  return $_primes_small[$n] if $n <= $#_primes_small;

  # Once past 10^12 or so, inverse_li gives better results.
  return Math::Prime::Util::inverse_li($n) if $n > 1e12;

  $n = _upgrade_to_float($n) if ref($n) || $n >= MPU_MAXPRIMEIDX;

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

  my $approx = $n * ( $flogn + $flog2n - 1
                      + (($flog2n - 2)/$flogn)
                      - ((($flog2n*$flog2n) - 6*$flog2n + 11) / (2*$flogn*$flogn))
                    );

  # Apply a correction to help keep values close.
  my $order = $flog2n/$flogn;
  $order = $order*$order*$order * $n;

  if    ($n <        259) { $approx += 10.4 * $order; }
  elsif ($n <        775) { $approx +=  6.3 * $order; }
  elsif ($n <       1271) { $approx +=  5.3 * $order; }
  elsif ($n <       2000) { $approx +=  4.7 * $order; }
  elsif ($n <       4000) { $approx +=  3.9 * $order; }
  elsif ($n <      12000) { $approx +=  2.8 * $order; }
  elsif ($n <     150000) { $approx +=  1.2 * $order; }
  elsif ($n <   20000000) { $approx +=  0.11 * $order; }
  elsif ($n <  100000000) { $approx +=  0.008 * $order; }
  elsif ($n <  500000000) { $approx += -0.038 * $order; }
  elsif ($n < 2000000000) { $approx += -0.054 * $order; }
  else                    { $approx += -0.058 * $order; }
  # If we want the asymptotic approximation to be >= actual, use -0.010.

  Mtoint($approx + 0.5);
}

#############################################################################

sub prime_count_approx {
  my($x) = @_;
  validate_integer_nonneg($x);
  #return (0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6)[$x] if $x < 16;
  return _tiny_prime_count($x) if $x < $_primes_small[-1];

  # Turn on high precision FP if they gave us a big number.
  $x = _upgrade_to_float($x) if ref($x) || $x >= 1e16;

  #    Method             10^10 %error  10^19 %error
  #    -----------------  ------------  ------------
  #    n/(log(n)-1)        .22%          .058%
  #    n/(ln(n)-1-1/ln(n)) .032%         .0041%
  #    average bounds      .0005%        .0000002%
  #    asymp               .0006%        .00000004%
  #    li(n)               .0007%        .00000004%
  #    li(n)-li(n^.5)/2    .0004%        .00000001%
  #    R(n)                .0004%        .00000001%
  #
  # Also consider: http://trac.sagemath.org/sage_trac/ticket/8135

  # Asymp:
  #   my $l1 = log($x);  my $l2 = $l1*$l1;  my $l4 = $l2*$l2;
  #   my $result = int( $x/$l1 + $x/$l2 + 2*$x/($l2*$l1) + 6*$x/($l4) + 24*$x/($l4*$l1) + 120*$x/($l4*$l2) + 720*$x/($l4*$l2*$l1) + 5040*$x/($l4*$l4) + 40320*$x/($l4*$l4*$l1) + 0.5 );
  # my $result = int( (prime_count_upper($x) + prime_count_lower($x)) / 2);
  # my $result = int( LogarithmicIntegral($x) );
  # my $result = int(LogarithmicIntegral($x) - LogarithmicIntegral(sqrt($x))/2);
  # my $result = RiemannR($x) + 0.5;

  # Make sure we get enough accuracy, and also not too much more than needed
  $x->accuracy(length($x->copy->as_int->bstr())+2) if ref($x) =~ /^Math::Big/;

  my $result;
  if ($Math::Prime::Util::_GMPfunc{"riemannr"} || !ref($x)) {
    # Fast if we have our GMP backend, and ok for native.
    $result = Math::Prime::Util::PP::RiemannR($x);
  } else {
    $x = _upgrade_to_float($x) unless ref($x) eq 'Math::BigFloat';
    $result = Math::BigFloat->new(0);
    $result->accuracy($x->accuracy) if ref($x) && $x->accuracy;
    $result += Math::BigFloat->new(MLi($x));
    $result -= Math::BigFloat->new(MLi(sqrt($x))/2);
    my $intx = ref($x) ? Math::BigInt->new($x->bfround(0)) : $x;
    for my $k (3 .. 1000) {
      my $m = moebius($k);
      next unless $m != 0;
      # With Math::BigFloat and the Calc backend, FP root is ungodly slow.
      # Use integer root instead.  For more accuracy (not useful here):
      # my $v = Math::BigFloat->new( "" . Mrootint($x->as_int,$k) );
      # $v->accuracy(length($v)+5);
      # $v = $v - Math::BigFloat->new(($v**$k - $x))->bdiv($k * $v**($k-1));
      # my $term = LogarithmicIntegral($v)/$k;
      my $term = MLi(Mrootint($intx,$k)) / $k;
      last if $term < .25;
      if ($m == 1) { $result->badd(Math::BigFloat->new($term)) }
      else         { $result->bsub(Math::BigFloat->new($term)) }
    }
  }

  Mtoint($result+0.5);
}

sub prime_count_lower {
  my($x) = @_;
  validate_integer_nonneg($x);

  return _tiny_prime_count($x) if $x < $_primes_small[-1];

  return reftyped($_[0], Math::Prime::Util::GMP::prime_count_lower($x))
    if $Math::Prime::Util::_GMPfunc{"prime_count_lower"};

  $x = _upgrade_to_float($x) if ref($x);

  my($result,$a);
  my $fl1 = log($x);
  my $fl2 = $fl1*$fl1;
  my $one = (ref($x) eq 'Math::BigFloat') ? $x->copy->bone : $x-$x+1.0;

  # Chebyshev            1*x/logx       x >= 17
  # Rosser & Schoenfeld  x/(logx-1/2)   x >= 67
  # Dusart 1999          x/logx*(1+1/logx+1.8/logxlogx)  x >= 32299
  # Dusart 2010          x/logx*(1+1/logx+2.0/logxlogx)  x >= 88783
  # Axler 2014 (1.2)     ""+...                          x >= 1332450001
  # Axler 2014 (1.2)     x/(logx-1-1/logx-...)           x >= 1332479531
  # Büthe 2015 (1.9)     li(x)-(sqrtx/logx)*(...)        x <= 10^19
  # Büthe 2014 Th 2      li(x)-logx*sqrtx/8Pi    x > 2657, x <= 1.4   * 10^25
  # Johnston 2021 Cor3.3 li(x)-logx*sqrtx/8Pi    x > 2657, x <= 1.101 * 10^26

  # Also see Dusart 2018: if RH and x >= 5639,
  #     |pi(x)-li(x)|<= x * (logx-loglogx)/(8*Pi*sqrtx)
  # TODO: evaluate this

  if ($x < 599) {                         # Decent for small numbers
    $result = $x / ($fl1 - 0.7);
  } elsif ($x < 52600000) {               # Dusart 2010 tweaked
    if    ($x <       2700) { $a = 0.30; }
    elsif ($x <       5500) { $a = 0.90; }
    elsif ($x <      19400) { $a = 1.30; }
    elsif ($x <      32299) { $a = 1.60; }
    elsif ($x <      88783) { $a = 1.83; }
    elsif ($x <     176000) { $a = 1.99; }
    elsif ($x <     315000) { $a = 2.11; }
    elsif ($x <    1100000) { $a = 2.19; }
    elsif ($x <    4500000) { $a = 2.31; }
    else                    { $a = 2.35; }
    $result = ($x/$fl1) * ($one + $one/$fl1 + $a/$fl2);
  } elsif ($x < 1.1e26 || getconfig()->{'assume_rh'}){
                                          # Büthe 2014/2015
    my $lix = MLi($x);
    my $sqx = sqrt($x);
    if ($x < 1e19) {
      $result = $lix - ($sqx/$fl1) * (1.94 + 3.88/$fl1 + 27.57/$fl2);
    } else {
      if (ref($x) eq 'Math::BigFloat') {
        my $xdigits = _find_big_acc($x);
        $result = $lix - ($fl1*$sqx / (Math::BigFloat->bpi($xdigits)*8));
      } else {
        $result = $lix - ($fl1*$sqx / PI_TIMES_8);
      }
    }
  } else {                                # Axler 2014 1.4
    my($fl3,$fl4) = ($fl2*$fl1,$fl2*$fl2);
    my($fl5,$fl6) = ($fl4*$fl1,$fl4*$fl2);
    $result = $x / ($fl1 - $one - $one/$fl1 - 2.65/$fl2 - 13.35/$fl3 - 70.3/$fl4 - 455.6275/$fl5 - 3404.4225/$fl6);
  }
  # This will truncate bigfloat or floats to native int or bigint class.
  Mtoint($result);
}

sub prime_count_upper {
  my($x) = @_;
  validate_integer_nonneg($x);

  # Give an exact answer for what we have in our little table.
  return _tiny_prime_count($x) if $x < $_primes_small[-1];

  return reftyped($_[0], Math::Prime::Util::GMP::prime_count_upper($x))
    if $Math::Prime::Util::_GMPfunc{"prime_count_upper"};

  $x = _upgrade_to_float($x) if ref($x);

  # Chebyshev:            1.25506*x/logx       x >= 17
  # Rosser & Schoenfeld:  x/(logx-3/2)         x >= 67
  # Panaitopol 1999:      x/(logx-1.112)       x >= 4
  # Dusart 1999:          x/logx*(1+1/logx+2.51/logxlogx)   x >= 355991
  # Dusart 2010:          x/logx*(1+1/logx+2.334/logxlogx)  x >= 2_953_652_287
  # Dusart 2018:          x/lx*(1+1/lx+2/lxlx+7.59/lxlxlx)  x > 1
  # Axler 2014:           x/(logx-1-1/logx-3.35/logxlogx...) x >= e^3.804
  # Büthe 2014 7.4        Schoenfeld bounds hold to x <= 1.4e25
  # Axler 2017 Prop 2.2   Schoenfeld bounds hold to x <= 5.5e25
  # Johnston 2021 Cor 3.3 Schoenfeld bounds hold to x <= 1.0e26
  # Skewes                li(x)                x < 1e14

  # TODO: Also look at these from Dusart (2018) [paywalled].
  # 1  If RH and x >= 5639, |pi(x)-li(x)|<= x * (logx-loglogx)/(8*Pi*sqrtx)
  # 2  pi(x) <= li(x) for all 2 <= x <= 10^20
  # 3  [li(x) - 2sqrt(x)/log(x)] <= pi(x) for 1090877 <= x <= 10^20
  #
  # See https://arxiv.org/pdf/2404.17165 page 9 for Mossinghoff and Trudgian.
  # Page 26 also points out the Dusart 2018 improvement to Schoenfeld.
  #   https://math.colgate.edu/~integers/y34/y34.pdf
  # Axler 2022:
  #   https://arxiv.org/pdf/2203.05917

  my($result,$a);
  my $fl1 = log($x);
  my $fl2 = $fl1 * $fl1;
  my $one = (ref($x) eq 'Math::BigFloat') ? $x->copy->bone : $x-$x+1.0;

  if ($x < 15900) {              # Tweaked Rosser-type
    $a = ($x < 1621) ? 1.048 : ($x < 5000) ? 1.071 : 1.098;
    $result = ($x / ($fl1 - $a)) + 1.0;
  } elsif ($x < 821800000) {     # Tweaked Dusart 2010
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
    $result = ($x/$fl1) * ($one + $one/$fl1 + $a/$fl2) + $one;
  } elsif ($x < 1e19) {                     # Skewes number lower limit
    $a = ($x < 110e7) ? 0.032 : ($x < 1001e7) ? 0.027 : ($x < 10126e7) ? 0.021 : 0.0;
    $result = MLi($x) - $a * $fl1*sqrt($x)/PI_TIMES_8;
  } elsif ($x < 1.1e26 || getconfig()->{'assume_rh'}) {
                                            # Schoenfeld / Büthe 2014 Th 7.4
    my $lix = MLi($x);
    my $sqx = sqrt($x);
    if (ref($x) eq 'Math::BigFloat') {
      my $xdigits = _find_big_acc($x);
      $result = $lix + ($fl1*$sqx / (Math::BigFloat->bpi($xdigits)*8));
    } else {
      $result = $lix + ($fl1*$sqx / PI_TIMES_8);
    }
  } else {                                  # Axler 2014 1.3
    my($fl3,$fl4) = ($fl2*$fl1,$fl2*$fl2);
    my($fl5,$fl6) = ($fl4*$fl1,$fl4*$fl2);
    $result = $x / ($fl1 - $one - $one/$fl1 - 3.35/$fl2 - 12.65/$fl3 - 71.7/$fl4 - 466.1275/$fl5 - 3489.8225/$fl6);
  }
  # This will truncate bigfloat or floats to native int or bigint class.
  Mtoint($result);
}

sub twin_prime_count {
  my($low,$high) = @_;
  if (defined $high) { validate_integer_nonneg($low); }
  else               { ($low,$high) = (2, $low);      }
  validate_integer_nonneg($high);
  my $sum = 0;
  while ($low <= $high) {
    my $seghigh = ($high-$high) + $low + 1e7 - 1;
    $seghigh = $high if $seghigh > $high;
    $sum += scalar(@{Math::Prime::Util::twin_primes($low,$seghigh)});
    $low = $seghigh + 1;
  }
  $sum;
}
sub _semiprime_count {
  my $n = shift;
  my($sum,$pc) = (0,0);
  Mforprimes( sub {
    $sum += Mprime_count(int($n/$_))-$pc++;
  }, Msqrtint($n));
  $sum;
}
sub semiprime_count {
  my($lo,$hi) = @_;
  if (defined $hi) { validate_integer_nonneg($lo); }
  else             { ($lo,$hi) = (2, $lo);         }
  validate_integer_nonneg($hi);
  # todo: threshold of fast count vs. walk
  if (($hi-$lo+1) < $hi / (sqrt($hi)/4)) {
    my $sum = 0;
    while ($lo <= $hi) {
      $sum++ if Mis_semiprime($lo);
      $lo++;
    }
    return $sum;
  }
  my $sum = _semiprime_count($hi) - (($lo < 4) ? 0 : semiprime_count($lo-1));
  $sum;
}

sub _kap_reduce_count {   # returns new k and n
  my($k, $n) = @_;

  my $pow3k = Mpowint(3, $k);
  while ($n < $pow3k) {
    $n = Mdivint($n, 2);
    $k--;
    $pow3k = Mdivint($pow3k, 3);
  }
  ($k, $n);
}
sub _kapc_final {              # k = 2
  my($n, $pdiv, $lo) = @_;
  my($sum, $hi, $pc) = (0,  Msqrtint(Mdivint($n,$pdiv)),  Mprime_count($lo)-1);
  my $nlim = int(INTMAX / $pdiv);
  Mforprimes( sub {
    my $npp = ($_<=$nlim) ? int($n/($pdiv*$_)) : Mdivint($n,Mmulint($pdiv,$_));
    $sum += Mprime_count($npp)-$pc++;
  }, $lo, $hi);
  $sum;
}
sub _kapc_count {
  my($n, $pdiv, $lo, $k) = @_;
  return _kapc_final($n, $pdiv, $lo)  if $k == 2;
  my($sum, $hi) = (0,  Mrootint(Mdivint($n,$pdiv),$k));
  Mforprimes(
    ($k == 3) ? sub { $sum += _kapc_final($n, Mmulint($pdiv,$_), $_); }
              : sub { $sum += _kapc_count($n, Mmulint($pdiv,$_), $_, $k-1); },
    $lo, $hi
  );
  $sum;
}
sub almost_prime_count {
  my($k,$n) = @_;
  validate_integer_nonneg($k);
  validate_integer_nonneg($n);
  return ($n >= 1) if $k == 0;
  my $ok = $k;
  ($k, $n) = _kap_reduce_count($k, $n);
  return $n if $k == 0;
  # If we reduced parameters, try again if XS might be able to do it.
  return Math::Prime::Util::almost_prime_count($k,$n) if $ok != $k && !ref($n) && getconfig()->{'xs'};
  return Mprime_count($n) if $k == 1;
  return Math::Prime::Util::semiprime_count($n) if $k == 2;
  return 0 if ($n >> $k) == 0;

  _kapc_count($n, 1, 2, $k);
}

sub _omega_prime_count_rec {
  my($k, $n,  $m, $p, $s, $j) = @_;
  $s = Mrootint(Mdivint($n,$m),$k) unless defined $s;
  $j = 1 unless defined $j;
  my $count = 0;

  if ($k == 2) {

    for (;  $p <= $s  ;  ++$j) {
      my $r = Mnext_prime($p);
      for (my $t = Mmulint($m, $p) ; $t <= $n ; $t = Mmulint($t, $p)) {
        my $w = Mdivint($n, $t);
        last if $r > $w;
        $count += Mprime_count($w) - $j;
        for (my $r2 = $r ; $r2 <= $w ; $r2 = Mnext_prime($r2)) {
          my $u = Mvecprod($t, $r2, $r2);
          last if $u > $n;
          for (; $u <= $n ; $u = Mmulint($u, $r2)) {
             ++$count;
          }
        }
      }
      $p = $r;
    }

  }  else  {

    for (;  $p <= $s  ;  ++$j) {
      my $r = Mnext_prime($p);
      for (my $t = Mmulint($m, $p) ; $t <= $n ; $t = Mmulint($t, $p)) {
        my $s = Mrootint(Mdivint($n, $t), $k - 1);
        last if $r > $s;
        $count += _omega_prime_count_rec($k-1, $n,  $t, $r, $s, $j+1);
      }
      $p = $r;
    }

  }
  $count;
}
sub omega_prime_count {
  my($k,$n) = @_;
  validate_integer_nonneg($k);
  validate_integer_nonneg($n);

  return ($n >= 1) ? 1 : 0 if $k == 0;
  return prime_power_count($n) if $k == 1;
  # find a simple formula for k=2.

  # Naive method
  # my ($sum, $low) = (0, Mpn_primorial($k));
  # for (my $i = $low; $i <= $n; $i++) {
  #   $sum++ if Mprime_omega($i) == $k;
  # }
  # return $sum;

  # Recursive method from trizen
  return _omega_prime_count_rec($k, $n,  1, 2);
}
sub ramanujan_prime_count {
  my($low,$high) = @_;
  if (defined $high) { validate_integer_nonneg($low); }
  else               { ($low,$high) = (2, $low);      }
  validate_integer_nonneg($high);
  my $sum = 0;
  while ($low <= $high) {
    my $seghigh = ($high-$high) + $low + 1e9 - 1;
    $seghigh = $high if $seghigh > $high;
    $sum += scalar(@{Math::Prime::Util::ramanujan_primes($low,$seghigh)});
    $low = $seghigh + 1;
  }
  $sum;
}

sub twin_prime_count_approx {
  my($n) = @_;
  validate_integer_nonneg($n);
  return twin_prime_count(3,$n) if $n < 2000;
  # Remove bigint / bigfloat.  Everything here will be done with native NV.
  $n = 0.0+"$n" if ref($n);
  my $logn = log($n);
  my $li2 = Math::Prime::Util::ExponentialIntegral($logn) + 2.8853900817779268147198494 - ($n/$logn);

  # Empirical correction factor
  my $fm;
  if    ($n <     4000) { $fm = 0.2952; }
  elsif ($n <     8000) { $fm = 0.3151; }
  elsif ($n <    16000) { $fm = 0.3090; }
  elsif ($n <    32000) { $fm = 0.3096; }
  elsif ($n <    64000) { $fm = 0.3100; }
  elsif ($n <   128000) { $fm = 0.3089; }
  elsif ($n <   256000) { $fm = 0.3099; }
  elsif ($n <   600000) { my($x0, $x1, $y0, $y1) = (1e6, 6e5, .3091, .3059);
                          $fm = $y0 + ($n - $x0) * ($y1-$y0) / ($x1 - $x0); }
  elsif ($n <  1000000) { my($x0, $x1, $y0, $y1) = (6e5, 1e6, .3062, .3042);
                          $fm = $y0 + ($n - $x0) * ($y1-$y0) / ($x1 - $x0); }
  elsif ($n <  4000000) { my($x0, $x1, $y0, $y1) = (1e6, 4e6, .3067, .3041);
                          $fm = $y0 + ($n - $x0) * ($y1-$y0) / ($x1 - $x0); }
  elsif ($n < 16000000) { my($x0, $x1, $y0, $y1) = (4e6, 16e6, .3033, .2983);
                          $fm = $y0 + ($n - $x0) * ($y1-$y0) / ($x1 - $x0); }
  elsif ($n < 32000000) { my($x0, $x1, $y0, $y1) = (16e6, 32e6, .2980, .2965);
                          $fm = $y0 + ($n - $x0) * ($y1-$y0) / ($x1 - $x0); }
  $li2 *= $fm * log(12+$logn)  if defined $fm;

  return int(1.32032363169373914785562422 * $li2 + 0.5);
}

sub semiprime_count_approx {
  my($n) = @_;
  validate_integer_nonneg($n);
  return 0 if $n < 4;
  $n = "$n" + 0.00000001;
  my $l1 = log($n);
  my $l2 = log($l1);
  #my $est = $n * $l2 / $l1;
  #my $est = $n * ($l2 + 0.302) / $l1;
  my $est = ($n/$l1) * (0.11147910114 + 0.00223801350*$l1 + 0.44233207922*$l2 + 1.65236647896*log($l2));
  Mtoint($est + 0.5);
}

sub almost_prime_count_approx {
  my($k,$n) = @_;
  validate_integer_nonneg($k);
  validate_integer_nonneg($n);
  return ($n >= 1) if $k == 0;
  return Math::Prime::Util::prime_count_approx($n) if $k == 1;
  return Math::Prime::Util::semiprime_count_approx($n) if $k == 2;
  return 0 if ($n >> $k) == 0;

  my $lo = Math::Prime::Util::almost_prime_count_lower($k, $n);
  my $hi = Math::Prime::Util::almost_prime_count_upper($k, $n);

  if ($k == 3) {
    my $x  = 0.0 + "$n";
    my $l1 = log($x);
    my $l2 = log($l1);
    my($a,$s) = (1.0,2.0);
    if    ($x <=      638) { $s = 1.554688; $a = 0.865814; }
    elsif ($x <=     1544) { $s = 1.050000; $a = 0.822256; }
    elsif ($x <=     1927) { $s = 0.625000; $a = 0.791747; }
    elsif ($x <=   486586) { $s = 2.865611; $a = 1.004090; }
    elsif ($x <=  1913680) { $s = 2.790963; $a = 0.999618; }
    elsif ($x <= 22347532) { $s = 2.719238; $a = 0.995635; }
    elsif ($x <= 2.95e8)   { $s = 2.584473; $a = 0.988802; }
    elsif ($x <= 4.20e9)   { $s = 2.457108; $a = 0.983098; }
    elsif ($x <= 7.07e10)  { $s = 2.352818; $a = 0.978931; }
    elsif ($x <= 1.36e12)  { $s = 2.269745; $a = 0.975953; }
    elsif ($x <= 4.1e13)   { $s = 2.203002; $a = 0.973796; }
    elsif ($x <= 9.2e14)   { $s = 2.148463; $a = 0.972213; }
    else                   { $s = 2.119279; $a = 0.971438; }
    my $est = 0.5*$a*$x*(($l2+0.26153)*($l2+0.26153)) / ($l1+$s) + 0.5;
    return $est < $lo ? $lo : $est > $hi ? $hi : Mtoint($est);
  }

  {
    my $mult = 0.5;
    if ($n < 2**32 && $k < 13) {
      $mult = 0.9;
    } elsif ($k > 11) {
      $mult = 0.20;
    } else {
      $mult = 0.76;
    }
    return Mtoint($lo + ($hi - $lo) * $mult + 0.5) unless ref($lo) || ref($hi);

    my $imult = int($mult * (1<<16));
    my $est = Maddint( Mlshiftint($lo,16), Mmulint(Msubint($hi,$lo),$imult) );
    return Mrshiftint($est,16);
  }
}

sub nth_twin_prime {
  my($n) = @_;
  return undef if $n < 0;  ## no critic qw(ProhibitExplicitReturnUndef)
  return (undef,3,5,11,17,29,41)[$n] if $n <= 6;

  my $p = Math::Prime::Util::nth_twin_prime_approx($n+200);
  my $tp = Math::Prime::Util::twin_primes($p);
  while ($n > scalar(@$tp)) {
    $n -= scalar(@$tp);
    $tp = Math::Prime::Util::twin_primes($p+1,$p+1e5);
    $p += 1e5;
  }
  return $tp->[$n-1];
}

sub nth_twin_prime_approx {
  my($n) = @_;
  validate_integer_nonneg($n);
  return nth_twin_prime($n) if $n < 6;
  $n = _upgrade_to_float($n) if ref($n) || $n > 127e14;   # TODO lower for 32-bit
  my $logn = log($n);
  my $nlogn2 = $n * $logn * $logn;

  return int(5.158 * $nlogn2/log(9+log($n*$n))) if $n > 59 && $n <= 1092;

  my $lo = int(0.7 * $nlogn2);
  my $hi = int( ($n > 1e16) ? 1.1 * $nlogn2
              : ($n >  480) ? 1.7 * $nlogn2
                            : 2.3 * $nlogn2 + 3 );

  _binary_search($n, $lo, $hi,
                 sub{Math::Prime::Util::twin_prime_count_approx(shift)},
                 sub{ ($_[2]-$_[1])/$_[1] < 1e-15 } );
}

sub nth_semiprime {
  my $n = shift;
  validate_integer_nonneg($n);
  return (undef,4,6,9,10,14,15,21,22)[$n] if $n <= 8;
  my $x = "$n" + 0.000000001; # Get rid of bigint so we can safely call log
  my $logx = log($x);
  my $loglogx = log($logx);
  my $a = ($n < 1000) ? 1.027 : ($n < 10000) ? 0.995 : 0.966;
  my $est = $a * $x * $logx / $loglogx;
  my $lo = ($n < 20000) ? int(0.97*$est)-1 : int(0.98*$est)-1;
  my $hi = ($n < 20000) ? int(1.07*$est)+1 : int(1.02*$est)+1;
  1+_binary_search($n,$lo,$hi, sub{Math::Prime::Util::semiprime_count(shift)});
}

sub nth_semiprime_approx {
  my $n = shift;
  validate_integer_nonneg($n);
  return (undef,4,6,9,10,14,15,21,22)[$n] if $n <= 8;
  $n = "$n" + 0.00000001;
  my $l1 = log($n);
  my $l2 = log($l1);
  my $est = 0.966 * $n * $l1 / $l2;
  Mtoint($est+0.5);
}

sub _almost_prime_count_asymptotic {
  my($k, $n) = @_;
  return 0 if ($n >> $k) == 0;
  return ($n >= 1) if $k == 0;

  my $x;
  if (ref($n) || $n > ~0) {
    $x = _upgrade_to_float($n);
  } else {
    $x = 0.0 + "$n";
  }
  my $logx = log($x);
  my $loglogx = log($logx);
  my $est = $x / $logx;
  my $numk = $k - ( ($k<7) ? 1 : ($k<12) ? 2 : ($k-6)>>2 );
  $est *= ($loglogx/$_) for 1 .. $numk;
  $est;  # Returns FP
}
sub _almost_prime_nth_asymptotic {
  my($k, $n) = @_;
  return 0 if $k == 0 || $n == 0;
  return Mpowint(2,$k) if $n == 1;

  my $x;
  if (ref($n) || $n > ~0) {
    require Math::BigFloat;
    Math::BigFloat->import();
    $x = Math::BigFloat->new($n);
  } else {
    $x = 0.0 + "$n";
  }
  my $logx = log($x);
  my $loglogx = log($logx);
  my $est = $x * $logx;
  my $numk = $k - ( ($k<7) ? 1 : ($k<12) ? 2 : ($k-6)>>2 );
  $est *= ($_/$loglogx) for 1 .. $numk;
  $est;  # Returns FP
}

sub almost_prime_count_lower {
  my($k, $n) = @_;
  validate_integer_nonneg($k);
  validate_integer_nonneg($n);


  return 0 if ($n >> $k) == 0;
  ($k, $n) = _kap_reduce_count($k, $n);
  return ($n >= 1) if $k == 0;
  return Math::Prime::Util::prime_count_lower($n) if $k == 1;

  my $bound = 0;
  my $x = 0.0 + "$n";
  my $logx = log($x);
  my $loglogx = log($logx);
  my $logplus = $loglogx + 0.26153;

  my @lower20 = (0,0, 0.8197, 0.8418, 0.5242, 0.5154,0.3053,0.1901,0.1253,0.0892,0.06551,0.05082,0.04101);
  my @lower32 = (0,0, 1.004,  0.7383, 0.6828, 0.5939,0.3594,0.2222,0.1438,0.09754,0.06981,0.05245,0.04151, 0.03461,0.03006,0.02709,0.02553,0.02502,0.02552,0.02697,0.02945);
  my @lower64 = (0,0,1.011,0.8093,0.7484,0.6465,0.3982,0.2463,0.1571,0.1048,0.07363,0.0545,0.0422, 0.0331,0.0270,0.0232,0.0208,0.0194,0.0190,0.0193,0.0203, 0.0222,0.0252,0.0295,0.0356,0.0444,0.0570,0.0753,0.102,0.14,0.20,0.297,0.44,0.68,1.07,1.71,2.8,4.7,8.0,13.89,23.98);
  # TODO: These are likely still too high
  my @lower   = (0,0,1.011,0.8093,0.7484,0.6465,0.3982,0.2463,0.1571,0.1048,0.07363,0.0545,0.0422, 0.0331,0.0270,
 0.0230,0.0200,0.0187,0.018,0.018,0.019,0.020,0.020,0.027,0.032,0.040,0.051,0.068,0.090,0.12,0.18,0.26,0.355);

  my $multl;
  my $isn64bit = Mrshiftint($n,64) == 0;
  if    ($n <=    1048575) { $multl = $lower20[$k]; }
  elsif ($n <= 4294967295) { $multl = $lower32[$k]; }
  elsif ($isn64bit)        { $multl = $lower64[$k]; }
  else {
    push @lower, 1.5 * $lower[$#lower] until defined $lower[$k];
    $multl = $lower[$k];
  }

  if ($k == 2) {
    if ($x <= 1e12) {
      $bound = $x * ($loglogx + 0.261536) / $logx;
    } else {
      # Bayless Theorem 5.2
      $bound = ($x * ($loglogx+0.1769)/$logx) * (1 + 0.4232/$logx);
      $multl = 1;
    }
  } elsif ($k == 3) {
    # Kinlaw Theorem 1, using custom multipliers for 64-bit n
    $bound = $x * $loglogx * $loglogx / (2*$logx);
    if ($n < 638) {
      $multl = 0.8418;
    } elsif ($n < 1927) {
      my $dist = ($x - 638) / (1926 - 638);
      $multl = (1.0-$dist) * 0.8939  +  $dist * 0.9233;
    } elsif ($n < 500194) {
      my $dist = ($x - 1927) / (500194 - 1927);
      $multl = (1.0-$dist) * 0.9233  +  $dist * 1.000;
    } elsif ($n <= 3184393786) {
      my $dist = ($x - 500194) / (3184393786 - 500194);
      $multl = (1.0-$dist) * 1.0000  +  $dist * 1.039;
    } else {
      $multl = $isn64bit ? 1.0004 : 1.0;
    }
  } elsif ($k == 4) {
    $bound = $x * $logplus*$logplus*$logplus / (6*$logx);
    $multl = 0.4999 if !$isn64bit;
  } else {
    $bound = $x / $logx;
    $logplus = $loglogx + (log("$k") * log(log("$k")) - 0.504377);
    $bound *= $logplus/$_ for 1 .. $k-1;
  }
  $bound *= $multl;
  $bound = 1 if $bound < 1;  # We would have returned zero earlier
  Mtoint($bound)
}
sub almost_prime_count_upper {
  my($k, $n) = @_;
  validate_integer_nonneg($k);
  validate_integer_nonneg($n);

  return 0 if ($n >> $k) == 0;
  ($k, $n) = _kap_reduce_count($k, $n);
  return ($n >= 1) if $k == 0;
  return Math::Prime::Util::prime_count_upper($n) if $k == 1;

  # In theory we might have reduced k/n to where XS can handle it.
  # We should consider handling that, especially for k >= 5.

  my $bound = 0;
  my $x = 0.0 + "$n";
  my $logx = log($x);
  my $loglogx = log($logx);
  my $logplus = $loglogx + 0.26153;

  my @upper20 = (0,0, 1.006,0.7385,0.6830,0.5940,0.3596,0.2227,0.1439, 0.09785,0.07016,0.05303,0.04202);
  my @upper32 = (0,0, 1.013,0.8094,0.7485, 0.6467,0.3984,0.2464,0.1572,0.1049,0.07364,0.05452,0.04266, 0.03542,0.03082,0.02798,0.02642,0.02585,0.02615,0.02808,0.03054);
  my @upper64 = (0,0, 1.028, 1.028, 1.3043,
  0.72208, 0.46609, 0.29340,0.18571,0.12063,0.0815,0.0575,0.0427,
  0.03490, 0.03007, 0.02710, 0.02554, 0.02504, 0.02554, 0.02699, 0.02954,
  0.03294, 0.03779, 0.04453, 0.05393, 0.06703, 0.08543, 0.1117, 0.1494,
  0.205,0.287,0.410,
  0.60,0.90,1.36,2.12,3.35,5.38,8.83,14.75,25.07);

  # TODO: These bounds are likely to not be accurate for large inputs

  my $multu;
  my $isn64bit = Mrshiftint($n,64) == 0;
  if    ($n <=    1048575) { $multu = $upper20[$k]; }
  elsif ($n <= 4294967295) { $multu = $upper32[$k]; }
  else {
    push @upper64, 2.1 * $upper64[$#upper64] until defined $upper64[$k];
    $multu = $upper64[$k];
  }

  if ($k == 2) {
    # Bayless Corollary 5.1
    $bound = 1.028 * $x * ($loglogx + 0.261536) / $logx;
  } elsif ($k == 3) {
    # Bayless Theorem 5.3
    $bound = $x * ($logplus * $logplus + 1.055852) / (2*$logx);
    $multu = 0.8711 if $n > 4294967295 && $isn64bit;
  } elsif ($k == 4) {
    # Bayless Theorem 5.4 part 1, multu = 1.3043
    $bound = $x * $logplus*$logplus*$logplus / (6*$logx);
    if ($x >= 1e12) {  # part 2
      $bound += + 0.511977 * $x * (log(log($x/4)) + 0.261536) / $logx;
      $multu = 1.028;
    }
    if ($isn64bit) {
      $multu = 0.780  if $n > 4294967295;
      $multu = 0.6921 if $x > 1e12;
    }
  } else {
    # We could use Bayless (2018) Theorem 3.5:
    #  # First we have Pi_k(x) -- the upper bound for the square free kaps.
    #   $bound = 1.028 * $x / $logx;
    #   $bound *= ($logplus/$_) for 1..$k-1;
    #  # Second, turn into Tau_k(x) using the paragraph before Theorem 5.4.
    #   my $sigmalim = Msqrtint(Mdivint($n, Mpowint(2,$k-2)));
    #   my $ix = Math::BigInt->new("$x");
    #   Mforprimes( sub {
    #     $bound += almost_prime_count_upper($k-2, Mdivint($ix,Mmulint($_,$_)));
    #   }, 2, $sigmalim);
    #  # This is incredibly slow. )
    #
    # Or use theorem 1 from:
    #   Erdős & Sárközy, "On the number of prime factors of integers", 1980.
    #
    # Or Hildebrand & Tenenbaum 1988:
    #   https://www.researchgate.net/publication/38333551_On_the_number_of_prime_factors_of_an_integer
    # Section 1 has lots of info.  Corollary 2 (page 476) is what we want.

    $bound = $x / $logx;
    $logplus = $loglogx + (log("$k") * log(log("$k")) - 0.504377);
    $bound *= $logplus/$_ for 1 .. $k-1;
  }

  $bound *= $multu;
  $bound = 1 if $bound < 1;  # We would have returned zero earlier
  Mtoint($bound + 1)
}

sub _kap_reduce_nth {   # returns reduction amount r
  my($k, $n) = @_;
  return 0 if $k <= 1;

  # We could calculate new values as needed.
  my @A078843 = (1, 2, 3, 5, 8, 14, 23, 39, 64, 103, 169, 269, 427, 676, 1065, 1669, 2628, 4104, 6414, 10023, 15608, 24281, 37733, 58503, 90616, 140187, 216625, 334527, 516126, 795632, 1225641, 1886570, 2901796, 4460359, 6851532, 10518476, 16138642, 24748319, 37932129, 58110457, 88981343, 136192537, 208364721, 318653143, 487128905, 744398307, 1137129971, 1736461477, 2650785552, 4045250962, 6171386419, 9412197641, 14350773978, 21874583987, 33334053149, 50783701654, 77348521640, 117780873397, 179306456282, 272909472119, 415284741506);
  my $r = 0;
  if ($k > $#A078843) {
    return 0 if $n >= $A078843[-1];
    $r = $k - $#A078843;
  }
  $r++ while $n < $A078843[$k-$r];
  $r;
}
sub _fast_small_nth_almost_prime {
  my($k,$n) = @_;
  croak "Internal kap out of range error" if $n >= 8 || $k < 2;
  return (0, 4,  6,  9, 10, 14, 15, 21)[$n] if $k == 2;
  return Mmulint((0, 8, 12, 18, 20, 27, 28, 30)[$n], Mlshiftint(1,$k-3));
}

sub nth_almost_prime_upper {
  my($k, $n) = @_;
  return undef if $n == 0;
  return (($n == 1) ? 1 : 0) if $k == 0;
  return Mnth_prime_upper($n) if $k == 1;
  return _fast_small_nth_almost_prime($k,$n) if $n < 8;

  my $r = _kap_reduce_nth($k,$n);
  if ($r > 0) {
    my $nth = Math::Prime::Util::nth_almost_prime_upper($k-$r, $n);
    return Mlshiftint($nth, $r);
  }

  my $lo = Mlshiftint(5,$k);   # $k >= 1, $n >= 8
  my $hi = Mtoint(1 + _almost_prime_nth_asymptotic($k, $n));
  # We just guessed at hi, so bump it up until it's in range
  my $rhi = almost_prime_count_lower($k, $hi);
  while ($rhi < $n) {
    $lo = Maddint($hi,1);
    $hi = Mvecsum($hi, int(1.02 * ("$hi"/"$rhi") * ("$n"-"$rhi")), 100);
    $rhi = almost_prime_count_lower($k, $hi);
  }
  while ($lo < $hi) {
    my $mid = $lo + (($hi-$lo) >> 1);
    if (almost_prime_count_lower($k,$mid) < $n) { $lo = $mid+1; }
    else                                        { $hi = $mid; }
  }
  $lo;
}
sub nth_almost_prime_lower {
  my($k, $n) = @_;
  return undef if $n == 0;
  return (($n == 1) ? 1 : 0) if $k == 0;
  return Math::Prime::Util::nth_prime_lower($n) if $k == 1;
  return _fast_small_nth_almost_prime($k,$n) if $n < 8;

  my $r = _kap_reduce_nth($k,$n);
  if ($r > 0) {
    my $nth = Math::Prime::Util::nth_almost_prime_lower($k-$r, $n);
    return Mlshiftint($nth, $r);
  }

  my $lo = Mlshiftint(5,$k);   # $k >= 1, $n >= 8
  my $hi = Mtoint(1 + _almost_prime_nth_asymptotic($k, $n));
  # We just guessed at hi, so bump it up until it's in range
  my $rhi = almost_prime_count_upper($k, $hi);
  while ($rhi < $n) {
    $lo = Maddint($hi,1);
    $hi = Mvecsum($hi, int(1.02 * ("$hi"/"$rhi") * ("$n"-"$rhi")), 100);
    $rhi = almost_prime_count_upper($k, $hi);
  }
  while ($lo < $hi) {
    my $mid = $lo + (($hi-$lo) >> 1);
    if (almost_prime_count_upper($k,$mid) < $n) { $lo = $mid+1; }
    else                                        { $hi = $mid; }
  }
  $lo;
}

sub nth_almost_prime_approx {
  my($k, $n) = @_;
  return undef if $n == 0;
  return Mlshiftint(1,$k) if $n == 1;
  return undef if $k == 0;  # n==1 already returned
  return Math::Prime::Util::nth_prime_approx($n) if $k == 1;
  return Math::Prime::Util::nth_semiprime_approx($n) if $k == 2;
  return _fast_small_nth_almost_prime($k,$n) if $n < 8;

  my $r = _kap_reduce_nth($k,$n);
  if ($r > 0) {
    my $nth = Math::Prime::Util::nth_almost_prime_approx($k-$r, $n);
    return Mmulint($nth, Mpowint(2,$r));
  }

  my $lo = Math::Prime::Util::nth_almost_prime_lower($k, $n);
  my $hi = Math::Prime::Util::nth_almost_prime_upper($k, $n);

  # TODO: Add interpolation speedup steps

  while ($lo < $hi) {
    my $mid = $lo + (($hi-$lo) >> 1);
    if (almost_prime_count_approx($k,$mid) < $n) { $lo = $mid+1; }
    else                                         { $hi = $mid; }
  }
  $lo;
}

sub _interp_linear {
  my($n, $rlo, $rhi, $lo, $hi) = @_;
  #return int( ($n-$rlo) * ($hi-$lo) / ($rhi-$rlo) );
  my $num = Mmulint( Msubint($n,$rlo), Msubint($hi,$lo) );
  my $den = Msubint($rhi, $rlo);
  return Mdivint(Maddint($num,$den>>1), $den);
  #return divint($num, $den);
}
sub _inverse_interpolate {
  my($lo, $hi, $n, $k, $callback) = @_;
  my($mid, $rmid, $rlo, $rhi);

  $rlo = $callback->($k, $lo);
  croak "interp: bad lower bound" if $rlo > $n;
  return $lo if $rlo == $n;   # If lo wasn't small enough, this could be wrong.

  # We have the exact value (rlo) at lo.
  #print "1   $lo  $hi   ",$hi-$lo,"\n";

  $rhi = $callback->($k, $hi) if $hi != 0;

  while ($hi == 0) {
    # Use lo/rlo to make an estimate

    # Make an estimate of where we will end up
    my $estf = ($rlo == 0) ? 1 : Mdivint(Mlshiftint($n,8),$rlo) - 1; # slightly lower
    $estf = 1+(1<<8) if $estf <= (1<<8);
    $estf = (8<<8) if $estf > (8<<8);
    $mid = Mrshiftint(Mmulint($estf,$lo),8);
    # rmid is the exact count at this estimate
    $rmid = $callback->($k, $mid);

    # Either we have a hi value, or we pull in lo and do it again.
    if ($rmid >= $n) { $hi = $mid;  $rhi = $rmid; }
    else             { $lo = $mid;  $rlo = $rmid; }
    #print "2   $lo  $hi   ",$hi-$lo,"\n";
  }
  croak "interp bad initial" unless $rlo <= $n && $rhi >= $n;
  return $lo if $rlo == $n;
  return (($rlo==$n || ($rlo<$n && $rhi>$n)) ? $lo : $hi) if $hi-$lo <= 1;

  # Step 1.  Linear interpolation while it centers.

  $mid = ($n == $rhi)  ?  $hi-1
                       :  Maddint($lo, _interp_linear($n,$rlo,$rhi,$lo,$hi));
  if ($mid == $lo) { $mid++; } elsif ($mid == $hi) { $mid--; }

  while ($rhi > $n && ($hi-$lo) > 1) {
    croak "interp: need 3 unique points" unless $lo < $mid && $mid < $hi;
    #print "I   $lo  $hi   ",$hi-$lo,"\n";
    $rmid = $callback->($k, $mid);
    if ($rmid >= $n) { ($hi,$rhi) = ($mid,$rmid); }
    else             { ($lo,$rlo) = ($mid,$rmid); }
    last if $rhi == $n;

    my $num = Mmulint(Msubint($n,$rmid),Msubint($hi,$lo));
    my $den = Msubint($rhi,$rlo);
    $mid = Maddint($mid, Mdivint($num, $den));
    # Fairly crude way of pulling in opposite side so we bracket.
    if    ($mid <= $lo) { $mid = Maddint($lo, Mdivint(Msubint($hi,$lo),100)); }
    elsif ($mid >= $hi) { $mid = Msubint($hi, Mdivint(Msubint($hi,$lo),100)); }
    if ($mid == $lo) { $mid++; } elsif ($mid == $hi) { $mid--; }
    croak "interp: range error" unless $lo <= $mid && $mid <= $hi;
  }

  return $lo if $rlo == $n;
  return (($rlo==$n || ($rlo<$n && $rhi>$n)) ? $lo : $hi) if $hi-$lo <= 1;

  croak "interp: bad step 1 interpolation" unless $rlo < $n && $rhi == $n;

  # Step 2.   Ridder's method until we're very close.

  croak "interp: Ridder initial assumption error" unless $rlo<$n && $rhi>=$n;
  #print "F   $lo  $hi   ",$hi-$lo,"\n";

  while (($hi-$lo > 8) && ($hi-$lo) > 1) {
    my($x0, $x2, $x1) = ($lo, $hi, Maddint($lo, Msubint($hi,$lo)>>1));
    my($rx1) = $callback->($k, $x1);
    my($fx0, $fx1, $fx2) = (Msubint($rlo,$n), Msubint($rx1,$n), Msubint($rhi,$n)+1);

    # Calculate new point using false position method
    #my $pos = (($x1-$x0) * "$fx1") / sqrt( "$fx1"*"$fx1" - "$fx0"*"$fx2" );
    #my $x3 = $x1 - int($pos+0.5);
    # Rather convoluted so it's all in integer.
    my $num = Mmulint($fx1, Msubint($x1,$x0));
    my $d1  = Msubint(Mmulint($fx1,$fx1),Mmulint($fx0,$fx2));
    my $den = Msqrtint(Mlshiftint($d1,64));
       $num = Mlshiftint($num, 32);
    my $pos = Mdivint(Maddint($num,$den>>1), $den);
    my $x3 = Msubint($x1, $pos);

    # print "    Ridder mid = $x1 - $pos = $x3\n";
    # print " $lo $x1 $x3 $hi\n";

    if ($x3 >= $hi || $x3 <= $lo || $x3 == $x1) {

      # The new point hasn't given us anything.  Just bisect.
      if ($rx1 >= $n) { $hi = $x1; $rhi = $rx1; }
      else            { $lo = $x1; $rlo = $rx1; }

    } else {

      my $rx3 = $callback->($k,$x3);
      if ($rx1 > $rx3) { ($x1,$x3,$rx1,$rx3) = ($x3,$x1,$rx3,$rx1); }
      if    ($rx1 >= $n) {                           $hi = $x1;  $rhi = $rx1; }
      elsif ($rx3 >= $n) { $lo = $x1;  $rlo = $rx1;  $hi = $x3;  $rhi = $rx3; }
      else               { $lo = $x3;  $rlo = $rx3;                           }

    }
    #print "R   $lo  $hi   ",$hi-$lo,"\n";
    croak "interp: Ridder step error" unless $rlo < $n && $rhi >= $n;
  }

  # Step 3.  Binary search.  Invariant:  f(lo) < n, f(hi) >= n

  while ($hi-$lo > 1) {
    $mid = $lo + (($hi-$lo) >> 1);
    $rmid = $callback->($k, $mid);
    if ($rmid < $n) { $lo = $mid; }
    else            { $hi = $mid; }
    #print "B  $lo  $hi   ",$hi-$lo,"\n";
  }
  $hi;
}

sub nth_almost_prime {
  my($k, $n) = @_;
  return undef if $n == 0;
  return Mlshiftint(1,$k) if $n == 1;
  return undef if $k == 0;  # n==1 already returned
  return Math::Prime::Util::nth_prime($n) if $k == 1;
  return Math::Prime::Util::nth_semiprime($n) if $k == 2;
  return _fast_small_nth_almost_prime($k,$n) if $n < 8;

  my $r = _kap_reduce_nth($k,$n);
  if ($r > 0) {
    my $nth = Math::Prime::Util::nth_almost_prime($k-$r, $n);
    return Mmulint($nth, Mpowint(2,$r));
  }

  my $lo = Math::Prime::Util::nth_almost_prime_lower($k, $n);

  return _inverse_interpolate($lo, 0, $n, $k, sub { Math::Prime::Util::almost_prime_count($_[0],$_[1]); });
  #my $ncalls = 0;
  #my $res = _inverse_interpolate($lo, 0, $n, $k, sub { $ncalls++; Math::Prime::Util::almost_prime_count($_[0],$_[1]); });
  #print "ncalls:  $ncalls\n";
  #return $ncalls;
  #return $res;
}

sub nth_omega_prime {
  my($k, $n) = @_;
  return undef if $n == 0;
  return Mpn_primorial($k) if $n == 1;
  return undef if $k == 0;  # n==1 already returned

  # Very inefficient algorithm.
  my $i = Mpn_primorial($k);
  while (1) {
    $i++ while Mprime_omega($i) != $k;
    return $i if --$n == 0;
    $i++;
  }
}

sub nth_ramanujan_prime_upper {
  my $n = shift;
  validate_integer_nonneg($n);
  return (0,2,11)[$n] if $n <= 2;

  if ($n < 50) {
    return Mnth_prime_upper(int(2.6*$n)) if $n <= 20;
    return 33+((310*Mnth_prime_upper(2*$n))>>8);
  }

  my $nth = Mnth_prime_upper(Mmulint($n,3));

  return  115+((727*$nth) >> 10) if $n < 647;

  # TODO: Ideally these would all be adjusted to make smooth transitions.

  my($add,$mul) = $n <      16000 ? ( 271,358)
                : $n <    1200000 ? (9450,350)
                : $n <    7000000 ? (5000,349)
                : $n <   90000000 ? (   0,348)
                : $n < 3100000000 ? (   0,347)
                :                   (   0,346);

  my $ret = Mrshiftint(Mmulint($mul,$nth),9);
  $ret = Maddint($ret,$add) if $add != 0;
  $ret;
}
sub nth_ramanujan_prime_lower {
  my $n = shift;
  validate_integer_nonneg($n);
  return (0,2,11)[$n] if $n <= 2;
  my $nth = Math::Prime::Util::nth_prime_lower(Mmulint($n,2));
  return Mdivint(Mmulint(275,$nth),256) if $n < 10000;
  return Mdivint(Mmulint(262,$nth),256) if $n < 1e10;
  $nth;
}
sub nth_ramanujan_prime_approx {
  my $n = shift;
  validate_integer_nonneg($n);
  return (0,2,11)[$n] if $n <= 2;
  my($lo,$hi) = (nth_ramanujan_prime_lower($n),nth_ramanujan_prime_upper($n));
  $lo + (($hi-$lo)>>1);
}
sub ramanujan_prime_count_upper {
  my $n = shift;
  validate_integer_nonneg($n);
  return (($n < 2) ? 0 : 1) if $n < 11;
  my $lo = Mdivint(prime_count_lower($n),3);
  my $hi = Mrshiftint(prime_count_upper($n),1);
  1+_binary_search($n, $lo, $hi,
                   sub{Math::Prime::Util::nth_ramanujan_prime_lower(shift)});
}
sub ramanujan_prime_count_lower {
  my $n = shift;
  validate_integer_nonneg($n);
  return (($n < 2) ? 0 : 1) if $n < 11;
  my $lo = int(prime_count_lower($n) / 3);
  my $hi = prime_count_upper($n) >> 1;
  _binary_search($n, $lo, $hi,
                 sub{Math::Prime::Util::nth_ramanujan_prime_upper(shift)});
}
sub ramanujan_prime_count_approx {
  my $n = shift;
  validate_integer_nonneg($n);
  return (($n < 2) ? 0 : 1) if $n < 11;
  #$n = _upgrade_to_float($n) if ref($n) || $n > 2e16;
  my $lo = ramanujan_prime_count_lower($n);
  my $hi = ramanujan_prime_count_upper($n);
  _binary_search($n, $lo, $hi,
                 sub{Math::Prime::Util::nth_ramanujan_prime_approx(shift)},
                 sub{ ($_[2]-$_[1])/$_[1] < 1e-15 } );
}

sub _sum_primes_n {
  my $n = shift;
  return (0,0,2,5,5)[$n] if $n < 5;
  my $r = Msqrtint($n);
  my $r2 = $r + Mdivint($n, $r+1);
  my(@V,@S);
  for my $k (0 .. $r2) {
    my $v = ($k <= $r) ? $k : Mdivint($n,($r2-$k+1));
    $V[$k] = $v;
    $S[$k] = Maddint(
              Mrshiftint(Mmulint($v, $v-1)),
              $v-1);
  }
  for my $p (2 .. $r) {
    next unless $S[$p] > $S[$p-1];
    my $sp = $S[$p-1];
    my $p2 = Mmulint($p,$p);
    for my $v (reverse @V) {
      last if $v < $p2;
      my($a,$b) = ($v,Mdivint($v,$p));
      $a = $r2 - Mdivint($n,$a) + 1 if $a > $r;
      $b = $r2 - Mdivint($n,$b) + 1 if $b > $r;
      $S[$a] -= Mmulint($p, $S[$b]-$sp);
      #$S[$a] = Msubint($S[$a], Mmulint($p, Msubint($S[$b],$sp)));
    }
  }
  $S[$r2];
}
sub sum_primes {
  my($low,$high) = @_;
  if (defined $high) { validate_integer_nonneg($low); }
  else               { ($low,$high) = (2, $low);      }
  validate_integer_nonneg($high);
  my $sum = 0;

  return $sum if $high < $low;

  # It's very possible we're here because they've counted too high.  Skip fwd.
  if ($low <= 2 && $high >= 29505444491) {
    ($low, $sum) = (29505444503, tobigint("18446744087046669523"));
  }

  return $sum if $low > $high;

  # Easy, not unreasonable, but seems slower than the windowed sum.
  # return _sum_primes_n($high) if $low <= 2;

  # Performance decision, which to use.
  if ( $high <= ~0 &&
       $high > (MPU_64BIT ? 2000000 : 320000) &&
       ($high-$low) > $high/50 &&
       !getconfig()->{'xs'}) {
    my $hsum = _sum_primes_n($high);
    my $lsum = ($low <= 2) ? 0 : _sum_primes_n($low - 1);
    return $hsum - $lsum;
  }

  # Sum in windows.
  # TODO: consider some skipping forward with small tables.
  my $xssum = (MPU_64BIT && $high < 6e14 && getconfig()->{'xs'});
  my $step = ($xssum && $high > 5e13) ? 1_000_000 : 11_000_000;
  Math::Prime::Util::prime_precalc(Msqrtint($high));
  while ($low <= $high) {
    my $next = Maddint($low, $step) - 1;
    $next = $high if $next > $high;
    $sum = Maddint($sum,
            ($xssum) ? Math::Prime::Util::sum_primes($low,$next)
                     : Mvecsum( @{Mprimes($low,$next)} ));
    last if $next == $high;
    $low = Maddint($next,1);
  }
  $sum;
}

sub print_primes {
  my($low,$high,$fd) = @_;
  if (defined $high) { validate_integer_nonneg($low); }
  else               { ($low,$high) = (2, $low);      }
  validate_integer_nonneg($high);

  $fd = fileno(STDOUT) unless defined $fd;
  open(my $fh, ">>&=", $fd);  # TODO .... or die

  if ($high >= $low) {
    my $p1 = $low;
    while ($p1 <= $high) {
      my $p2 = $p1 + 15_000_000 - 1;
      $p2 = $high if $p2 > $high;
      if ($Math::Prime::Util::_GMPfunc{"sieve_primes"}) {
        print $fh "$_\n" for Math::Prime::Util::GMP::sieve_primes($p1,$p2,0);
      } else {
        print $fh "$_\n" for @{Mprimes($p1,$p2)};
      }
      $p1 = $p2+1;
    }
  }
  close($fh);
}


#############################################################################

sub _mulmod {
  my($x, $y, $n) = @_;
  return (($x * $y) % $n) if ($x|$y) < MPU_HALFWORD;
  #return (($x * $y) % $n) if ($x|$y) < MPU_HALFWORD || $y == 0 || $x < int(~0/$y);
  my $r = 0;
  $x %= $n if $x >= $n;
  $y %= $n if $y >= $n;
  ($x,$y) = ($y,$x) if $x < $y;
  if ($n <= (~0 >> 1)) {
    while ($y > 1) {
      if ($y & 1) { $r += $x;  $r -= $n if $r >= $n; }
      $y >>= 1;
      $x += $x;  $x -= $n if $x >= $n;
    }
    if ($y & 1) { $r += $x;  $r -= $n if $r >= $n; }
  } else {
    while ($y > 1) {
      if ($y & 1) { $r = $n-$r;  $r = ($x >= $r) ? $x-$r : $n-$r+$x; }
      $y >>= 1;
      $x = ($x > ($n - $x))  ?  ($x - $n) + $x  :  $x + $x;
    }
    if ($y & 1) { $r = $n-$r;  $r = ($x >= $r) ? $x-$r : $n-$r+$x; }
  }
  $r;
}
sub _addmod {
  my($x, $y, $n) = @_;
  $x %= $n if $x >= $n;
  $y %= $n if $y >= $n;
  if (($n-$x) <= $y) {
    ($x,$y) = ($y,$x) if $y > $x;
    $x -= $n;
  }
  $x + $y;
}

# Note that Perl 5.6.2 with largish 64-bit numbers will break.  As usual.
sub _native_powmod {
  my($n, $power, $m) = @_;
  my $t = 1;
  $n = $n % $m;
  while ($power) {
    $t = ($t * $n) % $m if ($power & 1);
    $power >>= 1;
    $n = ($n * $n) % $m if $power;
  }
  $t;
}

sub _powmod {
  my($n, $power, $m) = @_;
  my $t = 1;

  $n %= $m if $n >= $m;
  if ($m < MPU_HALFWORD) {
    while ($power) {
      $t = ($t * $n) % $m if ($power & 1);
      $power >>= 1;
      $n = ($n * $n) % $m if $power;
    }
  } else {
    while ($power) {
      $t = _mulmod($t, $n, $m) if ($power & 1);
      $power >>= 1;
      $n = _mulmod($n, $n, $m) if $power;
    }
  }
  $t;
}

sub powint {
  my($a, $b) = @_;
  validate_integer($a);
  validate_integer($b);
  return reftyped($_[0], Math::Prime::Util::GMP::powint($a,$b))
    if $Math::Prime::Util::_GMPfunc{"powint"};
  croak "powint: exponent must be >= 0" if $b < 0;

  # Special cases for small a and b
  if ($a >= -1 && $a <= 4) {
    return ($b == 0) ? 1 : 0 if $a == 0;
    return 1 if $a == 1;
    return ($b % 2) ? -1 : 1 if $a == -1;
    if ($a == 2) {
      return ($b < MPU_MAXBITS)  ?  1<<$b  :  Mlshiftint(1,$b);
    }
    if ($a == 4) {
      return 1 << (2*$b) if $b < MPU_MAXBITS/2;
      return Mlshiftint(1,2*$b) if $b < 4000000000;
    }
  }

  return 1 if $b == 0;
  return $a if $b == 1;
  if ($b == 2) {
    return int("$a")*int("$a") if abs($a) <= (MPU_32BIT ? 65535 : 4294967295);
    return Mmulint($a,$a);
  }

  if (!ref($a) && !ref($b) && $b < MPU_MAXBITS) {
    if ($b == 3) {
      return int($a*$a*$a) if $a <= 99999;
      return Mmulint(int($a*$a), $a) if $a <= 31622776;
    } else {
      # Check if inside limit of int on 32-bit
      my $r = $a ** $b;
      return int($r) if $r < 1000000000000000 && $r > -1000000000000000;
      # Try to complete using a single mulint if we can
      $r = $a ** (($b+1)>>1);
      if ($r < 1000000000000000 && $r > -1000000000000000) {
        return Mmulint(int($r), $b&1 ? int($a**($b>>1)) : int($r));
      }
    }
    # Fall through
  }

  return Mmulint(Mmulint($a,$a),$a) if $b == 3;

  my $r = tobigint($a) ** tobigint($b);
  return $r <= INTMAX && $r >= INTMIN  ?  _bigint_to_int($r)  :  $r;
}

sub mulint {
  my($a, $b) = @_;
  validate_integer($a);
  validate_integer($b);
  return 0 if $a == 0 || $b == 0;
  return reftyped($_[0], Math::Prime::Util::GMP::mulint($a,$b))
    if $Math::Prime::Util::_GMPfunc{"mulint"};

  my $r = $a * $b;

  if (!ref($r)) {
    return $r if $r < INTMAX && $r > INTMIN;
    $r = tobigint($a) * $b;
  }
  return $r <= INTMAX && $r >= INTMIN  ?  _bigint_to_int($r)  :  $r;
}
sub addint {
  my($a, $b) = @_;
  validate_integer($a);
  validate_integer($b);
  return reftyped($_[0], Math::Prime::Util::GMP::addint($a,$b))
    if $Math::Prime::Util::_GMPfunc{"addint"};

  my $r = $a + $b;

  if (!ref($r)) {
    return $r if $r < INTMAX && $r > INTMIN;
    $r = tobigint($a) + $b;
  }
  return $r <= INTMAX && $r >= INTMIN  ?  _bigint_to_int($r)  :  $r;
}
sub subint {
  my($a, $b) = @_;
  validate_integer($a);
  validate_integer($b);
  return reftyped($_[0], Math::Prime::Util::GMP::subint($a,$b))
    if $Math::Prime::Util::_GMPfunc{"subint"};

  my $r = $a - $b;

  if (!ref($r)) {
    return $r if $r < INTMAX && $r > INTMIN;
    $r = tobigint($a) - $b;
  }
  return $r <= INTMAX && $r >= INTMIN  ?  _bigint_to_int($r)  :  $r;
}
sub add1int {
  my($a) = @_;
  validate_integer($a);
  my $r = $a+1;
  if (!ref($r)) {
    return $r if $r < INTMAX;
    $r = tobigint($a) + 1;
  }
  return $r <= INTMAX && $r >= INTMIN  ?  _bigint_to_int($r)  :  $r;
}
sub sub1int {
  my($a) = @_;
  validate_integer($a);
  my $r = $a-1;
  if (!ref($r)) {
    return $r if $r < INTMAX;
    $r = tobigint($a) - 1;
  }
  return $r <= INTMAX && $r >= INTMIN  ?  _bigint_to_int($r)  :  $r;
}

# For division / modulo, see:
#
# https://www.researchgate.net/publication/234829884_The_Euclidean_definition_of_the_functions_div_and_mod
#
# https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/divmodnote-letter.pdf

sub _tquotient {
  my($a,$b) = @_;
  return $a if $b == 1;

  $a = tobigint($a) if ($a >= SINTMAX || $a <= INTMIN) && !ref($a);
  $b = tobigint($b) if ($b >= SINTMAX || $b <= INTMIN) && !ref($b);
  my($refa,$refb) = (ref($a),ref($b));

  if (!$refa && !$refb) {
    # Both numbers are in signed range, so we can negate them.
    use integer;  # This is >>> SIGNED <<< integer.
    # Signed division is implementation defined in C89.
    return -(-$a /  $b)  if $a < 0 && $b > 0;
    return -( $a / -$b)  if $b < 0 && $a > 0;
    return  (-$a / -$b)  if $a < 0 && $b < 0;
    return  ( $a /  $b);
  }

  my $q;  # set this, turn into int and return at end
  if ($refa eq 'Math::GMPz' || $refb eq 'Math::GMPz') {
    $q = Math::GMPz->new();
    $a = Math::GMPz->new($a) unless $refa eq 'Math::GMPz';
    $b = Math::GMPz->new($b) unless $refb eq 'Math::GMPz';
    Math::GMPz::Rmpz_tdiv_q($q,$a,$b);
  } elsif ($refa eq 'Math::GMP' || $refb eq 'Math::GMP') {
    $a = Math::GMP->new($a) unless $refa eq 'Math::GMP';
    $b = Math::GMP->new($b) unless $refb eq 'Math::GMP';
    # op_div => mpz_div function (obsolete!).  bdiv => tdiv_qr
    ($q) = $a->bdiv($b);
  } else {
    # Force no upgrade so 'use bignum' won't screw us over.
    my $A = Math::BigInt->new("$a")->upgrade(undef);
    my $B = Math::BigInt->new("$b")->upgrade(undef);
    $q = abs($a) / abs($b);
    $q = -$q if ($a < 0) != ($b < 0);
    $q = $refa->new($q) if $refa ne 'Math::BigInt' && $refb ne 'Math::BigInt';
  }
  $q;
  #return $q <= INTMAX && $q >= INTMIN  ?  _bigint_to_int($q)  :  $q;
}

# Truncated Division
sub tdivrem {
  my($a,$b) = @_;
  validate_integer($a);
  validate_integer($b);
  croak "tdivrem: divide by zero" if $b == 0;
  my($q,$r);
  if (!ref($a) && !ref($b) && $a>=0 && $b>=0 && $a<SINTMAX && $b<SINTMAX) {
    { use integer;  $q = $a / $b; }
    $r = $a - $b * $q;
  } else {
    $q = _tquotient($a, $b);
    $r = $a - $b * $q;
    $q = _bigint_to_int($q) if ref($q) && $q <= INTMAX && $q >= INTMIN;
    $r = _bigint_to_int($r) if ref($r) && $r <= INTMAX && $r >= INTMIN;
  }
  ($q,$r);
}
# Floored Division
sub fdivrem {
  my($a,$b) = @_;
  validate_integer($a);
  validate_integer($b);
  croak "fdivrem: divide by zero" if $b == 0;
  my($q,$r);
  if (!ref($a) && !ref($b) && $a>=0 && $b>=0 && $a<SINTMAX && $b<SINTMAX) {
    use integer; $q = $a / $b;
  } else {
    $q = _tquotient($a, $b);
  }
  $r = $a - $b * $q;
  # qe = qt-I     re = rt+I*d    I = (rt >= 0) ? 0 : (b>0) ? 1 : -1;
  # qf = qt-I     rf = rt+I*d    I = (signum(rt) = -signum(b)) 1 : 0
  if ( ($r < 0 && $b > 0) || ($r > 0 && $b < 0) )
    { $q--; $r += $b; }
  $q = _bigint_to_int($q) if ref($q) && $q <= INTMAX && $q >= INTMIN;
  $r = _bigint_to_int($r) if ref($r) && $r <= INTMAX && $r >= INTMIN;
  ($q,$r);
}
# Ceiling Division
sub cdivrem {
  my($a,$b) = @_;
  validate_integer($a);
  validate_integer($b);
  croak "cdivrem: divide by zero" if $b == 0;
  my($q,$r);
  if (!ref($a) && !ref($b) && $a>=0 && $b>=0 && $a<SINTMAX && $b<SINTMAX) {
    use integer; $q = $a / $b;
  } else {
    $q = _tquotient($a, $b);
  }
  $r = $a - $b * $q;
  if ($r != 0 && (($a >= 0) == ($b >= 0)))
    { $q++; $r -= $b; }
  $q = _bigint_to_int($q) if ref($q) && $q <= INTMAX && $q >= INTMIN;
  $r = _bigint_to_int($r) if ref($r) && $r <= INTMAX && $r >= INTMIN;
  ($q,$r);
}
# Euclidean Division
sub divrem {
  my($a,$b) = @_;
  validate_integer($a);
  validate_integer($b);
  croak "divrem: divide by zero" if $b == 0;
  my($q,$r);
  if (!ref($a) && !ref($b) && $a>=0 && $b>=0 && $a<SINTMAX && $b<SINTMAX) {
    use integer; $q = $a / $b;
  } else {
    $q = _tquotient($a, $b);
  }
  $r = $a - $b * $q;
  if ($r <0) {
    if ($b > 0) { $q--; $r += $b; }
    else        { $q++; $r -= $b; }
  }
  $q = _bigint_to_int($q) if ref($q) && $q <= INTMAX && $q >= INTMIN;
  $r = _bigint_to_int($r) if ref($r) && $r <= INTMAX && $r >= INTMIN;
  ($q,$r);
}

sub divint {
  if ($_[1] > 0 && $_[0] >= 0) {   # Simple no-error all positive case
    my($a,$b) = @_;
    my $q;
    if (!ref($a) && !ref($b) && $a<SINTMAX && $b<SINTMAX) {
      use integer; $q = $a / $b;
    } else {
      $q = _tquotient($a, $b);
      $q = _bigint_to_int($q) if ref($q) && $q <= INTMAX;
    }
    return $q;
  }
  (fdivrem(@_))[0];
}
sub modint {
  if ($_[1] > 0 && $_[0] >= 0) {   # Simple no-error all positive case
    my $r;
    if (ref($_[1]) || ref($_[0])) {
      $r = $_[0] % $_[1];
      $r = _bigint_to_int($r) if $r <= INTMAX;
    } elsif ($_[1] < INTMAX && $_[0] < INTMAX) {
      $r = $_[0] % $_[1];
    } else {
      $r = tobigint($_[0]) % tobigint($_[1]);
      $r = _bigint_to_int($r) if $r <= INTMAX;
    }
    return $r;
  }
  (fdivrem(@_))[1];
}
sub cdivint {
  if ($_[1] > 0 && $_[0] >= 0) {   # Simple no-error all positive case
    my($a,$b) = @_;
    my $q;
    if (!ref($a) && !ref($b) && $a<SINTMAX && $b<SINTMAX) {
      use integer; $q = $a / $b;
      $q++ if $a != $b*$q;
    } else {
      $q = _tquotient($a, $b);
      $q++ if $a != $b*$q;
      $q = _bigint_to_int($q) if ref($q) && $q <= INTMAX;
    }
    return $q;
  }
  (cdivrem(@_))[0];
}

sub absint {
  my($n) = @_;
  validate_integer($n);
  return (($n >= 0) ? $n : -$n) if ref($n);
  $n =~ s/^-// if $n <= 0;
  reftyped($_[0], $n);
}
sub negint {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n == 0;  # Perl 5.6 has to have this: if $n=0 => -$n = -0
  return -$n if ref($n) || $n < (~0 >> 1);
  if ($n > 0) { $n = "-$n"; }
  else        { $n =~ s/^-//; }
  reftyped($_[0], $n);
}
sub signint {
  my($n) = @_;
  validate_integer($n);
  # Native ints, Math::BigInt, and Math::GMP work as-is, but some others don't.
  my $r = $n <=> 0;
  return $r < 0 ? -1 : $r > 0 ? 1 : 0;
}
sub cmpint {
  my($a, $b) = @_;
  validate_integer($a);
  validate_integer($b);
  my $r = $a <=> $b;
  return $r < 0 ? -1 : $r > 0 ? 1 : 0;
}

sub lshiftint {
  my($n, $k) = @_;
  validate_integer($n);
  if (!defined $k) { $k = 1; } else { validate_integer($k); }

  return rshiftint($n, Mnegint($k)) if $k < 0;  # Technically not supported
  return Mnegint(lshiftint(Mnegint($n),$k)) if $n < 0;

  if (!ref($n)) {
    return $n << $k if $n < INTMAX && $k < MPU_MAXBITS && $n == ($n<<$k)>>$k;
    $n = tobigint($n);
  }
  $n = $n << $k;
  return $n <= INTMAX ? _bigint_to_int($n) : $n;

  #my $k2 = (!defined $k) ? 2 : ($k < MPU_MAXBITS) ? (1<<$k) : Mpowint(2,$k);
  #Mmulint($n, $k2);
}
sub rshiftint {
  my($n, $k) = @_;
  validate_integer($n);
  if (!defined $k) { $k = 1; } else { validate_integer($k); }

  return lshiftint($n, Mnegint($k)) if $k < 0;  # Technically not supported
  return Mnegint(rshiftint(Mnegint($n),$k)) if $n < 0;

  if (!ref($n)) {
    # Pre 5.24.0, large right shifts were undefined.
    return $k < MPU_MAXBITS ? $n >> $k : 0  if $n < INTMAX;
    $n = tobigint($n);
  }
  $n = $n >> $k;
  return $n <= INTMAX ? _bigint_to_int($n) : $n;

  #my $k2 = (!defined $k) ? 2 : ($k < MPU_MAXBITS) ? (1<<$k) : Mpowint(2,$k);
  #(Mtdivrem($n, $k2))[0];
}

sub rashiftint {
  my($n, $k) = @_;
  validate_integer($n);
  if (!defined $k) { $k = 1; } else { validate_integer($k); }
  my $k2 = (!defined $k) ? 2 : ($k < MPU_MAXBITS) ? (1<<$k) : Mpowint(2,$k);
  Mdivint($n, $k2);
}

sub powersum {
  my($n, $k) = @_;
  validate_integer_nonneg($n);
  validate_integer_nonneg($k);

  return $n if $n <= 1 || $k == 0;

  return Mdivint(Mvecprod($n, Maddint($n,1), Maddint(Mmulint($n,2),1)),6) if $k==2;
  return Mdivint(Mvecprod(
                          $n, Maddint($n,1), Maddint(Mmulint($n,2),1),
                          Mvecsum( Mmulint(3,Mpowint($n,2)), Mmulint(3,$n), -1 )
                         ),30) if $k==4;

  my $a = Mrshiftint(Mmulint($n,Maddint($n,1)),1);
  return $a if $k == 1;
  return Mmulint($a,$a) if $k == 3;
  return Mdivint(Msubint(Mmulint(4,Mpowint($a,3)),Mmulint($a,$a)),3) if $k == 5;

  if ($k < $n) {
    return Mvecsum( map { Mvecprod( Mfactorial($_),
                                    Mbinomial($n+1,$_+1),
                                    Mstirling($k,$_,2)     ) } 1..$k );
  }

  Mvecsum(map { Mpowint($_,$k) } 1..$n);
}

# Make sure to work around RT71548, Math::BigInt::Lite,
# and use correct lcm semantics.
sub gcd {
  my $REF = undef;
  for my $n (@_) {
    my $refn = ref($n);
    if ($refn) { $REF = $refn; last; }
  }

  # Try all-native if all inputs are native ints.
  if (!$REF) {
    my($x,$y) = (shift || 0, 0);
    $x = -$x if $x < 0;
    while (@_) {
      $y = shift;
      while ($y) {  ($x,$y) = ($y, $x % $y);  }
      $x = -$x if $x < 0;
    }
    return $x;
  }

  my @N = map { ref($_) eq $REF ? $_ : $REF->new("$_") } @_;
  my $gcd;

  if ($REF eq 'Math::BigInt') {
    $gcd = Math::BigInt::bgcd(@N);
  } elsif ($REF eq 'Math::GMPz') {
    $gcd = Math::GMPz->new(shift(@N));
    Math::GMPz::Rmpz_gcd($gcd,$gcd,$_) for @N;
  } elsif ($REF eq 'Math::GMP') {
    $gcd = Math::GMP->new(shift(@N));
    $gcd = Math::GMP::gcd($gcd,$_) for @N;
  } else {
    $gcd = Math::BigInt::bgcd(map { Math::BigInt->new("$_") } @N);
    $gcd = tobigint($gcd);
  }
  $gcd = _bigint_to_int($gcd) if $gcd <= INTMAX;
  $gcd;
}
sub lcm {
  return 1 unless @_;
  return Mabsint($_[0]) if @_==1;

  my $lcm = Mabsint(shift);
  while (@_) {
    my $y = Mabsint(shift) || return 0;
    $lcm = Mmulint($lcm, $y / Mgcd($lcm,$y));
  }
  return $lcm;
}
sub gcdext {
  my($x,$y) = @_;
  if ($x == 0) { return (0, (-1,0,1)[($y>=0)+($y>0)], abs($y)); }
  if ($y == 0) { return ((-1,0,1)[($x>=0)+($x>0)], 0, abs($x)); }

  if ($Math::Prime::Util::_GMPfunc{"gcdext"}) {
    my($a,$b,$g) = Math::Prime::Util::GMP::gcdext($x,$y);
    $a = reftyped($_[0], $a);
    $b = reftyped($_[0], $b);
    $g = reftyped($_[0], $g);
    return ($a,$b,$g);
  }

  my($a,$b,$g,$u,$v,$w);
  if (abs($x) < (~0>>1) && abs($y) < (~0>>1)) {
    $x = _bigint_to_int($x) if ref($x);
    $y = _bigint_to_int($y) if ref($y);
    ($a,$b,$g,$u,$v,$w) = (1,0,$x,0,1,$y);
    while ($w != 0) {
      my $r = $g % $w;
      my $q = int(($g-$r)/$w);
      ($a,$b,$g,$u,$v,$w) = ($u,$v,$w,$a-$q*$u,$b-$q*$v,$r);
    }
  } else {
    ($a,$b,$g,$u,$v,$w) = (1,0,$x,0,1,$y);
    while ($w != 0) {
      my($q,$r) = Mdivrem($g,$w);
      ($a,$b,$g,$u,$v,$w) = ($u, $v, $w,
                             Msubint($a,Mmulint($q,$u)),
                             Msubint($b,Mmulint($q,$v)), $r);
    }
  }
  if ($g < 0) { ($a,$b,$g) = (-$a,-$b,-$g); }
  return ($a,$b,$g);
}

sub chinese2 {
  return (0,0) unless scalar @_;
  my($lcm, $sum);

  if ($Math::Prime::Util::_GMPfunc{"chinese2"} && $Math::Prime::Util::GMP::VERSION >= 0.53) {
    ($sum,$lcm) = Math::Prime::Util::GMP::chinese2(@_);
    if (defined $sum) {
      $sum = tobigint($sum);
      $sum = _bigint_to_int($sum) if ref($sum) && $sum <= INTMAX;
    }
    if (defined $lcm) {
      $lcm = tobigint($lcm);
      $lcm = _bigint_to_int($lcm) if ref($lcm) && $lcm <= INTMAX;
    }
    return ($sum,$lcm);
  }

  # Validate, copy, and do abs on the inputs.
  my @items;
  foreach my $aref (@_) {
    die "chinese arguments are two-element array references"
      unless ref($aref) eq 'ARRAY' && scalar @$aref == 2;
    my($a,$n) = @$aref;
    validate_integer($a);
    validate_integer($n);
    return (undef,undef) if $n == 0;
    $n = Mabsint($n);
    $a = Mmodint($a,$n);
    push @items, [$a,$n];
  }
  return @{$items[0]} if scalar @items == 1;
  @items = sort { $b->[1] <=> $a->[1] } @items;

  ($sum, $lcm) = @{shift @items};

  foreach my $aref (@items) {
    my($ai, $ni) = @$aref;
    # gcdext
    my($u,$v,$g,$s,$t,$w) = (1,0,$lcm,0,1,$ni);
    while ($w != 0) {
      my($q,$r) = Mdivrem($g,$w);
      ($u,$v,$g,$s,$t,$w) = ($s, $t, $w,
                             Msubint($u,Mmulint($q,$s)),
                             Msubint($v,Mmulint($q,$t)), $r);
    }
    ($u,$v,$g) = (-$u,-$v,-$g)  if $g < 0;
    return (undef,undef) if $g != 1 && ($sum % $g) != ($ai % $g); # Not co-prime
    $s = -$s if $s < 0;
    $t = -$t if $t < 0;
    $lcm = Mmulint($lcm, $s);
    $u = Maddint($u, $lcm) if $u < 0;
    $v = Maddint($v, $lcm) if $v < 0;
    my $vs = Mmulmod($v,$s,$lcm);
    my $ut = Mmulmod($u,$t,$lcm);
    my $m1 = Mmulmod($sum,$vs,$lcm);
    my $m2 = Mmulmod($ut,$ai,$lcm);
    $sum = Maddmod($m1, $m2, $lcm);
  }
  ($sum,$lcm);
}

sub chinese {
  (chinese2(@_))[0];
}

sub _from_128 {
  my($hi, $lo) = @_;
  return 0 unless defined $hi && defined $lo;
  Maddint(tobigint("$hi") << MPU_MAXBITS, $lo);
}

sub vecsum {
  return reftyped($_[0], @_ ? $_[0] : 0)  if @_ <= 1;

  return reftyped($_[0], Math::Prime::Util::GMP::vecsum(@_))
    if $Math::Prime::Util::_GMPfunc{"vecsum"};
  my $sum = 0;
  foreach my $v (@_) {
    $sum += $v;
    if ($sum > (INTMAX-250) || $sum < (INTMIN+250)) {
      # Sum again from the start using bigint sum
      $sum = 0;
      $sum = Maddint($sum,$_) for @_;
      return $sum;
    }
  }
  $sum;
}

sub _product_mulint {
  my($a, $b, $r) = @_;
  return $r->[$a] if $b <= $a;
  return Mmulint($r->[$a], $r->[$b]) if $b == $a+1;
  return Mmulint(Mmulint($r->[$a], $r->[$a+1]), $r->[$a+2]) if $b == $a+2;
  my $c = $a + (($b-$a+1)>>1);
  Mmulint( _product_mulint($a, $c-1, $r),  _product_mulint($c, $b, $r) );
}
sub _product_mult {
  my($a, $b, $r) = @_;
  return $r->[$a] if $b <= $a;
  return $r->[$a] * $r->[$a+1] if $b == $a+1;
  return $r->[$a] * $r->[$a+1] * $r->[$a+2] if $b == $a+2;
  my $c = $a + (($b-$a+1)>>1);
  _product_mult($a, $c-1, $r) * _product_mult($c, $b, $r);
}

sub vecprod {
  return 1 unless @_;
  return reftyped($_[0], Math::Prime::Util::GMP::vecprod(@_))
    if $Math::Prime::Util::_GMPfunc{"vecprod"};

  return $_[0] if @_ == 1;

  # Try native for non-negative/non-zero inputs
  if ($_[0] > 0 && $_[0] <= INTMAX && $_[1] > 0 && $_[1] <= INTMAX) {
    my $prod = shift @_;
    $prod *= shift @_
      while @_ && $_[0] > 0 && $_[0] <= INTMAX && int(INTMAX/$prod) > $_[0];
    return $prod if @_ == 0;
    unshift @_, $prod if $prod > 1;
  }

  return mulint($_[0], $_[1]) if @_ == 2;

  # Product tree:
  #
  # my $prod = _product_mulint(0, $#_, \@_);
  # $prod = _bigint_to_int($prod) if ref($prod) && $prod <= INTMAX && $prod >= INTMIN;
  # Or faster:
  my $prod = _product_mult(0, $#_, [map { tobigint($_) } @_]);
  $prod = _bigint_to_int($prod) if $prod <= INTMAX && $prod >= INTMIN;
  $prod;
}

sub vecmin {
  return unless @_;
  my $min = shift;
  for (@_) { $min = $_ if $_ < $min; }
  $min;
}
sub vecmax {
  return unless @_;
  my $max = shift;
  for (@_) { $max = $_ if $_ > $max; }
  $max;
}

sub vecextract {
  my($aref, $mask) = @_;

  return @$aref[@$mask] if ref($mask) eq 'ARRAY';

  # This is concise but very slow.
  # map { $aref->[$_] }  grep { $mask & (1 << $_) }  0 .. $#$aref;

  my($i, @v) = (0);
  while ($mask) {
    push @v, $i if $mask & 1;
    $mask >>= 1;
    $i++;
  }
  @$aref[@v];
}

sub vecequal {
  my($aref, $bref) = @_;
  croak "vecequal element not scalar or array reference"
    unless ref($aref) eq 'ARRAY' && ref($bref) eq 'ARRAY';
  return 0 unless $#$aref == $#$bref;
  my $i = 0;
  for my $av (@$aref) {
    my $bv = $bref->[$i++];
    next if !defined $av && !defined $bv;
    return 0 if !defined $av || !defined $bv;
    if ( ref($av) && ref($bv) &&
         (ref($av) =~ /^(ARRAY|HASH|CODE|FORMAT|IO|REGEXP)$/i) ||
         (ref($bv) =~ /^(ARRAY|HASH|CODE|FORMAT|IO|REGEXP)$/i) ) {
      next if (ref($av) eq ref($bv)) && vecequal($av, $bv);
      return 0;
    }
    # About 7x faster if we skip the validates.
    # _validate_integer($av);
    # _validate_integer($bv);
    return 0 unless "$av" eq "$bv";
  }
  1;
}

sub vecmex {
  my $items = scalar(@_);
  my @seen;
  for (@_) {
    $seen[$_] = 0 if $_ < $items;
  }
  for (0 .. $items-1) {
    return $_ unless defined $seen[$_];
  }
  return $items;
}

sub vecpmex {
  my $items = scalar(@_);
  my @seen;
  for (@_) {
    $seen[$_] = 0 if $_ <= $items;
  }
  for (1 .. $items) {
    return $_ unless defined $seen[$_];
  }
  return $items+1;
}

sub sumdigits {
  my($n,$base) = @_;
  my $sum = 0;
  $base =  2 if !defined $base && $n =~ s/^0b//;
  $base = 16 if !defined $base && $n =~ s/^0x//;
  if (!defined $base || $base == 10) {
    $n =~ tr/0123456789//cd;
    $sum += $_ for (split(//,$n));
  } else {
    croak "sumdigits: invalid base $base" if $base < 2;
    my $cmap = substr("0123456789abcdefghijklmnopqrstuvwxyz",0,$base);
    for my $c (split(//,lc($n))) {
      my $p = index($cmap,$c);
      $sum += $p if $p > 0;
    }
  }
  $sum;
}

sub is_happy {
  my($n, $base, $k) = @_;
  validate_integer_nonneg($n);

  my $h = 1;

  if (!defined $base && !defined $k) {   # default base 10 exponent 2
    while ($n > 1 && $n != 4) {
      my $sum = 0;
      $sum += $_*$_ for (split(//,$n));
      $n = $sum;
      $h++;
    }
    return ($n == 1) ? $h : 0;
  }

  if (defined $base) {
    validate_integer_nonneg($base);
    croak "is_happy: invalid base $base" if $base < 2 || $base > 36;
  } else {
    $base = 10;
  }
  if (defined $k) {
    validate_integer_nonneg($k);
    croak "is_happy: invalid exponent $k" if $k > 10;
  } else {
    $k = 2;
  }

  my %seen;
  while ($n > 1 && !exists $seen{$n}) {
    $seen{$n} = undef;
    if ($base == 10) {
      my $sum = 0;
      $sum += $_ ** $k for (split(//,$n));
      $n = $sum;
    } else {
      my @d;
      while ($n >= 1) {
        my $rem = $n % $base;
        push @d, ($k <= 6) ? int($rem ** $k) : Mpowint($rem,$k);
        #push @d, Mpowint($rem,$k);
        $n = ($n-$rem)/$base;    # Always an exact division
      }
      $n = Mvecsum(@d);
    }
    $h++;
  }
  return ($n == 1) ? $h : 0;
}



# Tonelli-Shanks
sub _sqrtmod_prime {
  my($a, $p) = @_;
  my($x, $q, $e, $t, $z, $r, $m, $b);
  my $Q = Msubint($p,1);

  if (($p % 4) == 3) {
    $r = Mpowmod($a, Mrshiftint(Maddint($p,1),2), $p);
    return undef unless Mmulmod($r,$r,$p) == $a;
    return $r;
  }
  if (($p % 8) == 5) {
    $m = Maddmod($a,$a,$p);
    $t = Mpowmod($m, Mrshiftint(Msubint($p,5),3), $p);
    $z = Mmulmod($m, Mmulmod($t,$t,$p), $p);
    $r = Mmulmod($t, Mmulmod($a, Msubmod($z,1,$p), $p), $p);
    return undef unless Mmulmod($r,$r,$p) == $a;
    return $r;
  }

  # Verify Euler's criterion for odd p
  return undef if $p != 2 && Mpowmod($a, Mrshiftint($Q,1), $p) != 1;

  # Cohen Algorithm 1.5.1.  Tonelli-Shanks.
  $e = Mvaluation($Q, 2);
  $q = Mdivint($Q, Mpowint(2,$e));
  $t = 3;
  while (Mkronecker($t,$p) != -1) {
    $t += 2;
    return undef if $t == 201 && !Mis_prime($p);
  }
  $z = Mpowmod($t, $q, $p);
  $b = Mpowmod($a, $q, $p);
  $r = $e;
  $q = ($q+1) >> 1;
  $x = Mpowmod($a, $q, $p);
  while ($b != 1) {
    $t = $b;
    for ($m = 0;  $m < $r && $t != 1;  $m++) {
      $t = Mmulmod($t, $t, $p);
    }
    $t = Mpowmod($z, Mlshiftint(1, $r-$m-1), $p);
    $x = Mmulmod($x, $t, $p);
    $z = Mmulmod($t, $t, $p);
    $b = Mmulmod($b, $z, $p);
    $r = $m;
  }
  # Expected to always be true.
  return undef unless Mmulmod($x,$x,$p) == $a;
  return $x;
}

sub _sqrtmod_prime_power {
  my($a,$p,$e) = @_;
  my($r,$s);

  if ($e == 1) {
    $a %= $p if $a >= $p;
    return $a if $p == 2 || $a == 0;
    $r = _sqrtmod_prime($a,$p);
    return (defined $r && (Mmulmod($r,$r,$p) == $a) ? $r : undef);
  }

  my $n = Mpowint($p,$e);
  my $pk = Mmulint($p,$p);

  return 0 if ($a % $n) == 0;

  if (($a % $pk) == 0) {
    my $apk = Mdivint($a, $pk);
    $s = _sqrtmod_prime_power($apk, $p, $e-2);
    return undef unless defined $s;
    return Mmulint($s,$p);
  }

  return undef if ($a % $p) == 0;

  my $ered = ($p > 2 || $e < 5)  ?  ($e+1) >> 1  :  ($e+3) >> 1;
  $s = _sqrtmod_prime_power($a,$p,$ered);
  return undef unless defined $s;

  my $np  = ($p == 2)  ?  Mmulint($n,$p)  :  $n;
  my $t1  = Msubmod($a, Mmulmod($s,$s,$np), $np);
  my $t2  = Maddmod($s, $s, $np);
  my $gcd = Mgcd($t1, $t2);
  $r = Maddmod($s, Mdivmod(Mdivint($t1,$gcd),Mdivint($t2,$gcd),$n), $n);
  return ((Mmulmod($r,$r,$n) == ($a % $n)) ? $r : undef);
}

sub _sqrtmod_composite {
  my($a,$n) = @_;

  return undef if $n <= 0;
  $a %= $n if $a >= $n;
  return $a if $n <= 2 || $a <= 1;
  return Msqrtint($a) if _is_perfect_square($a);

  my $N = 1;
  my $r = 0;
  foreach my $F (Mfactor_exp($n)) {
    my($f,$e) = @$F;
    my $fe = Mpowint($f, $e);
    my $s = _sqrtmod_prime_power($a, $f, $e);
    return undef unless defined $s;
    my $inv = Minvmod($N, $fe);
    my $t = Mmulmod($inv, Msubmod($s % $fe, $r % $fe, $fe), $fe);
    $r = Mmuladdmod($N, $t, $r, $n);
    $N = Mmulint($N, $fe);
  }
  #croak "Bad _sqrtmod_composite root $a,$n" unless Mmulmod($r,$r,$n) == $a;
  $r;
}

sub sqrtmod {
  my($a,$n) = @_;
  validate_integer($a);
  validate_integer_abs($n);
  return (undef,0)[$n] if $n <= 1;
  #return Mmodint(Msqrtint($a),$n) if _is_perfect_square($a);
  $a = Mmodint($a,$n);

  my $r = _sqrtmod_composite($a,$n);
  if (defined $r) {
    $r = $n-$r if $n-$r < $r;
    #croak "Bad _sqrtmod_composite root $a,$n" unless Mmulmod($r,$r,$n) == $a;
  }
  return $r;
}




# helper function for allsqrtmod() - return list of all square roots of
# a (mod p^k), assuming a integer, p prime, k positive integer.
sub _allsqrtmodpk {
  my($a,$p,$k) = @_;
  my $pk = Mpowint($p,$k);
  unless ($a % $p) {
    unless ($a % ($pk)) {
      # if p^k divides a, we need the square roots of zero, satisfied by
      # ip^j with 0 <= i < p^{floor(k/2)}, j = p^{ceil(k/2)}
      my $low = Mpowint($p,$k >> 1);
      my $high = ($k % 2)  ?  Mmulint($low, $p)  :  $low;
      return map Mmulint($high, $_), 0 .. $low - 1;
    }
    # p divides a, p^2 does not
    my $a2 = Mdivint($a,$p);
    return () if $a2 % $p;
    my $pj = Mdivint($pk, $p);
    return map {
      my $qp = Mmulint($_,$p);
      map Maddint($qp,Mmulint($_,$pj)), 0 .. $p - 1;
    } _allsqrtmodpk(Mdivint($a2,$p), $p, $k - 2);
  }
  my $q = _sqrtmod_prime_power($a,$p,$k);
  return () unless defined $q;
  return ($q, $pk - $q) if $p != 2;
  return ($q) if $k == 1;
  return ($q, $pk - $q) if $k == 2;
  my $pj = Mdivint($pk,$p);
  my $q2 = ($q * ($pj - 1)) % $pk;
  return ($q, $pk - $q, $q2, $pk - $q2);
}

# helper function for allsqrtmod() - return list of all square roots of
# a (mod p^k), assuming a integer, n positive integer > 1, f arrayref
# of [ p, k ] pairs representing factorization of n. Destroys f.
sub _allsqrtmodfact {
  my($a,$n,$f) = @_;
  my($p,$k) = @{ shift @$f };
  my @q = _allsqrtmodpk($a, $p, $k);
  return @q unless @$f;
  my $pk = Mpowint($p, $k);
  my $n2 = Mdivint($n, $pk);
  return map {
    my $q2 = $_;
    map Mchinese([ $q2, $n2 ], [ $_, $pk ]), @q;
  } _allsqrtmodfact($a, $n2, $f);
}

sub allsqrtmod {
  my($A,$n) = @_;
  validate_integer($A);
  validate_integer_abs($n);
  return $n ? (0) : () if $n <= 1;
  $A = Mmodint($A,$n);
  Mvecsort(  _allsqrtmodfact($A, $n, [ Mfactor_exp($n) ])  );
}


###############################################################################
#       Tonelli-Shanks kth roots
###############################################################################

# Algorithm 3.3, step 2 "Find generator"
sub _find_ts_generator {
  my ($a, $k, $p) = @_;
  # Assume:  k > 2,  1 < a < p,  p > 2,  k prime,  p prime

  my($e,$r) = (0, $p-1);
  while (!($r % $k)) {
    $e++;
    $r /= $k;
  }
  my $ke1 = Mpowint($k, $e-1);
  my($x,$m,$y) = (2,1);
  while ($m == 1) {
    $y = Mpowmod($x, $r, $p);
    $m = Mpowmod($y, $ke1, $p) if $y != 1;
    croak "bad T-S input" if $x >= $p;
    $x++;
  }
  ($y, $m);
}

sub _ts_rootmod {
  my($a, $k, $p, $y, $m) = @_;

  my($e,$r) = (0, $p-1);
  while (!($r % $k)) {
    $e++;
    $r /= $k;
  }
  # p-1 = r * k^e
  my $x = Mpowmod($a, Minvmod($k % $r, $r), $p);
  my $A = ($a == 0) ? 0 : Mmulmod(Mpowmod($x,$k,$p), Minvmod($a,$p), $p);

  ($y,$m) = _find_ts_generator($a,$k,$p) if $y == 0 && $A != 1;

  while ($A != 1) {
    my ($l,$T,$z) = (1,$A);
    while ($T != 1) {
      return 0 if $l >= $e;
      $z = $T;
      $T = Mpowmod($T, $k, $p);
      $l++;
    }
    # We want a znlog that takes gorder as well (k=znorder(m,p))
    my $kz = _negmod(znlog($z, $m, $p), $k);
    $m = Mpowmod($m, $kz, $p);
    $T = Mpowmod($y, Mmulint($kz,Mpowint($k,$e-$l)), $p);
    # In the loop we always end with l < e, so e always gets smaller
    $e = $l-1;
    $x = Mmulmod($x, $T, $p);
    $y = Mpowmod($T, $k, $p);
    return 0 if $y <= 1;  # In theory this will never be hit.
    $A = Mmulmod($A, $y, $p);
  }
  $x;
}

sub _compute_generator {
  my($l, $e, $r, $p) = @_;
  my($m, $lem1, $y) = (1, Mpowint($l, $e-1));
  for (my $x = 2; $m == 1; $x++) {
    $y = Mpowmod($x, $r, $p);
    next if $y == 1;
    $m = Mpowmod($y, $lem1, $p);
  }
  $y;
}

sub _rootmod_prime_splitk {
  my($a, $k, $p, $refzeta) = @_;

  $$refzeta = 1 if defined $refzeta;
  $a = Mmodint($a, $p) if $a >= $p;
  return $a if $a == 0 || ($a == 1 && !defined $refzeta);
  my $p1 = Msubint($p,1);

  if ($k == 2) {
    my $r = _sqrtmod_prime($a,$p);
    $$refzeta = (defined $r) ? $p1 : 0     if defined $refzeta;
    return $r;
  }

  # See Algorithm 2.1 of van de Woestijne (2006), or Lindhurst (1997).
  # The latter's proposition 7 generalizes to composite p.

  my $g = Mgcd($k, $p1);
  my $r = $a;

  if ($g != 1) {
    foreach my $fac (Mfactor_exp($g)) {
      my($F,$E) = @$fac;
      last if $r == 0;
      if (defined $refzeta) {
        my $V   = Mvaluation($p1, $F);
        my $REM = Mdivint($p1, Mpowint($F,$V));
        my $Y   = _compute_generator($F, $V, $REM, $p);
        $$refzeta = Mmulmod($$refzeta, Mpowmod($Y, Mpowint($F, $V-$E), $p), $p);
      }
      my ($y,$m) = _find_ts_generator($r, $F, $p);
      while ($E-- > 0) {
        $r = _ts_rootmod($r, $F, $p,  $y, $m);
      }
    }
  }
  if ($g != $k) {
    my($kg, $pg) = (Mdivint($k,$g), Mdivint($p1,$g));
    $r = Mpowmod($r, Minvmod($kg % $pg, $pg), $p);
  }
  return $r if Mpowmod($r, $k, $p) == $a;
  $$refzeta = 0 if defined $refzeta;
  undef;
}

sub _rootmod_composite1 {
  my($a,$k,$n) = @_;
  my $r;

  croak "_rootmod_composite1 bad parameters" if $a < 1 || $k < 2 || $n < 2;

  if (Mis_power($a, $k, \$r)) {
    return $r;
  }

  if (Mis_prime($n)) {
    return _rootmod_prime_splitk($a,$k,$n,undef);
  }

  # We should do this iteratively using cprod
  my @rootmap;
  foreach my $fac (Mfactor_exp($n)) {
    my($F,$E) = @$fac;
    my $FE = Mpowint($F,$E);
    my $A = $a % $FE;
    if ($E == 1) {
      $r = _rootmod_prime_splitk($A,$k,$F,undef)
    } else {
      # TODO: Fix this.  We should do this directly.
      $r = (allrootmod($A, $k, $FE))[0];
    }
    return undef unless defined $r && Mpowmod($r, $k, $FE) == $A;
    push @rootmap, [ $r, $FE ];
  }
  $r = Mchinese(@rootmap) if @rootmap > 1;

  #return (defined $r && Mpowmod($r, $k, $n) == ($a % $n))  ?  $r  :  undef;
  croak "Bad _rootmod_composite1 root $a,$k,$n" unless defined $r && Mpowmod($r,$k,$n) == ($a % $n);
  $r;
}

###############################################################################
#       Tonelli-Shanks kth roots  alternate version
###############################################################################

sub _ts_prime {
  my($a, $k, $p, $refzeta) = @_;

  my($e,$r) = (0, $p-1);
  while (!($r % $k)) {
    $e++;
    $r /= $k;
  }
  my $ke = Mdivint($p-1, $r);

  my $x = Mpowmod($a, Minvmod($k % $r, $r), $p);
  my $B = Mmulmod(Mpowmod($x, $k, $p), Minvmod($a, $p), $p);

  my($T,$y,$t,$A) = (2,1);
  while ($y == 1) {
    $t = Mpowmod($T, $r, $p);
    $y = Mpowmod($t, Mdivint($ke,$k), $p);
    $T++;
  }
  while ($ke != $k) {
    $ke = Mdivint($ke, $k);
    $T = $t;
    $t = Mpowmod($t, $k, $p);
    $A = Mpowmod($B, Mdivint($ke,$k), $p);
    while ($A != 1) {
      $x = Mmulmod($x, $T, $p);
      $B = Mmulmod($B, $t, $p);
      $A = Mmulmod($A, $y, $p);
    }
  }
  $$refzeta = $t if defined $refzeta;
  $x;
}

sub _rootmod_prime {
  my($a, $k, $p) = @_;

  # p must be a prime, k must be a prime.  Otherwise UNDEFINED.
  $a %= $p if $a >= $p;

  return $a if $p == 2 || $a == 0;
  return _sqrtmod_prime($a, $p) if $k == 2;

  # If co-prime, there is exactly one root.
  my $g = Mgcd($k, $p-1);
  return Mpowmod($a, Minvmod($k % ($p-1), $p-1), $p)  if $g == 1;
  # Check generalized Euler's criterion
  return undef if Mpowmod($a, Mdivint($p-1, $g), $p) != 1;

  _ts_prime($a, $k, $p);
}

sub _rootmod_prime_power {
  my($a,$k,$p,$e) = @_;        # prime k, prime p

  return _sqrtmod_prime_power($a, $p, $e) if $k == 2;
  return _rootmod_prime($a, $k, $p)       if $e == 1;

  my $n = Mpowint($p,$e);
  my $pk = Mpowint($p,$k);

  return 0 if ($a % $n) == 0;

  if (($a % $pk) == 0) {
    my $apk = Mdivint($a, $pk);
    my $s = _rootmod_prime_power($apk, $k, $p, $e-$k);
    return (defined $s)  ?  Mmulint($s,$p)  :  undef;
  }

  return undef if ($a % $p) == 0;

  my $ered = ($p > 2 || $e < 5)  ?  ($e+1) >> 1  :  ($e+3) >> 1;
  my $s = _rootmod_prime_power($a, $k, $p, $ered);
  return undef if !defined $s;

  my $np  = ($p == $k)  ?  Mmulint($n,$p)  :  $n;
  my $t = Mpowmod($s, $k-1, $np);
  my $t1  = Msubmod($a, Mmulmod($t,$s,$np), $np);
  my $t2  = Mmulmod($k, $t, $np);
  my $gcd = Mgcd($t1, $t2);
  my $r   = Maddmod($s,Mdivmod(Mdivint($t1,$gcd),Mdivint($t2,$gcd),$n),$n);
  return ((Mpowmod($r,$k,$n) == ($a % $n)) ? $r : undef);
}

sub _rootmod_kprime {
  my($a,$k,$n,@nf) = @_;       # k prime, n factored into f^e,f^e,...

  my($N,$r) = (1,0);
  foreach my $F (@nf) {
    my($f,$e) = @$F;
    my $fe = Mpowint($f, $e);
    my $s = _rootmod_prime_power($a, $k, $f, $e);
    return undef unless defined $s;
    my $inv = Minvmod($N, $fe);
    my $t = Mmulmod($inv, Msubmod($s % $fe, $r % $fe, $fe), $fe);
    $r = Mmuladdmod($N, $t, $r, $n);
    $N = Mmulint($N, $fe);
  }
  $r;
}

sub _rootmod_composite2 {
  my($a,$k,$n) = @_;

  croak "_rootmod_composite2 bad parameters" if $a < 1 || $k < 2 || $n < 2;

  my @nf = Mfactor_exp($n);

  return _rootmod_kprime($a, $k, $n, @nf) if Mis_prime($k);

  my $r = $a;
  foreach my $kf (Mfactor($k)) {
    $r = _rootmod_kprime($r, $kf, $n, @nf);
    if (!defined $r) {
      # Choose one.  The former is faster but makes more intertwined code.
      return _rootmod_composite1($a,$k,$n);
      #return (allrootmod($a,$k,$n))[0];
    }
  }
  croak "Bad _rootmod_composite2 root $a,$k,$n" unless defined $r && Mpowmod($r,$k,$n) == ($a % $n);
  $r;
}


###############################################################################
#       Modular k-th root
###############################################################################

sub rootmod {
  my($a,$k,$n) = @_;
  validate_integer($a);
  validate_integer($k);
  validate_integer_abs($n);
  return (undef,0)[$n] if $n <= 1;
  $a = Mmodint($a,$n);

  # Be careful with zeros, as we can't divide or invert them.
  if ($a == 0) {
    return ($k <= 0) ? undef : 0;
  }
  if ($k < 0) {
    $a = Minvmod($a, $n);
    return undef unless defined $a && $a > 0;
    $k = -$k;
  }
  return undef if $k == 0 && $a != 1;
  return 1 if $k == 0 || $a == 1;
  return $a if $k == 1;

  # Choose either one based on performance.
  my $r = _rootmod_composite1($a, $k, $n);
  #my $r = _rootmod_composite2($a, $k, $n);
  $r = $n-$r if defined $r && $k == 2 && ($n-$r) < $r; # Select smallest root
  $r;
}

###############################################################################
#       All modular k-th roots
###############################################################################

sub _allrootmod_cprod {
  my($aroots1, $p1, $aroots2, $p2) = @_;
  my($t, $n, $inv);

  $n = mulint($p1, $p2);
  $inv = Minvmod($p1, $p2);
  croak("CRT has undefined inverse") unless defined $inv;

  my @roots;
  for my $q1 (@$aroots1) {
    for my $q2 (@$aroots2) {
      $t = Mmulmod($inv, Msubmod($q2, $q1, $p2), $p2);
      $t = Mmuladdmod($p1, $t, $q1, $n);
      push @roots, $t;
    }
  }
  return @roots;
}

sub _allrootmod_prime {
  my($a,$k,$p) = @_;        # prime k, prime p
  $a %= $p if $a >= $p;     #$a = Mmodint($a,$p) if $a >= $p;

  return ($a) if $p == 2 || $a == 0;

  # If co-prime, there is exactly one root.
  my $g = Mgcd($k, $p-1);
  if ($g == 1) {
    my $r = Mpowmod($a, Minvmod($k % ($p-1), $p-1), $p);
    return ($r);
  }

  # Check generalized Euler's criterion
  return () if Mpowmod($a, Mdivint($p-1, $g), $p) != 1;

  # Special case for p=3 for performance
  return (1,2) if $p == 3;

  # A trivial brute force search:
  # return grep { Mpowmod($_,$k,$p) == $a } 0 .. $p-1;

  # Call one of the general TS solvers that also allow us to get all the roots.
  my $z;
  #my $r = _rootmod_prime_splitk($a, $k, $p, \$z);
  my $r = _ts_prime($a, $k, $p, \$z);
  croak "allrootmod: failed to find root" if $z==0 || Mpowmod($r,$k,$p) != $a;
  my @roots = ($r);
  my $r2 = Mmulmod($r,$z,$p);
  while ($r2 != $r && @roots < $k) {
    push @roots, $r2;
    $r2 = Mmulmod($r2, $z, $p);
  }
  croak "allrootmod: excess roots found" if $r2 != $r;
  return @roots;
}

sub _allrootmod_prime_power {
  my($a,$k,$p,$e) = @_;        # prime k, prime p

  return _allrootmod_prime($a, $k, $p) if $e == 1;

  my $n  = ($e<=13 && $p<=13)||($e<=5 && $p<=1000) ?int($p**$e):Mpowint($p,$e);
  my $pk = ($k<=13 && $p<=13)||($k<=5 && $p<=1000) ?int($p**$k):Mpowint($p,$k);
  my @roots;

  if (($a % $n) == 0) {
    my $t = Mdivint($e-1, $k) + 1;
    my $nt = Mpowint($p, $t);
    my $nr = Mpowint($p, $e-$t);
    @roots = map { Mmulmod($_, $nt, $n) } 0 .. $nr-1;
    return @roots;
  }

  if (($a % $pk) == 0) {
    my $apk = Mdivint($a, $pk);
    my $pe1 = Mpowint($p, $k-1);
    my $pek = Mpowint($p, $e-$k+1);
    my @roots2 = _allrootmod_prime_power($apk, $k, $p, $e-$k);
    for my $r (@roots2) {
      my $rp = Mmulmod($r, $p, $n);
      push @roots, Mmuladdmod($_, $pek, $rp, $n) for 0 .. $pe1-1;
    }
    return @roots;
  }

  return () if ($a % $p) == 0;

  my $np  = Mmulint($n,$p);
  my $ered = ($p > 2 || $e < 5)  ?  ($e+1) >> 1  :  ($e+3) >> 1;
  my @roots2 = _allrootmod_prime_power($a, $k, $p, $ered);

  if ($k != $p) {
    for my $s (@roots2) {
      my $t = Mpowmod($s, $k-1, $n);
      my $t1  = Msubmod($a, Mmulmod($t,$s,$n), $n);
      my $t2  = Mmulmod($k, $t, $n);
      my $gcd = Mgcd($t1, $t2);
      my $r   = Maddmod($s,Mdivmod(Mdivint($t1,$gcd),Mdivint($t2,$gcd),$n),$n);
      push @roots, $r;
    }
  } else {
    my @rootst;
    for my $s (@roots2) {
      my $t  = Mpowmod($s, $k-1, $np);
      my $t1  = Msubmod($a, Mmulmod($t,$s,$np), $np);
      my $t2  = Mmulmod($k, $t, $np);
      my $gcd = Mgcd($t1, $t2);
      my $r   = Maddmod($s,Mdivmod(Mdivint($t1,$gcd), Mdivint($t2,$gcd),$n),$n);
      push @rootst, $r if Mpowmod($r, $k, $n) == ($a % $n);
    }
    my $ndivp = Mdivint($n,$p);
    my $rset = [];
    for my $r (@rootst) {
      Msetinsert($rset, Mmulmod($r, Mmuladdmod($_, $ndivp, 1, $n), $n))
        for 0 .. $k-1;
    }
    @roots = @$rset;
  }
  return @roots;
}

sub _allrootmod_kprime {
  my($a,$k,$n,@nf) = @_;       # k prime, n factored into f^e,f^e,...

  return _allsqrtmodfact($a, $n, \@nf) if $k == 2;

  my $N = 1;
  my @roots;
  foreach my $F (@nf) {
    my($f,$e) = @$F;
    my @roots2 = ($e==1) ? _allrootmod_prime($a, $k, $f)
                         : _allrootmod_prime_power($a, $k, $f, $e);
    return () unless @roots2;
    my $fe = ($e <= 13 && $f <= 13) ? int($f**$e) : Mpowint($f, $e);
    if (scalar(@roots) == 0) {
      @roots = @roots2;
    } else {
      @roots = _allrootmod_cprod(\@roots, $N, \@roots2, $fe);
    }
    $N = Mmulint($N, $fe);
  }

  return @roots;
}

sub allrootmod {
  my($A,$k,$n) = @_;
  validate_integer($A);
  validate_integer($k);
  validate_integer_abs($n);

  return () if $n == 0;
  $A = Mmodint($A,$n);

  return () if $k <= 0 && $A == 0;

  if ($k < 0) {
    $A = Minvmod($A, $n);
    return () unless defined $A && $A > 0;
    $k = -$k;
  }

  # TODO: For testing
  #my @roots = sort { $a <=> $b }
  #            grep { Mpowmod($_,$k,$n) == $A } 0 .. $n-1;
  #return @roots;

  return ($A) if $n <= 2 || $k == 1;
  return ($A == 1) ? (0..$n-1) : ()  if $k == 0;

  my @roots;
  my @nf = Mfactor_exp($n);

  if (Mis_prime($k)) {
    @roots = _allrootmod_kprime($A, $k, $n, @nf);
  } else {
    @roots = ($A);
    for my $primek (Mfactor($k)) {
      my @rootsnew = ();
      for my $r (@roots) {
        push @rootsnew, _allrootmod_kprime($r, $primek, $n, @nf);
      }
      @roots = @rootsnew;
    }
  }

  Mvecsort(@roots);
}

################################################################################
################################################################################

sub _modabsint {
  my($a, $n) = @_;
  if ($n <= 1) {
    if ($n < 0) { $n = tobigint($n) if $n <= INTMIN && !ref($n);  $n = -$n; }
    return (undef,0)[$n] if $n <= 1;
  }
  if ($n < INTMAX && $a < INTMAX && $a > INTMIN) {
    $a = $n - ((-$a) % $n) if $a < 0;
    $a %= $n if $a >= $n;
  } else {
    $a = tobigint($a) % $n;
    $a = _bigint_to_int($a) if $a <= INTMAX;
  }
  $a;
}

sub addmod {
  my($a, $b, $n) = @_;
  if ($n <= 1) {
    if ($n < 0) { $n = tobigint($n) if $n <= INTMIN && !ref($n);  $n = -$n; }
    return (undef,0)[$n] if $n <= 1;
  }
  if ($n <= INTMAX && $a <= INTMAX && $b <= INTMAX && $a >= INTMIN && $b >= INTMIN) {
    $a = $n - ((-$a) % $n) if $a < 0;
    $b = $n - ((-$b) % $n) if $b < 0;
    $a %= $n if $a >= $n;
    $b %= $n if $b >= $n;
    return  $n-$a > $b  ?  $a+$b
          : $a > $b     ?  ($a-$n)+$b
          :                ($b-$n)+$a;
  }
  # Impl 1.  Make $a a bigint and let things promote.  Fastest.
  my $r = (tobigint($a) + $b) % $n;
  return $r <= INTMAX ? _bigint_to_int($r) : $r;

  # Impl 2.  Use Maddint but mod with a $n as a bigint.
  #my $r = Maddint($a,$b) % tobigint($n);
  #return $r <= INTMAX ? _bigint_to_int($r) : $r;

  # Impl 3.  Prefered, but slowest.  Probably fine when we use amagic in XS.
  #Mmodint(Maddint($a,$b),$n);
}
sub submod {
  my($a, $b, $n) = @_;
  if ($n <= 1) {
    if ($n < 0) { $n = tobigint($n) if $n <= INTMIN && !ref($n);  $n = -$n; }
    return (undef,0)[$n] if $n <= 1;
  }
  if ($n <= INTMAX && $a <= INTMAX && $b <= INTMAX && $a >= INTMIN && $b >= INTMIN) {
    $a = $n - ((-$a) % $n) if $a < 0;
    $b = $n - ((-$b) % $n) if $b < 0;
    $a %= $n if $a >= $n;
    $b %= $n if $b >= $n;
    $b = $n-$b;  # negate b then we add as above
    return  $n-$a > $b  ?  $a+$b
          : $a > $b     ?  ($a-$n)+$b
          :                ($b-$n)+$a;
  }
  my $r = (tobigint($a) - $b) % $n;
  return $r <= INTMAX ? _bigint_to_int($r) : $r;
}

sub mulmod {
  my($a, $b, $n) = @_;
  #if ($n <= 1) { # ABS(n) and handle mod 0 | mod 1.
  #  if ($n < 0) { $n = tobigint($n) if $n <= INTMIN && !ref($n);  $n = -$n; }
  #  return (undef,0)[$n] if $n <= 1;
  #}
  if ($n <= 1) {
    return (undef,0)[$n] if $n >= 0;
    $n = tobigint($n) if $n <= INTMIN && !ref($n);
    $n = -$n;
    return 0 if $n == 1;
  }

  # If n is a native int, we can reduce a and b then do everything native
  if ($n < INTMAX) {
    if ($a >= INTMAX || $a < 0 || $b >= INTMAX || $b < 0) {
      $a = _bigint_to_int(tobigint($a) % $n) if $a >= INTMAX || $a < 0;
      $b = _bigint_to_int(tobigint($b) % $n) if $b >= INTMAX || $b < 0;
    }
    return _mulmod($a,$b,$n);
  }

  # Try GMP
  return reftyped($_[0], Math::Prime::Util::GMP::mulmod($a,$b,$n))
    if $Math::Prime::Util::_GMPfunc{"mulmod"};

  my $refn = ref($n);
  if (!$refn) {
    $n = tobigint($n);
    $refn = ref($n);
  }
  $a = $refn->new("$a") unless ref($a) eq $refn;
  $b = $refn->new("$b") unless ref($b) eq $refn;
  my $r = ($a * $b) % $n;
  return $r <= INTMAX ? _bigint_to_int($r) : $r;
}

sub _bi_powmod {
  my($a, $b, $n) = @_;
  croak "_bi_powmod must have positive exponent" if $b < 0;
  croak "_bi_powmod must have n > 1" if $n <= 1;

  my $refn = ref($n);
  if (!$refn) {
    $n = tobigint($n);
    $refn = ref($n);
  }
  $b = $refn->new($b) unless ref($b) eq $refn;

  my $r = $refn->new($a);

  if ($refn eq 'Math::GMPz') {
    Math::GMPz::Rmpz_powm($r, $r, $b, $n);
  } elsif ($refn eq 'Math::GMP') {
    $r = $r->powm_gmp($b,$n);
  } elsif ($refn eq 'Math::BigInt') {
    $r->bmodpow($b,$n);
  } else {
    $r->bmodpow("$b","$n");
  }
  $r;
}

sub powmod {
  my($a, $b, $n) = @_;
  if ($n <= 1) {
    if ($n < 0) { $n = tobigint($n) if $n <= INTMIN && !ref($n);  $n = -$n; }
    return (undef,0)[$n] if $n <= 1;
  }
  return ($b > 0) ? 0 : 1  if $a == 0;

  if ($Math::Prime::Util::_GMPfunc{"powmod"}) {
    my $r = Math::Prime::Util::GMP::powmod($a,$b,$n);
    return (defined $r) ? reftyped($_[0], $r) : undef;
  }

  # If the exponent is negative: a=1/a ; b=-b
  if ($b < 0) {
    $a = Minvmod($a,$n);
    return undef unless defined $a;
    $b = tobigint($b) if $b <= INTMIN && !ref($b);
    $b = -$b;
  }

  if ($b <= 8) {
    return 1 if $b == 0;
    return _modabsint($a,$n) if $b == 1;
    return Mmulmod($a,$a,$n) if $b == 2;
    # For exponents 3-8, this can be 20x faster for native n
    if (!ref($n) && $a <= 31622776 && $a >= -31622776) {
      my $a2 = int($a*$a);
      return Mmulmod($a2,$a,$n) if $b == 3;
      my $a4 = Mmulmod($a2,$a2,$n);
      return $a4 if $b == 4;
      return Mmulmod($a4,$a,$n) if $b == 5;
      return Mmulmod($a4,$a2,$n) if $b == 6;
      return Mmulmod($a4,Mmulmod($a2,$a,$n),$n) if $b == 7;
      return Mmulmod($a4,$a4,$n) if $b == 8;
    }
  }

  my $r = _bi_powmod($a,$b,$n);
  return $r <= INTMAX ? _bigint_to_int($r) : $r;
}

sub muladdmod {
  my($a, $b, $c, $n) = @_;
  if ($n <= 1) {
    $n = Mnegint($n) if $n < 0;
    return (undef,0)[$n] if $n <= 1;
  }

  if (!ref($n) && $n <= INTMAX
               && $a <= INTMAX && $b <= INTMAX && $c <= INTMAX
               && $a >= INTMIN && $b >= INTMIN && $c >= INTMIN) {
    $a = $n - ((-$a) % $n) if $a < 0;
    $b = $n - ((-$b) % $n) if $b < 0;
    $c = $n - ((-$c) % $n) if $c < 0;
    #$c %= $n if $c >= $n;   # For mulsubmod
    return _addmod(_mulmod($a,$b,$n),$c,$n);
  }
  return reftyped($_[0], Math::Prime::Util::GMP::muladdmod($a,$b,$c,$n))
    if $Math::Prime::Util::_GMPfunc{"muladdmod"};

  $n = tobigint($n) unless ref($n);
  $a = tobigint($a) unless ref($a) || ref($b);
  $c = tobigint($c) unless ref($c);
  my $r = (($a * $b) + $c) % $n;
  return $r <= INTMAX ? _bigint_to_int($r) : $r;
}
sub mulsubmod {
  my($a, $b, $c, $n) = @_;
  if ($n <= 1) {
    $n = Mnegint($n) if $n < 0;
    return (undef,0)[$n] if $n <= 1;
  }

  if (!ref($n) && $n <= INTMAX
               && $a <= INTMAX && $b <= INTMAX && $c <= INTMAX
               && $a >= INTMIN && $b >= INTMIN && $c >= INTMIN) {
    $a = $n - ((-$a) % $n) if $a < 0;
    $b = $n - ((-$b) % $n) if $b < 0;
    $c = $n - ((-$c) % $n) if $c < 0;
    $c = ($c < $n)  ?  $n-$c  :  $n-($c % $n);  # $c = -$c (mod n)
    return _addmod(_mulmod($a,$b,$n),$c,$n);
  }
  return reftyped($_[0], Math::Prime::Util::GMP::mulsubmod($a,$b,$c,$n))
    if $Math::Prime::Util::_GMPfunc{"mulsubmod"};

  # return Msubmod(Mmulmod($a,$b,$n),$c,$n);

  $n = tobigint($n) unless ref($n);
  $a = tobigint($a) unless ref($a) || ref($b);
  $c = tobigint($c) unless ref($c);
  my $r = (($a * $b) - $c) % $n;
  return $r <= INTMAX ? _bigint_to_int($r) : $r;
}

sub invmod {
  my($a,$n) = @_;
  if ($n <= 1) {
    if ($n < 0) { $n = tobigint($n) if $n <= INTMIN && !ref($n);  $n = -$n; }
    return (undef,0)[$n] if $n <= 1;
  }
  return if $a == 0;

  if ($n < INTMAX) {  # Fast all native math
    my($t,$nt,$r,$nr) = (0, 1, $n, _modabsint($a,$n));
    while ($nr != 0) {
      # Use mod before divide to force correct behavior with high bit set
      my $quot = int( ($r-($r % $nr))/$nr );
      ($nt,$t) = ($t-$quot*$nt,$nt);
      ($nr,$r) = ($r-$quot*$nr,$nr);
    }
    return $r > 1  ?  undef
        :  $t < 0  ?  $t+$n
                   :  $t;
  }

  $n = tobigint($n);
  $a = tobigint($a) % $n;
  my $refn = ref($n);
  my $I;

  if ($refn eq 'Math::BigInt') {
    $I = $a->copy->bmodinv($n);
    $I = undef if defined $I && !$I->is_int();
  } elsif ($refn eq 'Math::GMPz') {
    $I = Math::GMPz->new();
    Math::GMPz::Rmpz_invert($I, $a, $n);
    $I = undef if defined $I && $I == 0;
  } elsif ($refn eq 'Math::GMP') {
    $I = $a->gmp_copy->bmodinv($n);
    $I = undef if defined $I && $I == 0;
  } else {
    $I = Math::BigInt->new("$a")->bmodinv("$n");
    $I = undef if defined $I && !$I->is_int();
    $I = tobigint("$I") if defined $I;
  }
  $I = _bigint_to_int($I) if defined $I && $I <= INTMAX;
  return $I;
}

sub divmod {
  my($a, $b, $n) = @_;
  if ($n <= 1) {
    if ($n < 0) { $n = tobigint($n) if $n <= INTMIN && !ref($n);  $n = -$n; }
    return (undef,0)[$n] if $n <= 1;
  }

  my $invb = Minvmod($b,$n);
  return undef unless defined $invb;
  return Mmulmod($a,$invb,$n);
}

sub negmod {
  my($a,$n) = @_;
  validate_integer($a);
  validate_integer($n);

  if ($n <= 0) {
    return undef if $n == 0;   # standard mod behavior with n = 0
    $n = tobigint($n) if $n <= INTMIN && !ref($n);
    $n = -$n;                  # we use |n|, unlike modint
  }
  # Easy:
  # Msubmod(0, $a, $n);

  $a = Mmodint($a,$n) if $a >= $n || $a < 0;
  return $a ? $n-$a : 0;
}

# No validation.
sub _negmod {
  my($a,$n) = @_;
  if ($n <= 0) {
    return undef if $n == 0;
    $n = tobigint($n) if $n <= INTMIN && !ref($n);
    $n = -$n;
  }
  $a = Mmodint($a,$n) if $a >= $n || $a < 0;
  return $a ? $n-$a : 0;
}

################################################################################
################################################################################

# no validation, x is allowed to be negative, y must be >= 0
sub _gcd_ui {
  my($x, $y) = @_;
  if ($y < $x) { ($x, $y) = ($y, $x); }
  elsif ($x < 0) { $x = -$x; }
  while ($y > 0) {
    ($x, $y) = ($y, $x % $y);
  }
  $x;
}

sub _powerof_ret {
  my($n, $refp) = @_;

  my $k = 2;
  while (1) {
    my $rk;
    my $r = Mrootint($n, $k, \$rk);
    return 0 if $r == 1;
    if ($rk == $n) {
      my $next = _powerof_ret($r, $refp);
      $$refp = $r if !$next && defined $refp;
      $k *= $next if $next != 0;
      return $k;
    }
    $k = Mnext_prime($k);
  }
  0;
}

sub is_power {
  my ($n, $a, $refp) = @_;
  validate_integer($n);
  if (!defined $a) { $a = 0; } else { validate_integer_nonneg($a); }
  croak("is_power third argument not a scalar reference") if defined($refp) && !ref($refp);
  return 0 if abs($n) <= 3 && !$a;

  if ($Math::Prime::Util::_GMPfunc{"is_power"} &&
      ($Math::Prime::Util::GMP::VERSION >= 0.42 ||
       ($Math::Prime::Util::GMP::VERSION >= 0.28 && $n > 0))) {
    my $k = Math::Prime::Util::GMP::is_power($n,$a);
    return 0 unless $k > 0;
    if (defined $refp) {
      $a = $k unless $a;
      my $isneg = ($n < 0);
      $n =~ s/^-// if $isneg;
      $$refp = Mrootint($n, $a);
      $$refp = reftyped($_[0], $$refp) if $$refp > INTMAX;
      $$refp = Mnegint($$refp) if $isneg;
    }
    return $k;
  }

  if ($a != 0) {
    return 1 if $a == 1;                  # Everything is a 1st power
    return 0 if $n < 0 && $a % 2 == 0;    # Negative n never an even power
    if ($a == 2) {
      if (_is_perfect_square($n)) {
        $$refp = Msqrtint($n) if defined $refp;
        return 1;
      }
    } else {

      my @rootmask = (
  0x00000000,0x00000000,0xfdfcfdec,0x54555454,0xfffcfffc,           # 0-4
  0x55555554,0xfdfdfdfc,0x55555554,0xfffffffc,0x55555554,0xfdfdfdfc,# 5-10
  0x55555554,0xfffdfffc,0xd5555556,0xfdfdfdfc,0xf57d57d6,0xfffffffc,# 11-16
  0xffffd556,0xfdfdfdfe,0xd57ffffe,0xfffdfffc,0xffd7ff7e,0xfdfdfdfe,# 17-22
  0xffffd7fe,0xfffffffc,0xffffffd6,0xfdfffdfe,0xd7fffffe,0xfffdfffe,# 23-28
  0xfff7fffe,0xfdfffffe,0xfffff7fe,0xfffffffc,0xfffffff6,0xfffffdfe,# 29-34
  0xf7fffffe,0xfffdfffe,0xfff7fffe,0xfdfffffe,0xfffff7fe,0xfffffffc # 35-40
      );
      return 0 if $a <= 40 && (1 << ($n & 31)) & $rootmask[$a];

      my $RK;
      if ($n >= 0) {
        my $root = Mrootint($n, $a, \$RK);
        if ($RK == $n) { $$refp = $root if defined $refp;  return 1; }
      } else {
        my $N = Mnegint($n);
        my $root = Mrootint($N, $a, \$RK);
        if ($RK == $N) { $$refp = Mnegint($root) if defined $refp;  return 1; }
      }
    }
    return 0;
  }

  my $negn = $n < 0;
  $n = Mnegint($n) if $negn;
  my $k = _powerof_ret($n, $refp);
  return 0 if $k < 2;
  if ($negn && $k % 2 == 0) {
    my $v = Mvaluation($k, 2);
    $k >>= $v;
    return 0 if $k < 2;
    $$refp = Mpowint($$refp, Mpowint(2,$v)) if defined $refp;
  }
  $$refp = Mnegint($$refp) if $negn && defined $refp;
  $k;
}

sub is_square {
  my($n) = @_;
  return 0 if $n < 0;
  #Mis_power($n,2);
  validate_integer($n);
  _is_perfect_square($n);
}

sub is_prime_power {
  my ($n, $refp) = @_;
  validate_integer($n);
  croak("is_prime_power second argument not a scalar reference") if defined($refp) && !ref($refp);
  return 0 if $n <= 1;

  if (Mis_prime($n)) { $$refp = $n if defined $refp; return 1; }
  my $r;
  my $k = Mis_power($n,0,\$r);
  if ($k) {
    $r = _bigint_to_int($r) if ref($r) && $r <= INTMAX;
    return 0 unless Mis_prime($r);
    $$refp = $r if defined $refp;
  }
  $k;
}

sub is_gaussian_prime {
  my($a,$b) = @_;
  validate_integer_abs($a);
  validate_integer_abs($b);
  return ((($b % 4) == 3) ? Mis_prime($b) : 0) if $a == 0;
  return ((($a % 4) == 3) ? Mis_prime($a) : 0) if $b == 0;
  Mis_prime( Maddint( Mmulint($a,$a), Mmulint($b,$b) ) );
}

sub is_polygonal {
  my ($n, $k, $refp) = @_;
  validate_integer($n);
  validate_integer_nonneg($k);
  croak("is_polygonal third argument not a scalar reference") if defined($refp) && !ref($refp);
  croak("is_polygonal: k must be >= 3") if $k < 3;
  return 0 if $n < 0;
  if ($n <= 1) { $$refp = $n if defined $refp; return 1; }

  if ($Math::Prime::Util::_GMPfunc{"polygonal_nth"}) {
    my $nth = Math::Prime::Util::GMP::polygonal_nth($n, $k);
    return 0 unless $nth;
    $$refp = reftyped($_[0], $nth) if defined $refp;
    return 1;
  }

  my($D,$R);
  if ($k == 4) {
    return 0 unless _is_perfect_square($n);
    $$refp = Msqrtint($n) if defined $refp;
    return 1;
  }
  if ($n <= MPU_HALFWORD && $k <= MPU_HALFWORD) {
    $D = ($k==3) ? 1+($n<<3) : (8*$k-16)*$n + ($k-4)*($k-4);
    return 0 unless _is_perfect_square($D);
    $D = $k-4 + Msqrtint($D);
    $R = 2*$k-4;
  } else {
    if ($k == 3) {
      $D = Maddint(1, Mmulint($n, 8));
    } else {
      $D = Maddint(Mmulint($n, Mmulint(8, $k) - 16), Mmulint($k-4,$k-4));
    }
    return 0 unless _is_perfect_square($D);
    $D = Maddint( Msqrtint($D), $k-4 );
    $R = Mmulint(2, $k) - 4;
  }
  return 0 if ($D % $R) != 0;
  $$refp = $D / $R if defined $refp;
  1;
}

sub is_sum_of_squares {
  my($n, $k) = @_;
  validate_integer_abs($n);
  if (defined $k) { validate_integer_nonneg($k); }
  else            { $k = 2; }

  return ($n == 0) ? 1 : 0 if $k == 0;
  return 1 if $k > 3;
  return _is_perfect_square($n) if $k == 1;

  return 1 if $n < 3;

  if ($k == 3) {
    my $tz = Mvaluation($n,2);
    return 1 if ($tz & 1) == 1;
    return 1 unless Mis_congruent(Mrshiftint($n,$tz), 7, 8);
    return 0;
  }

  # k = 2
  while (($n % 2) == 0) { $n >>= 1; }
  return 0 if ($n % 4) == 3;

  foreach my $F (Mfactor_exp($n)) {
    my($f,$e) = @$F;
    return 0 if ($e & 1) == 1 && ($f % 4) == 3;
  }
  1;
}

sub cornacchia {
  my($d, $n) = @_;
  validate_integer_nonneg($d);
  validate_integer_nonneg($n);

  return (0,0) if $n == 0;
  if ($d == 0) {
    return undef unless _is_perfect_square($n);
    return (Msqrtint($n), 0);
  }

  if (Mis_prime($n)) {
    my ($u,$rk);
    my $negd = _negmod($d,$n);
    return undef if Mkronecker($negd, $n) == -1;
    $u = _sqrtmod_prime($negd, $n);
    return undef unless defined $u;
    $u = $n-$u if $u > ($n>>1);
    {
      my $l = Msqrtint($n);
      my($a, $b) = ($n, $u);
      while ($a > $l) {
        ($a,$b) = ($b, $a % $b);
      }
      $rk = $a;
    }
    $u = _negmod(Mmulmod($rk,$rk,$n),$n);
    $u = (($u % $d) == 0) ? Mdivint($u,$d) : 0;
    return ($rk, Msqrtint($u)) if $u && _is_perfect_square($u);
    return undef;
  }

  my $limu = Msqrtint(Mdivint($n,$d));
  for my $u (0 .. $limu) {
    my $t = $n - Mvecprod($d,$u,$u);
    return (Msqrtint($t), $u) if _is_perfect_square($t);
  }
  undef;
}

sub is_congruent_number {
  my($n) = @_;
  validate_integer_nonneg($n);

  return ($n >= 5 && $n <= 7) if $n < 13;

  my $n8 = $n % 8;
  return 1 if $n8 == 5 || $n8 == 6 || $n8 == 7;

  if (!Mis_square_free($n)) {
    my $N = 1;
    foreach my $f (Mfactor_exp($n)) {
      my($p,$e) = @$f;
      $N = Mmulint($N,$p) if ($e % 2) == 1;
    }
    return is_congruent_number($N);
  }

  my $ndiv2 = Mrshiftint($n);

  if (Mis_even($n) && Mis_prime($ndiv2)) {
    my $p = $ndiv2;
    my $p8 = $p % 8;
    return 1 if $p8 == 3 || $p8 == 7;
    return 0 if $p8 == 5 || ($p % 16) == 9;
  } elsif (Mis_prime($n)) {
    return 0 if $n8 == 3;
    return 1 if $n8 == 5 || $n8 == 7;

    my $r = _sqrtmod_prime(2, $n);
    return 0 if defined $r && Mkronecker(1+$r, $n) == -1;
  } elsif (1) {
    my @factors = Mfactor($n);
    if (scalar(@factors) == 2) {
      my($p, $q) = ($factors[0], $factors[1]);
      my($p8, $q8) = ($p % 8, $q %8);
      return 0 if $p8 == 3 && $q8 == 3;
      return 0 if $p8 == 1 && $q8 == 3 && kronecker($p,$q) == -1;
      return 0 if $p8 == 3 && $q8 == 1 && kronecker($q,$p) == -1;
    } elsif (scalar(@factors) == 3 && $factors[0] == 2) {
      my($p, $q) = ($factors[1], $factors[2]);
      my($p8, $q8) = ($p % 8, $q %8);
      return 0 if $p8 == 5 && $q8 == 5;
      return 0 if $p8 == 1 && $q8 == 5 && kronecker($p,$q) == -1;
      return 0 if $p8 == 5 && $q8 == 1 && kronecker($q,$p) == -1;
    }
  }

  # General test
  my @sols = (0,0);
  if (Mis_odd($n)) {
    my $limz = Msqrtint($n >> 3);
    foreach my $z (0 .. $limz) {
      my $zsols = 0;
      my $n8z = $n - 8*$z*$z;
      my $limy = Msqrtint($n8z >> 1);
      foreach my $y (0 .. $limy) {
        my $x = $n8z - 2*$y*$y;
        $zsols += 1 << (1 + ($y>0) + ($z>0))
          if _is_perfect_square($x);
      }
      $sols[$z % 2] += $zsols;
    }
  } else {
    my $limz = Msqrtint($ndiv2 >> 3);
    foreach my $z (0 .. $limz) {
      my $zsols = 0;
      my $n8z = $ndiv2 - 8*$z*$z;   # ndiv2 odd => n8z is odd
      my $limx = Msqrtint($n8z);
      for (my $x = 1; $x <= $limx; $x += 2) {
        my $y = $n8z - $x*$x;
        $zsols += 1 << (1 + ($y>0) + ($z>0))
          if $y == 0 || _is_perfect_square($y);
      }
      $sols[$z % 2] += $zsols;
    }
  }
  return ($sols[0] == $sols[1]) ? 1 : 0;
}

sub is_perfect_number {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n <= 0;

  if (Mis_even($n)) {
    my $v = Mvaluation($n,2);
    my $m = Mrshiftint($n, $v);
    return 0 if Mrshiftint($m,$v) != 1;
    return 0 if Math::Prime::Util::hammingweight($m) != $v+1;
    return Math::Prime::Util::is_mersenne_prime($v+1);
  }

  # N is odd.  See https://www.lirmm.fr/~ochem/opn/
  return 0 if length($n) <= 2200;
  return 0 unless Mis_divisible($n, 105);
  return 0 unless Mis_congruent($n,  1, 12)
               || Mis_congruent($n,117,468)
               || Mis_congruent($n, 81, 324);
  Mcmpint($n,Msubint(Mdivisor_sum($n),$n)) == 0;
}

sub valuation {
  my($n, $k) = @_;
  # The validation in PP is 2x more time than our actual work.
  validate_integer_abs($n);
  validate_integer_positive($k);
  croak "valuation: k must be > 1" if $k <= 1;

  return if $k < 2;
  return (undef,0)[$n] if $n <= 1;
  my $v = 0;
  if ($k == 2) { # Accelerate power of 2
    my $s;
    if (!ref($n)) {
      return 0 if $n & 1;
      return 1 if $n & 2;
      return 2 if $n & 4;
      $s = sprintf("%b",$n);
    } elsif (ref($n) eq 'Math::BigInt') {
      $s = substr($n->as_bin,2);
    } else {
      $s = substr(Math::BigInt->new("$n")->as_bin,2);
    }
    return length($s) - rindex($s,'1') - 1;
  }
  while ( !($n % $k) ) {
    $n /= $k;
    $v++;
  }
  $v;
}

sub hammingweight {
  my $n = shift;
  return 0 + (Mtodigitstring($n,2) =~ tr/1//);
}

my @_digitmap = (0..9, 'a'..'z');
my %_mapdigit = map { $_digitmap[$_] => $_ } 0 .. $#_digitmap;
sub _splitdigits {
  my($n, $base, $len) = @_;    # n is num or bigint, base is in range
  validate_integer_nonneg($n);
  my @d;
  if ($base == 10) {
    @d = split(//,"$n");
  } elsif ($base == 2) {
    @d = split(//,substr(Math::BigInt->new("$n")->as_bin,2));
  } elsif ($base == 16) {
    @d = map { $_mapdigit{$_} } split(//,substr(Math::BigInt->new("$n")->as_hex,2));
  } else {
    # The validation turned n into a bigint if necessary
    while ($n >= 1) {
      my $rem = $n % $base;
      unshift @d, $rem;
      $n = ($n-$rem)/$base;    # Always an exact division
    }
  }
  if ($len >= 0 && $len != scalar(@d)) {
    while (@d < $len) { unshift @d, 0; }
    while (@d > $len) { shift @d; }
  }
  @d;
}

sub todigits {
  my($n,$base,$len) = @_;
  validate_integer_abs($n);
  $base = 10 unless defined $base;
  $len = -1 unless defined $len;
  die "Invalid base: $base" if $base < 2;
  return if $n == 0;
  _splitdigits($n, $base, $len);
}

sub todigitstring {
  my($n,$base,$len) = @_;
  validate_integer($n);
  $base = 10 unless defined $base;
  croak "Invalid base for string: $base" if $base < 2 || $base > 36;
  $len = -1 unless defined $len;
  $n =~ s/^-//;

  return "" if $len == 0 || $n == 0;

  if ($n < INTMAX) {
    if ($base != 2 && $base != 8 && $base != 16) {
      return join "", _splitdigits($n, $base, $len)  if $base <= 10;
      return join "", map { $_digitmap[$_] } _splitdigits($n, $base, $len);
    }
    my $s;
    $s = sprintf("%b",$n)  if $base ==  2;
    $s = sprintf("%o",$n)  if $base ==  8;
    $s = sprintf("%x",$n)  if $base == 16;
    if ($len > 0) {
      $s = substr($s,0,$len);
      $s = '0' x ($len-length($s)) . $s if length($s) < $len;
    }
    return $s;
  }

  $n = tobigint($n) unless ref($n);
  my $refn = ref($n);
  my $s;

  if ($refn eq 'Math::GMPz') {
    $s = Math::GMPz::Rmpz_get_str($n,$base);
  } elsif ($refn eq 'Math::GMP') {
    $s = Math::GMP::get_str_gmp($n,$base);
  } else {
    $n = Math::BigInt->new("$n") if $refn ne 'Math::BigInt';
    $s = $n->to_base($base);
  }
  if ($len > 0) {
    $s = substr($s,0,$len);
    $s = '0' x ($len-length($s)) . $s if length($s) < $len;
  }
  return lc($s);
}

sub _FastIntegerInput {
  my($digits, $B) = @_;
  return 0 if scalar(@$digits) == 0;
  return $digits->[0] if scalar(@$digits) == 1;
  my $L = [reverse @$digits];
  my $k = scalar(@$L);
  while ($k > 1) {
    my @T;
    for my $i (1 .. $k>>1) {
      my $x = $L->[2*$i-2];
      my $y = $L->[2*$i-1];
      push(@T, Maddint($x, Mmulint($B, $y)));
    }
    push(@T, $L->[$k-1]) if ($k&1);
    $L = \@T;
    $B = Mmulint($B, $B);
    $k = ($k+1) >> 1;
  }
  $L->[0];
}

sub fromdigits {
  my($r, $base) = @_;
  $base = 10 unless defined $base;
  my $refr = ref($r);

  if ($refr && $refr !~ /^Math::/) {
    croak "fromdigits: first argument must be a string or array reference"
      unless $refr eq 'ARRAY';
    return _FastIntegerInput($r,$base);
  }

  my $n;
  $r =~ s/^0*//;
  return 0 if $r eq "";
  { # Validate string
    my $cmap = substr("0123456789abcdefghijklmnopqrstuvwxyz",0,$base);
    croak "Invalid digit for base $base" if $r =~ /[^$cmap]/i;
  }
  if (defined $_BIGINT && $_BIGINT =~ /^Math::(GMPz|GMP)$/) {
    $n = $_BIGINT->new($r, $base);
  } else {
    # from_base is 2x slower than calling the method directly (TODO file an RT)
    if    ($base ==  2) { $n = Math::BigInt->from_bin($r); }
    elsif ($base ==  8) { $n = Math::BigInt->from_oct($r); }
    elsif ($base == 10) { $n = Math::BigInt->new($r); }
    elsif ($base == 16) { $n = Math::BigInt->from_hex($r); }
    else                { $n = Math::BigInt->from_base($r,$base); }
    $n = tobigint($n) if defined $_BIGINT && $_BIGINT ne 'Math::BigInt';
  }
  return $n <= INTMAX ? _bigint_to_int($n) : $n;
}

sub _validate_zeckendorf {
  my $s = shift;
  if ($s ne '0') {
    croak "fromzeckendorf: expected binary string"
      unless $s =~ /^1[01]*\z/;
    croak "fromzeckendorf: expected binary string in canonical Zeckendorf form"
      if $s =~ /11/;
  }
  1;
}

sub fromzeckendorf {
  my($s) = @_;
  _validate_zeckendorf($s);

  my($n, $fb, $fc) = (0, 1, 1);
  for my $c (split(//,reverse $s)) {
    $n = Maddint($n,$fc) if $c eq '1';
    ($fb, $fc) = ($fc, Maddint($fb,$fc));
  }
  $n;
}

sub tozeckendorf {
  my($n) = @_;
  validate_integer_nonneg($n);
  return '0' if $n == 0;

  my($rn, $s, $fa, $fb, $fc) = ($n, '', 0, 1, 1);
  my($i, $k);
  for ($k = 2; $fc <= $rn; $k++) {
    ($fa, $fb, $fc) = ($fb, $fc, Maddint($fb,$fc));
  }
  for ($i = $k-1; $i >= 2; $i--) {
    ($fc, $fb, $fa) = ($fb, $fa, Msubint($fb,$fa));
    if ($fc <= $rn) {
      $rn = Msubint($rn, $fc);
      $s .= '1';
    } else {
      $s .= '0';
    }
  }
  # croak "wrong tozeckendorf $n" unless $n == fromzeckendorf($s);
  $s;
}


sub sqrtint {
  my($n) = @_;
  validate_integer_nonneg($n);
  return int(sqrt("$n")) if $n <= 562949953421312;  # 2^49 safe everywhere

  my $refn = ref($n);
  my $R;

  if ($refn eq 'Math::BigInt') {
    $R = $n->copy->bsqrt;
  } elsif ($refn eq 'Math::GMPz') {
    $R = Math::GMPz->new();
    Math::GMPz::Rmpz_sqrt($R, $n);
  } elsif ($refn eq 'Math::GMP') {
    $R = $n->bsqrt();
  } else {
    $R = Math::BigInt->new("$n")->bsqrt;
  }
  $R = _bigint_to_int($R) if $R <= INTMAX;
  $R;
}

sub rootint {
  my ($n, $k, $refp) = @_;
  validate_integer_nonneg($n);
  validate_integer_positive($k);
  croak("rootint: third argument not a scalar reference") if defined $refp && !ref($refp);

  if ($k == 1) {
    $$refp = $n if defined $refp;
    return $n;
  }
  if (!ref($n)) {  # native integer
    if ($n == 0) {
      $$refp = 0 if defined $refp;
      return 0;
    }
    if ($k == 2 && $n <= 562949953421312) {
      my $R = int(sqrt($n));
      $$refp = $R*$R if defined $refp;
      return $R;
    }
    if ($k >= MPU_MAXBITS || $n >> $k == 0) {
      $$refp = 1 if defined $refp;
      return 1;
    }
    my $R = int($n ** (1/$k));  # Could be off by +/-1.
    my $F = $n <= 562949953421312 ? $R**$k : powint($R,$k);
    if ($F > $n) {
      $R--;
      $F = $n <= 562949953421312 ? $R**$k : powint($R,$k);
    } else {
      my $F1 = $n <= 562949953421312 ? ($R+1)**$k : powint($R+1,$k);
      if ($F1 <= $n) {
        $R++;
        $F = $F1;
      }
    }
    $$refp = $F if defined $refp;
    return $R;
  }

  # It's unclear whether we should add GMPfunc here.  We want it in logint
  # because it's slow or not included in Perl bigint classes.

  my $refn = ref($n);
  my $R;
  if ($refn eq 'Math::BigInt') {
    $R = $n->copy->broot($k);
  } elsif ($refn eq 'Math::GMPz') {
    $R = Math::GMPz->new();
    Math::GMPz::Rmpz_root($R, $n, $k);
  } elsif ($refn eq 'Math::GMP') {
    $R = $n->broot($k);
  } else {
    $R = Math::BigInt->new("$n")->broot($k);
  }
  $R = _bigint_to_int($R) if $R <= INTMAX;
  $$refp = Mpowint($R,$k) if defined $refp;
  $R;
}

sub _logint {
  my($n,$b) = @_;
  return 0 if $n < $b;
  return length("$n")-1 if $b == 10;
  if ($n < INTMAX) {
    return length(sprintf("%b",$n))-1 if $b ==  2;
    return length(sprintf("%o",$n))-1 if $b ==  8;
    return length(sprintf("%x",$n))-1 if $b == 16;
  }
  my $l;
  if (length("$n") > 150) {
    # Reduce size so native log works
    my $N = substr($n,0,80);
    my $reddigits = length("$n") - length($N);
    $l = log($N) + 2.302585092994045684*$reddigits;
  } else {
    $l = log("$n");
  }
  $l /= log($b);

  # Just in case something failed, escape via using Math::BigInt's blog
  if ($l == MPU_INFINITY || !defined($l<=>MPU_INFINITY)) {
    my $R = Math::BigInt->new("$n")>copy->blog($b);
    $R = _bigint_to_int($R) if $R <= INTMAX;
    return $R;
  }

  my $R = int($l);
  if ($R != int($l+1e-7) || $R != int($l-1e-7)) {
    my $BR = Mpowint($b,$R);
    if ($BR > $n) {
      $R--;
    } elsif ($BR < $n) {
      my $BRB = Mmulint($BR, $b);
      $R++ if $BRB <= $n;
    }
  }
  $R;
}

sub logint {
  my ($n, $b, $refp) = @_;
  croak("logint third argument not a scalar reference") if defined($refp) && !ref($refp);

  if ($Math::Prime::Util::_GMPfunc{"logint"}) {
    my $e = Math::Prime::Util::GMP::logint($n, $b);
    if (defined $refp) {
      # logint in 0.47, powmod in 0.36, powint in 0.52
      my $r = Math::Prime::Util::GMP::powmod($b, $e, $n);
      $r = $n if $r == 0;
      $$refp = reftyped($_[0], $r);
    }
    return reftyped($_[0], $e);
  }

  validate_integer_positive($n);
  validate_integer_nonneg($b);

  my $log = _logint($n,$b);
  $$refp = Mpowint($b,$log) if defined $refp;
  return $log;
}

# Seidel (Luschny), core using Trizen's simplications from Math::AnyNum.
# http://oeis.org/wiki/User:Peter_Luschny/ComputationAndAsymptoticsOfBernoulliNumbers#Bernoulli_numbers__after_Seidel
sub _bernoulli_seidel {
  my($n) = @_;
  return (1,1) if $n == 0;
  return (0,1) if $n > 1 && $n % 2;

  my @D = (0, 1, map { 0} 1 .. ($n>>1)-1);
  my ($h, $w) = (1, 1);

  foreach my $i (0 .. $n-1) {
    if ($w ^= 1) {
      $D[$_] = Maddint($D[$_],$D[$_-1]) for 1.. $h-1;
    } else {
      $w = $h++;
      $D[$w] = Maddint($D[$w],$D[$w+1]) while --$w;
    }
  }
  my $num = $D[$h-1];
  my $den = Msubint(Mpowint(2,$n+1),2);
  my $gcd = Mgcd($num,$den);
  ($num,$den) = map { Mdivint($_,$gcd) } ($num,$den) if $gcd > 1;
  $num = Mnegint($num) if ($n % 4) == 0;
  ($num,$den);
}

sub bernfrac {
  my $n = shift;
  return (1,1) if $n == 0;
  return (1,2) if $n == 1;    # We're choosing 1/2 instead of -1/2
  return (0,1) if $n < 0 || $n & 1;

  # We should have used one of the GMP functions before coming here.

  _bernoulli_seidel($n);
}

sub stirling {
  my($n, $m, $type) = @_;
  return 1 if $m == $n;
  return 0 if $n == 0 || $m == 0 || $m > $n;
  $type = 1 unless defined $type;
  croak "stirling type must be 1, 2, or 3" unless $type == 1 || $type == 2 || $type == 3;
  if ($m == 1) {
    return 1 if $type == 2;
    return Mfactorial($n) if $type == 3;
    return Mfactorial($n-1) if $n & 1;
    return Mvecprod(-1, Mfactorial($n-1));
  }
  return reftyped($_[0], Math::Prime::Util::GMP::stirling($n,$m,$type))
    if $Math::Prime::Util::_GMPfunc{"stirling"};
  # Go through vecsum with quoted negatives to make sure we don't overflow.
  my $s;
  if ($type == 3) {
    $s = Mvecprod( Mbinomial($n,$m), Mbinomial($n-1,$m-1), Mfactorial($n-$m) );
  } elsif ($type == 2) {
    my @terms;
    for my $j (1 .. $m) {
      my $t = Mmulint(
                Mpowint($j,$n),
                Mbinomial($m,$j)
              );
      $t = Mnegint($t) if ($m-$j) & 1;
      push @terms, $t;
    }
    $s = Mvecsum(@terms) / Mfactorial($m);
  } else {
    my @terms;
    for my $k (1 .. $n-$m) {
      my $t = Mvecprod(
        Mbinomial($k + $n - 1, $k + $n - $m),
        Mbinomial(2 * $n - $m, $n - $k - $m),
        Mstirling($k - $m + $n, $k, 2),
      );
      $t = Mnegint($t) if $k & 1;
      push @terms, $t;
    }
    $s = Mvecsum(@terms);
  }
  $s;
}

sub _harmonic_split { # From Fredrik Johansson
  my($a,$b) = @_;
  return (1, $a) if $b-$a == 1;
  return (Mvecsum($a,$a,1), Maddint(Mmulint($a,$a),$a)) if $b-$a == 2;
  my $m = Mrshiftint(Maddint($a,$b),1);
  my ($p,$q) = _harmonic_split($a, $m);
  my ($r,$s) = _harmonic_split($m, $b);
  (Maddint(Mmulint($p,$s),Mmulint($q,$r)), Mmulint($q,$s));
}

sub harmfrac {
  my($n) = @_;
  return (0,1) if $n <= 0;
  my($p,$q) = _harmonic_split(1, Maddint($n,1));
  my $gcd = Mgcd($p,$q);
  ($p,$q) = map { Mdivint($_,$gcd) } ($p,$q) if $gcd > 1;
  ($p,$q);
}

sub harmreal {
  my($n, $precision) = @_;

  do { require Math::BigFloat; Math::BigFloat->import(); } unless defined $Math::BigFloat::VERSION;
  return Math::BigFloat->bzero if $n <= 0;

  # Use asymptotic formula for larger $n if possible.  Saves lots of time if
  # the default Calc backend is being used.
  {
    my $sprec = $precision;
    $sprec = Math::BigFloat->precision unless defined $sprec;
    $sprec = 40 unless defined $sprec;
    if ( ($sprec <= 23 && $n >    54) ||
         ($sprec <= 30 && $n >   348) ||
         ($sprec <= 40 && $n >  2002) ||
         ($sprec <= 50 && $n > 12644) ) {
      $n = Math::BigFloat->new($n, $sprec+5);
      my($n2, $one, $h) = ($n*$n, Math::BigFloat->bone, Math::BigFloat->bzero);
      my $nt = $n2;
      my $eps = Math::BigFloat->new(10)->bpow(-$sprec-4);
      foreach my $d (-12, 120, -252, 240, -132, 32760, -12, 8160, -14364, 6600, -276, 65520, -12) { # OEIS A006593
        my $term = $one/($d * $nt);
        last if $term->bacmp($eps) < 0;
        $h += $term;
        $nt *= $n2;
      }
      $h->badd(scalar $one->copy->bdiv(2*$n));
      $h->badd(_Euler($sprec));
      $h->badd($n->copy->blog);
      $h->round($sprec);
      return $h;
    }
  }

  my($num,$den) = Math::Prime::Util::harmfrac($n);
  # Note, with Calc backend this can be very, very slow
  scalar Math::BigFloat->new($num)->bdiv($den, $precision);
}

sub is_pseudoprime {
  my($n, @bases) = @_;
  validate_integer($n);
  return 0 if $n < 0;
  @bases = (2) if scalar(@bases) == 0;
  return 0+($n >= 2) if $n < 3;

  foreach my $a (@bases) {
    croak "Base $a is invalid" if $a < 2;
    $a = $a % $n if $a >= $n;
    return 0 unless $a == 1 || Mpowmod($a, $n-1, $n) == 1;
  }
  1;
}

sub is_euler_pseudoprime {
  my($n, @bases) = @_;
  validate_integer($n);
  return 0 if $n < 0;
  @bases = (2) if scalar(@bases) == 0;
  return 0+($n >= 2) if $n < 3;
  return 0 if ($n % 2) == 0;

  foreach my $a (@bases) {
    croak "Base $a is invalid" if $a < 2;
    $a = $a % $n if $a >= $n;
    my $j = Mkronecker($a, $n);
    return 0 if $j == 0;   # gcd(a,n) != 1
    $j = ($j > 0) ? 1 : $n-1;
    return 0 unless Mpowmod($a, ($n-1)>>1, $n) == $j;
  }
  1;
}

sub is_euler_plumb_pseudoprime {
  my($n) = @_;
  validate_integer($n);
  return 0 if $n < 0;
  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0;
  my $nmod8 = $n % 8;
  my $exp = 1 + ($nmod8 == 1);
  my $ap = Mpowmod(2, ($n-1) >> $exp, $n);
  if ($ap ==    1) { return ($nmod8 == 1 || $nmod8 == 7); }
  if ($ap == $n-1) { return ($nmod8 == 1 || $nmod8 == 3 || $nmod8 == 5); }
  0;
}

sub _miller_rabin_2 {
  my($n, $nm1, $s, $d) = @_;
  return 0 if $n < 0;
  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0;

  if (ref($n)) {

    if (!defined $nm1) {
      $nm1 = Msubint($n,1);
      $s = Mvaluation($nm1,2);
      $d = Mrshiftint($nm1,$s);
    }
    my $x = _bi_powmod(2,$d,$n);
    return 1 if $x == 1 || $x == $nm1;
    foreach my $r (1 .. $s-1) {
      $x = Mmulmod($x,$x,$n);
      last if $x == 1;
      return 1 if $x == $nm1;
    }

  } else {

    if (!defined $nm1) {
      $nm1 = $n-1;
      $s = 0;
      $d = $nm1;
      while ( ($d & 1) == 0 ) {
        $s++;
        $d >>= 1;
      }
    }

    if ($n < MPU_HALFWORD) {
      my $x = _native_powmod(2, $d, $n);
      return 1 if $x == 1 || $x == $nm1;
      foreach my $r (1 .. $s-1) {
        $x = ($x*$x) % $n;
        last if $x == 1;
        return 1 if $x == $n-1;
      }
    } else {
      my $x = _powmod(2, $d, $n);
      return 1 if $x == 1 || $x == $nm1;
      foreach my $r (1 .. $s-1) {
        $x = ($x < MPU_HALFWORD) ? ($x*$x) % $n : _mulmod($x, $x, $n);
        last if $x == 1;
        return 1 if $x == $n-1;
      }
    }
  }
  0;
}

sub is_strong_pseudoprime {
  my($n, @bases) = @_;
  validate_integer($n);
  return 0 if $n < 0;
  return _miller_rabin_2($n) if scalar(@bases) == 0;

  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0;

  my @newbases;
  for my $a (@bases) {
    croak "Base $a is invalid" if $a < 2;
    $a %= $n if $a >= $n;
    next if $a <= 1 || $a == $n-1;
    if ($a == 2) {
      return 0 unless _miller_rabin_2($n);
      next;
    }
    push @newbases, $a;
  }
  return 1 if scalar(@newbases) == 0;
  @bases = @newbases;

  if (ref($n)) {

    my $nm1 = Msubint($n,1);
    my $s = Mvaluation($nm1,2);
    my $d = Mrshiftint($nm1,$s);

    foreach my $ma (@bases) {
      my $x = Mpowmod($ma,$d,$n);
      next if $x == 1 || $x == $nm1;
      foreach my $r (1 .. $s-1) {
        $x = Mmulmod($x,$x,$n);
        return 0 if $x == 1;
        last if $x == $nm1;
      }
      return 0 if $x != $nm1;
    }

  } else {

   my $s = 0;
   my $d = $n - 1;
   while ( ($d & 1) == 0 ) {
     $s++;
     $d >>= 1;
   }

   if ($n < MPU_HALFWORD) {
    foreach my $ma (@bases) {
      my $x = _native_powmod($ma, $d, $n);
      next if ($x == 1) || ($x == ($n-1));
      foreach my $r (1 .. $s-1) {
        $x = ($x*$x) % $n;
        return 0 if $x == 1;
        last if $x == $n-1;
      }
      return 0 if $x != $n-1;
    }
   } else {
    foreach my $ma (@bases) {
      my $x = _powmod($ma, $d, $n);
      next if ($x == 1) || ($x == ($n-1));
      foreach my $r (1 .. $s-1) {
        $x = ($x < MPU_HALFWORD) ? ($x*$x) % $n : _mulmod($x, $x, $n);
        return 0 if $x == 1;
        last if $x == $n-1;
      }
      return 0 if $x != $n-1;
    }
   }

  }
  1;
}


# Calculate Kronecker symbol (a|b).  Cohen Algorithm 1.4.10.
# Extension of the Jacobi symbol, itself an extension of the Legendre symbol.
sub kronecker {
  my($a, $b) = @_;
  return (abs($a) == 1) ? 1 : 0  if $b == 0;
  my $k = 1;
  if ($b % 2 == 0) {
    return 0 if $a % 2 == 0;
    my $v = 0;
    do { $v++; $b /= 2; } while $b % 2 == 0;
    $k = -$k if $v % 2 == 1 && ($a % 8 == 3 || $a % 8 == 5);
  }
  if ($b < 0) {
    $b = -$b;
    $k = -$k if $a < 0;
  }
  if ($a < 0) { $a = -$a; $k = -$k if $b % 4 == 3; }
  $b = _bigint_to_int($b) if ref($b) && $b <= INTMAX;
  $a = _bigint_to_int($a) if ref($a) && $a <= INTMAX;
  # Now:  b > 0, b odd, a >= 0
  while ($a != 0) {
    if ($a % 2 == 0) {
      my $v = 0;
      do { $v++; $a /= 2; } while $a % 2 == 0;
      $k = -$k if $v % 2 == 1 && ($b % 8 == 3 || $b % 8 == 5);
    }
    $k = -$k if $a % 4 == 3 && $b % 4 == 3;
    ($a, $b) = ($b % $a, $a);
    # If a,b are bigints and now small enough, finish as native.
    return $k * kronecker(_bigint_to_int($a),_bigint_to_int($b))
      if $a <= INTMAX && $b <= INTMAX && ref($a) && ref($b);
  }
  return ($b == 1) ? $k : 0;
}

sub is_qr {
  my($a, $n) = @_;
  validate_integer($a);
  validate_integer_abs($n);

  # return (defined Math::Prime::Util::sqrtmod($a,$n)) ? 1 : 0;

  return (undef,1,1)[$n] if $n <= 2;
  $a = Mmodint($a,$n);
  return 1 if $a <= 1;

  return 0+(Mkronecker($a,$n) == 1) if Mis_prime($n);

  foreach my $f (Mfactor_exp($n)) {
    my($p,$e) = @$f;
    next if $e == 1 && Mkronecker($a,$p) == 1;
    return 0 unless defined _sqrtmod_prime_power($a,$p,$e);
  }
  1;
}

sub _binomialu {
  my($r, $n, $k) = (1, @_);
  return ($k == $n) ? 1 : 0 if $k >= $n;
  $k = $n - $k if $k > ($n >> 1);
  foreach my $d (1 .. $k) {
    if ($r >= int(INTMAX/$n)) {
      my($g, $nr, $dr);
      $g = _gcd_ui($n, $d);   $nr = int($n/$g);   $dr = int($d/$g);
      $g = _gcd_ui($r, $dr);  $r  = int($r/$g);   $dr = int($dr/$g);
      return 0 if $r >= int(INTMAX/$nr);
      $r *= $nr;
      $r = int($r/$dr);
    } else {
      $r *= $n;
      $r = int($r/$d);
    }
    $n--;
  }
  $r;
}

sub binomial {
  my($n, $k) = @_;
  _validate_integer($n);
  _validate_integer($k);

  # 1. Try GMP
  return reftyped($_[0], Math::Prime::Util::GMP::binomial($n,$k))
    if $Math::Prime::Util::_GMPfunc{"binomial"} &&
       ($Math::Prime::Util::GMP::VERSION >= 0.53 || ($n >= 0 && $k >= 0 && $n < 4294967296 && $k < 4294967296));

  # 2. Exit early for known 0 cases, and adjust k to be positive.
  if ($n >= 0) {  return 0 if $k < 0 || $k > $n;  }
  else         {  return 0 if $k < 0 && $k > $n;  }
  $k = $n - $k if $k < 0;

  # TODO: consider reflection for large k (e.g. k=n-2 => k=2)
  # Also, be careful with large n and k with bigints.

  my $r;

  # 3. Try to do in integer Perl
  if (!ref($n)) {
    if ($n >= 0) {
      $r = _binomialu($n, $k);
      return $r  if $r > 0 && $r eq int($r);
    } else {
      $r = _binomialu(-$n+$k-1, $k);
      if ($r > 0 && $r eq int($r)) {
        return $r   if !($k & 1);
        return Mnegint($r);
      }
    }
  }

  # 4. Overflow.  Solve using Math::BigInt
  return 1 if $k == 0;                   # Work around bug in old
  return $n if $k == 1 || $k == $n-1;    # Math::BigInt (fixed in 1.90)

  my $R;
  $n = tobigint($n) unless ref($n);

  # Older Math::BigInt isn't right for negative n.  Adjust now.
  my $negate = 0;
  if ($n < 0) {
    $n = -$n + ($k-1);
    $negate = 1 if $k & 1;
  }

  if (defined $Math::GMPz::VERSION) {
    $R = Math::GMPz->new();
    Math::GMPz::Rmpz_bin_ui($R, Math::GMPz->new($n), $k);
  } elsif (defined $Math::GMP::VERSION && $n < 4294967296) {
    # This will silently coerce inputs to C 'long' type.
    $R = Math::GMP::bnok("$n","$k");
  } elsif ($n > INTMAX && $k < 100) {
    # Incomplete work around problem with Math::BigInt not liking bigint n.
    # Fixed in 2.003003.
    $R = Mdivint(Mfalling_factorial($n,$k),Mfactorial($k));
  } else {
    $R = Math::BigInt::bnok("$n","$k");
  }
  $R = -$R if $negate;
  return $R <= INTMAX && $R <= INTMIN            ?  _bigint_to_int($R)
       : defined $_BIGINT && $_BIGINT eq ref($R) ?  $R
       :                                            tobigint($R);
}

sub binomialmod {
  my($n,$k,$m) = @_;
  validate_integer($n);
  validate_integer($k);
  validate_integer_abs($m);
  return (undef,0)[$m] if $m <= 1;

  # Best if we have it.
  return reftyped($_[2], _gmpcall("binomialmod",$n,$k,$m))
    if $Math::Prime::Util::_GMPfunc{"binomialmod"};

  # Avoid the possible enormously slow bigint creation.
  if ($Math::Prime::Util::_GMPfunc{"binomial"} && $Math::Prime::Util::_GMPfunc{"modint"}) {
    return reftyped($_[2], Math::Prime::Util::GMP::modint(Math::Prime::Util::GMP::binomial($n,$k),$m));
  }

  return 1 if $k == 0 || $k == $n;
  return 0 if $n >= 0 && ($k < 0 || $k > $n);
  return 0 if $n  < 0 && ($k < 0 && $k > $n);
  return 0+!(($n-$k) & $k) if $m == 2;

  # TODO: Lucas split, etc.
  # 1. factorexp
  # 2.   bin[i] = _binomial_lucas_mod_prime_power(n, k, $f, $e)
  # 2a.            _factorialmod_without_prime
  # 3.   chinese(bin, p^e)
  # we can just run the more general code path.

  # Give up.
  return Mmodint(Mbinomial($n,$k),$m);
}

sub _falling_factorial {
  my($n,$m) = @_;
  if ($m <= 1) { return ($m == 0) ? 1 : $n }
  return 0 if $n >= 0 && $m > $n;
  return Mvecprod($n,map { Msubint($n,$_) } 1 .. Msubint($m,1))  if $m < 250;
  Mmulint(Mbinomial($n,$m),Mfactorial($m));
}
sub falling_factorial {
  my($n,$m) = @_;
  validate_integer($n);
  validate_integer_nonneg($m);
  _falling_factorial($n,$m);
}
sub rising_factorial {
  my($n,$m) = @_;
  validate_integer($n);
  validate_integer_nonneg($m);
  _falling_factorial(Mvecsum($n,$m,-1),$m);
}

sub factorial {
  my($n) = @_;
  return (1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600)[$n] if $n <= 12;
  return Math::GMP::bfac($n) if ref($n) eq 'Math::GMP';
  do { my $r = Math::GMPz->new(); Math::GMPz::Rmpz_fac_ui($r,$n); return $r; }
    if ref($n) eq 'Math::GMPz';
  if (Math::BigInt->config()->{lib} !~ /GMP|Pari/) {
    # It's not a GMP or GMPz object, and we have a slow bigint library.
    my $r;
    if (defined $Math::GMPz::VERSION) {
      $r = Math::GMPz->new(); Math::GMPz::Rmpz_fac_ui($r,$n);
    } elsif (defined $Math::GMP::VERSION) {
      $r = Math::GMP::bfac($n);
    } elsif (defined &Math::Prime::Util::GMP::factorial && getconfig()->{'gmp'}) {
      $r = Math::Prime::Util::GMP::factorial($n);
    }
    return reftyped($_[0], $r)    if defined $r;
  }
  # maybe roll our own: https://oeis.org/A000142/a000142.pdf
  my $r = Math::BigInt->new($n)->bfac();
  $r = _bigint_to_int($r) if $r <= INTMAX;
  $r;
}

sub factorialmod {
  my($n,$m) = @_;
  validate_integer($n);
  validate_integer_abs($m);
  return (undef,0)[$m] if $m <= 1;

  return reftyped($_[1], Math::Prime::Util::GMP::factorialmod($n,$m))
    if $Math::Prime::Util::_GMPfunc{"factorialmod"} && $n < ~0;

  return 0 if $n >= $m || $m == 1;

  return factorial($n) % $m if $n <= 10;

  my($F, $N, $m_prime) = (1, $n, Mis_prime($m));

  # Check for Wilson's theorem letting us go backwards
  $n = $m-$n-1 if $m_prime && $n > Mrshiftint($m);
  return ($n == 0) ? ($m-1) : 1  if $n < 2;

  if ($n > 100 && !$m_prime) {   # Check for a composite that leads to zero
    my $maxpk = 0;
    foreach my $f (Mfactor_exp($m)) {
      my $pk = Mmulint($f->[0],$f->[1]);
      $maxpk = $pk if $pk > $maxpk;
    }
    return 0 if $n >= $maxpk;
  }

  my($t,$e);
  Mforprimes( sub {
    ($t,$e) = ($n,0);
    while ($t > 0) {
      $t = int($t/$_);
      $e += $t;
    }
    $F = Mmulmod($F,Mpowmod($_,$e,$m),$m);
  }, 2, $n >> 1);
  Mforprimes( sub {
    $F = Mmulmod($F, $_, $m);
  }, ($n >> 1)+1, $n);

  # Adjust for Wilson's theorem if we used it
  if ($n != $N && $F != 0) {
    $F = Msubmod($m, $F, $m) if !($n & 1);
    $F = Minvmod($F, $m);
  }

  $F;
}

sub subfactorial {
  my($n) = @_;
  validate_integer_nonneg($n);
  if ($n <= 3) { return ($n == 0) ? 1 : $n-1; }
  Mvecsum(map{ Mvecprod((-1)**($n-$_),Mbinomial($n,$_),Mfactorial($_)) }0..$n);
}

sub fubini {
  my($n) = @_;
  validate_integer_nonneg($n);
  return 1 if $n <= 1;
  Mvecsum(map{ Mmulint(Mfactorial($_),Mstirling($n,$_,2)) }1..$n);
}


# Rational maps

sub _rational_cfrac {
  my($num,$den,$non_reduce_ok) = @_;
  my @CF;
  while ($den > 0) {
    my($quo,$rem) = Mtdivrem($num,$den);
    ($num,$den) = ($den,$rem);
    push @CF, $quo;
  }
  croak "Rational must be reduced" unless $num == 1 || $non_reduce_ok;
  @CF;
}

  # https://kconrad.math.uconn.edu/blurbs/ugradnumthy/contfrac-neg-invert.pdf
sub _negcfrac {
  my(@CF) = @_;
  my $neg0 = Mnegint($CF[0]);
  if (@CF == 1) {
    $CF[0] = $neg0;
  } elsif ($CF[1] == 1) {
    splice(@CF, 0, 3, Msubint($neg0,1), Maddint($CF[2],1));
  } else {
    splice(@CF, 0, 2, Msubint($neg0,1), 1, Msubint($CF[1],1));
  }
  @CF;
}

sub contfrac {
  my($num,$den) = @_;
  validate_integer($num);
  validate_integer_positive($den);

  my @CF = _rational_cfrac(Mabsint($num),$den,1);
  return ($num >= 0)  ?  @CF  :  _negcfrac(@CF);
}

sub from_contfrac {
  return (0,1) unless @_;

  my $b0 = shift @_;
  validate_integer($b0);

  my($A0,$A1,$B0,$B1) = (1,$b0,0,1);

  while (@_) {
    my $bi = shift @_;
    validate_integer_positive($bi);
    ($A0,$A1) = ($A1, Maddint(Mmulint($bi,$A1),$A0));
    ($B0,$B1) = ($B1, Maddint(Mmulint($bi,$B1),$B0));
  }
  return ($A1,$B1);
}

sub next_calkin_wilf {
  my($num,$den) = @_;
  validate_integer_positive($num);
  validate_integer_positive($den);
  # Check gcd to ensure a valid CW entry?
  ($den, Mvecprod(2,$den,Mdivint($num,$den)) + $den - $num);
}
sub next_stern_brocot {
  my($num,$den) = @_;
  validate_integer_positive($num);
  validate_integer_positive($den);
  # There should be a better solution
  nth_stern_brocot(Maddint(stern_brocot_n($num,$den),1));
}

sub calkin_wilf_n {
  my($num,$den) = @_;
  validate_integer_positive($num);
  validate_integer_positive($den);

  my @CF = _rational_cfrac($num,$den);
  # Note:  vecsum(@CF) gives the number of bits in the output

  $CF[-1]--;
  my $bitstr = '1';
  $bitstr .= (1-($_%2)) x $CF[$_]  for reverse 0 .. $#CF;
  return Mfromdigits($bitstr,2);
}
sub stern_brocot_n {
  my($num,$den) = @_;
  validate_integer_positive($num);
  validate_integer_positive($den);
  my @CF = _rational_cfrac($num,$den);
  $CF[-1]--;
  my $bitstr = '1';
  $bitstr .= (1-($_%2)) x $CF[$_]  for 0 .. $#CF;
  return Mfromdigits($bitstr,2);
}

sub nth_calkin_wilf {
  my($n) = @_;
  validate_integer_positive($n);

  my @M = (1,0);
  $M[$_] = Mvecsum(@M) for split(//, Mtodigitstring($n,2));
  ($M[1],$M[0]);
}
sub nth_stern_brocot {
  my($n) = @_;
  validate_integer_positive($n);

  my @M = (1,0);
  my @bits = split(//,Mtodigitstring($n,2));
  $M[$_] = Mvecsum(@M) for 1,reverse(@bits[1..$#bits]);
  ($M[1],$M[0]);
}

sub nth_stern_diatomic {
  my ($n) = @_;
  validate_integer_nonneg($n);
  my @M = (1,0);
  $M[$_] = Mvecsum(@M) for split(//, Mtodigitstring($n,2));
  $M[1];
}

sub farey {
  my($n,$k) = @_;
  validate_integer_positive($n);
  my $len = Maddint(Math::Prime::Util::sumtotient($n),1);

  my($p0, $q0, $p1, $q1, $p2, $q2, $j) = (0,1,1,$n);

  if (defined $k) {
    validate_integer_nonneg($k);
    return undef if $k >= $len;
    for (1 .. $k) {
      $j = Mdivint(($q0 + $n), $q1);
      $p2 = Mmulint($j, $p1) - $p0;
      $q2 = Mmulint($j, $q1) - $q0;
      ($p0, $q0, $p1, $q1) = ($p1, $q1, $p2, $q2);
    }
    return [$p0,$q0];
  }

  return $len unless wantarray;

  my @V;
  for (1 .. $len) {
    push @V, [$p0, $q0];
    $j = Mdivint(($q0 + $n), $q1);
    $p2 = Mmulint($j, $p1) - $p0;
    $q2 = Mmulint($j, $q1) - $q0;
    ($p0, $q0, $p1, $q1) = ($p1, $q1, $p2, $q2);
  }
  @V;
}

# Uses gcdext to find next entry with only one point.
sub next_farey {
  my($n,$frac) = @_;
  validate_integer_positive($n);
  croak "next_farey second argument not an array reference" unless ref($frac) eq 'ARRAY';
  my($p,$q) = @$frac;
  validate_integer_nonneg($p);
  validate_integer_positive($q);
  return undef if $p >= $q;
  my($u,$v,$g) = Mgcdext($p,$q);
  ($p,$q) = (Mdivint($p,$g),Mdivint($q,$g)) if $g != 1;
  my $d = Mmulint(Mdivint(($n+$u),$q),$q) - $u;
  my $c = Mdivint((Mmulint($d,$p)+1),$q);
  [$c,$d];
}

sub farey_rank {
  my($n,$frac) = @_;
  validate_integer_positive($n);
  croak "next_farey second argument not an array reference" unless ref($frac) eq 'ARRAY';
  my($p,$q) = @$frac;
  validate_integer_nonneg($p);
  validate_integer_positive($q);

  return 0 if $p == 0;

  my $g = Mgcd($p,$q);
  ($p,$q) = (Mdivint($p,$g),Mdivint($q,$g)) if $g != 1;

  my @count = (0,0,map { Mdivint(Mmulint($p,$_)-1,$q); } 2..$n);
  my $sum = 1;
  for my $i (2 .. $n) {
    my $icount = $count[$i];
    for (my $j = Mmulint($i,2); $j <= $n; $j = Maddint($j,$i)) {
      $count[$j] -= $icount;
    }
    $sum += $icount;
  }
  $sum;
}

# End of Rational maps


sub _is_perfect_square {
  my($n) = @_;
  return (1,1,0,0,1)[$n] if $n <= 4;

  if (ref($n)) {
    return 0 if ((1 << Mmodint($n,32)) & 0xfdfcfdec);
    my $sq = Msqrtint($n);
    return 1 if Mmulint($sq,$sq) == $n;
  } else {
    return 0 if (1 << ($n & 31)) & 0xfdfcfdec;
    my $sq = int(sqrt($n));
    return 1 if ($sq*$sq) == $n;
  }
  0;
}

sub is_primitive_root {
  my($a, $n) = @_;
  validate_integer($a);
  validate_integer_abs($n);

  return (undef,1)[$n] if $n <= 1;
  $a = Mmodint($a, $n) if $a < 0 || $a >= $n;
  return 0+($a == $n-1) if $n <= 4;
  return 0 if $a <= 1;

  return Math::Prime::Util::GMP::is_primitive_root($a,$n)
    if $Math::Prime::Util::_GMPfunc{"is_primitive_root"};

  # my $order = Mznorder($a,$n);  return 0 unless defined $order;  return 0+($order == Mtotient($n));

  if (Mis_even($n)) {
    return 0 if ($n % 4) == 0; # n can't still be even after we shift it
    return 0 if Mis_even($a);  # n and a cannot both be even
    $n = Mrshiftint($n,1);     # a is odd, so it is a primroot of p^k also
  }
  return 0 if Mgcd($a, $n) != 1;
  return 0 if _is_perfect_square($a);

  my ($p,$k,$phi);
  $k = Mis_prime_power($n,\$p);
  return 0 if !$k;
  $n = $p;
  $phi = Msubint($n,1);
  return 0 if $k > 1 && Mpowmod($a, $phi, Mmulint($p,$p)) == 1;

  return 0 if Mkronecker($a,$n) != -1;
  return 0 if ($phi % 3) == 0 && Mpowmod($a,Mdivint($phi,3),$n) == 1;
  return 0 if ($phi % 5) == 0 && Mpowmod($a,Mdivint($phi,5),$n) == 1;
  foreach my $f (Mfactor_exp($phi)) {
    my $fp = $f->[0];
    return 0 if $fp > 5 && Mpowmod($a, Mdivint($phi,$fp), $n) == 1;
  }
  1;
}

sub znorder {
  my($a, $n) = @_;
  validate_integer_abs($n);
  return (undef,1)[$n] if $n <= 1;
  $a = Mmodint($a, $n);
  return undef if $a <= 0;
  return 1 if $a == 1;

  return reftyped($_[0], Math::Prime::Util::GMP::znorder($a,$n))
    if $Math::Prime::Util::_GMPfunc{"znorder"};

  # Sadly, Calc/FastCalc are horrendously slow for this function.
  return undef if Mgcd($a, $n) > 1;

  # The answer is one of the divisors of phi(n) and lambda(n).
  my $lambda = Math::Prime::Util::carmichael_lambda($n);
  $a = tobigint($a);

  # This is easy and usually fast, but can bog down with too many divisors.
  if ($lambda <= 2**64) {
    foreach my $k (Mdivisors($lambda)) {
      return $k if Mpowmod($a,$k,$n) == 1;
    }
    return undef;
  }

  # Algorithm 1.7 from A. Das applied to Carmichael Lambda.
  my $k = 1;
  foreach my $f (Mfactor_exp($lambda)) {
    my($pi, $ei, $enum) = ($f->[0],$f->[1], 0);
    my $phidiv = Mdivint($lambda, Mpowint($pi,$ei));
    my $b = Mpowmod($a, $phidiv, $n);
    while ($b != 1) {
      return undef if $enum++ >= $ei;
      $b = Mpowmod($b, $pi, $n);
      $k = Mmulint($k, $pi);
    }
  }
  $k;
}

sub _dlp_trial {
  my ($a,$g,$p,$limit) = @_;
  $limit = $p if !defined $limit || $limit > $p;

  if ($limit < 1_000_000_000) {
    my $t = $g;
    for my $k (1 .. $limit) {
      return $k if $t == $a;
      $t = Mmulmod($t, $g, $p);
    }
    return 0;
  }

  ($a, $g, $p, $limit) = map { tobigint($_) } ($a, $g, $p, $limit);
  my $t = tobigint($g);
  for (my $k = tobigint(1); $k < $limit; $k++) {
    return Maddint($k,0) if $t == $a;
    $t *= $g;
    $t %= $p;
  }
  0;
}
sub _dlp_bsgs {
  my ($a,$g,$p,$_verbose) = @_;
  my $invg = Minvmod($g, $p);
  return 0 unless defined $invg;
  my $N = Maddint(Msqrtint($p-1),1);
  # Limit for time and space.
  my $b = $N > 4_000_000 ? 4_000_000 : $N;

  my %hash;
  my $am = 1;
  my $gm = Mpowmod($invg, $N, $p);
  my $key = $a;
  my $r;

  print "  BSGS starting $b loops\n" if $_verbose > 1;
  foreach my $m (0 .. $b) {
    # Baby Step
    if ($m <= $N) {
      $r = $hash{"$am"};
      if (defined $r) {
        print "  bsgs found in stage 1 after $m tries\n" if $_verbose;
        $r = Mmuladdmod($r, $N, $m, $p);
        return $r;
      }
      $hash{"$am"} = $m;
      $am = Mmulmod($am,$g,$p);
      if ($am == $a) {
        print "  bsgs found during bs\n" if $_verbose;
        return $m+1;
      }
    }

    # Giant Step
    $r = $hash{"$key"};
    if (defined $r) {
      print "  bsgs found in stage 2 after $m tries\n" if $_verbose;
      $r = Mmuladdmod($m, $N, $r, $p);
      return $r;
    }
    $hash{"$key"} = $m;
    $key = Mmulmod($key,$gm,$p);
  }
  0;
}

sub znlog {
  my($a, $g, $n) = @_;
  validate_integer($a);
  validate_integer($g);
  validate_integer_abs($n);
  return (undef,0,1)[$n] if $n <= 1;
  $a = Mmodint($a, $n);
  $g = Mmodint($g, $n);
  return 0 if $a == 1 || $g == 0 || $n < 2;

  my $_verbose = getconfig()->{'verbose'};

  # For large p, znorder can be very slow.  Do a small trial test first.
  my $x = _dlp_trial($a, $g, $n, 200);

  if ($x == 0) {
    ($a,$g,$n) = map { tobigint($_) } ($a,$g,$n);
    $x = _dlp_bsgs($a, $g, $n, $_verbose);
    $x = _bigint_to_int($x) if ref($x) && $x <= INTMAX;
    return $x if $x > 0 && Mpowmod($g,$x,$n) == $a;
    print "  BSGS giving up\n" if $x == 0 && $_verbose;
    print "  BSGS incorrect answer $x\n" if $x > 0 && $_verbose > 1;
    $x = _dlp_trial($a,$g,$n);
  }
  $x = _bigint_to_int($x) if ref($x) && $x <= INTMAX;
  return ($x == 0) ? undef : $x;
}

sub znprimroot {
  my($n) = @_;
  validate_integer_abs($n);
  return (undef,0,1,2,3)[$n] if $n <= 4;
  return if $n % 4 == 0;

  my $iseven = Mis_even($n);
  $n = Mrshiftint($n,1) if $iseven;

  my($k,$p);
  $k = Mis_prime_power($n, \$p);
  return if $k < 1;
  return 5 if $p == 3 && $iseven;
  my $ispow = ($k > 1);

  my $phi = $p-1;
  my $psquared = $ispow ? Mmulint($p,$p) : 0;

  my @phidivfac = map  { Mdivint($phi, $_) }
                  grep { $_ > 2 }
                  map  { $_->[0] }  Mfactor_exp($phi);
  my $a = 1;
  while (1) {
    $a += $iseven ? 2 : 1;
    return if $a >= $p;
    next if $a == 4 || $a == 8 || $a == 9;
    next if Mkronecker($a,$p) != -1;
    next if Mvecany(sub { Mpowmod($a,$_,$p) == 1 }, @phidivfac);
    return $a unless $ispow && Mpowmod($a,$phi,$psquared) == 1;
  }
}

sub qnr {
  my($n) = @_;
  validate_integer_abs($n);
  return (undef,1,2)[$n] if $n <= 2;

  return 2 if Mkronecker(2,$n) == -1;

  if (Mis_prime($n)) {
    for (my $a = 3; $a < $n; $a = Mnext_prime($a)) {
      return $a if Mkronecker($a,$n) == -1;
    }
  } else {
    if ($n % 2 == 0) {
      my $e = Mvaluation($n, 2);
      $n >>= $e;
      return 2 if $n == 1 || $e >= 2;
    }
    return 2 if !($n%3) || !($n%5) || !($n%11) || !($n%13) || !($n%19);
    my @F = Mfactor_exp($n);
    for (my $a = 2; $a < $n; $a = Mnext_prime($a)) {
      for my $pe (@F) {
        my $p = $pe->[0];
        return $a if $a < $p && Mkronecker($a,$p) == -1;
      }
    }
  }
  0;
}


# Find first D in sequence (5,-7,9,-11,13,-15,...) where (D|N) == -1
sub _lucas_selfridge_params {
  my($n) = @_;

  # D is typically quite small: 67 max for N < 10^19.  However, it is
  # theoretically possible D could grow unreasonably.  I'm giving up at 4000M.
  my $d = 5;
  my $sign = 1;
  while (1) {
    my $gcd = Mgcd($d, $n);
    return (0,0,0) if $gcd > 1 && $gcd != $n;  # Found divisor $d
    my $j = Mkronecker($d * $sign, $n);
    last if $j == -1;
    $d += 2;
    croak "Could not find Jacobi sequence for $n" if $d > 4_000_000_000;
    $sign = -$sign;
  }
  my $D = $sign * $d;
  my $P = 1;
  my $Q = int( (1 - $D) / 4 );
  ($P, $Q, $D)
}

sub _lucas_extrastrong_params {
  my($n, $increment) = @_;
  $increment = 1 unless defined $increment;

  my ($P, $Q, $D) = (3, 1, 5);
  while (1) {
    my $gcd = Mgcd($D, $n);
    return (0,0,0) if $gcd > 1 && $gcd != $n;  # Found divisor $d
    last if Mkronecker($D, $n) == -1;
    $P += $increment;
    croak "Could not find Jacobi sequence for $n" if $P > 65535;
    $D = $P*$P - 4;
  }
  ($P, $Q, $D);
}

# returns U_k, V_k, Q_k all mod n
sub lucas_sequence {
  my($n, $P, $Q, $k) = @_;

  croak "lucas_sequence: n must be > 0" if $n < 1;
  croak "lucas_sequence: k must be >= 0" if $k < 0;
  return (0,0,0) if $n == 1;

  if ($Math::Prime::Util::_GMPfunc{"lucas_sequence"} && $Math::Prime::Util::GMP::VERSION >= 0.30 && !ref($P) && !ref($Q)) {
    return maybetobigintall(
             Math::Prime::Util::GMP::lucas_sequence($n, $P, $Q, $k)
           );
  }

  return (lucasuvmod($P,$Q,$k,$n), Mpowmod($Q,$k,$n));
}

sub lucasuv {
  my($P, $Q, $k) = @_;

  croak "lucasuv: k must be >= 0" if $k < 0;
  return (0,2) if $k == 0;

  if ($Math::Prime::Util::_GMPfunc{"lucasuv"} && $Math::Prime::Util::GMP::VERSION >= 0.53) {
    return maybetobigintall(
             Math::Prime::Util::GMP::lucasuv($P, $Q, $k)
           );
  }

  # Do this very generic.  Optimize later if needed (D=0,Q=1,Q=-1,n odd).

  ($P,$Q) = map { tobigint($_) } ($P,$Q);
  my($Uh, $Vl, $Vh, $Ql, $Qh) = map { tobigint($_) } (1, 2, $P, 1, 1);

  my $s = 0;
  my @kbits = Mtodigits($k, 2);
  while ($kbits[-1] == 0) { $s++; pop @kbits; }  # Remove trailing zeros.
  pop @kbits;                                    # Remove trailing 1.

  foreach my $bit (@kbits) {
    $Ql *= $Qh;
    if ($bit) {
      $Qh = $Ql * $Q;
      $Uh = $Uh * $Vh;
      $Vl = $Vh * $Vl - $P * $Ql;
      $Vh = $Vh * $Vh - ($Qh+$Qh);
    } else {
      $Qh = $Ql;
      $Uh = $Uh * $Vl - $Ql;
      $Vh = $Vh * $Vl - $P * $Ql;
      $Vl = $Vl * $Vl - ($Ql+$Ql);
    }
  }
  $Ql *= $Qh;
  $Qh = $Ql * $Q;
  $Uh = $Uh * $Vl - $Ql;
  $Vl = $Vh * $Vl - $P * $Ql;
  $Ql *= $Qh;
  for (1 .. $s) {
    $Uh *= $Vl;
    $Vl = $Vl * $Vl - ($Ql+$Ql);
    $Ql *= $Ql;
  }
  $Uh = _bigint_to_int($Uh) if $Uh <= INTMAX && $Uh >= INTMIN;
  $Vl = _bigint_to_int($Vl) if $Vl <= INTMAX && $Vl >= INTMIN;
  ($Uh, $Vl);
}

sub lucasuvmod {
  my($P, $Q, $k, $n) = @_;
  validate_integer($P);
  validate_integer($Q);
  validate_integer_nonneg($k);
  validate_integer_abs($n);
  return if $n == 0;
  return (0,0) if $n == 1;
  return (0, Mmodint(2,$n)) if $k == 0;

  if ($Math::Prime::Util::_GMPfunc{"lucasuvmod"} && $Math::Prime::Util::GMP::VERSION >= 0.53) {
    return maybetobigintall(
             Math::Prime::Util::GMP::lucasuvmod($P, $Q, $k, $n)
           );
  }

  $P = Mmodint($P,$n) if $P < 0 || $P >= $n;
  $Q = Mmodint($Q,$n) if $Q < 0 || $Q >= $n;
  my $D = Mmulsubmod($P, $P, Mmulmod(4,$Q,$n), $n);

  if ($D == 0) {
    my $S = Mdivmod($P, 2, $n);
    if ($S) {
      my $U = Mmulmod($k, Mpowmod($S, $k-1, $n), $n);
      my $V = Mmulmod(2,  Mpowmod($S, $k,   $n), $n);
      return ($U, $V);
    }
  }

  my @kbits = Mtodigits($k, 2);
  shift @kbits;  # Remove leading 1
  my $U = 1;
  my $V = $P;
  my $invD = Minvmod($D, $n);
  my $nisodd = Mis_odd($n);

  if ($Q == 1 && $invD) {
    $U = Mmulsubmod($P, $P, 2, $n);
    foreach my $bit (@kbits) {
      my $T = Mmulsubmod($U, $V, $P, $n);
      if ($bit) {
        $V = $T;
        $U = Mmulsubmod($U, $U, 2, $n);
      } else {
        $U = $T;
        $V = Mmulsubmod($V, $V, 2, $n);
      }
    }
    $V = Mmodint($V,$n);
    $U = Maddmod($U, $U, $n);
    $U = Msubmod($U, Mmulmod($V, $P, $n), $n);
    $U = Mmulmod($U, $invD, $n);
  } elsif ($nisodd && ($Q == 1 || $Q == ($n-1))) {
    my $ps = ($P == 1);
    my $qs = ($Q == 1);
    my $halfn = Maddint(Mrshiftint($n,1),1);
    foreach my $bit (@kbits) {
      $U = Mmulmod($U, $V, $n);
      $V = ($qs) ? Mmulsubmod($V,$V,2,$n) : Mmuladdmod($V,$V,2,$n);
      $qs = 1;
      if ($bit) {
        my $t = Mmulmod($U, $D, $n);
        $U = (!$ps) ? Mmuladdmod($U,$P,$V,$n) : Maddmod($U,$V,$n);
        if (Mis_odd($U)) {
          $U = Maddint(Mrshiftint($U, 1), $halfn);
        } else {
          $U = Mrshiftint($U, 1);
        }
        $V = (!$ps) ? Mmuladdmod($V,$P,$t,$n) : Maddmod($V,$t,$n);
        if (Mis_odd($V)) {
          $V = Maddint(Mrshiftint($V, 1), $halfn);
        } else {
          $V = Mrshiftint($V, 1);
        }
        $qs = ($Q==1);
      }
    }
  } elsif ($nisodd) {
    my $Qk = $Q;
    my $halfn = Maddint(Mrshiftint($n,1),1);
    foreach my $bit (@kbits) {
      $U = Mmulmod($U, $V, $n);
      $V = Mmulsubmod($V, $V, Maddmod($Qk, $Qk, $n), $n);
      $Qk = Mmulmod($Qk, $Qk, $n);
      if ($bit) {
        my $t = Mmulmod($U, $D, $n);
        $U = Mmuladdmod($U, $P, $V, $n);
        if (Mis_odd($U)) {
          $U = Maddint(Mrshiftint($U, 1), $halfn);
        } else {
          $U = Mrshiftint($U, 1);
        }
        $V = Mmuladdmod($V, $P, $t, $n);
        if (Mis_odd($V)) {
          $V = Maddint(Mrshiftint($V, 1), $halfn);
        } else {
          $V = Mrshiftint($V, 1);
        }
        $Qk = Mmulmod($Qk, $Q, $n);
      }
    }
  } else {
    my ($s, $Uh, $Vl, $Vh, $Ql, $Qh) = (0, 1, 2, $P, 1, 1);
    unshift @kbits, 1;                             # Add back leading 1.
    while ($kbits[-1] == 0) { $s++; pop @kbits; }  # Remove trailing zeros.
    pop @kbits;                                    # Remove trailing 1.
    foreach my $bit (@kbits) {
      $Ql = Mmulmod($Ql, $Qh, $n);
      if ($bit) {
        $Qh = Mmulmod($Ql, $Q, $n);
        $Uh = Mmulmod($Uh, $Vh, $n);
        $Vl = Mmulsubmod($Vh, $Vl, Mmulmod($P,  $Ql, $n), $n);
        $Vh = Mmulsubmod($Vh, $Vh, Maddmod($Qh, $Qh, $n), $n);
      } else {
        $Qh = $Ql;
        $Uh = Mmulsubmod($Uh, $Vl, $Ql, $n);
        $Vh = Mmulsubmod($Vh, $Vl, Mmulmod($P,  $Ql, $n), $n);
        $Vl = Mmulsubmod($Vl, $Vl, Maddmod($Ql, $Ql, $n), $n);
      }
    }
    $Ql = Mmulmod($Ql, $Qh, $n);
    $Qh = Mmulmod($Ql, $Q, $n);
    $Uh = Mmulsubmod($Uh, $Vl, $Ql, $n);
    $Vl = Mmulsubmod($Vh, $Vl, Mmulmod($P, $Ql, $n), $n);
    $Ql = Mmulmod($Ql, $Qh, $n);
    for (1 .. $s) {
      $Uh = Mmulmod($Uh, $Vl, $n);
      $Vl = Mmulsubmod($Vl, $Vl, Maddmod($Ql, $Ql, $n), $n);
      $Ql = Mmulmod($Ql, $Ql, $n);
    }
    ($U, $V) = ($Uh, $Vl);
  }
  ($U,$V);
}

sub lucasu {
  return maybetobigint( Math::Prime::Util::GMP::lucasu($_[0], $_[1], $_[2]) )
    if $Math::Prime::Util::_GMPfunc{"lucasu"};
  (lucasuv(@_))[0];
}
sub lucasv {
  return maybetobigint( Math::Prime::Util::GMP::lucasv($_[0], $_[1], $_[2]) )
    if $Math::Prime::Util::_GMPfunc{"lucasv"};
  (lucasuv(@_))[1];
}

sub lucasumod {
  return maybetobigint( Math::Prime::Util::GMP::lucasumod($_[0], $_[1], $_[2], $_[3]) )
    if $Math::Prime::Util::_GMPfunc{"lucasumod"};
  (lucasuvmod(@_))[0];
}
sub lucasvmod {
  my($P, $Q, $k, $n) = @_;
  return maybetobigint( Math::Prime::Util::GMP::lucasvmod($P, $Q, $k, $n) )
    if $Math::Prime::Util::_GMPfunc{"lucasvmod"};
  validate_integer($P);
  validate_integer($Q);
  validate_integer_nonneg($k);
  validate_integer_abs($n);
  return if $n == 0;

  return (lucasuvmod($P, $Q, $k, $n))[1] if $Q != 1;

  # Fast algorithm for Q=1
  $P = Mmodint($P, $n);
  my $V = 2;
  my $U = $P;
  foreach my $bit (Mtodigits($k, 2)) {
    my $T = Mmulsubmod($U, $V, $P, $n);
    if ($bit) {
      $V = $T;
      $U = Mmulsubmod($U, $U, 2, $n);
    } else {
      $U = $T;
      $V = Mmulsubmod($V, $V, 2, $n);
    }
  }
  return $V;
}

my %_ppc = (3 => 8, 5 => 20, 7 => 16, 11 => 10, 13 => 28, 17 => 36, 19 => 18);
sub _pisano_pp {
  my($p,$e) = @_;
  return 1 if $e == 0;
  return 3 << ($e-1) if $p == 2 && $e < 32;
  return Mlshiftint(3,$e-1) if $p == 2;
  my $k = $_ppc{$p};

  if (!defined $k) {
    $k = Msubint($p, Mkronecker(5,$p));
    for my $f (Mfactor_exp($k)) {
      my($fac,$exp) = @$f;
      for my $j (1 .. $exp) {
        my $rk = Mdivint($k,$fac);
        last if Mlucasumod(1, $p-1, $rk, $p) != 0;
        $k = $rk;
      }
    }
    $_ppc{$p} = $k;
  }
  $k = Mmulint($k, Mpowint($p, $e-1)) if $e > 1;
  $k;
}
sub pisano_period {
  my($n) = @_;
  validate_integer_nonneg($n);
  return 0 if $n < 0;
  return (0,1,3,8,6,20,24,16,12,24,60)[$n] if $n <= 10;

  my $k = Mlcm(map { _pisano_pp($_->[0],$_->[1]) } Mfactor_exp($n));

  my $lim = Mmulint(6,$n);
  for (my $ret = $k;  $ret <= $lim;  $ret = Maddint($ret,$k)) {
    return $ret if Mlucasumod(1, -1, Msubint($ret,1), $n) == 1;
  }
  undef;
}

sub is_lucas_pseudoprime {
  my($n) = @_;

  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0 || _is_perfect_square($n);

  my ($P, $Q, $D) = _lucas_selfridge_params($n);
  return 0 if $D == 0;  # We found a divisor in the sequence
  die "Lucas parameter error: $D, $P, $Q\n" if ($D != $P*$P - 4*$Q);

  my($U, $V) = lucasuvmod($P, $Q, $n+1, $n);
  return ($U == 0) ? 1 : 0;
}

sub is_strong_lucas_pseudoprime {
  my($n) = @_;

  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0 || _is_perfect_square($n);

  my ($P, $Q, $D) = _lucas_selfridge_params($n);
  return 0 if $D == 0;  # We found a divisor in the sequence
  die "Lucas parameter error: $D, $P, $Q\n" if ($D != $P*$P - 4*$Q);

  my $m = $n+1;
  my($s, $k) = (0, $m);
  while ( $k > 0 && !($k % 2) ) {
    $s++;
    $k >>= 1;
  }
  my($U, $V) = lucasuvmod($P, $Q, $k, $n);
  return 1 if $U == 0;

  my $Qk = Mpowmod($Q,$k,$n);
  foreach my $r (0 .. $s-1) {
    return 1 if $V == 0;
    if ($r < ($s-1)) {
      $V = Mmulsubmod($V, $V, Maddmod($Qk,$Qk,$n), $n);
      $Qk = Mmulmod($Qk, $Qk, $n);
    }
  }
  return 0;
}

sub is_extra_strong_lucas_pseudoprime {
  my($n) = @_;

  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0 || _is_perfect_square($n);

  my ($P, $Q, $D) = _lucas_extrastrong_params($n);
  return 0 if $D == 0;  # We found a divisor in the sequence
  die "Lucas parameter error: $D, $P, $Q\n" if ($D != $P*$P - 4*$Q);

  # This would be a great place to use a factor remove function
  my($s, $k) = (0, Maddint($n,1));
  while (Mis_even($k) && $k != 0) {
    $s++;
    $k = Mrshiftint($k,1);
  }

  my($U, $V) = lucasuvmod($P, $Q, $k, $n);
  return 1 if $U == 0 && ($V == 2 || $V == Msubint($n,2));
  foreach my $r (0 .. $s-2) {
    return 1 if $V == 0;
    $V = Mmulsubmod($V, $V, 2, $n);
  }
  return 0;
}

sub is_almost_extra_strong_lucas_pseudoprime {
  my($n, $increment) = @_;
  $increment = 1 unless defined $increment;

  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0 || _is_perfect_square($n);

  my ($P, $Q, $D) = _lucas_extrastrong_params($n, $increment);
  return 0 if $D == 0;  # We found a divisor in the sequence
  die "Lucas parameter error: $D, $P, $Q\n" if ($D != $P*$P - 4*$Q);

  my($TWO, $V, $W) = map { tobigint($_) } (2, $P, $P*$P-2);
  $n = tobigint($n);

  my $kstr = todigitstring($n + 1, 2);
  $kstr =~ s/(0*)\z//;   # Remove trailing zeros
  my $s = length($1);
  my $bpos = 0;
  while (++$bpos < length($kstr)) {
    if (substr($kstr,$bpos,1)) {
      $V = Mmulsubmod($V, $W, $P,   $n);
      $W = Mmulsubmod($W, $W, $TWO, $n);
    } else {
      $W = Mmulsubmod($W, $V, $P,   $n);
      $V = Mmulsubmod($V, $V, $TWO, $n);
    }
  }
  return 1 if $V == 2 || $V == ($n-$TWO);
  foreach my $r (0 .. $s-2) {
    return 1 if $V == 0;
    $V = Mmulsubmod($V, $V, $TWO, $n);
  }
  0;
}

sub is_frobenius_khashin_pseudoprime {
  my($n) = @_;
  return 0+($n >= 2) if $n < 4;
  return 0 unless $n % 2;
  return 0 if _is_perfect_square($n);

  $n = tobigint($n);

  my($k,$c) = (2,1);
  if    ($n % 4 == 3) { $c = $n-1; }
  elsif ($n % 8 == 5) { $c = 2; }
  else {
    do {
      $c += 2;
      $k = Mkronecker($c, $n);
    } while $k == 1;
  }
  return 0 if $k == 0 || ($k == 2 && !($n % 3));

  my $ea = ($k == 2) ? 2 : 1;
  my($ra,$rb,$a,$b,$d) = ($ea,1,$ea,1,$n-1);
  while ($d != 0) {
    if ($d % 2 == 1) {
      ($ra, $rb) = ( (($ra*$a)%$n + ((($rb*$b)%$n)*$c)%$n) % $n,
                     (($rb*$a)%$n + ($ra*$b)%$n) % $n );
    }
    $d >>= 1;
    if ($d != 0) {
      ($a, $b) = ( (($a*$a)%$n + ((($b*$b)%$n)*$c)%$n) % $n,
                   (($b*$a)%$n + ($a*$b)%$n) % $n );
    }
  }
  return ($ra == $ea && $rb == $n-1) ? 1 : 0;
}

sub is_frobenius_underwood_pseudoprime {
  my($n) = @_;
  return 0+($n >= 2) if $n < 4;
  return 0 unless $n % 2;

  my($a, $temp1, $temp2);
  if ($n % 4 == 3) {
    $a = 0;
  } else {
    for ($a = 1; $a < 1000000; $a++) {
      next if $a==2 || $a==4 || $a==7 || $a==8 || $a==10 || $a==14 || $a==16 || $a==18;
      my $j = Mkronecker($a*$a - 4, $n);
      last if $j == -1;
      return 0 if $j == 0 || ($a == 20 && _is_perfect_square($n));
    }
  }
  $temp1 = Mgcd(($a+4)*(2*$a+5), $n);
  return 0 if $temp1 != 1 && $temp1 != $n;

  $n = tobigint($n);
  my($s, $t, $ap2) = map { tobigint($_) } (1, 2, $a+2);
  my $np1string = todigitstring($n+1,2);
  my $np1len = length($np1string);

  foreach my $bit (1 .. $np1len-1) {
    $temp2 = $t+$t;
    $temp2 += ($s * $a)  if $a != 0;
    $temp1 = $temp2 * $s;
    $temp2 = $t - $s;
    $s += $t;
    $t = ($s * $temp2) % $n;
    $s = $temp1 % $n;
    if ( substr( $np1string, $bit, 1 ) ) {
      if ($a == 0)  { $temp1 = $s + $s; }
      else          { $temp1 = $s * $ap2; }
      $temp1 += $t;
      $t = $t + $t - $s;
      $s = $temp1;
    }
  }
  $temp1 = (2*$a+5) % $n;
  return ($s == 0 && $t == $temp1) ? 1 : 0;
}

sub _perrin_signature {
  my($n) = @_;
  my @S = (1,$n-1,3, 3,0,2);
  return @S if $n <= 1;

  my @nbin = Mtodigits($n,2);
  shift @nbin;

  while (@nbin) {
    my @SUB = map { Maddmod($n-$S[5-$_], $n-$S[5-$_],$n) } 0..5;
    my @T = map { Mmuladdmod($S[$_], $S[$_], $SUB[$_], $n); } 0..5;

    my $T01 = Msubmod($T[2], $T[1], $n);
    my $T34 = Msubmod($T[5], $T[4], $n);
    my $T45 = Maddmod($T34, $T[3], $n);
    if (shift @nbin) {
      @S = ($T[0], $T01, $T[1], $T[4], $T45, $T[5]);
    } else {
      @S = ($T01, $T[1], Maddmod($T01,$T[0],$n), $T34, $T[4], $T45);
    }
  }
  @S;
}

sub is_perrin_pseudoprime {
  my($n, $restrict) = @_;
  _validate_integer($n);
  if (defined $restrict) { _validate_integer_nonneg($restrict); }
  else                   { $restrict = 0; }
  return 0+($n >= 2) if $n < 4;
  return 0 if $restrict > 2 && ($n % 2) == 0;

  my @S = _perrin_signature($n);
  return 0 unless $S[4] == 0;
  return 1 if $restrict == 0;
  return 0 unless $S[1] == Msubint($n,1);
  return 1 if $restrict == 1;
  my $j = Mkronecker(-23,$n);
  if ($j == -1) {
    my $B = $S[2];
    my $B2 = Mmulmod($B,$B,$n);
    my $A = Msubmod(Mmuladdmod(3, $B, 1, $n), $B2, $n);
    my $C = Mmulsubmod(3,$B2,2,$n);
    return 1 if $S[0] == $A && $S[2] == $B && $S[3] == $B && $S[5] == $C && $B != 3 && Mmulsubmod($B2,$B,$B,$n) == 1;
  } else {
    return 0 if $j == 0 && $n != 23 && $restrict > 2;
    return 1 if $S[0] == 1 && $S[2] == 3 && $S[3] == 3 && $S[5] == 2;
    return 1 if $S[0] == 0 && $S[5] == $n-1 && $S[2] != $S[3] && Maddmod($S[2],$S[3],$n) == $n-3 && Mmulmod(Msubmod($S[2],$S[3],$n),Msubmod($S[2],$S[3],$n),$n) == $n-(23%$n);
  }
  0;
}

# Aebi and Cairns (2008)
sub _catgamma {
  my($n,$mod) = @_;

  # Theorem 6, allowing us to possibly reduce n
  if ($mod < $n) {
    my $NP = Mdivint($n,$mod);
    if ($NP & 1) {  # odd
      return $mod*$NP == $n  ?  _catgamma($NP,$mod)  :  0;
    } else {
      return Mmulmod(_catgamma($NP+1,$mod),_catgamma($n-$mod*$NP,$mod),$mod);
    }
  }
  # Section 5 rephrases Theorem 2 into the middle binomial.
  my $N = Msubint($n,1);
  my $m = Mrshiftint($N,1);
  my $r = Math::Prime::Util::binomialmod($N, $m, $mod);
  return ($m & 1) ? $mod-$r : $r;
}
sub _catvtest {
  my($n,$p) = @_;
  while ($n = int($n/$p)) { return 1 if $n % 2; }
  0;
}
sub is_catalan_pseudoprime {
  my($n) = @_;
  return 0+($n >= 2) if $n < 4;
  return 0 unless $n & 1;

  {
    my @f = Mtrial_factor($n, 10000);
    if (@f == 2 && is_prime($f[1]) && $f[0] != $f[1]) {
      my($p,$q) = ($f[0],$f[1]);  # two primes, q > p
      return 0 if 2*$p+1 >= $q; # by Theorem 6(a)
      # Proposition 3 (semiprimes)
      return 0 unless _catgamma($q,$p) == 1 && _catgamma($p,$q) == 1;
    }
    if (is_prime($f[-1])) {  # fully factored
      for my $F (vecuniq(@f)) {
        return 0 if _catvtest($n-1,$F);
      }
    }
  }
  return _catgamma($n,$n) == 1 ? 1 : 0;
}

sub is_frobenius_pseudoprime {
  my($n, $P, $Q) = @_;
  ($P,$Q) = (0,0) unless defined $P && defined $Q;
  return 0+($n >= 2) if $n < 4;

  $n = tobigint($n);
  return 0 if Mis_even($n);

  my($k, $Vcomp, $D, $Du) = (0, 4);
  if ($P == 0 && $Q == 0) {
    ($P,$Q) = (-1,2);
    while ($k != -1) {
      $P += 2;
      $P = 5 if $P == 3;  # Skip 3
      $D = $P*$P-4*$Q;
      $Du = ($D >= 0) ? $D : -$D;
      last if $P >= $n || $Du >= $n;   # TODO: remove?
      $k = Mkronecker($D, $n);
      return 0 if $k == 0;
      return 0 if $P == 10001 && _is_perfect_square($n);
    }
  } else {
    $D = $P*$P-4*$Q;
    $Du = ($D >= 0) ? $D : -$D;
    croak "Frobenius invalid P,Q: ($P,$Q)" if _is_perfect_square($Du);
  }
  return (Mis_prime($n) ? 1 : 0) if $n <= $Du || $n <= abs($Q) || $n <= abs($P);
  return 0 if Mgcd(abs($P*$Q*$D), $n) > 1;

  if ($k == 0) {
    $k = Mkronecker($D, $n);
    return 0 if $k == 0;
    my $Q2 = (2*abs($Q)) % $n;
    $Vcomp = ($k == 1) ? 2 : ($Q >= 0) ? $Q2 : $n-$Q2;
  }

  my($U, $V) = lucasuvmod($P, $Q, $n-$k, $n);
  return 1 if $U == 0 && $V == $Vcomp;
  0;
}

# Since people have graciously donated millions of CPU years to doing these
# tests, it would be rude of us not to use the results.  This means we don't
# actually use the pretest and Lucas-Lehmer test coded below for any reasonable
# size number.
# See: http://www.mersenne.org/report_milestones/
my %_mersenne_primes;
undef @_mersenne_primes{2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841};

sub is_mersenne_prime {
  my $p = shift;

  # Use the known Mersenne primes
  return 1 if exists $_mersenne_primes{$p};
  return 0 if $p < 79711549; # GIMPS has tested and verified all below
  # Past this we do a generic Mersenne prime test

  return 1 if $p == 2;
  return 0 unless is_prob_prime($p);
  return 0 if $p > 3 && $p % 4 == 3 && $p < ((~0)>>1) && is_prob_prime($p*2+1);
  my $mp = Msubint(Mlshiftint(1,$p), 1);

  # Definitely faster than using Math::BigInt that doesn't have GMP.
  return (0 == (Math::Prime::Util::GMP::lucasuvmod(4, 1, $mp+1, $mp))[0])
    if $Math::Prime::Util::_GMPfunc{"lucasuvmod"};

  my $V = 4;
  for my $k (3 .. $p) {
    $V = Mmulsubmod($V, $V, 2, $mp);
  }
  return $V == 0;
}


sub _poly_new {
  my($refn, @poly) = @_;
  push @poly, 0 unless scalar @poly;
  @poly = map { tobigint("$_") } @poly if $refn;
  return \@poly;
}

#sub _poly_print {
#  my($poly) = @_;
#  carp "poly has null top degree" if $#$poly > 0 && !$poly->[-1];
#  foreach my $d (reverse 1 .. $#$poly) {
#    my $coef = $poly->[$d];
#    print "", ($coef != 1) ? $coef : "", ($d > 1) ? "x^$d" : "x", " + "
#      if $coef;
#  }
#  my $p0 = $poly->[0] || 0;
#  print "$p0\n";
#}

sub _poly_mod_mul {
  my($px, $py, $r, $n) = @_;

  my $px_degree = $#$px;
  my $py_degree = $#$py;
  my @res = ref($n)  ?  map { tobigint(0) } 0..$r-1  :  map { 0 } 0..$r-1;

  # convolve(px, py) mod (X^r-1,n)
  my @indices_y = grep { $py->[$_] } (0 .. $py_degree);
  foreach my $ix (0 .. $px_degree) {
    my $px_at_ix = $px->[$ix];
    next unless $px_at_ix;
    foreach my $iy (@indices_y) {
      my $rindex = ($ix + $iy) % $r;  # reduce mod X^r-1
      $res[$rindex] = ($res[$rindex] + $px_at_ix * $py->[$iy]) % $n;
    }
  }
  # In case we had upper terms go to zero after modulo, reduce the degree.
  pop @res while !$res[-1];
  return \@res;
}

sub _poly_mod_pow {
  my($pn, $power, $r, $mod) = @_;
  my $res = _poly_new(ref($mod), 1);
  my $p = $power;

  while ($p) {
    $res = _poly_mod_mul($res, $pn, $r, $mod) if ($p % 2) != 0;
    $p >>= 1;
    $pn  = _poly_mod_mul($pn,  $pn, $r, $mod) if $p;
  }
  return $res;
}

sub _test_anr {
  my($a, $n, $r) = @_;
  my $pp = _poly_mod_pow(_poly_new(ref($n), $a, 1), $n, $r, $n);
  my $nr = $n % $r;
  $pp->[$nr] = (($pp->[$nr] || 0) -  1) % $n;  # subtract X^(n%r)
  $pp->[  0] = (($pp->[  0] || 0) - $a) % $n;  # subtract a
  return 0 if scalar grep { $_ } @$pp;
  1;
}

sub _log_gamma {
  my $x = shift;
  my @lanczos = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
      771.32342877765313, -176.61502916214059, 12.507343278686905,
      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7);
  my($base,$sum) = ($x+7.5, 0);
  $sum += $lanczos[$_] / ($x + $_)  for (8,7,6,5,4,3,2,1);
  $sum += $lanczos[0];
  return 0.91893853320467274178 + log($sum/$x) + (($x+0.5)*log($base)-$base);
}
sub _log_binomial {
  my($n,$k) = @_;
  return 0 if $n < $k;
  return _log_gamma($n+1) - _log_gamma($k+1) - _log_gamma($n-$k+1);
}
sub _log_bern41_binomial {
  my($r,$d,$i,$j,$s) = @_;
  return _log_binomial( 2*$s,    $i)
       + _log_binomial( $d,      $i)
       + _log_binomial( 2*$s-$i, $j)
       + _log_binomial( $r-2-$d, $j);
}
sub _bern41_acceptable {
  my($n,$r,$s) = @_;
  my $scmp = int(sqrt(($r-1)/3.0) + 0.99999) * log($n);
  my $d = int(0.5 * ($r-1));
  my $i = int(0.475 * ($r-1));
  my $j = $i;
  $d = $r-2 if $d > $r-2;
  $i = $d if $i > $d;
  $j = $r-2-$d if $j > ($r-2-$d);
  return _log_bern41_binomial($r,$d,$i,$j,$s) >= $scmp;
}

sub is_aks_prime {
  my $n = shift;
  _validate_integer($n);
  return 0 if $n < 2 || Mis_power($n);
  return 1 if $n == 2;

  if ($n > 11) {
    return 0 if Mis_divisible($n,2,3,5,7,11);
  }

  my($starta, $s);
  my $_verbose = getconfig()->{'verbose'};

  my $log2n = log($n)/log(2) + 0.0001;      # Error on large side.
  my $r0    = ($log2n > 32 ? 0.010 : 0.003) * $log2n * $log2n;
  my $rmult =  $log2n > 32 ? 6     : 30;

  my $r = Mnext_prime($r0 < 2 ? 2 : Mtoint($r0));
  while (   !Math::Prime::Util::is_primitive_root($n,$r)
         || !_bern41_acceptable($n,$r,$rmult * ($r-1))) {
    $r = next_prime($r);
  }

  {
    my $bi = 1;
    my $bj = $rmult * ($r-1);
    while ($bi < $bj) {
      $s = $bi + (($bj-$bi) >> 1);
      if (!_bern41_acceptable($n, $r, $s)) { $bi = $s+1; }
      else                                 { $bj = $s; }
    }
    $s = $bj;
    croak "AKS: internal error bad s" unless _bern41_acceptable($n, $r, $s);
    # S will range from 2 to s+1
    $starta = 2;
    $s = $s+1;
  }
  my $slim = $s * ($s-1);
  print "# aks trial to $slim\n" if $_verbose >= 2;
  {
    my @f = Mtrial_factor($n, $slim);
    return 0 if @f >= 2;
  }
  return 1 if Mmulint($slim,$slim) >= $n;
  # Check b^(n-1) = 1 mod n for b in [2..s]
  for my $a (2 .. $s) {
    return 0 if Mpowmod($a, $n-1, $n) != 1;
  }

  if ($n < (MPU_HALFWORD-1) ) {
    $n = _bigint_to_int($n) if ref($n);
  } else {
    $n = tobigint($n);
  }

  print "# aks r = $r  s = $s\n" if $_verbose;
  local $| = 1 if $_verbose > 1;
  for (my $a = $starta; $a <= $s; $a++) {
    return 0 unless _test_anr($a, $n, $r);
    print "." if $_verbose > 1;
  }
  print "\n" if $_verbose > 1;

  return 1;
}


################################################################################

sub factor_exp {
  my($n) = @_;
  _validate_integer_nonneg($n);

  my %exponents;
  my @factors = grep { !$exponents{$_}++ } Mfactor($n);
  return scalar @factors unless wantarray;
  return (map { [$_, $exponents{$_}] } @factors);
}

sub _basic_factor {
  # MODIFIES INPUT SCALAR
  return ($_[0] == 1) ? () : ($_[0])   if $_[0] < 4;

  my @factors;
  if (!ref($_[0])) {
    while ( !($_[0] % 2) ) { push @factors, 2;  $_[0] = int($_[0] / 2); }
    while ( !($_[0] % 3) ) { push @factors, 3;  $_[0] = int($_[0] / 3); }
    while ( !($_[0] % 5) ) { push @factors, 5;  $_[0] = int($_[0] / 5); }
  } else {
    if (Mgcd($_[0], 30) != 1) {
      while ($_[0] % 2 == 0) { push @factors, 2;  $_[0] >>= 1; }
      while ($_[0] % 3 == 0) { push @factors, 3;  $_[0] = Mdivint($_[0],3); }
      while ($_[0] % 5 == 0) { push @factors, 5;  $_[0] = Mdivint($_[0],5); }
    }
  }

  if ($_[0] > 1 && _is_prime7($_[0])) {
    push @factors, $_[0];
    $_[0] = 1;
  }
  @factors;
}

sub trial_factor {
  my($n, $limit) = @_;
  _validate_integer_nonneg($n);
  _validate_integer_nonneg($limit) if defined $limit;

  return ($n==1) ? () : ($n)  if $n < 4;

  if ($Math::Prime::Util::_GMPfunc{"trial_factor"}) {
    my @f = defined $limit ? Math::Prime::Util::GMP::trial_factor($n,$limit)
                           : Math::Prime::Util::GMP::trial_factor($n);
    return ref($_[0]) ? maybetobigintall(@f) : @f;
  }

  my @factors;
  # Don't use _basic_factor here -- they want a trial forced.

  # For 32-bit n, we can simplify things a lot.
  if ($n <= 4294967295) {
    my $sqrtn = int(sqrt($n));
    $limit = $sqrtn if !defined $limit || $limit > $sqrtn;

    if ($limit >= 2 && ($n % 2 == 0)) {
      do { push @factors, 2;  $n >>= 1; } while ($n % 2) == 0;
      $sqrtn = int(sqrt($n));
      $limit = $sqrtn if $sqrtn < $limit;
    }
    for my $p (3,5,7,11,13,17,19,23,29,31,37,41,43,47,53) {
      last if $n == 1 || $p > $limit;
      if ($n % $p == 0) {
        do { push @factors, $p;  $n = int($n/$p); } while ($n % $p) == 0;
        $sqrtn = int(sqrt($n));
        $limit = $sqrtn if $sqrtn < $limit;
      }
    }
    return @factors if $n == 1;
    return (@factors,$n) if $limit < 59;

    push @_primes_small, @{Mprimes($_primes_small[-1]+1, $limit+72)}
      if $limit > $_primes_small[-1];

    for my $i (17 .. $#_primes_small) {
      my $p = $_primes_small[$i];
      last if $p > $limit;
      if (($n % $p) == 0) {
        do { push @factors, $p;  $n = int($n/$p); } while ($n % $p) == 0;
        last if $n == 1;
        $sqrtn = int(sqrt($n));
        $limit = $sqrtn if $sqrtn < $limit;
      }
    }
    push @factors, $n if $n > 1;
    return @factors;
  }

  my $start_idx = 1;
  # Expand small primes if it would help.
  push @_primes_small, @{Mprimes($_primes_small[-1]+1, 100_003)}
    if $n > 400_000_000
    && $_primes_small[-1] < 99_000
    && (!defined $limit || $limit > $_primes_small[-1]);

  my $sqrtn = Msqrtint($n);
  $limit = $sqrtn if !defined $limit || $limit > $sqrtn;

  # Do initial bigint reduction.  Hopefully reducing it to native int.
  if (ref($n)) {
    my $ismbi = ref($n) eq 'Math::BigInt';
    while ($start_idx <= $#_primes_small) {
      # Math::BigInt is *terribly* slow doing mods.  Use GCDs => 2-3x faster.
      if ($ismbi && $start_idx <= $#_primes_small-2 && $_primes_small[$start_idx+2] <= $limit && $_primes_small[$start_idx] <= 99989) {
        my $g = $_primes_small[$start_idx+0] * $_primes_small[$start_idx+1] * $_primes_small[$start_idx+2];
        if ($n->bgcd($g)->is_one) {
          $start_idx += 3;
          next;
        }
      }
      my $f = $_primes_small[$start_idx++];
      last if $f > $limit;
      if ($n % $f == 0) {
        do {
          push @factors, $f;
          $n = Mdivint($n,$f);
        } while $n % $f == 0;
        last if $n < INTMAX;
        $sqrtn = Msqrtint($n);
        $limit = $sqrtn if $limit > $sqrtn;
      }
    }
    return @factors if $n == 1;
    return (@factors,$n) if $start_idx <= $#_primes_small && $_primes_small[$start_idx] > $limit;
  }

  if (!ref($n)) {
    for my $i ($start_idx .. $#_primes_small) {
      my $p = $_primes_small[$i];
      last if $p > $limit;
      if (($n % $p) == 0) {
        do { push @factors, $p;  $n = int($n/$p); } while ($n % $p) == 0;
        last if $n == 1;
        my $newlim = int( sqrt($n) + 0.001);
        $limit = $newlim if $newlim < $limit;
      }
    }
    if ($_primes_small[-1] < $limit) {
      my $inc = (($_primes_small[-1] % 6) == 1) ? 4 : 2;
      my $p = $_primes_small[-1] + $inc;
      while ($p <= $limit) {
        if (($n % $p) == 0) {
          do { push @factors, $p;  $n = int($n/$p); } while ($n % $p) == 0;
          last if $n == 1;
          my $newlim = int( sqrt($n) + 0.001);
          $limit = $newlim if $newlim < $limit;
        }
        $p += ($inc ^= 6);
      }
    }
  } else {   # n is a bigint.  Use mod-210 wheel trial division.
    # Generating a wheel mod $w starting at $s:
    # mpu 'my($s,$w,$t)=(11,2*3*5); say join ",",map { ($t,$s)=($_-$s,$_); $t; } grep { gcd($_,$w)==1 } $s+1..$s+$w;'
    # Should start at $_primes_small[$start_idx], do 11 + next multiple of 210.

    my @incs = (2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4,6,8,4,2,4,2,4,8,6,4,6,2,4,6,2,6,6,4,2,4,6,2,6,4,2,4,2,10,2,10);
    my $f = 11; while ($f <= $_primes_small[$start_idx-1]-210) { $f += 210; }
    SEARCH: while ($f <= $limit) {
      foreach my $finc (@incs) {
        if ($n % $f == 0 && $f <= $limit) {
          do {
            push @factors, $f;
            $n = Mdivint($n,$f);
          } while $n % $f == 0;
          last SEARCH if $n == 1;
          $sqrtn = Msqrtint($n);
          $limit = $sqrtn if $limit > $sqrtn;
        }
        $f += $finc;
      }
    }
  }
  push @factors, $n  if $n > 1;
  @factors;
}

my $_holf_r;
my @_fsublist = (
  [ "power",      sub { _power_factor (shift) } ],
  [ "pbrent 32k", sub { pbrent_factor (shift,   32*1024, 1, 1) } ],
  [ "p-1 1M",     sub { pminus1_factor(shift, 1_000_000, undef, 1); } ],
  [ "ECM 1k",     sub { ecm_factor    (shift,     1_000,   5_000, 15) } ],
  [ "pbrent 512k",sub { pbrent_factor (shift,  512*1024, 7, 1) } ],
  [ "p-1 4M",     sub { pminus1_factor(shift, 4_000_000, undef, 1); } ],
  [ "ECM 10k",    sub { ecm_factor    (shift,    10_000,  50_000, 10) } ],
  [ "pbrent 512k",sub { pbrent_factor (shift,  512*1024, 11, 1) } ],
  [ "HOLF 256k",  sub { holf_factor   (shift, 256*1024, $_holf_r); $_holf_r += 256*1024; } ],
  [ "p-1 20M",    sub { pminus1_factor(shift,20_000_000); } ],
  [ "ECM 100k",   sub { ecm_factor    (shift,   100_000, 800_000, 10) } ],
  [ "HOLF 512k",  sub { holf_factor   (shift, 512*1024, $_holf_r); $_holf_r += 512*1024; } ],
  [ "pbrent 2M",  sub { pbrent_factor (shift, 2048*1024, 13, 1) } ],
  [ "HOLF 2M",    sub { holf_factor   (shift, 2048*1024, $_holf_r); $_holf_r += 2048*1024; } ],
  [ "ECM 1M",     sub { ecm_factor    (shift, 1_000_000, 1_000_000, 10) } ],
  [ "p-1 100M",   sub { pminus1_factor(shift, 100_000_000, 500_000_000); } ],
);

sub factor {
  my($n) = @_;
  _validate_integer_nonneg($n);

  my @factors;
  if ($n < 4) {
    @factors = ($n == 1) ? () : ($n);
    return @factors;
  }

  if ($Math::Prime::Util::_GMPfunc{"factor"}) {
    my @factors = Math::Prime::Util::GMP::factor($n);
    return ref($_[0]) ? maybetobigintall(@factors) : @factors;
  }

  $n = Maddint($n,0) if ref($n);  # Ensure we have a copy
  my $lim = 4999;  # How much trial factoring to do

  # For native integers, we could save a little time by doing hardcoded trials
  # by 2-29 here.  Skipping it.

  push @factors, Mtrial_factor($n, $lim);
  return @factors if $factors[-1] < $lim*$lim;
  $n = pop(@factors);

  my @nstack = ($n);
  while (@nstack) {
    $n = pop @nstack;
    # Don't use bignum on $n if it has gotten small enough.
    $n = _bigint_to_int($n) if ref($n) && $n <= INTMAX;
    #print "Looking at $n with stack ", join(",",@nstack), "\n";
    while ( ($n >= ($lim*$lim)) && !_is_prime7($n) ) {
      my @ftry;
      $_holf_r = 1;
      foreach my $sub (@_fsublist) {
        last if scalar @ftry >= 2;
        print "  starting $sub->[0]\n" if getconfig()->{'verbose'} > 1;
        @ftry = $sub->[1]->($n);
      }
      if (scalar @ftry > 1) {
        #print "  split into ", join(",",@ftry), "\n";
        $n = shift @ftry;
        $n = _bigint_to_int($n) if ref($n) && $n <= INTMAX;
        push @nstack, @ftry;
      } else {
        #warn "trial factor $n\n";
        push @factors, Mtrial_factor($n);
        #print "  trial into ", join(",",@factors), "\n";
        $n = 1;
        last;
      }
    }
    push @factors, $n  if $n != 1;
  }
  Mvecsort(@factors);
}

sub _found_factor {
  my($f, $n, $what, @factors) = @_;
  if ($f == 1 || $f == $n) {
    push @factors, $n;
  } else {
    my $f2 = Mdivint($n,$f);
    croak "internal error in $what" unless Mmulint($f,$f2) == $n;
    ($f,$f2) = ($f2,$f) if $f > $f2;
    push @factors, $f, $f2;
    # MPU::GMP prints this type of message if verbose, so do the same.
    print "$what found factor $f\n" if getconfig()->{'verbose'} > 0;
  }
  @factors;
}

################################################################################

# TODO:
sub squfof_factor { Mtrial_factor(@_) }
sub lehman_factor { Mtrial_factor(@_) }
sub pplus1_factor { pminus1_factor(@_) }

sub _power_factor {
  my $r;
  my $k = Mis_power($_[0],0,\$r);
  return ($_[0]) unless $k > 1;
  print "power found factor $r\n" if getconfig()->{'verbose'} > 0;
  map { $r } 1..$k;
}

sub prho_factor {
  my($n, $rounds, $pa, $skipbasic) = @_;
  validate_integer_nonneg($n);
  if (defined $rounds) { validate_integer_nonneg($rounds); }
  else                 { $rounds = 4*1024*1024; }
  if (defined $pa)     { validate_integer_nonneg($pa); }
  else                 { $pa = 3; }

  my @factors;
  if (!$skipbasic) {
    @factors = _basic_factor($n);
    return @factors if $n < 4;
  }

  my($U,$V) = (7,7);

  if (ref($n) || $n >= MPU_HALFWORD) {

    my $inner = 32;
    $rounds = int( ($rounds + $inner-1) / $inner );
    while ($rounds-- > 0) {
      my($m, $oldU, $oldV, $f) = (1, $U, $V);
      for my $i (1 .. $inner) {
        $U = Mmuladdmod($U, $U, $pa, $n);
        $V = Mmuladdmod($V, $V, $pa, $n);
        $V = Mmuladdmod($V, $V, $pa, $n);
        $f = ($U > $V) ? Msubint($U,$V) : Msubint($V,$U);
        $m = Mmulmod($m,$f,$n);
      }
      $f = Mgcd($m,$n);
      next if $f == 1;
      if ($f == $n) {
        ($U, $V) = ($oldU, $oldV);
        for my $i (1 .. $inner) {
          $U = Mmuladdmod($U, $U, $pa, $n);
          $V = Mmuladdmod($V, $V, $pa, $n);
          $V = Mmuladdmod($V, $V, $pa, $n);
          $f = ($U > $V) ? Msubint($U,$V) : Msubint($V,$U);
          $f = Mgcd($f, $n);
          last if $f != 1;
        }
        last if $f == 1 || $f == $n;
      }
      return _found_factor($f, $n, "prho-bigint", @factors);
    }

  } else {

    my $inner = 32;
    $rounds = int( ($rounds + $inner-1) / $inner );
    while ($rounds-- > 0) {
      my($m, $oldU, $oldV, $f) = (1, $U, $V);
      for my $i (1 .. $inner) {
        $U = ($U * $U + $pa) % $n;
        $V = ($V * $V + $pa) % $n;
        $V = ($V * $V + $pa) % $n;
        $f = ($U > $V) ? $U-$V : $V-$U;
        $m = ($m * $f) % $n;
      }
      $f = _gcd_ui( $m, $n );
      next if $f == 1;
      if ($f == $n) {
        ($U, $V) = ($oldU, $oldV);
        for my $i (1 .. $inner) {
          $U = ($U * $U + $pa) % $n;
          $V = ($V * $V + $pa) % $n;
          $V = ($V * $V + $pa) % $n;
          $f = ($U > $V) ? $U-$V : $V-$U;
          $f = _gcd_ui( $f, $n);
          last if $f != 1;
        }
        last if $f == 1 || $f == $n;
      }
      return _found_factor($f, $n, "prho-32", @factors);
    }

  }
  push @factors, $n;
  @factors;
}

sub pbrent_factor {
  my($n, $rounds, $pa, $skipbasic) = @_;
  validate_integer_nonneg($n);
  if (defined $rounds) { validate_integer_nonneg($rounds); }
  else                 { $rounds = 4*1024*1024; }
  if (defined $pa)     { validate_integer_nonneg($pa); }
  else                 { $pa = 3; }

  my @factors;
  if (!$skipbasic) {
    @factors = _basic_factor($n);
    return @factors if $n < 4;
  }

  my($Xi,$Xm) = (2,2);

  if (ref($n) || $n >= MPU_HALFWORD) {

    # Same code as the GMP version, but runs *much* slower.  Even with
    # Math::BigInt::GMP it's >200x slower.  With the default Calc backend
    # it's thousands of times slower.
    my($inner,$r,$saveXi,$f) = (32,1);

    while ($rounds > 0) {
      my $rleft = ($r > $rounds) ? $rounds : $r;
      while ($rleft > 0) {
        my $dorounds = ($rleft > $inner) ? $inner : $rleft;
        my $m = 1;
        $saveXi = Maddint($Xi,0);
        foreach my $i (1 .. $dorounds) {
          $Xi = Mmuladdmod($Xi, $Xi, $pa, $n);
          $m = Mmulmod($m, $Xi > $Xm ? $Xi-$Xm : $Xm-$Xi,$n);
        }
        $rleft -= $dorounds;
        $rounds -= $dorounds;
        $f = Mgcd($m,$n);
        last unless $f == 1;
      }
      if ($f == 1) {
        $r *= 2;
        $Xm = Maddint($Xi,0);
        next;
      }
      if ($f == $n) {  # back up to determine the factor
        $Xi = Maddint($saveXi,0);
        do {
          $Xi = Mmuladdmod($Xi, $Xi, $pa, $n);
          $f = Mgcd($Xi > $Xm ? $Xi-$Xm : $Xm-$Xi, $n);
        } while ($f != 1 && $r-- != 0);
        last if $f == 1 || $f == $n;
      }
      return _found_factor($f, $n, "pbrent", @factors);
    }

  } else {

    # Doing the gcd batching as above works pretty well here, but it's a lot
    # of code for not much gain for general users.
    for my $i (1 .. $rounds) {
      $Xi = ($Xi * $Xi) % $n;
      $Xi += $pa; $Xi -= $n if $Xi >= $n;
      my $f = _gcd_ui( ($Xi>$Xm) ? $Xi-$Xm : $Xm-$Xi, $n);
      return _found_factor($f, $n, "pbrent-32",@factors) if $f != 1 && $f != $n;
      $Xm = $Xi if ($i & ($i-1)) == 0;  # i is a power of 2
    }

  }
  push @factors, $n;
  @factors;
}

sub pminus1_factor {
  my($n, $B1, $B2, $skipbasic) = @_;
  validate_integer_nonneg($n);
  validate_integer_nonneg($B1) if defined $B1;
  validate_integer_nonneg($B2) if defined $B2;

  my @factors;
  if (!$skipbasic) {
    @factors = _basic_factor($n);
    return @factors if $n < 4;
  }

  if (!ref($n)) {
    # Stage 1 only
    my $sqrtn = Msqrtint($n);
    $B1 = !defined $B1 || $B1 > $sqrtn ? $sqrtn : $sqrtn;
    my $sqrtb1 = int(sqrt($B1));
    my($pc_beg, $pc_end) = (2, 6_000-1);
    my $pa = 2;

    while (1) {
      $pc_end = $B1 if $pc_end > $B1;
      foreach my $q (@{Mprimes($pc_beg, $pc_end)}) {
        my $k = $q;
        if ($q <= $sqrtb1) {
          my $kmin = int($B1 / $q);
          while ($k <= $kmin) { $k *= $q; }
        }
        $pa = Mpowmod($pa, $k, $n);
        if ($pa == 0) { push @factors, $n; return @factors; }
        my $f = Mgcd($pa-1, $n);
        return _found_factor($f, $n, "pminus1-64", @factors) if $f != 1;
      }
      last if $pc_end >= $B1;
      ($pc_beg, $pc_end) = ($pc_end+1, $pc_end+18000);
    }
    push @factors, $n;
    return @factors;
  }

  if (!defined $B1) {
    for my $mul (1, 100, 1000, 10_000, 100_000, 1_000_000) {
      $B1 = 1000 * $mul;
      $B2 = 1*$B1;
      #warn "Trying p-1 with $B1 / $B2\n";
      my @nf = pminus1_factor($n, $B1, $B2);
      if (scalar @nf > 1) {
        push @factors, @nf;
        return @factors;
      }
    }
    push @factors, $n;
    return @factors;
  }
  $B2 = 1*$B1 unless defined $B2;

  $n = tobigint($n) if !ref($n) || (defined $_BIGINT && $_BIGINT ne ref($n));
  # bigints:  n, pa, t, savea, [stage2] b, bm

  my ($j, $q, $saveq) = (32, 2, 2);
  my $pa = tobigint(2);
  my $t  = tobigint(1);
  my $savea = $pa+0;
  my $f = 1;
  my($pc_beg, $pc_end) = (2, 2+100_000);

  while (1) {
    $pc_end = $B1 if $pc_end > $B1;
    my @bprimes = @{ Mprimes($pc_beg, $pc_end) };
    foreach my $q (@bprimes) {
      my($k, $kmin) = ($q, int($B1 / $q));
      while ($k <= $kmin) { $k *= $q; }
      $t *= $k;                         # accumulate powers for a
      if ( ($j++ % 64) == 0) {
        next if $pc_beg > 2 && ($j-1) % 256;
        $pa = _bi_powmod($pa, $t, $n);
        $t = tobigint(1);
        if ($pa == 0) { push @factors, $n; return @factors; }
        $f = Mgcd($pa-1, $n);
        last if $f == $n;
        return _found_factor($f, $n, "pminus1-bigint $B1", @factors) unless $f == 1;
        $saveq = $q;
        $savea = $pa+0;
      }
    }
    $q = $bprimes[-1];
    last if $f != 1 || $pc_end >= $B1;
    ($pc_beg, $pc_end) = (Maddint($pc_end,1), Maddint($pc_end,500_000));
  }
  $pa = _bi_powmod($pa, $t, $n);
  if ($pa == 0) { push @factors, $n; return @factors; }
  $f = Mgcd($pa-1, $n);
  if ($f == $n) {
    $q = $saveq;
    $pa = $savea+0;
    while ($q <= $B1) {
      my ($k, $kmin) = ($q, int($B1 / $q));
      while ($k <= $kmin) { $k *= $q; }
      $pa = _bi_powmod($pa, $k, $n);
      my $f = Mgcd($pa-1, $n);
      if ($f == $n) { push @factors, $n; return @factors; }
      last if $f != 1;
      $q = Mnext_prime($q);
    }
  }
  # STAGE 2
  if ($f == 1 && $B2 > $B1) {
    my $bm = $pa + 0;
    my $b = tobigint(1);
    my @precomp_bm;
    $precomp_bm[0] = ($bm * $bm) % $n;
    $precomp_bm[$_] = ($precomp_bm[$_-1] * $bm * $bm) % $n for 1..19;
    $pa = _bi_powmod($pa, $q, $n);

    my $j = 1;
    $pc_beg = $q+1;
    $pc_end = Maddint($pc_beg, 100_000);
    while (1) {
      $pc_end = $B2 if $pc_end > $B2;
      my @bprimes = @{ Mprimes($pc_beg, $pc_end) };
      foreach my $i (0 .. $#bprimes) {
        my $diff = $bprimes[$i] - $q;
        $q = $bprimes[$i];
        my $qdiff = ($diff >> 1) - 1;
        $precomp_bm[$qdiff] = _bi_powmod($bm, $diff, $n)
          unless defined $precomp_bm[$qdiff];
        $pa = ($pa * $precomp_bm[$qdiff]) % $n;
        if ($pa == 0) { push @factors, $n; return @factors; }
        $b *= ($pa-1);
        if (($j++ % 128) == 0) {
          $b %= $n;
          $f = Mgcd($b, $n);
          last if $f != 1;
        }
      }
      last if $f != 1 || $pc_end >= $B2;
      ($pc_beg, $pc_end) = (Maddint($pc_end,1), Maddint($pc_end,500_000));
    }
    $f = Mgcd($b, $n);
  }
  return _found_factor($f, $n, "pminus1-bigint $B1/$B2", @factors);
}

sub cheb_factor {
  my($n, $B1, $initx, $skipbasic) = @_;
  validate_integer_nonneg($n);
  validate_integer_nonneg($B1) if defined $B1;
  validate_integer_nonneg($initx) if defined $initx;

  my @factors;
  if (!$skipbasic) {
    @factors = _basic_factor($n);
    return @factors if $n < 4;
  }

  my $x = (defined $initx && $initx > 0)  ?  $initx  :  72;  # Arbitrary
  my $B = (defined $B1 && $B1 > 0) ? $B1 : Mmulint(Mpowint(Mlogint($n,2),2),8);
  $B = Msqrtint($n) if $B > Msqrtint($n);
  my $sqrtB = Msqrtint($B);
  my $inv = Minvmod(2,$n);
  my $f = 1;

  my @bprimes = @{ Mprimes(2, $B) };
  foreach my $p (@bprimes) {
    my $xx = Maddmod($x,$x,$n);
    if ($p <= $sqrtB) {
      my $plgp = Mpowint($p, Mlogint($B, $p));
      $x = Mmulmod(Math::Prime::Util::lucasvmod($xx, 1, $plgp, $n), $inv, $n);
    } else {
      $x = Mmulmod(Math::Prime::Util::lucasvmod($xx, 1, $p, $n), $inv, $n);
    }
    $f = Mgcd($x-1, $n);
    last if $f != 1;
  }
  return _found_factor($f, $n, "cheb", @factors);
}

sub holf_factor {
  my($n, $rounds, $startrounds) = @_;
  validate_integer_nonneg($n);
  if (defined $rounds) { validate_integer_nonneg($rounds); }
  else                 { $rounds = 64*1024*1024; }
  $startrounds = 1 if (!defined $startrounds) || ($startrounds < 1);

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if (ref($n)) {
    for my $i ($startrounds .. $rounds) {
      my $ni = Mmulint($n,$i);
      my $s = Msqrtint($ni);
      if (Mmulint($s,$s) == $ni) {
        # s^2 = n*i, so m = s^2 mod n = 0.  Hence f = GCD(n, s) = GCD(n, n*i)
        my $f = Mgcd($ni, $n);
        return _found_factor($f, $n, "HOLF", @factors);
      }
      $s = Maddint($s,1);
      my $m = Msubint(Mmulint($s,$s),$ni);
      if (Mis_power($m, 2, \my $f)) {
        $f = Mgcd($n, $s > $f ? $s-$f : $f-$s);
        return _found_factor($f, $n, "HOLF ($i rounds)", @factors);
      }
    }
  } else {
    for my $i ($startrounds .. $rounds) {
      my $s = int(sqrt($n * $i));
      $s++ if ($s * $s) != ($n * $i);
      my $m = ($s < MPU_HALFWORD) ? ($s*$s) % $n : _mulmod($s, $s, $n);
      # Check for perfect square
      my $mc = $m & 31;
      next unless $mc==0||$mc==1||$mc==4||$mc==9||$mc==16||$mc==17||$mc==25;
      my $f = int(sqrt($m));
      next unless $f*$f == $m;
      $f = _gcd_ui($s - $f,  $n);
      return _found_factor($f, $n, "HOLF ($i rounds)", @factors);
    }
  }
  push @factors, $n;
  @factors;
}

sub fermat_factor {
  my($n, $rounds) = @_;
  validate_integer_nonneg($n);
  if (defined $rounds) { validate_integer_nonneg($rounds); }
  else                 { $rounds = 64*1024*1024; }

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if (ref($n)) {
    my $pa = Msqrtint($n);
    return _found_factor($pa, $n, "Fermat", @factors) if Mmulint($pa,$pa) == $n;
    $pa = Maddint($pa,1);
    my $b2 = Msubint(Mmulint($pa,$pa),$n);
    my $lasta = Maddint($pa,$rounds);
    while ($pa <= $lasta) {
      if (Mis_power($b2, 2, \my $s)) {
        my $i = Msubint($pa,($lasta-$rounds))+1;
        return _found_factor(Msubint($pa,$s), $n, "Fermat ($i rounds)", @factors);
      }
      $pa = Maddint($pa,1);
      $b2 = Msubint(Mmulint($pa,$pa),$n);
    }
  } else {
    my $pa = int(sqrt($n));
    return _found_factor($pa, $n, "Fermat", @factors) if $pa*$pa == $n;
    $pa++;
    my $b2 = $pa*$pa - $n;
    my $lasta = $pa + $rounds;
    while ($pa <= $lasta) {
      my $mc = $b2 & 31;
      if ($mc==0||$mc==1||$mc==4||$mc==9||$mc==16||$mc==17||$mc==25) {
        my $s = int(sqrt($b2));
        if ($s*$s == $b2) {
          my $i = $pa-($lasta-$rounds)+1;
          return _found_factor($pa - $s, $n, "Fermat ($i rounds)", @factors);
        }
      }
      $pa++;
      $b2 = $pa*$pa-$n;
    }
  }
  push @factors, $n;
  @factors;
}


sub ecm_factor {
  my($n, $B1, $B2, $ncurves) = @_;
  validate_integer_nonneg($n);

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if ($Math::Prime::Util::_GMPfunc{"ecm_factor"}) {
    $B1 = 0 if !defined $B1;
    $ncurves = 0 if !defined $ncurves;
    my @ef = Math::Prime::Util::GMP::ecm_factor($n, $B1, $ncurves);
    if (@ef > 1) {
      my $ecmfac = reftyped($n, $ef[-1]);
      return _found_factor($ecmfac, $n, "ECM (GMP) B1=$B1 curves $ncurves", @factors);
    }
    push @factors, $n;
    return @factors;
  }

  $ncurves = 10 unless defined $ncurves;

  if (!defined $B1) {
    for my $mul (1, 10, 100, 1000, 10_000, 100_000, 1_000_000) {
      $B1 = 100 * $mul;
      $B2 = 10*$B1;
      #warn "Trying ecm with $B1 / $B2\n";
      my @nf = ecm_factor($n, $B1, $B2, $ncurves);
      if (scalar @nf > 1) {
        push @factors, @nf;
        return @factors;
      }
    }
    push @factors, $n;
    return @factors;
  }

  $B2 = 10*$B1 unless defined $B2;
  my $sqrt_b1 = int(sqrt($B1)+1);

  # Affine code.  About 3x slower than the projective, and no stage 2.
  #
  #if (!defined $Math::Prime::Util::ECAffinePoint::VERSION) {
  #  eval { require Math::Prime::Util::ECAffinePoint; 1; }
  #  or do { croak "Cannot load Math::Prime::Util::ECAffinePoint"; };
  #}
  #my @bprimes = @{ primes(2, $B1) };
  #my $irandf = Math::Prime::Util::_get_rand_func();
  #foreach my $curve (1 .. $ncurves) {
  #  my $a = $irandf->($n-1);
  #  my $b = 1;
  #  my $ECP = Math::Prime::Util::ECAffinePoint->new($a, $b, $n, 0, 1);
  #  foreach my $q (@bprimes) {
  #    my $k = $q;
  #    if ($k < $sqrt_b1) {
  #      my $kmin = int($B1 / $q);
  #      while ($k <= $kmin) { $k *= $q; }
  #    }
  #    $ECP->mul($k);
  #    my $f = $ECP->f;
  #    if ($f != 1) {
  #      last if $f == $n;
  #      warn "ECM found factors with B1 = $B1 in curve $curve\n";
  #      return _found_factor($f, $n, "ECM B1=$B1 curve $curve", @factors);
  #    }
  #    last if $ECP->is_infinity;
  #  }
  #}

  require Math::Prime::Util::ECProjectivePoint;
  require Math::Prime::Util::RandomPrimes;

  # With multiple curves, it's better to get all the primes at once.
  # The downside is this can kill memory with a very large B1.
  my @bprimes = @{ Mprimes(3, $B1) };
  foreach my $q (@bprimes) {
    last if $q > $sqrt_b1;
    my($k,$kmin) = ($q, int($B1/$q));
    while ($k <= $kmin) { $k *= $q; }
    $q = $k;
  }
  my @b2primes = ($B2 > $B1) ? @{Mprimes($B1+1, $B2)} : ();

  foreach my $curve (1 .. $ncurves) {
    my $sigma = Murandomm($n-6) + 6;
    my ($u, $v) = ( ($sigma*$sigma - 5) % $n, (4 * $sigma) % $n );
    my ($x, $z) = ( ($u*$u*$u) % $n,  ($v*$v*$v) % $n );
    my $cb = (4 * $x * $v) % $n;
    my $ca = ( (($v-$u)**3) * (3*$u + $v) ) % $n;
    my $f = Math::BigInt::bgcd( $cb, $n );
    $f = Math::BigInt::bgcd( $z, $n ) if $f == 1;
    next if $f == $n;
    return _found_factor($f,$n, "ECM B1=$B1 curve $curve", @factors) if $f != 1;
    $cb = Math::BigInt->new("$cb") unless ref($cb) eq 'Math::BigInt';
    $u = $cb->copy->bmodinv($n);
    $ca = (($ca*$u) - 2) % $n;

    my $ECP = Math::Prime::Util::ECProjectivePoint->new($ca, $n, $x, $z);
    my $fm = $n-$n+1;
    my $i = 15;

    for (my $q = 2; $q < $B1; $q *= 2) { $ECP->double(); }
    foreach my $k (@bprimes) {
      $ECP->mul($k);
      $fm = ($fm * $ECP->x() ) % $n;
      if ($i++ % 32 == 0) {
        $f = Math::BigInt::bgcd($fm, $n);
        last if $f != 1;
      }
    }
    $f = Math::BigInt::bgcd($fm, $n);
    next if $f == $n;

    if ($f == 1 && $B2 > $B1) { # BEGIN STAGE 2
      my $D = int(sqrt($B2/2));  $D++ if $D % 2;
      my $one = $n - $n + 1;
      my $g = $one;

      my $S2P = $ECP->copy->normalize;
      $f = $S2P->f;
      if ($f != 1) {
        next if $f == $n;
        #warn "ECM S2 normalize f=$f\n" if $f != 1;
        return _found_factor($f, $n, "ECM S2 B1=$B1 curve $curve");
      }
      my $S2x = $S2P->x;
      my $S2d = $S2P->d;
      my @nqx = ($n-$n, $S2x);

      foreach my $i (2 .. 2*$D) {
        my($x2, $z2);
        if ($i % 2) {
          ($x2, $z2) = Math::Prime::Util::ECProjectivePoint::_addx($nqx[($i-1)/2], $nqx[($i+1)/2], $S2x, $n);
        } else {
          ($x2, $z2) = Math::Prime::Util::ECProjectivePoint::_double($nqx[$i/2], $one, $n, $S2d);
        }
        $nqx[$i] = $x2;
        #($f, $u, undef) = _extended_gcd($z2, $n);
        $f = Math::BigInt::bgcd( $z2, $n );
        last if $f != 1;
        $u = $z2->copy->bmodinv($n);
        $nqx[$i] = ($x2 * $u) % $n;
      }
      if ($f != 1) {
        next if $f == $n;
        #warn "ECM S2 1: B1 $B1 B2 $B2 curve $curve f=$f\n";
        return _found_factor($f, $n, "ECM S2 B1=$B1 curve $curve", @factors);
      }

      $x = $nqx[2*$D-1];
      my $m = 1;
      while ($m < ($B2+$D)) {
        if ($m != 1) {
          my $oldx = $S2x;
          my ($x1, $z1) = Math::Prime::Util::ECProjectivePoint::_addx($nqx[2*$D], $S2x, $x, $n);
          $f = Math::BigInt::bgcd( $z1, $n );
          last if $f != 1;
          $u = $z1->copy->bmodinv($n);
          $S2x = ($x1 * $u) % $n;
          $x = $oldx;
          last if $f != 1;
        }
        if ($m+$D > $B1) {
          my @p = grep { $_ >= $m-$D && $_ <= $m+$D } @b2primes;
          foreach my $i (@p) {
            last if $i >= $m;
            $g = ($g * ($S2x - $nqx[$m+$D-$i])) % $n;
          }
          foreach my $i (@p) {
            next unless $i > $m;
            next if $i > ($m+$m) || is_prime($m+$m-$i);
            $g = ($g * ($S2x - $nqx[$i-$m])) % $n;
          }
          $f = Math::BigInt::bgcd($g, $n);
          #warn "ECM S2 3: found $f in stage 2\n" if $f != 1;
          last if $f != 1;
        }
        $m += 2*$D;
      }
    } # END STAGE 2

    next if $f == $n;
    if ($f != 1) {
      #warn "ECM found factors with B1 = $B1 in curve $curve\n";
      return _found_factor($f, $n, "ECM B1=$B1 curve $curve", @factors);
    }
    # end of curve loop
  }
  push @factors, $n;
  @factors;
}

sub divisors {
  my($n,$k) = @_;
  validate_integer_nonneg($n);
  if (defined $k) {
    validate_integer_nonneg($k);
    $k = $n if $k > $n;
  } else {
    $k = $n;
  }

  if (!wantarray) {
    # In scalar context, returns sigma_0(n).  Very fast.
    return Mdivisor_sum($n,0) if $k >= $n;
    my @div = divisors($n,$k);
    return scalar(@div);
  }

  return ()  if $n == 0 || $k == 0;
  return (1) if $n == 1 || $k == 1;

  my @d;
  if ($Math::Prime::Util::_GMPfunc{"divisors"}) {
    # This trips an erroneous compile time error without the eval.
    if ($k < $n && $Math::Prime::Util::GMP::VERSION >= 0.53) {
      eval "\@d = Math::Prime::Util::GMP::divisors(\"$n\",\"$k\"); ";  ## no critic qw(ProhibitStringyEval)
    } else {
      eval "\@d = Math::Prime::Util::GMP::divisors(\"$n\"); ";  ## no critic qw(ProhibitStringyEval)
      @d = grep { $_ <= $k } @d  if $k < $n;
    }
    @d = map { $_ <= ~0 ? $_ : ref($n)->new($_) } @d   if ref($n);
    return @d;
  }

  my @pe = Mfactor_exp($n);
  return (1,$n) if @pe == 1 && $pe[0]->[1] == 1 && $n <= $k;

  @d = (1);
  for my $pe (@pe) {
    my($p,$e) = @$pe;
    last if $p > $k;
    my @t;
    push @d,  @t = map { Mmulint($_,$p) } @d;              # multiply through
    push @d,  @t = map { Mmulint($_,$p) } @t  for 2 .. $e; # repeat
  }

  @d = grep { $_ <= $k } @d  if $k < $n;
  Mvecsort(@d);
}


################################################################################


sub _chebyshev_theta {
  my($n,$low) = @_;
  my($sum,$high) = (0.0, 0);
  while ($low <= $n) {
    $high = $low + 1e6;
    $high = $n if $high > $n;
    $sum += log($_) for @{Mprimes($low,$high)};
    $low = $high+1;
  }
  $sum;
}
sub chebyshev_theta {
  my($n) = @_;
  validate_integer_nonneg($n);
  _chebyshev_theta($n,2);
}

sub chebyshev_psi {
  my($n) = @_;
  validate_integer_nonneg($n);
  return 0 if $n <= 1;
  my ($sum, $logn, $sqrtn) = (0.0, log($n), Msqrtint($n));

  # Sum the log of prime powers first
  for my $p (@{Mprimes($sqrtn)}) {
    my $logp = log($p);
    $sum += $logp * int($logn/$logp+1e-15);
  }
  # The rest all have exponent 1: add them in using the segmenting theta code
  $sum += _chebyshev_theta($n, $sqrtn+1);

  $sum;
}

sub hclassno {
  my $n = shift;
  validate_integer($n);

  return -1 if $n == 0;
  return 0 if $n < 0 || ($n % 4) == 1 || ($n % 4) == 2;
  return 2 * (2,3,6,6,6,8,12,9,6,12,18,12,8,12,18,18,12,15,24,12,6,24,30,20,12,12,24,24,18,24)[($n>>1)-1] if $n <= 60;

  my ($h, $square, $b, $b2) = (0, 0, $n & 1, ($n+1) >> 2);

  if ($b == 0) {
    my $lim = Msqrtint($b2);
    if (_is_perfect_square($b2)) {
      $square = 1;
      $lim--;
    }
    #$h += scalar(grep { $_ <= $lim } divisors($b2));
    for my $i (1 .. $lim) { $h++ unless $b2 % $i; }
    ($b,$b2) = (2, ($n+4) >> 2);
  }
  while ($b2 * 3 < $n) {
    $h++ unless $b2 % $b;
    my $lim = Msqrtint($b2);
    if (_is_perfect_square($b2)) {
      $h++;
      $lim--;
    }
    #$h += 2 * scalar(grep { $_ > $b && $_ <= $lim } divisors($b2));
    for my $i ($b+1 .. $lim) { $h += 2 unless $b2 % $i; }
    $b += 2;
    $b2 = ($n+$b*$b) >> 2;
  }
  return (($b2*3 == $n) ? 2*(3*$h+1) : $square ? 3*(2*$h+1) : 6*$h) << 1;
}

# Ramanujan Tau using Cohen's method with Hurwitz class numbers.
# Also see Lygeros (2010).
# The two hclassno calls could be collapsed with some work
sub _tauprime {
  my $p = shift;
  return -24 if $p == 2;

  my $sum = 0;
  my($p9,$pp7) = (Mmulint(9,$p), Mvecprod(7,$p,$p));
  for my $t (1 .. Msqrtint($p)) {
    my $t2 = Mpowint($t,2);
    my $v = Msubint($p,$t2);
    my $T1 = Mpowint($t2,3);
    my $T2 = Maddint( Msubint(Mvecprod(4,$t2,$t2), Mmulint($p9,$t2)), $pp7);
    my $T3;
    my $v4 = $v % 4;
    if ($v4 == 0) {
      $T3 = Maddint(Mmulint(2,Mhclassno($v)), Mhclassno(Mmulint(4,$v)) );
    } elsif ($v4 == 3) {
      $T3 = Mmulint( $v%8 == 3 ? 6 : 4, Mhclassno($v) );
    } else {
      $T3 = Mhclassno(Mmulint(4,$v));
    }

    $sum = Maddint($sum, Mvecprod($T1, $T2, $T3));
  }
  Mvecsum( Mmulint( 28, Mpowint($p,6)),
           Mmulint(-28, Mpowint($p,5)),
           Mmulint(-90, Mpowint($p,4)),
           Mmulint(-35, Mpowint($p,3)),
           -1,
           Mmulint(-32,Mdivint($sum,3)) );
}

# Recursive method for handling prime powers
sub _taupower {
  my($p, $e, $tp) = @_;
  return 1 if $e <= 0;
  $tp = _tauprime($p) unless defined $tp;

  return $tp if $e == 1;

  my $p11 = Mpowint($p,11);
  return Msubint(Mpowint($tp,2), $p11) if $e == 2;
  return Msubint(Mpowint($tp,3), Mvecprod(2,$tp,$p11)) if $e == 3;
  return Mvecsum(Mpowint($tp,4), Mvecprod(-3,Mpowint($tp,2),$p11), Mpowint($p11,2)) if $e == 4;

  # Recurse -3
  my $F3 = Msubint(Mpowint($tp,3),Mvecprod(2,$tp,$p11));
  my $F4 = Msubint(Mmulint($p11,$p11),Mvecprod($tp,$tp,$p11));
  Maddint( Mmulint($F3,_taupower($p,$e-3,$tp)),
           Mmulint($F4,_taupower($p,$e-4,$tp)) );
}

sub ramanujan_tau {
  my $n = shift;
  validate_integer_nonneg($n);
  return 0 if $n <= 0;

  # Use GMP if we have no XS or if size is small
  if ($n < 100000 || !getconfig()->{'xs'}) {
    if ($Math::Prime::Util::_GMPfunc{"ramanujan_tau"}) {
      return reftyped($_[0], Math::Prime::Util::GMP::ramanujan_tau($n));
    }
  }
  Mvecprod(map { _taupower($_->[0],$_->[1]) } Mfactor_exp($n));
}

sub _Euler {
 my($dig) = @_;
 return Math::Prime::Util::GMP::Euler($dig)
   if $dig > 70 && $Math::Prime::Util::_GMPfunc{"Euler"};
 '0.57721566490153286060651209008240243104215933593992359880576723488486772677766467';
}
sub _Li2 {
 my($dig) = @_;
 return Math::Prime::Util::GMP::li(2,$dig)
   if $dig > 70 && $Math::Prime::Util::_GMPfunc{"li"};
 '1.04516378011749278484458888919461313652261557815120157583290914407501320521';
}

sub ExponentialIntegral {
  my($x) = @_;
  return - MPU_INFINITY if $x == 0;
  return 0              if $x == - MPU_INFINITY;
  return MPU_INFINITY   if $x == MPU_INFINITY;

  if ($Math::Prime::Util::_GMPfunc{"ei"}) {
    $x = Math::BigFloat->new("$x") if defined $bignum::VERSION && ref($x) ne 'Math::BigFloat';
    return 0.0 + Math::Prime::Util::GMP::ei($x,40) if !ref($x);
    my $str = Math::Prime::Util::GMP::ei($x, _find_big_acc($x));
    return $x->copy->bzero->badd($str);
  }

  $x = Math::BigFloat->new("$x") if defined $bignum::VERSION && ref($x) ne 'Math::BigFloat';

  my $tol = 1e-16;
  my $sum = 0.0;
  my($y, $t);
  my $c = 0.0;
  my $val; # The result from one of the four methods

  if ($x < -1) {
    # Continued fraction
    my $lc = 0;
    my $ld = 1 / (1 - $x);
    $val = $ld * (-exp($x));
    for my $n (1 .. 100000) {
      $lc = 1 / (2*$n + 1 - $x - $n*$n*$lc);
      $ld = 1 / (2*$n + 1 - $x - $n*$n*$ld);
      my $old = $val;
      $val *= $ld/$lc;
      last if abs($val - $old) <= ($tol * abs($val));
    }
  } elsif ($x < 0) {
    # Rational Chebyshev approximation
    my @C6p = ( -148151.02102575750838086,
                 150260.59476436982420737,
                  89904.972007457256553251,
                  15924.175980637303639884,
                   2150.0672908092918123209,
                    116.69552669734461083368,
                      5.0196785185439843791020);
    my @C6q = (  256664.93484897117319268,
                 184340.70063353677359298,
                  52440.529172056355429883,
                   8125.8035174768735759866,
                    750.43163907103936624165,
                     40.205465640027706061433,
                      1.0000000000000000000000);
    my $sumn = $C6p[0]-$x*($C6p[1]-$x*($C6p[2]-$x*($C6p[3]-$x*($C6p[4]-$x*($C6p[5]-$x*$C6p[6])))));
    my $sumd = $C6q[0]-$x*($C6q[1]-$x*($C6q[2]-$x*($C6q[3]-$x*($C6q[4]-$x*($C6q[5]-$x*$C6q[6])))));
    $val = log(-$x) - ($sumn / $sumd);
  } elsif ($x < -log($tol)) {
    # Convergent series
    my $fact_n = 1;
    $y = _Euler(18)-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    $y = log($x)-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    for my $n (1 .. 200) {
      $fact_n *= $x/$n;
      my $term = $fact_n / $n;
      $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      last if $term < $tol;
    }
    $val = $sum;
  } else {
    # Asymptotic divergent series
    my $invx = 1.0 / $x;
    my $term = $invx;
    $sum = 1.0 + $term;
    for my $n (2 .. 200) {
      my $last_term = $term;
      $term *= $n * $invx;
      last if $term < $tol;
      if ($term < $last_term) {
        $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      } else {
        $y = (-$last_term/3)-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
        last;
      }
    }
    $val = exp($x) * $invx * $sum;
  }
  $val;
}

sub LogarithmicIntegral {
  my($x,$opt) = @_;
  return 0              if $x == 0;
  return - MPU_INFINITY if $x == 1;
  return MPU_INFINITY   if $x == MPU_INFINITY;
  croak "Invalid input to LogarithmicIntegral:  x must be > 0" if $x <= 0;
  $opt = 0 unless defined $opt;

  if ($Math::Prime::Util::_GMPfunc{"li"}) {
    $x = Math::BigFloat->new("$x") if defined $bignum::VERSION && ref($x) ne 'Math::BigFloat';
    return 0.0 + Math::Prime::Util::GMP::li($x,40) if !ref($x);
    my $str = Math::Prime::Util::GMP::li($x, _find_big_acc($x));
    return $x->copy->bzero->badd($str);
  }

  if ($x == 2) {
    my $li2const = (ref($x) eq 'Math::BigFloat') ? Math::BigFloat->new(_Li2(_find_big_acc($x))) : 0.0+_Li2(30);
    return $li2const;
  }

  if (defined $bignum::VERSION) {
    # If bignum is on, always use Math::BigFloat.
    $x = Math::BigFloat->new("$x") if ref($x) ne 'Math::BigFloat';
  } elsif (ref($x)) {
    # bignum is off, use native if small, BigFloat otherwise.
    if ($x <= 1e16) {
      $x = _bigint_to_int($x);
    } else {
      $x = _upgrade_to_float($x) if ref($x) ne 'Math::BigFloat';
    }
  }
  # Make sure we preserve whatever accuracy setting the input was using.
  $x->accuracy($_[0]->accuracy) if ref($x) && ref($_[0]) =~ /^Math::Big/ && $_[0]->accuracy;

  # Do divergent series here for big inputs.  Common for big pc approximations.
  # Why is this here?
  #   1) exp(log(x)) results in a lot of lost precision
  #   2) exp(x) with lots of precision turns out to be really slow, and in
  #      this case it was unnecessary.
  my $tol = 1e-16;
  my $xdigits = 0;
  my $finalacc = 0;
  if (ref($x) =~ /^Math::Big/) {
    $xdigits = _find_big_acc($x);
    my $xlen = length($x->copy->bfloor->bstr());
    $xdigits = $xlen if $xdigits < $xlen;
    $finalacc = $xdigits;
    $xdigits += length(int(log(0.0+"$x"))) + 1;
    $tol = Math::BigFloat->new(10)->bpow(-$xdigits);
    $x->accuracy($xdigits);
  }
  my $logx = $xdigits ? $x->copy->blog(undef,$xdigits) : log($x);

  # TODO: See if we can tune this
  if (0 && $x >= 1) {
    _upgrade_to_float();
    my $sum = Math::BigFloat->new(0);
    my $inner_sum = Math::BigFloat->new(0);
    my $p = Math::BigFloat->new(-1);
    my $factorial = 1;
    my $power2 = 1;
    my $q;
    my $k = 0;
    my $neglogx = -$logx;
    for my $n (1 .. 1000) {
      $factorial = mulint($factorial, $n);
      $q = mulint($factorial, $power2);
      $power2 = mulint(2, $power2);
      while ($k <= ($n-1)>>1) {
        $inner_sum += Math::BigFloat->new(1) / (2*$k+1);
        $k++;
      }
      $p *= $neglogx;
      my $term = ($p / $q) * $inner_sum;
      $sum += $term;
      last if abs($term) < $tol;
    }
    $sum *= sqrt($x);
    return 0.0+_Euler(18) + log($logx) + $sum unless ref($x)=~/^Math::Big/;
    my $val = Math::BigFloat->new(_Euler(40))->badd("".log($logx))->badd("$sum");
    $val->accuracy($finalacc) if $xdigits;
    return $val;
  }

  if ($x > 1e16) {
    my $invx = ref($logx) ? Math::BigFloat->bone / $logx : 1.0/$logx;
    # n = 0  =>  0!/(logx)^0 = 1/1 = 1
    # n = 1  =>  1!/(logx)^1 = 1/logx
    my $term = $invx;
    my $sum = 1.0 + $term;
    for my $n (2 .. 1000) {
      my $last_term = $term;
      $term *= $n * $invx;
      last if $term < $tol;
      if ($term < $last_term) {
        $sum += $term;
      } else {
        $sum -= ($last_term/3);
        last;
      }
      $term->bround($xdigits) if $xdigits;
    }
    $invx *= $sum;
    $invx *= $x;
    $invx->accuracy($finalacc) if ref($invx) && $xdigits;
    return $invx;
  }
  # Convergent series.
  if ($x >= 1) {
    my $fact_n = 1.0;
    my $nfac = 1.0;
    my $sum  = 0.0;
    for my $n (1 .. 200) {
      $fact_n *= $logx/$n;
      my $term = $fact_n / $n;
      $sum += $term;
      last if $term < $tol;
      $term->bround($xdigits) if $xdigits;
    }

    return 0.0+_Euler(18) + log($logx) + $sum unless ref($x) =~ /^Math::Big/;

    my $val = Math::BigFloat->new(_Euler(40))->badd("".log($logx))->badd("$sum");
    $val->accuracy($finalacc) if $xdigits;
    return $val;
  }

  ExponentialIntegral($logx);
}

# Riemann Zeta function for native integers.
my @_Riemann_Zeta_Table = (
  0.6449340668482264364724151666460251892,  # zeta(2) - 1
  0.2020569031595942853997381615114499908,
  0.0823232337111381915160036965411679028,
  0.0369277551433699263313654864570341681,
  0.0173430619844491397145179297909205279,
  0.0083492773819228268397975498497967596,
  0.0040773561979443393786852385086524653,
  0.0020083928260822144178527692324120605,
  0.0009945751278180853371459589003190170,
  0.0004941886041194645587022825264699365,
  0.0002460865533080482986379980477396710,
  0.0001227133475784891467518365263573957,
  0.0000612481350587048292585451051353337,
  0.0000305882363070204935517285106450626,
  0.0000152822594086518717325714876367220,
  0.0000076371976378997622736002935630292,
  0.0000038172932649998398564616446219397,
  0.0000019082127165539389256569577951013,
  0.0000009539620338727961131520386834493,
  0.0000004769329867878064631167196043730,
  0.0000002384505027277329900036481867530,
  0.0000001192199259653110730677887188823,
  0.0000000596081890512594796124402079358,
  0.0000000298035035146522801860637050694,
  0.0000000149015548283650412346585066307,
  0.0000000074507117898354294919810041706,
  0.0000000037253340247884570548192040184,
  0.0000000018626597235130490064039099454,
  0.0000000009313274324196681828717647350,
  0.0000000004656629065033784072989233251,
  0.0000000002328311833676505492001455976,
  0.0000000001164155017270051977592973835,
  0.0000000000582077208790270088924368599,
  0.0000000000291038504449709968692942523,
  0.0000000000145519218910419842359296322,
  0.0000000000072759598350574810145208690,
  0.0000000000036379795473786511902372363,
  0.0000000000018189896503070659475848321,
  0.0000000000009094947840263889282533118,
);


sub RiemannZeta {
  my($x) = @_;

  # Try our GMP code if possible.
  if ($Math::Prime::Util::_GMPfunc{"zeta"}) {
    my($wantbf,$xdigits) = _bfdigits($x);
    # If we knew the *exact* number of zero digits, we could let GMP zeta
    # handle the correct rounding.  But we don't, so we have to go over.
    my $zero_dig = "".int($x / 3) - 1;
    my $strval = Math::Prime::Util::GMP::zeta($x, $xdigits + 8 + $zero_dig);
    if ($strval =~ s/^(1\.0*)/./) {
      $strval .= "e-".(length($1)-2) if length($1) > 2;
    } else {
      $strval =~ s/^(-?\d+)/$1-1/e;
    }

    return ($wantbf)  ?  Math::BigFloat->new($strval,$wantbf)  : 0.0 + $strval;
  }

  # If we need a bigfloat result, then call our PP routine.
  if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
    require Math::Prime::Util::ZetaBigFloat;
    return Math::Prime::Util::ZetaBigFloat::RiemannZeta($x);
  }

  # Native float results
  return 0.0 + $_Riemann_Zeta_Table[int($x)-2]
    if $x == int($x) && defined $_Riemann_Zeta_Table[int($x)-2];
  my $tol = 1.11e-16;

  # Series based on (2n)! / B_2n.
  # This is a simplification of the Cephes zeta function.
  my @A = (
      12.0,
     -720.0,
      30240.0,
     -1209600.0,
      47900160.0,
     -1892437580.3183791606367583212735166426,
      74724249600.0,
     -2950130727918.1642244954382084600497650,
      116467828143500.67248729113000661089202,
     -4597978722407472.6105457273596737891657,
      181521054019435467.73425331153534235290,
     -7166165256175667011.3346447367083352776,
      282908877253042996618.18640556532523927,
  );
  my $s = 0.0;
  my $rb = 0.0;
  foreach my $i (2 .. 10) {
    $rb = $i ** -$x;
    $s += $rb;
    return $s if abs($rb/$s) < $tol;
  }
  my $w = 10.0;
  $s = $s  +  $rb*$w/($x-1.0)  -  0.5*$rb;
  my $ra = 1.0;
  foreach my $i (0 .. 12) {
    my $k = 2*$i;
    $ra *= $x + $k;
    $rb /= $w;
    my $t = $ra*$rb/$A[$i];
    $s += $t;
    $t = abs($t/$s);
    last if $t < $tol;
    $ra *= $x + $k + 1.0;
    $rb /= $w;
  }
  return $s;
}

# Riemann R function
sub RiemannR {
  my($x) = @_;

  croak "Invalid input to ReimannR:  x must be > 0" if $x <= 0;

  # With MPU::GMP v0.49 this is fast.
  if ($Math::Prime::Util::_GMPfunc{"riemannr"}) {
    my($wantbf,$xdigits) = _bfdigits($x);
    my $strval = Math::Prime::Util::GMP::riemannr($x, $xdigits);
    return ($wantbf)  ?  Math::BigFloat->new($strval,$wantbf)  :  0.0 + $strval;
  }


# TODO: look into this as a generic solution
if (0 && $Math::Prime::Util::_GMPfunc{"zeta"}) {
  my($wantbf,$xdigits) = _bfdigits($x);
  $x = _upgrade_to_float($x);

  my $extra_acc = 4;
  $xdigits += $extra_acc;
  $x->accuracy($xdigits);

  my $logx = log($x);
  my $part_term = $x->copy->bone;
  my $sum = $x->copy->bone;
  my $tol = $x->copy->bone->brsft($xdigits-1, 10);
  my $bigk = $x->copy->bone;
  my $term;
  for my $k (1 .. 10000) {
    $part_term *= $logx / $bigk;
    my $zarg = $bigk->copy->binc;
    my $zeta = (RiemannZeta($zarg) * $bigk) + $bigk;
    #my $strval = Math::Prime::Util::GMP::zeta($k+1, $xdigits + int(($k+1) / 3));
    #my $zeta = Math::BigFloat->new($strval)->bdec->bmul($bigk)->badd($bigk);
    $term = $part_term / $zeta;
    $sum += $term;
    last if $term < ($tol * $sum);
    $bigk->binc;
  }
  $sum->bround($xdigits-$extra_acc);
  my $strval = "$sum";
  return ($wantbf)  ?  Math::BigFloat->new($strval,$wantbf)  :  0.0 + $strval;
}

  if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
    require Math::Prime::Util::ZetaBigFloat;
    return Math::Prime::Util::ZetaBigFloat::RiemannR($x);
  }

  my $sum = 0.0;
  my $tol = 1e-18;
  my($c, $y, $t) = (0.0);
  if ($x > 10**17) {
    my @mob = Mmoebius(0,300);
    for my $k (1 .. 300) {
      next if $mob[$k] == 0;
      my $term = $mob[$k] / $k * MLi($x**(1.0/$k));
      $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      last if abs($term) < ($tol * abs($sum));
    }
  } else {
    $y = 1.0-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    my $flogx = log($x);
    my $part_term = 1.0;
    for my $k (1 .. 10000) {
      my $zeta = ($k <= $#_Riemann_Zeta_Table)
                 ? $_Riemann_Zeta_Table[$k+1-2]    # Small k from table
                 : RiemannZeta($k+1);              # Large k from function
      $part_term *= $flogx / $k;
      my $term = $part_term / ($k + $k * $zeta);
      $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      last if $term < ($tol * $sum);
    }
  }
  return $sum;
}

sub LambertW {
  my $x = shift;
  croak "Invalid input to LambertW:  x must be >= -1/e" if $x < -0.36787944118;
  $x = _upgrade_to_float($x) if ref($x) eq 'Math::BigInt';
  my $xacc = ref($x) ? _find_big_acc($x) : 0;
  my $w;

  if ($Math::Prime::Util::_GMPfunc{"lambertw"}) {
    my $w = (!$xacc)
          ? 0.0 + Math::Prime::Util::GMP::lambertw($x)
          : $x->copy->bzero->badd(Math::Prime::Util::GMP::lambertw($x, $xacc));
    return $w;
  }

  # Approximation
  if ($x < -0.06) {
    my $ti = $x * 2 * exp($x-$x+1) + 2;
    return -1 if $ti <= 0;
    my $t  = sqrt($ti);
    $w = (-1 + 1/6*$t + (257/720)*$t*$t + (13/720)*$t*$t*$t) / (1 + (5/6)*$t + (103/720)*$t*$t);
  } elsif ($x < 1.363) {
    my $l1 = log($x + 1);
    $w = $l1 * (1 - log(1+$l1) / (2+$l1));
  } elsif ($x < 3.7) {
    my $l1 = log($x);
    my $l2 = log($l1);
    $w = $l1 - $l2 - log(1 - $l2/$l1)/2.0;
  } else {
    my $l1 = log($x);
    my $l2 = log($l1);
    my $d1 = 2 * $l1 * $l1;
    my $d2 = 3 * $l1 * $d1;
    my $d3 = 2 * $l1 * $d2;
    my $d4 = 5 * $l1 * $d3;
    $w = $l1 - $l2 + $l2/$l1 + $l2*($l2-2)/$d1
       + $l2*(6+$l2*(-9+2*$l2))/$d2
       + $l2*(-12+$l2*(36+$l2*(-22+3*$l2)))/$d3
       + $l2*(60+$l2*(-300+$l2*(350+$l2*(-125+12*$l2))))/$d4;
  }

  # Now iterate to get the answer
  #
  # Newton:
  #   $w = $w*(log($x) - log($w) + 1) / ($w+1);
  # Halley:
  #   my $e = exp($w);
  #   my $f = $w * $e - $x;
  #   $w -= $f / ($w*$e+$e - ($w+2)*$f/(2*$w+2));
  #
  # Also see https://arxiv.org/pdf/2008.06122

  if (!$xacc) {
    my $eps = 1.054e-8;  # sqrt(double eps)
    for (1 .. 200) {
      last if $w == 0;
      my $w1 = $w + 1;
      my $zn = log($x/$w) - $w;
      my $qn = $w1 * 2 * ($w1+(2*$zn/3));
      my $en = ($zn/$w1) * ($qn-$zn)/($qn-$zn*2);
      my $wen = $w * $en;
      $w += $wen;
      last if abs($en) < $eps;
    }
    return $w;
  }

  # Fritsch converges quadratically, so tolerance could be 4x smaller.  Use 2x.
  my $tol = 10**(-int(1+$xacc/2));
  $w->accuracy($xacc+15);
  for (1 .. 200) {
    last if $w == 0;
    my $w1 = $w + 1;
    my $zn = log($x/$w) - $w;
    my $qn = $w1 * 2 * ($w1+(2*$zn/3));
    my $en = ($zn/$w1) * ($qn-$zn)/($qn-$zn*2);
    my $wen = $w * $en;
    $w += $wen;
    last if abs($wen) < $tol;
  }
  $w->accuracy($xacc);
  $w;
}

my $_Pi = "3.141592653589793238462643383279503";
sub Pi {
  my $digits = shift;
  return 0.0+$_Pi unless $digits;
  return 0.0+sprintf("%.*lf", $digits-1, $_Pi) if $digits < 15;
  return _upgrade_to_float($_Pi, $digits) if $digits < 30;

  # Performance ranking:
  #   MPU::GMP         Uses AGM or Ramanujan/Chudnosky with binary splitting
  #   MPFR             Uses AGM, from 1x to 1/4x the above
  #   Perl AGM w/GMP   also AGM, nice growth rate, but slower than above
  #   C pidigits       much worse than above, but faster than the others
  #   Perl AGM         without Math::BigInt::GMP, it's sluggish
  #   Math::BigFloat   new versions use AGM, old ones are *very* slow
  #
  # With a few thousand digits, any of the top 4 are fine.
  # At 10k digits, the first two are pulling away.
  # At 50k digits, the first three are 5-20x faster than C pidigits, and
  #   pray you're not having to the Perl BigFloat methods without GMP.
  # At 100k digits, the first two are 15x faster than the third, C pidigits
  #   is 200x slower, and the rest thousands of times slower.
  # At 1M digits, the first is under 1 second, MPFR under 2 seconds,
  #   Perl AGM (Math::BigInt::GMP) is over a minute, and C piigits at 1.5 hours.
  #
  # Interestingly, Math::BigInt::Pari, while greatly faster than Calc, is
  # *much* slower than GMP for these operations (both AGM and Machin).  While
  # Perl AGM with the Math::BigInt::GMP backend will pull away from C pidigits,
  # using it with the other backends doesn't do so.
  #
  # The GMP program at https://gmplib.org/download/misc/gmp-chudnovsky.c
  # will run ~4x faster than MPFR and ~1.5x faster than MPU::GMP.

  my $have_bigint_gmp = Math::BigInt->config()->{lib} =~ /GMP/;
  my $have_xdigits    = getconfig()->{'xs'};
  my $_verbose        = getconfig()->{'verbose'};

  if ($Math::Prime::Util::_GMPfunc{"Pi"}) {
    print "  using MPUGMP for Pi($digits)\n" if $_verbose;
    return _upgrade_to_float( Math::Prime::Util::GMP::Pi($digits) );
  }

  # We could consider looking for Math::MPFR or Math::Pari

  # This has a *much* better growth rate than the later solutions.
  if ( !$have_xdigits || ($have_bigint_gmp && $digits > 100) ) {
    print "  using Perl AGM for Pi($digits)\n" if $_verbose;
    # Brent-Salamin (aka AGM or Gauss-Legendre)
    $digits += 8;
    my $HALF = _upgrade_to_float(0.5);
    my ($an, $bn, $tn, $pn) = ($HALF->copy->bone, $HALF->copy->bsqrt($digits),
                               $HALF->copy->bmul($HALF), $HALF->copy->bone);
    while ($pn < $digits) {
      my $prev_an = $an->copy;
      $an->badd($bn)->bmul($HALF, $digits);
      $bn->bmul($prev_an)->bsqrt($digits);
      $prev_an->bsub($an);
      $tn->bsub($pn * $prev_an * $prev_an);
      $pn->badd($pn);
    }
    $an->badd($bn);
    $an->bmul($an,$digits)->bdiv(4*$tn, $digits-8);
    return $an;
  }

  # Spigot method in C.  Low overhead but not good growth rate.
  if ($have_xdigits) {
    print "  using XS spigot for Pi($digits)\n" if $_verbose;
    return _upgrade_to_float(Math::Prime::Util::_pidigits($digits));
  }

  # We're going to have to use the Math::BigFloat code.
  # 1) it rounds incorrectly (e.g. 761, 1372, 1509,...).
  #    Fix by adding some digits and rounding.
  # 2) AGM is *much* faster once past ~2000 digits
  # 3) It is very slow without the GMP backend.  The Pari backend helps
  #    but it still pretty bad.  With Calc it's glacial for large inputs.

  #           Math::BigFloat                AGM              spigot   AGM
  # Size     GMP    Pari  Calc        GMP    Pari  Calc        C      C+GMP
  #   500   0.04    0.60   0.30      0.08    0.10   0.47      0.09    0.06
  #  1000   0.04    0.11   1.82      0.09    0.14   1.82      0.09    0.06
  #  2000   0.07    0.37  13.5       0.09    0.34   9.16      0.10    0.06
  #  4000   0.14    2.17 107.8       0.12    1.14  39.7       0.20    0.06
  #  8000   0.52   15.7              0.22    4.63 186.2       0.56    0.08
  # 16000   2.73  121.8              0.52   19.2              2.00    0.08
  # 32000  15.4                      1.42                     7.78    0.12
  #                                   ^                        ^       ^
  #                                   |      use this THIRD ---+       |
  #                use this SECOND ---+                                |
  #                                                  use this FIRST ---+
  # approx
  # growth  5.6x    7.6x   8.0x      2.7x    4.1x   4.7x      3.9x    2.0x

  print "  using BigFloat for Pi($digits)\n" if $_verbose;
  _upgrade_to_float(0);
  return Math::BigFloat::bpi($digits+10)->round($digits);
}

################################################################################

sub forprimes {
  my($sub, $beg, $end) = @_;
  if (defined $end) { _validate_integer_nonneg($beg); }
  else              { ($beg,$end) = (2, $beg);        }
  _validate_integer_nonneg($end);
  $beg = 2 if $beg < 2;

  my $oldforexit = Math::Prime::Util::_start_for_loop();
  {
    my $pp;
    local *_ = \$pp;
    for (my $p = Mnext_prime($beg-1);  $p <= $end;  $p = Mnext_prime($p)) {
      $pp = $p;
      $sub->();
      last if Math::Prime::Util::_get_forexit();
    }
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}


sub _forcomp_sub {
  my($what, $sub, $beg, $end) = @_;
  if (defined $end) { _validate_integer_nonneg($beg); }
  else              { ($beg,$end) = (0, $beg);        }
  _validate_integer_nonneg($end);

  my $cinc = 1;
  if ($what eq 'oddcomposites') {
    $beg = 9 if $beg < 9;
    $beg++ unless $beg & 1;
    $cinc = 2;
  } else {
    $beg = 4 if $beg < 4;
  }
  $end = tobigint(~0) if $end == ~0 && !ref($end);
  my $oldforexit = Math::Prime::Util::_start_for_loop();
  {
    my $pp;
    local *_ = \$pp;
    for (my $p = Mnext_prime($beg-1);  $beg <= $end;  $p = Mnext_prime($p)) {
      for ( ; $beg < $p && $beg <= $end ; $beg += $cinc ) {
        $pp = $beg;
        $sub->();
        last if Math::Prime::Util::_get_forexit();
      }
      $beg += $cinc;
      last if Math::Prime::Util::_get_forexit();
    }
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}
sub forcomposites {
  _forcomp_sub('composites', @_);
}
sub foroddcomposites {
  _forcomp_sub('oddcomposites', @_);
}
sub forsemiprimes {
  foralmostprimes($_[0], 2, $_[1], $_[2]);
}

sub _forfac_sub {
  my($sf, $sub, $beg, $end) = @_;
  if (defined $end) { _validate_integer_nonneg($beg); }
  else              { ($beg,$end) = (1, $beg);        }
  _validate_integer_nonneg($end);
  $beg = 1 if $beg < 1;

  my $oldforexit = Math::Prime::Util::_start_for_loop();
  {
    my $pp;
    local *_ = \$pp;
    while ($beg <= $end) {
      if (!$sf || Mis_square_free($beg)) {
        $pp = $beg;
        if ($sf == 2) {
          $sub->();
        } else {
          my @f = Mfactor($beg);
          $sub->(@f);
        }
        last if Math::Prime::Util::_get_forexit();
      }
      $beg++;
    }
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}
sub forfactored {
  _forfac_sub(0, @_);
}
sub forsquarefree {
  _forfac_sub(1, @_);
}
sub forsquarefreeint {
  _forfac_sub(2, @_);
}

sub fordivisors {
  my($sub, $n) = @_;
  _validate_integer_nonneg($n);
  my @divisors = Mdivisors($n);
  my $oldforexit = Math::Prime::Util::_start_for_loop();
  {
    my $pp;
    local *_ = \$pp;
    foreach my $d (@divisors) {
      $pp = $d;
      $sub->();
      last if Math::Prime::Util::_get_forexit();
    }
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}

sub forpart {
  my($sub, $n, $rhash) = @_;
  _forcompositions(1, $sub, $n, $rhash);
}
sub forcomp {
  my($sub, $n, $rhash) = @_;
  _forcompositions(0, $sub, $n, $rhash);
}
sub _forcompositions {
  my($ispart, $sub, $n, $rhash) = @_;
  validate_integer_nonneg($n);
  my($mina, $maxa, $minn, $maxn, $primeq) = (1,$n,1,$n,-1);
  if (defined $rhash) {
    croak "forpart second argument must be a hash reference"
      unless ref($rhash) eq 'HASH';
    if (defined $rhash->{amin}) {
      $mina = $rhash->{amin};
      validate_integer_nonneg($mina);
    }
    if (defined $rhash->{amax}) {
      $maxa = $rhash->{amax};
      validate_integer_nonneg($maxa);
    }
    $minn = $maxn = $rhash->{n} if defined $rhash->{n};
    $minn = $rhash->{nmin} if defined $rhash->{nmin};
    $maxn = $rhash->{nmax} if defined $rhash->{nmax};
    validate_integer_nonneg($minn);
    validate_integer_nonneg($maxn);
    if (defined $rhash->{prime}) {
      $primeq = $rhash->{prime};
      validate_integer_nonneg($primeq);
    }
   $mina = 1 if $mina < 1;
   $maxa = $n if $maxa > $n;
   $minn = 1 if $minn < 1;
   $maxn = $n if $maxn > $n;
   $primeq = 2 if $primeq != -1 && $primeq != 0;
  }

  $sub->() if $n == 0 && $minn <= 1;
  return if $n < $minn || $minn > $maxn || $mina > $maxa || $maxn <= 0 || $maxa <= 0;

  my $oldforexit = Math::Prime::Util::_start_for_loop();
  my ($x, $y, $r, $k);
  my @a = (0) x ($n);
  $k = 1;
  $a[0] = $mina - 1;
  $a[1] = $n - $mina + 1;
  while ($k != 0) {
    $x = $a[$k-1]+1;
    $y = $a[$k]-1;
    $k--;
    $r = $ispart ? $x : 1;
    while ($r <= $y) {
      $a[$k] = $x;
      $x = $r;
      $y -= $x;
      $k++;
    }
    $a[$k] = $x + $y;
    # Restrict size
    while ($k+1 > $maxn) {
      $a[$k-1] += $a[$k];
      $k--;
    }
    next if $k+1 < $minn;
    # Restrict values
    if ($mina > 1 || $maxa < $n) {
      last if $a[0] > $maxa;
      if ($ispart) {
        next if $a[$k] > $maxa;
      } else {
        next if Mvecany(sub{ $_ < $mina || $_ > $maxa }, @a[0..$k]);
      }
    }
    next if $primeq == 0 && Mvecany(sub{ Mis_prime($_) }, @a[0..$k]);
    next if $primeq == 2 && Mvecany(sub{ !Mis_prime($_) }, @a[0..$k]);
    last if Math::Prime::Util::_get_forexit();
    $sub->(@a[0 .. $k]);
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}
sub forcomb {
  my($sub, $n, $k) = @_;
  validate_integer_nonneg($n);

  my($begk, $endk);
  if (defined $k) {
    validate_integer_nonneg($k);
    return if $k > $n;
    $begk = $endk = $k;
  } else {
    $begk = 0;
    $endk = $n;
  }

  my $oldforexit = Math::Prime::Util::_start_for_loop();
  for my $k ($begk .. $endk) {
    if ($k == 0) {
      $sub->();
    } else {
      my @c = 0 .. $k-1;
      while (1) {
        $sub->(@c);
        last if Math::Prime::Util::_get_forexit();
        next if $c[-1]++ < $n-1;
        my $i = $k-2;
        $i-- while $i >= 0 && $c[$i] >= $n-($k-$i);
        last if $i < 0;
        $c[$i]++;
        while (++$i < $k) { $c[$i] = $c[$i-1] + 1; }
      }
    }
    last if Math::Prime::Util::_get_forexit();
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}
sub _forperm {
  my($sub, $n, $all_perm) = @_;
  my $k = $n;
  my @c = reverse 0 .. $k-1;
  my $inc = 0;
  my $send = 1;
  my $oldforexit = Math::Prime::Util::_start_for_loop();
  while (1) {
    if (!$all_perm) {   # Derangements via simple filtering.
      $send = 1;
      for my $p (0 .. $#c) {
        if ($c[$p] == $k-$p-1) {
          $send = 0;
          last;
        }
      }
    }
    if ($send) {
      $sub->(reverse @c);
      last if Math::Prime::Util::_get_forexit();
    }
    if (++$inc & 1) {
      @c[0,1] = @c[1,0];
      next;
    }
    my $j = 2;
    $j++ while $j < $k && $c[$j] > $c[$j-1];
    last if $j >= $k;
    my $m = 0;
    $m++ while $c[$j] > $c[$m];
    @c[$j,$m] = @c[$m,$j];
    @c[0..$j-1] = reverse @c[0..$j-1];
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}
sub forperm {
  my($sub, $n, $k) = @_;
  validate_integer_nonneg($n);
  croak "Too many arguments for forperm" if defined $k;
  return $sub->() if $n == 0;
  return $sub->(0) if $n == 1;
  _forperm($sub, $n, 1);
}
sub forderange {
  my($sub, $n, $k) = @_;
  validate_integer_nonneg($n);
  croak "Too many arguments for forderange" if defined $k;
  return $sub->() if $n == 0;
  return if $n == 1;
  _forperm($sub, $n, 0);
}

sub _multiset_permutations {
  my($sub, $prefix, $ar, $sum) = @_;

  return if $sum == 0;

  # Remove any values with 0 occurances
  my @n = grep { $_->[1] > 0 } @$ar;

  if ($sum == 1) {                       # A single value
    $sub->(@$prefix, $n[0]->[0]);
  } elsif ($sum == 2) {                  # Optimize the leaf case
    my($n0,$n1) = map { $_->[0] } @n;
    if (@n == 1) {
      $sub->(@$prefix, $n0, $n0);
    } else {
      $sub->(@$prefix, $n0, $n1);
      $sub->(@$prefix, $n1, $n0) unless Math::Prime::Util::_get_forexit();
    }
  } elsif (0 && $sum == scalar(@n)) {         # All entries have 1 occurance
    # TODO:  Figure out a way to use this safely.  We need to capture any
    #        lastfor that was seen in the forperm.
    my @i = map { $_->[0] } @n;
    Math::Prime::Util::forperm(sub { $sub->(@$prefix, @i[@_]) }, 1+$#i);
  } else {                               # Recurse over each leading value
    for my $v (@n) {
      $v->[1]--;
      push @$prefix, $v->[0];
      no warnings 'recursion';
      _multiset_permutations($sub, $prefix, \@n, $sum-1);
      pop @$prefix;
      $v->[1]++;
      last if Math::Prime::Util::_get_forexit();
    }
  }
}

sub numtoperm {
  my($n,$k) = @_;
  validate_integer_nonneg($n);
  validate_integer($k);
  return () if $n == 0;
  return (0) if $n == 1;
  my $f = Mfactorial($n-1);
  $k %= Mmulint($f,$n) if $k < 0 || int($k/$f) >= $n;
  my @S = map { $_ } 0 .. $n-1;
  my @V;
  while ($n-- > 0) {
    my $i = int($k/$f);
    push @V, splice(@S,$i,1);
    last if $n == 0;
    $k -= $i*$f;
    $f /= $n;
  }
  @V;
}

sub permtonum {
  my $A = shift;
  croak "permtonum argument must be an array reference"
    unless ref($A) eq 'ARRAY';
  my $n = scalar(@$A);
  return 0 if $n == 0;
  {
    my %S;
    for my $v (@$A) {
      croak "permtonum invalid permutation array"
        if !defined $v || $v < 0 || $v >= $n || $S{$v}++;
    }
  }
  my $f = factorial($n-1);
  my $rank = 0;
  for my $i (0 .. $n-2) {
    my $k = 0;
    for my $j ($i+1 .. $n-1) {
      $k++ if $A->[$j] < $A->[$i];
    }
    $rank = Maddint($rank, Mmulint($k,$f));
    $f /= $n-$i-1;
  }
  $rank;
}

sub randperm {
  my($n,$k) = @_;
  validate_integer_nonneg($n);
  if (defined $k) {
    validate_integer_nonneg($k);
  }
  $k = $n if !defined($k) || $k > $n;
  return () if $k == 0;

  my @S;
  if ("$k"/"$n" <= 0.30) {
    my %seen;
    my $v;
    for my $i (1 .. $k) {
      do { $v = Murandomm($n); } while $seen{$v}++;
      push @S,$v;
    }
  } else {
    @S = (0..$n-1);
    for my $i (0 .. $n-2) {
      last if $i >= $k;
      my $j = Murandomm($n-$i);
      @S[$i,$i+$j] = @S[$i+$j,$i];
    }
    $#S = $k-1;
  }
  return @S;
}

sub shuffle {
  my @S=@_;
  # Note: almost all the time is spent in urandomm.
  for (my $i = $#S; $i >= 1; $i--) {
    my $j = Murandomm($i+1);
    @S[$i,$j] = @S[$j,$i];
  }
  @S;
}

sub vecsample {
  my $k = shift;
  return () if $k == 0 || @_ == 0;
  my $R = $_[0];
  my $isarr = (@_ > 1 || !ref($R) || ref($R) ne 'ARRAY');
  my $len = $isarr  ?  scalar(@_)  :  scalar(@$R);

  $k = $len if $k > $len;
  my @I = ($len-1, 0 .. $len-2);
  my $j;
  my @O = map { $j = Murandomm(scalar(@I));    # random index from remaining
                @I[0,$j] = @I[$j,0];           # move to front
                shift @I;                      # take it off
              } 1 .. $k;
  return $isarr  ?  @_[@O]  :  @$R[@O];
}

###############################################################################

sub vecsort {
  my(@s) = @_;
  # If we have a single array reference, unpack it.
  @s = @{$s[0]} if scalar(@s) == 1 && (ref($s[0]) || '') eq 'ARRAY';

  # Validate and convert everything into a native int or bigint
  validate_integer($_) for @s;

  # See https://github.com/perl/perl5/issues/12803 for various discussion.
  # Optimize to skip the sorting.
  return scalar(@s) unless wantarray;

  # Before Perl 5.26, numerical sort used doubles (sigh).
  if ($] < 5.026) {
    @s = sort { 0+($a<=>$b) } @s;  # Prevent sort from using built-in compare
  } else {
    @s = sort { $a<=>$b } @s;
  }
  return @s;
}

# In-place sort.
sub vecsorti {
  my($r) = @_;
  croak 'Not an array reference' unless (ref($r) || '') eq 'ARRAY';
  validate_integer($_) for @$r;
  if ($] < 5.026) { @$r = sort { 0+($a<=>$b) } @$r; }
  else            { @$r = sort {    $a<=>$b  } @$r; }
  return $r;
}

sub setbinop (&$;$) {   ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, $ra, $rb) = @_;
  croak 'Not a subroutine reference' unless (ref($sub) || '') eq 'CODE';
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY';
  if (defined $rb) {
    croak 'Not an array reference' unless (ref($rb) || '') eq 'ARRAY';
  } else {
    $rb = $ra;
  }

  my $caller = caller();
  no strict 'refs'; ## no critic(strict)
  local(*{$caller.'::a'}) = \my $a;
  local(*{$caller.'::b'}) = \my $b;

  # Typically faster and less memory to push them all instead of hashing here.
  my @set;
  for my $ia (@$ra) {
    $a = $ia;
    for my $ib (@$rb) {
      $b = $ib;
      push @set, $sub->();
    }
  }
  Mtoset(\@set);
}

sub sumset {
  my($ra,$rb) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY';
  if (defined $rb) {
    croak 'Not an array reference' unless (ref($rb) || '') eq 'ARRAY';
  } else {
    $rb = $ra;
  }
  return () if scalar(@$ra) == 0 || scalar(@$rb) == 0;

  validate_integer($_) for @$ra;
  if ($ra != $rb) { validate_integer($_) for @$rb; }

  my @set;
  for my $x (@$ra) {
    for my $y (@$rb) {
      push @set, Maddint($x,$y);
    }
  }
  Mtoset(\@set);
}

sub vecuniq {
  my %seen = ();
  my $k;
  # Validation means about 1.4x slower.
  #my @T = @_; return grep { validate_integer($_) && not $seen{$k = $_}++; } @T;
  # We have decided to skip validation and not support undefined values.
  return grep { not $seen{$k = $_}++; } @_;
}

sub vecfreq {
  my %count = ();
  my $countundef = 0;
  my $k;
  for (@_) {
    if (defined $_) { $count{$k = $_}++; } else { $countundef++; }
  }
  return wantarray ? %count : scalar(keys %count)   if !$countundef;
  return 1 + scalar(keys %count)                    if !wantarray;
  undef $k;
  return (%count, (\$k => $countundef));
}

sub vecsingleton {
  my %count = ();
  my ($countundef,$k) = (0);
  # Filter later duplicates during the count stage for a ~10% speedup.
  # Idea from List::MoreUtil.
  return grep { (defined $_ ? $count{$k=$_} : $countundef) == 1 }
         grep { ! (defined $_ ? $count{$k = $_}++ : $countundef++) }
         @_;
}

# SET/VEC generic.

sub setunion {
  my($ra,$rb) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY'
                                     && (ref($rb) || '') eq 'ARRAY';
  # return toset([@$ra,@$rb]);
  my(%seen,$k);
  my @res = Mtoset([grep { not $seen{$k = $_}++ } @$ra,@$rb]);
  return wantarray ? @res : scalar(@res);
}
sub setintersect {
  my($ra,$rb) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY'
                                     && (ref($rb) || '') eq 'ARRAY';
  ($ra,$rb) = ($rb,$ra) if scalar(@$ra) > scalar(@$rb);  # Performance
  if (scalar(@$ra) == 0) { return wantarray ? () : 0; }
  my %ina;
  undef @ina{@$ra};
  my @res = Mtoset([grep { exists $ina{$_} } @$rb]);
  return wantarray ? @res : scalar(@res);
}
sub setminus {
  my($ra,$rb) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY'
                                     && (ref($rb) || '') eq 'ARRAY';
  return @$ra if scalar(@$rb) == 0;
  my %inb;
  undef @inb{@$rb};
  my @res = Mtoset([grep { !exists $inb{$_} } @$ra]);
  return wantarray ? @res : scalar(@res);
}
sub setdelta {
  my($ra,$rb) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY'
                                     && (ref($rb) || '') eq 'ARRAY';
  return @$ra if scalar(@$rb) == 0;
  return @$rb if scalar(@$ra) == 0;
  my(%ina, %inb);
  undef @ina{@$ra};
  undef @inb{@$rb};
  my @s =  grep { !exists $inb{$_} } @$ra;
  push @s, grep { !exists $ina{$_} } @$rb;
  my @res =  Mtoset(\@s);
  return wantarray ? @res : scalar(@res);
}

# Can do setminus([$min..$max],\@L) albeit 2x slower
sub _setcomplement {
  my($ra, $min, $max) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY';
  validate_integer($min);
  validate_integer($max);
  my %ina;
  $ina{$_} = undef for @$ra;
  my @s;
  if ((ref($min) && !ref($max)) || (!ref($min) && ref($max))) {
    while ($min <= $max) {
      push @s, $min unless exists $ina{$min};
      $min = Maddint($min,1);
    }
  } else {
    while ($min <= $max) {
      push @s, $min unless exists $ina{$min};
      $min++;
    }
  }
  @s;
}

sub toset {
  my($ra) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY';

  validate_integer($_) for @$ra;
  return @$ra if scalar(@$ra) <= 1;
  my($k,%seen);
  Mvecsort( grep { not $seen{$k = $_}++; } @$ra );
}

# Is the second set a subset of the first set?
sub setcontains {
  my($set, $sub) = @_;
  my @newset;
  if (!ref($sub) || ref($sub) ne 'ARRAY') {
    validate_integer($sub);
    @newset = ($sub);
  } else {
    @newset = Mtoset($sub);
  }
  return 1 if @newset == 0;
  return 0 if @$set == 0 || @newset > @$set || $newset[-1] > $set->[-1] || $newset[0] < $set->[0];

  if (@$set <= 150 || (@$set <= 250 && @newset > 2)) {   # Linear search
    for my $sv (@$set) {
      if ($sv >= $newset[0]) {
        return 0 if $sv > $newset[0];
        shift @newset;
        return 1 if @newset == 0;
      }
    }
    return 0;
  }

  my $newlo = 0;
  # The next value is probably in this range.  Can save a lot of steps.
  my $range = Mcdivint(scalar(@$set),(scalar(@newset)+1) >> 1);
  for my $v (@newset) {
    my($lo,$hi) = ($newlo,$#$set);
    $hi = $lo + $range if $hi-$lo > $range && $set->[$lo+$range] >= $v;
    while ($lo < $hi) {
      my $mid = $lo + (($hi-$lo) >> 1);
      if ($set->[$mid] < $v) { $lo = $mid+1; }
      else                   { $hi = $mid; }
    }
    return 0 if $set->[$hi] != $v;
    $newlo = $hi+1;
  }
  1;
}

# Are any members of sub also a member of set?
sub setcontainsany {
  my($set, $sub) = @_;
  $sub = [$sub] if !ref($sub) || ref($sub) ne 'ARRAY';
  croak 'Not an array reference' unless (ref($set) || '') eq 'ARRAY'
                                     && (ref($sub) || '') eq 'ARRAY';

  # For better performance, make sub the larger
  ($set,$sub) = ($sub,$set) if scalar(@$set) > scalar(@$sub);
  return 0 if @$set == 0;

  my %ina;
  undef @ina{@$set};
  return 0 + (scalar(grep { exists $ina{$_} } @$sub) > 0);
}

sub _setinsert1 {
  my($rset, $v) = @_;
  validate_integer($v);

  if (scalar(@$rset) == 0 || $v > $rset->[-1]) {
    push @$rset, $v;
  } elsif ($v < $rset->[0]) {
    unshift @$rset, $v;
  } elsif (scalar(@$rset) > 1) {
    my($lo,$hi) = (0,$#$rset);
    while ($lo < $hi) {
      my $mid = $lo + (($hi-$lo) >> 1);
      if ($rset->[$mid] < $v) { $lo = $mid+1; }
      else                    { $hi = $mid; }
    }
    return 0 if $rset->[$hi] == $v;
    croak "internal too high" if $hi > 0 && $rset->[$hi-1] >= $v;
    croak "internal too low"  if $rset->[$hi] <= $v;
    splice @$rset, $hi, 0, $v;
  } else {
    return 0;   # Single element already in list.
  }
  1;
}

sub setinsert {
  my($set, $in) = @_;
  my @newset;
  if (!ref($in) || ref($in) ne 'ARRAY') {
    validate_integer($in);
    @newset = ($in);
  } else {
    @newset = Mtoset($in);
  }
  my $setsize = scalar(@$set);
  if ($setsize == 0 || $newset[0] > $set->[-1]) {
    push @$set, @newset;
  } elsif ($newset[-1] < $set->[0]) {
    unshift @$set, @newset;
  } elsif (@newset > 100) {
    @$set = Msetunion($set,\@newset);
  } else {
    # 1. values in front and back.
    my($nbeg,$nend) = (0,0);
    $nend++ while $nend < scalar(@newset) && $newset[-1 - $nend] > $set->[-1];
    push @$set, splice(@newset,-$nend) if $nend > 0;
    $nbeg++ while $nbeg < scalar(@newset) && $newset[$nbeg] < $set->[0];
    unshift @$set, splice(@newset,0,$nbeg) if $nbeg > 0;
    # 2. values in the middle.
    my $start = 0;
    my $range = Mcdivint(scalar(@$set),(scalar(@newset)+2) >> 1);
    for my $v (@newset) {
      my($lo,$hi) = ($start,$#$set);
      $hi = $lo + $range if $hi-$lo > $range && $set->[$lo+$range] >= $v;
      while ($lo < $hi) {
        my $mid = $lo + (($hi-$lo) >> 1);
        if ($set->[$mid] < $v) { $lo = $mid+1; }
        else                   { $hi = $mid; }
      }
      splice @$set, $hi, 0, $v if $set->[$hi] != $v;
      $start = $hi+1;
    }
  }
  return scalar(@$set) - $setsize;
}

sub _setremove1 {
  my($rset, $v) = @_;
  validate_integer($v);

  return 0 if scalar(@$rset) == 0 || $v > $rset->[-1] || $v < $rset->[0];

  my($lo,$hi) = (0,$#$rset);
  while ($lo < $hi) {
    my $mid = $lo + (($hi-$lo) >> 1);
    if ($rset->[$mid] < $v) { $lo = $mid+1; }
    else                    { $hi = $mid; }
  }
  return 0 if $rset->[$hi] != $v;
  splice @$rset, $hi, 1;
  1;
}

sub setremove {
  my($set, $in) = @_;
  my @newset;
  if (!ref($in) || ref($in) ne 'ARRAY') {
    validate_integer($in);
    @newset = ($in);
  } else {
    @newset = Mtoset($in);
  }
  my $setsize = scalar(@$set);
  return 0 if $setsize == 0 || scalar(@newset) == 0 || $newset[0] > $set->[-1] || $newset[-1] < $set->[0];

  if (@newset > 100) {
    @$set = Math::Prime::Util::setminus($set,\@newset);
  } else {
    #_setremove1($set, $_) for @newset;
    for my $v (@newset) {
      next if $v < $set->[0];
      last if $v > $set->[-1];
      my($lo,$hi) = (0,$#$set);
      while ($lo < $hi) {
        my $mid = $lo + (($hi-$lo) >> 1);
        if ($set->[$mid] < $v) { $lo = $mid+1; }
        else                   { $hi = $mid; }
      }
      if ($set->[$hi] == $v) {
        splice @$set, $hi, 1;
        last if @$set == 0;
      }
    }
  }
  return $setsize - scalar(@$set);
}

sub setinvert {
  my($set, $in) = @_;
  my @newset;
  if (!ref($in) || ref($in) ne 'ARRAY') {
    validate_integer($in);
    @newset = ($in);
  } else {
    @newset = Mtoset($in);
    return 0 if scalar(@newset) == 0;
  }
  # newset is in integer set form.  No duplicates.
  my $setsize = scalar(@$set);
  if ($setsize == 0) {
    @$set = @newset;
    return scalar(@$set);
  }
  if (@newset > 100) {
    @$set = Math::Prime::Util::setdelta($set,\@newset);
  } else {
    for my $v (@newset) {
      #_setremove1($set,$v) or _setinsert1($set,$v);
      my($lo,$hi) = (0,$#$set);
      while ($lo < $hi) {
        my $mid = $lo + (($hi-$lo) >> 1);
        if ($set->[$mid] < $v) { $lo = $mid+1; }
        else                   { $hi = $mid; }
      }
      if ($set->[$hi] == $v) {
        splice @$set, $hi, 1;
      } else {
        splice @$set, $hi, 0, $v;
      }
    }
  }
  return scalar(@$set) - $setsize;
}

sub set_is_disjoint {
  my($s,$t) = @_;
  croak 'Not an array reference' unless (ref($s) || '') eq 'ARRAY'
                                     && (ref($t) || '') eq 'ARRAY';
  return 0 + (scalar(Msetintersect($s,$t) == 0));
}
sub set_is_equal {
  my($s,$t) = @_;
  croak 'Not an array reference' unless (ref($s) || '') eq 'ARRAY'
                                     && (ref($t) || '') eq 'ARRAY';
  return 0 + (@$s == @$t && scalar(Msetintersect($s,$t)) == @$t);
}
sub set_is_subset {
  my($s,$t) = @_;
  croak 'Not an array reference' unless (ref($s) || '') eq 'ARRAY'
                                     && (ref($t) || '') eq 'ARRAY';
  return 0 + (@$s >= @$t && scalar(Msetintersect($s,$t)) == @$t);
}
sub set_is_proper_subset {
  my($s,$t) = @_;
  croak 'Not an array reference' unless (ref($s) || '') eq 'ARRAY'
                                     && (ref($t) || '') eq 'ARRAY';
  return 0 + (@$s > @$t && scalar(Msetintersect($s,$t)) == @$t);
}
sub set_is_superset {
  my($s,$t) = @_;
  croak 'Not an array reference' unless (ref($s) || '') eq 'ARRAY'
                                     && (ref($t) || '') eq 'ARRAY';
  return 0 + (@$s <= @$t && scalar(Msetintersect($s,$t)) == @$s);
}
sub set_is_proper_superset {
  my($s,$t) = @_;
  croak 'Not an array reference' unless (ref($s) || '') eq 'ARRAY'
                                     && (ref($t) || '') eq 'ARRAY';
  return 0 + (@$s < @$t && scalar(Msetintersect($s,$t)) == @$s);
}
sub set_is_proper_intersection {
  my($s,$t) = @_;
  croak 'Not an array reference' unless (ref($s) || '') eq 'ARRAY'
                                     && (ref($t) || '') eq 'ARRAY';
  my $minsize = (scalar(@$s) < scalar(@$t)) ? scalar(@$s) : scalar(@$t);
  my $intersize = scalar(Msetintersect($s,$t));
  return ($intersize > 0 && $intersize < $minsize) ? 1 : 0;
}

sub is_sidon_set {
  my($ra) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY';

  my %sums;
  my @S = Mtoset($ra);  # Validated, sorted, deduped.
  while (@S) {
    my $x = pop @S;
    return 0 if $x < 0;
    for my $y ($x, @S) {
      my $s = Maddint($x, $y);
      return 0 if exists $sums{$s};
      $sums{$s} = undef;
    }
  }
  1;
}

sub is_sumfree_set {
  my($ra) = @_;
  croak 'Not an array reference' unless (ref($ra) || '') eq 'ARRAY';

  my %ina;
  my @S = Mtoset($ra);  # Validated, sorted, deduped.
  $ina{$_}=undef for @S;
  while (@S) {
    my $x = pop @S;
    for my $y ($x, @S) {
      return 0 if exists $ina{Maddint($x,$y)};
    }
  }
  1;
}


###############################################################################

sub foralmostprimes {
  my($sub, $k, $lo, $hi) = @_;
  validate_integer_nonneg($k);
  return if $k == 0;
  if (defined $hi) { validate_integer_nonneg($lo); }
  else             { ($lo,$hi) = (1, $lo);         }
  validate_integer_nonneg($hi);

  $lo = Mvecmax($lo, Mpowint(2, $k));
  return if $lo > $hi;

  #return Math::Prime::Util::forprimes($sub,$lo,$hi) if $k == 1;

  my $estcount = almost_prime_count_approx($k,$hi) - almost_prime_count_approx($k,$lo);
  my $nsegs = "$estcount" / 1e6;
  my $len = Maddint(Msubint($hi,$lo),1);
  my $segsize = ($nsegs <= 1.1) ? $len : int("$len"/$nsegs);
  if ($segsize < 5*1e6) { $segsize = 5e6; }
  # warn "  estcount $estcount   nsegs $nsegs   segsize $segsize\n";

  my $oldforexit = Math::Prime::Util::_start_for_loop();
  while ($lo <= $hi) {
    my $seghi = Mvecmin($hi, Maddint($lo,$segsize)-1);
    my $ap = Math::Prime::Util::almost_primes($k, $lo, $seghi);
    #my $ap = [];  _genkap($lo, $seghi, $k, 1, 2, sub { push @$ap,$_[0]; });
    # warn "  from $lo to $seghi found ",scalar(@$ap), " $k-almost-primes\n";
    {
      my $pp;
      local *_ = \$pp;
      for my $kap (@$ap) {
        $pp = $kap;
        $sub->();
        last if Math::Prime::Util::_get_forexit();
      }
    }
    $lo = Maddint($seghi,1);
    last if Math::Prime::Util::_get_forexit();
  }
  Math::Prime::Util::_end_for_loop($oldforexit);
}



###############################################################################
#       Random numbers
###############################################################################

# PPFE:  irand irand64 drand random_bytes csrand srand _is_csprng_well_seeded
sub urandomb {
  my($n) = @_;
  return 0 if $n <= 0;
  return ( Math::Prime::Util::irand() >> (32-$n) ) if $n <= 32;
  return ( Math::Prime::Util::irand64() >> (64-$n) ) if MPU_MAXBITS >= 64 && $n <= 64;
  my $bytes = Math::Prime::Util::random_bytes(($n+7)>>3);
  return _frombinary( substr(unpack("B*",$bytes),0,$n) );
}
sub urandomm {
  my($n) = @_;
  # validate_integer_nonneg($n);
  return reftyped($_[0], Math::Prime::Util::GMP::urandomm($n))
    if $Math::Prime::Util::_GMPfunc{"urandomm"};
  return 0 if $n <= 1;
  my $r;
  if ($n <= 4294967295) {
    my $rmin = (4294967295 - ($n-1)) % $n;
    do { $r = Math::Prime::Util::irand(); } while $r < $rmin;
  } elsif (!ref($n)) {
    my $rmin = (~0 - ($n-1)) % $n;
    do { $r = Math::Prime::Util::irand64(); } while $r < $rmin;
  } else {
    # TODO: verify and try to optimize this
    my $bytes = 1 + length(todigitstring($n,16));
    my $rmax = Msubint(Mpowint(2,$bytes*8),1);
    my $overflow = $rmax - ($rmax % $n);
    do { $r = Murandomb($bytes*8); } while $r >= $overflow;
  }
  return $r % $n;
}

sub random_prime {
  my($low, $high) = @_;
  if (scalar(@_) == 1) { ($low,$high) = (2,$low);       }
  else                 { validate_integer_nonneg($low); }
  validate_integer_nonneg($high);

  return reftyped($_[0], Math::Prime::Util::GMP::random_prime($low, $high))
    if $Math::Prime::Util::_GMPfunc{"random_prime"};

  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_prime($low,$high);
}

sub random_ndigit_prime {
  my($digits) = @_;
  validate_integer_nonneg($digits);
  croak "random_ndigit_prime digits must be >= 1" unless $digits >= 1;
  return reftyped($_[0], Math::Prime::Util::GMP::random_ndigit_prime($digits))
    if $Math::Prime::Util::_GMPfunc{"random_ndigit_prime"} && !getconfig()->{'nobigint'};
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_ndigit_prime($digits);
}
sub random_nbit_prime {
  my($bits) = @_;
  validate_integer_nonneg($bits);
  croak "random_nbit_prime bits must be >= 2" unless $bits >= 2;
  return reftyped($_[0], Math::Prime::Util::GMP::random_nbit_prime($bits))
    if $Math::Prime::Util::_GMPfunc{"random_nbit_prime"};
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_nbit_prime($bits);
}
sub random_safe_prime {
  my($bits) = @_;
  validate_integer_nonneg($bits);
  croak "random_safe_prime bits must be >= 3" unless $bits >= 3;
  return reftyped($_[0], eval "Math::Prime::Util::GMP::random_safe_prime($bits)")  ## no critic qw(ProhibitStringyEval)
    if $Math::Prime::Util::_GMPfunc{"random_safe_prime"};
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_safe_prime($bits);
}
sub random_strong_prime {
  my($bits) = @_;
  validate_integer_nonneg($bits);
  croak "random_strong_prime bits must be >= 128" unless $bits >= 128;
  return reftyped($_[0], eval "Math::Prime::Util::GMP::random_strong_prime($bits)")  ## no critic qw(ProhibitStringyEval)
    if $Math::Prime::Util::_GMPfunc{"random_strong_prime"};
  require Math::Prime::Util::RandomPrimes;
  return Math::Prime::Util::RandomPrimes::random_strong_prime($bits);
}

sub random_proven_prime {
  random_maurer_prime(@_);
}

sub random_maurer_prime {
  my($bits) = @_;
  validate_integer_nonneg($bits);
  croak "random_maurer_prime bits must be >= 2" unless $bits >= 2;

  return reftyped($_[0], Math::Prime::Util::GMP::random_maurer_prime($bits))
    if $Math::Prime::Util::_GMPfunc{"random_maurer_prime"};

  require Math::Prime::Util::RandomPrimes;
  my ($n, $cert) = Math::Prime::Util::RandomPrimes::random_maurer_prime_with_cert($bits);
  croak "maurer prime $n failed certificate verification!"
        unless Math::Prime::Util::verify_prime($cert);

  return $n;
}

sub random_shawe_taylor_prime {
  my($bits) = @_;
  validate_integer_nonneg($bits);
  croak "random_shawe_taylor_prime bits must be >= 2" unless $bits >= 2;

  return reftyped($_[0], Math::Prime::Util::GMP::random_shawe_taylor_prime($bits))
    if $Math::Prime::Util::_GMPfunc{"random_shawe_taylor_prime"};

  require Math::Prime::Util::RandomPrimes;
  my ($n, $cert) = Math::Prime::Util::RandomPrimes::random_shawe_taylor_prime_with_cert($bits);
  croak "shawe-taylor prime $n failed certificate verification!"
        unless Math::Prime::Util::verify_prime($cert);

  return $n;
}

sub miller_rabin_random {
  my($n, $k, $seed) = @_;
  validate_integer_nonneg($n);
  if (scalar(@_) == 1 ) { $k = 1; } else { validate_integer_nonneg($k); }

  return 1 if $k <= 0;

  if ($Math::Prime::Util::_GMPfunc{"miller_rabin_random"}) {
    return Math::Prime::Util::GMP::miller_rabin_random($n, $k, $seed) if defined $seed;
    return Math::Prime::Util::GMP::miller_rabin_random($n, $k);
  }

  # getconfig()->{'assume_rh'})  ==>  2*log(n)^2
  if ($k >= int(3*$n/4) ) {
    for (2 .. int(3*$n/4)+2) {
      return 0 unless Math::Prime::Util::is_strong_pseudoprime($n, $_);
    }
    return 1;
  }
  my $brange = $n-2;
  return 0 unless Math::Prime::Util::is_strong_pseudoprime($n, Murandomm($brange)+2 );
  $k--;
  while ($k > 0) {
    my $nbases = ($k >= 20) ? 20 : $k;
    return 0 unless is_strong_pseudoprime($n, map { Murandomm($brange)+2 } 1 .. $nbases);
    $k -= $nbases;
  }
  1;
}

sub random_semiprime {
  my($b) = @_;
  validate_integer_nonneg($b);
  croak "random_semiprime bits must be >= 4" unless $b >= 4;

  my $n;
  my $min = Mpowint(2,$b-1);
  my $max = $min + ($min - 1);
  my $L = $b >> 1;
  my $N = $b - $L;
  do {
    $n = Mmulint(random_nbit_prime($L), random_nbit_prime($N));
  } while $n < $min || $n > $max;
  $n;
}

sub random_unrestricted_semiprime {
  my($b) = @_;
  validate_integer_nonneg($b);
  croak "random_unrestricted_semiprime bits must be >= 3" unless $b >= 3;

  my $n;
  my $min = Mpowint(2,$b-1);
  my $max = Maddint($min, $min - 1);

  if ($b <= MPU_MAXBITS) {
    do {
      $n = $min + Murandomb($b-1);
    } while !Mis_semiprime($n);
  } else {
    # Try to get probabilities right for small divisors
    my %M = (
      2 => 1.91218397452243,
      3 => 1.33954826555021,
      5 => 0.854756717114822,
      7 => 0.635492301836862,
      11 => 0.426616792046787,
      13 => 0.368193843118344,
      17 => 0.290512701603111,
      19 => 0.263359264658156,
      23 => 0.222406328935102,
      29 => 0.181229250520242,
      31 => 0.170874199059434,
      37 => 0.146112155735473,
      41 => 0.133427839963585,
      43 => 0.127929010905662,
      47 => 0.118254609086782,
      53 => 0.106316418106489,
      59 => 0.0966989675438643,
      61 => 0.0938833658008547,
      67 => 0.0864151823151671,
      71 => 0.0820822953188297,
      73 => 0.0800964416340746,
      79 => 0.0747060914833344,
      83 => 0.0714973706654851,
      89 => 0.0672115468436284,
      97 => 0.0622818892486191,
      101 => 0.0600855891549939,
      103 => 0.0590613570015407,
      107 => 0.0570921135626976,
      109 => 0.0561691667641485,
      113 => 0.0544330141081874,
      127 => 0.0490620204315701,
    );
    my ($p,$r);
    $r = Math::Prime::Util::drand();
    for my $prime (2..113,127) {
      next unless defined $M{$prime};
      my $PR = $M{$prime} / $b  +  0.19556 / $prime;
      if ($r <= $PR) {
        $p = $prime;
        last;
      }
      $r -= $PR;
    }
    if (!defined $p) {
      # Idea from Charles Greathouse IV, 2010.  The distribution is right
      # at the high level (small primes weighted more and not far off what
      # we get with the uniform selection), but there is a noticeable skew
      # toward primes with a large gap after them.  For instance 3 ends up
      # being weighted as much as 2, and 7 more than 5.
      #
      # Since we handled small divisors earlier, this is less bothersome.
      my $M = 0.26149721284764278375542683860869585905;
      my $weight = $M + log($b * log(2)/2);
      my $minr = log(log(131));
      do {
        $r  = Math::Prime::Util::drand($weight) - $M;
      } while $r < $minr;
      my $a;
      if ($r <= 3.54) {
        # result under 10^15, can do directly
        $a = int( exp(exp($r)) + 0.5 );
      } elsif ($Math::Prime::Util::_GMPfunc{"expreal"}) {
        # Use our fast arbitrary precision expreal.
        my $digits = $r < 4.45 ? 40 : int(exp($r)/2.2 + 2);  # overestimate
        my $re = Math::Prime::Util::GMP::expreal($r,$digits);
        $a = Math::Prime::Util::GMP::expreal($re,$digits);
        $a = Mtoint($a);  #_upgrade_to_float($a)->as_int;
      } else {
        # exp(x)=exp(x/n)^n
        # We could use Math::BigFloat but it's sooooooooooo slow.
        my $re = exp($r);
        my $redd = 1+int($re/34.5);
        $a = Mpowint(int(exp($re/$redd)+0.5), $redd);
      }
      $p = $a < 2 ? 2 : Mprev_prime($a+1);
    }
    my $ranmin = Mcdivint($min, $p);
    my $ranmax = Mdivint($max, $p);
    my $q = random_prime($ranmin, $ranmax);
    $n = Mmulint($p,$q);
  }
  $n;
}

sub random_factored_integer {
  my($n) = @_;
  validate_integer_positive($n);

  while (1) {
    my @S = ($n);
    # make s_i chain
    push @S, 1 + Murandomm($S[-1])  while $S[-1] > 1;
    # first is n, last is 1
    @S = grep { Mis_prime($_) } @S[1 .. $#S-1];
    my $r = Mvecprod(@S);
    return ($r, [@S]) if $r <= $n && (1+Murandomm($n)) <= $r;
  }
}



1;

__END__


# ABSTRACT: Pure Perl version of Math::Prime::Util

=pod

=encoding utf8


=head1 NAME

Math::Prime::Util::PP - Pure Perl version of Math::Prime::Util


=head1 VERSION

Version 0.73


=head1 SYNOPSIS

The functionality is basically identical to L<Math::Prime::Util>, as this
module is just the Pure Perl implementation.  This documentation will only
note differences.

  # Normally you would just import the functions you are using.
  # Nothing is exported by default.
  use Math::Prime::Util ':all';


=head1 DESCRIPTION

Pure Perl implementations of prime number utilities that are normally
handled with XS or GMP.  Having the Perl implementations (1) provides examples,
(2) allows the functions to run even if XS isn't available, and (3) gives
big number support if L<Math::Prime::Util::GMP> isn't available.

All routines should work with native integers or multi-precision numbers.  To
enable big numbers, use bigint:

    use bigint;
    say prime_count_approx(1000000000000000000000000)'
    # says 18435599767347541878147

Or string inputs:

    say prime_count_approx("1000000000000000000000000")'
    # identical output.

Some functions will be very slow.
L<Math::Prime::Util::GMP> has much faster versions of many of these functions.
Alternately, L<Math::Pari> has a lot of these types of functions.


=head1 FUNCTIONS

=head2 euler_phi

Takes a I<single> integer input and returns the Euler totient.

=head2 euler_phi_range

Takes two values defining a range C<low> to C<high> and returns an array
with the totient of each value in the range, inclusive.

=head2 moebius

Takes a I<single> integer input and returns the Moebius function.

=head2 moebius_range

Takes two values defining a range C<low> to C<high> and returns an array
with the Moebius function of each value in the range, inclusive.


=head1 LIMITATIONS

The SQUFOF and Fermat factoring algorithms are not implemented yet.

Some of the prime methods use more memory than they should, as the segmented
sieve is not properly used in C<primes> and C<prime_count>.


=head1 PERFORMANCE

Performance compared to the XS/C code is quite poor for many operations.  Some
operations that are relatively close for small and medium-size values:

  next_prime / prev_prime
  is_prime / is_prob_prime
  is_strong_pseudoprime
  ExponentialIntegral / LogarithmicIntegral / RiemannR
  primearray

Operations that are slower include:

  primes
  random_prime / random_ndigit_prime
  factor / factor_exp / divisors
  nth_prime
  prime_count
  is_aks_prime

Performance improvement in this code is still possible.  The prime sieve is
over 2x faster than anything I was able to find online, but it is still has
room for improvement.

L<Math::Prime::Util::GMP> offers C<C+XS+GMP> support for most of the important
functions, and will be vastly faster for most operations.  If you install that
module, L<Math::Prime::Util> will load it automatically, meaning you should
not have to think about what code is actually being used (C, GMP, or Perl).

Memory use will generally be higher for the PP code, and in some cases B<much>
higher.  Some of this may be addressed in a later release.

For small values (e.g. primes and prime counts under 10M) most of this will
not matter.


=head1 SEE ALSO

L<Math::Prime::Util>

L<Math::Prime::Util::GMP>


=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012-2026 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
