# The PP front end, only loaded if XS is not used.  It is intended to load
# directly into the package namespace.

use strict;
use warnings;
use Math::Prime::Util::PP;

*_validate_num = \&Math::Prime::Util::PP::_validate_num;
*_prime_memfreeall = \&Math::Prime::Util::PP::_prime_memfreeall;
*prime_memfree  = \&Math::Prime::Util::PP::prime_memfree;
*prime_precalc  = \&Math::Prime::Util::PP::prime_precalc;


sub mertens {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::mertens(@_);
}

sub moebius {
  if (scalar @_ <= 1) {
    my($n) = @_;
    return 0 if defined $n && $n < 0;
    _validate_num($n) || _validate_positive_integer($n);
    return Math::Prime::Util::PP::moebius($n);
  }
  my($lo, $hi) = @_;
  _validate_num($lo) || _validate_positive_integer($lo);
  _validate_num($hi) || _validate_positive_integer($hi);
  return Math::Prime::Util::PP::moebius_range($lo, $hi);
}

sub euler_phi {
  if (scalar @_ <= 1) {
    my($n) = @_;
    return 0 if defined $n && $n < 0;
    _validate_num($n) || _validate_positive_integer($n);
    return Math::Prime::Util::PP::euler_phi($n);
  }
  my($lo, $hi) = @_;
  _validate_num($lo) || _validate_positive_integer($lo);
  _validate_num($hi) || _validate_positive_integer($hi);
  return Math::Prime::Util::PP::euler_phi_range($lo, $hi);
}
sub jordan_totient {
  my($k, $n) = @_;
  _validate_positive_integer($k);
  return 0 if defined $n && $n < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::jordan_totient(@_);
}
sub carmichael_lambda {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::carmichael_lambda(@_);
}

sub nth_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_prime(@_);
}
sub nth_prime_lower {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_prime_lower($n);
}
sub nth_prime_upper {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_prime_upper($n);
}
sub nth_prime_approx {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_prime_approx($n);
}
sub prime_count_lower {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::prime_count_lower($n);
}
sub prime_count_upper {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::prime_count_upper($n);
}
sub prime_count_approx {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::prime_count_approx($n);
}
  

sub is_prime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_prime($n);
}
sub is_prob_prime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_prob_prime($n);
}
sub is_pseudoprime {
  my($n, $base) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  _validate_positive_integer($base);
  return Math::Prime::Util::PP::is_pseudoprime($n, $base);
}
sub is_strong_pseudoprime {
  my($n, @bases) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  croak "No bases given to miller_rabin" unless @bases;
  return Math::Prime::Util::PP::is_strong_pseudoprime($n, @bases);
}
sub is_lucas_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_lucas_pseudoprime($n);
}
sub is_strong_lucas_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_strong_lucas_pseudoprime($n);
}
sub is_extra_strong_lucas_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_extra_strong_lucas_pseudoprime($n);
}
sub is_almost_extra_strong_lucas_pseudoprime {
  my($n, $increment) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  if (defined $increment) { _validate_positive_integer($increment, 1, 256);
  } else                  { $increment = 1; }
  return Math::Prime::Util::PP::is_almost_extra_strong_lucas_pseudoprime($n, $increment);
}
sub is_frobenius_underwood_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_frobenius_underwood_pseudoprime($n);
}
sub is_aks_prime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_aks_prime($n);
}


sub kronecker {
  my($a, $b) = @_;
  my ($va, $vb) = ($a, $b);
  $va = -$va if defined $va && $va < 0;
  $vb = -$vb if defined $vb && $vb < 0;
  _validate_positive_integer($va);
  _validate_positive_integer($vb);
  return Math::Prime::Util::PP::kronecker(@_);
}

sub znorder {
  my($a, $n) = @_;
  _validate_positive_integer($a);
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::znorder($a, $n);
}

sub znlog {
  my($a, $g, $p) = @_;
  _validate_positive_integer($a);
  _validate_positive_integer($g);
  _validate_positive_integer($p);
  return Math::Prime::Util::PP::znlog($a, $g, $p);
}

sub znprimroot {
  my($n) = @_;
  $n = -$n if defined $n && $n =~ /^-\d+/;   # TODO: fix this for string bigints
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::znprimroot($n);
}

sub trial_factor {
  my($n, $maxlim) = @_;
  _validate_positive_integer($n);
  if (defined $maxlim) {
    _validate_positive_integer($maxlim);
    return Math::Prime::Util::PP::trial_factor($n, $maxlim);
  }
  return Math::Prime::Util::PP::trial_factor($n);
}
sub fermat_factor {
  my($n, $rounds) = @_;
  _validate_positive_integer($n);
  if (defined $rounds) {
    _validate_positive_integer($rounds);
    return Math::Prime::Util::PP::fermat_factor($n, $rounds);
  }
  return Math::Prime::Util::PP::fermat_factor($n);
}
sub holf_factor {
  my($n, $rounds) = @_;
  _validate_positive_integer($n);
  if (defined $rounds) {
    _validate_positive_integer($rounds);
    return Math::Prime::Util::PP::holf_factor($n, $rounds);
  }
  return Math::Prime::Util::PP::holf_factor($n);
}
sub squfof_factor {
  my($n, $rounds) = @_;
  _validate_positive_integer($n);
  if (defined $rounds) {
    _validate_positive_integer($rounds);
    return Math::Prime::Util::PP::squfof_factor($n, $rounds);
  }
  return Math::Prime::Util::PP::squfof_factor($n);
}
sub pbrent_factor {
  my($n, $rounds, $pa) = @_;
  _validate_positive_integer($n);
  if (defined $rounds) { _validate_positive_integer($rounds);
  } else               { $rounds = 4*1024*1024; }
  if (defined $pa    ) { _validate_positive_integer($pa);
  } else               { $pa = 3; }
  return Math::Prime::Util::PP::pbrent_factor($n, $rounds, $pa);
}
sub prho_factor {
  my($n, $rounds, $pa) = @_;
  _validate_positive_integer($n);
  if (defined $rounds) { _validate_positive_integer($rounds);
  } else               { $rounds = 4*1024*1024; }
  if (defined $pa    ) { _validate_positive_integer($pa);
  } else               { $pa = 3; }
  return Math::Prime::Util::PP::prho_factor($n, $rounds, $pa);
}
sub pminus1_factor {
  my($n, $B1, $B2) = @_;
  _validate_positive_integer($n);
  _validate_positive_integer($B1) if defined $B1;
  _validate_positive_integer($B2) if defined $B2;
  Math::Prime::Util::PP::pminus1_factor($n, $B1, $B2);
}
*pplus1_factor = \&pminus1_factor;
sub ecm_factor {
  my($n, $B1, $B2, $ncurves) = @_;
  _validate_positive_integer($n);
  _validate_positive_integer($B1) if defined $B1;
  _validate_positive_integer($B2) if defined $B2;
  _validate_positive_integer($ncurves) if defined $ncurves;
  Math::Prime::Util::PP::ecm_factor($n, $B1, $B2, $ncurves);
}

sub gcd {
  return Math::Prime::Util::PP::gcd(@_);
}
sub lcm {
  return Math::Prime::Util::PP::lcm(@_);
}

sub legendre_phi {
  my($x, $a) = @_;
  _validate_positive_integer($x);
  _validate_positive_integer($a);
  return Math::Prime::Util::PP::legendre_phi($x, $a);
}

1;
