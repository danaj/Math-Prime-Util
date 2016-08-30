package Math::Prime::Util::PPFE;
use strict;
use warnings;
use Math::Prime::Util::PP;

# The PP front end, only loaded if XS is not used.
# It is intended to load directly into the MPU namespace.

package Math::Prime::Util;
use Carp qw/carp croak confess/;

*_validate_num = \&Math::Prime::Util::PP::_validate_num;
*_validate_integer = \&Math::Prime::Util::PP::_validate_integer;
*_prime_memfreeall = \&Math::Prime::Util::PP::_prime_memfreeall;
*prime_memfree  = \&Math::Prime::Util::PP::prime_memfree;
*prime_precalc  = \&Math::Prime::Util::PP::prime_precalc;


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
  return Math::Prime::Util::PP::jordan_totient($k, $n);
}
sub ramanujan_sum {
  my($k,$n) = @_;
  _validate_positive_integer($k);
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::ramanujan_sum($k, $n);
}
sub carmichael_lambda {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::carmichael_lambda($n);
}
sub mertens {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::mertens($n);
}
sub liouville {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::liouville($n);
}
sub exp_mangoldt {
  my($n) = @_;
  return 1 if defined $n && $n <= 1;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::exp_mangoldt($n);
}
sub hclassno {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::hclassno($n);
}


sub next_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::next_prime($n);
}
sub prev_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::prev_prime($n);
}
sub nth_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_prime($n);
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
sub twin_prime_count {
  my($low,$high) = @_;
  return Math::Prime::Util::PP::twin_prime_count($low,$high);
}
sub ramanujan_prime_count {
  my($low,$high) = @_;
  return Math::Prime::Util::PP::ramanujan_prime_count($low,$high);
}
sub sum_primes {
  my($low,$high) = @_;
  return Math::Prime::Util::PP::sum_primes($low,$high);
}
sub print_primes {
  my($low,$high,$fd) = @_;
  return Math::Prime::Util::PP::print_primes($low,$high,$fd);
}
sub sieve_range {
  my($n,$width,$depth) = @_;
  return Math::Prime::Util::PP::sieve_range($n,$width,$depth);
}
sub sieve_prime_cluster {
  my($low,$high,@cl) = @_;
  return Math::Prime::Util::PP::sieve_prime_cluster($low,$high,@cl);
}
sub twin_prime_count_approx {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::twin_prime_count_approx($n);
}
sub nth_twin_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_twin_prime($n);
}
sub nth_twin_prime_approx {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_twin_prime_approx($n);
}
sub nth_ramanujan_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::nth_ramanujan_prime($n);
}


*is_prime          = \&Math::Prime::Util::PP::is_prime;
*is_prob_prime     = \&Math::Prime::Util::PP::is_prob_prime;
*is_provable_prime = \&Math::Prime::Util::PP::is_provable_prime;
*is_bpsw_prime     = \&Math::Prime::Util::PP::is_bpsw_prime;

sub is_pseudoprime {
  my($n, @bases) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  croak "No bases given to is_strong_pseudoprime" unless @bases;
  return Math::Prime::Util::PP::is_pseudoprime($n, @bases);
}
sub is_euler_pseudoprime {
  my($n, @bases) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  croak "No bases given to is_euler_pseudoprime" unless @bases;
  return Math::Prime::Util::PP::is_euler_pseudoprime($n, @bases);
}
sub is_strong_pseudoprime {
  my($n, @bases) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  croak "No bases given to is_strong_pseudoprime" unless @bases;
  return Math::Prime::Util::PP::is_strong_pseudoprime($n, @bases);
}
sub is_euler_plumb_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_euler_plumb_pseudoprime($n);
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
sub is_perrin_pseudoprime {
  my($n,$restrict) = @_;
  return 0 if defined $n && int($n) < 0;
  $restrict = 0 unless defined $restrict;
  _validate_positive_integer($n);
  _validate_positive_integer($restrict);
  return Math::Prime::Util::PP::is_perrin_pseudoprime($n, $restrict);
}
sub is_catalan_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_catalan_pseudoprime($n);
}
sub is_frobenius_pseudoprime {
  my($n, $P, $Q) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  # TODO: validate P & Q
  return Math::Prime::Util::PP::is_frobenius_pseudoprime($n, $P, $Q);
}
sub is_frobenius_underwood_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_frobenius_underwood_pseudoprime($n);
}
sub is_frobenius_khashin_pseudoprime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_frobenius_khashin_pseudoprime($n);
}
sub is_aks_prime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_aks_prime($n);
}
sub is_ramanujan_prime {
  my($n) = @_;
  return 0 if defined $n && int($n) < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_ramanujan_prime($n);
}
sub is_mersenne_prime {
  my($p) = @_;
  _validate_positive_integer($p);
  return Math::Prime::Util::PP::is_mersenne_prime($p);
}
sub is_square_free {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_square_free($n);
}
sub is_carmichael {
  my($n) = @_;
  #_validate_positive_integer($n);
  return Math::Prime::Util::PP::is_carmichael($n);
}
sub is_quasi_carmichael {
  my($n) = @_;
  #_validate_positive_integer($n);
  return Math::Prime::Util::PP::is_quasi_carmichael($n);
}
sub is_primitive_root {
  my($a,$n) = @_;
  $n = -$n if defined $n && $n < 0;
  $a %= $n if defined $a && $a < 0;
  _validate_positive_integer($a);
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::is_primitive_root($a,$n);
}


sub lucas_sequence {
  my($n, $P, $Q, $k) = @_;
  my ($vp, $vq) = ($P, $Q);
  $vp = -$vp if defined $vp && $vp < 0;
  $vq = -$vq if defined $vq && $vq < 0;
  _validate_positive_integer($n);
  _validate_positive_integer($vp);
  _validate_positive_integer($vq);
  _validate_positive_integer($k);
  return Math::Prime::Util::PP::lucas_sequence(@_);
}
sub lucasu {
  my($P, $Q, $k) = @_;
  my ($vp, $vq) = ($P, $Q);
  $vp = -$vp if defined $vp && $vp < 0;
  $vq = -$vq if defined $vq && $vq < 0;
  _validate_positive_integer($vp);
  _validate_positive_integer($vq);
  _validate_positive_integer($k);
  return Math::Prime::Util::PP::lucasu(@_);
}
sub lucasv {
  my($P, $Q, $k) = @_;
  my ($vp, $vq) = ($P, $Q);
  $vp = -$vp if defined $vp && $vp < 0;
  $vq = -$vq if defined $vq && $vq < 0;
  _validate_positive_integer($vp);
  _validate_positive_integer($vq);
  _validate_positive_integer($k);
  return Math::Prime::Util::PP::lucasv(@_);
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

sub factorial {
  my($n) = @_;
  _validate_integer($n);
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
  _validate_positive_integer($n);
  _validate_positive_integer($k);
  _validate_positive_integer($type) if defined $type;
  return Math::Prime::Util::PP::stirling($n, $k, $type);
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
  $n =~ s/^-(\d+)/$1/ if defined $n;
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

sub divisors {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::divisors($n);
}

sub divisor_sum {
  my($n, $k) = @_;
  _validate_positive_integer($n);
  _validate_positive_integer($k) if defined $k && ref($k) ne 'CODE';
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
sub chinese {
  # TODO: make sure we're not modding their data
  foreach my $aref (@_) {
    die "chinese arguments are two-element array references"
      unless ref($aref) eq 'ARRAY' && scalar @$aref == 2;
    _validate_integer($aref->[0]);
    _validate_integer($aref->[1]);
  }
  return Math::Prime::Util::PP::chinese(@_);
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
sub invmod {
  my ($a, $n) = @_;
  _validate_integer($a);
  _validate_integer($n);
  return Math::Prime::Util::PP::invmod($a,$n);
}
sub sqrtmod {
  my ($a, $n) = @_;
  _validate_integer($a);
  _validate_integer($n);
  return Math::Prime::Util::PP::sqrtmod($a,$n);
}
sub addmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b>=0?$b:-$b); _validate_integer($n);
  return Math::Prime::Util::PP::addmod($a,$b, $n);
}
sub mulmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b>=0?$b:-$b); _validate_integer($n);
  return Math::Prime::Util::PP::mulmod($a,$b, $n);
}
sub divmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b>=0?$b:-$b); _validate_integer($n);
  return Math::Prime::Util::PP::divmod($a,$b, $n);
}
sub powmod {
  my ($a, $b, $n) = @_;
  _validate_integer($a); _validate_integer($b>=0?$b:-$b); _validate_integer($n);
  return Math::Prime::Util::PP::powmod($a,$b, $n);
}
sub sqrtint {
  my($n) = @_;
  _validate_integer($n);
  return Math::Prime::Util::PP::sqrtint($n);
}
sub rootint {
  my($n, $k, $refp) = @_;
  _validate_positive_integer($n);
  _validate_positive_integer($k);
  return Math::Prime::Util::PP::rootint($n, $k, $refp);
}
sub logint {
  my($n, $b, $refp) = @_;
  _validate_positive_integer($n);
  _validate_positive_integer($b);
  return Math::Prime::Util::PP::logint($n, $b, $refp);
}

sub legendre_phi {
  my($x, $a) = @_;
  _validate_positive_integer($x);
  _validate_positive_integer($a);
  return Math::Prime::Util::PP::legendre_phi($x, $a);
}

sub chebyshev_theta {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::chebyshev_theta($n);
}
sub chebyshev_psi {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::chebyshev_psi($n);
}
sub ramanujan_tau {
  my($n) = @_;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::ramanujan_tau($n);
}

sub is_power {
  my($n, $a, $refp) = @_;
  my $vn = "$n";  $vn =~ s/^-//;
  _validate_positive_integer($vn);
  _validate_positive_integer($a) if defined $a;
  $vn = '-'.$vn if $n < 0;
  return Math::Prime::Util::PP::is_power($vn, $a, $refp);
}
sub is_prime_power {
  my($n, $refp) = @_;
  my $vn = "$n";  $vn =~ s/^-//;
  _validate_positive_integer($vn);
  $vn = '-'.$vn if $n < 0;
  return Math::Prime::Util::PP::is_prime_power($vn, $refp);
}
sub valuation {
  my($n, $k) = @_;
  $n = -$n if defined $n && $n < 0;
  $k = -$k if defined $k && $k < 0;
  _validate_positive_integer($n);
  _validate_positive_integer($k);
  return Math::Prime::Util::PP::valuation($n, $k);
}
sub todigits {
  my($n,$base,$len) = @_;
  return Math::Prime::Util::PP::todigits($n,$base,$len);
}
sub todigitstring {
  my($n,$base,$len) = @_;
  return Math::Prime::Util::PP::todigitstring($n,$base,$len);
}
sub fromdigits {
  my($n,$base) = @_;
  return Math::Prime::Util::PP::fromdigits($n,$base);
}
sub hammingweight {
  my($n) = @_;
  $n = -$n if defined $n && $n < 0;
  _validate_positive_integer($n);
  return Math::Prime::Util::PP::hammingweight($n);
}

sub Pi {
  my($digits) = @_;
  _validate_positive_integer($digits) if defined $digits;
  return Math::Prime::Util::PP::Pi($digits);
}

#############################################################################

sub forprimes (&$;$) {    ## no critic qw(ProhibitSubroutinePrototypes)
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

sub forcomposites(&$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, $beg, $end) = @_;
  if (!defined $end) { $end = $beg; $beg = 4; }
  _validate_num($beg) || _validate_positive_integer($beg);
  _validate_num($end) || _validate_positive_integer($end);
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

sub foroddcomposites(&$;$) { ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, $beg, $end) = @_;
  if (!defined $end) { $end = $beg; $beg = 9; }
  _validate_num($beg) || _validate_positive_integer($beg);
  _validate_num($end) || _validate_positive_integer($end);
  $beg = 9 if $beg < 9;
  $beg++ unless $beg & 1;
  $end = Math::BigInt->new(''.~0) if ref($end) ne 'Math::BigInt' && $end == ~0;
  {
    my $pp;
    local *_ = \$pp;
    for ( ; $beg <= $end ; $beg += 2 ) {
      if (!is_prime($beg)) {
        $pp = $beg;
        $sub->();
      }
    }
  }
}

sub fordivisors (&$) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my($sub, $n) = @_;
  _validate_num($n) || _validate_positive_integer($n);
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
  { my $pp; local *_ = \$pp;
    for my $v (@_) { $pp = $v; return 1 if $sub->(); }
  }
  undef;
}
sub vecall (&@) {       ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  { my $pp; local *_ = \$pp;
    for my $v (@_) { $pp = $v; return if !$sub->(); }
  }
  1;
}
sub vecnone (&@) {      ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  { my $pp; local *_ = \$pp;
    for my $v (@_) { $pp = $v; return if $sub->(); }
  }
  1;
}
sub vecnotall (&@) {    ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  { my $pp; local *_ = \$pp;
    for my $v (@_) { $pp = $v; return 1 if !$sub->(); }
  }
  undef;
}

sub vecfirst (&@) {     ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  #for (@_) { return $_ if &{$sub}(); }  return undef;
  { my $pp; local *_ = \$pp;
    for my $v (@_) { $pp = $v; return $v if $sub->(); }
  }
  undef;
}

sub vecfirstidx (&@) {     ## no critic qw(ProhibitSubroutinePrototypes)
  my $sub = shift;
  { my $pp; local *_ = \$pp; my $i = 0;
    for my $v (@_) { $pp = $v; return $i if $sub->(); $i++; }
  }
  -1;
}

sub vecextract {
  my($aref, $mask) = @_;
  croak "vecextract first argument must be an array reference"
    unless ref($aref) eq 'ARRAY';
  return Math::Prime::Util::PP::vecextract(@_);
}

sub sumdigits {
  return Math::Prime::Util::PP::sumdigits(@_);
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

Copyright 2014-2016 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
