package Math::Prime::Util::PP;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::PP::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PP::VERSION = '0.41';
}

BEGIN {
  do { require Math::BigInt;  Math::BigInt->import(try=>"GMP,Pari"); }
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
  use constant CONST_EULER     => '0.577215664901532860606512090082402431042159335939923598805767';
  use constant CONST_LI2       => '1.04516378011749278484458888919461313652261557815120157583290914407501320521';
  use constant BZERO           => Math::BigInt->bzero;
  use constant BONE            => Math::BigInt->bone;
  use constant BTWO            => Math::BigInt->new(2);
  use constant B_PRIM759       => Math::BigInt->new("64092011671807087969");
  use constant B_PRIM235       => Math::BigInt->new("30");
  use constant PI_TIMES_8      => 25.13274122871834590770114707;
}

{
  my $_have_MPFR = -1;
  sub _MPFR_available {
    if ($_have_MPFR < 0) {
      $_have_MPFR = 0;
      $_have_MPFR = 1 if (!defined $ENV{MPU_NO_MPFR} || $ENV{MPU_NO_MPFR} != 1)
                      && eval { require Math::MPFR; $Math::MPFR::VERSION>=2.03; };
    }
    return $_have_MPFR;
  }
}

my $_precalc_size = 0;
sub prime_precalc {
  my($n) = @_;
  croak "Input must be a positive integer" unless _is_positive_int($n);
  $_precalc_size = $n if $n > $_precalc_size;
}
sub prime_memfree {
  $_precalc_size = 0;
}
sub _get_prime_cache_size { $_precalc_size }
sub _prime_memfreeall { prime_memfree; }


sub _is_positive_int {
  ((defined $_[0]) && $_[0] ne '' && ($_[0] !~ tr/0123456789//c));
}

sub _bigint_to_int {
  return (OLD_PERL_VERSION) ? unpack(UVPACKLET,pack(UVPACKLET,"$_[0]"))
                            : int("$_[0]");
}

sub _upgrade_to_float {
  do { require Math::BigFloat; Math::BigFloat->import(); }
    if !defined $Math::BigFloat::VERSION;
  return Math::BigFloat->new($_[0]);
}

# Get the accuracy of variable x, or the max default from BigInt/BigFloat
# One might think to use ref($x)->accuracy() but numbers get upgraded and
# downgraded willy-nilly, and it will do the wrong thing from the user's
# perspective.
sub _find_big_acc {
  my($x) = @_;

  $b = $x->accuracy();
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


sub _validate_num {
  my($n, $min, $max) = @_;
  croak "Parameter must be defined" if !defined $n;
  return 0 if ref($n);
  croak "Parameter must be a positive integer" if $n eq '';
  croak "Parameter '$n' must be a positive integer"
          if $n =~ tr/0123456789//c && $n !~ /^\+\d+$/;
  croak "Parameter '$n' must be >= $min" if defined $min && $n < $min;
  croak "Parameter '$n' must be <= $max" if defined $max && $n > $max;
  substr($_[0],0,1,'') if substr($n,0,1) eq '+';
  return 0 unless $n < ~0 || int($n) eq ''.~0;
  1;
}

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
      $_[0] = Math::BigInt->new($strn)
    }
  }
  $_[0]->upgrade(undef) if ref($_[0]) && $_[0]->upgrade();
  croak "Parameter '$_[0]' must be >= $min" if defined $min && $_[0] < $min;
  croak "Parameter '$_[0]' must be <= $max" if defined $max && $_[0] > $max;
  1;
}

sub _validate_integer {
  my($n) = @_;
  croak "Parameter must be defined" if !defined $n;
  if (ref($n) eq 'CODE') {
    $_[0] = $_[0]->();
    $n = $_[0];
  }
  my $poscmp = OLD_PERL_VERSION ?  562949953421312 : ''.~0;
  my $negcmp = OLD_PERL_VERSION ? -562949953421312 : -(~0 >> 1);
  if (ref($n) eq 'Math::BigInt') {
    croak "Parameter '$n' must be an integer" if !$n->is_int();
    $_[0] = _bigint_to_int($_[0]) if $n <= $poscmp && $n >= $negcmp;
  } else {
    my $strn = "$n";
    croak "Parameter '$strn' must be an integer"
      if $strn =~ tr/-0123456789//c && $strn !~ /^[-+]?\d+$/;
    if ($n <= $poscmp && $n >= $negcmp) {
      $_[0] = $strn if ref($n);
    } else {
      $_[0] = Math::BigInt->new($strn)
    }
  }
  $_[0]->upgrade(undef) if ref($_[0]) && $_[0]->upgrade();
  1;
}


my @_primes_small = (0,2,3,5);
{
  my $sieveref = _sieve_erat_string(5000);
  my $n = 5;
  foreach my $s (split("0", substr($$sieveref, 3), -1)) {
    $n += 2 + 2 * length($s);
    push @_primes_small, $n if $n <= 5000;
  }
}
my @_prime_next_small = (
   2,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,
   29,29,29,29,29,29,31,31,37,37,37,37,37,37,41,41,41,41,43,43,47,
   47,47,47,53,53,53,53,53,53,59,59,59,59,59,59,61,61,67,67,67,67,67,67,71);

# For wheel-30
my @_prime_indices = (1, 7, 11, 13, 17, 19, 23, 29);
my @_nextwheel30 = (1,7,7,7,7,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,29,29,29,29,29,29,1);
my @_prevwheel30 = (29,29,1,1,1,1,1,1,7,7,7,7,11,11,13,13,13,13,17,17,19,19,19,19,23,23,23,23,23,23);

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

  if (ref($n) eq 'Math::BigInt') {
    return 0 unless Math::BigInt::bgcd($n, B_PRIM759)->is_one;
    return 0 unless _miller_rabin_2($n);
    my $is_esl_prime = is_extra_strong_lucas_pseudoprime($n);
    return ($is_esl_prime)  ?  (($n <= "18446744073709551615") ? 2 : 1)  :  0;
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

  if ($n <= 1_000_000) {
    $n = _bigint_to_int($n) if ref($n) eq 'Math::BigInt';
    my $limit = int(sqrt($n));
    my $i = 61;
    while (($i+30) <= $limit) {
      return 0 if !($n % $i);  $i += 6;
      return 0 if !($n % $i);  $i += 4;
      return 0 if !($n % $i);  $i += 2;
      return 0 if !($n % $i);  $i += 4;
      return 0 if !($n % $i);  $i += 2;
      return 0 if !($n % $i);  $i += 4;
      return 0 if !($n % $i);  $i += 6;
      return 0 if !($n % $i);  $i += 2;
    }
    while (1) {
      last if $i > $limit;  return 0 if !($n % $i);  $i += 6;
      last if $i > $limit;  return 0 if !($n % $i);  $i += 4;
      last if $i > $limit;  return 0 if !($n % $i);  $i += 2;
      last if $i > $limit;  return 0 if !($n % $i);  $i += 4;
      last if $i > $limit;  return 0 if !($n % $i);  $i += 2;
      last if $i > $limit;  return 0 if !($n % $i);  $i += 4;
      last if $i > $limit;  return 0 if !($n % $i);  $i += 6;
      last if $i > $limit;  return 0 if !($n % $i);  $i += 2;
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
  return 0 if int($n) < 0;
  _validate_positive_integer($n);

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
  _validate_positive_integer($n);
  return 0 unless _miller_rabin_2($n);
  if ($n <= 18446744073709551615) {
    return is_almost_extra_strong_lucas_pseudoprime($n) ? 2 : 0;
  }
  return is_extra_strong_lucas_pseudoprime($n) ? 1 : 0;
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
# of naive sieves like found on RosettaCode.

sub _sieve_erat_string {
  my($end) = @_;
  $end-- if ($end & 1) == 0;
  my $s_end = $end >> 1;

  my $whole = int( $s_end / 15);   # Prefill with 3 and 5 already marked.
  croak "Sieve too large" if $whole > 1_145_324_612;  # ~32 GB string
  my $sieve = "100010010010110" . "011010010010110" x $whole;
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
  my($beg,$end) = @_;
  croak "Internal error: segment beg is even" if ($beg % 2) == 0;
  croak "Internal error: segment end is even" if ($end % 2) == 0;
  croak "Internal error: segment end < beg" if $end < $beg;
  croak "Internal error: segment beg should be >= 3" if $beg < 3;
  my $range = int( ($end - $beg) / 2 ) + 1;

  # Prefill with 3 and 5 already marked, and offset to the segment start.
  my $whole = int( ($range+14) / 15);
  my $startp = ($beg % 30) >> 1;
  my $sieve = substr("011010010010110", $startp) . "011010010010110" x $whole;
  # Set 3 and 5 to prime if we're sieving them.
  substr($sieve,0,2) = "00" if $beg == 3;
  substr($sieve,0,1) = "0"  if $beg == 5;
  # Get rid of any extra we added.
  substr($sieve, $range) = '';

  # If the end value is below 7^2, then the pre-sieve is all we needed.
  return \$sieve if $end < 49;

  my $limit = int(sqrt($end)) + 1;
  # For large value of end, it's a huge win to just walk primes.
  my $primesieveref = _sieve_erat($limit);
  my $p = 7-2;
  foreach my $s (split("0", substr($$primesieveref, 3), -1)) {
    $p += 2 + 2 * length($s);
    my $p2 = $p*$p;
    if ($p2 < $beg) {
      my $f = 1+int(($beg-1)/$p);
      $p2 = $p * ($f + !($f & 1));
    } elsif ($p2 > $end) { last; }
    # With large bases and small segments, it's common to find we don't hit
    # the segment at all.  Skip all the setup if we find this now.
    if ($p2 <= $end) {
      # Inner loop marking multiples of p
      # (everything is divided by 2 to keep inner loop simpler)
      my $filter_end = ($end - $beg) >> 1;
      my $filter_p2  = ($p2  - $beg) >> 1;
      while ($filter_p2 <= $filter_end) {
        substr($sieve, $filter_p2, 1) = "1";
        $filter_p2 += $p;
      }
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
  _validate_positive_integer($low);
  _validate_positive_integer($high);

  return if $low > $high;

  my @primes;
  $low-- if $low >= 2;
  my $curprime = next_prime($low);
  while ($curprime <= $high) {
    push @primes, $curprime;
    $curprime = next_prime($curprime);
  }
  return \@primes;
}

sub primes {
  my($low,$high) = @_;
  if (scalar @_ > 1) {
    _validate_positive_integer($low);
    _validate_positive_integer($high);
  } else {
    ($low,$high) = (2, $low);
    _validate_positive_integer($high);
  }
  my $sref = [];
  return $sref if ($low > $high) || ($high < 2);

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

  if ($low == 7) {
    my $sieveref = _sieve_erat($high);
    my $n = $low - 2;
    foreach my $s (split("0", substr($$sieveref, 3), -1)) {
      $n += 2 + 2 * length($s);
      push @$sref, $n if $n <= $high;
    }
  } else {
    my $sieveref = _sieve_segment($low,$high);
    my $n = $low - 2;
    foreach my $s (split("0", $$sieveref, -1)) {
      $n += 2 + 2 * length($s);
      push @$sref, $n if $n <= $high;
    }
  }
  $sref;
}

sub next_prime {
  my($n) = @_;
  _validate_positive_integer($n);

  return $_prime_next_small[$n] if $n <= $#_prime_next_small;
  if ($n < $_primes_small[-1]) {
    my $i = _tiny_prime_count($n);  # Binary search
    return $_primes_small[$i+1];
  }

  return Math::BigInt->new(MPU_32BIT ? "4294967311" : "18446744073709551629")
    if ref($n) ne 'Math::BigInt' && $n >= MPU_MAXPRIME;
  # n is now either 1) not bigint and < maxprime, or (2) bigint and >= uvmax

  #$n = ($n+1) | 1;
  #while (    !($n%3) || !($n%5) || !($n%7) || !($n%11) || !($n%13)
  #        || !_is_prime7($n) ) {
  #  $n += 2;
  #}
  my $m = $n % 30;
  my $d = ($n - $m) / 30;
  if ($m == 29) { $d++;  $m = 1;} else { $m = $_nextwheel30[$m]; }
  $n = $d*30+$m;
  while ( !($n%7) || !_is_prime7($n) ) {
    $m = $_nextwheel30[$m];
    $d++ if $m == 1;
    $n = $d*30+$m;
  }
  return $n;
}

sub prev_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  if ($n <= 11) {
    return ($n <= 2) ? 0 : ($n <= 3) ? 2 : ($n <= 5) ? 3 : ($n <= 7) ? 5 : 7;
  }

  #$n++ if ($n % 2) == 0;
  #do {
  #  $n -= 2;
  #} while ( (($n % 3) == 0) || (($n % 5) == 0) || (!_is_prime7($n)) );
  #return $n;

  $n -= ($n & 1) ? 2 : 1;
  my $nmod6 = $n % 6;
  if ($nmod6 == 5) {
    return $n if ($n % 5) != 0 && ($n % 7) != 0 && _is_prime7($n);
    $n -= 4;
  } elsif ($nmod6 == 3) {
    $n -= 2;
  }

  while (1) {
    return $n if ($n % 5) != 0 && ($n % 7) != 0 && _is_prime7($n);
    $n -= 2;
    return $n if ($n % 5) != 0 && ($n % 7) != 0 && _is_prime7($n);
    $n -= 4;
  }
  return $n;

  # This is faster for larger intervals, slower for short ones.
  #my $base = 30 * int($n/30);
  #my $in = 0;  $in++ while ($n - $base) > $_prime_indices[$in];
  #if (--$in < 0) {  $base -= 30; $in = 7;  }
  #$n = $base + $_prime_indices[$in];
  #while (!_is_prime7($n)) {
  #  if (--$in < 0) {  $base -= 30; $in = 7;  }
  #  $n = $base + $_prime_indices[$in];
  #}
  #$n;

  #my $m = $n % 30;
  #my $d = int( ($n - $m) / 30 );
  #do {
  #  $m = $_prevwheel30[$m];
  #  $d-- if $m == 29;
  #} while (!_is_prime7($d*30+$m));
  #$d*30+$m;
}

sub partitions {
  my $n = shift;

  my $d = int(sqrt($n+1));
  my @pent = (1, map { (($_*(3*$_+1))>>1, (($_+1)*(3*$_+2))>>1) } 1 .. $d);
  my $ZERO = ($n >= ((~0 > 4294967295) ? 400 : 270)) ? BZERO : 0;
  my @part = ($ZERO+1);
  foreach my $j (scalar @part .. $n) {
    my ($psum1, $psum2, $k) = ($ZERO, $ZERO, 1);
    foreach my $p (@pent) {
      last if $p > $j;
      if ((++$k) & 2) { $psum1 += $part[ $j - $p ] }
      else            { $psum2 += $part[ $j - $p ] }
    }
    $part[$j] = $psum1 - $psum2;
  }
  return $part[$n];
}

sub primorial {
  my $n = shift;

  my $max = (MPU_32BIT) ? 29 : (OLD_PERL_VERSION) ? 43 : 53;
  my $pn = (ref($_[0]) eq 'Math::BigInt') ? $_[0]->copy->bone()
         : ($n >= $max) ? Math::BigInt->bone()
         : 1;
  if (ref($pn) eq 'Math::BigInt') {
    my $start = 2;
    if ($n >= 97) {
      $start = 101;
      $pn->bdec->badd(Math::BigInt->new("2305567963945518424753102147331756070"));
    }
    my @plist = @{primes($start,$n)};
    while (@plist > 2 && $plist[2] < 1625) {
      $pn->bmul( Math::BigInt->new(shift(@plist)*shift(@plist)*shift(@plist)) );
    }
    while (@plist > 1 && $plist[1] < 65536) {
      $pn->bmul( Math::BigInt->new(shift(@plist)*shift(@plist)) );
    }
    $pn->bmul($_) for @plist;
  } else {
    foreach my $p (@{primes($n)}) {  $pn *= $p;  }
  }
  return $pn;
}

sub consecutive_integer_lcm {
  my $n = shift;

  my $max = (MPU_32BIT) ? 22 : (OLD_PERL_VERSION) ? 37 : 46;
  my $pn = ref($n) ? ref($n)->new(1) : ($n >= $max) ? Math::BigInt->bone() : 1;
  for (my $p = 2; $p <= $n; $p = next_prime($p)) {
    my($p_power, $pmin) = ($p, int($n/$p));
    $p_power *= $p while $p_power <= $pmin;
    $pn *= $p_power;
  }
  $pn = _bigint_to_int($pn) if $pn <= ''.~0;
  return $pn;
}

sub jordan_totient {
  my($k, $n) = @_;
  return ($n == 1) ? 1 : 0  if $k == 0;
  return euler_phi($n)      if $k == 1;
  return 0 if $n < 0;
  return ($n == 1) ? 1 : 0  if $n <= 1;

  my @pe = Math::Prime::Util::factor_exp($n);
  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';
  my $totient = BONE->copy;
  foreach my $f (@pe) {
    my ($p, $e) = @$f;
    $p = Math::BigInt->new("$p")->bpow($k);
    $totient->bmul($p->copy->bdec());
    $totient->bmul($p) for 2 .. $e;
  }
  $totient = _bigint_to_int($totient) if $totient->bacmp(''.~0) <= 0;
  return $totient;
}

sub euler_phi {
  return euler_phi_range(@_) if scalar @_ > 1;
  my($n) = @_;
  return 0 if $n < 0;
  return $n if $n <= 1;

  my @pe = Math::Prime::Util::factor_exp($n);
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

sub euler_phi_range {
  my($n, $nend) = @_;
  return () if $nend < $n;
  return euler_phi($n) if $n == $nend;
  my @totients;
  if ($nend > 2**30) {
    while ($n < $nend) {
      push @totients, euler_phi($n++);
    }
  } else {
    @totients = (0 .. $nend);
    foreach my $i (2 .. $nend) {
      next unless $totients[$i] == $i;
      $totients[$i] = $i-1;
      foreach my $j (2 .. int($nend / $i)) {
        $totients[$i*$j] -= $totients[$i*$j]/$i;
      }
    }
    splice(@totients, 0, $n) if $n > 0;
  }
  return @totients;
}

sub moebius {
  return moebius_range(@_) if scalar @_ > 1;
  my($n) = @_;
  return ($n == 1) ? 1 : 0  if $n <= 1;
  return 0 if ($n >= 49) && (!($n % 4) || !($n % 9) || !($n % 25) || !($n%49) );
  my @factors = Math::Prime::Util::factor($n);
  foreach my $i (1 .. $#factors) {
    return 0 if $factors[$i] == $factors[$i-1];
  }
  return ((scalar @factors) % 2) ? -1 : 1;
}

sub moebius_range {
  my($lo, $hi) = @_;
  return () if $hi < $lo;
  return moebius($lo) if $lo == $hi;
  if ($hi > 2**32) {
    my @mu;
    while ($lo < $hi) {
      push @mu, moebius($lo++);
    }
    return @mu;
  }
  my @mu = map { 1 } $lo .. $hi;
  $mu[0] = 0 if $lo == 0;
  my($p, $sqrtn) = (2, int(sqrt($hi)+0.5));
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
    $p = next_prime($p);
  }
  foreach my $i ($lo .. $hi) {
    my $m = $mu[$i-$lo];
    $m *= -1 if abs($m) != $i;
    $mu[$i-$lo] = ($m>0) - ($m<0);
  }
  return @mu;
}

sub mertens {
  my($n) = @_;
    # This is the most basic Deléglise and Rivat algorithm.  u = n^1/2
  # and no segmenting is done.  Their algorithm uses u = n^1/3, breaks
  # the summation into two parts, and calculates those in segments.  Their
  # computation time growth is half of this code.
  return $n if $n <= 1;
  my $u = int(sqrt($n));
  my @mu = (0, Math::Prime::Util::moebius(1, $u)); # Hold values of mu for 0-u
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

sub liouville {
  my($n) = @_;
  my $l = (-1) ** scalar factor($n);
  return $l;
}

# Exponential of Mangoldt function (A014963).
# Return p if n = p^m [p prime, m >= 1], 1 otherwise.
sub exp_mangoldt {
  my($n) = @_;
  return 1 if defined $n && $n <= 1;  # n <= 1
  return 2 if ($n & ($n-1)) == 0;     # n power of 2
  return 1 unless $n & 1;             # even n can't be p^m
  return $n if Math::Prime::Util::is_prob_prime($n); # prime n returns n

  my $k = Math::Prime::Util::is_power($n);
  if ($k >= 2) {
    my $root = Math::BigInt->new("$n")->broot($k);
    return $root if Math::Prime::Util::is_prob_prime($root);
  }
  1;
}

sub carmichael_lambda {
  my($n) = @_;
  return euler_phi($n) if $n < 8;                # = phi(n) for n < 8
  return euler_phi($n)/2 if ($n & ($n-1)) == 0;  # = phi(n)/2 for 2^k, k>2

  my @pe = Math::Prime::Util::factor_exp($n);
  $pe[0]->[1]-- if $pe[0]->[0] == 2 && $pe[0]->[1] > 2;

  my $lcm = Math::BigInt::blcm(
    map { $_->[0]->copy->bpow($_->[1]->copy->bdec)->bmul($_->[0]->copy->bdec) }
    map { [ map { Math::BigInt->new("$_") } @$_ ] }
    @pe
  );
  $lcm = _bigint_to_int($lcm) if $lcm->bacmp(''.~0) <= 0;
  return $lcm;
}


my @_ds_overflow =  # We'll use BigInt math if the input is larger than this.
  (~0 > 4294967295)
   ? (124, 3000000000000000000, 3000000000, 2487240, 64260, 7026)
   : ( 50,           845404560,      52560,    1548,   252,   84);
sub divisor_sum {
  my($n, $k) = @_;
  return 1 if $n == 1;

  if (defined $k && ref($k) eq 'CODE') {
    my $sum = $n-$n;
    my $refn = ref($n);
    foreach my $d (Math::Prime::Util::divisors($n)) {
      $sum += $k->( $refn ? $refn->new("$d") : $d );
    }
    return $sum;
  }

  croak "Second argument must be a code ref or number"
    unless !defined $k || _validate_num($k) || _validate_positive_integer($k);
  $k = 1 if !defined $k;

  my $will_overflow = ($k == 0) ? (length($n) >= $_ds_overflow[0])
                    : ($k <= 5) ? ($n >= $_ds_overflow[$k])
                    : 1;

  # The standard way is:
  #    my $pk = $f ** $k;  $product *= ($pk ** ($e+1) - 1) / ($pk - 1);
  # But we get less overflow using:
  #    my $pk = $f ** $k;  $product *= $pk**E for E in 0 .. e
  # Also separate BigInt and do fiddly bits for better performance.

  my $product = 1;
  if (!$will_overflow) {
    foreach my $f (Math::Prime::Util::factor_exp($n)) {
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
    foreach my $f (Math::Prime::Util::factor_exp($n)) {
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

#############################################################################
#                       Lehmer prime count
#
#my @_s0 = (0);
#my @_s1 = (0,1);
#my @_s2 = (0,1,1,1,1,2);
my @_s3 = (0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7,8);
my @_s4 = (0,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7,7,8,8,8,8,8,8,9,9,9,9,10,10,11,11,11,11,12,12,12,12,12,12,13,13,13,13,13,13,14,14,15,15,15,15,15,15,16,16,16,16,17,17,18,18,18,18,18,18,19,19,19,19,20,20,20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,22,23,23,24,24,24,24,25,25,26,26,26,26,27,27,27,27,27,27,27,27,28,28,28,28,28,28,29,29,29,29,30,30,30,30,30,30,31,31,32,32,32,32,33,33,33,33,33,33,34,34,35,35,35,35,35,35,36,36,36,36,36,36,37,37,37,37,38,38,39,39,39,39,40,40,40,40,40,40,41,41,42,42,42,42,42,42,43,43,43,43,44,44,45,45,45,45,46,46,47,47,47,47,47,47,47,47,47,47,48);
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
  return _tablephi($x,$a) if $a <= 6;
  $primes = primes(Math::Prime::Util::nth_prime_upper($a+1)) unless defined $primes;
  return ($x > 0 ? 1 : 0) if $x < $primes->[$a];

  my $sum = 0;
  my %vals = ( $x => 1 );
  while ($a > 6) {
    my $primea = $primes->[$a-1];
    my %newvals;
    while (my($v,$c) = each %vals) {
      my $sval = int($v / $primea);
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
  $high-- unless ($high & 1);
  return 1 + ${_sieve_erat($high)} =~ tr/0//;
}

sub _count_with_sieve {
  my ($sref, $low, $high) = @_;
  ($low, $high) = (2, $low) if !defined $high;
  my $count = 0;
  if   ($low < 3) { $low = 3; $count++; }
  else            { $low |= 1; }
  $high-- unless ($high & 1);
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
  my $z = (ref($x) ne 'Math::BigInt')
        ? int(sqrt($x+0.5))
        : int(Math::BigFloat->new($x)->badd(0.5)->bsqrt->bfloor->bstr);
  my $a = _lehmer_pi(int(sqrt($z)+0.5));
  my $b = _lehmer_pi($z);
  my $c = _lehmer_pi(int( (ref($x) ne 'Math::BigInt')
                          ? $x**(1/3)+0.5
                          : Math::BigFloat->new($x)->broot(3)->badd(0.5)->bfloor
                     ));
  ($z, $a, $b, $c) = map { (ref($_) =~ /^Math::Big/) ? _bigint_to_int($_) : $_ }
                     ($z, $a, $b, $c);

  # Generate at least b primes.
  my $bth_prime_upper = ($b <= 10) ? 29 : int($b*(log($b) + log(log($b)))) + 1;
  my $primes = primes( $bth_prime_upper );

  my $sum = int(($b + $a - 2) * ($b - $a + 1) / 2);
  $sum += legendre_phi($x, $a, $primes);

  # Get a big sieve for our primecounts.  The C code compromises with either
  # b*10 or x^3/5, as that cuts out all the inner loop sieves and about half
  # of the big outer loop counts.
  # Our sieve count isn't nearly as optimized here, so error on the side of
  # more primes.  This uses a lot more memory but saves a lot of time.
  my $sref = _sieve_erat( int($x / $primes->[$a] / 5) );

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
  if (!defined $high) {
    $high = $low;
    $low = 2;
  }
  _validate_positive_integer($low);
  _validate_positive_integer($high);

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
    my $curprime = next_prime($low-1);
    while ($curprime <= $high) {
      $count++;
      $curprime = next_prime($curprime);
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
  _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  if ($n > MPU_MAXPRIMEIDX && ref($n) ne 'Math::BigFloat') {
    do { require Math::BigFloat; Math::BigFloat->import(); }
      if !defined $Math::BigFloat::VERSION;
    $n = Math::BigFloat->new("$n")
  }

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
  _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if $n > MPU_MAXPRIMEIDX || $n > 2**45;

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

# The nth prime will be greater than or equal to this number
sub nth_prime_lower {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n) if $n > MPU_MAXPRIMEIDX || $n > 2**45;

  my $flogn  = log($n);
  my $flog2n = log($flogn);  # Note distinction between log_2(n) and log^2(n)

  # Dusart 1999 page 14, for all n >= 2
  #my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.25)/$flogn));
  # Dusart 2010 page 2, for all n >= 3
  my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.10)/$flogn));

  return int($lower);
}

sub nth_prime_approx {
  my($n) = @_;
  _validate_num($n) || _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  $n = _upgrade_to_float($n)
    if ref($n) eq 'Math::BigInt' || $n >= MPU_MAXPRIMEIDX;

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
  elsif ($n <        775) { $approx +=  7.52* $order; }
  elsif ($n <       1271) { $approx +=  5.6 * $order; }
  elsif ($n <       2000) { $approx +=  5.2 * $order; }
  elsif ($n <       4000) { $approx +=  4.3 * $order; }
  elsif ($n <      12000) { $approx +=  3.0 * $order; }
  elsif ($n <     150000) { $approx +=  2.1 * $order; }
  elsif ($n <  200000000) { $approx +=  0.0 * $order; }
  else                    { $approx += -0.025 * $order; }
  # If we want the asymptotic approximation to be >= actual, use -0.010.

  return int($approx + 0.5);
}

#############################################################################

sub prime_count_approx {
  my($x) = @_;
  _validate_num($x) || _validate_positive_integer($x);

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
  if ( $x < 1e36 || _MPFR_available() ) {
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
  if ($x > 1000_000_000_000 && Math::Prime::Util::prime_get_config()->{'assume_rh'}) {
    # Schoenfeld bound
    my $lix = LogarithmicIntegral($x);
    my $sqx = sqrt($x);
    if (ref($x) eq 'Math::BigFloat') {
      my $xdigits = _find_big_acc($x);
      $result = $lix - (($sqx*$flogx) / (Math::BigFloat->bpi($xdigits)*8));
    } else {
      $result = $lix - (($sqx*$flogx) / PI_TIMES_8);
    }
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
  if ($x > 10000_000_000_000 && Math::Prime::Util::prime_get_config()->{'assume_rh'}) {
    # Schoenfeld bound
    my $lix = LogarithmicIntegral($x);
    my $sqx = sqrt($x);
    if (ref($x) eq 'Math::BigFloat') {
      my $xdigits = _find_big_acc($x);
      $result = $lix + (($sqx*$flogx) / (Math::BigFloat->bpi($xdigits)*8));
    } else {
      $result = $lix + (($sqx*$flogx) / PI_TIMES_8);
    }
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

sub twin_prime_count {
  my($low,$high) = @_;
  if (defined $high) { _validate_positive_integer($low); }
  else               { ($low,$high) = (2, $low);         }
  _validate_positive_integer($high);
  my $sum = 0;
  # TODO: I suspect calling primes() on segments would be faster in most cases.
  if ($high >= $low) {
    my $p = next_prime($low-1);
    my $p2 = next_prime($p);
    while ($p <= $high) {
      $sum++ if ($p2-$p) == 2;
      ($p, $p2) = ($p2, next_prime($p2));
    }
  }
  $sum;
}

sub twin_prime_count_approx {
  my($n) = @_;
  return twin_prime_count(3,$n) if $n < 2000;
  $n = _upgrade_to_float($n) if ref($n);
  my $logn = log($n);
  my $li2 = ExponentialIntegral($logn) + 2.8853900817779268147198494 - ($n/$logn);
  if    ($n <     4000) { $li2 *= 1.0005 * log(log(log($n*4000))); }
  elsif ($n <     8000) { $li2 *= 0.9734 * log(log(log($n*4000))); }
  elsif ($n <    32000) { $li2 *= 0.8967 * log(log(log($n*4000))); }
  elsif ($n <   200000) { $li2 *= 0.8937 * log(log(log($n*4000))); }
  elsif ($n <  1000000) { $li2 *= 0.8793 * log(log(log($n*4000))); }
  elsif ($n <  4000000) { $li2 *= 0.8766 * log(log(log($n*4000))); }
  elsif ($n < 10000000) { $li2 *= 0.8664 * log(log(log($n*4000))); }
  return int(1.32032363169373914785562422 * $li2 + 0.5);
}

sub nth_twin_prime {
  my($n) = @_;
  my($nth, $p, $p2) = (0, 0, 3);
  while ($n > 0) {
    if (($p2-$p) == 2) { $nth = $p; $n--; }
    ($p, $p2) = ($p2, next_prime($p2));
  }
  $nth;
}

sub nth_twin_prime_approx {
  my($n) = @_;
  return nth_twin_prime($n) if $n < 6;
  $n = _upgrade_to_float($n) if ref($n);
  my $nlogn2 = $n * log($n) * log($n);
  return int(5.023 * $nlogn2/log(log(1600*$n*$n))) if $n > 59 && $n < 1200;

  my $lo = int(1.0 * $nlogn2);
  my $hi = int(3.0 * $nlogn2 + 3);
  while ($lo < $hi) {
    my $mid = $lo + (($hi-$lo) >> 1);
    if (twin_prime_count_approx($mid) < $n) { $lo = $mid+1; }
    else                                    { $hi = $mid;   }
  }
  return $lo;
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

# Make sure to work around RT71548 here, and correct lcm semantics.
sub gcd {
  my $gcd = Math::BigInt::bgcd( map {
    my $v = ($_ < 2147483647 || ref($_) eq 'Math::BigInt') ? $_ : "$_";
    $v;
  } @_ );
  $gcd = _bigint_to_int($gcd) if $gcd->bacmp(''.~0) <= 0;
  return $gcd;
}
sub lcm {
  my $lcm = Math::BigInt::blcm( map {
    my $v = ($_ < 2147483647 || ref($_) eq 'Math::BigInt') ? $_ : "$_";
    return 0 if $v == 0;
    $v = -$v if $v < 0;
    $v;
  } @_ );
  $lcm = _bigint_to_int($lcm) if $lcm->bacmp(''.~0) <= 0;
  return $lcm;
}
sub vecsum {
  my $sum = 0;
  my $neglim = -(~0 >> 1) - 1;
  foreach my $v (@_) {
    $sum += $v;
    if ($sum > ~0 || $sum < $neglim) {
      $sum = BZERO->copy;
      $sum->badd("$_") for @_;
      return $sum;
    }
  }
  $sum;
}
sub invmod {
  my($a,$n) = @_;
  return if $n == 0 || $a == 0;
  return 0 if $n == 1;
  $n = -$n if $n < 0;  # Pari semantics
  if ($n > ~0) {
    my $invmod = Math::BigInt->new("$a")->bmodinv("$n");
    return if !defined $invmod || $invmod->is_nan;
    $invmod = _bigint_to_int($invmod) if $invmod->bacmp(''.~0) <= 0;
    return $invmod;
  }
  my($t,$nt,$r,$nr) = (0, 1, $n, $a % $n);
  while ($nr != 0) {
    # Use mod before divide to force correct behavior with high bit set
    my $quot = int( ($r-($r % $nr))/$nr );
    ($nt,$t) = ($t-$quot*$nt,$nt);
    ($nr,$r) = ($r-$quot*$nr,$nr);
  }
  return if $r > 1;
  $t += $n if $t < 0;
  $t;
}

# no validation, x is allowed to be negative, y must be >= 0
sub _gcd_ui {
  my($x, $y) = @_;
  if ($y < $x) { ($x, $y) = ($y, $x); }
  elsif ($x < 0) { $x = -$x; }
  while ($y > 0) {
    # y1 <- x0 % y0 ; x1 <- y0
    my $t = $y;
    $y = $x % $y;
    $x = $t;
  }
  $x;
}

sub is_power {
  my ($n, $a) = @_;
  return 0 if $n <= 3;
  if (defined $a && $a != 0) {
    return _is_perfect_square($n) if $a == 2;
    $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';
    return $n->copy->broot($a)->bint->bpow($a) == $n;
  }
  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';
  my $e = 2;
  while (1) {
    my $root = $n->copy()->broot($e);
    last if $root->is_one();
    if ($root->copy->bpow($e) == $n) {
      my $next = is_power($root);
      return ($next == 0) ? $e : $e * $next;
    }
    $e = next_prime($e);
  }
  0;
}

sub valuation {
  my($n, $k) = @_;
  return 0 if $n < 2 || $k < 2;
  my $v = 0;
  if ($k == 2) { # Accelerate power of 2
    if (ref($n) eq 'Math::BigInt') {   # This can pay off for big inputs
      return 0 unless $n->is_even;
      my $s = $n->as_bin;              # We could do same for k=10
      return length($s) - rindex($s,'1') - 1;
    }
    while (!($n & 0xFFFF) ) {  $n >>=16;  $v +=16;  }
    while (!($n & 0x000F) ) {  $n >>= 4;  $v += 4;  }
  }
  while ( !($n % $k) ) {
    $n /= $k;
    $v++;
  }
  $v;
}


sub is_pseudoprime {
  my($n, $base) = @_;
  return 0 if int($n) < 0;
  _validate_positive_integer($n);

  if ($n < 5) { return ($n == 2) || ($n == 3) ? 1 : 0; }
  $base = 2 if !defined $base;
  croak "Base $base is invalid" if $base < 2;
  if ($base >= $n) {
    $base = $base % $n;
    return 1 if $base <= 1 || $base == $n-1;
  }
  my $x = (ref($n) eq 'Math::BigInt')
        ? $n->copy->bzero->badd($base)->bmodpow($n-1,$n)
        : _powmod($base, $n-1, $n);
  return ($x == 1) ? 1 : 0;
}

sub _miller_rabin_2 {
  my($n, $nm1, $s, $d) = @_;

  if ( ref($n) eq 'Math::BigInt' ) {

    if (!defined $nm1) {
      $nm1 = $n->copy->bdec();
      $s = 0;
      $d = $nm1->copy;
      do {
        $s++;
        $d->brsft(BONE);
      } while $d->is_even;
    }
    my $x = BTWO->copy->bmodpow($d,$n);
    return 1 if $x->is_one || $x->bcmp($nm1) == 0;
    foreach my $r (1 .. $s-1) {
      $x->bmul($x)->bmod($n);
      last if $x->is_one;
      return 1 if $x->bcmp($nm1) == 0;
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
  return 0 if int($n) < 0;
  _validate_positive_integer($n);

  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0;

  if ($bases[0] == 2) {
    return 0 unless _miller_rabin_2($n);
    shift @bases;
    return 1 unless @bases;
  }

  # Die on invalid bases
  foreach my $base (@bases) { croak "Base $base is invalid" if $base < 2 }
  # Make sure we handle big bases ok.
  @bases = grep { $_ > 1 }  map { ($_ >= $n) ? $_ % $n : $_ }  @bases;

  if ( ref($n) eq 'Math::BigInt' ) {

    my $nminus1 = $n->copy->bdec();
    my $s = 0;
    my $d = $nminus1->copy;
    do {  # n is > 3 and odd, so n-1 must be even
      $s++;
      $d->brsft(BONE);
    } while $d->is_even;
    # Different way of doing the above.  Fewer function calls, slower on ave.
    #my $dbin = $nminus1->as_bin;
    #my $last1 = rindex($dbin, '1');
    #my $s = length($dbin)-2-$last1+1;
    #my $d = $nminus1->copy->brsft($s);

    foreach my $ma (@bases) {
      my $x = $n->copy->bzero->badd($ma)->bmodpow($d,$n);
      next if $x->is_one || $x->bcmp($nminus1) == 0;
      foreach my $r (1 .. $s-1) {
        $x->bmul($x); $x->bmod($n);
        return 0 if $x->is_one;
        do { $ma = 0; last; } if $x->bcmp($nminus1) == 0;
      }
      return 0 if $ma != 0;
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
    if (   ref($a) eq 'Math::BigInt' && $a <= ''.~0
        && ref($b) eq 'Math::BigInt' && $b <= ''.~0) {
      return $k * kronecker(_bigint_to_int($a),_bigint_to_int($b));
    }
  }
  return ($b == 1) ? $k : 0;
}

sub _binomialu {
  my($r, $n, $k) = (1, @_);
  return ($k == $n) ? 1 : 0 if $k >= $n;
  $k = $n - $k if $k > ($n >> 1);
  foreach my $d (1 .. $k) {
    if ($r >= int(~0/$n)) {
      my($g, $nr, $dr);
      $g = _gcd_ui($n, $d);   $nr = int($n/$g);   $dr = int($d/$g);
      $g = _gcd_ui($r, $dr);  $r  = int($r/$g);   $dr = int($dr/$g);
      return 0 if $r >= int(~0/$nr);
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

  # 1. Try GMP
  if (defined &Math::Prime::Util::GMP::binomial && Math::Prime::Util::prime_get_config()->{'gmp'}) {
    return Math::Prime::Util::_reftyped($_[0], Math::Prime::Util::GMP::binomial($n,$k));
  }

  # 2. Exit early for known 0 cases, and adjust k to be positive.
  if ($n >= 0) {  return 0 if $k < 0 || $k > $n;  }
  else         {  return 0 if $k < 0 && $k > $n;  }
  $k = $n - $k if $k < 0;

  # 3. Try to do in integer Perl
  my $r;
  if ($n >= 0) {
    $r = _binomialu($n, $k);
    return $r  if $r > 0;
  } else {
    $r = _binomialu(-$n+$k-1, $k);
    return $r   if $r > 0 && !($k & 1);
    return -$r  if $r > 0 && $r <= (~0>>1);
  }

  # 4. Overflow.  Solve using Math::BigInt
  return 1 if $k == 0;        # Work around bug in old
  return $n if $k == $n-1;    # Math::BigInt (fixed in 1.90)
  if ($n >= 0) {
    $r = Math::BigInt->new($n)->bnok($k);
    $r = _bigint_to_int($r) if $r->bacmp(''.~0) <= 0;
  } else { # Math::BigInt is incorrect for negative n
    $r = Math::BigInt->new(-$n+$k-1)->bnok($k);
    if ($k & 1) {
      $r->bneg;
      $r = _bigint_to_int($r) if $r->bacmp(''.(~0>>1)) <= 0;
    } else {
      $r = _bigint_to_int($r) if $r->bacmp(''.~0) <= 0;
    }
  }
  $r;
}

sub _is_perfect_square {
  my($n) = @_;

  if (ref($n) eq 'Math::BigInt') {
    my $mc = _bigint_to_int($n & 31);
    if ($mc==0||$mc==1||$mc==4||$mc==9||$mc==16||$mc==17||$mc==25) {
      my $sq = $n->copy->bsqrt->bfloor;
      $sq->bmul($sq);
      return 1 if $sq == $n;
    }
  } else {
    my $mc = $n & 31;
    if ($mc==0||$mc==1||$mc==4||$mc==9||$mc==16||$mc==17||$mc==25) {
      my $sq = int(sqrt($n));
      return 1 if ($sq*$sq) == $n;
    }
  }
  0;
}

sub znorder {
  my($a, $n) = @_;
  return if $n <= 0;
  return (undef,1)[$a] if $a <= 1;
  return 1 if $n == 1;

  # Sadly, Calc/FastCalc are horrendously slow for this function.
  return if Math::BigInt::bgcd($a, $n) > 1;

  # The answer is one of the divisors of phi(n) and lambda(n).
  my $lambda = Math::Prime::Util::carmichael_lambda($n);
  $a = Math::BigInt->new("$a") unless ref($a) eq 'Math::BigInt';

  # This is easy and usually fast, but can bog down with too many divisors.
  if ($lambda <= 2**64) {
    foreach my $k (Math::Prime::Util::divisors($lambda)) {
      return $k if $a->copy->bmodpow("$k", $n)->is_one;
    }
    return;
  }

  # Algorithm 1.7 from A. Das applied to Carmichael Lambda.
  $lambda = Math::BigInt->new("$lambda") unless ref($lambda) eq 'Math::BigInt';
  my $k = Math::BigInt->bone;
  foreach my $f (Math::Prime::Util::factor_exp($lambda)) {
    my($pi, $ei, $enum) = (Math::BigInt->new("$f->[0]"), $f->[1], 0);
    my $phidiv = $lambda / ($pi**$ei);
    my $b = $a->copy->bmodpow($phidiv, $n);
    while ($b != 1) {
      return if $enum++ >= $ei;
      $b = $b->copy->bmodpow($pi, $n);
      $k *= $pi;
    }
  }
  $k = _bigint_to_int($k) if $k->bacmp(''.~0) <= 0;
  return $k;
}

# This is just a stupid brute force search.
sub znlog {
  my ($a,$g,$p) =
    map { ref($_) eq 'Math::BigInt' ? $_ : Math::BigInt->new("$_") } @_;
  for (my $k = BONE->copy; $k < $p; $k->binc) {
    my $t = $g->copy->bmodpow($k, $p);
    if ($t == $a) {
      $k = _bigint_to_int($k) if $k->bacmp(''.~0) <= 0;
      return $k;
    }
  }
  return;
}

sub znprimroot {
  my($n) = @_;
  $n = -$n if $n < 0;
  if ($n <= 4) {
    return if $n == 0;
    return $n-1;
  }
  return if $n % 4 == 0;
  my $a = 1;
  my $phi = $n-1;
  if (!is_prob_prime($n)) {
    $phi = euler_phi($n);
    # Check that a primitive root exists.
    return if $phi != Math::Prime::Util::carmichael_lambda($n);
  }
  my @exp = map { Math::BigInt->new("$_") }
            map { int($phi/$_->[0]) }
            Math::Prime::Util::factor_exp($phi);
  #print "phi: $phi  factors: ", join(",",factor($phi)), "\n";
  #print "  exponents: ", join(",", @exp), "\n";
  my $bign = (ref($n) eq 'Math::BigInt') ? $n : Math::BigInt->new("$n");
  while (1) {
    my $fail = 0;
    do { $a++ } while kronecker($a,$n) == 0;
    return if $a >= $n;
    foreach my $f (@exp) {
      if ( Math::BigInt->new($a)->bmodpow($f, $bign)->is_one ) {
        $fail = 1;
        last;
      }
    }
    return $a if !$fail;
  }
}


# Find first D in sequence (5,-7,9,-11,13,-15,...) where (D|N) == -1
sub _lucas_selfridge_params {
  my($n) = @_;

  # D is typically quite small: 67 max for N < 10^19.  However, it is
  # theoretically possible D could grow unreasonably.  I'm giving up at 4000M.
  my $d = 5;
  my $sign = 1;
  while (1) {
    my $gcd = (ref($n) eq 'Math::BigInt') ? Math::BigInt::bgcd($d, $n)
                                          : _gcd_ui($d, $n);
    return (0,0,0) if $gcd > 1 && $gcd != $n;  # Found divisor $d
    my $j = kronecker($d * $sign, $n);
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
    my $gcd = (ref($n) eq 'Math::BigInt') ? Math::BigInt::bgcd($D, $n)
                                          : _gcd_ui($D, $n);
    return (0,0,0) if $gcd > 1 && $gcd != $n;  # Found divisor $d
    last if kronecker($D, $n) == -1;
    $P += $increment;
    croak "Could not find Jacobi sequence for $n" if $P > 65535;
    $D = $P*$P - 4;
  }
  ($P, $Q, $D);
}

# returns U_k, V_k, Q_k all mod n
sub lucas_sequence {
  my($n, $P, $Q, $k) = @_;

  croak "lucas_sequence: n must be >= 2" if $n < 2;
  croak "lucas_sequence: k must be >= 0" if $k < 0;
  croak "lucas_sequence: P out of range" if abs($P) >= $n;
  croak "lucas_sequence: Q out of range" if abs($Q) >= $n;

  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';

  my $ZERO = $n->copy->bzero;
  $P = $ZERO+$P unless ref($P) eq 'Math::BigInt';
  $Q = $ZERO+$Q unless ref($Q) eq 'Math::BigInt';
  my $D = $P*$P - BTWO*BTWO*$Q;
  if ($D->is_zero) {
    my $S = ($ZERO+$P) >> 1;
    my $U = $S->copy->bmodpow($k-1,$n)->bmul($k)->bmod($n);
    my $V = $S->copy->bmodpow($k,$n)->bmul(BTWO)->bmod($n);
    my $Qk = ($ZERO+$Q)->bmodpow($k, $n);
    return ($U, $V, $Qk);
  }
  my $U = BONE->copy;
  my $V = $P->copy;
  my $Qk = $Q->copy;

  return (BZERO->copy, BTWO->copy) if $k == 0;
  $k = Math::BigInt->new("$k") unless ref($k) eq 'Math::BigInt';
  my $kstr = substr($k->as_bin, 2);
  my $bpos = 0;

  if ($Q->is_one) {
    my $Dinverse = $D->copy->bmodinv($n);
    if ($P > BTWO && !$Dinverse->is_nan) {
      # Calculate V_k with U=V_{k+1}
      $U = $P->copy->bmul($P)->bsub(BTWO)->bmod($n);
      while (++$bpos < length($kstr)) {
        if (substr($kstr,$bpos,1)) {
          $V->bmul($U)->bsub($P  )->bmod($n);
          $U->bmul($U)->bsub(BTWO)->bmod($n);
        } else {
          $U->bmul($V)->bsub($P  )->bmod($n);
          $V->bmul($V)->bsub(BTWO)->bmod($n);
        }
      }
      # Crandall and Pomerance eq 3.13: U_n = D^-1 (2V_{n+1} - PV_n)
      $U = $Dinverse * (BTWO*$U - $P*$V);
    } else {
      while (++$bpos < length($kstr)) {
        $U->bmul($V)->bmod($n);
        $V->bmul($V)->bsub(BTWO)->bmod($n);
        if (substr($kstr,$bpos,1)) {
          my $T1 = $U->copy->bmul($D);
          $U->bmul($P)->badd( $V);
          $U->badd($n) if $U->is_odd;
          $U->brsft(BONE);
          $V->bmul($P)->badd($T1);
          $V->badd($n) if $V->is_odd;
          $V->brsft(BONE);
        }
      }
    }
  } else {
    my $qsign = ($Q == -1) ? -1 : 0;
    while (++$bpos < length($kstr)) {
      $U->bmul($V)->bmod($n);
      if    ($qsign ==  1) { $V->bmul($V)->bsub(BTWO)->bmod($n); }
      elsif ($qsign == -1) { $V->bmul($V)->badd(BTWO)->bmod($n); }
      else { $V->bmul($V)->bsub($Qk->copy->blsft(BONE))->bmod($n); }
      if (substr($kstr,$bpos,1)) {
        my $T1 = $U->copy->bmul($D);
        $U->bmul($P)->badd( $V);
        $U->badd($n) if $U->is_odd;
        $U->brsft(BONE);

        $V->bmul($P)->badd($T1);
        $V->badd($n) if $V->is_odd;
        $V->brsft(BONE);

        if ($qsign != 0) { $qsign = -1; }
        else             { $Qk->bmul($Qk)->bmul($Q)->bmod($n); }
      } else {
        if ($qsign != 0) { $qsign = 1; }
        else             { $Qk->bmul($Qk)->bmod($n); }
      }
    }
    if    ($qsign ==  1) { $Qk->bneg; }
    elsif ($qsign == -1) { $Qk = $n->copy->bdec; }
  }
  $U->bmod($n);
  $V->bmod($n);
  return ($U, $V, $Qk);
}

sub is_lucas_pseudoprime {
  my($n) = @_;

  return 0+($n >= 2) if $n < 4;
  return 0 if ($n % 2) == 0 || _is_perfect_square($n);

  my ($P, $Q, $D) = _lucas_selfridge_params($n);
  return 0 if $D == 0;  # We found a divisor in the sequence
  die "Lucas parameter error: $D, $P, $Q\n" if ($D != $P*$P - 4*$Q);

  my($U, $V, $Qk) = lucas_sequence($n, $P, $Q, $n+1);
  return $U->is_zero ? 1 : 0;
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
  my($U, $V, $Qk) = lucas_sequence($n, $P, $Q, $k);

  return 1 if $U->is_zero;
  foreach my $r (0 .. $s-1) {
    return 1 if $V->is_zero;
    if ($r < ($s-1)) {
      $V->bmul($V)->bsub(BTWO*$Qk)->bmod($n);
      $Qk->bmul($Qk)->bmod($n);
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

  # We have to convert n to a bigint or Math::BigInt::GMP's stupid set_si bug
  # (RT 71548) will hit us and make the test $V == $n-2 always return false.
  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';

  my($s, $k) = (0, $n->copy->binc);
  while ($k->is_even && !$k->is_zero) {
    $s++;
    $k->brsft(BONE);
  }

  my($U, $V, $Qk) = lucas_sequence($n, $P, $Q, $k);

  return 1 if $U->is_zero && ($V == BTWO || $V == ($n - BTWO));
  foreach my $r (0 .. $s-2) {
    return 1 if $V->is_zero;
    $V->bmul($V)->bsub(BTWO)->bmod($n);
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

  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';

  my $ZERO = $n->copy->bzero;
  my $TWO = $ZERO->copy->binc->binc;
  my $V = $ZERO + $P;           # V_{k}
  my $W = $ZERO + $P*$P-$TWO;   # V_{k+1}
  my $kstr = substr($n->copy->binc()->as_bin, 2);
  $kstr =~ s/(0*)$//;
  my $s = length($1);
  my $bpos = 0;
  while (++$bpos < length($kstr)) {
    if (substr($kstr,$bpos,1)) {
      $V->bmul($W)->bsub($P  )->bmod($n);
      $W->bmul($W)->bsub($TWO)->bmod($n);
    } else {
      $W->bmul($V)->bsub($P  )->bmod($n);
      $V->bmul($V)->bsub($TWO)->bmod($n);
    }
  }

  return 1 if $V == 2 || $V == ($n-$TWO);
  foreach my $r (0 .. $s-2) {
    return 1 if $V->is_zero;
    $V->bmul($V)->bsub($TWO)->bmod($n);
  }
  return 0;
}

sub is_frobenius_underwood_pseudoprime {
  my($n) = @_;
  return 0+($n >= 2) if $n < 4;

  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';
  return 0 if $n->is_even || _is_perfect_square($n);

  my $ZERO = $n->copy->bzero;
  my $ONE = $ZERO->copy->binc;
  my $TWO = $ONE->copy->binc;
  my $fa = $ZERO + 1;
  my $fb = $ZERO + 2;
  my $x;

  if ($n % 4 == 3) {
    $x = 0;
  } elsif ($n % 6 == 5) {
    $x = 1;
  } else {
    $x = 3;
    $x++ while kronecker($x*$x-4, $n) == 1;
    return 0 if kronecker($x*$x-4, $n) == 0; # gcd(x*x-4,n) is a factor
  }
  my $np1string = substr( $n->copy->binc->as_bin, 2);
  my $np1len = length($np1string);
  my $multiplier = $x+2;
  $multiplier %= $n if $multiplier > $n;
  foreach my $bit (1 .. $np1len-1) {
    my $t = ( $fa * ( $fa * $x + $TWO * $fb ) ) % $n;
    $fb   = ( ( $fb - $fa ) * ( $fb + $fa ) ) % $n;
    $fa   = $t;
    if ( substr( $np1string, $bit, 1 ) ) {
      $t  = $fa * $multiplier + $fb;
      $fb = $fb * $TWO - $fa;
      $fa = $t;
    }
  }
  $fa %= $n;
  $fb %= $n;
  return ($fa == 0 && $fb == (($TWO * $x + 5) % $n)) ? 1 : 0;
}


my $_poly_bignum;
sub _poly_new {
  my @poly = @_;
  push @poly, 0 unless scalar @poly;
  if ($_poly_bignum) {
    @poly = map { (ref $_ eq 'Math::BigInt')
                  ?  $_->copy
                  :  Math::BigInt->new("$_"); } @poly;
  }
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
  my @res = map { $_poly_bignum ? Math::BigInt->bzero : 0 } 0 .. $r-1;

  # convolve(px, py) mod (X^r-1,n)
  my @indices_y = grep { $py->[$_] } (0 .. $py_degree);
  foreach my $ix (0 .. $px_degree) {
    my $px_at_ix = $px->[$ix];
    next unless $px_at_ix;
    if ($_poly_bignum) {
      foreach my $iy (@indices_y) {
        my $rindex = ($ix + $iy) % $r;  # reduce mod X^r-1
        $res[$rindex]->badd($px_at_ix->copy->bmul($py->[$iy]))->bmod($n);
      }
    } else {
      foreach my $iy (@indices_y) {
        my $rindex = ($ix + $iy) % $r;  # reduce mod X^r-1
        $res[$rindex] = ($res[$rindex] + $px_at_ix * $py->[$iy]) % $n;
      }
    }
  }
  # In case we had upper terms go to zero after modulo, reduce the degree.
  pop @res while !$res[-1];
  return \@res;
}

sub _poly_mod_pow {
  my($pn, $power, $r, $mod) = @_;
  my $res = _poly_new(1);
  my $p = $power;

  while ($p) {
    $res = _poly_mod_mul($res, $pn, $r, $mod) if ($p & 1);
    $p >>= 1;
    $pn  = _poly_mod_mul($pn,  $pn, $r, $mod) if $p;
  }
  return $res;
}

sub _test_anr {
  my($a, $n, $r) = @_;
  my $pp = _poly_mod_pow(_poly_new($a, 1), $n, $r, $n);
  $pp->[$n % $r] = (($pp->[$n % $r] || 0) -  1) % $n;  # subtract X^(n%r)
  $pp->[      0] = (($pp->[      0] || 0) - $a) % $n;  # subtract a
  return 0 if scalar grep { $_ } @$pp;
  1;
}

sub is_aks_prime {
  my $n = shift;
  return 0 if $n < 2 || is_power($n);

  my($log2n, $limit);
  if ($n > 2**48) {
    do { require Math::BigFloat; Math::BigFloat->import(); }
      if !defined $Math::BigFloat::VERSION;
    # limit = floor( log2(n) * log2(n) ).  o_r(n) must be larger than this
    my $floatn = Math::BigFloat->new("$n");
    #my $sqrtn = _bigint_to_int($floatn->copy->bsqrt->bfloor);
    # The following line seems to trigger a memory leak in Math::BigFloat::blog
    # (the part where $MBI is copied to $int) if $n is a Math::BigInt::GMP.
    $log2n = $floatn->copy->blog(2);
    $limit = _bigint_to_int( ($log2n * $log2n)->bfloor );
  } else {
    $log2n = log($n)/log(2) + 0.0001;      # Error on large side.
    $limit = int( $log2n*$log2n + 0.0001 );
  }

  my $r = next_prime($limit);
  foreach my $f (@{primes(0,$r-1)}) {
    return 1 if $f == $n;
    return 0 if !($n % $f);
  }

  while ($r < $n) {
    return 0 if !($n % $r);
    #return 1 if $r >= $sqrtn;
    last if znorder($n, $r) > $limit;  # Note the arguments!
    $r = next_prime($r);
  }

  return 1 if $r >= $n;

  # Since r is a prime, phi(r) = r-1
  my $rlimit = (ref($log2n) eq 'Math::BigFloat')
             ? _bigint_to_int( Math::BigFloat->new("$r")->bdec()
                                           ->bsqrt->bmul($log2n)->bfloor)
             : int( (sqrt(($r-1)) * $log2n) + 0.001 );

  $_poly_bignum = 1;
  if ( $n < (MPU_HALFWORD-1) ) {
    $_poly_bignum = 0;
    #$n = _bigint_to_int($n) if ref($n) eq 'Math::BigInt';
  } else {
    $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';
  }

  my $_verbose = Math::Prime::Util::prime_get_config()->{'verbose'};
  print "# aks r = $r  s = $rlimit\n" if $_verbose;
  local $| = 1 if $_verbose > 1;
  for (my $a = 1; $a <= $rlimit; $a++) {
    return 0 unless _test_anr($a, $n, $r);
    print "." if $_verbose > 1;
  }
  print "\n" if $_verbose > 1;

  return 1;
}


sub _basic_factor {
  # MODIFIES INPUT SCALAR
  return ($_[0] == 1) ? () : ($_[0])   if $_[0] < 4;

  my @factors;
  if (ref($_[0]) ne 'Math::BigInt') {
    while ( !($_[0] % 2) ) { push @factors, 2;  $_[0] = int($_[0] / 2); }
    while ( !($_[0] % 3) ) { push @factors, 3;  $_[0] = int($_[0] / 3); }
    while ( !($_[0] % 5) ) { push @factors, 5;  $_[0] = int($_[0] / 5); }
  } else {
    # Without this, the bdivs will try to convert the results to BigFloat
    # and lose precision.
    $_[0]->upgrade(undef) if ref($_[0]) && $_[0]->upgrade();
    if (!Math::BigInt::bgcd($_[0], B_PRIM235)->is_one) {
      while ( $_[0]->is_even)   { push @factors, 2;  $_[0]->brsft(BONE); }
      foreach my $div (3, 5) {
        my ($q, $r) = $_[0]->copy->bdiv($div);
        while ($r->is_zero) {
          push @factors, $div;
          $_[0] = $q;
          ($q, $r) = $_[0]->copy->bdiv($div);
        }
      }
    }
    $_[0] = _bigint_to_int($_[0]) if $_[0] <= ''.~0;
  }

  if ( ($_[0] > 1) && _is_prime7($_[0]) ) {
    push @factors, $_[0];
    $_[0] = 1;
  }
  @factors;
}

sub trial_factor {
  my($n, $maxlim) = @_;
  $maxlim = $n unless defined $maxlim;

  # Don't use _basic_factor here -- they want a trial forced.
  my @factors;
  if ($n < 4) {
    @factors = ($n == 1) ? () : ($n);
    return @factors;
  }
  if (ref($n) ne 'Math::BigInt' || !Math::BigInt::bgcd($n, 30030)->is_one) {
    while ( !($n % 2) ) { push @factors, 2;  $n = int($n / 2); }
    while ( !($n % 3) ) { push @factors, 3;  $n = int($n / 3); }
    while ( !($n % 5) ) { push @factors, 5;  $n = int($n / 5); }
    while ( !($n % 7) ) { push @factors, 7;  $n = int($n / 7); }
    while ( !($n %11) ) { push @factors,11;  $n = int($n /11); }
    while ( !($n %13) ) { push @factors,13;  $n = int($n /13); }
  }
  $n = _bigint_to_int($n) if ref($n) eq 'Math::BigInt' && $n <= ''.~0;
  return @factors if $n < 4;

  my $limit = int(sqrt($n) + 0.001);
  $limit = $maxlim if $limit > $maxlim;

  if (ref($n) eq 'Math::BigInt') {
    my $f = Math::BigInt->new(17);
    $limit = Math::BigInt->new("$limit");
    my @incs = map { Math::BigInt->new($_) } (2, 4, 6, 2, 6, 4, 2, 4);
    SEARCH: while ($f <= $limit) {
      foreach my $finc (@incs) {
        if ($n->copy->bmod($f)->is_zero && $f->bacmp($limit) <= 0) {
          my $sf = ($f <= ''.~0) ? _bigint_to_int($f) : $f;
          do { push @factors, $sf; $n = ($n/$f)->bfloor(); } while (($n % $f) == 0);
          last SEARCH if $n->is_one;
          $limit = int( sqrt($n) + 0.001);
          $limit = $maxlim if $limit > $maxlim;
          $limit = Math::BigInt->new("$limit");
        }
        $f->badd($finc);
      }
    }
  } else {
    my $f = 17;
    SEARCH: while ($f <= $limit) {
      foreach my $finc (2, 4, 6, 2, 6, 4, 2, 4) {
        if ( (($n % $f) == 0) && ($f <= $limit) ) {
          do { push @factors, $f; $n = int($n/$f); } while (($n % $f) == 0);
          last SEARCH if $n == 1;
          $limit = int( sqrt($n) + 0.001);
          $limit = $maxlim if $limit > $maxlim;
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
  sub { prho_factor   (shift,    8*1024, 3) },
  sub { pminus1_factor(shift,    10_000); },
  sub { pbrent_factor (shift,   32*1024, 1) },
  sub { pminus1_factor(shift, 1_000_000); },
  sub { pbrent_factor (shift,  512*1024, 7) },
  sub { ecm_factor    (shift,     1_000,   5_000, 10) },
  sub { pminus1_factor(shift, 4_000_000); },
  sub { pbrent_factor (shift,  512*1024, 11) },
  sub { ecm_factor    (shift,    10_000,  50_000, 10) },
  sub { holf_factor   (shift, 256*1024, $_holf_r); $_holf_r += 256*1024; },
  sub { pminus1_factor(shift,20_000_000); },
  sub { ecm_factor    (shift,   100_000, 800_000, 10) },
  sub { holf_factor   (shift, 512*1024, $_holf_r); $_holf_r += 512*1024; },
  sub { pbrent_factor (shift, 2048*1024, 13) },
  sub { holf_factor   (shift, 2048*1024, $_holf_r); $_holf_r += 2048*1024; },
  sub { ecm_factor    (shift, 1_000_000, 1_000_000, 10) },
  sub { pminus1_factor(shift, 100_000_000, 500_000_000); },
);

sub factor {
  my($n) = @_;
  _validate_positive_integer($n);

  return trial_factor($n) if $n < 1_000_000;

  $n = $n->copy if ref($n) eq 'Math::BigInt';

  my @factors;
  # Use 'n=int($n/7)' instead of 'n/=7' to not "upgrade" n to a Math::BigFloat.
  if (ref($n) eq 'Math::BigInt') {
    while ($n->is_even) { push @factors,  2;  $n->brsft(BONE); }
    if (!Math::BigInt::bgcd($n, "3234846615")->is_one) {
      foreach my $div (3, 5, 7, 11, 13, 17, 19, 23, 29) {
        my ($q, $r) = $n->copy->bdiv($div);
        while ($r->is_zero) {
          push @factors, $div;
          $n = $q;
          ($q, $r) = $n->copy->bdiv($div);
        }
      }
    }
  } else {
    while (($n %  2) == 0) { push @factors,  2;  $n = int($n /  2); }
    while (($n %  3) == 0) { push @factors,  3;  $n = int($n /  3); }
    while (($n %  5) == 0) { push @factors,  5;  $n = int($n /  5); }
    while (($n %  7) == 0) { push @factors,  7;  $n = int($n /  7); }
    while (($n % 11) == 0) { push @factors, 11;  $n = int($n / 11); }
    while (($n % 13) == 0) { push @factors, 13;  $n = int($n / 13); }
    while (($n % 17) == 0) { push @factors, 17;  $n = int($n / 17); }
    while (($n % 19) == 0) { push @factors, 19;  $n = int($n / 19); }
    while (($n % 23) == 0) { push @factors, 23;  $n = int($n / 23); }
    while (($n % 29) == 0) { push @factors, 29;  $n = int($n / 29); }
  }
  if ($n < (31*31)) {
    push @factors, $n  if $n != 1;
    return @factors;
  }

  my @nstack = ($n);
  while (@nstack) {
    $n = pop @nstack;
    # Don't use bignum on $n if it has gotten small enough.
    $n = _bigint_to_int($n) if ref($n) eq 'Math::BigInt' && $n <= ''.~0;
    #print "Looking at $n with stack ", join(",",@nstack), "\n";
    while ( ($n >= (31*31)) && !_is_prime7($n) ) {
      my @ftry;
      $_holf_r = 1;
      foreach my $sub (@_fsublist) {
        last if scalar @ftry >= 2;
        @ftry = $sub->($n);
      }
      if (scalar @ftry > 1) {
        #print "  split into ", join(",",@ftry), "\n";
        $n = shift @ftry;
        $n = _bigint_to_int($n) if ref($n) eq 'Math::BigInt' && $n <= ''.~0;
        push @nstack, @ftry;
      } else {
        #warn "trial factor $n\n";
        push @factors, trial_factor($n);
        #print "  trial into ", join(",",@factors), "\n";
        $n = 1;
        last;
      }
    }
    push @factors, $n  if $n != 1;
  }
  @factors = sort {$a<=>$b} @factors;
  return @factors;
}

sub _found_factor {
  my($f, $n, $what, @factors) = @_;
  if ($f == 1 || $f == $n) {
    push @factors, $n;
  } else {
    # Perl 5.6.2 needs things spelled out for it.
    my $f2 = (ref($n) eq 'Math::BigInt') ? $n->copy->bdiv($f)->as_int
                                         : int($n/$f);
    push @factors, $f;
    push @factors, $f2;
    croak "internal error in $what" unless $f * $f2 == $n;
    # MPU::GMP prints this type of message if verbose, so do the same.
    print "$what found factor $f\n" if Math::Prime::Util::prime_get_config()->{'verbose'} > 0;
  }
  @factors;
}

# TODO:
sub squfof_factor { trial_factor(@_) }

sub prho_factor {
  my($n, $rounds, $pa) = @_;
  $rounds = 4*1024*1024 unless defined $rounds;
  $pa = 3 unless defined $pa;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  my $inloop = 0;
  my $U = 7;
  my $V = 7;

  if ( ref($n) eq 'Math::BigInt' ) {

    $pa = Math::BigInt->new("$pa");
    $U = $n->copy->bzero->badd($U);
    $V = $n->copy->bzero->badd($V);
    for my $i (1 .. $rounds) {
      # Would use bmuladd here, but old Math::BigInt's barf with scalar $pa.
      $U->bmul($U)->badd($pa)->bmod($n);
      $V->bmul($V)->badd($pa);
      $V->bmul($V)->badd($pa)->bmod($n);
      my $f = Math::BigInt::bgcd($U-$V, $n);
      if ($f->bacmp($n) == 0) {
        last if $inloop++;  # We've been here before
      } elsif (!$f->is_one) {
        return _found_factor($f, $n, "prho", @factors);
      }
    }

  } elsif ($n < MPU_HALFWORD) {

    for my $i (1 .. $rounds) {
      $U = ($U * $U + $pa) % $n;
      $V = ($V * $V + $pa) % $n;
      $V = ($V * $V + $pa) % $n;
      my $f = _gcd_ui( $U-$V,  $n );
      if ($f == $n) {
        last if $inloop++;  # We've been here before
      } elsif ($f != 1) {
        return _found_factor($f, $n, "prho", @factors);
      }
    }

  } else {

    for my $i (1 .. $rounds) {
      if ($n <= (~0 >> 1)) {
       $U = _mulmod($U, $U, $n);  $U += $pa;  $U -= $n if $U >= $n;
       $V = _mulmod($V, $V, $n);  $V += $pa;  # Let the mulmod handle it
       $V = _mulmod($V, $V, $n);  $V += $pa;  $V -= $n if $V >= $n;
      } else {
       $U = _mulmod($U, $U, $n); $U=$n-$U;  $U = ($pa>=$U) ? $pa-$U : $n-$U+$pa;
       $V = _mulmod($V, $V, $n); $V=$n-$V;  $V = ($pa>=$V) ? $pa-$V : $n-$V+$pa;
       $V = _mulmod($V, $V, $n); $V=$n-$V;  $V = ($pa>=$V) ? $pa-$V : $n-$V+$pa;
      }
      my $f = _gcd_ui( $U-$V,  $n );
      if ($f == $n) {
        last if $inloop++;  # We've been here before
      } elsif ($f != 1) {
        return _found_factor($f, $n, "prho", @factors);
      }
    }

  }
  push @factors, $n;
  @factors;
}

sub pbrent_factor {
  my($n, $rounds, $pa) = @_;
  $rounds = 4*1024*1024 unless defined $rounds;
  $pa = 3 unless defined $pa;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  my $Xi = 2;
  my $Xm = 2;

  if ( ref($n) eq 'Math::BigInt' ) {

    # Same code as the GMP version, but runs *much* slower.  Even with
    # Math::BigInt::GMP it's >200x slower.  With the default Calc backend
    # it's thousands of times slower.
    my $inner = 256;
    my $zero = $n->copy->bzero;
    my $saveXi;
    my $f;
    $Xi = $zero->copy->badd($Xi);
    $Xm = $zero->copy->badd($Xm);
    my $r = 1;
    while ($rounds > 0) {
      my $rleft = ($r > $rounds) ? $rounds : $r;
      while ($rleft > 0) {
        my $dorounds = ($rleft > $inner) ? $inner : $rleft;
        my $m = $zero->copy->bone;
        $saveXi = $Xi->copy;
        foreach my $i (1 .. $dorounds) {
          $Xi->bmul($Xi)->badd($pa)->bmod($n);
          $m->bmul($Xi - $Xm);
        }
        $rleft -= $dorounds;
        $rounds -= $dorounds;
        $m->bmod($n);
        $f = Math::BigInt::bgcd( $m,  $n);
        last if $f != 1;
      }
      if ($f == 1) {
        $r *= 2;
        $Xm = $Xi->copy;
        next;
      }
      if ($f == $n) {  # back up to determine the factor
        $Xi = $saveXi->copy;
        do {
          $Xi->bmul($Xi)->badd($pa)->bmod($n);
          $f = Math::BigInt::bgcd($Xm-$Xi, $n);
        } while ($f != 1 && $r-- != 0);
        last if $f == 1 || $f == $n;
      }
      return _found_factor($f, $n, "pbrent", @factors);
    }

  } elsif ($n < MPU_HALFWORD) {

    for my $i (1 .. $rounds) {
      $Xi = ($Xi * $Xi + $pa) % $n;
      my $f = _gcd_ui($Xm-$Xi, $n);
      return _found_factor($f, $n, "pbrent", @factors) if $f != 1 && $f != $n;
      $Xm = $Xi if ($i & ($i-1)) == 0;  # i is a power of 2
    }

  } else {

    for my $i (1 .. $rounds) {
      # Xi^2+a % n
      $Xi = _mulmod($Xi, $Xi, $n);
      $Xi = (($n-$Xi) > $pa)  ?  $Xi+$pa  :  $Xi+$pa-$n;
      my $f = _gcd_ui($Xm-$Xi, $n);
      return _found_factor($f, $n, "pbrent", @factors) if $f != 1 && $f != $n;
      $Xm = $Xi if ($i & ($i-1)) == 0;  # i is a power of 2
    }

  }
  push @factors, $n;
  @factors;
}

sub pminus1_factor {
  my($n, $B1, $B2) = @_;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if ( ref($n) ne 'Math::BigInt' ) {
    # Stage 1 only
    $B1 = 10_000_000 unless defined $B1;
    my $pa = 2;
    my $f = 1;
    my($pc_beg, $pc_end, @bprimes);
    $pc_beg = 2;
    $pc_end = $pc_beg + 100_000;
    my $sqrtb1 = int(sqrt($B1));
    while (1) {
      $pc_end = $B1 if $pc_end > $B1;
      @bprimes = @{ primes($pc_beg, $pc_end) };
      foreach my $q (@bprimes) {
        my $k = $q;
        if ($q <= $sqrtb1) {
          my $kmin = int($B1 / $q);
          while ($k <= $kmin) { $k *= $q; }
        }
        $pa = _powmod($pa, $k, $n);
        if ($pa == 0) { push @factors, $n; return @factors; }
        my $f = _gcd_ui( $pa-1, $n );
        return _found_factor($f, $n, "pminus1", @factors) if $f != 1;
      }
      last if $pc_end >= $B1;
      $pc_beg = $pc_end+1;
      $pc_end += 500_000;
    }
    push @factors, $n;
    return @factors;
  }

  # Stage 2 isn't really any faster than stage 1 for the examples I've tried.
  # Perl's overhead is greater than the savings of multiply vs. powmod

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

  my $one = $n->copy->bone;
  my ($j, $q, $saveq) = (32, 2, 2);
  my $t = $one->copy;
  my $pa = $one->copy->binc();
  my $savea = $pa->copy;
  my $f = $one->copy;
  my($pc_beg, $pc_end, @bprimes);

  $pc_beg = 2;
  $pc_end = $pc_beg + 100_000;
  while (1) {
    $pc_end = $B1 if $pc_end > $B1;
    @bprimes = @{ primes($pc_beg, $pc_end) };
    foreach my $q (@bprimes) {
      my($k, $kmin) = ($q, int($B1 / $q));
      while ($k <= $kmin) { $k *= $q; }
      $t *= $k;                         # accumulate powers for a
      if ( ($j++ % 64) == 0) {
        next if $pc_beg > 2 && ($j-1) % 256;
        $pa->bmodpow($t, $n);
        $t = $one->copy;
        if ($pa == 0) { push @factors, $n; return @factors; }
        $f = Math::BigInt::bgcd( $pa-1, $n );
        last if $f == $n;
        return _found_factor($f, $n, "pminus1", @factors) unless $f->is_one;
        $saveq = $q;
        $savea = $pa->copy;
      }
    }
    $q = $bprimes[-1];
    last if !$f->is_one || $pc_end >= $B1;
    $pc_beg = $pc_end+1;
    $pc_end += 500_000;
  }
  undef @bprimes;
  $pa->bmodpow($t, $n);
  if ($pa == 0) { push @factors, $n; return @factors; }
  $f = Math::BigInt::bgcd( $pa-1, $n );
  if ($f == $n) {
    $q = $saveq;
    $pa = $savea->copy;
    while ($q <= $B1) {
      my ($k, $kmin) = ($q, int($B1 / $q));
      while ($k <= $kmin) { $k *= $q; }
      $pa->bmodpow($k, $n);
      my $f = Math::BigInt::bgcd( $pa-1, $n );
      if ($f == $n) { push @factors, $n; return @factors; }
      last if !$f->is_one;
      $q = next_prime($q);
    }
  }
  # STAGE 2
  if ($f->is_one && $B2 > $B1) {
    my $bm = $pa->copy;
    my $b = $one->copy;
    my @precomp_bm;
    $precomp_bm[0] = ($bm * $bm) % $n;
    foreach my $j (1..19) {
      $precomp_bm[$j] = ($precomp_bm[$j-1] * $bm * $bm) % $n;
    }
    $pa->bmodpow($q, $n);
    my $j = 1;
    $pc_beg = $q+1;
    $pc_end = $pc_beg + 100_000;
    while (1) {
      $pc_end = $B2 if $pc_end > $B2;
      @bprimes = @{ primes($pc_beg, $pc_end) };
      foreach my $i (0 .. $#bprimes) {
        my $diff = $bprimes[$i] - $q;
        $q = $bprimes[$i];
        my $qdiff = ($diff >> 1) - 1;
        if (!defined $precomp_bm[$qdiff]) {
          $precomp_bm[$qdiff] = $bm->copy->bmodpow($diff, $n);
        }
        $pa->bmul($precomp_bm[$qdiff])->bmod($n);
        if ($pa == 0) { push @factors, $n; return @factors; }
        $b->bmul($pa-1);
        if (($j++ % 128) == 0) {
          $b->bmod($n);
          $f = Math::BigInt::bgcd( $b, $n );
          last if !$f->is_one;
        }
      }
      last if !$f->is_one || $pc_end >= $B2;
      $pc_beg = $pc_end+1;
      $pc_end += 500_000;
    }
    $f = Math::BigInt::bgcd( $b, $n );
  }
  return _found_factor($f, $n, "pminus1", @factors);
}

sub holf_factor {
  my($n, $rounds, $startrounds) = @_;
  $rounds = 64*1024*1024 unless defined $rounds;
  $startrounds = 1 unless defined $startrounds;
  $startrounds = 1 if $startrounds < 1;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if ( ref($n) eq 'Math::BigInt' ) {
    for my $i ($startrounds .. $rounds) {
      my $ni = $n->copy->bmul($i);
      my $s = $ni->copy->bsqrt->bfloor->as_int;
      if ($s * $s == $ni) {
        # s^2 = n*i, so m = s^2 mod n = 0.  Hence f = GCD(n, s) = GCD(n, n*i)
        my $f = Math::BigInt::bgcd($ni, $n);
        return _found_factor($f, $n, "HOLF", @factors);
      }
      $s->binc;
      my $m = ($s * $s) - $ni;
      # Check for perfect square
      my $mc = _bigint_to_int($m & 31);
      next unless $mc==0||$mc==1||$mc==4||$mc==9||$mc==16||$mc==17||$mc==25;
      my $f = $m->copy->bsqrt->bfloor->as_int;
      next unless ($f*$f) == $m;
      $f = Math::BigInt::bgcd( ($s > $f) ? $s-$f : $f-$s,  $n);
      return _found_factor($f, $n, "HOLF ($i rounds)", @factors);
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
  $rounds = 64*1024*1024 unless defined $rounds;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if ( ref($n) eq 'Math::BigInt' ) {
    my $pa = $n->copy->bsqrt->bfloor->as_int;
    return _found_factor($pa, $n, "Fermat", @factors) if $pa*$pa == $n;
    $pa++;
    my $b2 = $pa*$pa - $n;
    my $lasta = $pa + $rounds;
    while ($pa <= $lasta) {
      my $mc = _bigint_to_int($b2 & 31);
      if ($mc==0||$mc==1||$mc==4||$mc==9||$mc==16||$mc==17||$mc==25) {
        my $s = $b2->copy->bsqrt->bfloor->as_int;
        if ($s*$s == $b2) {
          my $i = $pa-($lasta-$rounds)+1;
          return _found_factor($pa - $s, $n, "Fermat ($i rounds)", @factors);
        }
      }
      $pa++;
      $b2 = $pa*$pa-$n;
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

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

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
  my @bprimes = @{ primes(3, $B1) };
  foreach my $q (@bprimes) {
    last if $q > $sqrt_b1;
    my($k,$kmin) = ($q, int($B1/$q));
    while ($k <= $kmin) { $k *= $q; }
    $q = $k;
  }
  my @b2primes = ($B2 > $B1) ? @{primes($B1+1, $B2)} : ();
  my $irandf = Math::Prime::Util::RandomPrimes::get_randf();

  foreach my $curve (1 .. $ncurves) {
    my $sigma = $irandf->($n-1-6) + 6;
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
  my($n) = @_;
  _validate_positive_integer($n);

  # In scalar context, returns sigma_0(n).  Very fast.
  return Math::Prime::Util::divisor_sum($n,0) unless wantarray;
  return ($n == 0) ? (0,1) : (1)  if $n <= 1;

  my %all_factors;
  my @factors = Math::Prime::Util::factor($n);
  return (1,$n) if scalar @factors == 1;

  if (ref($n) eq 'Math::BigInt') {
    foreach my $f1 (@factors) {
      my $big_f1 = Math::BigInt->new("$f1");
      my @to_add = map { ($_ <= ''.~0) ? _bigint_to_int($_) : $_ }
                   grep { $_ < $n }
                   map { $big_f1 * $_ }
                   keys %all_factors;
      undef @all_factors{ $f1, @to_add };
    }
  } else {
    foreach my $f1 (@factors) {
      my @to_add = grep { $_ < $n }
                   map { $f1 * $_ }
                   keys %all_factors;
      undef @all_factors{ $f1, @to_add };
    }
  }
  # Add 1 and n
  undef $all_factors{1};
  undef $all_factors{$n};
  my @divisors = sort {$a<=>$b} keys %all_factors;
  return @divisors;
}


sub chebyshev_theta {
  my($n) = @_;
  my $sum = 0.0;
  for (my $p = 2; $p <= $n; $p = next_prime($p)) {
    $sum += log($p);
  }
  return $sum;
}

sub chebyshev_psi {
  my($n) = @_;
  return 0 if $n <= 1;

  my ($sum, $p, $logn, $sqrtn) = (0.0, 2, log($n), int(sqrt($n)));

  for ( ; $p <= $sqrtn; $p = next_prime($p)) {
    my $logp = log($p);
    $sum += $logp * int($logn/$logp+1e-15);
  }

  for ( ; $p <= $n; $p = next_prime($p)) {
    $sum += log($p);
  }

  return $sum;
}

sub ExponentialIntegral {
  my($x) = @_;
  return - MPU_INFINITY if $x == 0;
  return 0              if $x == - MPU_INFINITY;
  return MPU_INFINITY   if $x == MPU_INFINITY;

  # Gotcha -- MPFR decided to make negative inputs return NaN.  Grrr.
  if ($x > 0 && _MPFR_available()) {
    my $wantbf = 0;
    my $xdigits = 17;
    if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
      do { require Math::BigFloat; Math::BigFloat->import(); }
        if !defined $Math::BigFloat::VERSION;
      $x = Math::BigFloat->new("$x") if ref($x) ne 'Math::BigFloat';
      $wantbf = _find_big_acc($x);
      $xdigits = $wantbf;
    }
    my $rnd = 0;  # MPFR_RNDN;
    my $bit_precision = int($xdigits * 3.322) + 4;
    my $rx = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rx, $bit_precision);
    Math::MPFR::Rmpfr_set_str($rx, "$x", 10, $rnd);
    my $eix = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($eix, $bit_precision);
    Math::MPFR::Rmpfr_eint($eix, $rx, $rnd);
    my $strval = Math::MPFR::Rmpfr_get_str($eix, 10, 0, $rnd);
    Math::MPFR::Rmpfr_free_cache();
    return ($wantbf)  ?  Math::BigFloat->new($strval,$wantbf)  :  0.0 + $strval;
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
    $y = CONST_EULER-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
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
  my($x) = @_;
  return 0              if $x == 0;
  return - MPU_INFINITY if $x == 1;
  return MPU_INFINITY   if $x == MPU_INFINITY;
  croak "Invalid input to LogarithmicIntegral:  x must be > 0" if $x <= 0;

  # Remember MPFR eint doesn't handle negative inputs
  if ($x >= 1 && _MPFR_available()) {
    my $wantbf = 0;
    my $xdigits = 18;
    if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
      $wantbf = _find_big_acc($x);
      $xdigits = $wantbf;
    }
    $xdigits += length(int(log(0.0+"$x"))) + 1;
    my $rnd = 0;  # MPFR_RNDN;
    my $bit_precision = int($xdigits * 3.322) + 4;
    my $rx = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rx, $bit_precision);
    Math::MPFR::Rmpfr_set_str($rx, "$x", 10, $rnd);
    Math::MPFR::Rmpfr_log($rx, $rx, $rnd);
    my $lix = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($lix, $bit_precision);
    Math::MPFR::Rmpfr_eint($lix, $rx, $rnd);
    my $strval = Math::MPFR::Rmpfr_get_str($lix, 10, 0, $rnd);
    Math::MPFR::Rmpfr_free_cache();
    return ($wantbf)  ?  Math::BigFloat->new($strval,$wantbf)  :  0.0 + $strval;
  }

  if ($x == 2) {
    my $li2const = (ref($x) eq 'Math::BigFloat') ? Math::BigFloat->new(CONST_LI2) : 0.0+CONST_LI2;
    return $li2const;
  }

  $x = Math::BigFloat->new("$x") if defined $bignum::VERSION && ref($x) ne 'Math::BigFloat';

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
    my $xlen = length($x->bfloor->bstr());
    $xdigits = $xlen if $xdigits < $xlen;
    $finalacc = $xdigits;
    $xdigits += length(int(log(0.0+"$x"))) + 1;
    $tol = Math::BigFloat->new(10)->bpow(-$xdigits);
    $x->accuracy($xdigits);
  }
  my $logx = $xdigits ? $x->copy->blog(undef,$xdigits) : log($x);

  if ($x > 1e16) {
    my $invx = 1.0 / $logx;
    # n = 0  =>  0!/(logx)^0 = 1/1 = 1
    # n = 1  =>  1!/(logx)^1 = 1/logx
    my $term = $invx;
    my $sum = 1.0 + $term;
    for my $n (2 .. 200) {
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
    my $val = $x * $invx * $sum;
    $val->accuracy($finalacc) if $xdigits;
    return $val;
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
    my $eulerconst = (ref($x) eq 'Math::BigFloat') ? Math::BigFloat->new(CONST_EULER) : 0.0+CONST_EULER;
    my $val = $eulerconst + log($logx) + $sum;
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

  # Use MPFR if possible.
  if (_MPFR_available()) {
    my $wantbf = 0;
    my $xdigits = 17;
    if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
      do { require Math::BigFloat; Math::BigFloat->import(); }
        if !defined $Math::BigFloat::VERSION;
      if (ref($x) eq 'Math::BigInt') {
        my $xacc = $x->accuracy();
        $x = Math::BigFloat->new($x);
        $x->accuracy($xacc) if $xacc;
      }
      $x = Math::BigFloat->new("$x") if ref($x) ne 'Math::BigFloat';
      $wantbf = _find_big_acc($x);
      $xdigits = $wantbf;
    }
    my $rnd = 0;  # MPFR_RNDN;
    my $bit_precision = int($xdigits * 3.322) + 7;
    my $rx = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rx, $bit_precision);
    Math::MPFR::Rmpfr_set_str($rx, "$x", 10, $rnd);
    my $zetax = Math::MPFR->new();
    # Add more bits to account for the leading zeros.
    my $extra_bits = int(abs($x));
    Math::MPFR::Rmpfr_set_prec($zetax, $bit_precision + $extra_bits);
    Math::MPFR::Rmpfr_zeta($zetax, $rx, $rnd);
    Math::MPFR::Rmpfr_sub_ui($zetax, $zetax, 1, $rnd);
    my $strval = Math::MPFR::Rmpfr_get_str($zetax, 10, $xdigits, $rnd);
    return ($wantbf)  ?  Math::BigFloat->new($strval,$wantbf)  :  0.0 + $strval;
  }

  if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
    # No MPFR, BigFloat
    require Math::Prime::Util::ZetaBigFloat;
    return Math::Prime::Util::ZetaBigFloat::RiemannZeta($x);
  }

  # No MPFR, no BigFloat.
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

  # Use MPFR if possible.
  if (_MPFR_available()) {
    my $wantbf = 0;
    my $xdigits = 17;
    if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
      do { require Math::BigFloat; Math::BigFloat->import(); }
        if !defined $Math::BigFloat::VERSION;
      if (ref($x) eq 'Math::BigInt') {
        my $xacc = $x->accuracy();
        $x = Math::BigFloat->new($x);
        $x->accuracy($xacc) if $xacc;
      }
      $x = Math::BigFloat->new("$x") if ref($x) ne 'Math::BigFloat';
      $wantbf = _find_big_acc($x);
      $xdigits = $wantbf;
    }
    my $rnd = 0;  # MPFR_RNDN;
    my $bit_precision = int($xdigits * 3.322) + 8;  # Add some extra

    my $rlogx = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rlogx, $bit_precision);
    Math::MPFR::Rmpfr_set_str($rlogx, "$x", 10, $rnd);
    Math::MPFR::Rmpfr_log($rlogx, $rlogx, $rnd);

    my $rpart_term = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rpart_term, $bit_precision);
    Math::MPFR::Rmpfr_set_str($rpart_term, "1", 10, $rnd);

    my $rzeta = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rzeta, $bit_precision);
    my $rterm = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rterm, $bit_precision);

    my $rsum = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rsum, $bit_precision);
    Math::MPFR::Rmpfr_set_str($rsum, "1", 10, $rnd);

    my $rstop = Math::MPFR->new();
    Math::MPFR::Rmpfr_set_prec($rstop, $bit_precision);
    Math::MPFR::Rmpfr_set_str($rstop, "1e-$xdigits", 10, $rnd);

    for my $k (1 .. 10000) {
      Math::MPFR::Rmpfr_mul($rpart_term, $rpart_term, $rlogx, $rnd);
      Math::MPFR::Rmpfr_div_ui($rpart_term, $rpart_term, $k, $rnd);

      Math::MPFR::Rmpfr_zeta_ui($rzeta, $k+1, $rnd);
      Math::MPFR::Rmpfr_sub_ui($rzeta, $rzeta, 1, $rnd);
      Math::MPFR::Rmpfr_mul_ui($rzeta, $rzeta, $k, $rnd);
      Math::MPFR::Rmpfr_add_ui($rzeta, $rzeta, $k, $rnd);
      Math::MPFR::Rmpfr_div($rterm, $rpart_term, $rzeta, $rnd);

      last if Math::MPFR::Rmpfr_less_p($rterm, $rstop);
      Math::MPFR::Rmpfr_add($rsum, $rsum, $rterm, $rnd);
    }
    my $strval = Math::MPFR::Rmpfr_get_str($rsum, 10, $xdigits, $rnd);
    return ($wantbf)  ?  Math::BigFloat->new($strval,$wantbf)  :  0.0 + $strval;
  }

  if (defined $bignum::VERSION || ref($x) =~ /^Math::Big/) {
    require Math::Prime::Util::ZetaBigFloat;
    return Math::Prime::Util::ZetaBigFloat::RiemannR($x);
  }


  my $tol = 1e-16;
  my $sum = 0.0;
  my($y, $t);
  my $c = 0.0;

  $y = 1.0-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
  my $flogx = log($x);
  my $part_term = 1.0;
  for my $k (1 .. 10000) {
    # Small k from table, larger k from function
    my $zeta = ($k <= $#_Riemann_Zeta_Table)
               ? $_Riemann_Zeta_Table[$k+1-2]
               : RiemannZeta($k+1);
    $part_term *= $flogx / $k;
    my $term = $part_term / ($k + $k * $zeta);
    $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    last if $term < ($tol * $sum);
  }
  return $sum;
}

sub forpart {
  my($sub, $n, $rhash) = @_;
  _validate_positive_integer($n);
  my($mina, $maxa, $minn, $maxn) = (0,$n,0,$n);
  if (defined $rhash) {
    croak "forpart second argument must be a hash reference"
      unless ref($rhash) eq 'HASH';
    $minn = $maxn = $rhash->{n} if defined $rhash->{n};
    $mina = $rhash->{amin} if defined $rhash->{amin};
    $maxa = $rhash->{amax} if defined $rhash->{amax};
    $minn = $rhash->{nmin} if defined $rhash->{nmin};
    $maxn = $rhash->{nmax} if defined $rhash->{nmax};
   _validate_positive_integer($mina);
   _validate_positive_integer($maxa);
   _validate_positive_integer($minn);
   _validate_positive_integer($maxn);
  }
  my @x = (1) x ($n+1);
  my $m = ($n > 0) ? 1 : 0;
  my $h = 1;
  $x[1] = $n;
  while (1) {
    $sub->(@x[1 .. $m])
      if $m >= $minn && $m <= $maxn && $x[$m] >= $mina && $x[1] <= $maxa;
    last if $x[1] <= 1;
    if ($x[$h] == 2) {
      $m++; $x[$h] = 1; $h--; }
    else {
      my $r = $x[$h] - 1;
      my $t = $m - $h + 1;
      $x[$h] = $r;
      while ($t >= $r) { $h++; $x[$h] = $r; $t -= $r; }
      if ($t == 0) { $m = $h }
      else         { $m = $h+1; if ($t > 1) { $h++; $x[$h] = $t; } }
    }
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

Version 0.41


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
big number support if L<Math::Prime::Util::GMP> isn't available.  This is a
subset of L<Math::Prime::Util>'s functionality.

All routines should work with native integers or multi-precision numbers.  To
enable big numbers, use bigint or bignum:

    use bigint;
    say prime_count_approx(1000000000000000000000000)'
    # says 18435599767347543283712

This is still experimental, and some functions will be very slow.  The
L<Math::Prime::Util::GMP> module has much faster versions of many of these
functions.  Alternately, L<Math::Pari> has a lot of these types of functions.


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

Copyright 2012-2014 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
