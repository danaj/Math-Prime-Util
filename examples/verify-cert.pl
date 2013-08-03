#!/usr/bin/env perl
use warnings;
use strict;
use Math::BigInt lib=>"GMP,Pari";
use Math::Prime::Util qw/:all/;
use Time::HiRes qw(gettimeofday tv_interval);
use Getopt::Long;
$|++;

# MPU and PRIMO certificate verification.
# Written by Dana Jacobsen, 2013.

# Requires Math::Prime::Util v0.30 or later.
# Will be very slow without Math:::Prime::Util::GMP for EC operations.

# Exits with:
#  0  all numbers verified prime
#  1  at least one number verified composite
#  2  incorrect or incomplete conditions.  Cannot verify.
#  3  certificate file cannot be parsed or no number found

# The candidate number is always checked against is_prime first.  That
# performs an extra-strong Lucas pseudoprime test followed by at least
# one additional M-R test using a random base.

my $verbose = 2;
my $quiet;
my $verb;
my $timing;
GetOptions("verbose+" => \$verb,
           "quiet" => \$quiet,
           "timing" => \$timing,
          ) or die "Error in option parsing\n";
$verbose = $verb if defined $verb;
$verbose = 0 if $quiet;

sub error ($) {
  my $message = shift;
  warn "\n$message\n" if $verbose;
  exit(3);  # error in certificate
}

sub fail ($) {
  my $message = shift;
  warn "\n$message\n" if $verbose;
  exit(2);  # Failed a condition
}

my $orig_N;
my $N;
my %parts;  # Map of "N is prime if Q is prime"
my %proof_funcs = (
  ECPP        =>  \&prove_ecpp,    # Standard ECPP proof
  ECPP3       =>  \&prove_ecpp3,   # Primo type 3
  ECPP4       =>  \&prove_ecpp4,   # Primo type 4
  BLS15       =>  \&prove_bls15,   # basic n+1, includes Primo type 2
  BLS3        =>  \&prove_bls3,    # basic n-1
  BLS5        =>  \&prove_bls5,    # much better n-1
  SMALL       =>  \&prove_small,   # n <= 2^64
  POCKLINGTON =>  \&prove_pock,    # simple n-1, Primo type 1
  LUCAS       =>  \&prove_lucas,   # n-1 completely factored
);
my $smallval = Math::BigInt->new(2)->bpow(64);
my $step = 1;
my $base = 10;
my $cert_type = 'Unknown';
my $start_time;

while (<>) {
  next if /^\s*#/ or /^\s*$/;  # Skip comments and blank lines

  chomp;

  if (/^\[(\S+) - Primality Certificate\]/) {
    error "Unknown certificate type: $1" unless $1 eq 'MPU' || $1 eq 'PRIMO';
    $cert_type = $1;
    next;
  }

  if ( ($cert_type eq 'PRIMO' && /^\[Candidate\]/) || ($cert_type eq 'MPU' && /^Proof for:/) ) {
    if (defined $N) {
      # Done with this number, starting the next.
      print " " x 60, "\r" if $verbose == 2;
      if (final_verify($N)) {
        print "PRIME\n" if $verbose;
      } else {
        print "NOT PROVEN\n" if $verbose;
        exit(2);
      }
      undef $N;
      undef %parts;
      $step = 1;
    }
    if ($cert_type eq 'PRIMO') {
      ($N) = primo_read_vars('Candidate', qw/N/);
    } else {
      ($N) = read_vars('Proof for', qw/N/);
    }
    $start_time = [gettimeofday];
    $orig_N = $N;
    if    ($verbose == 1) {  print "N $N";  }
    elsif ($verbose == 2) {  print "$N\n";  }
    if (!is_prime($N)) {
      print "COMPOSITE\n" if $verbose;
      exit(1);
    }
    next;
  }

  if ($cert_type eq 'PRIMO') {
    if (/^Type\s*=\s*(\d+)/) {
      my $type = $1;
      error("Starting type without telling me the N value!") unless defined $N;
      if ($type == 4) {
        my ($n, $f) = verify_ecpp4( $N, primo_read_vars('4', qw/S R J T/) );
        $N = $f;
      } elsif ($type == 3) {
        my ($n, $f) = verify_ecpp3( $N, primo_read_vars('3', qw/S R A B T/) );
        $N = $f;
      } elsif ($type == 2) {
        my ($s, $r, $q) = primo_read_vars('2', qw/S R Q/);
        my $p = ($q->is_odd()) ? 2 : 1;
        my ($n, $f) = verify_bls15( $N, $r, $p, $q );
        $N = $f;
      } elsif ($type == 1) {
        my ($s, $r, $b) = primo_read_vars('1', qw/S R B/);
        fail "Type 1: $N failed SR + 1 = N" unless $s*$r+1 == $N;
        my ($n, $f) = verify_pock( $N, $r, $b );   # S = (N-1)/r
        $N = $f;
      } elsif ($type == 0) {
        # Final
      } else {
        error "Unknown type: $type";
      }
      if    ($verbose == 1) {  print ".";  }
      elsif ($verbose == 2) {  printf "step %2d: %4d digits  type %d\r", $step++, length($N), $type; }
    }
  } elsif ($cert_type eq 'MPU') {
    if (/^Base (\d+)/) {
      $base = $1;
      error "Invalid base: $base" unless $base == 10 || $base == 16 || $base == 62;
      error "Sorry, only base 10 implemented in this version" unless $base == 10;
    } elsif (/^Type (.*?)\s*$/) {
      error("Starting type without telling me the N value!") unless defined $N;
      my $type = $1;
      $type =~ tr/a-z/A-Z/;
      error("Unknown type: $type") unless defined $proof_funcs{$type};
      my ($n, @q) = $proof_funcs{$type}->();
      $parts{$n} = [@q];
      if    ($verbose == 1) {  print ".";  }
      elsif ($verbose == 2) {  printf "step %2d: %4d digits  type %-12s\r", $step++, length($n), $type; }
    }
  }
}
error("No N found") unless defined $N;
print " " x 60, "\r" if $verbose == 2;
if (final_verify($N)) {
  print "PRIME\n" if $verbose;
  exit(0);
} else {
  print "NOT PROVEN\n" if $verbose;
  exit(2);
}

sub final_verify {
  my $n = shift;
  die "Internal error: argument not defined" unless defined $n;

  if ($timing) {
    my $seconds = tv_interval($start_time);
    printf "%7.6f seconds for verification of %d digit number\n", $seconds, length($orig_N);
  }

  if ($cert_type eq 'PRIMO') {
    fail "Type 0: $n failed N > 18" unless $n > 18;
    fail "Type 0: $n failed N < 34 * 10^13" unless $n < (34*10**13);
    fail "Type 0: $n failed spsp(2,3,5,7,11,13,17)" unless is_strong_pseudoprime($n,2,3,5,7,11,13,17);
    return 1;
  }

  my @qs = ($n);
  while (@qs) {
    my $q = shift @qs;
    # Check that this q has a chain
    if (!defined $parts{$q}) {
      # Auto-small: handle small q right here.
      if ($q <= $smallval) {
        fail "Small n $q does not pass BPSW" unless is_prime($q);
        next;
      } else {
        error "q value $q has no proof\n";
      }
    }
    die "Internal error: Invalid parts entry" unless ref($parts{$q}) eq 'ARRAY';
    # q is prime if all it's chains are prime.
    push @qs, @{$parts{$q}};
  }
  1;
}


##############################################################################
# MPU Proof handlers
##############################################################################

sub prove_ecpp {
  verify_ecpp( read_vars('ECPP', qw/N A B M Q X Y/) );
}
sub prove_ecpp3 {
  verify_ecpp3( read_vars('ECPP3', qw/N S R A B T/) );
}
sub prove_ecpp4 {
  verify_ecpp4( read_vars('ECPP4', qw/N S R J T/) );
}
sub prove_bls15 {
  verify_bls15( read_vars('BLS15', qw/N Q LP LQ/) );
}
sub prove_bls3 {
  verify_bls3( read_vars('BLS3', qw/N Q A/) );
}
sub prove_pock {
  verify_pock( read_vars('POCKLINGTON', qw/N Q A/) );
}
sub prove_small {
  verify_small( read_vars('Small', qw/N/) );
}
sub prove_bls5 {
  # No good way to do this using read_vars
  my ($n, @Q, @A);
  my $index = 0;
  $Q[0] = Math::BigInt->new(2);  # 2 is implicit
  while (1) {
    my $line = <>;
    error("end of file during type BLS5") unless defined $line;
    # Skip comments and blank lines
    next if $line =~ /^\s*#/ or $line =~ /^\s*$/;
    # Stop when we see a line starting with -.
    last if $line =~ /^-/;
    chomp($line);
    if ($line =~ /^N\s+(\d+)/) {
      error("BLS5: N redefined") if defined $n;
      $n = Math::BigInt->new("$1");
    } elsif ($line =~ /^Q\[(\d+)\]\s+(\d+)/) {
      $index++;
      error("BLS5: Invalid index: $1") unless $1 == $index;
      $Q[$1] = Math::BigInt->new("$2");
    } elsif ($line =~ /^A\[(\d+)\]\s+(\d+)/) {
      error("BLS5: Invalid index: A[$1]") unless $1 >= 0 && $1 <= $index;
      $A[$1] = Math::BigInt->new("$2");
    } else {
      error("Unrecognized line: $line");
    }
  }
  verify_bls5($n, \@Q, \@A);
}

sub prove_lucas {
  # No good way to do this using read_vars
  my ($n, @Q, $a);
  my $index = 0;
  while (1) {
    my $line = <>;
    error("end of file during type Lucas") unless defined $line;
    # Skip comments and blank lines
    next if $line =~ /^\s*#/ or $line =~ /^\s*$/;
    chomp($line);
    if ($line =~ /^N\s+(\d+)/) {
      error("Lucas: N redefined") if defined $n;
      $n = Math::BigInt->new("$1");
    } elsif ($line =~ /^Q\[(\d+)\]\s+(\d+)/) {
      $index++;
      error("Lucas: Invalid index: $1") unless $1 == $index;
      $Q[$1] = Math::BigInt->new("$2");
    } elsif ($line =~ /^A\s+(\d+)/) {
      $a = Math::BigInt->new("$1");
      last;
    } else {
      error("Unrecognized line: $line");
    }
  }
  verify_lucas($n, \@Q, $a);
}

##############################################################################
# Proof verifications
##############################################################################

sub verify_ecpp {
  my ($n, $a, $b, $m, $q, $x, $y) = @_;
  $a %= $n if $a < 0;
  $b %= $n if $b < 0;
  fail "ECPP: $n failed N > 0" unless $n > 0;
  fail "ECPP: $n failed gcd(N, 6) = 1" unless Math::BigInt::bgcd($n, 6) == 1;
  fail "ECPP: $n failed gcd(4*a^3 + 27*b^2, N) = 1"
    unless Math::BigInt::bgcd(4*$a*$a*$a+27*$b*$b,$n) == 1;
  fail "ECPP: $n failed Y^2 = X^3 + A*X + B mod N"
    unless ($y*$y) % $n == ($x*$x*$x + $a*$x + $b) % $n;
  fail "ECPP: $n failed M >= N - 2*sqrt(N) + 1" unless $m >= $n - 2*$n->copy->bsqrt() + 1;
  fail "ECPP: $n failed M <= N + 2*sqrt(N) + 1" unless $m <= $n + 2*$n->copy->bsqrt() + 1;
  fail "ECPP: $n failed Q > (N^(1/4)+1)^2" unless $q > $n->copy->broot(4)->badd(1)->bpow(2);
  fail "ECPP: $n failed Q < N" unless $q < $n;
  fail "ECPP: $n failed M != Q" unless $m != $q;
  my ($mdivq, $rem) = $m->copy->bdiv($q);
  fail "ECPP: $n failed Q divides M" unless $rem == 0;

  # Now verify the elliptic curve
  my $correct_point = 0;
  if (prime_get_config->{'gmp'} && defined &Math::Prime::Util::GMP::_validate_ecpp_curve) {
    $correct_point = Math::Prime::Util::GMP::_validate_ecpp_curve($a, $b, $n, $x, $y, $m, $q);
  } else {
    if (!defined $Math::Prime::Util::ECAffinePoint::VERSION) {
      eval { require Math::Prime::Util::ECAffinePoint; 1; }
      or do { die "Cannot load Math::Prime::Util::ECAffinePoint"; };
    }
    my $ECP = Math::Prime::Util::ECAffinePoint->new($a, $b, $n, $x, $y);
    # Compute U = (m/q)P, check U != point at infinity
    $ECP->mul( $m->copy->bdiv($q)->as_int );
    if (!$ECP->is_infinity) {
      # Compute V = qU, check V = point at infinity
      $ECP->mul( $q );
      $correct_point = 1 if $ECP->is_infinity;
    }
  }
  fail "ECPP: $n failed elliptic curve conditions" unless $correct_point;
  ($n, $q);
}
sub verify_ecpp3 {
  my ($n, $s, $r, $a, $b, $t) = @_;
  fail "ECPP3: $n failed |A| <= N/2" unless 2*abs($a) <= $n;
  fail "ECPP3: $n failed |B| <= N/2" unless 2*abs($b) <= $n;
  fail "ECPP3: $n failed T >= 0" unless $t >= 0;
  fail "ECPP3: $n failed T < N" unless $t < $n;
  my $l = ($t*$t*$t + $a*$t + $b) % $n;
  verify_ecpp( $n,
               ($a * $l*$l) % $n,
               ($b * $l*$l*$l) % $n,
               $r*$s,
               $r,
               ($t*$l) % $n,
               ($l*$l) % $n    );
}
sub verify_ecpp4 {
  my ($n, $s, $r, $j, $t) = @_;
  fail "ECPP4: $n failed |J| <= N/2" unless 2*abs($j) <= $n;
  fail "ECPP4: $n failed T >= 0" unless $t >= 0;
  fail "ECPP4: $n failed T < N" unless $t < $n;
  my $a = 3 * $j * (1728 - $j);
  my $b = 2 * $j * (1728 - $j) * (1728 - $j);
  my $l = ($t*$t*$t + $a*$t + $b) % $n;
  verify_ecpp( $n,
               ($a * $l*$l) % $n,
               ($b * $l*$l*$l) % $n,
               $r*$s,
               $r,
               ($t*$l) % $n,
               ($l*$l) % $n    );
}

sub verify_bls15 {
  my ($n, $q, $lp, $lq) = @_;
  fail "BLS15: $n failed Q odd" unless $q->is_odd();
  fail "BLS15: $n failed Q > 2" unless $q > 2;
  my ($m, $rem) = ($n+1)->copy->bdiv($q);
  fail "BLS15: $n failed Q divides N+1" unless $rem == 0;
  fail "BLS15: $n failed MQ-1 = N" unless $m*$q-1 == $n;
  fail "BLS15: $n failed M > 0" unless $m > 0;
  fail "BLS15: $n failed 2Q-1 > sqrt(N)" unless 2*$q-1 > $n->copy->bsqrt();
  my $D = $lp*$lp - 4*$lq;
  fail "BLS15: $n failed D != 0" unless $D != 0;
  fail "BLS15: $n failed jacobi(D,N) = -1" unless _jacobi($D,$n) == -1;
  fail "BLS15: $n failed V_{m/2} mod N != 0"
      unless (lucas_sequence($n, $lp, $lq, $m/2))[1] != 0;
  fail "BLS15: $n failed V_{(N+1)/2} mod N == 0"
      unless (lucas_sequence($n, $lp, $lq, ($n+1)/2))[1] == 0;
  ($n, $q);
}

sub verify_bls3 {
  my ($n, $q, $a) = @_;
  fail "BLS3: $n failed Q odd" unless $q->is_odd();
  fail "BLS3: $n failed Q > 2" unless $q > 2;
  my ($m, $rem) = ($n-1)->copy->bdiv($q);
  fail "BLS3: $n failed Q divides N-1" unless $rem == 0;
  fail "BLS3: $n failed MQ+1 = N" unless $m*$q+1 == $n;
  fail "BLS3: $n failed M > 0" unless $m > 0;
  fail "BLS3: $n failed 2Q+1 > sqrt(n)" unless 2*$q+1 > $n->copy->bsqrt();
  fail "BLS3: $n failed A^((N-1)/2) = N-1 mod N" unless $a->copy->bmodpow(($n-1)/2, $n) == $n-1;
  fail "BLS3: $n failed A^(M/2) != N-1 mod N" unless $a->copy->bmodpow($m/2,$n) != $n-1;
  ($n, $q);
}
sub verify_pock {
  my ($n, $q, $a) = @_;
  my ($m, $rem) = ($n-1)->copy->bdiv($q);
  fail "Pocklington: $n failed Q divides N-1" unless $rem == 0;
  fail "Pocklington: $n failed M is even" unless $m->is_even();
  fail "Pocklington: $n failed M > 0" unless $m > 0;
  fail "Pocklington: $n failed M < Q" unless $m < $q;
  fail "Pocklington: $n failed MQ+1 = N" unless $m*$q+1 == $n;
  fail "Pocklington: $n failed A > 1" unless $a > 1;
  fail "Pocklington: $n failed A^(N-1) mod N = 1"
    unless $a->copy->bmodpow($n-1, $n) == 1;
  fail "Pocklington: $n failed gcd(A^M - 1, N) = 1"
    unless Math::BigInt::bgcd($a->copy->bmodpow($m, $n)-1, $n) == 1;
  ($n, $q);
}
sub verify_small {
  my ($n) = @_;
  fail "Small n $n is > 2^64\n" unless $n <= $smallval;
  fail "Small n $n does not pass BPSW" unless is_prime($n);
  ($n);
}

sub verify_bls5 {
  my ($n, $Qr, $Ar) = @_;
  my @Q = @{$Qr};
  my @A = @{$Ar};
  my $nm1 = $n - 1;
  my $F = Math::BigInt->bone;
  my $R = $nm1->copy;
  my $index = $#Q;
  foreach my $i (0 .. $index) {
    error "BLS5: $n failed Q[$i] doesn't exist" unless defined $Q[$i];
    $A[$i] = Math::BigInt->new(2) unless defined $A[$i];
    fail "BLS5: $n failed Q[$i] > 1" unless $Q[$i] > 1;
    fail "BLS5: $n failed Q[$i] < N-1" unless $Q[$i] < $nm1;
    fail "BLS5: $n failed A[$i] > 1" unless $A[$i] > 1;
    fail "BLS5: $n failed A[$i] < N" unless $A[$i] < $n;
    fail "BLS5: $n failed Q[$i] divides N-1" unless ($nm1 % $Q[$i]) == 0;
    while (($R % $Q[$i]) == 0) {
      $F *= $Q[$i];
      $R /= $Q[$i];
    }
  }
  die "BLS5: Internal error R != (N-1)/F\n" unless $R == $nm1/$F;
  fail "BLS5: $n failed F is even" unless $F->is_even();
  fail "BLS5: $n failed gcd(F, R) = 1\n" unless Math::BigInt::bgcd($F,$R) == 1;
  my ($s, $r) = $R->copy->bdiv(2*$F);
  my $P = ($F+1) * (2 * $F * $F + ($r-1)*$F + 1);
  fail "BLS5: $n failed n < P" unless $n < $P;
  fail "BLS5: $n failed s=0 OR r^2-8s not a perfect square"
    unless $s == 0 or !_is_perfect_square($r*$r - 8*$s);
  foreach my $i (0 .. $index) {
    my $a = $A[$i];
    my $q = $Q[$i];
    fail "BLS5: $n failed A[i]^(N-1) mod N = 1"
      unless $a->copy->bmodpow($nm1, $n) == 1;
    fail "BLS5: $n failed gcd(A[i]^((N-1)/Q[i])-1, N) = 1"
      unless Math::BigInt::bgcd($a->copy->bmodpow($nm1/$q, $n)-1, $n) == 1;
  }
  ($n, @Q);
}

sub verify_lucas {
  my ($n, $Qr, $a) = @_;
  my @Q = @{$Qr};
  my $index = $#Q;
  fail "Lucas: $n failed A > 1" unless $a > 1;
  fail "Lucas: $n failed A < N" unless $a < $n;
  my $nm1 = $n - 1;
  my $F = Math::BigInt->bone;
  my $R = $nm1->copy;
  fail "Lucas: $n failed A^(N-1) mod N = 1"
    unless $a->copy->bmodpow($nm1, $n) == 1;
  foreach my $i (1 .. $index) {
    error "Lucas: $n failed Q[$i] doesn't exist" unless defined $Q[$i];
    fail "Lucas: $n failed Q[$i] > 1" unless $Q[$i] > 1;
    fail "Lucas: $n failed Q[$i] < N-1" unless $Q[$i] < $nm1;
    fail "Lucas: $n failed Q[$i] divides N-1" unless ($nm1 % $Q[$i]) == 0;
    fail "Lucas: $n failed A^((N-1)/Q[$i]) mod N != 1"
      unless $a->copy->bmodpow($nm1/$Q[$i], $n) != 1;
    while (($R % $Q[$i]) == 0) {
      $F *= $Q[$i];
      $R /= $Q[$i];
    }
  }
  fail("Lucas: $n failed N-1 has only factors Q") unless $R == 1 && $F == $nm1;
  shift @Q;  # Remove Q[0]
  ($n, @Q);
}


##############################################################################
# Utility functions
##############################################################################


sub read_vars {
  my $type = shift;
  my %vars = map { $_ => 1 } @_;
  my %return;
  while (scalar keys %vars) {
    my $line = <>;
    error("end of file during type $type") unless defined $line;
    # Skip comments and blank lines
    next if $line =~ /^\s*#/ or $line =~ /^\s*$/;
    chomp($line);
    error("Still missing values in type $type") if $line =~ /^Type /;
    if ($line =~ /^(\S+)\s+(-?\d+)/) {
      my ($var, $val) = ($1, $2);
      $var =~ tr/a-z/A-Z/;
      error("Type $type: repeated or inappropriate var: $line") unless defined $vars{$var};
      $return{$var} = $val;
      delete $vars{$var};
    } else {
      error("Unrecognized line: $line");
    }
  }
  # Now return them in the order given, turned into bigints.
  return map { Math::BigInt->new("$return{$_}") } @_;
}

sub primo_read_vars {
  my $type = shift;
  my %vars = map { $_ => 1 } @_;
  my %return;
  while (scalar keys %vars) {
    my $line = <>;
    error("end of file during type $type") unless defined $line;
    error("blank line during type $type") if $line =~ /^\s*$/;
    chomp($line);
    error("Still missing values in type $type") if $line =~ /^Type=/;
    if ($line =~ /^(\S+)\s*=\s*(\S+)/) {
      my ($var, $val) = ($1, $2);
      $var =~ tr/a-z/A-Z/;
      $val = "0x$val" if $var =~ s/\$$//;
      # For Primo, just skip things we don't understand.
      next unless defined $vars{$var};
      $return{$var} = $val;
      delete $vars{$var};
    } else {
      error("Unrecognized line: $line");
    }
  }
  # Now return them in the order given, turned into bigints.
  my @ret;
  foreach my $var (@_) {
    my $sign = 1;
    $sign = -1 if $return{$var} =~ s/^(0x)?-/$1/;
    push @ret, Math::BigInt->new($return{$var}) * $sign;
  }
  return @ret;
}




sub _is_perfect_square {
  my($n) = @_;

  if (ref($n) eq 'Math::BigInt') {
    my $mc = int(($n & 31)->bstr);
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

# Calculate Jacobi symbol (M|N)
sub _jacobi {
  my($n, $m) = @_;
  return 0 if $m <= 0 || ($m % 2) == 0;
  my $j = 1;
  if ($n < 0) {
    $n = -$n;
    $j = -$j if ($m % 4) == 3;
  }
  # Split loop so we can reduce n/m to non-bigints after first iteration.
  if ($n != 0) {
    while (($n % 2) == 0) {
      $n >>= 1;
      $j = -$j if ($m % 8) == 3 || ($m % 8) == 5;
    }
    ($n, $m) = ($m, $n);
    $j = -$j if ($n % 4) == 3 && ($m % 4) == 3;
    $n = $n % $m;
    $n = int($n->bstr) if ref($n) eq 'Math::BigInt' && $n <= ''.~0;
    $m = int($m->bstr) if ref($m) eq 'Math::BigInt' && $m <= ''.~0;
  }
  while ($n != 0) {
    while (($n % 2) == 0) {
      $n >>= 1;
      $j = -$j if ($m % 8) == 3 || ($m % 8) == 5;
    }
    ($n, $m) = ($m, $n);
    $j = -$j if ($n % 4) == 3 && ($m % 4) == 3;
    $n = $n % $m;
  }
  return ($m == 1) ? $j : 0;
}
