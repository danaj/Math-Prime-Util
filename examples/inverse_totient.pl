use warnings;
use strict;
use Math::Prime::Util qw/:all/;
use Getopt::Long;

my %opts;
GetOptions(\%opts,
           'count',
           'help',
          ) || die_usage();
die_usage() if exists $opts{'help'};
my $n = shift;
die_usage() unless defined $n && length($n) > 0 && $n !~ tr/0123456789//c;
if (exists $opts{'count'}) {
  print scalar inverse_euler_phi($n), "\n";
} else {
  print join("\n", inverse_euler_phi($n)), "\n";
}

sub die_usage {
  die "Usage: $0 [-count] <m>\n\nPrint all n such that euler_phi(n) = m.\nIf -count is used, just prints the number of such n.\n";
}


sub inverse_euler_phi {
  my $N = shift;

  my $do_bigint = ($N > 2**49);
  if ($do_bigint) {
    # Math::GMPz and Math::GMP are fast.  Math::BigInt::GMP is 10x slower.
    eval { use Math::GMPz; $do_bigint = "Math::GMPz"; 1; } ||
    eval { use Math::GMP; $do_bigint = "Math::GMP"; 1; } ||
    eval { use Math::BigInt try=>"GMP,Pari"; $do_bigint = "Math::BigInt"; 1; };
    $N = $do_bigint->new("$N");
  }

  return wantarray ? (1,2) : 2 if $N == 1;
  return wantarray ? () : 0 if $N < 1 || ($N & 1);
  if (is_prime($N >> 1)) {   # Coleman Remark 3.3 (Thm 3.1) and Prop 6.2
    return wantarray ? () : 0             if !is_prime($N+1);
    return wantarray ? ($N+1, 2*$N+2) : 2 if $N >= 10;
  }

  #if (!wantarray) { return a014197($N) }

  # Based on invphi.gp v1.3 by Max Alekseyev
  my @L;
  fordivisors { $n=$_;
    $n = $do_bigint->new("$n") if $do_bigint;
    my $p = $n+1;
    if (is_prime($p)) {
      if ( ($N % $p) != 0 ) {
        push @L, [ [$n, $p] ];
      } else {
        my $v = valuation($N, $p);
        my $t = $N / $p**$v;
        push @L, [ [$n,$p], map { [$n*$p**($_-1), $p**$_] } 2..$v+1 ];
      }
    }
  } $N;

  if (!wantarray) {   # Just count.  Much less memory.
    my %r = ( 1 => 1 );
    my($l0, $l1);
    foreach my $Li (@L) {
      my %t;
      foreach my $Lij (@$Li) {
        my($l0, $l1) = @$Lij;
        fordivisors {
          $t{$_*$l0} += $r{$_} if defined $r{$_};
        } $N / $l0;
      }
      while (my($i,$vec) = each(%t)) { $r{$i} += $t{$i}; }
    }
    return (defined $r{$N}) ? $r{$N} : 0;
  }

  my %r = ( 1 => [1] );
  my($l0, $l1);
  foreach my $Li (@L) {
    my %t;
    foreach my $Lij (@$Li) {
      my($l0, $l1) = @$Lij;
      foreach my $n (divisors($N / $l0)) {
        push @{ $t{$n*$l0} }, map { $_ * $l1 } @{ $r{$n} }
          if defined $r{$n};
      }
    }
    while (my($i,$vec) = each(%t)) { push @{$r{$i}}, @$vec; }
  }
  return () unless defined $r{$N};
  delete @r{ grep { $_ != $N } keys %r };  # Delete all intermediate results
  my @result = sort { $a <=> $b } @{$r{$N}};
  return @result;
}

# Simple recursive count translated from Pari.
sub a014197 {
  my($n,$m) = @_;
  $m=1 unless defined $m;
  return 1+($m<2) if $n == 1;
  # TODO: divisor_sum with sub ought to be faster
  #divisor_sum( $n, sub { my $d=shift;
  #  return 0 if $d < $m || !is_prime($d+1);
  #  my($p, $q) = ($d+1, $n/$d);
  #  vecsum( map { a014197($q/($p**$_), $p) } 0 .. valuation($q,$p) );
  #} );
  my($sum,$p,$q) = (0);
  fordivisors {
    if ($_ >= $m && is_prime($_+1)) {
      ($p,$q)=($_+1,$n/$_);
      $sum += vecsum( map { a014197($q/($p**$_), $p) } 0 .. valuation($q,$p) );
    }
  } $n;
  $sum;
}
