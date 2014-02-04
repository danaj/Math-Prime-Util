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

  # Based on invphi.gp v1.3 by Max Alekseyev
  my @L;
  foreach my $n (divisors($N)) {
    $n = $do_bigint->new("$n") if $do_bigint;
    my $p = $n+1;
    next unless is_prime($p);
    if ( ($N % $p) != 0 ) {
      push @L, [ [$n, $p] ];
    } else {
      my ($v, $t) = (2, int($N/$p));
      while (!($t % $p)) { $v++; $t = int($t/$p); }
      push @L, [ [$n,$p], map { [$n*$p**($_-1), $p**$_] } 2..$v ];
    }
  }

  if (!wantarray) {   # Just count.  Much less memory.
    my %r = ( 1 => 1 );
    my($l0, $l1);
    foreach my $Li (@L) {
      my %t;
      foreach my $Lij (@$Li) {
        my($l0, $l1) = @$Lij;
        foreach my $n (divisors($N / $l0)) {
          $t{$n*$l0} += $r{$n} if defined $r{$n};
        }
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
