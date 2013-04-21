#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use bigint try => 'GMP';
use Math::Prime::Util qw/factor nth_prime/;
$| = 1;
no bigint;

# Allow execution of any of these functions in the command line
my @mpu_funcs = (qw/next_prime prev_prime prime_count nth_prime random_prime
                    random_ndigit_prime random_nbit_prime random_strong_prime
                    random_maurer_prime primorial pn_primorial moebius mertens
                    euler_phi jordan_totient exp_mangoldt divisor_sum
                    consecutive_integer_lcm/);
my %mpu_func_map;

my %opts;
GetOptions(\%opts,
           'version',   # turn off MPU::GMP for debugging
           'help',
          ) || die_usage();
if (exists $opts{'version'}) {
  my $version_str =
   "factor.pl version 1.1 using Math::Prime::Util $Math::Prime::Util::VERSION";
  $version_str .= " and MPU::GMP $Math::Prime::Util::GMP::VERSION"
    if Math::Prime::Util::prime_get_config->{'gmp'};
  $version_str .= "\nWritten by Dana Jacobsen.\n";
  die "$version_str";
}
die_usage() if exists $opts{'help'};

if (@ARGV) {
  foreach my $n (@ARGV) {
    $n = eval_expr($n) unless $n =~ /^\d+$/;
    print "$n: ", join(" ", factor($n)), "\n";
  }
} else {
  while (<>) {
    chomp;
    foreach my $n (split / /) {
      $n = eval_expr($n) unless $n =~ /^\d+$/;
      print "$n: ", join(" ", factor($n)), "\n";
    }
  }
}

# This is rather braindead.  We're going to eval their input so they can give
# arbitrary expressions.  But we only want to allow math-like strings.
sub eval_expr {
  my $expr = shift;
  die "$expr cannot be evaluated" if $expr =~ /:/;  # Use : for escape
  if (scalar(keys %mpu_func_map) == 0) {
    my $n = 10;
    foreach my $func (@mpu_funcs) {
      $mpu_func_map{$func} = sprintf("%03d", $n++);
    }
  }
  $expr =~ s/\blog\(/:001(/g;
  foreach my $func (@mpu_funcs) {
    $expr =~ s/\b$func\(/:$mpu_func_map{$func}(/g;
  }
  die "$expr cannot be evaluated" if $expr =~ tr|-0123456789+*/() :||c;
  $expr =~ s/:001/log/g;
  foreach my $func (@mpu_funcs) {
    $expr =~ s/:$mpu_func_map{$func}\(/Math::Prime::Util::$func(/g;
  }
  $expr =~ s/(\d+)/ Math::BigInt->new("$1") /g;
  my $res = eval $expr; ## no critic
  die "Cannot eval: $expr\n" if !defined $res;
  $res = int($res->bstr) if ref($res) eq 'Math::BigInt' && $res <= ~0;
  $res;
}


sub die_usage {
  die <<EOU;
Usage: $0 [options] [number] ...

Print the prime factors of each positive integer given on the command line,
or reads numbers from standard input if called without arguments.

Math expressions may be given as arguments, which will be evaluated before
factoring.  This includes most Math::Prime::Util functions including things
like prime_count(#), nth_prime(#), primorial(#), random_nbit_prime(#), etc.

  --help       displays this help message
  --version    displays the version information

Part of the Math::Prime::Util $Math::Prime::Util::VERSION package, wrapping
the factor() function.  See 'man Math::Prime::Util' for more information.
EOU
}
