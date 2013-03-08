#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use bigint try => 'GMP';
use Math::Prime::Util qw/factor nth_prime/;
$| = 1;
no bigint;

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
  $expr =~ s/nth_prime\(/:1(/g;
  $expr =~ s/log\(/:2(/g;
  die "$expr cannot be evaluated" if $expr =~ tr|-0123456789+*/() :||c;
  $expr =~ s/:1/nth_prime/g;
  $expr =~ s/:2/log/g;
  $expr =~ s/(\d+)/ Math::BigInt->new($1) /g;
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

  --help       displays this help message
  --version    displays the version information

Part of the Math::Prime::Util $Math::Prime::Util::VERSION package, wrapping
the factor() function.  See 'man Math::Prime::Util' for more information.
EOU
}
