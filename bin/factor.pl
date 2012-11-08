#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use bigint try => 'GMP';
use Math::Prime::Util qw/factor/;
$| = 1;
no bigint;

my %opts;
GetOptions(\%opts,
           'version',   # turn off MPU::GMP for debugging
           'help',
          ) || die_usage();
if (exists $opts{'version'}) {
  my $version_str =
   "factor.pl version 1.0 using Math::Prime::Util $Math::Prime::Util::VERSION";
  $version_str .= " and MPU::GMP $Math::Prime::Util::GMP::VERSION"
    if Math::Prime::Util::prime_get_config->{'gmp'};
  $version_str .= "\nWritten by Dana Jacobsen.\n";
  die "$version_str";
}
die_usage() if exists $opts{'help'};

if (@ARGV) {
  foreach my $n (@ARGV) {
    print "$n: ", join(" ", factor($n)), "\n";
  }
} else {
  while (<>) {
    chomp;
    foreach my $n (split / /) {
      print "$n: ", join(" ", factor($n)), "\n";
    }
  }
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
