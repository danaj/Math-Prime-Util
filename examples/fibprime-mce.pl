#!/usr/bin/env perl
use strict;
use warnings;

# Overkill, but let's try to select a good bigint module.
my $bigint_class;
if      (eval { require Math::GMPz; 1; }) {
  $bigint_class = "Math::GMPz";
} elsif (eval { require Math::GMP; 1; }) {
  $bigint_class = "Math::GMP";
} else {
  require Math::BigInt;
  Math::BigInt->import(try=>"GMP,Pari");
  $bigint_class = "Math::BigInt";
}

use Math::Prime::Util ':all';
use Time::HiRes qw(gettimeofday tv_interval);
use MCE::Util qw(get_ncpu);
use MCE;
$| = 1;

# Find Fibonacci primes in parallel, using Math::Prime::Util and MCE.
#
# Dana Jacobsen, 2012.
# Mario Roy, 2014.
#
# Runs about the same speed as the threads version, but doesn't need
# a threaded Perl.
#
# n32 ( F50833) in  4776s on 3930k 4.2GHz, SERIAL
# n32 ( F50833) in   754s on 3930k 4.2GHz, 12 CPU
# n32 ( F50833) in   472s on EC2 c3.8xlarge, 32 CPU
# n32 ( F50833) in   323s on EC2 c4.8xlarge, 36 CPU
#
# n36 (F148091) in 26245s on 3930k 4.2GHz, 12 CPU
# n36 (F148091) in 14380s on EC2 c3.8xlarge, 32 CPU
#

my $time_start = [gettimeofday];
my $nworkers = get_ncpu();
warn "Using $nworkers CPUs\n";
prime_precalc(10_000_000);

sub fib_n {
   my ($n, $fibstate) = @_;
   @$fibstate = (1, $bigint_class->new(0), $bigint_class->new(1))
      unless defined $fibstate->[0];
   my ($curn, $a, $b) = @$fibstate;
   die "fib_n only increases" if $n < $curn;
   do { ($a, $b) = ($b, $a+$b); } for (1 .. $n-$curn);
   @$fibstate = ($n, $a, $b);
   $b;
}

sub nth_iter {
   my $n = 0; my $order_id = 1; my %tmp;
   return sub {
      $tmp{$_[0]} = $_[1];   ## @_ = ( $nth, [ $k, $time_int ] )
      while (1) {
         last if not exists $tmp{$order_id};
         if (defined $tmp{$order_id}) {
            my ($k, $time_int) = @{ $tmp{$order_id} };
            printf "%3d %7d %20.5f\n", ++$n, $k, $time_int;
         }
         delete $tmp{$order_id++};
      }
   }
}

my $mce = MCE->new(
   max_workers => $nworkers, gather => nth_iter,

   user_func => sub {
      my @fibstate; my $nth = MCE->wid();

      while (1) {
         # Exploit knowledge that excepting k=4, all prime F_k have a prime k.
         my $k = ($nth <= 2) ?  2 + $nth  :  nth_prime($nth);
         my $Fk = fib_n($k, \@fibstate);
         if (is_prob_prime($Fk)) {
            MCE->gather($nth, [ $k, tv_interval($time_start) ]);
         } else {
            MCE->gather($nth, undef);
         }
         $nth += $nworkers;
      }
   }

)->run;

