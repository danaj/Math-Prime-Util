#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util ":all";
$| = 1;  # fast pipes

my %nthprimes = (
                  1 => 2,
                 10 => 29,
                100 => 541,
               1000 => 7919,
              10000 => 104729,
             100000 => 1299709,
            1000000 => 15485863,
           10000000 => 179424673,
          100000000 => 2038074743,
         1000000000 => 22801763489,
        10000000000 => 252097800623,
       100000000000 => 2760727302517,
      1000000000000 => 29996224275833,
     10000000000000 => 323780508946331,
    100000000000000 => 3475385758524527,
   1000000000000000 => 37124508045065437,
  10000000000000000 => 394906913903735329,
 100000000000000000 => 4185296581467695669,
);

printf("  N    %12s  %12s\n", "nth_approx", "percent");
printf("-----  %12s  %12s\n", '-'x12, '-'x12);
foreach my $n (sort {$a<=>$b} keys %nthprimes) {
  my $nth  = $nthprimes{$n};
  my $ntha = nth_prime_approx($n);

  printf "10^%2d %13lu  %12.7f\n", length($n)-1, abs($nth-$ntha), 100*($ntha-$nth)/$nth;
}

print "\n";
print "Lower / Upper bounds.  Percentages.\n";
print "\n";

printf("  N    %12s  %12s\n", "lower", "upper");
printf("-----  %12s  %12s\n", '-'x12,'-'x12);
foreach my $n (sort {$a<=>$b} keys %nthprimes) {
  my $nth  = $nthprimes{$n};
  my $nthl  = nth_prime_lower($n);
  my $nthu  = nth_prime_upper($n);

  printf "10^%2d  %12.7f  %12.7f\n",
         length($n)-1, 100.0*($nth-$nthl)/$nth, 100.0*($nthu-$nth)/$nth;
}
