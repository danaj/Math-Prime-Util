#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/:all/;
use Time::HiRes qw(gettimeofday tv_interval);
use feature 'say';
use Math::GMPz;
use Math::GMP;
use Math::BigInt;
$| = 1;  # fast pipes

my $class = 'Math::GMPz';
prime_set_config(bigint=>$class);

my($r,$x);

$r = is_semiprime(165);
die "is_semiprime $r" unless $r == 0 && !ref($r);
print ".";
$r = is_semiprime(166);
die "is_semiprime $r" unless $r == 1 && !ref($r);
print ".";
$r = is_semiprime($class->new("27273137616939507011"));
die "is_semiprime $r" unless $r == 1 && !ref($r);
print ".";

$r = random_semiprime(32);
die "random_semiprime 32" unless is_semiprime($r) && !ref($r);
print ".";
$r = random_semiprime(65);
die "random_semiprime 65" unless is_semiprime($r) && ref($r) eq $class;
print ".";

$r = random_unrestricted_semiprime(32);
die "random_unrestricted_semiprime 32" unless is_semiprime($r) && !ref($r);
print ".";
$r = random_unrestricted_semiprime(65);
die "random_unrestricted_semiprime 65" unless is_semiprime($r) && ref($r) eq $class;
print ".";

$r = semi_primes(1000000000,1000000010);
die "semiprimes small" unless $r->[0] == 1000000006 && !ref($r->[0]);
print ".";

$x = powint(2,70);
die "powint" unless ref($x) eq $class;
print ".";

$r = semi_primes($x,$x+10);
die "semiprimes large" unless $r->[0] == $x+9 && ref($r->[0]) eq $class;
print ".";

$r = semiprime_count(1000000000);
die "semiprime_count" unless $r == 160788536 && !ref($r);
print ".";

$r = semiprime_count_approx($x);
die "semiprime_count_approx" unless ref($r) eq $class;
print ".";

$r = nth_semiprime(1000000);
die "nth_semiprime" unless $r == 5109839 && !ref($r);
print ".";

$r = nth_semiprime_approx($x);
die "nth_semiprime_approx" unless ref($r) eq $class;
print ".";

my @S;
forsemiprimes { push @S, $_; } 3000000000,3000000005;
die "forsemiprimes small" unless vecequal(\@S,[3000000001,3000000002,3000000005]) && vecall { !ref($_) } @S;
print ".";

@S=();
forsemiprimes { push @S, $_; } $x,$x+10;
die "forsemiprimes large" unless vecequal(\@S,["1180591620717411303433"]) && vecall { ref($_) eq $class } @S;
print ".";

print " PASS\n";
