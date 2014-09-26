package ntheory;
use strict;
use warnings;

BEGIN {
  $ntheory::AUTHORITY = 'cpan:DANAJ';
  $ntheory::VERSION = '0.45';
}

BEGIN {
  require Math::Prime::Util; 
  *ntheory:: = *Math::Prime::Util::;
}

1;

__END__


# ABSTRACT: Number theory utilities

=pod

=encoding utf8

=for stopwords ntheory

=head1 NAME

ntheory - Number theory utilities

=head1 SEE

See L<Math::Prime::Util> for complete documentation.

=head1 QUICK REFERENCE

=head2 PRIMALITY

  is_prob_prime(n)                    primality test (BPSW)
  is_prime(n)                         primality test (BPSW + extra)
  is_provable_prime(n)                primality test with proof
  is_provable_prime_with_cert(n)      primality test: (isprime,cert)
  prime_certificate(n)                as above with just certificate
  verify_prime(cert)                  verify a primality certificate
  is_aks_prime                        AKS deterministic test (slow)

=head2 PROBABLE PRIME TESTS

  is_pseudoprime(n,bases)                  Fermat probable prime tests
  is_strong_pseudoprime(n,bases)           Miller-Rabin tests to bases
  is_lucas_pseudoprime(n)                  Lucas test
  is_strong_lucas_pseudoprime(n)           strong Lucas test
  is_almost_extra_strong_lucas_pseudoprime(n, [incr])   AES Lucas test
  is_extra_strong_lucas_pseudoprime(n)     extra strong Lucas test
  is_frobenius_pseudoprime(n, [a,b])       Frobenius quadratic test
  is_perrin_pseudoprime(n)                 Perrin test
  is_frobenius_underwood_pseudoprime(n)    combined PSP and Lucas
  is_bpsw_prime(n)                         combined SPSP-2 and ES Lucas
  miller_rabin_random(n, ntests)           perform random-base MR tests

=head2 PRIMES

  primes([start,] end)                array ref of primes
  twin_primes([start,] end)           array ref of twin primes
  next_prime(n)                       next prime > n
  prev_prime(n)                       previous prime < n
  prime_count(n)                      count of primes <= n
  prime_count(start, end)             count of primes in range
  prime_count_lower(n)                fast lower bound for prime count
  prime_count_upper(n)                fast upper bound for prime count
  prime_count_approx(n)               fast approximate count of primes
  nth_prime(n)                        the nth prime (n=1 returns 2)
  nth_prime_lower(n)                  fast lower bound for nth prime
  nth_prime_upper(n)                  fast upper bound for nth prime
  nth_prime_approx(n)                 fast approximate nth prime
  twin_prime_count(n)                 count of twin primes <= n
  twin_prime_count(start, end)        count of twin primes in range
  twin_prime_count_approx(n)          fast approx count of twin primes
  nth_twin_prime(n)                   the nth twin prime (n=1 returns 3)
  nth_twin_prime_approx(n)            fast approximate nth twin prime
  legendre_phi(n,a)                   # below n not div by first a primes
  prime_precalc(n)                    precalculate primes to n

=head2 FACTORING

  factor(n)                           array of prime factors of n
  factor_exp(n)                       array of [p,k] factors p^k
  divisors(n)                         array of divisors of n
  divisor_sum(n)                      sum of divisors
  divisor_sum(n,k)                    sum of k-th power of divisors
  divisor_sum(n,sub{...})             sum of code run for each divisor
  znlog(a, g, p)                      solve k in a = g^k mod p

=head2 ITERATORS

  forprimes { ... } [start,] end      loop over primes in range
  forcomposites { ... } [start,] end  loop over composites in range
  foroddcomposites {...} [start,] end loop over odd composites in range
  fordivisors { ... } n               loop over the divisors of n
  forpart { ... } n [,{...}]          loop over integer partitions
  forcomb { ... } n, k                loop over combinations
  forperm { ... } n                   loop over permutations
  prime_iterator                      returns a simple prime iterator
  prime_iterator_object               returns a prime iterator object

=head2 RANDOM PRIMES

  random_prime([start,] end)          random prime in a range
  random_ndigit_prime(n)              random prime with n digits
  random_nbit_prime(n)                random prime with n bits
  random_strong_prime(n)              random strong prime with n bits
  random_proven_prime(n)              random n-bit prime with proof
  random_proven_prime_with_cert(n)    as above and include certificate
  random_maurer_prime(n)              random n-bit prime w/ Maurer's alg.
  random_maurer_prime_with_cert(n)    as above and include certificate
  random_shawe_taylor_prime(n)        random n-bit prime with S-T alg.
  random_shawe_taylor_prime_with_cert(n) as above including certificate

=head2 MATH

  vecsum(@list)                       integer sum of list
  vecprod(@list)                      integer product of list
  vecmin(@list)                       minimum of list of integers
  vecmax(@list)                       maximum of list of integers
  is_power(n)                         return k if n = p^k for integer p
  gcd(@list)                          greatest common divisor
  lcm(@list)                          least common multiple
  gcdext(x,y)                         return (u,v,d) where u*x+v*y=d
  chinese([a,mod1],[b,mod2],...)      Chinese Remainder Theorem
  primorial(n)                        product of primes below n
  pn_primorial(n)                     product of first n primes
  factorial(n)                        product of first n integers: n!
  binomial(n,k)                       binomial coefficient
  partitions(n)                       number of integer partitions
  valuation(n,k)                      number of times n is divisible by k
  kronecker(a,b)                      Kronecker (Jacobi) symbol
  invmod(a,n)                         inverse of a modulo n
  moebius(n)                          Moebius function of n
  moebius(beg, end)                   array of Moebius in range
  mertens(n)                          sum of Moebius for 1 to n
  euler_phi(n)                        Euler totient of n
  euler_phi(beg, end)                 Euler totient for a range
  jordan_totient(n,k)                 Jordan's totient
  carmichael_lambda(n)                Carmichael's Lambda function
  exp_mangoldt                        exponential of Mangoldt function
  liouville(n)                        Liouville function
  znorder(a,n)                        multiplicative order of a mod n
  znprimroot(n)                       smallest primitive root
  chebyshev_theta(n)                  first Chebyshev function
  chebyshev_psi(n)                    second Chebyshev function
  consecutive_integer_lcm(n)          lcm(1 .. n)
  lucas_sequence(n, P, Q, k)          (U_k,V_k,Q_k) for Lucas(P,Q) mod n
  bernfrac(n)                         Bernoulli number as (num,den)
  bernreal(n)                         Bernoulli number as BigFloat
  
=head2 NON-INTEGER MATH

  ExponentialIntegral(x)              Ei(x)
  LogarithmicIntegral(x)              li(x)
  RiemannZeta(x)                      ζ(s)-1, real-valued Riemann Zeta
  RiemannR(x)                         Riemann's R function
  LambertW(k)                         Lambert W: W for C<k = W exp(W)>
  Pi([n])                             The constant π (NV or n digits)
  
=head2 SUPPORT

  prime_get_config                    gets hash ref of current settings
  prime_set_config(%hash)             sets parameters
  prime_memfree                       frees any cached memory


=head1 COPYRIGHT

Copyright 2011-2014 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
