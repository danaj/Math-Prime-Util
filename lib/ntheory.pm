package ntheory;
use strict;
use warnings;

BEGIN {
  $ntheory::AUTHORITY = 'cpan:DANAJ';
  $ntheory::VERSION = '0.73';
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

Tags:
  :all         to import almost all functions
  :rand        to import rand, srand, irand, irand64

=head2 PRIMALITY

  is_prob_prime(n)                    primality test (BPSW)
  is_prime(n)                         primality test (BPSW + extra)
  is_provable_prime(n)                primality test with proof
  is_provable_prime_with_cert(n)      primality test: (isprime,cert)
  prime_certificate(n)                as above with just certificate
  verify_prime(cert)                  verify a primality certificate
  is_mersenne_prime(p)                is 2^p-1 prime or composite
  is_aks_prime(n)                     AKS deterministic test (slow)
  is_ramanujan_prime(n)               is n a Ramanujan prime
  is_gaussian_prime(a,b)              is a+bi a Gaussian prime

=head2 PROBABLE PRIME TESTS

  is_pseudoprime(n,bases)                  Fermat probable prime test
  is_euler_pseudoprime(n,bases)            Euler test to bases
  is_euler_plumb_pseudoprime(n)            Euler Criterion test
  is_strong_pseudoprime(n,bases)           Miller-Rabin test to bases
  is_lucas_pseudoprime(n)                  Lucas test
  is_strong_lucas_pseudoprime(n)           strong Lucas test
  is_almost_extra_strong_lucas_pseudoprime(n, [incr])   AES Lucas test
  is_extra_strong_lucas_pseudoprime(n)     extra strong Lucas test
  is_frobenius_pseudoprime(n, [a,b])       Frobenius quadratic test
  is_frobenius_underwood_pseudoprime(n)    combined PSP and Lucas
  is_frobenius_khashin_pseudoprime(n)      Khashin's 2013 Frobenius test
  is_perrin_pseudoprime(n [,r])            Perrin test
  is_catalan_pseudoprime(n)                Catalan test
  is_bpsw_prime(n)                         combined SPSP-2 and ES Lucas
  miller_rabin_random(n, ntests)           perform random-base MR tests

=head2 PRIMES

  primes([start,] end)                array ref of primes
  prime_powers([start,] end)          array ref of prime powers
  twin_primes([start,] end)           array ref of twin primes
  semi_primes([start,] end)           array ref of semiprimes
  almost_primes(k, [start,] end)      array ref of k-almost-primes
  omega_primes(k, [start,] end)       array ref of k-omega-primes
  ramanujan_primes([start,] end)      array ref of Ramanujan primes
  sieve_prime_cluster(start, end, @C) list of prime k-tuples
  sieve_range(n, width, depth)        sieve out small factors to depth
  next_prime(n)                       next prime > n
  prev_prime(n)                       previous prime < n
  next_prime_power(n)                 next prime power > n
  prev_prime_power(n)                 previous prime power < n
  prime_count(n)                      count of primes <= n
  prime_count(start, end)             count of primes in range
  prime_count_lower(n)                fast lower bound for prime count
  prime_count_upper(n)                fast upper bound for prime count
  prime_count_approx(n)               fast approximate prime count
  prime_power_count(n)                count of prime powers <= n
  prime_power_count(start, end)       count of prime powers in range
  prime_power_count_lower(n)          fast lower bound for prime power count
  prime_power_count_upper(n)          fast upper bound for prime power count
  prime_power_count_approx(n)         fast approximate prime power count
  nth_prime(n)                        the nth prime (n=1 returns 2)
  nth_prime_lower(n)                  fast lower bound for nth prime
  nth_prime_upper(n)                  fast upper bound for nth prime
  nth_prime_approx(n)                 fast approximate nth prime
  nth_prime_power(n)                  the nth prime power (n=1 returns 2)
  nth_prime_power_lower(n)            fast lower bound for nth prime power
  nth_prime_power_upper(n)            fast upper bound for nth prime power
  nth_prime_power_approx(n)           fast approximate nth prime power
  twin_prime_count(n)                 count of twin primes <= n
  twin_prime_count(start, end)        count of twin primes in range
  twin_prime_count_approx(n)          fast approximate twin prime count
  nth_twin_prime(n)                   the nth twin prime (n=1 returns 3)
  nth_twin_prime_approx(n)            fast approximate nth twin prime
  semiprime_count(n)                  count of semiprimes <= n
  semiprime_count(start, end)         count of semiprimes in range
  semiprime_count_approx(n)           fast approximate semiprime count
  nth_semiprime(n)                    the nth semiprime
  nth_semiprime_approx(n)             fast approximate nth semiprime
  almost_prime_count(k,n)             count of k-almost-primes
  almost_prime_count_approx(k,n)      fast approximate k-almost-prime count
  almost_prime_count_lower(k,n)       fast k-almost-prime count lower bound
  almost_prime_count_upper(k,n)       fast k-almost-prime count upper bound
  nth_almost_prime(k,n)               the nth number with exactly k factors
  nth_almost_prime_approx(k,n)        fast approximate nth k-almost prime
  nth_almost_prime_lower(k,n)         fast nth k-almost prime lower bound
  nth_almost_prime_upper(k,n)         fast nth k-almost prime upper bound
  omega_prime_count(k,n)              count divisible by exactly k primes
  nth_omega_prime(k,n)                the nth number div by exactly k primes
  ramanujan_prime_count(n)            count of Ramanujan primes <= n
  ramanujan_prime_count(start, end)   count of Ramanujan primes in range
  ramanujan_prime_count_lower(n)      fast lower bound for Ramanujan count
  ramanujan_prime_count_upper(n)      fast upper bound for Ramanujan count
  ramanujan_prime_count_approx(n)     fast approximate Ramanujan count
  nth_ramanujan_prime(n)              the nth Ramanujan prime (Rn)
  nth_ramanujan_prime_lower(n)        fast lower bound for Rn
  nth_ramanujan_prime_upper(n)        fast upper bound for Rn
  nth_ramanujan_prime_approx(n)       fast approximate Rn
  legendre_phi(n,a)                   # below n not div by first a primes
  inverse_li(n)                       integer inverse logarithmic integral
  inverse_li_nv(x)                    float inverse logarithmic integral
  prime_precalc(n)                    precalculate primes to n
  sum_primes([start,] end)            return summation of primes in range
  print_primes(start,end[,fd])        print primes to stdout or fd

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
  forsemiprimes {...} [start,] end    loop over semiprimes in range
  foralmostprimes {...} k,[beg,],end  loop over k-almost-primes in range
  forfactored {...} [start,] end      loop with factors
  forsquarefree {...} [start,] end    loop with factors of square-free n
  forsquarefreeint {...} [start,] end loop over square-free n
  fordivisors { ... } n               loop over the divisors of n
  forpart { ... } n [,{...}]          loop over integer partitions
  forcomp { ... } n [,{...}]          loop over integer compositions
  forcomb { ... } n, k                loop over combinations
  forperm { ... } n                   loop over permutations
  formultiperm { ... } \@n            loop over multiset permutations
  forderange { ... } n                loop over derangements
  forsetproduct { ... } \@a[,...]     loop over Cartesian product of lists
  prime_iterator                      returns a simple prime iterator
  prime_iterator_object               returns a prime iterator object
  lastfor                             stop iteration of for.... loop

=head2 RANDOM NUMBERS

  irand                               random 32-bit integer
  irand64                             random 64-bit integer
  drand([limit])                      random NV in [0,1) or [0,limit)
  random_bytes(n)                     string with n random bytes
  entropy_bytes(n)                    string with n entropy-source bytes
  urandomb(n)                         random integer less than 2^n
  urandomm(n)                         random integer less than n
  csrand(data)                        seed the CSPRNG with binary data
  srand([seed])                       simple seed (exported with :rand)
  rand([limit])                       alias for drand (exported with :rand)
  random_factored_integer(n)          random [1..n] and array ref of factors

=head2 RANDOM PRIMES

  random_prime([start,] end)          random prime in a range
  random_ndigit_prime(n)              random prime with n digits
  random_nbit_prime(n)                random prime with n bits
  random_safe_prime(n)                random safe prime with n bits
  random_strong_prime(n)              random strong prime with n bits
  random_proven_prime(n)              random n-bit prime with proof
  random_proven_prime_with_cert(n)    as above and include certificate
  random_maurer_prime(n)              random n-bit prime w/ Maurer's alg.
  random_maurer_prime_with_cert(n)    as above and include certificate
  random_shawe_taylor_prime(n)        random n-bit prime with S-T alg.
  random_shawe_taylor_prime_with_cert(n) as above including certificate
  random_unrestricted_semiprime(n)    random n-bit semiprime
  random_semiprime(n)                 as above with equal size factors

=head2 LISTS

  vecsum(@list)                       integer sum of list
  vecprod(@list)                      integer product of list
  vecmin(@list)                       minimum of list of integers
  vecmax(@list)                       maximum of list of integers
  vecuniq(@list)                      remove duplicates from list of integers
  vecsort(@list)                      numerically sort a list of integers
  vecsorti(\@list)                    in-place numeric sort a list ref
  vecextract(\@list, mask)            select from list based on mask
  vecequal(\@list1, \@list2)          compare equality of two array refs
  vecreduce { ... } @list             reduce / left fold applied to list
  vecall { ... } @list                return true if all are true
  vecany { ... } @list                return true if any are true
  vecnone { ... } @list               return true if none are true
  vecnotall { ... } @list             return true if not all are true
  vecfirst { ... } @list              return first value that evals true
  vecfirstidx { ... } @list           return first index that evals true
  vecmex(@list)                       return least non-neg value not in list
  vecpmex(@list)                      return least positive value not in list
  vecsample(k,@list)                  return k random elements of list

  toset(\@A)                          convert to unique sorted integer list
  setinsert(\@A,$v)                   insert integer v into integer set A
  setinsert(\@A,\@B)                  insert list B values into integer set A
  setremove(\@A,$v)                   remove integer v from integer set A
  setremove(\@A,\@B)                  remove list B values from integer set A
  setinvert(\@A,$v)                   if v is in set A, remove, otherwise add
  setinvert(\@A,\@B)                  invert for all values in integer list B
  setcontains(\@A,$v)                 is integer v in integer set A
  setcontains(\@A,\@B)                is int set B a subset of int set A
  setcontainsany(\@A,\@B)             is any value in B in int set A
  setbinop { ... } \@A[,\@B]          apply operation to all a,b [a:A,b:B]
  sumset \@A[,\@B]                    apply a+b to all a,b [a:A,b:B]
  setunion(\@A,\@B)                   union of two integer lists
  setintersect(\@A,\@B)               intersection of two integer lists
  setminus(\@A,\@B)                   difference of two integer lists
  setdelta(\@A,\@B)                   symmetric difference of two int lists
  is_sidon_set(\@L)                   is integer list L a Sidon set
  is_sumfree_set(\@L)                 is integer list L a sum-free set
  set_is_disjoint(\@A,\@B)            is set B disjoint from set A
  set_is_equal(\@A,\@B)               is set B equal to set A
  set_is_subset(\@A,\@B)              is set B a subset of set A
  set_is_proper_subset(\@A,\@B)       is set B a proper subset of set A
  set_is_superset(\@A,\@B)            is set B a superet of set A
  set_is_proper_superset(\@A,\@B)     is set B a proper superet of set A
  set_is_proper_intersection(\@A,\@B) is set B a proper intersection of set A

=head2 MATH

  todigits(n[,base[,len]])            convert n to digit array in base
  todigitstring(n[,base[,len]])       convert n to string in base
  fromdigits(\@d,[,base])             convert base digit vector to number
  fromdigits(str,[,base])             convert base digit string to number
  sumdigits(n)                        sum of digits, with optional base
  tozeckendorf(n)                     convert n to Zeckendorf/Fibbinary
  fromzeckendorf(str)                 convert Zeckendorf binary str to num
  is_odd(n)                           return 1 if n is odd, 0 otherwise
  is_even(n)                          return 1 if n is even, 0 otherwise
  is_divisible(n,d)                   return 1 if n divisible by d
  is_congruent(n,c,d)                 return 1 if n is congruent to c mod d
  is_qr(a,n)                          return 1 if a is quadratic non-res mod n
  is_square(n)                        return 1 if n is a perfect square
  is_power(n)                         return k if n = c^k for integer c
  is_power(n,k)                       return 1 if n = c^k for integer c, k
  is_power(n,k,\$root)                as above but also set $root to c
  is_perfect_power(n)                 return 1 if n = c^k for c != 0, k > 1
  is_prime_power(n)                   return k if n = p^k for prime p, k > 0
  is_prime_power(n,\$p)               as above but also set $p to p
  is_square_free(n)                   return true if no repeated factors
  is_powerfree(n[,k])                 is n free of any k-th powers
  is_cyclic(n)                        does n have only one group of order n
  is_carmichael(n)                    is n a Carmichael number
  is_quasi_carmichael(n)              is n a quasi-Carmichael number
  is_primitive_root(r,n)              is r a primitive root mod n
  is_pillai(n)                        v where  v! % n == n-1  and  n % v != 1
  is_semiprime(n)                     does n have exactly 2 prime factors
  is_almost_prime(k,n)                does n have exactly k prime factors
  is_omega_prime(k,n)                 is n divisible by exactly k primes
  is_chen_prime(n)                    is n prime and n+2 prime or semiprime
  is_polygonal(n,k)                   is n a k-polygonal number
  is_polygonal(n,k,\$root)            as above but also set $root
  is_sum_of_squares(n[,k])            is n a sum of k (def 2) squares
  is_congruent_number(n)              is n a congruent number
  is_perfect_number(n)                is n equal to sum of its proper divisors
  is_fundamental(d)                   is d a fundamental discriminant
  is_totient(n)                       is n = euler_phi(x) for some x
  is_lucky(n)                         is n a lucky number
  is_happy(n)                         if n a happy number, returns height
  is_happy(n,base,exponent)           if n a S_b_e happy number, returns height
  is_smooth(n,k)                      is n a k-smooth number
  is_rough(n,k)                       is n a k-rough number
  is_powerful(n[,k])                  is n a k-powerful number
  is_practical(n)                     is n a practical number
  is_delicate_prime(n)                is n a digitally delicate prime
  powint(a,b)                         signed integer a^b
  mulint(a,b)                         signed integer a * b
  addint(a,b)                         signed integer a + b
  subint(a,b)                         signed integer a - b
  add1int(n)                          signed integer n + 1
  sub1int(n)                          signed integer n - 1
  divint(a,b)                         signed integer a / b     (floor)
  modint(a,b)                         signed integer a % b     (floor)
  cdivint(a,b)                        signed integer a / b     (ceilint)
  divrem(a,b)                         return (quot,rem) of a/b (Euclidian)
  fdivrem(a,b)                        return (quot,rem) of a/b (floored)
  cdivrem(a,b)                        return (quot,rem) of a/b (ceiling)
  tdivrem(a,b)                        return (quot,rem) of a/b (truncated)
  lshiftint(n,k)                      left shift n by k bits
  rshiftint(n,k)                      right shift n by k bits (truncate)
  rashiftint(n,k)                     right shift n by k bits (floor)
  absint(n)                           integer absolute value
  negint(n)                           integer negation
  cmpint(a,b)                         integer comparison (like <=>)
  signint(n)                          integer sign (-1,0,1)
  sqrtint(n)                          integer square root
  rootint(n,k)                        integer k-th root
  rootint(n,k,\$rk)                   as above but also set $rk to r^k
  logint(n,b)                         integer logarithm
  logint(n,b,\$be)                    as above but also set $be to b^e
  gcd(@list)                          greatest common divisor
  lcm(@list)                          least common multiple
  gcdext(x,y)                         return (u,v,d) where u*x+v*y=d
  chinese([a,mod1],[b,mod2],...)      Chinese Remainder Theorem
  chinese2([a,mod1],[b,mod2],...)     Chinese Remainder Theorem
  frobenius_number(@list)             Frobenius Number of a set
  primorial(n)                        product of primes below n
  pn_primorial(n)                     product of first n primes
  factorial(n)                        product of first n integers: n!
  factorialmod(n,m)                   factorial mod m
  subfactorial(n)                     count of derangements of n objects
  binomial(n,k)                       binomial coefficient
  binomialmod(n,k,m)                  binomial(n,k) mod m
  falling_factorial(x,n)              falling factorial
  rising_factorial(x,n)               rising factorial
  partitions(n)                       number of integer partitions
  valuation(n,k)                      number of times n is divisible by k
  hammingweight(n)                    population count (# of binary 1s)
  kronecker(a,b)                      Kronecker (Jacobi) symbol
  negmod(a,n)                         -a mod n
  addmod(a,b,n)                       a + b mod n
  submod(a,b,n)                       a - b mod n
  mulmod(a,b,n)                       a * b mod n
  muladdmod(a,b,c,n)                  a * b + c mod n
  mulsubmod(a,b,c,n)                  a * b - c mod n
  divmod(a,b,n)                       a / b mod n
  powmod(a,b,n)                       a ^ b mod n
  invmod(a,n)                         inverse of a modulo n
  sqrtmod(a,n)                        modular square root
  rootmod(a,k,n)                      modular k-th root
  allsqrtmod(a,n)                     list of all modular square roots
  allrootmod(a,k,n)                   list of all modular k-th roots
  cornacchia(d,n)                     find x,y for x^2 + d * y^2 = n
  prime_bigomega(n)                   number of prime factors
  prime_omega(n)                      number of distinct prime factors
  moebius(n)                          Moebius function of n
  moebius(beg, end)                   list of Moebius in range
  mertens(n)                          sum of Moebius for 1 to n
  euler_phi(n)                        Euler totient of n
  euler_phi(beg, end)                 Euler totient for a range
  inverse_totient(n)                  image of Euler totient
  jordan_totient(k,n)                 Jordan's totient
  sumtotient(n)                       sum of Euler totient for 1 to n
  carmichael_lambda(n)                Carmichael's Lambda function
  ramanujan_sum(k,n)                  Ramanujan's sum
  exp_mangoldt                        exponential of Mangoldt function
  liouville(n)                        Liouville function
  sumliouville(n)                     sum of Liouville for 1 to n
  znorder(a,n)                        multiplicative order of a mod n
  znprimroot(n)                       smallest primitive root
  qnr(n)                              least quadratic non-residue
  chebyshev_theta(n)                  first Chebyshev function
  chebyshev_psi(n)                    second Chebyshev function
  hclassno(n)                         Hurwitz class number H(n) * 12
  ramanujan_tau(n)                    Ramanujan's Tau function
  consecutive_integer_lcm(n)          lcm(1 .. n)
  lucasu(P, Q, k)                     U_k for Lucas(P,Q)
  lucasv(P, Q, k)                     V_k for Lucas(P,Q)
  lucasuv(P, Q, k)                    (U_k,V_k) for Lucas(P,Q)
  lucasumod(P, Q, k, n)               U_k for Lucas(P,Q) mod n
  lucasvmod(P, Q, k, n)               V_k for Lucas(P,Q) mod n
  lucasuvmod(P, Q, k, n)              (U_k,V_k,Q^k) for Lucas(P,Q) mod n
  pisano_period(n)                    The period of Fibonacci numbers mod n
  bernfrac(n)                         Bernoulli number as (num,den)
  bernreal(n)                         Bernoulli number as BigFloat
  harmfrac(n)                         Harmonic number as (num,den)
  harmreal(n)                         Harmonic number as BigFloat
  stirling(n,m,[type])                Stirling numbers of 1st or 2nd type
  fubini(n)                           Fubini (Ordered Bell) number
  numtoperm(n,k)                      kth lexico permutation of n elems
  permtonum([a,b,...])                permutation number of given perm
  randperm(n,[k])                     random permutation of n elems
  shuffle(...)                        random permutation of an array
  lucky_numbers(n)                    array ref of lucky sieve up to n
  lucky_count(n)                      count of lucky numbers <= n
  lucky_count(start, end)             count of lucky numbers in range
  lucky_count_lower(n)                fast lower bound for lucky count
  lucky_count_upper(n)                fast upper bound for lucky count
  lucky_count_approx(n)               fast approximate lucky count
  nth_lucky(n)                        nth entry in lucky sieve
  nth_lucky_lower(n)                  fast lower bound for nth lucky number
  nth_lucky_upper(n)                  fast upper bound for nth lucky number
  nth_lucky_approx(n)                 fast approximate nth lucky number
  powerful_numbers([lo,]hi[,k])       array ref of k-powerful lo to hi
  powerful_count(n[,k])               count of k-powerful numbers <= n
  sumpowerful(n[,k])                  sum of k-powerful numbers <= n
  nth_powerful(n[,k])                 the nth k-powerful number
  next_perfect_power(n)               the next perfect power > n
  prev_perfect_power(n)               the previous perfect power < n
  perfect_power_count(n)              count of perfect powers <= n
  perfect_power_count(start, end)     count of perfect powers in range
  perfect_power_count_lower(n)        fast lower bound for perf power count
  perfect_power_count_upper(n)        fast upper bound for perf power count
  perfect_power_count_approx(n)       fast approximate perfect power count
  nth_perfect_power(n)                the nth perfect power
  nth_perfect_power_lower(n)          fast lower bound for nth perfect power
  nth_perfect_power_upper(n)          fast upper bound for nth perfect power
  nth_perfect_power_approx(n)         fast approximate nth perfect power
  next_chen_prime(n)                  next Chen prime > n
  smooth_count(n,k)                   count of k-smooth numbers <= n
  rough_count(n,k)                    count of k-rough numbers <= n
  powerfree_count(n[,k])              count of k-powerfree numbers <= n
  nth_powerfree(n[,k])                the nth k-powerfree number
  powerfree_sum(n[,k])                sum of k-powerfree numbers <= n
  powerfree_part(n[,k])               remove excess powers so n is k-free
  powerfree_part_sum(n[,k])           sum of k-powerfree parts for 1 to n
  squarefree_kernel(n)                integer radical of |n|
  powersum(n,k)                       sum of kth powers from 1 to n

=head2 RATIONALS

  contfrac(n,d)                       list of continued fraction for n/d
  next_calkin_wilf(n,d)               next breadth-first CW rational
  next_stern_brocot(n,d)              next breadth-first SB rational
  calkin_wilf_n(n,d)                  index of breadth-first CW rational
  stern_brocot_n(n,d)                 index of breadth-first SB rational
  nth_calkin_wilf(n)                  CW rational at breadth-first index n
  nth_stern_brocot(n)                 SB rational at breadth-first index n
  nth_stern_diatomic(n)               Stern's Diatomic series; fusc(n)
  farey(n)                            list of Farey sequence order n
  farey(n,k)                          k'th entry of Farey sequence order n
  next_farey(n,[p,q])                 next order-n rational after p/q
  farey_rank(n,[p,q])                 number of F_n less than p/q

=head2 NON-INTEGER MATH

  ExponentialIntegral(x)              Ei(x)
  LogarithmicIntegral(x)              li(x)
  RiemannZeta(x)                      ζ(s)-1, real-valued Riemann Zeta
  RiemannR(x)                         Riemann's R function
  LambertW(k)                         Lambert W: solve for W in k = W exp(W)
  Pi([n])                             The constant π (NV or n digits)

=head2 SUPPORT

  prime_get_config                    gets hash ref of current settings
  prime_set_config(%hash)             sets parameters
  prime_memfree                       frees any cached memory


=head1 COPYRIGHT

Copyright 2011-2025 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
