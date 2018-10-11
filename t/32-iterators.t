#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/primes prev_prime next_prime
                         forprimes forcomposites foroddcomposites fordivisors
                         forpart forcomp forcomb forperm forderange formultiperm
                         forfactored forsquarefree forsemiprimes
                         lastfor
                         is_power is_semiprime vecsum sqrtint
                         prime_iterator prime_iterator_object/;
use Math::BigInt try => "GMP,Pari";
use Math::BigFloat;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $broken64 = (18446744073709550592 == ~0);

plan tests => 8        # forprimes errors
            + 12 + 7   # forprimes simple
            + 3        # forcomposites simple
            + 2        # fordivisors simple
            + 3        # iterator errors
            + 7        # iterator simple
            + 1        # other forprimes
            + 2        # forprimes/iterator nesting
            + 3        # forprimes BigInt/BigFloat
            + 3        # oo iterator errors
            + 7        # oo iterator simple
            + 25       # oo iterator methods
            + 12       # lastfor
            + 5        # forfactored and forsquarefree
            + 1        # forsemiprimes
            + 0;

ok(!eval { forprimes { 1 } undef; },   "forprimes undef");
ok(!eval { forprimes { 1 } 2, undef; },   "forprimes 2,undef");
ok(!eval { forprimes { 1 } undef, 2; },   "forprimes 2,undef");
# This is caught at compile type because of the prototype
#ok(!eval { forprimes { 1 } 2, 3, 4; },   "forprimes 2,3,4");
ok(!eval { forprimes { 1 } -2, 3; },   "forprimes -2,3");
ok(!eval { forprimes { 1 } 2, -3; },   "forprimes 2,-3");
ok(!eval { forprimes { 1 } "abc"; },   "forprimes abc");
ok(!eval { forprimes { 1 } 2, "abc"; },   "forprimes 2, abc");
ok(!eval { forprimes { 1 } 5.6; },   "forprimes abc");

{my @t; forprimes {push @t,$_} 1; is_deeply( [@t], [], "forprimes 1" ); }
{my @t; forprimes {push @t,$_} 2; is_deeply( [@t], [2], "forprimes 3" ); }
{my @t; forprimes {push @t,$_} 3; is_deeply( [@t], [2,3], "forprimes 3" ); }
{my @t; forprimes {push @t,$_} 4; is_deeply( [@t], [2,3], "forprimes 4" ); }
{my @t; forprimes {push @t,$_} 5; is_deeply( [@t], [2,3,5], "forprimes 5" ); }
{my @t; forprimes {push @t,$_} 3,5; is_deeply( [@t], [3,5], "forprimes 3,5" ); }
{my @t; forprimes {push @t,$_} 3,6; is_deeply( [@t], [3,5], "forprimes 3,6" ); }
{my @t; forprimes {push @t,$_} 3,7; is_deeply( [@t], [3,5,7], "forprimes 3,7" ); }
{my @t; forprimes {push @t,$_} 5,7; is_deeply( [@t], [5,7], "forprimes 5,7" ); }
{my @t; forprimes {push @t,$_} 6,7; is_deeply( [@t], [7], "forprimes 6,7" ); }
{my @t; forprimes {push @t,$_} 5,11; is_deeply( [@t], [5,7,11], "forprimes 5,11" ); }
{my @t; forprimes {push @t,$_} 7,11; is_deeply( [@t], [7,11], "forprimes 7,11" ); }

{ my @t; forprimes { push @t, $_ } 50;
  is_deeply( [@t], [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47], "forprimes 50" );
}
{ my @t; forprimes { push @t, $_ } 2,20;
  is_deeply( [@t], [2,3,5,7,11,13,17,19], "forprimes 2,20" );
}
{ my @t; forprimes { push @t, $_ } 20,30;
  is_deeply( [@t], [23,29], "forprimes 20,30" );
}
{ my @t; forprimes { push @t, $_ } 199, 223;
  is_deeply( [@t], [199,211,223], "forprimes 199,223" );
}
{ my @t; forprimes { push @t, $_ } 31398, 31468;
  is_deeply( [@t], [], "forprimes 31398,31468 (empty region)" );
}
{ my @t; forprimes { push @t, $_ } 2147483647,2147483659;
  is_deeply( [@t], [2147483647,2147483659], "forprimes 2147483647,2147483659" );
}
{ my @t; forprimes { push @t, $_ } 3842610774,3842611326;
  is_deeply( [@t], [3842611109,3842611139,3842611163,3842611181,3842611211,3842611229,3842611249,3842611259,3842611261,3842611291,3842611301], "forprimes 3842610774,3842611326" );
}
{ my @t; forcomposites { push @t, $_ } 2147483647,2147483659;
  is_deeply( [@t], [qw/2147483648 2147483649 2147483650 2147483651 2147483652 2147483653 2147483654 2147483655 2147483656 2147483657 2147483658/], "forcomposites 2147483647,2147483659" );
}
{ my @t; forcomposites { push @t, $_ } 50;
  is_deeply( [@t], [qw/4 6 8 9 10 12 14 15 16 18 20 21 22 24 25 26 27 28 30 32 33 34 35 36 38 39 40 42 44 45 46 48 49 50/], "forcomposites 50" );
}
{ my @t; forcomposites { push @t, $_ } 200,410;
  is_deeply( [@t], [qw/200 201 202 203 204 205 206 207 208 209 210 212 213 214 215 216 217 218 219 220 221 222 224 225 226 228 230 231 232 234 235 236 237 238 240 242 243 244 245 246 247 248 249 250 252 253 254 255 256 258 259 260 261 262 264 265 266 267 268 270 272 273 274 275 276 278 279 280 282 284 285 286 287 288 289 290 291 292 294 295 296 297 298 299 300 301 302 303 304 305 306 308 309 310 312 314 315 316 318 319 320 321 322 323 324 325 326 327 328 329 330 332 333 334 335 336 338 339 340 341 342 343 344 345 346 348 350 351 352 354 355 356 357 358 360 361 362 363 364 365 366 368 369 370 371 372 374 375 376 377 378 380 381 382 384 385 386 387 388 390 391 392 393 394 395 396 398 399 400 402 403 404 405 406 407 408 410/], "forcomposites 200,410" );
}
{
  my $a = 0;
  fordivisors { $a += $_ + $_*$_ } 54321;
  is($a, 3287796520, "fordivisors: d|54321: a+=d+d^2");
  # Matches Math::Pari:
  #   my $a = PARI(0); my $j; fordiv(54321,$j,sub { $a += $j + $j**2 });
}
{
  # Pari: v=List(); for(n=1, 50, fordiv(n, d, listput(v, d))); Vec(v)
  my @A027750 = (1,1,2,1,3,1,2,4,1,5,1,2,3,6,1,7,1,2,4,8,1,3,9,1,2,5,10,1,11,1,2,3,4,6,12,1,13,1,2,7,14,1,3,5,15,1,2,4,8,16,1,17,1,2,3,6,9,18,1,19,1,2,4,5,10,20,1,3,7,21,1,2,11,22,1,23,1,2,3,4,6,8,12,24,1,5,25,1,2,13,26,1,3,9,27,1,2,4,7,14,28,1,29,1,2,3,5,6,10,15,30,1,31,1,2,4,8,16,32,1,3,11,33,1,2,17,34,1,5,7,35,1,2,3,4,6,9,12,18,36,1,37,1,2,19,38,1,3,13,39,1,2,4,5,8,10,20,40,1,41,1,2,3,6,7,14,21,42,1,43,1,2,4,11,22,44,1,3,5,9,15,45,1,2,23,46,1,47,1,2,3,4,6,8,12,16,24,48,1,7,49,1,2,5,10,25,50);
  my @a;
  do { fordivisors { push @a, $_ } $_ } for 1..50;
  is_deeply(\@a, \@A027750, "A027750 using fordivisors");
}

ok(!eval { prime_iterator(-2); }, "iterator -2");
ok(!eval { prime_iterator("abc"); }, "iterator abc");
ok(!eval { prime_iterator(4.5); }, "iterator 4.5");

{ my $it = prime_iterator();
  is_deeply( [map { $it->() } 1..10], [2,3,5,7,11,13,17,19,23,29], "iterator first 10 primes" );
}
{my $it = prime_iterator(47);
  is_deeply( [map { $it->() } 1..5], [47,53,59,61,67], "iterator 5 primes starting at 47" );
}
{my $it = prime_iterator(199);
  is_deeply( [map { $it->() } 1..3], [199,211,223], "iterator 3 primes starting at 199" );
}
{my $it = prime_iterator(200);
  is_deeply( [map { $it->() } 1..3], [211,223,227], "iterator 3 primes starting at 200" );
}
{my $it = prime_iterator(31397);
  is_deeply( [map { $it->() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31397" );
}
{my $it = prime_iterator(31396);
  is_deeply( [map { $it->() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31396" );
}
{my $it = prime_iterator(31398);
  is_deeply( [map { $it->() } 1..3], [31469,31477,31481], "iterator 3 primes starting at 31398" );
}

# Make sure things work when the type of $_ changes
{ my $sum = 0;
  forprimes { $sum += int(12345/$_) } 1000;
  is(27053, $sum, "forprimes handles \$_ type changes");
}

# For fun, nest them.
{
  my @t;
  forprimes {
    forprimes {
      forprimes { push @t, $_ } $_,$_+10;
    } 10*$_,10*$_+10;
  } 10;
  is_deeply( [@t], [qw/23 29 31 29 31 37 31 37 41 37 41 43 47 53 59 61 59 61 67 71 73 79 73 79 83 79 83 89/], "triple nested forprimes" );
}
{
  my @t;
  my $ita = prime_iterator();
  while ((my $a = $ita->()) <= 10) {
    my $itb = prime_iterator(10*$a);
    while ((my $b = $itb->()) <= 10*$a+10) {
      my $itc = prime_iterator($b);
      while ((my $c = $itc->()) <= $b+10) {
        push @t, $c;
      }
    }
  }
  is_deeply( [@t], [qw/23 29 31 29 31 37 31 37 41 37 41 43 47 53 59 61 59 61 67 71 73 79 73 79 83 79 83 89/], "triple nested iterator" );
}

# With BigInt and BigFloat objects
{ my @t;
  forprimes { push @t, $_ } Math::BigInt->new("5"), Math::BigInt->new("11");
  is_deeply( [@t], [5,7,11], "forprimes with BigInt range" );
}
{ my @t;
  forprimes { push @t, $_ } Math::BigFloat->new("5"), Math::BigFloat->new("11");
  is_deeply( [@t], [5,7,11], "forprimes with BigFloat range" );
}
{my $it = prime_iterator(Math::BigInt->new("68719476736"));
  is_deeply( [map { $it->() } 1..3],
             [68719476767,68719476851,68719476853],
             "iterator 3 primes with BigInt start" );
}

# Test new object iterator
ok(!eval { prime_iterator_object(-2); }, "iterator -2");
ok(!eval { prime_iterator_object("abc"); }, "iterator abc");
ok(!eval { prime_iterator_object(4.5); }, "iterator 4.5");

{ my $it = prime_iterator_object();
  is_deeply( [map { $it->iterate() } 1..10], [2,3,5,7,11,13,17,19,23,29], "iterator first 10 primes" );
}
{my $it = prime_iterator_object(47);
  is_deeply( [map { $it->iterate() } 1..5], [47,53,59,61,67], "iterator 5 primes starting at 47" );
}
{my $it = prime_iterator_object(199);
  is_deeply( [map { $it->iterate() } 1..3], [199,211,223], "iterator 3 primes starting at 199" );
}
{my $it = prime_iterator_object(200);
  is_deeply( [map { $it->iterate() } 1..3], [211,223,227], "iterator 3 primes starting at 200" );
}
{my $it = prime_iterator_object(31397);
  is_deeply( [map { $it->iterate() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31397" );
}
{my $it = prime_iterator_object(31396);
  is_deeply( [map { $it->iterate() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31396" );
}
{my $it = prime_iterator_object(31398);
  is_deeply( [map { $it->iterate() } 1..3], [31469,31477,31481], "iterator 3 primes starting at 31398" );
}

{
  my $it = prime_iterator_object;
  do { $it->next } for 1..10;
  is( $it->value(), 31, "iterator object moved forward 10 now returns 31");
  $it->prev;
  is( $it->value(), 29, "iterator object moved back now returns 29");
  is( $it->iterate(), 29, "iterator object iterates to 29");
  is( $it->iterate(), 31, "iterator object iterates to 31");
  $it->rewind->next->next->next->prev;
  is( $it->value(), 5, "iterator object rewind and move returns 5");
  # Validate that it automatically handles bigint range traversal.
  SKIP: {
    skip "Skipping bigint traversals on a Perl that can't add correctly",5 if $broken64;
    my $top_prime = prev_prime(~0);
    my $big_prime = next_prime(Math::BigInt->new(''.~0));
    ok( $big_prime > ~0, "internal check, next_prime on big int works");
    $it->rewind($top_prime);
    is( $it->value(), $top_prime, "iterator object can rewind to $top_prime");
    $it->next;
    is( $it->value(), $big_prime, "iterator object next is $big_prime");
    $it->rewind(~0);
    is( $it->value(), $big_prime, "iterator object rewound to ~0 is $big_prime");
    $it->prev;
    is( $it->value(), $top_prime, "iterator object prev goes back to $top_prime");
  }

  # Validation for the Math::NumSeq compatiblity stuff
  $it->rewind;
  do { $it->next } for 1..200;
  is( $it->tell_i(), 201, "iterator object tell_i");
  is( $it->i_start, 1, "iterator object i_start = 1");
  like( $it->description, qr/prime numbers/, "iterator object description");
  is( $it->values_min, 2, "iterator object values_min = 2");
  is( $it->values_max, undef, "iterator object values_max = undef");
  # missing: characteristic
  is( $it->oeis_anum, "A000040", "iterator object oeis_anum = A000040");
  # missing: parameter_info_array / parameter_info_list
  is( $it->seek_to_i(156)->value, 911, "iterator object seek_to_i goes to nth prime");
  is( $it->seek_to_value(156)->value, 157, "iterator object seek_to_value goes to value");
  is( $it->ith(589), 4289, "iterator object ith returns nth prime");
  ok( $it->pred(577), "iterator object pred returns true if is_prime");
  is( $it->value_to_i(4289), 589, "iterator object value_to_i works");
  is( $it->value_to_i(4290), undef, "iterator object value_to_i for non-prime returns undef");
  is( $it->value_to_i_floor(4290), 589, "iterator object value_to_i_floor");
  is( $it->value_to_i_ceil(4290), 590, "iterator object value_to_i_ceil");
  my $est = $it->value_to_i_estimate( 4171510507 );
  my $act = 197710788;
  # We will get an estimate that is much, much closer than Math::NumSeq
  ok( ($est > ($act-500)) && ($est < ($act+500)),
      "iterator object value_to_i_estimage is in range");
}

{ my @zn;
  forprimes {
    my $p=$_;
    forprimes {
      lastfor, push @zn,$_ if $_ % $p == 1;
    } 1000;
  } 100;
  is_deeply( \@zn,
             [3,7,11,29,23,53,103,191,47,59,311,149,83,173,283,107,709,367,269,569,293,317,167,179,389],
             "lastfor works in forprimes" );
}
{ my @zn;
  forprimes {
    my $p=$_;
    forcomposites {
      lastfor, push @zn,$_ if $_ % $p == 1;
    } 1000;
  } 100;
  is_deeply( \@zn,
             [9,4,6,8,12,14,18,20,24,30,32,38,42,44,48,54,60,62,68,72,74,80,84,90,98],
             "lastfor works in forcomposites" );
}
{ my @zn;
  forprimes {
    my $p=$_;
    foroddcomposites {
      lastfor, push @zn,$_ if $_ % $p == 1;
    } 1000;
  } 100;
  is_deeply( \@zn,
             [9,25,21,15,45,27,35,39,93,117,63,75,165,87,95,213,119,123,135,143,147,159,333,357,195],
             "lastfor works in foroddcomposites" );
}
{ my @powers;
  for my $n (1..20) {
    fordivisors { lastfor,push @powers,$_ if is_power($_) } $n;
  }
  is_deeply( \@powers, [4,4,9,4,4,9,4], "lastfor works in fordivisors" );
}
{ my $firstpart;
  forpart { lastfor,return if @_ < 4; $firstpart++; } 7;
  is($firstpart, 6, "lastfor works in forpart");
}
{ my $firstcomp;
  forcomp { lastfor,return if @_ < 4; $firstcomp++; } 7;
  is($firstcomp, 15, "lastfor works in forcomp");
}
{ my $smallcomb;
  forcomb { lastfor,return if vecsum(@_) > 11; $smallcomb++; } 7,4;
  is($smallcomb, 9, "lastfor works in forcomb");
}
{ my $t;
  forperm { lastfor,return if $_[3]==5; $t++; } 7;
  is($t, 12, "lastfor works in forperm");
}
{ my $t;
  forderange { lastfor,return if $_[3]==5; $t++; } 7;
  is($t, 5, "lastfor works in forderange");
}
{ my $t;
  formultiperm { lastfor if "miles" eq join("",@_); $t++; } [split(//,"smile")];
  is($t, 81, "lastfor works in formultiperm");
}
{ my @ps;
  forprimes {
    lastfor if $_ >= 7;
    # Note we keep going, unlike "last".
    push @ps, $_;
    forcomposites {
      push @ps,$_;
    } $_;
    # Our lastfor indicator is separate from the inside loop.
  } 20;
  is_deeply( \@ps, [2,3,5,4,7,4,6], "nested lastfor semantics" );
}
{
  my $t;
  forcomposites { $t=$_; lastfor if $_ > 2000; } 20000;
  is($t, 2001, "lastfor in forcomposites stops appropriately");
}

sub a053462 {
  my($s,$n)=(0,10**$_[0]-1);
  forsquarefree { $s += int($n / ($_*$_)) * ((scalar(@_) & 1)?-1:1); } sqrtint($n);
  $s;
}

{
  my $s;
  $s=0; forfactored { $s += $_ } 1; is($s, 1, "forfactored {} 1");
  $s=0; forfactored { $s += vecsum($_,@_) } 100; is($s, 7330, "forfactored {} 100");
  $s=0; forsquarefree { $s += vecsum($_,@_) } 100; is($s, 4763, "forsquarefree {} 100");
  $s=0; forfactored { $s += vecsum($_,@_) } 1e8,1e8+10; is($s, 1208835222, "forfactored {} 10^8,10^8+10");
  is( a053462(6), 607926, "A053462 using forsquarefree");
}

{
  my @got;
  forsemiprimes { push @got, $_; } 1000;
  is_deeply(\@got, [grep { is_semiprime($_) } 0 .. 1000], "forsemiprimes 1000");
}
