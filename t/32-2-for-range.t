#!/usr/bin/env perl
use strict;
use warnings;

# On 64-bit, without GMP or a good Math::BigInt backend, bigint range tests for
# forprimes, forcomposites, and foroddcomposites are relatively slow.

use Test::More;
use Math::Prime::Util qw/forprimes forcomposites foroddcomposites
                         forsemiprimes foralmostprimes
                         is_semiprime is_almost_prime divisors/;

our $iter_scope_probe;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

subtest 'forprimes invalid arguments' => sub {
  ok(!eval { forprimes { 1 } undef; },       "forprimes undef");
  ok(!eval { forprimes { 1 } 2, undef; },    "forprimes 2,undef");
  ok(!eval { forprimes { 1 } undef, 2; },    "forprimes undef,2");
  ok(!eval { forprimes { 1 } -2, 3; },       "forprimes -2,3");
  ok(!eval { forprimes { 1 } 2, -3; },       "forprimes 2,-3");
  ok(!eval { forprimes { 1 } "abc"; },       "forprimes abc");
  ok(!eval { forprimes { 1 } 2, "abc"; },    "forprimes 2,abc");
  ok(!eval { forprimes { 1 } 5.6; },         "forprimes 5.6");
};

subtest 'forprimes ranges' => sub {
  my @T1 = (
    [1, []],
    [2, [2]],
    [3, [2,3]],
    [4, [2,3]],
    [5, [2,3,5]],
    [50, [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]],
  );
  for my $i (@T1) {
    my($n,$exp) = @$i;
    my @got;
    forprimes { push @got, $_; } $n;
    is_deeply(\@got, $exp, "forprimes $n");
  }

  my @T2 = (
    [  0,  0, []],
    [  0,  1, []],
    [  3,  5, [3,5]],
    [  3,  6, [3,5]],
    [  3,  7, [3,5,7]],
    [  5,  7, [5,7]],
    [  6,  7, [7]],
    [  5, 11, [5,7,11]],
    [  7, 11, [7,11]],
    [  2, 20, [2,3,5,7,11,13,17,19]],
    [ 20, 30, [23,29]],
    [199,223, [199,211,223]],
    [31398, 31468, []],
    [2147483647,2147483659, [2147483647,2147483659]],
    [3842610774,3842611326, [3842611109,3842611139,3842611163,3842611181,3842611211,3842611229,3842611249,3842611259,3842611261,3842611291,3842611301]],
  );
  for my $i (@T2) {
    my($lo,$hi,$exp) = @$i;
    my @got;
    forprimes { push @got, $_; } $lo,$hi;
    is_deeply(\@got, $exp, "forprimes $lo,$hi");
  }

  my $sum = 0;
  forprimes { $sum += int(12345/$_) } 1000;
  is(27053, $sum, "forprimes handles \$_ type changes");

  my @t;
  forprimes {
    forprimes {
      forprimes { push @t, $_ } $_,$_+10;
    } 10*$_,10*$_+10;
  } 10;
  is_deeply( \@t, [qw/23 29 31 29 31 37 31 37 41 37 41 43 47 53 59 61 59 61 67 71 73 79 73 79 83 79 83 89/], "triple nested forprimes" );
};

subtest 'forprimes callback scope' => sub {
  my @before;
  local $iter_scope_probe = "seed";
  forprimes {
    push @before, $iter_scope_probe;
    local $iter_scope_probe = $_;
  } 2,11;
  is_deeply(\@before, [("seed") x 5], "forprimes local scalar restored per callback");

  my @refs;
  forprimes {
    my $x = $_;
    push @refs, \$x;
  } 2,11;
  is_deeply([map { $$_ } @refs], [2,3,5,7,11], "forprimes lexical refs preserve iteration values");

  @refs = ();
  forprimes {
    my $x = $_;
    push @refs, \$x;
  } 2,7;
  ${$refs[0]} = 42;
  is_deeply([map { $$_ } @refs], [42,3,5,7], "forprimes lexical refs are distinct scalars");
};

subtest 'forprimes nested calls' => sub {
  my $sum = 0;
  forprimes {
    my @d = divisors(6486480 * ($_-1));
    $sum += 1+$#d;
  } 2,10;
  is($sum, 2016, "Nested call to large divisors inside forprimes");
};

subtest 'forcomposites' => sub {
  my @t;

  @t=(); forcomposites { push @t, $_ } 2147483647,2147483659;
  is_deeply( \@t, [qw/2147483648 2147483649 2147483650 2147483651 2147483652 2147483653 2147483654 2147483655 2147483656 2147483657 2147483658/], "forcomposites 2147483647,2147483659" );

  @t=(); forcomposites { push @t, $_ } 50;
  is_deeply( \@t, [qw/4 6 8 9 10 12 14 15 16 18 20 21 22 24 25 26 27 28 30 32 33 34 35 36 38 39 40 42 44 45 46 48 49 50/], "forcomposites 50" );

  @t=(); forcomposites { push @t, $_ } 200,410;
  is_deeply( \@t, [qw/200 201 202 203 204 205 206 207 208 209 210 212 213 214 215 216 217 218 219 220 221 222 224 225 226 228 230 231 232 234 235 236 237 238 240 242 243 244 245 246 247 248 249 250 252 253 254 255 256 258 259 260 261 262 264 265 266 267 268 270 272 273 274 275 276 278 279 280 282 284 285 286 287 288 289 290 291 292 294 295 296 297 298 299 300 301 302 303 304 305 306 308 309 310 312 314 315 316 318 319 320 321 322 323 324 325 326 327 328 329 330 332 333 334 335 336 338 339 340 341 342 343 344 345 346 348 350 351 352 354 355 356 357 358 360 361 362 363 364 365 366 368 369 370 371 372 374 375 376 377 378 380 381 382 384 385 386 387 388 390 391 392 393 394 395 396 398 399 400 402 403 404 405 406 407 408 410/], "forcomposites 200,410" );
};

subtest 'forsemiprimes' => sub {
  my @got;
  forsemiprimes { push @got, $_; } 1000;
  is_deeply(\@got, [grep { is_semiprime($_) } 0 .. 1000], "forsemiprimes 1000");
};

subtest 'foralmostprimes' => sub {
  my $num = 0;
  foralmostprimes { $num++; } 0,1,1000;
  is($num, 0, "foralmostprimes 0,1000 is empty");

  for my $k (1 .. 10) {
    my @got;
    foralmostprimes { push @got, $_; } $k,1000;
    is_deeply(\@got, [grep { is_almost_prime($k,$_) } 0 .. 1000], "foralmostprimes $k,1000");
  }
};

subtest 'bigint ranges' => sub {
  my($E,%d);
  if ($use64) {
    $E = 66;
    $d{primes} = [166,199,169];
    $d{semi} = [98,99,99];
    $d{almost} = [30,32,31];
    $d{comp} = [1506,1508,1506,1508];
    $d{oddcomp} = [1506,1511,1509,1511];
  } else {
    $E = 36;
    $d{primes} = [11,71,31];
    $d{semi} = [710,711,711];
    $d{almost} = [246,248,247];
    $d{comp} = [2166,2168,2166,2168];
    $d{oddcomp} = [1800,1805,1803,1805];
  }

  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{primes});
    forprimes { push @r,"$_" } $arg1, $arg2;
    is_deeply(\@r, \@res, "forprimes {} 2^$E+$a1, 2^$E+$a2");
  }
  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{semi});
    forsemiprimes { push @r,"$_" } $arg1, $arg2;
    is_deeply(\@r, \@res, "forsemiprimes {} 2^$E+$a1, 2^$E+$a2");
  }
  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{almost});
    foralmostprimes { push @r,"$_" } 3, $arg1, $arg2;
    is_deeply(\@r, \@res, "foralmostprimes {} 3, 2^$E+$a1, 2^$E+$a2");
  }
  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{comp});
    forcomposites { push @r,"$_" } $arg1, $arg2;
    is_deeply(\@r, \@res, "forcomposites {} 2^$E+$a1, 2^$E+$a2");
  }
  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{oddcomp});
    foroddcomposites { push @r,"$_" } $arg1, $arg2;
    is_deeply(\@r, \@res, "foroddcomposites {} 2^$E+$a1, 2^$E+$a2");
  }
};

sub split_d {
  my($E,$arr) = @_;
  my $a1 = $arr->[0];
  my $a2 = $arr->[1];
  return ($a1,$a2,map { ref($_) ? $_ : 2**$E + $_    } @$arr) if $E < 48;
  return ($a1,$a2,map { ref($_) ? $_ : plus_2_66($_) } @$arr) if $E == 66;
  die "unsupported test exponent $E";
}

sub plus_2_66 {
  my $add = shift;
  my $final = $add + 6464;
  die "too large" if $final > 99999;
  return "737869762948382" . sprintf("%05d",$final);
}

done_testing();
