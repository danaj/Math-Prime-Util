#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Test::Warn;
use Math::Prime::Util qw/is_prime is_provable_prime is_provable_prime_with_cert
                         prime_certificate verify_prime
                         prime_get_config
                        /;
use Math::BigInt try => 'GMP';

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $broken64 = (18446744073709550592 == ~0);

my @plist = qw/20907001 809120722675364249 677826928624294778921
               980098182126316404630169387/;

## This is too slow without Math::Prime::Util::GMP.
#push @plist, '3364125245431456304736426076174232972735419017865223025179282077503701'
#             if prime_get_config->{'gmp'};
#
#push @plist, '6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151'
#             if $extra
#             && (prime_get_config->{'gmp'} || Math::BigInt->config()->{lib} eq 'Math::BigInt::GMP');
#
#push @plist, '531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127'
#             if $extra && prime_get_config->{'gmp'};

@plist = sort { $a<=>$b }
         map { $_ > ~0 ? Math::BigInt->new("$_") : $_ }
         @plist;

plan tests => 0
            + 2  # is_provable_prime
            + 6 * scalar(@plist)
            + 6  # hand-done proofs
            + 28 # borked up certificates generating warnings
            + 6  # verification failures (tiny/BPSW)
            + 8  # verification failures (Lucas/Pratt)
            + 12 # verification failures (n-1)
            + 7  # verification failures (ECPP)
            + 0;

is( is_provable_prime(871139809), 0, "871139809 is composite" );
is( is_provable_prime(1490266103), 2, "1490266103 is provably prime" );

foreach my $p (@plist) {
 
  ok( is_prime($p), "$p is prime" );
  SKIP: {
    skip "Broken 64-bit causes trial factor to barf", 5
      if $broken64 && $p > 2**48;
    my($isp, $cert_ref) = is_provable_prime_with_cert($p);
    is( $isp, 2, "   is_provable_prime_with_cert returns 2" );
    ok( defined($cert_ref) && ref($cert_ref) eq 'ARRAY' && scalar(@$cert_ref) >= 1,
        "   certificate is non-null" );
    ok( verify_prime($cert_ref), "   verification of certificate reference done" );
    # Note, in some cases the two certs could be non-equal (but both must be valid!)
    my @cert = prime_certificate($p);
    ok( scalar(@cert) >= 1, "   prime_certificate is also non-null" );
    # TODO: compare certificates and skip if equal
    ok( verify_prime(@cert), "   verification of prime_certificate done" );
  }
}

# Some hand-done proofs
SKIP: {
  skip "Skipping 2**521-1 verification without Math::BigInt::GMP", 1
       unless Math::BigInt->config()->{lib} eq 'Math::BigInt::GMP';
  my @proof = ('6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151', 'n-1',
  [ 2,3,5,11,17,31,41,53,131,157,521,1613,61681,8191,42641,858001,51481, '7623851', '308761441' ],
  [ 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3 ] );
  ok( verify_prime(@proof), "2**521-1 primality proof verified" );
}
SKIP: {
  skip "Skipping 2**607-1 verification without Math::BigInt::GMP", 1
       unless Math::BigInt->config()->{lib} eq 'Math::BigInt::GMP';
  my @proof = ('531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127', 'n-1',
  [ 2,3,7,607,'112102729', '341117531003194129', '7432339208719',
    ['845100400152152934331135470251', 'n-1',
      [2,5,11,31,41,101,251,601,1801],
      [2,3,3,3,3,2,3,3,3]
    ] ],
  [ 3,5,3,2,3,3,3,3 ] );
  ok( verify_prime(@proof), "2**607-1 primality proof verified" );
}
{
  my @proof = ('809120722675364249', "n-1", ["B", 20000, '22233477760919', 2], [2, 4549], [3, 2]);
  ok( verify_prime(@proof), "n-1 T7 primality proof of 809120722675364249 verified" );
}
{
  my @proof = (20907001, "Pratt", [ 2,
                                 3,
                                 5,
                                 [23,"Pratt",[2,[11,"Pratt",[2,5],2]],5],
                                 [101, "Pratt", [2, 5], 2],
                               ], 14 );
  ok( verify_prime(@proof), "simple Lucas/Pratt proof verified" );
}
{
  my @proof = ('3364125245431456304736426076174232972735419017865223025179282077503701', 'n-1',
    [2,5,127, ['28432789963853652887491983185920687231739655787', 'n-1',
                [ 2,3,163,650933, [ '44662634059309451871488121651101494489', 'n-1',
                                    [ 2,3,23,4021,2321273 ],
                                    [ 11, 2, 2, 2, 2 ]
                                  ] ],
                [ 2, 2, 2, 2, 2 ]
              ],
              '9316417838190714313' ],
    [ 2, 2, 2, 2, 2 ]);
  ok( verify_prime(@proof), "simple n-1 proof verified" );
}
{
  my @proof = ('677826928624294778921',"AGKM", ['677826928624294778921', '404277700094248015180', '599134911995823048257', '677826928656744857936', '104088901820753203', ['2293544533', '356794037129589115041']], ['104088901820753203', '0', '73704321689372825', '104088902465395836', '1112795797', ['3380482019', '53320146243107032']], ['1112795797', '0', '638297481', '1112860899', '39019', ['166385704', '356512285']]);
  ok( verify_prime(@proof), "ECPP primality proof of 677826928624294778921 verified" );
}

# Failures for verify_prime

# First, let's get the borked up formats, which is should warn about.
{
  my $result;
  warning_like { $result = verify_prime([1490266103, 'INVALID', 1, 2, 3]) }
               { carped => qr/^verify_prime: / },
               "warning for unknown method";
  is( $result, 0, "   ...and returns 0" );

  warning_like { $result = verify_prime([1490266103, 'Pratt', 1, 2, 3]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid Lucas/Pratt";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'Pratt', 1, [2], 3]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid Lucas/Pratt";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'Pratt', [1], 2, 3]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid Lucas/Pratt";
  is( $result, 0, "   ...and returns 0" );

  warning_like { $result = verify_prime([1490266103, 'n-1', 1, 2, 3]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid n-1 (too many arguments)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'n-1', 1, 2]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid n-1 (non-array f,a)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'n-1', [1], 2]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid n-1 (non-array a)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'n-1', [2, 13, 19, 1597, 1889], [2, 2, 2, 2]]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid n-1 (too few a values)";
  is( $result, 0, "   ...and returns 0" );

  warning_like { $result = verify_prime([1490266103, 'ECPP']) }
               { carped => qr/^verify_prime: / },
               "warning for invalid ECPP (no n-certs)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'ECPP', 15]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid ECPP (non-array block)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'ECPP', [15,16,17]]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid ECPP (wrong size block)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'ECPP', [694361, 694358, 0, 695162, 26737, [348008, 638945]]]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid ECPP (block n != q)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'ECPP', [1490266103, 1442956066, 1025050760, 1490277784, 2780369, 531078754]]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid ECPP (block point wrong format)";
  is( $result, 0, "   ...and returns 0" );
  warning_like { $result = verify_prime([1490266103, 'ECPP', [1490266103, 1442956066, 1025050760, 1490277784, 2780369, [531078754, 0, 195830554]]]) }
               { carped => qr/^verify_prime: / },
               "warning for invalid ECPP (block point wrong format)";
  is( $result, 0, "   ...and returns 0" );
}

is( verify_prime([]), 0, "verify null is composite" );
is( verify_prime([2]), 1, "verify [2] is prime" );
is( verify_prime([9]), 0, "verify [9] is composite" );
is( verify_prime([14]), 0, "verify [14] is composite" );
is( verify_prime(['28446744073709551615']), 0, "verify BPSW with n > 2^64 fails" );
is( verify_prime([871139809]), 0, "verify BPSW with composite fails" );

is( verify_prime([1490266103, 'Pratt', [2,13,19,1597,1889], 5]), 1, "Lucas/Pratt proper" );
is( verify_prime([1490266103, 'Pratt', [4,13,19,1597,1889], 5]), 0, "Pratt with non-prime factors" );
is( verify_prime([1490266103, 'Pratt', [[4],13,19,1597,1889], 5]), 0, "Pratt with non-prime factors" );
is( verify_prime([1490266103, 'Pratt', [2,13,29,1597,1889], 5]), 0, "Pratt with wrong factors" );
is( verify_prime([1490266103, 'Pratt', [2,13,19,1597], 5]), 0, "Pratt with not enough factors" );
is( verify_prime([1490266103, 'Pratt', [2,13,19,1597,1889], 1490266103]), 0, "Pratt with coprime a" );
is( verify_prime([185156263, 'Pratt', [2,3,3,10286459], 2]), 0, "Pratt with non-psp a" );
is( verify_prime([1490266103, 'Pratt', [2,13,19,1597,1889], 3]), 0, "Pratt with a not valid for all f" );

is( verify_prime([1490266103, 'n-1', [2, 13, 19, 1597, 1889], [5, 2, 2, 2, 2]]), 1, "n-1 proper" );
is( verify_prime([1490266103, 'n-1', [2, 23, 19, 1597, 1889], [5, 2, 2, 2, 2]]), 0, "n-1 with wrong factors" );
is( verify_prime([1490266103, 'n-1', [13, 19, 1597, 1889], [2, 2, 2, 2]]), 0, "n-1 without 2 as a factor" );
is( verify_prime([1490266103, 'n-1', [2, 13, 1889, 30343], [5, 2, 2, 2]]), 0, "n-1 with a non-prime factor" );
is( verify_prime([1490266103, 'n-1', [2, 13, 1889, [30343]], [5, 2, 2, 2]]), 0, "n-1 with a non-prime array factor" );
# I don't know how to make F and R (A and B) to not be coprime
is( verify_prime(['9848131514359', 'n-1', ["B", 20000, 890588851, 2], [2, 3, 19, 97], [3, 5, 2, 2]]), 1, "n-1 T7 proper" );
is( verify_prime(['9848131514359', 'n-1', ["B", 20000, 890588951, 2], [2, 3, 19, 97], [3, 5, 2, 2]]), 0, "n-1 T7 with misfactor" );
is( verify_prime(['9848131514359', 'n-1', ["B", 0, 890588851, 2], [2, 3, 19, 97], [3, 5, 2, 2]]), 0, "n-1 T7 with B < 1" );
is( verify_prime(['9848131514359', 'n-1', ["B", 20000, 16921188169, 2], [2, 3, 97], [3, 5, 2]]), 0, "n-1 T7 with wrong B" );
is( verify_prime([1490266103, 'n-1', [2, 13], [5, 2]]), 0, "n-1 without enough factors" );
is( verify_prime([914144252447488195, 'n-1', [2, 3, 11, 17, 1531], [2, 2, 2, 2, 2]]), 0, "n-1 with bad BLS75 r/s" );
is( verify_prime([1490266103, 'n-1', [2, 13, 19, 1597, 1889], [3, 2, 2, 2, 2]]), 0, "n-1 with bad a value" );

is( verify_prime([1490266103, "ECPP",
                 [1490266103, 1442956066, 1025050760, 1490277784, 2780369, [531078754, 195830554]],
                 [2780369, 2780360, 0, 2777444, 694361, [2481811, 1317449]],
                 [694361, 694358, 0, 695162, 26737, [348008, 638945]]]),
                 1, "ECPP proper" );
is( verify_prime([1490266103, "ECPP",
                 [1490266103, 1442956066, 1025050760, 1490277784, 5560738, [531078754, 195830554]],
                 [5560738, 2780360, 0, 2777444, 694361, [2481811, 1317449]]]),
                 0, "ECPP q is divisible by 2" );
is( verify_prime([74468183, "ECPP",
                 [74468183, 89, 1629, 74475075, 993001, [47943960, 8832604]],
                 [993001, 0, 992984, 994825, 3061, [407531, 231114]]]),
                 0, "ECPP a/b invalid" );
is( verify_prime([1490266103, "ECPP",
                 [1490266103, 1442956066, 1025050760, 1490277784, 536, [531078754, 195830554]],
                 [536, 2780360, 0, 2777444, 694361, [2481811, 1317449]]]),
                 0, "ECPP q is too small" );
is( verify_prime([694361, "ECPP",
                 [694361, 694358, 0, 30, 26737, [264399, 59977]]]),
                 0, "ECPP multiplication wrong (infinity)" );
is( verify_prime([694361, "ECPP",
                 [694361, 694358, 0, 695161, 26737, [264399, 59977]]]),
                 0, "ECPP multiplication wrong (not infinity)" );
is( verify_prime([1490266103, "ECPP",
                 [1490266103, 1442956066, 1025050760, 1490277784, 2780369, [531078754, 195830554]],
                 [2780369, 2780360, 0, 2777444, 694361, [2481811, 1317449]],
                 [694361, 694358, 0, 695162, [26737, "n-1", [2],[2]], [348008, 638945]]]),
                 0, "ECPP non-prime last q" );
