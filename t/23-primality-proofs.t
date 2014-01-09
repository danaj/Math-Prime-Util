#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_prime is_provable_prime is_provable_prime_with_cert
                         prime_certificate verify_prime
                         prime_get_config prime_set_config
                        /;
use Math::BigInt try => 'GMP';

my $use_test_warn;
BEGIN {
  eval "use Test::Warn";
  $use_test_warn = $@ ? 0 : 1;
}

my $extra = 0+(defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING});
my $use64 = ~0 > 4294967295;
my $broken64 = (18446744073709550592 == ~0);
# Do some tests only if:
#   EXTENDED_TESTING is on OR we have the GMP backend
# Note that with Calc, these things are incredibly slow.
my $doexpensive = 0 + ($extra || ( (!$use64 || !$broken64) && Math::BigInt->config()->{lib} eq 'Math::BigInt::GMP' ));

my @plist = qw/20907001 809120722675364249/;
if ($extra || $use64) {
  push @plist, "677826928624294778921";
}
if ($extra || prime_get_config->{'gmp'}) {
  push @plist, "980098182126316404630169387";
}

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
                                   #  hand-done proofs
            + 1*$doexpensive       #  n-1 for 2^521-1
            + 1*$extra             #  n-1 for 2^607-1
            #+ (($doexpensive && !$broken64) ? 1 : 0)  # n-1 proof
            + (($doexpensive) ? 1 : 0)  # n-1 proof
            + 2                    #  Pratt and ECPP
            + 28 # borked up certificates generate warnings
            + 6  # verification failures (tiny/BPSW)
            + 8  # verification failures (Lucas/Pratt)
            + 8  # verification failures (n-1)
            + 7  # verification failures (ECPP)
            + 3  # Verious other types
            + 0;

is( is_provable_prime(871139809), 0, "871139809 is composite" );
is( is_provable_prime(1490266103), 2, "1490266103 is provably prime" );

foreach my $p (@plist) {
 
  SKIP: {
    skip "Broken 64-bit causes trial factor to barf", 6
      if $broken64 && $p > 2**60;
    ok( is_prime($p), "$p is prime" );
    my($isp, $cert) = is_provable_prime_with_cert($p);
    is( $isp, 2, "   is_provable_prime_with_cert returns 2" );
    ok( defined($cert) && $cert =~ /^Type/m,
        "   certificate is non-null" );
    prime_set_config(verbose=>1);
    ok( verify_prime($cert), "   verification of certificate for $p done" );
    prime_set_config(verbose=>0);
    # Note, in some cases the certs could be non-equal (but both must be valid!)
    my $cert2 = prime_certificate($p);
    ok( defined($cert2) && $cert2 =~ /^Type/m,
        "   prime_certificate is also non-null" );
    if ($cert2 eq $cert) {
      ok(1, "   certificate is identical to first");
    } else {
      ok( verify_prime($cert2), "   different cert, verified" );
    }
  }
}

# Some hand-done proofs
if ($doexpensive) {
  my $proof = <<EOPROOF;
[MPU - Primality Certificate]
Version 1.0

Proof for:
N 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151

Type BLS5
N  6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151
Q[1]  65427463921
Q[2]  308761441
Q[3]  7623851
Q[4]  409891
Q[5]  61681
Q[6]  1613
Q[7]  521
Q[8]  157
Q[9]  131
Q[10]  53
Q[11]  41
Q[12]  31
Q[13]  17
A[0]  3
A[1]  3
A[2]  3
A[3]  3
A[4]  3
A[5]  3
A[6]  3
A[8]  3
A[9]  3
A[10]  3
A[11]  3
A[12]  3
A[13]  3
----
EOPROOF
  ok( verify_prime($proof), "2**521-1 primality proof verified" );
}
if ($extra) {
  my $proof = <<EOPROOF;
[MPU - Primality Certificate]
Version 1.0

Proof for:
N 531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127

Type BLS5
N  531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127
Q[1]  845100400152152934331135470251
Q[2]  341117531003194129
Q[3]  7432339208719
Q[4]  607
A[0]  3
A[1]  3
A[2]  3
A[3]  3
----

Type BLS5
N  845100400152152934331135470251
Q[1]  1801
Q[2]  601
Q[3]  251
Q[4]  101
A[1]  3
A[2]  3
A[3]  3
----
EOPROOF
  ok( verify_prime($proof), "2**607-1 primality proof verified" );
}
if ($doexpensive) {
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
  my @proof = (20907001, "Pratt", [ 2,
                                 3,
                                 5,
                                 [23,"Pratt",[2,[11,"Pratt",[2,5],2]],5],
                                 [101, "Pratt", [2, 5], 2],
                               ], 14 );
  ok( verify_prime(@proof), "simple Lucas/Pratt proof verified" );
}
{
  my $proof = <<EOPROOF;
[MPU - Primality Certificate]
Version 1.0

Proof for:
N 1030291136596639351761062717

Type ECPP
N  1030291136596639351761062717
A  1030291136596639351761062709
B  0
M  1030291136596575744618987466
Q  317851433704525489
X  4215121326
Y  246323849244309081587435955
EOPROOF

  ok( verify_prime($proof), "ECPP primality proof of 1030291136596639351761062717 verified" );
}
#{
#  my @proof = ('809120722675364249', "n-1", ["B", 20000, '22233477760919', 2], [2, 4549], [3, 2]);
#  ok( verify_prime(@proof), "n-1 T7 primality proof of 809120722675364249 verified" );
#}

# Failures for verify_prime

# First, let's get the borked up formats, which is should warn about.
SKIP: {
  skip "No Test::Warn", 28 unless $use_test_warn;
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
#is( verify_prime(['9848131514359', 'n-1', ["B", 20000, 890588851, 2], [2, 3, 19, 97], [3, 5, 2, 2]]), 1, "n-1 T7 proper" );
#is( verify_prime(['9848131514359', 'n-1', ["B", 20000, 890588951, 2], [2, 3, 19, 97], [3, 5, 2, 2]]), 0, "n-1 T7 with misfactor" );
#is( verify_prime(['9848131514359', 'n-1', ["B", 0, 890588851, 2], [2, 3, 19, 97], [3, 5, 2, 2]]), 0, "n-1 T7 with B < 1" );
#is( verify_prime(['9848131514359', 'n-1', ["B", 20000, 16921188169, 2], [2, 3, 97], [3, 5, 2]]), 0, "n-1 T7 with wrong B" );
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

{
  my $header = "[MPU - Primality Certificate]\nVersion 1.0\nProof for:";
  {
  my $cert = join "\n", $header,
                       "N 2297612322987260054928384863",
                       "Type Pocklington",
                       "N  2297612322987260054928384863",
                       "Q  16501461106821092981",
                       "A  5";
  is( verify_prime($cert), 1, "Verify Pocklington");
  }
  {
  my $cert = join "\n", $header,
                       "N 5659942549665396263282978117",
                       "Type BLS15",
                       "N  5659942549665396263282978117",
                       "Q  42941814754495493",
                       "LP 2",
                       "LQ 3";
  is( verify_prime($cert), 1, "Verify BLS15");
  }
  {
  my $cert = join "\n", $header,
                       "N 43055019307158602560279",
                       "Type ECPP3",
                       "N 43055019307158602560279",
                       "S 106563369",
                       "R 404032076977387",
                       "A 0",
                       "B 4",
                       "T 1";
  is( verify_prime($cert), 1, "Verify ECPP3");
  }
}
