#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_prime is_provable_prime is_provable_prime_with_cert
                         prime_certificate verify_prime
                         prime_get_config
                        /;
use Math::BigInt try => 'GMP';

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

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
            + 2 # is_provable_prime
            + 6 * scalar(@plist)
            + 6 # hand-done proofs
            + 17 # verification failures
            + 0;

is( is_provable_prime(871139809), 0, "871139809 is composite" );
is( is_provable_prime(1490266103), 2, "1490266103 is provably prime" );

foreach my $p (@plist) {
 
  ok( is_prime($p), "$p is prime" );
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
is( verify_prime([]), 0, "verify null is composite" );
is( verify_prime([2]), 1, "verify [2] is prime" );
is( verify_prime([9]), 0, "verify [9] is composite" );
is( verify_prime([14]), 0, "verify [14] is composite" );
is( verify_prime(['28446744073709551615']), 0, "verify BPSW with n > 2^64 fails" );
is( verify_prime([871139809]), 0, "verify BPSW with composite fails" );
is( verify_prime([1490266103, 'INVALID', 1, 2, 3]), 0, "unknown method" );
is( verify_prime([1490266103, 'Pratt', 1, 2, 3]), 0, "Pratt with wrong count" );
is( verify_prime([1490266103, 'Pratt', 1, [2]]), 0, "Pratt with non-array arguments" );
is( verify_prime([1490266103, 'Pratt', [1], [2]]), 0, "Pratt with non-array arguments" );
is( verify_prime([1490266103, 'Pratt', [4,13,19,1597,1889], 5]), 0, "Pratt with non-prime factors" );
is( verify_prime([1490266103, 'Pratt', [[4],13,19,1597,1889], 5]), 0, "Pratt with non-prime factors" );
is( verify_prime([1490266103, 'Pratt', [2,13,29,1597,1889], 5]), 0, "Pratt with wrong factors" );
is( verify_prime([1490266103, 'Pratt', [2,13,19,1597], 5]), 0, "Pratt with not enough factors" );
is( verify_prime([1490266103, 'Pratt', [2,13,19,1597,1889], 1490266103]), 0, "Pratt with coprime a" );
is( verify_prime([185156263, 'Pratt', [2,3,3,10286459], 2]), 0, "Pratt with non-psp a" );
is( verify_prime([1490266103, 'Pratt', [2,13,19,1597,1889], 3]), 0, "Pratt with a not valid for all f" );
