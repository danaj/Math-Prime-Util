requires 'ExtUtils::MakeMaker';


requires 'Exporter', '5.562';
requires 'XSLoader', '0.01';
requires 'Carp';
requires 'Tie::Array';
requires 'base';
requires 'Config';
requires 'Math::BigInt', '1.88';
requires 'Math::BigFloat', '1.59';

requires 'Bytes::Random::Secure', '0.23';

recommends 'Math::Prime::Util::GMP', '0.15';
recommends 'Math::BigInt::GMP';
recommends 'Math::MPFR', '2.03';


on test => sub {
  requires 'Test::More', '0.45';
  requires 'bignum', '0.22';
  recommends 'Test::Warn';
};
