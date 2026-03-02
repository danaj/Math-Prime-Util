requires 'ExtUtils::MakeMaker';

requires 'Exporter', '5.57';
requires 'XSLoader', '0.01';
requires 'Carp';
requires 'Tie::Array';
requires 'base';
requires 'constant';
requires 'Config';
requires 'Math::BigInt', '1.999814';
requires 'Math::BigFloat', '1.59';

recommends 'Math::Prime::Util::GMP', '0.53';
recommends 'Math::BigInt::GMP';
recommends 'Math::GMPz', '0.68';
recommends 'Digest::SHA', '5.87';

on test => sub {
  requires 'Test::More', '0.96';
  requires 'bignum', '0.65';
  recommends 'Test::Warn';
};
