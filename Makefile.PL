use ExtUtils::MakeMaker;

WriteMakefile1(
    NAME         => 'Math::Prime::Util',
    ABSTRACT     => 'Utilities related to prime numbers, including fast sieves and factoring',
    VERSION_FROM => 'lib/Math/Prime/Util.pm',
    LICENSE      => 'perl',
    AUTHOR       => 'Dana A Jacobsen <dana@acm.org>',

    OBJECT       => 'cache.o '  .
                    'factor.o ' .
                    'aks.o '    .
                    'lehmer.o ' .
                    'sieve.o '  .
                    'util.o '   .
                    'XS.o',
    LIBS         => ['-lm'],

    EXE_FILES    => ['bin/primes.pl', 'bin/factor.pl'],

    BUILD_REQUIRES=>{
                      'Test::More'       => '0.45',
                      'bignum'           => '0.22',   # 'use bigint' in tests
                    },
    PREREQ_PM    => {
                      'Exporter'         => '5.562',
                      'XSLoader'         => '0.01',
                      'Carp'             => 0,
                      'Tie::Array'       => 0,
                      'base'             => 0,
                      'Config'           => 0,
                      'Math::BigInt'     => '1.88',
                      'Math::BigFloat'   => '1.59',
                      'Bytes::Random::Secure' => '0.23',
                    },
    META_MERGE   => {
               resources => {
                  homepage   => 'https://github.com/danaj/Math-Prime-Util',
                  repository => 'https://github.com/danaj/Math-Prime-Util',
               },
               recommends => {
                 'Math::Prime::Util::GMP' => 0.06,
                 'Math::BigInt::GMP'      => 0,
                 'Math::MPFR'             => 2.03,
               },
             },

    MIN_PERL_VERSION => 5.006002,
);

sub WriteMakefile1 {   # Cribbed from eumm-upgrade by Alexandr Ciornii
  my %params = @_;
  my $eumm_version = $ExtUtils::MakeMaker::VERSION;
  $eumm_version = eval $eumm_version;

  if ($params{BUILD_REQUIRES} and $eumm_version < 6.5503) {
      #EUMM 6.5502 has problems with BUILD_REQUIRES
      $params{PREREQ_PM}={ %{$params{PREREQ_PM} || {}} , %{$params{BUILD_REQUIRES}} };
      delete $params{BUILD_REQUIRES};
  }
  delete $params{CONFIGURE_REQUIRES} if $eumm_version < 6.52;
  delete $params{MIN_PERL_VERSION} if $eumm_version < 6.48;
  delete $params{META_MERGE} if $eumm_version < 6.46;
  delete $params{META_ADD} if $eumm_version < 6.46;
  delete $params{LICENSE} if $eumm_version < 6.31;
  delete $params{AUTHOR} if $] < 5.005;
  delete $params{ABSTRACT_FROM} if $] < 5.005;
  delete $params{BINARY_LOCATION} if $] < 5.005;

  WriteMakefile(%params);
}
