package Math::Prime::Util::ChaCha;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::ChaCha::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::ChaCha::VERSION = '0.61';
}

###############################################################################
# Begin ChaCha core, reference RFC 7539
# with change to make blockcount/nonce be 64/64 from 32/96
# Dana Jacobsen, 9 Apr 2017

BEGIN {
  use constant ROUNDS => 20;
}

sub _quarterround {
  my($a,$b,$c,$d) = @_;
  $a=($a+$b)&0xFFFFFFFF;  $d^=$a;  $d=(($d << 16) | ($d >>(32-16)))&0xFFFFFFFF;
  $c=($c+$d)&0xFFFFFFFF;  $b^=$c;  $b=(($b << 12) | ($b >>(32-12)))&0xFFFFFFFF;
  $a=($a+$b)&0xFFFFFFFF;  $d^=$a;  $d=(($d <<  8) | ($d >>(32- 8)))&0xFFFFFFFF;
  $c=($c+$d)&0xFFFFFFFF;  $b^=$c;  $b=(($b <<  7) | ($b >>(32- 7)))&0xFFFFFFFF;
  ($a,$b,$c,$d);
}

sub _test_qr {
  return unless ROUNDS == 20;
  my($a,$b,$c,$d) = _quarterround(0x11111111,0x01020304,0x9b8d6f43,0x01234567);
  #printf "  %08x  %08x  %08x  %08x\n", $a,$b,$c,$d;
  die "QR test 2.1.1 fail 1" unless $a == 0xea2a92f4 && $b == 0xcb1cf8ce && $c == 0x4581472e && $d == 0x5881c4bb;

  my @x = map { hex } (qw/879531e0  c5ecf37d  516461b1  c9a62f8a
                          44c20ef3  3390af7f  d9fc690b  2a5f714c
                          53372767  b00a5631  974c541a  359e9963
                          5c971061  3d631689  2098d9d6  91dbd320/);
  my @e = @x;
  @e[ 2, 7, 8,13] = (0xbdb886dc, 0xcfacafd2, 0xe46bea80, 0xccc07c79);
  @x[ 2, 7, 8,13] = _quarterround(@x[ 2, 7, 8,13]);
  #printf "  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n", @x;
  for (0..15) {
    die "QR test 2.2.1 fail 2 at $_" unless $x[$_] == $e[$_];
  }
}
_test_qr;

#  State is:
#       cccccccc  cccccccc  cccccccc  cccccccc
#       kkkkkkkk  kkkkkkkk  kkkkkkkk  kkkkkkkk
#       kkkkkkkk  kkkkkkkk  kkkkkkkk  kkkkkkkk
#       bbbbbbbb  nnnnnnnn  nnnnnnnn  nnnnnnnn
#
#     c=constant k=key b=blockcount n=nonce

sub _core {
  my(@x) = @_;
  die "Invalid ChaCha state" unless scalar(@x) == 16;
  my @j = @x;
  for (1 .. ROUNDS/2) {

    #@x[ 0, 4, 8,12] = _quarterround(@x[ 0, 4, 8,12]);
    $x[ 0]=($x[ 0]+$x[ 4])&0xFFFFFFFF;  $x[12]^=$x[ 0];  $x[12]=(($x[12] << 16) | ($x[12] >> 16))&0xFFFFFFFF;
    $x[ 8]=($x[ 8]+$x[12])&0xFFFFFFFF;  $x[ 4]^=$x[ 8];  $x[ 4]=(($x[ 4] << 12) | ($x[ 4] >> 20))&0xFFFFFFFF;
    $x[ 0]=($x[ 0]+$x[ 4])&0xFFFFFFFF;  $x[12]^=$x[ 0];  $x[12]=(($x[12] <<  8) | ($x[12] >> 24))&0xFFFFFFFF;
    $x[ 8]=($x[ 8]+$x[12])&0xFFFFFFFF;  $x[ 4]^=$x[ 8];  $x[ 4]=(($x[ 4] <<  7) | ($x[ 4] >> 25))&0xFFFFFFFF;

    #@x[ 1, 5, 9,13] = _quarterround(@x[ 1, 5, 9,13]);
    $x[ 1]=($x[ 1]+$x[ 5])&0xFFFFFFFF;  $x[13]^=$x[ 1];  $x[13]=(($x[13] << 16) | ($x[13] >> 16))&0xFFFFFFFF;
    $x[ 9]=($x[ 9]+$x[13])&0xFFFFFFFF;  $x[ 5]^=$x[ 9];  $x[ 5]=(($x[ 5] << 12) | ($x[ 5] >> 20))&0xFFFFFFFF;
    $x[ 1]=($x[ 1]+$x[ 5])&0xFFFFFFFF;  $x[13]^=$x[ 1];  $x[13]=(($x[13] <<  8) | ($x[13] >> 24))&0xFFFFFFFF;
    $x[ 9]=($x[ 9]+$x[13])&0xFFFFFFFF;  $x[ 5]^=$x[ 9];  $x[ 5]=(($x[ 5] <<  7) | ($x[ 5] >> 25))&0xFFFFFFFF;

    #@x[ 2, 6,10,14] = _quarterround(@x[ 2, 6,10,14]);
    $x[ 2]=($x[ 2]+$x[ 6])&0xFFFFFFFF;  $x[14]^=$x[ 2];  $x[14]=(($x[14] << 16) | ($x[14] >> 16))&0xFFFFFFFF;
    $x[10]=($x[10]+$x[14])&0xFFFFFFFF;  $x[ 6]^=$x[10];  $x[ 6]=(($x[ 6] << 12) | ($x[ 6] >> 20))&0xFFFFFFFF;
    $x[ 2]=($x[ 2]+$x[ 6])&0xFFFFFFFF;  $x[14]^=$x[ 2];  $x[14]=(($x[14] <<  8) | ($x[14] >> 24))&0xFFFFFFFF;
    $x[10]=($x[10]+$x[14])&0xFFFFFFFF;  $x[ 6]^=$x[10];  $x[ 6]=(($x[ 6] <<  7) | ($x[ 6] >> 25))&0xFFFFFFFF;

    #@x[ 3, 7,11,15] = _quarterround(@x[ 3, 7,11,15]);
    $x[ 3]=($x[ 3]+$x[ 7])&0xFFFFFFFF;  $x[15]^=$x[ 3];  $x[15]=(($x[15] << 16) | ($x[15] >> 16))&0xFFFFFFFF;
    $x[11]=($x[11]+$x[15])&0xFFFFFFFF;  $x[ 7]^=$x[11];  $x[ 7]=(($x[ 7] << 12) | ($x[ 7] >> 20))&0xFFFFFFFF;
    $x[ 3]=($x[ 3]+$x[ 7])&0xFFFFFFFF;  $x[15]^=$x[ 3];  $x[15]=(($x[15] <<  8) | ($x[15] >> 24))&0xFFFFFFFF;
    $x[11]=($x[11]+$x[15])&0xFFFFFFFF;  $x[ 7]^=$x[11];  $x[ 7]=(($x[ 7] <<  7) | ($x[ 7] >> 25))&0xFFFFFFFF;

    #@x[ 0, 5,10,15] = _quarterround(@x[ 0, 5,10,15]);
    $x[ 0]=($x[ 0]+$x[ 5])&0xFFFFFFFF;  $x[15]^=$x[ 0];  $x[15]=(($x[15] << 16) | ($x[15] >> 16))&0xFFFFFFFF;
    $x[10]=($x[10]+$x[15])&0xFFFFFFFF;  $x[ 5]^=$x[10];  $x[ 5]=(($x[ 5] << 12) | ($x[ 5] >> 20))&0xFFFFFFFF;
    $x[ 0]=($x[ 0]+$x[ 5])&0xFFFFFFFF;  $x[15]^=$x[ 0];  $x[15]=(($x[15] <<  8) | ($x[15] >> 24))&0xFFFFFFFF;
    $x[10]=($x[10]+$x[15])&0xFFFFFFFF;  $x[ 5]^=$x[10];  $x[ 5]=(($x[ 5] <<  7) | ($x[ 5] >> 25))&0xFFFFFFFF;

    #@x[ 1, 6,11,12] = _quarterround(@x[ 1, 6,11,12]);
    $x[ 1]=($x[ 1]+$x[ 6])&0xFFFFFFFF;  $x[12]^=$x[ 1];  $x[12]=(($x[12] << 16) | ($x[12] >> 16))&0xFFFFFFFF;
    $x[11]=($x[11]+$x[12])&0xFFFFFFFF;  $x[ 6]^=$x[11];  $x[ 6]=(($x[ 6] << 12) | ($x[ 6] >> 20))&0xFFFFFFFF;
    $x[ 1]=($x[ 1]+$x[ 6])&0xFFFFFFFF;  $x[12]^=$x[ 1];  $x[12]=(($x[12] <<  8) | ($x[12] >> 24))&0xFFFFFFFF;
    $x[11]=($x[11]+$x[12])&0xFFFFFFFF;  $x[ 6]^=$x[11];  $x[ 6]=(($x[ 6] <<  7) | ($x[ 6] >> 25))&0xFFFFFFFF;

    #@x[ 2, 7, 8,13] = _quarterround(@x[ 2, 7, 8,13]);
    $x[ 2]=($x[ 2]+$x[ 7])&0xFFFFFFFF;  $x[13]^=$x[ 2];  $x[13]=(($x[13] << 16) | ($x[13] >> 16))&0xFFFFFFFF;
    $x[ 8]=($x[ 8]+$x[13])&0xFFFFFFFF;  $x[ 7]^=$x[ 8];  $x[ 7]=(($x[ 7] << 12) | ($x[ 7] >> 20))&0xFFFFFFFF;
    $x[ 2]=($x[ 2]+$x[ 7])&0xFFFFFFFF;  $x[13]^=$x[ 2];  $x[13]=(($x[13] <<  8) | ($x[13] >> 24))&0xFFFFFFFF;
    $x[ 8]=($x[ 8]+$x[13])&0xFFFFFFFF;  $x[ 7]^=$x[ 8];  $x[ 7]=(($x[ 7] <<  7) | ($x[ 7] >> 25))&0xFFFFFFFF;

    #@x[ 3, 4, 9,14] = _quarterround(@x[ 3, 4, 9,14]);
    $x[ 3]=($x[ 3]+$x[ 4])&0xFFFFFFFF;  $x[14]^=$x[ 3];  $x[14]=(($x[14] << 16) | ($x[14] >> 16))&0xFFFFFFFF;
    $x[ 9]=($x[ 9]+$x[14])&0xFFFFFFFF;  $x[ 4]^=$x[ 9];  $x[ 4]=(($x[ 4] << 12) | ($x[ 4] >> 20))&0xFFFFFFFF;
    $x[ 3]=($x[ 3]+$x[ 4])&0xFFFFFFFF;  $x[14]^=$x[ 3];  $x[14]=(($x[14] <<  8) | ($x[14] >> 24))&0xFFFFFFFF;
    $x[ 9]=($x[ 9]+$x[14])&0xFFFFFFFF;  $x[ 4]^=$x[ 9];  $x[ 4]=(($x[ 4] <<  7) | ($x[ 4] >> 25))&0xFFFFFFFF;
  }
  map { ($x[$_] + $j[$_]) & 0xFFFFFFFF } 0 .. 15;
}
sub _test_core {
  return unless ROUNDS == 20;
  my $init_state = '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001090000004a00000000000000';
  my @state = map { hex } unpack "(a8)*", $init_state;
    #my @c = (0x61707865, 0x3320646e, 0x79622d32, 0x6b206574);
    #my $key = join("", map { chr } 0 .. 31);
    #my @key = unpack("V8",$key);
    #my @nonce = unpack("V3",pack("N3",0x09,0x04a,0x00));
    #my $counter = 1;
    #my @state = (@c, @key, $counter, @nonce);
  my $instr = join("",map { sprintf("%08x",$_) } @state);
  die "Block function fail test 2.3.2 input" unless $instr eq '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001090000004a00000000000000';
  @state = _core(@state);
  my $outstr = join("",map { sprintf("%08x",$_) } @state);
  #printf "  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n", @state;
  #print "outstr: $outstr\n";
  die "Block function fail test 2.3.2 output" unless $outstr eq 'e4e7f11015593bd11fdd0f50c47120a3c7f4d1c70368c0339aaa22044e6cd4c3466482d209aa9f0705d7c214a2028bd9d19c12b5b94e16dee883d0cb4e3c50a2';
}
_test_core();

sub _keystream {
  my($nbytes, $rstate) = @_;
  croak "Keystream invalid state" unless scalar(@$rstate) == 16;
  my $stream = '';
  while ($nbytes > 0) {
    my @block = _core(@$rstate);
    $stream .= pack("L16",@block);
    #print "block $rstate->[12]: ",join(" ",map { sprintf("%08x",$_) } @block),"\n";
    $nbytes -= 64;
    $rstate->[12]++;
    if ($rstate->[12] > 4294967295) {
      $rstate->[12] = 0;  $rstate->[13]++;
    }
  }
  return $stream;
}
sub _test_keystream {
  return unless ROUNDS == 20;
  my $init_state = '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001000000004a00000000000000';
  my @state = map { hex } unpack "(a8)*", $init_state;
  my $instr = join("",map { sprintf("%08x",$_) } @state);
  croak "Block function fail test 2.4.2 input" unless $instr eq '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001000000004a00000000000000';
  my $keystream = _keystream(114, \@state);
  # Verify new state
  my $outstr = join("",map { sprintf("%08x",$_) } @state);
  #print "outstr: $outstr\n";
  croak "Block function fail test 2.4.2 output" unless $outstr eq '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000003000000004a00000000000000';
  # Rather tediaus way of doing it, but very explicit
  my $ksstr = join("",map { sprintf("%02x",$_) } map { ord } split(//,$keystream));
  croak "Block function fail test 2.4.2 keystream" unless substr($ksstr,0,2*114) eq '224f51f3401bd9e12fde276fb8631ded8c131f823d2c06e27e4fcaec9ef3cf788a3b0aa372600a92b57974cded2b9334794cba40c63e34cdea212c4cf07d41b769a6749f3f630f4122cafe28ec4dc47e26d4346d70b98c73f3e9c53ac40c5945398b6eda1a832c89c167eacd901d7e2bf363';
}
_test_keystream();

# End ChaCha core
###############################################################################

my $_goodseed;
my $_state;     # This should be far better hidden from other modules
my $_stream;
my $_have;
my $_sptr;

sub _is_csprng_well_seeded { $_goodseed }

sub seed_csprng {
  my($seed) = @_;
  $_goodseed = length($seed) >= 16;
  my(@seed) = unpack("L10",$seed);
  # Ensure there are exactly 10 defined seed values
  $#seed = 9;
  for (0..9) { $seed[$_] = 0 unless defined $seed[$_]; }
  $_state = [0x61707865, 0x3320646e, 0x79622d32, 0x6b206574, @seed[0..7], 0, 0, @seed[8..9]];
  $_stream = '';
  $_have = 0;
  $_sptr = 0;
}
# TODO replace this
sub srand {
  my $seed = shift;
  $seed = CORE::rand unless defined $seed;
  seed_csprng(pack("L2", ($seed >> 32) & 0xFFFFFFFF, $seed & 0xFFFFFFFF));
  $seed;
}
sub irand {
  if ($_have < 4) {
    my $keystream = _keystream(64, $_state);
    $_stream = substr($_stream,$_sptr,$_have) . $keystream;
    $_sptr = 0;
    $_have = length($_stream);
  }
  $_have -= 4;
  $_sptr += 4;
  return unpack("L",substr($_stream, $_sptr-4, 4));
}
sub irand64 {
  if ($_have < 8) {
    my $keystream = _keystream(64, $_state);
    $_stream = substr($_stream,$_sptr,$_have) . $keystream;
    $_sptr = 0;
    $_have = length($_stream);
  }
  $_have -= 8;
  $_sptr += 8;
  return unpack("Q", substr($_stream, $_sptr-8, 8));
}
sub random_bytes {
  my($bytes) = @_;
  $bytes = (defined $bytes) ? int abs $bytes : 0;

  if ($_have < $bytes) {
    my $keystream = _keystream($bytes, $_state);
    $_stream = substr($_stream,$_sptr,$_have) . $keystream;
    $_sptr = 0;
    $_have = length($_stream);
  }
  $_have -= $bytes;
  $_sptr += $bytes;
  return substr($_stream, $_sptr-$bytes, $bytes);
}

1;
