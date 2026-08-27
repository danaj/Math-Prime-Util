\\ Compare PARI/GP with bench/bitops-2048.pl using the same 2048-bit value,
\\ masks, bit positions, and cyclic 56-operation workload.
\\
\\   gp -q bench/bitops-2048.gp
\\   MPU_BITOPS_SECONDS=2 gp -q bench/bitops-2048.gp

WIDTH = 2048;
POSITION_COUNT = 8;
OPS_PER_ROUND = 56;

initial_bit(k) =
{
  if(k == WIDTH-1, 1, ((k * 73 + k \ 8 + 19) % 29) < 14)
};
and_bit(k) =
{
  if(k == WIDTH-1, 1, ((k * 37 + k \ 5 + 7) % 31) < 23)
};
or_bit(k) = ((k * 43 + k \ 9 + 11) % 37) < 15;
xor_bit(k) =
{
  if(k == WIDTH-1, 0, ((k * 53 + k \ 7 + 3) % 41) < 19)
};
andnot_bit(k) =
{
  if(k == WIDTH-1, 0, ((k * 61 + k \ 6 + 17) % 43) < 17)
};

build_values() =
{
  my(initial = 0, and_mask = 0, and_restore = 0,
     or_mask = 0, or_keep = 0, xor_mask = 0,
     andnot_mask = 0, andnot_restore = 0);
  my(ib, ab, ob, xb, nb);

  forstep(k = WIDTH-1, 0, -1,
    ib = initial_bit(k);
    ab = and_bit(k);
    ob = or_bit(k);
    xb = xor_bit(k);
    nb = andnot_bit(k);
    initial        = 2 * initial        + ib;
    and_mask       = 2 * and_mask       + ab;
    and_restore    = 2 * and_restore    + (ib && !ab);
    or_mask        = 2 * or_mask        + ob;
    or_keep        = 2 * or_keep        + (ib || !ob);
    xor_mask       = 2 * xor_mask       + xb;
    andnot_mask    = 2 * andnot_mask    + nb;
    andnot_restore = 2 * andnot_restore + (ib && nb);
  );
  [initial, and_mask, and_restore, or_mask, or_keep, xor_mask,
   andnot_mask, andnot_restore]
};

pick_positions(wanted, offset) =
{
  my(v = vector(POSITION_COUNT), n = 0, k = offset);
  while(n < POSITION_COUNT,
    if(initial_bit(k) == wanted, n++; v[n] = k);
    k = (k + 251) % (WIDTH-1);
  );
  v
};

values = build_values();
initial        = values[1];
and_mask       = values[2];
and_restore    = values[3];
or_mask        = values[4];
or_keep        = values[5];
xor_mask       = values[6];
andnot_mask    = values[7];
andnot_restore = values[8];
set_pos   = pick_positions(0, 17);
clear_pos = pick_positions(1, 89);
flip_pos  = vector(POSITION_COUNT, i, (137 + 337 * (i-1)) % (WIDTH-1));
x = initial;

bit_work() =
{
  for(i = 1, #set_pos,   bitset(~x, set_pos[i]));
  for(i = 1, #set_pos,   bitclear(~x, set_pos[i]));
  for(i = 1, #clear_pos, bitclear(~x, clear_pos[i]));
  for(i = 1, #clear_pos, bitset(~x, clear_pos[i]));
  for(i = 1, #flip_pos,  bitflip(~x, flip_pos[i]));
  for(i = 1, #flip_pos,  bitflip(~x, flip_pos[i]));
  x = bitand(x, and_mask);
  x = bitor(x, and_restore);
  x = bitor(x, or_mask);
  x = bitand(x, or_keep);
  x = bitxor(x, xor_mask);
  x = bitxor(x, xor_mask);
  x = bitnegimply(x, andnot_mask);
  x = bitor(x, andnot_restore);
};

measure_rate(minimum_ms) =
{
  my(batch = 1, calls = 0, t, elapsed, scale);

  while(1,
    t = getwalltime();
    for(i = 1, batch, bit_work());
    elapsed = getwalltime() - t;
    if(elapsed >= 20, break);
    scale = if(elapsed > 0, floor(20.0 / elapsed) + 1, 10);
    if(scale > 10, scale = 10);
    batch *= scale;
  );

  t = getwalltime();
  while(getwalltime() - t < minimum_ms,
    for(i = 1, batch, bit_work());
    calls += batch;
  );
  elapsed = getwalltime() - t;
  [1000.0 * calls / elapsed, elapsed]
};

duration = getenv("MPU_BITOPS_SECONDS");
if(duration == 0, duration = 0.75, duration = eval(duration));
if(duration <= 0, error("MPU_BITOPS_SECONDS must be positive"));

bit_work();
if(x != initial, error("PARI/GP workload did not restore its input"));
result = measure_rate(1000.0 * duration);
if(x != initial, error("PARI/GP workload did not restore its input"));

printf("PARI/GP 2048-bit cyclic bit-operation benchmark\n");
printf("%d operations/round, %.2fs measurement\n", OPS_PER_ROUND, duration);
printf("rounds/s: %.1f\n", result[1]);
printf("ns/op:   %.1f\n", 10^9 / (result[1] * OPS_PER_ROUND));
