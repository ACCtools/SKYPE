# Path Checking

## Compact Checklist

1. Open both `n.paf` and `n.index.txt`.
2. Read start/end telomere labels from the tuple lines.
3. Read the first-column state using current code semantics:
   `0=DIR_BAK`, `1=DIR_FOR`, `2=DIR_OUT_FOR`, `3=DIR_OUT_BAK`, `4=DIR_IN_FOR`, `5=DIR_IN_BAK`.
4. For the first chunk, use `CTG_TELCON` and `CTG_DIR` to locate the telomere on contig front/back.
5. Use `DIR_OUT_FOR` or `DIR_OUT_BAK` to decide that the path starts at that telomeric side and moves inward.
6. Treat the first `*.paf` column, `CTG_NAME`, as the breakend-unit key. Consecutive rows with the same `CTG_NAME` are one local unit.
7. Every increase of `chr_change` or `dir_change` in `*.index.txt` is a local transition that must be checked with the neighboring `CTG_NAME` units.
8. Within a unit, `DIR_FOR` should follow stored chunk order and `DIR_BAK` should follow reverse chunk order.
9. If the last chunk has `CTG_TELCON=0`, use the terminal tuple telomere label from `*.index.txt` as the implied endpoint telomere.
10. For the last chunk, use `DIR_IN_FOR` or `DIR_IN_BAK` plus that telomere label to see whether the propagated orientation ends on the terminal telomeric side.
11. Remember that current code closes the terminal telomere on the flipped inbound layer relative to the outbound telomere-connected layer.
12. If any local transition or the terminal side is inconsistent, the path is contradictory.

## Current vs Legacy Output

Current code writes:

- terminal states `2/3/4/5`
- tuple lines `(chr_label, layer, chr_change, dir_change)`

Older output may still show:

- generic `2=DIR_OUT`, `3=DIR_IN`
- tuple lines `(chr_label, chr_change, dir_change)`

When reading old examples, keep the old file content as reference material only. When checking current output, prefer the current code semantics above.

## Breakend Unit Rule

- `CTG_NAME` is the original contig/unitig name and corresponds to the first `*.paf` column.
- Consecutive rows with the same `CTG_NAME` should be read together as one breakend unit.
- Do not skip directly from the start telomere to the final telomere.
- The important checks happen at each place where `*.index.txt` increments `chr_change` or `dir_change`.
- At each such step, decide whether the path moves into the next unit in stored order or reverse-complement order.
- Only after every local transition is consistent should the path be called globally consistent.

## Worked Example

This is a legacy-format contradictory example. It is still useful because the endpoint inconsistency is biologically real, but the terminal state encoding predates the current doubled terminal states.

Original file paths:

- `/home/hyunwoo/60g_skype/30_skype_pipe/U2OS_20_58_22/00_raw/chr2b_chr15b/1.paf`
- `/home/hyunwoo/60g_skype/30_skype_pipe/U2OS_20_58_22/00_raw/chr2b_chr15b/1.index.txt`

Inline copy of `1.paf`:

```text
('chr2b', 0, 0)
ptg000473l	27864	0	9209	+	chr2	242696752	135664279	135673488	60	3	2247	2247	0	0	chr2b	0	0	0	+	chr2	0.5936
utg058613l	31782	0	11797	-	chr2	242696752	31270	43066	60	1	5820	5821	0	0	0	0	0	0	-	chr15	1.70292
utg058613l	31782	11797	31782	-	chr15	99753195	17445128	17465114	0	1	5820	5821	0	0	0	chr15	rin	rin	-	chr15	1.70293
utg075486l	13248	931	13248	-	chr15	99753195	13034119	13046488	1	3	6228	6228	0	0	chr15b	chr15	rin	rin	-	chr15	1.89366
('chr15b', 1, 0)
[('chr2', 135642218), ('chr15', 4430995)]
```

Inline copy of `1.index.txt`:

```text
('chr2b', 0, 0)
2	2247	0	0
1	5820	0	0
1	5821	1	0
3	6228	1	0
('chr15b', 1, 0)
[('chr2', 135642218), ('chr15', 4430995)]
```

Index states:

- `('chr2b', 0, 0)`
- `(legacy DIR_OUT=2, 2247, 0, 0)`
- `(DIR_FOR, 5820, 0, 0)`
- `(DIR_FOR, 5821, 1, 0)`
- `(legacy DIR_IN=3, 6228, 1, 0)`
- `('chr15b', 1, 0)`

Interpretation:

- Chunk `2247` has `CTG_TELCON=chr2b`, `CTG_DIR='+'`.
- Therefore the telomere is on the contig back side.
- Because the path starts at `chr2b` and the first internal state is legacy `DIR_OUT`, the virtual chromosome must start from that back telomeric side.
- So the current path orientation begins as `2-`.
- The middle states are `DIR_FOR`, `DIR_FOR`, so the path orientation is propagated without flipping: `2- => 15-`.
- Chunk `6228` has `CTG_TELCON=chr15b`, `CTG_DIR='-'`.
- For a `-` chunk, `b` means the telomere is on the contig front side.
- But the propagated orientation is already `15-`, so `chr15b` lands on the wrong side for a terminal telomere.

Verdict:

- This path is `contradictory`, not a valid telomere-to-telomere virtual chromosome under this orientation rule.

## Correction Note

One easy mistake is to treat terminal `BAK/FOR` as if it directly matches the telomere suffix.

- Incorrect shortcut: `DIR_IN_BAK` + `chrXb` => valid
- Correct rule: check propagated orientation, `CTG_DIR`, telomere label, and the inbound flip rule together

That shortcut caused an earlier wrong answer for the historical `U2OS_20_39_04/00_raw/chr15b_chrXb/27` example. After the graph-side terminal fix, regenerated U2OS output no longer produces a `chr15b_chrXb` raw folder for that case.
