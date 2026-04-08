---
name: virtual-chromosome-path-reader
description: Read and validate SKYPE virtual chromosome path files under `00_raw/`. Use when Codex needs to interpret a `*.paf` and matching `*.index.txt`, determine telomere-side orientation from `CTG_TELCON`, `CTG_DIR`, and current terminal state, group consecutive chunks by `CTG_NAME` into breakend units, propagate path orientation through `DIR_FOR`/`DIR_BAK` across each `chr_change` or `dir_change` transition, and decide whether a telomere-to-telomere path is orientation-consistent or contradictory.
---

# Virtual Chromosome Path Reader

Interpret SKYPE `00_raw/<chrpair>/<n>.paf` and `<n>.index.txt` as one virtual chromosome path. Always read both files together.

## Workflow

### 1. Read the path structure from `*.index.txt`

Interpret the first column as graph state:

- `1` = `DIR_FOR`
- `0` = `DIR_BAK`
- `2` = `DIR_OUT_FOR`
- `3` = `DIR_OUT_BAK`
- `4` = `DIR_IN_FOR`
- `5` = `DIR_IN_BAK`

Current code also writes telomere tuple lines as:

- `(chr_label, layer, chr_change, dir_change)`

For example:

- `('chrXb', 0, 1, 1)` means `chrXb` on `DIR_BAK` layer, with `chr_change=1`, `dir_change=1`

Legacy outputs may still use:

- `2` = generic `DIR_OUT`
- `3` = generic `DIR_IN`
- tuple lines like `(chr_label, chr_change, dir_change)`

If a historical example and current code disagree, prefer current code semantics.

### 2. Anchor the path with the endpoint chunks

Do not start by reading middle chunk order. Start from the telomeres.

For the first and last chunk, read:

- `CTG_TELCON`
- `CTG_DIR`
- whether the state is `DIR_OUT_FOR`/`DIR_OUT_BAK` or `DIR_IN_FOR`/`DIR_IN_BAK`
- the layer stored in the start/end tuple lines

If an endpoint chunk has `CTG_TELCON=0`, use the corresponding tuple-line telomere label from `*.index.txt` as the implied telomere label for that endpoint chunk.

Telomere-side rule:

- If `CTG_DIR='+'`:
  - `*f` means the telomere is on the contig front
  - `*b` means the telomere is on the contig back
- If `CTG_DIR='-'`:
  - `*f` means the telomere is on the contig back
  - `*b` means the telomere is on the contig front

Endpoint state rule:

- `DIR_OUT_FOR` and `DIR_OUT_BAK` mean the path starts at that telomeric side and leaves the chunk toward the interior.
- `DIR_IN_FOR` and `DIR_IN_BAK` mean the path comes from the interior and ends at that telomeric side.

Current implementation note:

- `telomere_connected_list.txt` stores the telomere-connected chunk layer for the outbound direction.
- The inbound terminal state is the flipped layer, not the same layer.
- So do not assume `DIR_IN_BAK` automatically means "`*b` telomere is correct" or `DIR_IN_FOR` automatically means "`*f` telomere is correct".
- Always combine terminal state, propagated orientation, `CTG_DIR`, and the telomere label.

Use these two rules together to decide the local orientation of the endpoint chunk in the current virtual chromosome.

### 3. Propagate orientation through the middle of the path

After fixing the start orientation, propagate it through the internal chunks.

- `DIR_FOR` means keep the current path orientation.
- `DIR_BAK` means interpret that segment as reverse-complement relative to the stored chunk order.
- `CTG_NAME` is the first column of `*.paf`. Consecutive rows with the same `CTG_NAME` belong to one contig/unitig block and should be treated as one breakend unit.

For chunks from the same contig:

- increasing chunk index under `DIR_FOR` is consistent
- decreasing chunk index under `DIR_BAK` is consistent

Breakend-by-breakend rule:

- Do not validate only the two telomere ends.
- In `*.index.txt`, every increase of `chr_change` or `dir_change` marks a transition that must be checked locally.
- Around that transition, inspect the surrounding consecutive `CTG_NAME` block on each side and decide whether the current orientation changes consistently through that breakend unit.
- If the new node index is larger than the previous node index, the new local orientation follows the new chunk's `CTG_DIR`.
- If the new node index is smaller than the previous node index, the new local orientation is the flipped value of the new chunk's `CTG_DIR`.
- This means statements like `15- => X-` or `15- => X+` are not summary labels only. Each one must be justified by the specific breakend unit crossed at that step.

Do not infer endpoint orientation from `DIR_FOR`/`DIR_BAK` alone. Endpoint chunks are special and must stay anchored by `CTG_TELCON` plus terminal state and telomere tuple layer.

### 4. Check breakend and terminal consistency

A valid virtual chromosome must have telomeres at both ends in the current propagated orientation.

Run this checklist:

1. Confirm the starting telomere matches the first chunk's `CTG_TELCON`.
2. Determine the starting chunk orientation from `CTG_TELCON`, `CTG_DIR`, and `DIR_OUT_FOR`/`DIR_OUT_BAK`.
3. Group consecutive same-`CTG_NAME` rows into breakend units.
4. For each increase of `chr_change` or `dir_change`, verify that the orientation transition implied by the neighboring breakend units is locally consistent.
5. Within each unit, confirm that `DIR_FOR` follows stored order and `DIR_BAK` follows reverse order.
6. If the final chunk has `CTG_TELCON=0`, replace it with the final tuple telomere label from `*.index.txt`.
7. Confirm the final chunk's telomere label lands on the terminal side implied by `DIR_IN_FOR`/`DIR_IN_BAK`.
8. Apply the current-code terminal flip rule: inbound telomere closure uses the flipped layer relative to the outbound telomere-connected layer.
9. If any local breakend transition or the final telomere side is inconsistent, mark the path as contradictory.

## Common Pitfall

One earlier misread came from treating terminal `BAK/FOR` as if it directly matched the telomere label suffix.

- Wrong shortcut: `DIR_IN_BAK` + `chrXb` => automatically valid
- Correct rule: terminal `BAK/FOR` must still be checked against propagated orientation, `CTG_DIR`, and the inbound flip rule

The historical `U2OS_20_39_04/00_raw/chr15b_chrXb/27` case looked valid under the shortcut above, but it was actually contradictory. After the graph-side terminal fix in `02_Build_Breakend_Graph_Limited.py`, that `chr15b_chrXb` raw folder no longer appears in regenerated output.

## Output Format

When analyzing a path, report:

1. Start telomere and how the first chunk is oriented.
2. Breakend units in order, using `CTG_NAME`.
3. Middle path orientation transitions, one breakend crossing at a time.
4. End telomere and whether it matches the propagated orientation.
5. Final verdict: `consistent` or `contradictory`.

## Reference

For a worked example, a compact checklist, and notes on legacy vs current output formats, read [references/path-checking.md](references/path-checking.md).
