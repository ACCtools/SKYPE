# Depth NNLS convergence

By default, stage 23 fits every stage-22 feature with the unregularized objective
`min ||A w - B||²`, subject to `w >= 0`, without an intercept. It uses
`scipy.optimize.nnls` in Float64, with `adelie.solver.bvls` as the working-set
fallback. The historical `nnls_solver.jl` is not called by this stage.

For more than 2,048 feature columns, the fit grows a working set in batches of
32. Each subproblem uses the same SciPy NNLS routine. A full gradient over every
column determines which excluded columns could improve the objective. Previously
admitted columns stay in subsequent subproblems, including those whose current
coefficient is zero. This avoids cycling among nearly interchangeable profiles.
The working solve checks positive and zero coefficients against a scaled KKT
tolerance of `1e-12`; it does not stop just because no new column was selected.
Small matrices continue to use one full solve.

When a subproblem fails, returns invalid coefficients, increases the objective
beyond rounding tolerance, leaves an unresolved admitted-column violation, or
reaches the 128-round limit, the fit switches to Adelie on all aggregated
columns. Non-finite working-set KKT values also trigger this fallback. Adelie
uses nonnegative bounds, Float64, `tol=1e-7`, `max_iters=100000`, and one internal
thread. It is imported only if a fallback is needed. Its dense matrix wrapper
reuses the aggregated matrix and avoids a temporary squared copy of every entry.
The fallback starts a new Adelie solve and does not retry full SciPy NNLS.
No columns are permanently filtered, no penalty is added, and no earlier graph
solution is needed. Small problems use full SciPy NNLS and accept its returned
result without a separate KKT check.

Consecutive rows with exactly identical feature values are combined before
solving. For a run of length `k`, the equivalent row and target are
`sqrt(k) * A[i, :]` and `sum(B[run]) / sqrt(k)`. The original squared residual
equals the compressed squared residual plus the within-run variance of B,
which is independent of w. No feature columns, event identities, observations,
or depth information affecting the minimizer are removed. This is an exact
least-squares transformation; it does not denoise B. Only the numerical solve
uses these combined rows. Predictions and denoising still use the original bins.

The working-set loop is the only runtime KKT check. Its successful final
iteration provides the recorded KKT value, and there is no post-solve KKT
calculation or acceptance threshold. With `g = A.T @ (A @ w - B)`, the KKT
violation is `abs(g[j])` for positive coefficients and `max(-g[j], 0)` for zero
coefficients. Each violation is divided by `||A[:, j]|| * ||B_original||`.
The gradient and column norms are mathematically preserved by row aggregation;
using the original target norm preserves the KKT scale even when within-run
target noise is removed from the numerical solve. Zero columns contribute zero;
when the original B is zero, its normalizer is one. Both full SciPy results and
Adelie fallback results use the solver's own stopping criteria. Adelie's
tolerance measures coordinate-descent convergence
and is not a bound on this normalized KKT value. Non-finite or negative
coefficients, invalid predictions, and errors reported by Adelie still fail the
stage. The original-row gradient is never rechecked in the fit. The certificate
helper remains available for explicit numerical audits.

`nnls_diagnostics.json` records the backend, arithmetic precision, number of
original and combined rows, feature count, squared residual, relative error,
and the aggregated-matrix KKT value (`kkt_matrix: "aggregated"`,
`kkt_tolerance: 1e-12`) for successful working-set results.
An Adelie fallback records `strategy: "adelie_fallback"`,
`solver: "adelie.solver.bvls"`, the working-set fallback reason, and Adelie's
tolerance and iteration limit. Full SciPy and Adelie fallback results both use
`convergence_check: "solver_default"`; their `kkt_scaled_max`, `kkt_tolerance`,
`kkt_matrix`, and `converged` are `null` because strict KKT convergence was not
checked. Working-set results use `convergence_check: "aggregated_kkt"` and
`converged: true`. Squared residual
and relative error are still computed on the original bins; they require only
prediction with nonzero coefficients, not a gradient over all feature columns.
`weight.npy` retains one coefficient per original feature,
including zero columns, and `A_idx_list.pkl` retains the original feature order.
The diagnostics also record the strategy, working-set rounds, largest subproblem,
and any fallback reason. Prediction skips exactly zero coefficients and reads
nonzero columns in bounded blocks, avoiding a Float64 copy of the entire matrix.

The pipeline always passes one BLAS thread to stage 23 for fitting, optimality
checks, and prediction, then the stage restores the previous settings. The
`run.sh --threads N` option controls the other stages independently. Direct
stage-23 invocation accepts `-t`/`--thread`/`--threads` for benchmarking and
defaults to one thread. Public functions in `depth_nnls.py` inherit their
caller's BLAS settings. This uses the existing `threadpoolctl` dependency
without changing environment variables.
Diagnostics record both the requested limit and each loaded BLAS library's
effective limit, which can be lower if that library has a build-time cap.
More threads do not guarantee a faster fit: repeated matrix-vector products
and small subproblems can cost more to coordinate than they gain in parallelism.

This reduces repeated least-squares work on large graphs. Reading the input and
checking all columns still costs time proportional to matrix size. Working-set
overhead can be visible on medium-sized problems. The benchmark must compare
both implementations at the same precision and acceptance tolerance.

The previous Adelie defaults could stop with coordinate updates small relative
to total target energy while still having a material optimality residual.
Its screening check tests unscreened features, not convergence of every
already-screened coordinate. A successful return was therefore insufficient
evidence of the accuracy needed when comparing graph limits. The fallback now
deliberately accepts that lower convergence accuracy to avoid a slow full SciPy
solve; graph-limit raw-error monotonicity is not assured on fallback runs.

For fixed B and nested columns, the exact raw NNLS optimum cannot increase when
the graph adds columns: the previous solution, extended by zeros, remains
feasible. It can stay equal. The denoised-target score is a different objective
and has no such monotonicity guarantee. Neither score establishes SV accuracy.
Rank-deficient matrices can also have different nonnegative coefficients with
the same optimal prediction; changing solver may change individual SV weights.

A graph warm start is optional, not necessary for this solver. A future warm
start must match exact column profiles and identical target coordinates, sum
weights when matching duplicate columns, and retain a feasible incumbent if a
new numerical result has higher raw SSE. Copying coefficients by position is
invalid when graph enumeration changes feature order. If graph truncation
removes previous columns, the nesting guarantee no longer follows.

## Optional high-depth policies

Native runs accept `--depth-policy` through the existing `--option_skype` string:

```bash
bash run.sh --skype-start-at 22 \
  --option_skype '--depth-policy nclose_huber --robust-sigma-multiplier 3' \
  H2009 path/to/result
```

The default `legacy` policy retains the existing depth cutoff. `all` restores
amplitude-only exclusions, while `nclose_l2` and `nclose_huber` restore only
contiguous high-depth blocks intersecting an NClose's original aligned segments.
The NClose gate uses the current graph topology, including zero-weight paths;
it does not require raw-read support or a positive fitted NClose weight.
Existing CenSat exclusions remain in place. NClose-gated policies apply to the
native assembly route. Stage 22 builds the current structure model before the
gate, so an earlier fitting result is not required.

For `nclose_huber`, original clean bins retain their quadratic loss. Restored
bins use one-sided Huber loss: with residual `r = A w - B`, the loss is
`r² / 2` for `r >= -tau`, and `-tau*r - tau²/2` otherwise. `tau` is the positive
sigma multiplier times the chromosome noise estimated from original clean bins.
The solver uses L-BFGS-B followed by weighted NNLS and checks the actual robust
gradient against a scaled KKT tolerance of `1e-10`. It retains individual
high-depth observations while aggregating only quadratic clean rows. This is a
different objective and convergence certificate from the default NNLS above.

If there are no high-depth candidates, `high_depth_gate.tsv` is a header-only
table with a stable schema. If no candidate passes the NClose gate, the report
retains those candidates with `included=false`. In both cases `nclose_huber`
calls the existing raw NNLS directly at stage 23, records
`high_depth_fallback: "no_high_depth_rows"`, and writes an empty
`unexplained_high_depth.npy`. The robust optimizer and any cached initial weights
are unnecessary for this fallback.

With restored rows, `unexplained_high_depth.npy` records
`max(B - A w - tau, 0)` in the order of `high_depth_rows` in `23_input.pkl`.
These values are diagnostics, not extra path depth or SV support. Stage 31 uses
the same saved row coordinates for plots and VCF depth annotations; it never
reconstructs a separate amplitude mask. Changing the policy requires rerunning
from stage 22 so matrix rows and their metadata stay aligned.
