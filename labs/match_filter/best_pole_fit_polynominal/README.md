# Lowpass 1st Polynomial

This experiment approximates the pointwise minimax denominator from `best_pole_fit`
with piecewise degree-7 Chebyshev polynomials. The input is the normalized cutoff
`x = 2 fc / fs`, so the coefficients work across sample rates within
`0.0008333333 <= x <= 0.9958333`.

The default fit uses 32 intervals uniform in `sqrt(x)` and nine minimax reference
designs per interval. The polynomial coefficients use a local least-squares fit
to those minimax denominator parameters; the resulting filter response is then
checked against the pointwise minimax response.
Each runtime evaluation selects one interval, evaluates its `p` and `q`
polynomials with Clenshaw's recurrence, then computes the biquad numerator from
the existing DC / cutoff / Nyquist power match. The Python constants
`PIECE_COUNT`, `POLYNOMIAL_DEGREE`, `FIT_CUTOFF_GRID_SIZE`,
`VALIDATION_CUTOFF_GRID_SIZE`, and `RESPONSE_GRID_SIZE` control regeneration.

Validation uses 129 points uniform in `sqrt(x)` and checks each response on
400001 points. The generator
refuses to emit a header when any validation design is unstable, infeasible, or
exceeds the configured `0.14 dB` error ceiling.

The outer cutoff grid is independent of the fixed per-response sample count:
`FIT_CUTOFF_GRID_SIZE` controls reference designs per segment, while
`RESPONSE_GRID_SIZE` controls each response minimax check. Coefficients use the
Chebyshev basis and are evaluated with Clenshaw recurrence in the generated C++.

Run `python qwqdsp/labs/match_filter/best_pole_fit_polynominal/lowpass_1st.py`
to recompute the coefficients, plot, and C++ header. The generated header is
`lowpass_1st_polynomial.hpp` and exposes
`qwqdsp_filter::best_fit_detail::LowpassOnepolePolynomial(wc)`.
