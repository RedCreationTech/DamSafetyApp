# SP1 test-only fixtures

square.mtx/e833.mtx are real expert Abaqus2025 initial free8x8 CPE4R matrices, mm/N, section1mm. Tests use SI coordinates and unit width1m, so target K=1e6*native. No native nonlinear/FINITE or dynamic result is claimed.

P01..P06 and four material tables are frozen six real-input fixtures from7e0328e. Reused only for material-core replay; they are NOT the future one-QP center trajectory. CDP production implementation remains897.

All source hashes in source-identity.json. No test output/ODB/SIM is included.

`959_v2_3_accepted_path.txt` is the frozen v2-3 element 959/IP1 native accepted path: Gravity 1–24 plus Earthquake 1–163, 187 increments from zero to EQ2.80. Each row contains only stage/increment, adjacent native times, dt and the six D components. It contains no Abaqus plastic or damage state. SHA-256: `cf9e6bfe276c60d78b6927f0a3459dcebfebc735bfd7a38d2280fdfe42328315`. Preparation qualification verifies positive det(F), cond(A+I)≤1e12, finite values, positive dt and exact binary64 readback; EQ0 is not replayed.

`959_v2_3_mid_proxy.csv` is the existing M_mid rotation-bearing elastic proxy, Gravity 0–24 plus Earthquake 1–94 (119 rows, 118 increments, through EQ1.68). SHA-256: `535e066d5ecdf68740dde02626273143573cacdba53db75ce1ff63339cfb94f4`. `DS04GPeakReloadReplay.Native959AcceptedPath` uses only its four `midpoint_cayley_S_Pa_*` columns to report native-minus-proxy stress differences; this is not a correctness gate. Later frames emit null proxy fields rather than fabricated values.

The test runs native, dt_plus and dt_minus arms independently from zero old16/total state, committing only successful native updates. It emits 561 frame JSONL records plus arm/clock summaries. Hard gates cover input identity/order/time continuity, finite values, total-strain roundoff, immutable old_state, native exceptions, at most 256 substeps, state/damage ranges and irreversibility. The unrotated elastic recurrence uses only the technical budget `1e-5 Pa + 1e-12*stress_scale` before native plasticity. No physical-closeness threshold or strict-clock closure is asserted.
