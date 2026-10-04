# 🐛🧊 Isolated reciprocal lambda-bound correction

🌳📌 Code commit2070fbd8f919d459e64231c540eec4f3069ce1e4 on feat/lambda-bound-wave-a, based on3be573960e6429e7d7506580812cb406e7fa8de3; source change is only the reciprocal lower ice lambda expression.
📨🔒 Diego explicitly authorized this narrow experimental commit through COOR_ASSISTANT, relayed by COOR messy pumas; missing dependency compile/CESM BFB suites are recorded rather than represented as passed.

## 📜🔍 Upstream origin

🐛📐 Upstream5434e2b009295889f073ffb0df8649605199b09d, authored/committed by Jian Sun on2022-05-12 with title “Merge the OpenACC optimization based on the PUMAS kernel”, split the correct reciprocal vector expression into scalar assignments and lost denominator grouping.
📜📋 Git blame attributes the older size-limit comment/bound to Katetc,2019-10-07; commit metadata attributes introduction of the faulty scalar line to Jian Sun, without proving who personally typed every change.

```diff
-  ice_lambda_bounds = 1._r8/[2._r8*dcs, 1.e-6_r8]
+  ice_lambda_bounds(1) = 1._r8/2._r8*dcs
```

🛠️📐 Experimental correction restores the intended grouping:

```fortran
ice_lambda_bounds(1) = 1._r8/(2._r8*dcs)
```

## 🧪📦 Executed evidence and limits

✅🔖 Committed utilities bytes match mode1 source SHA2562a8c811597c9305ef5dfc5cc54b42c1bdf5459467c5aadf15f3c64edd656f185 exactly; Levante job27879900 completed0:0 with GCC11.2 and assertions for actual reciprocal initialization at dcs500/dcs350 plus scalar/vector/2D activation, boundary, interior, upper-limit and empty-cloud cases.
✅📋 Mode0 baseline job27879817 separately validates legacy initialization/nonactivation and a reciprocal property-copy oracle; its output matches mode1 fixture output, which is utility-oracle agreement rather than full-core/host BFB.
🔖📦 Candidate wrapper SHA256ee720f81564d41408773ecfe1ea885ec89d914fc18edf166d0408c1c1eb4b1c7, fixture SHA25609863083a5e379dc788695b48812fa6fd5ed82c6e4cb9876eb0d1c32c02b7a7a, binary SHA2561be61b7afe5e3f2e547dd8a4bf3b4f5c987058ed8f1339fc2d578c13ee394524; all22 candidate artifact hashes verify,740KB preserved.
🛑📋 Initial mode0 job27879770 failed before compiler execution because gcc/13.3.0 was unavailable; its16KB failed evidence remains preserved.
🚧📘 Dependency tests/compile_check.sh and tests/bfb.sh run were not executed: required built CESM case modules/compiler log, Euler modules and restart baselines are unavailable on Levante; no full-core/CAM/GPU/model validation, golden renewal or scientific acceptance is claimed.
📘🔒 The AGENTS.md:11 precommit exception is scoped only to this tested one-line experimental commit; it is not a general test waiver or broader default-policy approval.

```text
/scratch/b/b380602/pumas_lambda_wave_a/lambda-candidate.27879900.DJnhrr
/scratch/b/b380602/pumas_lambda_wave_a/lambda-kernel.27879817.LWURPg
/scratch/b/b380602/pumas_lambda_wave_a/lambda-kernel.27879770.L7KRfO
/work/bb1143/b380602/DEV/wt/stability-pumas-lambda/handoffs/handoff.pumas_lambda_wave_a.md
```

## 🔒📋 Current state and next gates

🔒🌳 This commit lives only in the isolated experimental checkout/branch; upstream mirrors, production parent gitlink and shared MESSy/ICON/EMAC sources/builds remain unchanged.
🔒📋 Production pin update, promotion, host runs, golden baselines, public upstream reporting and issue16 closure remain held; the upstream-default policy still requires explicit disposition before integration.
🧮✅ Cumulative lambda jobs0.000225694444NH, zero reserved; no new jobs were run for committing the already-tested bytes.
🪞📘 Reusable lesson: record baseline Git identity separately from a tested precommit source hash, and distinguish utilities kernel proof from full dependency/CESM/host suites.
