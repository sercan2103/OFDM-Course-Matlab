# Presentation Build Instructions
## OFDM Term Project — June 12, 2026 at 14:30
## Duration: ~10 min talk + 5 min Q&A | Tool: PowerPoint (.pptx)

---

## AESTHETIC GUIDELINES (apply to every slide)

- **Color palette:** Dark navy background (#0D1B2A) with white text, accent color electric blue (#00B4D8) for highlights, soft orange (#FF9E00) for warnings/emphasis, light gray (#C9D6DF) for secondary text. Alternatively use a clean white background with dark navy text and the same accent colors — whichever looks more professional to you.
- **Font:** Use a modern sans-serif throughout. Recommended: "Calibri", "Inter", or "Segoe UI". Headings bold 28–32pt, body 18–20pt, captions 14pt.
- **No walls of text.** Each slide should have at most 4–5 bullet points. If you feel the urge to write a paragraph, cut it in half.
- **Every slide needs a visual.** Equations, diagrams, plots, tables, or block diagrams. No slide should be text-only.
- **Consistent layout.** Use the same title bar position, same margin widths, same footer on every slide.
- **Footer on every slide (small, bottom center):** "Kemal Sercan İlhan — OFDM Spring 2026"
- **Slide numbers** bottom right.
- **Transitions:** None or very subtle fade only. No spinning, flying, or bouncing animations.
- **Figures:** All axis labels readable (min 14pt), grid lines present, legend present. Export at 150 dpi minimum.

---

## SLIDE-BY-SLIDE INSTRUCTIONS

---

### SLIDE 1 — Title Slide

**Content:**
- Title (large, centered): **"High-Doppler OFDM Receiver Design"**
- Subtitle: *"Beating the Tier-1 Baseline at f_d = 3000 Hz"*
- Name: Kemal Sercan İlhan — 22050211002
- Course: OFDM Spring 2026 — Term Project
- Date: June 12, 2026

**Visual:**
- Background: a subtle waveform or constellation diagram graphic (blurred, low opacity) to fill dead space. Or a clean gradient background.
- The title text should feel bold and confident — this is the first impression.

**Do NOT put:** bullet points, tables, or any technical content.

---

### SLIDE 2 — The Problem (30 seconds)

**Title:** "The Mission"

**Content (two columns):**
Left column — "What we are given":
- One OFDM frame (10 symbols, QPSK, uncoded)
- FPV_C drone channel (4 taps, Jakes, f_d up to 3000 Hz)
- Residual STO: ±32 samples
- Residual CFO: ±12 kHz
- AWGN sweep 0–30 dB Eb/N0

Right column — "What we must build":
- Receiver: ofdm_rx(rx_samples, params) → bits
- Beat the Tier-1 baseline at **f_d = 3000 Hz**
- Low operation count (profiler-graded)

**Visual:**
- Simple TX→Channel→Impairments→RX block diagram. Draw it as a horizontal flow:
  [Bits] → [QPSK+Pilots] → [IFFT+CP] → [FPV_C Channel] → [STO+CFO+Doppler+AWGN] → **[YOUR RX]** → [Bits]
  Highlight the "YOUR RX" box in accent blue. The impairments box in orange.

---

### SLIDE 3 — Why the Baseline Fails (1 minute — this is the most important slide)

**Title:** "Why Tier-1 Fails at f_d = 3000 Hz"

**Content:**
- dense_4sym: pilots only on symbols {0, 3, 6, 9}
- Gap between pilot symbols = 3 symbols = **130 µs**
- Jakes correlation: **J₀(2π × 3000 × 130µs) ≈ 0**
- Channel fully decorrelated — no interpolation can bridge this gap
- Result: ~7.5% BER floor at 30 dB regardless of SNR

**Visual (CRITICAL — this slide must have a strong visual):**
Draw a timeline of 10 OFDM symbols (boxes labeled 0–9). Color the pilot symbols (0,3,6,9) green and the rest gray. Draw red "?" marks on symbols 1,2,4,5,7,8 to show unknown channel. Below the timeline, show the Jakes J₀ curve: x-axis is time delay (0 to 200 µs), y-axis is correlation (0 to 1). Draw a vertical dashed line at 130 µs where J₀ ≈ 0. Label it "dense_4sym gap — correlation = 0". This is the killer visualization that explains everything.

**Bottom callout box (orange):** "Zero correlation between pilot snapshots = interpolation is blind = irreducible error floor"

---

### SLIDE 4 — Our Solution Overview (30 seconds)

**Title:** "Our Approach: M2 + M3 + M4 Combined"

**Content (table format, 4 rows):**
| Problem | Method Used | Key Idea |
|---|---|---|
| Residual STO | M3 — Phase-slope estimator | Linear phase ramp → fit slope across pilots |
| Residual CFO | M2 — Van de Beek CP correlation | CP is copy of tail → phase shift = CFO |
| Doppler tracking | Pilot rearrangement + M3 Wiener smoother | Comb pilots every symbol → J₀ ≈ 0.84 |
| ICI | M4 — Decision-directed cancellation | Subtract estimated ICI before re-equalizing |

**Visual:** The RX chain block diagram from the project description (STO→CFO→CP remove→FFT→Channel Est→Doppler→Equalizer→QPSK demap). Color each block differently: STO blue, CFO blue, Channel Est green, Equalizer green. Label each block with the method name (M2/M3/M4).

---

### SLIDE 5 — Pilot Pattern: The Big Win (1.5 minutes)

**Title:** "Pilot Pattern Rearrangement: comb_10sym_64"

**Content:**
Left side bullets:
- dense_4sym: 4 pilot symbols, gap = 130 µs, J₀ ≈ 0
- comb_10sym_64: **all 10 symbols** carry pilots, gap = 43 µs, **J₀ ≈ 0.84**
- 64-pilot cap respected: even RBs 3 pilots/RB, odd RBs 2 pilots/RB
- Average pilot spacing: 4.9 SCs (down from 6 SCs)

**Visual (side by side — two pilot grid diagrams):**
Left diagram: dense_4sym — show a 10×12 grid (symbols × SCs within one RB). Color pilot cells green on symbols 0,3,6,9 at even SCs. Rest is gray. Label "4 pilot symbols — 130 µs gap".
Right diagram: comb_10sym_64 — same grid but now pilot cells appear on ALL 10 symbols (at positions 1,5,9 for even RBs). Color these green. Rest gray. Label "10 pilot symbols — 43 µs gap".

Below both diagrams: two small J₀ bars or arrows: "J₀ = 0" (red, dense_4sym) vs "J₀ = 0.84" (green, comb_10sym_64).

**Bottom line (bold):** "This single change eliminated the Doppler dead zone."

---

### SLIDE 6 — STO + CFO Correction (1 minute)

**Title:** "Residual STO and CFO Correction"

**Content — two halves:**

Left half: **STO (M3 — Phase-slope)**
- STO of δ samples → phase ramp e^{-j2πkδ/N} across subcarriers
- Measure phase difference between adjacent pilots: Δφ = -2π·Δk·δ/N
- Solve for δ, average over all 64 pilot pairs × 10 symbols
- Range: ±128 samples >> spec ±32 ✓
- Why not CP correlator? → flat-top ambiguity for small offsets

Right half: **CFO (M2 — Van de Beek)**
- CFO rotates every sample → CP and tail of same symbol differ by e^{j2πΔf·N/Fs}
- Correlate CP against tail, extract angle, average over 10 symbols
- Range: ±15 kHz >> spec ±12 kHz ✓
- Unbiased under multipath (all taps within 288-sample CP)

**Visual:**
Left: a small diagram showing the subcarrier phase ramp — x axis is subcarrier index k, y axis is phase angle(H_LS[k]). Show a diagonal line (the ramp) with pilot points as dots. Arrow showing "slope = -2πδ/N".
Right: a diagram of one OFDM symbol time-domain: show CP block and data block. Draw an arrow from CP back to the tail of the data block with label "copy". Label the phase rotation Δφ between them.

---

### SLIDE 7 — Channel Estimation: Wiener Smoother (1.5 minutes)

**Title:** "Doppler-Aware Channel Estimation"

**Content:**
- Baseline: LS at pilots → spline frequency interpolation → nearest-neighbour in time
- Our upgrade: **Wiener temporal smoother** before frequency interpolation
- Model: H_LS[p,m] = H_true[p,m] + noise
- Temporal correlation: R_HH(i,j) = J₀(2π·f_d·|i−j|·T_sym) — Jakes model
- Smoother: **W = R_HH · (R_HH + σ²I)⁻¹**
- Applied: H_smooth = H_LS_pilots · W (64 pilots × 10 symbols matrix)
- σ² = 5×10⁻⁴ tuned to 30 dB grading point

**Visual (key visualization — 2x1 grid of subplots):**
Left plot: "Without Wiener smoother" — show 10 noisy LS estimates for one pilot subcarrier across symbols 1–10. Scatter plot (dots), true channel as smooth curve underneath. High noise scatter around the true curve.
Right plot: "With Wiener smoother" — same true curve, but now the smoothed estimates (connected line) hug the true curve much more closely. Show at f_d=3000 Hz the smoother follows the Jakes variation; at f_d=0 it averages to a flat line.

**Caption:** "At f_d=0: W averages all 10 → noise ↓√10. At f_d=3000: W follows Jakes variation while suppressing noise."

---

### SLIDE 8 — ICI Equalization (1 minute)

**Title:** "ICI-Aware Equalization (M4)"

**Content:**
- At f_d=3000 Hz: channel changes within one symbol (T_u ≈ 33 µs)
- Time-varying H → subcarriers leak into each other (ICI)
- ICI model (linearly varying H):
  Y[k] = H_avg[k]·X[k] + Σ dH[k−Δk]·w(Δk)·X[k−Δk] + noise
- ICI weights: w(Δk) = 1 / (N · (e^{−j2πΔk/N} − 1))
- Algorithm: ZF → hard decisions → estimate dH → subtract ICI → re-equalize
- Q=6 neighbours each side → captures ~94% of ICI energy

**Visual:**
A bar chart showing ICI weight magnitude |w(Δk)| vs Δk for Δk = 1 to 10. The bars decay rapidly (1/Δk shape). Draw a vertical dashed line at Δk=6 with label "Q=6: 94% captured". Shade bars 1–6 in blue (captured), bars 7–10 in gray (ignored). This makes the Q=6 choice visually obvious.

Small side note: "Weights decay as 1/(2π|Δk|) — most ICI comes from nearest neighbours"

---

### SLIDE 9 — BER Results (1.5 minutes — the payoff slide)

**Title:** "Results: 7.5× Better than Tier-1 at f_d = 3000 Hz"

**Content:**
- Use the actual BER PNG you generated (İlhan_ber.png) as the main visual
- The plot already has all 4 Doppler curves + Tier-1 baseline dashed black line

**Annotations to add ON TOP of the plot (use PowerPoint text boxes/arrows):**
- Arrow pointing to Tier-1 purple curve at 30 dB with label: "Tier-1: 7.53×10⁻²"
- Arrow pointing to our fd=3000 curve at 30 dB: "Ours: 9.99×10⁻³ → **7.5× better**"
- Arrow pointing to fd=0 curve: "fd=0: 1.47×10⁻³ (static channel floor)"

**Below the plot, a small summary table:**
| f_d (Hz) | Tier-1 @ 30 dB | Ours @ 30 dB | Gain |
|---|---|---|---|
| 1000 | 1.54×10⁻² | 3.26×10⁻³ | 4.7× |
| 2000 | 4.50×10⁻² | 5.53×10⁻³ | 8.1× |
| **3000** | **7.53×10⁻²** | **9.99×10⁻³** | **7.5×** |

**Design note:** Make this slide feel like a "results" moment. Use bold colors, make the 7.5× number large and prominent.

---

### SLIDE 10 — Complexity (45 seconds)

**Title:** "Complexity: ~5.7×10⁷ MACs/frame"

**Content — two columns:**

Left: Profiler table
| Function | Time (s) | % |
|---|---|---|
| sto_estimate_correct | 0.042 | 46% |
| channel_estimate | 0.030 | 33% |
| equalize | 0.010 | 11% |
| cfo + cp + demap | 0.008 | 9% |
| **Total** | **0.091** | 100% |

Right: Op-count methodology
- flops() removed from MATLAB R13
- Calibrated MAC rate: 1.02×10⁹ MACs/s
- Algorithmic time: 0.056 s (total − rng overhead)
- **Estimated: ~5.7×10⁷ MACs/frame**

**Visual:** A horizontal bar chart showing % breakdown by block (same data as profiler table). Bars colored by category: STO (blue), channel_estimate (green), equalize (orange), rest (gray). Note under STO bar: "35 ms = rng overhead, not algorithm".

---

### SLIDE 11 — What Worked / What Didn't (45 seconds)

**Title:** "What Worked and What Didn't"

**Content — two columns with icons:**

✅ **Worked:**
- Pilot rearrangement → comb_10sym_64 (biggest gain)
- Van de Beek CFO (M2) — simple, robust, sufficient range
- Wiener time smoother — +4% complexity, −13–19% BER
- ICI cancellation Q=6 — captures 94% ICI energy

❌ **Tried and Reverted:**
- MMSE equalizer → σ²→0 (LS absorbs noise, wrong bottleneck)
- Wiener freq. interpolation → DC gap (711 bins) breaks matrix
- Pilot-based CFO refinement → adds Doppler noise > residual it fixes
- 2nd ICI iteration → −1.3% BER at +40% op-count (bad tradeoff)
- DD channel refinement → residual is time-domain (not freq-domain) problem

**Visual:** A 2-column layout with green checkmarks on left, red X marks on right. Keep it clean — no extra decoration needed, the content speaks for itself.

---

### SLIDE 12 — Summary & Key Takeaways (30 seconds)

**Title:** "Summary"

**Content (5 bullets maximum, large text):**
- The Tier-1 baseline fails because dense_4sym has J₀ ≈ 0 at f_d=3000 Hz
- Pilot rearrangement to comb_10sym_64 was the single most impactful fix
- Combined M2 (CFO) + M3 (STO + Wiener) + M4 (ICI) beats Tier-1 by **7.5×**
- Estimated complexity: ~5.7×10⁷ MACs/frame
- Residual floor: non-linear Jakes variation within the symbol — linear ICI model's limit

**Visual:** Repeat the BER result plot thumbnail on the right side (small), or a simple "scorecard" box:
```
fd=3000 Hz, 30 dB
Tier-1:  7.53×10⁻²
Ours:    9.99×10⁻³
Gain:    7.5×
```
Make the scorecard look like a green badge/callout.

---

### SLIDE 13 — Q&A Preparation Slide (have this ready but don't show unless needed)

**Title:** "Anticipated Q&A"

**Content (6 rows — question + 1-line answer):**
| Question | One-line answer |
|---|---|
| Why did Tier-1 fail? | J₀≈0 at 130 µs gap — zero correlation, interpolation is blind |
| Why comb_10sym over dense_4sym? | J₀=0.84 at 43 µs — strong correlation every symbol |
| How does Wiener smoother work? | W = R_HH(R_HH+σ²I)⁻¹, Jakes-shaped low-pass in time |
| How does ICI cancellation work? | Estimate channel slope dH, subtract Σ dH·w(Δk)·X_dec |
| If CFO jumped to 20 kHz, what breaks? | Van de Beek range ±15 kHz → ambiguous estimate → all blocks fail |
| Why not MMSE equalizer? | σ²→0 because LS absorbs noise — dominant error is interpolation, not noise |

**Design note:** This slide is a safety net for Q&A. Keep it simple — just the table. You can flip to it quickly if needed.

---

## SLIDE COUNT SUMMARY

| # | Title | Time |
|---|---|---|
| 1 | Title | — |
| 2 | The Mission | 30s |
| 3 | Why Tier-1 Fails ⭐ | 60s |
| 4 | Our Approach Overview | 30s |
| 5 | Pilot Pattern ⭐ | 90s |
| 6 | STO + CFO | 60s |
| 7 | Channel Estimation ⭐ | 90s |
| 8 | ICI Equalization | 60s |
| 9 | BER Results ⭐ | 90s |
| 10 | Complexity | 45s |
| 11 | What Worked / Didn't | 45s |
| 12 | Summary | 30s |
| 13 | Q&A backup | (hidden) |

⭐ = most important slides — spend the most time on these.
Total talk time: ~10 minutes. Perfect.

---

## IMPORTANT REMINDERS

1. **Slide 3 (Why Tier-1 Fails) and Slide 5 (Pilot Pattern) are the core of the presentation.** If the grader only remembers two things, it should be these. Make these slides visually striking.

2. **Do not read from the slides.** The bullets are cues, not scripts. The slides should show visuals while you talk over them.

3. **The BER plot (Slide 9) is your proof.** Point to specific curves and numbers. Practice saying "at 30 dB our curve sits here at 9.99×10⁻³, versus Tier-1 here at 7.53×10⁻², that is a 7.5× improvement."

4. **Be ready to be asked "walk me through one line of your code physically."** Prepare to explain `W_t = R_HH / (R_HH + sigma2 * eye(N))` from channel_estimate.m in plain words.

5. **The Q&A grader will push on WHY, not WHAT.** Every statement you make should have a because. "We used van de Beek *because* its range covers ±15 kHz which exceeds the ±12 kHz spec, and it requires no pilots so it doesn't consume any of our 64-pilot budget."

6. **Aesthetic reminder:** consistent fonts, no clipart, no WordArt, no default Office themes. If using PowerPoint defaults, at minimum choose "Office Theme" or "Facet" — avoid "Ion", "Wisp", or anything with heavy background patterns that compete with the content.
