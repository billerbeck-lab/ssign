# Secretion-system substrate localization (per type)

Reference for how each secretion system that ssign handles delivers its
substrate, and therefore what each of ssign's three "secreted protein" signals
(DeepLocPro localization, SignalP signal peptide, DeepSecE substrate-type) should
and *does* report. Compiled from a literature sweep (citations below, with DOIs),
cross-checked against ssign's own predictions on a curated effector set.

## The two-family framework (biology)

Secretion systems split cleanly by *where the substrate ends up*:

**Family A — envelope-crossing, substrate ends up outside the producing cell.**
T1SS, T2SS, T5SS (a/b/c). The mature protein is extracellular or surface-attached.
Localization predictors are the *right* instrument here.

**Family B — injection, substrate is delivered into another cell.**
T3SS, T4SS, T6SS. The effector is a cytoplasmic protein in the producing cell
until it is fired/injected directly into a host or neighbouring cell; it never
free-floats outside. By strict producer-cell biology these are cytoplasmic, so
localization predictors "should" miss them. (But see the empirical section — in
practice DeepLocPro calls many of them extracellular anyway, because it is trained
on functional "secreted" annotations.)

The SignalP axis is orthogonal: a substrate is SignalP-positive only if it crosses
the inner membrane via Sec/Tat with a cleavable N-terminal signal peptide.

## Per-type summary

| Type | Mechanism | Sec signal peptide (SignalP) | Substrate destination | Family |
|---|---|---|---|---|
| **T1SS** | 1-step, ABC/MFP/TolC, no periplasm | **No** (uncleaved C-terminal signal) | Extracellular / surface (RTX toxins, serralysins, S-layer) | A |
| **T2SS** | 2-step: Sec/Tat → periplasm → OM | **Yes** (Sec SPI/SPII, or Tat) | Extracellular / surface (cholera toxin, cellulases, lipases) | A |
| **T3SS** | 1-step injectisome, cytoplasm → host cytosol | **No** (N-terminal + chaperone-binding signal, uncleaved) | Host-cell cytoplasm | B |
| **T4SS** | 1-step, VirD4 coupling protein, cytoplasm → recipient | **No** (C-terminal / internal signal) | Host/recipient cytoplasm | B |
| **T5aSS** | Sec → periplasm; β-barrel in OM, passenger to surface | **Yes** (often ESPR-extended) | Surface-retained OR released (β-barrel stays OM) | A |
| **T5bSS** | Two-partner: TpsA + TpsB, both Sec-exported | **Yes** (both partners) | TpsA released extracellular; TpsB stays OM | A |
| **T5cSS** | Trimeric autotransporter, Sec → periplasm → OM trimer | **Yes** | Surface-retained (uncleaved adhesin); β-barrel OM | A |
| **T6SS** | 1-step contractile injection, cytoplasm → target | **No** | Target-cell cytoplasm or periplasm | B |

Notable exception: **pertussis-toxin / Ptl-type T4SS** is genuinely two-step (Sec →
periplasm → Ptl across OM → extracellular), so its subunits are SignalP-positive and
extracellular, the opposite of VirD4-dependent host-injected effectors.

## What DeepLocPro / DeepSecE actually predict (empirical)

From ssign's own predictions on curated corpus effectors (ssign's default runs
apply DLP/DSE to the ±3-gene **neighbourhood** of detected systems only, so n is
the subset of known effectors that landed in a neighbourhood: small for some
types, and biased toward "near machinery". Indicative, not definitive;
whole-genome predictions would tighten these):

| Type | DLP ran on | % called **Extracellular** (≥0.8) | DeepSecE call (where it ran) |
|---|---|---|---|
| T1SS | 15 | **100%** (15/15), mean prob 0.97 | 15/15 → T1SS |
| T2SS | 4 | **100%** (4/4), mean 0.92 | 3/4 → T2SS |
| T3SS | 21 | **71%** (15/21), mean 0.66 | 13/21 → T3SS (we flag/exclude these) |
| T4SS | 1 | 0% (0/1) → Cytoplasmic Membrane | 1/1 → T4SS |
| T6SS | 40 | **68%** (27/40), mean 0.55 | 28/40 → T6SS |

**The load-bearing empirical finding:** DeepLocPro calls a *majority* of T1/T2/T3/T6
effectors "Extracellular", **including the injected Family-B types T3SS and T6SS**.
First-principles biology says injected effectors are cytoplasmic and a localization
tool "should" miss them — but DeepLocPro is trained on database localization
annotations, where virulence effectors are typically labelled secreted, so it flags
them extracellular regardless of the injection nuance. This means:

- DLP-extracellular is a **usable** signal for T3SS and T6SS effectors, not a
  mismatch. (T4SS is the one type that empirically looks cytoplasmic — but n=1 here,
  so it is unverified; T4SS deserves its own whole-genome check.)
- DeepSecE's sequence classifier independently flags T1/T6 well, and T3SS too — but
  T3SS-DSE is the flagellar-false-positive source (DeepSecE has no flagellum class, so
  flagellar proteins — T3SS homologs — funnel into its T3SS bin and over-call it),
  so ssign deliberately excludes DeepSecE for T3SS and relies on DLP + proximity.

## Implications for ssign

1. **Localization-based substrate calling works across both families in practice**,
   because DeepLocPro empirically calls most effectors extracellular even for the
   injected types. The strict "injected ⇒ cytoplasmic ⇒ unfindable" worry does not
   hold for DLP as actually trained.
2. **SignalP is informative only for Family A** (T2SS, T5SS, Ptl-T4SS). It is
   correctly silent for T1SS (C-terminal signal) and the injected types — a
   SignalP-negative is expected, not evidence against secretion, for those.
3. **DeepSecE T3SS stays excluded** (flagellar FPs); T3SS leans on DLP + proximity.
   DeepSecE is genuinely the right sequence signal for T6SS (and would be for T3SS if
   the flagellar problem were solved).
4. **To make these numbers definitive**, re-extract whole-genome DLP/DSE predictions
   for known effectors per type (the default runs are neighbourhood-only). A
   whole-genome run (e.g. Salmonella LT2) is the natural source for a clean T3SS check.

## Citations (per type; DOIs with PubMed verification)

**T1SS:** Holland et al. 2016 EcoSal Plus 10.1128/ecosalplus.ESP-0019-2015; Spitz et al.
2022 Front Physiol 10.3389/fphys.2021.804646; Zhang et al. 1995 Biochemistry
10.1021/bi00013a007 (C-terminal uncleaved signal).

**T2SS:** Korotkov & Sandkvist 2019 EcoSal Plus 10.1128/ecosalplus.ESP-0034-2018;
Nivaskumar & Francetic 2014 BBA 10.1016/j.bbamcr.2013.12.020; Ball et al. 2016 Sci Rep
10.1038/srep27675 (Tat subset).

**T3SS:** Lara-Tejero & Galán 2019 EcoSal Plus 10.1128/ecosalplus.ESP-0039-2018;
Stebbins & Galán 2001 Nature 10.1038/35102073 (chaperone holds effector unfolded);
Myeni et al. 2013 PLoS One 10.1371/journal.pone.0060499 (SipB/SipC translocon ≠ effector).

**T4SS:** Christie et al. 2014 BBA 10.1016/j.bbamcr.2013.12.019; Li, Hu & Christie 2019
Microbiol Spectr 10.1128/microbiolspec.PSIB-0012-2018; Vergunst et al. 2003 Plant Physiol
10.1104/pp.103.029223 (C-terminal signal); Rambow-Larsen & Weiss 2004 J Bacteriol
10.1128/JB.186.1.43-50.2004 (Ptl two-step exception).

**T5SS:** Leo, Grin & Linke 2012 Phil Trans R Soc B 10.1098/rstb.2011.0208 (T5a–e
classification, ESPR); Doyle & Bernstein 2021 Mol Cell 10.1016/j.molcel.2021.02.023
(T5a passenger via BamA); Jacob-Dubuisson et al. 2013 Res Microbiol
10.1016/j.resmic.2013.03.009 (T5b TPS); Cotter, Surana & St Geme 2005 Trends Microbiol
10.1016/j.tim.2005.03.004 (T5c trimeric).

**T6SS:** Cherrak et al. 2019 Microbiol Spectr 10.1128/microbiolspec.PSIB-0031-2019;
Durand et al. 2014 Trends Microbiol 10.1016/j.tim.2014.06.004; Russell et al. 2011 Nature
10.1038/nature10244 (effectors reach prey periplasm via apparatus, not donor periplasm).
