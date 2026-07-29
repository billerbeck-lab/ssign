# Secretion-system substrate localization

How each secretion system ssign handles delivers its substrate.

## The two-family framework

Secretion systems split by *where the substrate ends up*:

**Family A, envelope-crossing.** T1SS, T2SS, T5SS (a/b/c). The mature protein ends
up extracellular or surface-attached.

**Family B, injection.** T3SS, T4SS, T6SS. The effector stays cytoplasmic in the
producing cell until it is transported directly into a host or neighbouring cell; it
never free-floats outside. By strict producer-cell biology these are cytoplasmic,
so localization predictors "should" miss them. However, in practice DeepLocPro calls many
of them extracellular anyway, because it is trained on functional "secreted"
annotations.

## Per-type summary

| Type | Mechanism | Sec signal peptide (SignalP) | Substrate destination | Family |
|---|---|---|---|---|
| **T1SS** | 1-step, ABC/MFP/TolC, no periplasm | **No** (uncleaved C-terminal signal) | Extracellular / surface | A |
| **T2SS** | 2-step: Sec/Tat → periplasm → OM | **Yes** | Extracellular / surface | A |
| **T3SS** | 1-step injectisome, cytoplasm → host cytosol | **No** (N-terminal + chaperone-binding signal, uncleaved) | Host-cell cytoplasm | B |
| **T4SS** | 1-step, VirD4 coupling protein, cytoplasm → recipient | **No** | Host/recipient cytoplasm | B |
| **T5aSS** | Sec → periplasm; β-barrel in OM, passenger to surface | **Yes** | Surface-retained OR released | A |
| **T5bSS** | Two-partner: TpsA + TpsB, both Sec-exported | **Yes** | TpsA released extracellular; TpsB stays OM | A |
| **T5cSS** | Trimeric autotransporter, Sec → periplasm → OM trimer | **Yes** | Surface-retained | A |
| **T6SS** | 1-step contractile injection, cytoplasm → target | **No** | Target-cell cytoplasm or periplasm | B |

## What DeepLocPro & DeepSecE actually predict

| Type | DLP ran on | % called **Extracellular** (≥0.8) | DeepSecE call (where it ran) |
|---|---|---|---|
| T1SS | 15 | **100%** (15/15), mean prob 0.97 | 15/15 → T1SS |
| T2SS | 4 | **100%** (4/4), mean 0.92 | 3/4 → T2SS |
| T3SS | 21 | **71%** (15/21), mean 0.66 | 13/21 → T3SS (we flag/exclude these) |
| T4SS | 1 | 0% (0/1) → Cytoplasmic Membrane | 1/1 → T4SS |
| T6SS | 40 | **68%** (27/40), mean 0.55 | 28/40 → T6SS |

**The load-bearing empirical finding:** DeepLocPro calls a *majority* of T1/T2/T3/T6
effectors "Extracellular", **including the injected Family-B types T3SS and T6SS**.
Biology says injected effectors are cytoplasmic and a localization tool "should" miss
them, but DeepLocPro is trained on database annotations where virulence effectors are
typically labelled secreted, so it flags them extracellular regardless. This means:

- DLP-extracellular is a **usable** signal for T3SS and T6SS effectors, not a
  mismatch. (T4SS is the one type that empirically looks cytoplasmic, but n=1 here,
  so it is unverified and deserves its own whole-genome check.)
- DeepSecE's sequence classifier independently flags T1/T6 well, and T3SS too, but
  T3SS-DSE is the flagellar-false-positive source (DeepSecE has no flagellum class, so
  flagellar proteins, T3SS homologs, funnel into its T3SS bin and over-call it), so
  ssign excludes DeepSecE for T3SS and relies on DLP + proximity.

## Implications for ssign

1. **Localization-based calling works across both families in practice**: DeepLocPro
   calls most effectors extracellular even for the injected types, so the "injected ⇒
   cytoplasmic ⇒ unfindable" worry does not hold for DLP as trained.
2. **SignalP is informative only for Family A** (T2SS, T5SS, Ptl-T4SS). It is correctly
   silent for T1SS (C-terminal signal) and the injected types; there a SignalP-negative
   is expected, not evidence against secretion.
3. **DeepSecE T3SS stays excluded** (flagellar FPs); T3SS leans on DLP + proximity.
   DeepSecE is the right sequence signal for T6SS.
4. **To make these numbers definitive**, re-extract whole-genome DLP/DSE predictions
   per type (e.g. Salmonella LT2 for a clean T3SS check).

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
