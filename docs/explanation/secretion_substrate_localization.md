# Secretion-system substrate localization

How each secretion system ssign handles delivers its substrate.

## The two-family framework

Secretion systems split by *where the substrate ends up*:

**Family A, envelope-crossing.** T1SS, T2SS, T5SS (a/b/c). The mature protein ends
up extracellular or surface-attached.

**Family B, injection.** T3SS, T4SS, T6SS. The secreted protein stays cytoplasmic in the
producing cell until it is transported directly into a host or neighbouring cell.

## Per-type summary

| Type | Mechanism | Sec signal peptide (SignalP) | Substrate destination | Family |
|---|---|---|---|---|
| **T1SS** | 1-step, ABC/MFP/TolC, no periplasm | **No** | Extracellular / surface | A |
| **T2SS** | 2-step: Sec/Tat → periplasm → OM | **Yes** | Extracellular / surface | A |
| **T3SS** | 1-step injectisome, cytoplasm → host cytosol | **No** | Host-cell cytoplasm | B |
| **T4SS** | 1-step, VirD4 coupling protein, cytoplasm → recipient | **No** | Host/recipient cytoplasm | B |
| **T5aSS** | Sec → periplasm; β-barrel in OM, passenger to surface | **Yes** | Surface-retained OR released | A |
| **T5bSS** | Two-partner: TpsA + TpsB, both Sec-exported | **Yes** | TpsA released extracellular; TpsB stays OM | A |
| **T5cSS** | Trimeric autotransporter, Sec → periplasm → OM trimer | **Yes** | Surface-retained | A |
| **T6SS** | 1-step contractile injection, cytoplasm → target | **No** | Target-cell cytoplasm or periplasm | B |

## In context of how ssign's works

- **DeepLocPro** (localization) is one of the "is this secreted?" signals. In
  practice it flags many injected Family-B secreted proteins (T3SS, T6SS) as extracellular
  too, not just the envelope-crossing Family-A types, because it is trained on
  functional "secreted" annotations rather than strict producer-cell localization.
  So a localization call is usable across both families.
- **SignalP** is informative for T5SS self-detection, where a
  Sec signal is a positive trigger (see
  [`design_decisions.md` § 3.2](design_decisions.md#32-signalp-is-evidence-only-not-a-trigger)).
- **DeepSecE** is excluded for T3SS (flagellar proteins funnel into its T3SS bin
  and over-call it) but otherwise is used in a similar way to DeepLocPro.

How often each of these signals actually predicts a validated secreted protein, per
system type, is measured in [`benchmarks.md`](../benchmark/benchmarks.md).

## Citations (per type; DOIs with PubMed verification)

**T1SS:** Holland et al. 2016 EcoSal Plus 10.1128/ecosalplus.ESP-0019-2015; Spitz et al.
2022 Front Physiol 10.3389/fphys.2021.804646; Zhang et al. 1995 Biochemistry
10.1021/bi00013a007

**T2SS:** Korotkov & Sandkvist 2019 EcoSal Plus 10.1128/ecosalplus.ESP-0034-2018;
Nivaskumar & Francetic 2014 BBA 10.1016/j.bbamcr.2013.12.020; Ball et al. 2016 Sci Rep
10.1038/srep27675 

**T3SS:** Lara-Tejero & Galán 2019 EcoSal Plus 10.1128/ecosalplus.ESP-0039-2018;
Stebbins & Galán 2001 Nature 10.1038/35102073;
Myeni et al. 2013 PLoS One 10.1371/journal.pone.0060499

**T4SS:** Christie et al. 2014 BBA 10.1016/j.bbamcr.2013.12.019; Li, Hu & Christie 2019
Microbiol Spectr 10.1128/microbiolspec.PSIB-0012-2018; Vergunst et al. 2003 Plant Physiol
10.1104/pp.103.029223; Rambow-Larsen & Weiss 2004 J Bacteriol
10.1128/JB.186.1.43-50.2004 

**T5SS:** Leo, Grin & Linke 2012 Phil Trans R Soc B 10.1098/rstb.2011.0208; Doyle & Bernstein 2021 Mol Cell 10.1016/j.molcel.2021.02.023; Jacob-Dubuisson et al. 2013 Res Microbiol
10.1016/j.resmic.2013.03.009; Cotter, Surana & St Geme 2005 Trends Microbiol
10.1016/j.tim.2005.03.004

**T6SS:** Cherrak et al. 2019 Microbiol Spectr 10.1128/microbiolspec.PSIB-0031-2019;
Durand et al. 2014 Trends Microbiol 10.1016/j.tim.2014.06.004; Russell et al. 2011 Nature
10.1038/nature10244
