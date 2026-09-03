# References for the substrate figure (svg/05_substrate_signals.svg)

One section per bar. Each entry is tied to the specific claim it supports, so a
claim can be checked without reading the whole list.

**Provenance note.** The bars were built from literature research whose sources
were not recorded at the time. This list was assembled and verified afterwards
(2026-08-23): every entry below was retrieved from PubMed and the authors,
title, journal, year and DOI were checked against the record. It is the evidence
base for what the figure asserts, not a transcript of what was read while
drawing it. Where the figure makes a claim no single paper supports, that is
said explicitly rather than papered over with a review citation.

Reading order within each section is roughly: what the signal is → where it is →
what is *not* known.

---

## T1SS — HlyA (P08715)

Rating on the bar: *location known, no consensus motif*.

**The signal is C-terminal, uncleaved, and portable**

- Nicaud JM, Mackman N, Gray L, Holland IB (1986). The C-terminal, 23 kDa peptide
  of E. coli haemolysin 2001 contains all the information necessary for its
  secretion by the haemolysin (Hly) export machinery. *FEBS Lett* 204:331-5.
  doi:10.1016/0014-5793(86)80838-9
- Koronakis V, Stanley P, Koronakis E, Hughes C (1992). The HlyB/HlyD-dependent
  secretion of toxins by gram-negative bacteria. *FEMS Microbiol Immunol*
  5:45-53. (PMID 1419114) — states the signal is uncleaved and that secretion has
  no periplasmic intermediate.
- Chervaux C, Sauvonnet N, Le Clainche A, Kenny B, Hung AL, Broome-Smith JK,
  Holland IB (1995). Secretion of active beta-lactamase to the medium mediated by
  the E. coli haemolysin transport pathway. *Mol Gen Genet* 249:237-45.
  doi:10.1007/BF00290371 — the signal moves an unrelated, active enzyme.

**Mapped to the last ~50-60 residues, but with no sequence consensus**

- Stanley P, Koronakis V, Hughes C (1991). Mutational analysis supports a role for
  multiple structural features in the C-terminal secretion signal of E. coli
  haemolysin. *Mol Microbiol* 5:2391-403.
  doi:10.1111/j.1365-2958.1991.tb02085.x — signal within the last 48 aa; three
  structural features, not a motif.
- Kenny B, Chervaux C, Holland IB (1994). Evidence that residues -15 to -46 of the
  haemolysin secretion signal are involved in early steps in secretion, leading to
  recognition of the translocator. *Mol Microbiol* 11:99-109.
  doi:10.1111/j.1365-2958.1994.tb00293.x
- Chervaux C, Holland IB (1996). Random and directed mutagenesis to elucidate the
  functional importance of helix II and F-989 in the C-terminal secretion signal
  of E. coli hemolysin. *J Bacteriol* 178:1232-6.
  doi:10.1128/jb.178.4.1232-1236.1996 — the dissenting result: individual
  residues matter irrespective of a defined secondary structure.
- Holland IB, Schmitt L, Young J (2005). Type 1 protein secretion in bacteria, the
  ABC-transporter dependent pathway. *Mol Membr Biol* 22:29-39.
  doi:10.1080/09687860500042013 — review; "(uncleaved), poorly conserved
  secretion signal of approximately 50 residues".
- Jumpertz T, Chervaux C, Racher K, Zouhair M, Blight MA, Holland IB, Schmitt L
  (2010). Mutations affecting the extreme C terminus of E. coli haemolysin A
  reduce haemolytic activity by altering the folding of the toxin. *Microbiology*
  156:2495-2505. doi:10.1099/mic.0.038562-0
- Spitz O, Erenburg IN, Kanonenberg K, Peherstorfer S, Lenders MHH, Reiners J,
  Ma M, Luisi BF, Smits SHJ, Schmitt L (2022). Identity determinants of the
  translocation signal for a type 1 secretion system. *Front Physiol* 12:804646.
  doi:10.3389/fphys.2021.804646 — amphipathic helix 975-987 plus F990; states
  plainly that the RTX domain and the secretion signal are separate elements.

**The RTX repeats are NOT the signal, but they do have secretion roles**

- Pourhassan N Z, Hachani E, Spitz O, Smits SHJ, Schmitt L (2022). Investigations
  on the substrate binding sites of hemolysin B, an ABC transporter, of a type 1
  secretion system. *Front Microbiol* 13:1055032.
  doi:10.3389/fmicb.2022.1055032 — GG repeats bind HlyB's C39-like domain while
  the helix binds the NBD; the ~120-residue spacing between them is conserved and
  functionally required.
- Bumba L, Masin J, Macek P, Wald T, Motlova L, Bibova I, Klimova N, Bednarova L,
  Veverka V, Kachala M, Svergun DI, Barinka C, Sebo P (2016). Calcium-driven
  folding of RTX domain beta-rolls ratchets translocation of RTX proteins through
  type I secretion ducts. *Mol Cell* 62:47-62.
  doi:10.1016/j.molcel.2016.03.018 — the repeats are the motor, not the address.
- Sebo P, Ladant D (1993). Repeat sequences in the Bordetella pertussis adenylate
  cyclase toxin can be recognized as alternative carboxy-proximal secretion
  signals by the E. coli alpha-haemolysin translocator. *Mol Microbiol* 9:999-1009.
  doi:10.1111/j.1365-2958.1993.tb01229.x — the one result that blurs the line: a
  cryptic second signal maps inside the CyaA repeat region.
- Masi M, Wandersman C (2010). Multiple signals direct the assembly and function
  of a type 1 secretion system. *J Bacteriol* 192:3861-9.
  doi:10.1128/JB.00178-10 — HasA, a 188-residue T1SS substrate with no RTX
  repeats at all; also shows signals distributed through the polypeptide.

---

## T2SS — PulA (P07206)

Rating on the bar: *no linear signal* (this refers to step 2, the outer-membrane
step; step 1 has a textbook cleaved signal peptide).

**Two steps, and the T2SS only handles the second one**

- Korotkov KV, Sandkvist M (2019). Architecture, function, and substrates of the
  type II secretion system. *EcoSal Plus* 8(2).
  doi:10.1128/ecosalplus.ESP-0034-2018 — substrates engage the T2SS only after
  they have already reached the periplasm via Sec or Tat.
- Korotkov KV, Sandkvist M, Hol WGJ (2012). The type II secretion system:
  biogenesis, molecular architecture and mechanism. *Nat Rev Microbiol* 10:336-51.
  doi:10.1038/nrmicro2762
- Cianciotto NP, White RC (2017). Expanding role of type II secretion in bacterial
  pathogenesis and beyond. *Infect Immun* 85:e00014-17. doi:10.1128/IAI.00014-17
  — the breadth of the substrate list, which is the reason no common motif is
  expected.

**Step 2 has no linear signal: the central claim of this bar**

- Francetic O, Pugsley AP (2005). Towards the identification of type II secretion
  signals in a nonacylated variant of pullulanase from Klebsiella oxytoca.
  *J Bacteriol* 187:7045-55. doi:10.1128/JB.187.20.7045-7055.2005 — the single
  best citation here: at least three separate regions of PulA carry secretion
  information, none is individually necessary, and the effect depends on context.
  Also the source for PulA being a surface-anchored lipoprotein.
- Bortoli-German I, Brun E, Py B, Chippaux M, Barras F (1994). Periplasmic
  disulphide bond formation is essential for cellulase secretion by the plant
  pathogen Erwinia chrysanthemi. *Mol Microbiol* 11:545-53.
  doi:10.1111/j.1365-2958.1994.tb00335.x — recognition requires a folded
  substrate, i.e. it is conformational.
- D'Enfert C, Pugsley AP (1989). Klebsiella pneumoniae pulS gene encodes an outer
  membrane lipoprotein required for pullulanase secretion. *J Bacteriol*
  171:3673-9. doi:10.1128/jb.171.7.3673-3679.1989

**Step 1: Sec or Tat, and the caveats on that wording**

- Voulhoux R, Ball G, Ize B, Vasil ML, Lazdunski A, Wu LF, Filloux A (2001).
  Involvement of the twin-arginine translocation system in protein secretion via
  the type II pathway. *EMBO J* 20:6735-41. doi:10.1093/emboj/20.23.6735 — the
  paper that established Tat as a second route into the secreton.
- Rossier O, Cianciotto NP (2005). The Legionella pneumophila tatB gene
  facilitates secretion of phospholipase C, growth under iron-limiting
  conditions, and intracellular infection. *Infect Immun* 73:2020-32.
  doi:10.1128/IAI.73.4.2020-2032.2005 — the Tat-dependent Legionella substrate.
- DebRoy S, Dao J, Soderberg M, Rossier O, Cianciotto NP (2006). Legionella
  pneumophila type II secretome reveals unique exoproteins and a chitinase that
  promotes bacterial persistence in the lung. *PNAS* 103:19146-51.
  doi:10.1073/pnas.0608279103 — the size of the Lsp substrate set, against which
  the Tat-dependent fraction is small.
- Ferrandez Y, Condemine G (2008). Novel mechanism of outer membrane targeting of
  proteins in Gram-negative bacteria. *Mol Microbiol* 69:1349-57.
  doi:10.1111/j.1365-2958.2008.06366.x — PnlH, a NON-cleavable Tat signal anchor;
  the exception to "the signal peptide is always cleaved".
- Rondelet A, Condemine G (2012). SurA is involved in the targeting to the outer
  membrane of a Tat signal sequence-anchored protein. *J Bacteriol* 194:6131-42.
  doi:10.1128/JB.01419-12
- Tullman-Ercek D, DeLisa MP, Kawarasaki Y, Iranpour P, Ribnicky B, Palmer T,
  Georgiou G (2007). Export pathway selectivity of E. coli twin arginine
  translocation signal peptides. *J Biol Chem* 282:8309-16.
  doi:10.1074/jbc.M610507200 — 16 of 27 Tat signal peptides also work through
  Sec, so a twin-arginine motif is suggestive, not proof.

**Not citable, and deliberately so:** no systematic Sec:Tat census of T2SS
substrates has been published. Do not quote a percentage.

---

## T3SS — YopE (P31493)

Rating on the bar: *location known, no consensus motif*.

**N-terminal, uncleaved, no consensus**

- Michiels T, Wattiau P, Brasseur R, Ruysschaert JM, Cornelis G (1990). Secretion
  of Yop proteins by Yersiniae. *Infect Immun* 58:2840-9.
  doi:10.1128/iai.58.9.2840-2849.1990 — no cleaved signal sequence and no
  C-terminal domain.
- Michiels T, Cornelis GR (1991). Secretion of hybrid proteins by the Yersinia Yop
  export system. *J Bacteriol* 173:1677-85. doi:10.1128/jb.173.5.1677-1685.1991 —
  the signal lies in the N-terminal 98 residues of YopE; the authors conclude it
  is "conformational rather than sequential".
- Lloyd SA, Sjostrom M, Andersson S, Wolf-Watz H (2002). Molecular
  characterization of type III secretion signals via analysis of synthetic
  N-terminal amino acid sequences. *Mol Microbiol* 43:51-9.
  doi:10.1046/j.1365-2958.2002.02738.x — amphipathicity, not sequence, predicts a
  working signal.

**Bipartite: the chaperone-binding domain is a second, independent signal**

- Wattiau P, Bernier B, Deslee P, Michiels T, Cornelis GR (1994). Individual
  chaperones required for Yop secretion by Yersinia. *PNAS* 91:10493-7.
  doi:10.1073/pnas.91.22.10493
- Cheng LW, Anderson DM, Schneewind O (1997). Two independent type III secretion
  mechanisms for YopE in Yersinia enterocolitica. *Mol Microbiol* 24:757-65.
  doi:10.1046/j.1365-2958.1997.3831750.x — the source for drawing two blocks: a
  signal in residues 1-15 and a SycE-dependent one in residues 15-100.
- Cheng LW, Schneewind O (1999). Yersinia enterocolitica type III secretion: on
  the role of SycE in targeting YopE into HeLa cells. *J Biol Chem* 274:22102-8.
  doi:10.1074/jbc.274.31.22102

**The mRNA-signal question, which the figure deliberately leaves open**

- Anderson DM, Schneewind O (1997). A mRNA signal for the type III secretion of
  Yop proteins by Yersinia enterocolitica. *Science* 278:1140-3.
  doi:10.1126/science.278.5340.1140 — the original proposal.
- Lloyd SA, Norman M, Rosqvist R, Wolf-Watz H (2001). Yersinia YopE is targeted
  for type III secretion by N-terminal, not mRNA, signals. *Mol Microbiol*
  39:520-31. doi:10.1046/j.1365-2958.2001.02271.x — the direct rebuttal.
- Lloyd SA, Forsberg A, Wolf-Watz H, Francis MS (2001). Targeting exported
  substrates to the Yersinia TTSS: different functions for different signals?
  *Trends Microbiol* 9:367-71. doi:10.1016/s0966-842x(01)02100-x — summary of the
  controversy; the reason the figure says "remains unresolved".

**Why effector predictors are classifiers, not regex**

- Arnold R, Brandmaier S, Kleine F, Tischler P, Heinz E, Behrens S, Niinikoski A,
  Mewes HW, Horn M, Rattei T (2009). Sequence-based prediction of type III
  secreted proteins. *PLoS Pathog* 5:e1000376. doi:10.1371/journal.ppat.1000376 —
  ~71% sensitivity, ~85% selectivity from composition features alone.

---

## T4SS — VirE2 (P08062)

Rating on the bar: *location known, no consensus motif*. The bar's note claims
two incompatible flavours of signal; these are the two.

**Agrobacterium: C-terminal, positively charged**

- Vergunst AC, van Lier MCM, den Dulk-Ras A, Hooykaas PJJ (2003). Recognition of
  the Agrobacterium tumefaciens VirE2 translocation signal by the VirB/D4
  transport system does not require VirE1. *Plant Physiol* 133:978-88.
  doi:10.1104/pp.103.029223 — the C-terminal 50 residues of VirE2 are sufficient;
  also the source for the chaperone (VirE1) being dispensable for recognition.
- Vergunst AC, van Lier MCM, den Dulk-Ras A, Grosse Stuve TA, Ouwehand A,
  Hooykaas PJJ (2005). Positive charge is an important feature of the C-terminal
  transport signal of the VirB/D4-translocated proteins of Agrobacterium. *PNAS*
  102:832-7. doi:10.1073/pnas.0406241102 — consensus R-X(7)-R-X-R-X-R.

**Legionella: C-terminal, acidic**

- Nagai H, Cambronne ED, Kagan JC, Amor JC, Kahn RA, Roy CR (2005). A C-terminal
  translocation signal required for Dot/Icm-dependent delivery of the Legionella
  RalF protein to host cells. *PNAS* 102:826-31. doi:10.1073/pnas.0406239101 —
  20 C-terminal residues suffice; hydrophobic residue at position -3/-4.
- Huang L, Boyd D, Amyot WM, Hempstead AD, Luo ZQ, O'Connor TJ, Chen C, Machner M,
  Montminy T, Isberg RR (2011). The E Block motif is associated with Legionella
  pneumophila translocated substrates. *Cell Microbiol* 13:227-45.
  doi:10.1111/j.1462-5822.2010.01531.x — the glutamate-rich E-block.

Note that these two published in the same PNAS issue, back to back, and describe
signals with opposite charge. That is the point the bar's note makes.

---

## T5SS — Hbp (O88093)

Rating on the bar: *well defined*. The only bar with that rating.

- Otto BR, Sijbrandi R, Luirink J, Oudega B, Heddle JG, Mizutani K, Park SY,
  Tame JRH (2005). Crystal structure of hemoglobin protease, a heme binding
  autotransporter protein from pathogenic E. coli. *J Biol Chem* 280:17339-45.
  doi:10.1074/jbc.M412885200 — the passenger structure the bar's proportions
  are drawn from.
- Szabady RL, Peterson JH, Skillman KM, Bernstein HD (2005). An unusual signal
  peptide facilitates late steps in the biogenesis of a bacterial
  autotransporter. *PNAS* 102:221-6. doi:10.1073/pnas.0406055102 — the extended
  signal peptide region found in a subset of autotransporters; the reason the
  signal-peptide block is drawn slightly long.
- Oliver DC, Huang G, Nodel E, Pleasance S, Fernandez RC (2003). A conserved
  region within the Bordetella pertussis autotransporter BrkA is necessary for
  folding of its passenger domain. *Mol Microbiol* 47:1367-83.
  doi:10.1046/j.1365-2958.2003.03377.x — the autochaperone block.
- Jong WSP, ten Hagen-Jongman CM, den Blaauwen T, Slotboom DJ, Tame JRH,
  Wickstrom D, de Gier JW, Otto BR, Luirink J (2007). Limited tolerance towards
  folded elements during secretion of the autotransporter Hbp. *Mol Microbiol*
  63:1524-36. doi:10.1111/j.1365-2958.2007.05605.x
- Ieva R, Bernstein HD (2009). Interaction of an autotransporter passenger domain
  with BamA during its translocation across the bacterial outer membrane. *PNAS*
  106:19120-5. doi:10.1073/pnas.0907912106 — BAM-dependent, C-to-N translocation,
  SurA/Skp in the periplasm. This is the paper behind SurA appearing in the T5SS
  machine column of the other figure.

---

## T6SS — Tse2 (Q9I0E0)

Rating on the bar: *no linear signal*. One bar covers both effector classes; the
specialised case is named in the note rather than given a bar of its own (see
NOTES.md for why).

**Cargo effectors: the class drawn**

- Hood RD, Singh P, Hsu F, Guvener T, Carl MA, Trinidad RRS, Silverman JM,
  Ohlson BB, Hicks KG, Plemel RL, Li M, Schwarz S, Wang WY, Merz AJ, Goodlett DR,
  Mougous JD (2010). A type VI secretion system of Pseudomonas aeruginosa targets
  a toxin to bacteria. *Cell Host Microbe* 7:25-37.
  doi:10.1016/j.chom.2009.12.007 — Tse1-3 identified; Tse2 is the example on the
  bar.
- Silverman JM, Agnello DM, Zheng H, Andrews BT, Li M, Catalano CE, Gonen T,
  Mougous JD (2013). Haemolysin coregulated protein is an exported receptor and
  chaperone of type VI secretion substrates. *Mol Cell* 51:584-93.
  doi:10.1016/j.molcel.2013.07.025 — the mechanism the bar depicts: the effector
  binds the inside of the Hcp ring. Recognition is by a carrier, not a sequence.
- Salomon D, Kinch LN, Trudgian DC, Guo X, Klimko JA, Grishin NV, Mirzaei H,
  Orth K (2014). Marker for type VI secretion system effectors. *PNAS*
  111:9271-6. doi:10.1073/pnas.1406110111 — the MIX motif drawn as the optional
  N-terminal block. Note it is a *marker*, useful for finding effectors, and is
  explicitly shown not to be required for secretion.
- Unterweger D, Kostiuk B, Pukatzki S (2017). Adaptor proteins of type VI
  secretion system effectors. *Trends Microbiol* 25:8-10.
  doi:10.1016/j.tim.2016.10.003 — the three adaptor superfamilies.

**Specialised effectors: the class named in the note**

- Pukatzki S, Ma AT, Revel AT, Sturtevant D, Mekalanos JJ (2007). Type VI
  secretion system translocates a phage tail spike-like protein into target cells
  where it cross-links actin. *PNAS* 104:15508-13.
  doi:10.1073/pnas.0706532104 — VgrG-1 of V. cholerae with its actin-crosslinking
  extension: the first T6SS effector ever demonstrated, and the original
  specialised one.
- Shneider MM, Buth SA, Ho BT, Basler M, Mekalanos JJ, Leiman PG (2013).
  PAAR-repeat proteins sharpen and diversify the type VI secretion system spike.
  *Nature* 500:350-353. doi:10.1038/nature12453 — effector domains fused directly
  to a structural component of the machine.

**No universal signal, and the relative abundance of the two classes**

- Jurenas D, Journet L (2021). Activity, delivery, and diversity of type VI
  secretion effectors. *Mol Microbiol* 115:383-394. doi:10.1111/mmi.14648 — the
  review that supports the strong wording "no universal signal exists".
- Allsopp LP, Bernal P (2023). Killing in the name of: T6SS structure and
  effector diversity. *Microbiology (Reading)* 169:001367.
  doi:10.1099/mic.0.001367 — current terminology ("evolved" VgrG is the older
  term, "specialised" is now standard) and the point that Rhs effectors can fall
  in either class.
- Habich A, Chaves Vargas V, Robinson LA, Allsopp LP, Unterweger D (2025).
  Distribution of the four type VI secretion systems in Pseudomonas aeruginosa
  and classification of their core and accessory effectors. *Nat Commun* 16:888.
  doi:10.1038/s41467-024-54649-5 — the abundance data: cargo effectors dominate
  the conserved core repertoire, specialised effectors dominate the accessory
  repertoire, core:accessory is 2:3, and a genome carries 17-25 T6SS effectors.
  This is the source for the "neither class dominates cleanly" note in NOTES.md.

---

## Attribution

All bibliographic records above were retrieved and verified via PubMed.
