# Vendor / synthesizer templates

CaptureForge emits a **generic** order package (`baits.fasta` + `baits.csv`). Map it to a specific
synthesizer's template before ordering. CaptureForge does **not** contact any vendor.

Common targets (column-mapping notes to be added as needed):

| Vendor | Product | Notes |
|---|---|---|
| Twist Bioscience | Custom Panel / oligo pool | accepts FASTA + target BED |
| IDT | xGen Custom / Lockdown Probes | 120 nt probes, biotinylated |
| Agilent | SureSelect / custom baits | 120 nt RNA baits |
| Arbor Biosciences | myBaits | 80–120 nt, IVT RNA baits |
| GenScript / Twist | oligo pools | for `ivt_rna` format (T7 + handles) |

The `baits.csv` columns (`bait_id, sequence, copy_number, gene, class, source_region, strand, GC, Tm,
off_target_flag`) cover what every vendor template needs; reshape/rename per the chosen vendor.
