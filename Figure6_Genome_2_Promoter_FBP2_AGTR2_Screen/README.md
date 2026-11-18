# Promoter retention screen — AGTR2 and FBP2

**Question:** is the AGTR2 promoter still present at the AGTR2 locus, and is its absence in *Rhinolophus*
a real deletion or a failure of the method?

The second half of that question is why FBP2 is here. A screen that returns nothing proves nothing unless
the same screen, in the same genomes, returns something when the sequence is there.

## The two controls that make a negative mean something

**Sensitivity — FBP2, searched across the whole genome.** FBP2 is present in every species here, so its
promoter must be found. It is, in all nine bats, at 556–644 bp and 73–78% identity, *including all four
Rhinolophus*. The assay works in exactly the genomes where AGTR2 returns nothing.

**Specificity — cross-class searches.** Bats and birds are separated by roughly 300 million years and
non-coding sequence does not survive that. Every bat promoter is searched against every bird genome and
the reverse; all must be negative. They are. The assay does not return a positive for whatever it is fed.

**Self-hit.** Each query is searched against the genome it claims to come from and must recover ~100% of
itself at ~100% identity. This is what catches a sequence labelled as one species and taken from another —
it caught exactly that in the AGTR2 promoter supplied for this project, which self-hits at 96.7% against
*Myotis myotis* and is therefore a different *Myotis*. That query is kept in the run and reported as
failing, alongside two independently derived *M. myotis* queries that pass at 100%.

## What the analysis does

1. **Finds the locus in each genome.** Human proteins are searched with `tblastn`, chained into collinear
   gene models, and accepted only when a `blastp` search of the model against the complete human RefSeq
   proteome returns the same gene. AGTR2's promoter is sought upstream of AGTR2 where the gene exists, and
   in the PLS3–SLC6A14 interval where it does not. FBP2 is not anchored at all.
2. **Searches with `blastn`.** Collinear hits at one locus are chained; the total number of query bases
   aligned is reported, counted once each, with length-weighted identity. Low-complexity filtering is on:
   a promoter containing a poly-A tract will otherwise match poly-A runs throughout any genome.
3. **Measures a noise floor.** Each interval is shuffled 20 times preserving length and dinucleotide
   composition and rescanned. For a whole-genome search the null is inverted — the *query* is shuffled 10
   times instead, since a genome cannot be. The floor is the greater of the best shuffled score and the
   mean plus three standard deviations.
4. **Measures a detection limit.** The query is mutated from 5% to 60% divergence, embedded in shuffled
   sequence, and rescanned, so a negative states the divergence at which it would still have found
   something.
5. **Checks the assembly.** N content and longest N run of every interval are reported. An absence inside
   a gap is `NOT_READABLE`, not a deletion.

A promoter is called detected only when it clears the floor **and** has E < 1e-5 **and** aligns at least
100 query bases. The floor alone is not enough: with a floor of zero, any chance match would qualify.

## Queries

`queries.tsv` names, for each query, where to search and where the sequence comes from. A sequence is
either supplied as a FASTA or cut from a named genome by the program, which records the coordinates it
used. Nothing is typed, so nothing can be mislabelled without the self-hit control saying so.

`split_promoters.py` turns a plain text file of promoter sequences into one FASTA each, printing the
length and md5 of every file it writes.

## Running it

Edit the paths at the top of `run_PROMOTER.sh`, then:

```
bash run_PROMOTER.sh
python make_promoter_figure.py results
```

The figure is drawn from `PROMOTER.tsv` alone, so it cannot disagree with the table and redrawing needs
no re-run.

Requirements: Python 3.8+, NCBI BLAST+ 2.13+ (`blastn`, `tblastn`, `blastp`, `makeblastdb`, `blastdbcmd`)
on the `PATH`, and `matplotlib` for the figure.

## Output

| File | Contents |
|---|---|
| `PROMOTER.tsv` | one row per species per query: call, aligned bp, identity, bit score, E-value, noise floor, interval and its N content |
| `PROMOTER_blocks.tsv` | every aligned block in query and scaffold coordinates — what a locus diagram is drawn from |
| `PROMOTER_detection_limit.tsv` | the divergence at which each negative would still have been detected |
| `PROMOTER_controls.txt` | self-hit, cross-class specificity, noise floors, detection limits, assembly gaps, query and assembly provenance |
| `PROMOTER_methods.txt` | Methods generated from the parameters of the run |
| `PROMOTER_run_log.txt` | script and query checksums, random seed, every parameter, assembly provenance |
| `PROMOTER_figure.png` / `.pdf` | the figure |

## Notes

- BLAST databases must not sit in a cloud-synced folder; BLAST memory-maps them and placeholder files
  cannot be mapped (`Database memory map file error`).
- The genome assemblies are not redistributed here. Each is identified by its NCBI accession in
  `genomes.tsv` and again in the run log.

## Licence

MIT — see [LICENSE](LICENSE). The licence covers the code, not the third-party assemblies it analyses.
