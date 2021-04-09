#!/bin/bash
python -m jcvi.formats.chain fromagp groups.agp draft.asm.fasta groups.asm
liftOver -minMatch=0 draft.gff3 groups.chain asm.gff3 unmapped -gff