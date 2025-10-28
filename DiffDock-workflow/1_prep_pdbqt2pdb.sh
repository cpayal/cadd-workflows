#!/bin/bash


## This step could be performed from scratch using pdn-tools to extract relvant chains --> propka to determine the pKa --> tleap to protonate
## and then manually renaming the HIS to HID/HIP/HIE BUT, at position 160 --> he threonine is phosphorylated (TPO), and wasnt recognized by PROPKA
## Hence I am using the pdbqt prepared for autodock to convert to a pdb --> also keeps the 2 methods uniform.
## The pdbqt was prepared with meeko_prepare.py, which automatically detrmines the protonation states, and assigns HIS in the correct protonation states.
prot_prefix=$1
out_prot=$(echo $prot_prefix | awk -F'.' '{print $1}' )
cut -c-66 ${prot_prefix} > ${out_prot}.pdb
