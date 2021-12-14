# PP_genRescue
### Code used for "Genomic erosion in a demographically recovered bird species during conservation rescue"
The code is run in four steps:
1. PP_nonWF_burnin_v2_20_07_21.slim: used for generating a burnin stage; produces a treeseq recording without mutations later recapitated and mutated using pyslim (see https://pyslim.readthedocs.io/en/latest/tutorial.html)
2. PP_nonWF_burnin_tree2FullOutput.py: The resulting tree file is recapitated and mutated to incoporate neutral mutations
3. PP_nonWF_step2_v2_20_07_21.slim: takes mutated output from step #2 and runs a shorter burnin for deleterious mutations, it outputs genetic load metrics and slim output files (needed for step 4) at regular basis to track when equilibrium of masked load is reached. 
4. PP_nonWF_step3_v2_20_07_21.slim: Once equilibrium has been reached it uses the file "PP_nonWF_final_demoHistory_goodOne_YEAR2200.txt" to simulate the PP demographic trajectory and conservation managment scenarios (see paper)

### minimal example (unrealistic parameters, e.g. small population and few genes to run faster) 
```slim -d seed=1234 -d Nin=1000 -d g=100 -d geneLength=1000 -d "outPref_in='test_run_burnIn'" -d "chr_genes_in='CDS_chr_prop_len3400_genes13840_totalLen_47Mb.txt'" PP_nonWF_burnin_v2_20_07_21.slim```

```python PP_nonWF_burnin_tree2FullOutput.py --tree test_run_burnIn_genLen1000_genNo100_K1000_seed1234_treeSeq_gen5000.tree --genTime 3.3 --seed 1234 --gen 10 --mID m1000 --U 1 --neutP 0.3 --recRate test_run_burnIn_genLen1000_genNo100_K1000_seed1234_recRates.txt --recPos test_run_burnIn_genLen1000_genNo100_K1000_seed1234_recPos.txt --NeRatio 8 --Outprefix test_run_burnIn```

```slim -d U=1 -d "slim_in='test_run_burnIn_pyslimParams_gensRun5000_genTime3.3_genOut10_mIDm1000_Glen100099.0_U1.0_rho0.001_neutP0.3_pyslimResults_pi2.70e-03_Nc1017_Ne135_mutated_fullOutput.txt'" -d Nin=1000 -d g=100 -d geneLength=1000 -d "outPref_in='test_run_step2'" -d "chr_genes_in='CDS_chr_prop_len3400_genes13840_totalLen_47Mb.txt'" PP_nonWF_step2_v2_20_07_21.slim```

```slim -d "demo_file_path='PP_nonWF_final_demoHistory_goodOne_YEAR2200.txt'" -d "supMode='N'" -d mng_captive=2 -d U=1 -d "slim_in='test_run_step2_genLen1000_genNo100_totalLen0.100098_Nfinal1000_rho0.001_mu5e-06_U1_seed1789525100269_afterTreeSeqBurn_gen2000.slim'" -d Nin=1000 -d g=100 -d geneLength=1000 -d "outPref_in='test_run_step3'" -d "chr_genes_in='CDS_chr_prop_len3400_genes13840_totalLen_47Mb.txt'" PP_nonWF_step3_v2_20_07_21.slim```

Dendencies:
SLiM3, python and python libraries: msprime, pyslim, tskit, numpy, argparse, allel and re



For further information contact Dr Hernan Eduardo Morales Villegas hernanm@sund.ku.dk
