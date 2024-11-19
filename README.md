# ProSNP-DL

1. Introduction.
A deep learning based pipeline to identify functional non-coding germline variants mediating prostate cancer susceptibility 
This pipeline first generates a sequence-based prostate enhancer deep learning model. By applying the deep learning model to the 20 centeral [-200bp, +200bp] windows overlapping a high-Fst SNP, it evaluates the impact of the alternative allele by its essential window number (EWN), i.e. how many times that alternative allele change an non-enhancer to an enhnancer and vice versa according to the model. If the alternative allele with gained-enhancer effects have EWN ≥ 5, it will be considered as a gained eSNP; if the alternative allele with lost-enhancer effects have EWN ≥ 5, it will be considered as a lost eSNP.

![Fig1_github](https://github.com/user-attachments/assets/5d0a0740-08d4-44a3-8bab-6ac94572553c)


2. Usage
- The "model" folder contains the codes for training and testing. The data for training, validation, and testing can be downloaded from https://drive.google.com/drive/folders/1LV2_0kkPfAI6276tpE8g0GUn_Kw0zpxm?usp=drive_link 
and the trained deep learning model (hdf5 file) can be downloaded from https://drive.google.com/file/d/1ZZmbk_mIp27asQv6L44dlAEN8mgzOb5j/view?usp=drive_link

To run the python code, you need to impport tensorflow and Keras which could be loaded through conda environment or image file. For the training step, simply run: python train.py /path/to/your/input/data/ train_input_file validation_input_file. Similarly, to test the model, run: python test.py /path/to/your/input/data/ test_inpput

- In the folder "scoreSNP", there're three scripts corresponding to three steps. Before running the three steps, we need to partition the SNP files to 100 to 1000 bins, so that we can run the deep learning model to score the mutated sequences of all the sliding windows for sets of SNPs parallelly. Here we used the lost eSNP rs10095018 as an example to show the format of the partioned SNP file, named "1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.bed.bin_1".
  
a. In step1, run: perl step1_batch_extract_slidingWindowOverlap_highFstSNPs_seq.pl 1
where the input argument "1" here is the bin ID, which ranges from 1 to X (X=100 to 1000, as determined by users). We applied nibFrag to extract genomic   sequences of hg19. The commands of nibFrag can be downloaded here: https://drive.google.com/file/d/1gvV8OwgvnSi6CIhAAt31JLgMWUgUMo7W/view?usp=drive_link
Step 1 will then generate the sequences of the overlapping 1-KB windows of the SNPs with mutated sequences as well. For the example SNP rs10095018, the output file "1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.slidingWindowOverlapSNP.DLinput.bin_1" is the input for next step, which is the input of the deep learning enhancer model.

b. In step2, you need to impport tensorflow and Keras which could be loaded through conda environment or image file to run the python code:
python step2_score_permuatedSeq.py 1
where 1 is the bin ID as the input, which could range from 1 to 100 (or 1000).

c. Next, the output of step2 is a column of probablity generated by the model. We name it "DLscore_bin_1" for now. After catting all the 1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.slidingWindowOverlapSNP.DLinput.bin_* and DLscore_bin_* files. In the unix environment, run the following commands:

`paste 1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.slidingWindowOverlapSNP.DLinput DLscore > 1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.slidingWindowOverlapSNP.DLinput.withDLscore`
  
`cut -f2-4 1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.slidingWindowOverlapSNP.DLinput.withDLscore > slidingWindowOverlapHighFstSNP.withDLscore`
  
The file named "slidingWindowOverlapHighFstSNP.withDLscore.example" is the output of example SNP rs10095018 after this step.

d. We then run step3 to get the list of gained and lost eSNPs:
   perl step3_getEssentialMuts.pl

For more details, please refer to our preprint https://www.medrxiv.org/content/10.1101/2024.11.14.24317278v1.full-text
