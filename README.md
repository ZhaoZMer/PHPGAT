# PHPGAT
PHPGAT is a phage host prediction tool

## Installation
Download the package by 
```
git clone https://github.com/ZhaoZMer/PHPGAT
cd PHPGAT
```

## Required Dependencies
* [Diamond](https://github.com/bbuchfink/diamond)
* BLAST
* MCL
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Barrnap](https://github.com/tseemann/barrnap)
* [CRT](https://www.room220.com/)
* [PhageRBPdetection](https://github.com/dimiboeckaerts/PhageRBPdetection)
* [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)

### An easiler way to install
*Note*: we suggest you to install all the package using conda (both miniconda and [Anaconda](https://anaconda.org/) are ok)

After cloning this respository, you can use anaconda to install the **PHPGATv2.yaml**. The command is: `conda env create -f PHPGATv2.yaml -n PHPGAT`

### Prepare the databae
Due to the limited size of the GitHub,you can use ncbi-genome-download to download the dataset.

```
ncbi-genome-download -p 4 -F fasta,gff,genbank,protein-fasta  --assembly-accessions allhost.txt bacteria 
```
Move to the prokaryote folder after completing the command.

## Usage

**Options**


      --contigs INPUT_FA
                            input fasta file
      --len MINIMUM_LEN
                            predict only for sequence >= len bp (default 8000)
      --model MODEL (pretrain or retrain)
                            predicting host with pretrained parameters or retrained paramters (default pretrain)
      --topk TOPK_PRED
                            The host prediction with topk score (default 1)
```
conda activate PHPGAT
python run_First.py --contigs test.fasta --len 1000 --model pretrain --topk 20 
```  
+ The input of --contifs is the file name of the input file.
+ The input of --len is is the minimum input length limit, that is, predict only for sequence >= len bp (default 0).
+ The input of --model is predict host with pretrained parameters or retrained paramters (default pretrain).
+ The input of --topk is predict the top k hosts(default 1).

### Output
The format of the output file is a csv file ("final_prediction.csv") which contain the prediction of each phages.

## Extension
If you have more prokaryotic genomes databse,you can place your prokaryotic genomes into *prokaryote/* folder and add an entry of taxonomy information into *dataset/prokaryote.csv*.

If you have more phage-host interactions database,you can place your phage's genomes into thr *phages_train.fasta* file and add an entry of taxonomy information into *dataset/phages.csv*
```
python run_First.py --model retrain
```
## Contact
Please contact ZhaoZM(754506029@qq.com or GitHub Issues) with any questions, concerns or comments.

Thank you for using PHPGAT!
