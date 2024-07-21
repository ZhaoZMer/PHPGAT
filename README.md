# PHPGATv2
PHPGATv2 is a phage host prediction tool

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

After cloning this respository, you can use anaconda to install the **PHPGATv2.yaml**. The command is: `conda env create -f PHPGATv2.yaml -n PHPGATv2`

### Prepare the databae
Due to the limited size of the GitHub,you can use ncbi-genome-download to download the dataset.

```
ncbi-genome-download -p 4 -F fasta,gff,genbank,protein-fasta  --assembly-accessions allhost.txt bacteria 
```
Move to the prokaryote folder after completing the command.

## Usage
```
conda activate PHPGATv2
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

If you have more phage-host interactions database,you can place your phage's genomes into thr *nucl.fasta* file and add an entry of taxonomy information into *dataset/phages.csv*
```
python run_First.py --model retrain
``` 