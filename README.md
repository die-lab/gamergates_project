# gamergates_project
Progetto per il laboratorio di genomica comparata.

L'idea dietro questo progetto sarebbe quello di replicare uno studio condotto su una specie di formiche con struttura sociale divisa in caste. La specie in questione riesce a prolungare la vita della colonia anche alla morte della regina. Alcuni individui, definiti come gamergates, diventano ovodepositori e mostrano i tipici comportamenti da regina, molto differenti da quelli di raccolta delle altre operaie. È riconoscibile una trascrizione differenziale tra individui di caste differenti (gamergates/workers), in particolare di una proteina, la corazonina.

Inoltre, studi di caso e controllo tramite l'iniezione di questa proteina hanno dimostrato che le formiche iniettate manifestano comportamenti da operaia, quale anche una propensione per la ricerca del cibo. La corazonina quindi è in grado di fare lo shift verso la casta operaia. Al contrario, in individui con bassi livelli di corazonina (ma non solo), presentano alti livelli d i un altra proteina, la vitellogenina.

Il discorso è molto più lungo e complicato di così, per cui rimando all'articolo: ***https://doi.org/10.1016/j.cell.2017.07.014***

Quello che ho cercato di fare è un'analisi di trascrizione differenziale tra operaie e regina di Harpegnathos saltator, utilizzando gli stessi dati del paper di cui parlo sopra. 

#### SAMPLE

Sample | SRA |
------ | --- | 
workers,120d | SRR5909317 |
workers,120d | SRR5909319 | 
workers,120d | SRR5909321 | 
queens, 120d | SRR5909303 |
queens, 120d | SRR5909301 | 
queens, 120d | SRR5909299 |

#### Download and check fastq file
Ho scaricato le reads da NCBI prendendo i codici SRA dal paper (quelli in tabella qua sopra). Poi le ho validate utilizzando la cartella che viene scaricata automaticamente con `fastq-dump`. In più per ogni campione è stata fatta l'analisi utilizzando fastqc.
```download
fastq-dump --defline-seq '@$sn[_$rn]/$ri' <SRACODE>
```
```check
vdb-validate ~/ncbi/public/sra/*.sra
```
```create an apposite directory and move all the fastq file there
if [ -d fastq ]
  then rm -r fastq
  fi
  
mkdir fastq  
mv *.fastq fastq/.
```  
```create the directory and perform the fastqc analyses
mkdir fastqc
cd fastqc
fastqc ../fastq/*.fastq -o fastqc
```

#### Trimming
Le reads sono state trimmate e poi analizzate con fastqc, per vedere quanto erano migliorate con il trimming.
```create the directory for storing trimmed reads
if [ -d trimmomatic ]
  then rm -r trimmomatic
  fi
  
mkdir trimmomatic  
```
```trim
cd fastq/.

for fastq in *.fastq
  do /usr/local/anaconda3/share/trimmomatic-0.39-2/trimmomatic SE -threads 5 -phred33 $fastq ../trimmomatic/${fastq%.fastq}.trim.fastq ILLUMINACLIP:/usr/local/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 2>../trimmomatic/stats_trimmomatic
  done
  
cd ../.  
```
```fastqc of trimmed reads
cd ~/project
fastqc trimmomatic/*.trim.fastq -o fastqc
```

