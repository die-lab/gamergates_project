# gamergates_project
Progetto per il laboratorio di genomica comparata.

L'idea dietro questo progetto sarebbe quello di replicare uno studio condotto su una specie di formiche con struttura sociale divisa in caste. La specie in questione riesce a prolungare la vita della colonia anche alla morte della regina. Alcuni individui, definiti come gamergates, diventano ovodepositori e mostrano i tipici comportamenti da regina, molto differenti da quelli di raccolta delle altre operaie. È riconoscibile una trascrizione differenziale tra individui di caste differenti (gamergates/workers), in particolare di una proteina, la corazonina.

Inoltre, studi di caso e controllo tramite l'iniezione di questa proteina hanno dimostrato che le formiche iniettate manifestano comportamenti da operaia, quale anche una propensione per la ricerca del cibo. La corazonina quindi è in grado di fare lo shift verso la casta operaia. Al contrario, in individui con bassi livelli di corazonina (ma non solo), presentano alti livelli di un altra proteina, la vitellogenina.

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
In effeti le raw reads scaricate da NCBI non erano un granchè. In particolare erano quasi allarmanti i valori di *per base sequence content* e di *per sequence GC content*. Sono invece meno problematici i brutti score ottenuti nelle voci di *overrapresented sequence* e di *adapter content*. Essendo le raw reads ancora da trimmare e quindi con gli adapter, è normale trovare sequenze sovrarappresentate. Come dice pure fastqc nella sua analisi, la più probabile origine di quei brutti valori sono gli adapter, non ancora tagliati via da trimmomatic. Passando ai fastqc delle reads trimmate, i valori sono più buoni. Mentre in generale la voce *per sequence GC content* è migliorata, e sono spariti i picchi che le caratterizzavano (e che le rendevano "brutte"), il valore *per base sequence content* in alcuni casi sembra essere peggiorato. Potrebbe anche essere che sia stato fastqc a riconoscerlo come peggiorato quando non lo è. Graficamente infatti si vede una minore oscillazione dei valori. l'*adapter content* è positivo per tutti i campioni, una volta trimmati, il che significa che gli adapters sono stati riconosciuti e tagliati correttamente. Il *seqeunce duplication levels* da, per tutti i campioni (sia per i fastqc deille raw reads che per quelli trimmati), un warning o un failure. Questo potrebbe derivare dalla presenza di un contaminante o da un bias nel sequenziamento. Sembra che questo errore sia abbastanza comune, soprattutto nelle library da RNA-seq, quindi non ci ho dato tropppo peso. Dei campioni trimmati solamente uno (SRR5909319) mostra una *overrapresented sequence*. Anche in questo caso non mi sono preoccupato, essendo presente in un solo campione su sei. Fastqc non è riuscito ad individuare l'origine di questa sequenza sovrarappresentata, come invece riusciva a fare per gli adapter nell'analisi dei raw reads. Anche in questo caso però credo si tratti di un adapter, che probabilmente non è stato riconosciuto come tale perchè il match tra questo e uno di riferimento stava sotto ad una certa soglia. 


