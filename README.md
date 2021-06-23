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
Le reads sono state trimmate e poi analizzate con fastqc, per vedere quanto erano migliorate con il trimming. I risultati si possono vedere [qui.](https://github.com/die-lab/gamergates_project/tree/main/fastqc)
```create the directory for storing trimmed reads
if [ -d trimmomatic ]
  then rm -r trimmomatic
  fi
  
mkdir trimmomatic  
```
```trim
cd fastq/.

for fastq in *.fastq
  do /usr/local/anaconda3/share/trimmomatic-0.39-2/trimmomatic SE -threads 5 -phred33 $fastq ../trimmomatic/${fastq%.fastq}.trim.fastq ILLUMINACLIP:/usr/local/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 2>>../trimmomatic/stats_trimmomatic
  done
  
cd ../.  
```
```fastqc of trimmed reads
cd ~/project
fastqc trimmomatic/*.trim.fastq -o fastqc
```
In effeti le raw reads scaricate da NCBI non erano un granchè. In particolare erano quasi allarmanti i valori di *per base sequence content* e di *per sequence GC content*. Sono invece meno problematici i brutti score ottenuti nelle voci di *overrapresented sequence* e di *adapter content*. Essendo le raw reads ancora da trimmare e quindi con gli adapter, è normale trovare sequenze sovrarappresentate. Come dice pure fastqc nella sua analisi, la più probabile origine di quei brutti valori sono gli adapter, non ancora tagliati via da trimmomatic. Passando ai fastqc delle reads trimmate, i valori sono più buoni. Mentre in generale la voce *per sequence GC content* è migliorata, e sono spariti i picchi che le caratterizzavano (e che le rendevano "brutte"), il valore *per base sequence content* in alcuni casi sembra essere peggiorato. Potrebbe anche essere che sia stato fastqc a riconoscerlo come peggiorato quando non lo è. Graficamente infatti si vede una minore oscillazione dei valori. l'*adapter content* è positivo per tutti i campioni, una volta trimmati, il che significa che gli adapters sono stati riconosciuti e tagliati correttamente. Il *seqeunce duplication levels* da, per tutti i campioni (sia per i fastqc deille raw reads che per quelli trimmati), un warning o un failure. Questo potrebbe derivare dalla presenza di un contaminante o da un bias nel sequenziamento. Sembra che questo errore sia abbastanza comune, soprattutto nelle library da RNA-seq per trascrizione differenziale, ma anche per altri motivi, quindi non ci ho dato tropppo peso. Dei campioni trimmati solamente uno (SRR5909319) mostra una *overrapresented sequence*. Anche in questo caso non mi sono preoccupato, essendo presente in un solo campione su sei. Fastqc non è riuscito ad individuare l'origine di questa sequenza sovrarappresentata, come invece riusciva a fare per gli adapter nell'analisi dei raw reads. Anche in questo caso però credo si tratti di un adapter, che probabilmente non è stato riconosciuto come tale perchè il match tra questo e uno di riferimento stava sotto ad una certa soglia. 

Sample | SRA | RawReads | Trimmed | Surviving |
------ | --- | -------- | ------- | --------- | 
workers,120d | SRR5909317 | 21626186 | 14108660 | 0.65 | 
workers,120d | SRR5909319 | 16669963 | 11135596 | 0.67 |
workers,120d | SRR5909321 | 17764503 | 12087966 | 0.68 |
queens, 120d | SRR5909303 | 22630380 | 15362292 | 0.68 |
queens, 120d | SRR5909301 | 22593965 | 15206110 | 0.67 |
queens, 120d | SRR5909299 | 20493983 | 12306915 | 0.60 |

#### Trinity
Sono passato quindi alla fase di assemblaggio del trascrittoma. Essendoci relativamente poche sequenze, o perlomeno file fastq con una dimensione accettabile, ho fatto partire l'analisi di Trinity senza prima fare il subsampling. Sarebbe stato meglio avere più campioni e di questi prenderne solamente una parte, così da avere delle sequenze più rappresentative. Putroppo però i campioni erano solamente questi qui.
```create the directory for running, and storing, Trinity output.
if [ -d trinity ]
  then rm -r trinity
  fi

mkdir trinity
cd trinity 
```
```automatic trinity
list=$(ls ../trimmomatic/*.fastq)
echo $list > lista
sed -i 's/ /,/g' lista
file=$(cat lista)

Trinity --seqType fq --single $file --CPU 6 --max_memory 20G

#serve per pulire la cartella di trimmomatic da nuovi file.
mkdir trimmed_readcount
cd trimmed_readcount
mv ../../trimmomatic/*.readcount .
cd ..
```
 ##### Checking trinity output
 Prima di tutto è stata fatta una RNA-seq read rapresentation dell'assemblaggio. Le statistiche dell'allineamento, visibili qui, non sono male. Meno del 9% delle reads da cui si è ottenuto l'assemblaggio non mappano sull'assemblaggio stesso. 
 ```
  cd /home/STUDENTI/diego.carli/project/trinity/trinity_out_dir/.
  mkdir checking
  cd checking
  ```
  ```
  bowtie2-build ../Trinity.fasta trinity_fasta
  cat ../../../trimmomatic/*.fastq > all_reads.fastq
  bowtie2 -p 10 -q --no-unal -k 20 -x trinity_fasta -U all_reads.fastq 2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam 
  
  rm all_reads.fastq
  cd ..
 ```
  Poi si fanno un po' di statistiche. Dal file [trinity_check.txt](https://github.com/die-lab/gamergates_project/blob/main/trinity/trinity_out_dir/checking/trinity_check.txt) vengono mostrate alcune statistiche molto utili, come ad esempio l'N50, che raggiunge un valore abbastanza buono, sia considerando tutti i "geni", ovvero i contig trovati, che considerando solamente le isoforme più lunghe. 
 ```
 TrinityStats.pl Trinity.fasta > trinity_check.txt
 ```
 Poi gli assemblaggi sono stati caricati su gVolante  per l'analisi della qualità. LQuesta analisi è stata ripetuta più volte, settando diversamente i parametri.
 Per vedere i risultati andare alla pagina dedicata ed inserire il job ID (questo qui sotto, ma senza "ID").
 
 * ID 202106031822-6R7ADFXGD3F82VXR. sequence type: trans, program: Busco1, Reference_gene_set: Arthropoda
 * ID 202106031921-GV3JR6CZ26AAZKMX. sequence type: trans, program: Busco2-3, Reference_gene_set: Arthropoda

È evidente che l'analisi con Busco,v2-v3 è molto più soddisfacente rispetto a Busco,v1.
Probabilmente dipende dal tipo di set di geni ortologhi nella quale questi programmi fanno la ricerca. Busco,v2-v3 essendo stato rilasciato dopo, avrà probabilmente un dataset più grande e più aggiornato. Con solamente 26 gebi core mancanti secondo l'analisi di gVolante più promettente, penso si possa considera l'assemblaggio come buono. 

#### Isoform redundancy
Ho utilizzato Cd-hit perchè abbiamo visto che non necessariamente l'isoforma più lunga, che trovava lo script in perl di trinity per ridurre le isoforme, era quella più corretta . Vado allora con CD-HIT, che però vuole un file formattato diverso da quello che esce da trinity.
```modifico l'header e metto la sequenza su una sola riga usando lo script di Mariangela
sed -i 's/ len=.\+//g' Trinity.fasta
cp /home/STUDENTI/diego.carli/didattica-main/2_bash_scripting/scripts/oneline_fasta.sh .
sh oneline_fasta.sh Trinity.fasta
```
```
cd-hit-est -i Trinity.fasta -o output.fasta -T 12 -t 1 -c 0.9
```
Non sono sicuro di cosa ci sia nei due file di output che mi ha fatto cd-hit-est, magari chiedi. Inizialmente pensavo che da una parte, in <output.fasta>, ci fossero solamente le reads rappresentative, mentre in <output.fasta.clstr> ci fossero tutti i diversi cluster. Ciò non mi torna però perchè se greppo ">" in entrambi i file mi danno un numero molto diverso di righe. Dovrebbero infatti avere, questi due file, lo stesso numero di righe, con da una parte solamente la sequenza rappresentativa dopo l'header con >, e dall'altra l'header, ossia il numero di cluster, indicato con >, seguito da un numero variabile di righe che mi indicano quali reads sono state abbinate a quel cluster. Riassumendo, pensavo ci fosse una sequenza rappresentativa per ogni cluster, ma invece non è così (il numero di cluster è molto più alto). Il file trinity_results/cd_hit.log qui su github mi da le informazioni sulla corsa di cd-hit-est.


