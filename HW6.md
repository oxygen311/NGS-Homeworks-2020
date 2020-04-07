### Команды
Commands for assemlies:
```
canu -p ecoli -d ecoli_pacbio_10x genomeSize=4.7m stopOnLowCoverage=5 -pacbio-raw pacbio_10x.fq.gz
canu -p ecoli -d ecoli_pacbio_20x genomeSize=4.7m -pacbio-raw pacbio_20x.fq.gz
canu -p ecoli -d ecoli_pacbio_40x genomeSize=4.7m -pacbio-raw pacbio_40x.fq.gz
canu -p ecoli -d ecoli_pacbio_80x genomeSize=4.7m -pacbio-raw pacbio_80x.fq.gz
spades.py -o illumina -t 16 -m 16 --pe1-1 illumina.100x.1.fq.gz --pe1-2 illumina.100x.2.fq.gz
spades.py -o illumina_pacbio_10x -t 16 -m 16 --pe1-1 illumina.100x.1.fq.gz --pe1-2 illumina.100x.2.fq.gz --pacbio pacbio_10x.fq.gz
spades.py -o illumina_pacbio_20x -t 16 -m 16 --pe1-1 illumina.100x.1.fq.gz --pe1-2 illumina.100x.2.fq.gz --pacbio pacbio_20x.fq.gz
spades.py -o illumina_pacbio_40x -t 16 -m 16 --pe1-1 illumina.100x.1.fq.gz --pe1-2 illumina.100x.2.fq.gz --pacbio pacbio_40x.fq.gz
spades.py -o illumina_pacbio_80x -t 16 -m 16 --pe1-1 illumina.100x.1.fq.gz --pe1-2 illumina.100x.2.fq.gz --pacbio pacbio_80x.fq.gz
```

Command for QUAST:
```
quast.py -o quast_out -r reference.fasta pacbio_10x.fasta pacbio_20x.fasta pacbio_40x.fasta pacbio_80x.fasta illumina.fasta illumina_pacbio_10x.fasta illumina_pacbio_20x.fasta illumina_pacbio_40x.fasta illumina_pacbio_80x.fasta
```

### Результаты

Удобнее всего смотреть на результаты в файле `report.html`.

**Genome fraction (%)** получается маленький (80%) только если исопльзовать риды pacbioс покрытием 10х (а реальное покрытие после тримминга вообще 7х). В случае pacbio c покрытием 20х и только ридов иллюмина покрытие ~98%, что уже хорошо, но еще не идеал. А во всех остальных случаях покрытие за 99%. Видимо, ридов pacbio даже с маленьким покрытием хватает, чтобы разрешить порядок в случае сборки ридов илюмина.

**Duplication ratio** несоклько выше у Canu, чем у Spades, собрали лишнего.

По **Largest alignment** и **Total aligned length** уверенно лидируют гибридные сборки.

По **NGA50** и **LGA50** хорошо себя показывают гибридные сборки spades и Canu с 80х покрытием, собирая почти весь геном в один скафолд.

**Но** при это у гибридных сборок Spades определяется 6 **misassemblies**, а у Canu стабильно 8-9 **misassemblies**, тут идеально справилась только чистая illumina. Возможно, это всё таки Structural variation.

На удивление, количество **mismatches**, больше у гибридной сборки Spades, чем у сборки только на ридах pacbio через Canu. Количество ошибок Canu с pacbio даже меньше, чем на ридах Illumina.

### Итог:
Судя по картинкам и предыдущим результатам:
Pacbio с маленьким покрытие - плохо, вообще ничего не соберется.
Отдельно риды illumina - тоже плохо, не понятно как разрешать сложные места типа повторов и большие контиги не собираются.

Начиная с покрытия 40x риды pacbio показывают себя хорошо. Причем при покрытии 80х всё собирается вообще хорошо и с наименьшим количеством мисметчей.

С точки зрения итогового результата, гибридным сборкам хватает и pacbio с покрытием 10x, дальше качество/размер контигов не растет. Видимо, в своей основе spades использует скорее illumina риды, а не pacbio.

Кажется, что идеальный пайплайн - собирать pacbio/nanopore сколько собирается и на сколько хватает денег, а потом полировать ошибки через риды илюмина, но это не точно.