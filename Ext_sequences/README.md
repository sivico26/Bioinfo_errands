# Documentación `ext.sh`

`ext.sh` extrae secuencias de una secuencia referencia (genoma), basado en un alineamiento hecho a la referencia u otra secuencia que será alineada para obtener dicho alineamiento.

La secuencia a extraer será de la *referencia*. Algunas opciones permiten modificar qué secuencias se extraen como se explica a continuación.

### Instalación y dependencias

Al ser un script de `sh` ,`ext.sh` no requiere instalación, aunque debe ser utilizado en una terminal Unix. Sin embargo si requiere [`magicblast`](https://ncbi.github.io/magicblast/), [`bioawk`](https://github.com/lh3/bioawk) y [`bedtools`](https://bedtools.readthedocs.io/en/latest/) para poder ejecutarse. Estas pueden instalarse fácilmente usando conda.

Para crear un ambiente nuevo llamado *Ext* donde instalar las dependencias:

```bash
conda create -n Ext -c bioconda magicblast bioawk bedtools
```

Si se quiere instalar las dependencias en un ambiente existente llamado *my_env*

```bash
## Se activa el ambiente y se instalan las dependencias.
conda activate my_env
conda install -c bioconda magicblast bioawk bedtools

## Se instalan las dependencias sin activar el ambiente
conda install -n my_env -c bioconda magicblast bioawk bedtools
```

### Opciones

Para ver las opciones usar `ext.sh --help`:

```bash
Extract sequences from a reference given an alignment or another sequence that alings to it.

Usage: ./ext.sh [-q|--query <arg>] [-a|--alignment <arg>] [-r|--reference <arg>] [-o|--out-aln <arg>] [-s|--sequence <arg>] [-e|--(no-)exons] [-n|--(no-)no-splice] [-g|--(no-)gene-exclusion] [-t|--(no-)utrs] [-u|--upbases <arg>] [-d|--downbases <arg>] [-h|--help]

        -q, --query: Path to sequence to align to reference (.fasta). (no default)
        -a, --alignment: Path to alignment (.sam). If provided --query and --out-aln would be ignored. (no default)
        -r, --reference: Path to reference to align to (.fasta). (no default)
        -o, --out-aln: Path to save resulting alignment (.sam). (no default)
        -s, --sequence: Path to save extracted sequences (.fasta). (no default)
        -e, --exons, --no-exons: Split the intervals found in the alignment and extract those sequences. Not compatible with --utrs. [false]. (off by default)
        -n, --no-splice, --no-no-splice: Turn off the splice awareness of the alignment [false]. (off by default)
        -g, --gene-exclusion, --no-gene-exclusion: Exclude the main alignment from the extracted sequences. Use with --utrs to extract only UTRs [false]. (off by default)
        -t, --utrs, --no-utrs: Include upstream and downstream bases of the alignment in the extracted sequences [false]. (off by default)
        -u, --upbases: Number of upstream bases to extract [200]. (default: '200')
        -d, --downbases: Number of downstream bases to extract [200]. (default: '200')
        -h, --help: Prints help

Must have magicblast, bioawk, and bedtools available on the PATH.
```

#### Argumentos requeridos

Hay argumentos que `ext.sh` requiere en todos sus usos: requiere una de dos entradas 

- `-r`, `--reference`: Path a la referencia de la cual se extraerán las secuencias. Se esperan un archivo tipo fasta.
- `-s`, `--sequence`: Path donde escribir las secuencias extraídas. Las secuencias serán escritas en formato fasta.

Adicionalmente, `ext.sh` requiere que se especifique lo siguiente según como se utilice:

- **Si se tiene una secuencia a alinear con la referencia**: En esta modalidad, `ext.sh` alineará la secuencia a la referencia y luego extraerá la secuencia indicada por el alineamiento resultante. `-q`, `--query` indica el Path a la secuencia que se va a alinear; se espera un archivo tipo fasta.  Adicionalmente, se debe especificar `-o`, `--out-aln` , que indica el Path donde se escribirá el resultado del alineamiento (en formato [sam](https://en.wikipedia.org/wiki/SAM_(file_format))).  

  **Nota**: `ext.sh` se ha probado con secuencias individuales. Usar multifasta como *queries* no se ha probado pero se espera se incorpore pronto en la funcionalidad.

- **Si se tiene un alineamiento hecho a la referencia**: En esta modalidad, `ext.sh` omitirá la fase de alineamiento y extraerá la secuencia asociada al alineamiento provisto. `-a`, `--alignment` indica el Path a dicho alineamiento. Se espera que el alineamiento esté en formato [sam](https://en.wikipedia.org/wiki/SAM_(file_format)).

#### Comportamiento predeterminado y opciones

El comportamiento predeterminado de `ext.sh` es hacer un alineamiento *splice-aware* a la referencia (utilizando `magicblast`) y extraer las secuencia de la referencia directamente a dicho alineamiento (definida por el rango desde el principio hasta el fin del alineamiento). 

Un ejemplo de uso:

```bash
## Usando una secuencia (gen.fasta) como input.
./ext.sh -r genoma.fasta -q gen.fasta -o aln_out.sam -s sec_extraida.fasta

## Usando un alineamiento (aln_in.sam) como input.
./ext.sh -r genoma.fasta -a aln_in.sam -s sec_extraida.fasta
```

Varias opciones de `ext.sh` permiten modificar este comportamiento:

- `-n`, `--no-splice`: Deshabilita el *splice-aware* del alineamiento. Útil si no la secuencia usada para el query no es RNA. 

- `-e`,`--exons`:  De ser posible, el alineamiento será dividido en intervalos definidos por el CIGAR y las secuencias extraídas serán las que sean definidas por estos nuevos intervalos. Este escenario puede darse cuando se está alineado todo un gen y el alineamiento al genoma abarca los exones e intrones; utilizar esta opción en este caso extraería, de forma separada, las secuencias de todos los exones encontrados. Esa opción es incompatible con `-t`, `--utrs`.

  ```bash
  ## Usando una secuencia (gen.fasta) como input.
  ./ext.sh --exons -r genoma.fasta -q gen.fasta -o aln_out.sam -s exons.fasta
  ```

- `-t`, `--utrs`: Además de la secuencia especificada por el alineamiento, se incluirán las bases *upstream* y *downstream* en la secuencia extraída. El número de bases predeterminadas adicionales a extraer son 200. Este valor puede modificarse con las opciones `-u`, `--upbases` para las bases *upstream* y `-d`, `--downbases` para las bases *downstream*. 

  ```bash
  ## Usando una secuencia (gen.fasta) como input. Extrae +200 bp upstream y downstream
  ./ext.sh -t -r genoma.fasta -q gen.fasta -o aln_out.sam -s sec_extraida_con_utrs.fasta
  
  ## Usando una secuencia (gen.fasta) como input. Extrae +100 bp upstream y 0 downstream
  ./ext.sh --utrs -u 100 -d 0 -r genoma.fasta -q gen.fasta -o aln_out.sam -s sec_extraida_con_utrs.fasta
  ```

- `-g`, `--gene-exclusion`:  La secuencia extraida no incluirá el intervalo definido por el alineamiento. Esta opción está pensada para ser usada con `-t`, `--utrs` con la cual se extraerían sólo las secuencias *upstream* y *downstream* del intervalo definido por el alineamiento.

  ```bash
  ## Usando una secuencia (gen.fasta) como input. Extrae +200 bp upstream y downstream
  ./ext.sh -g -t -r genoma.fasta -q gen.fasta -o aln_out.sam -s utrs.fasta
  
  ## Usando una secuencia (gen.fasta) como input. Extrae +100 bp upstream y 0 downstream
  ./ext.sh --utrs -u 100 -d 0 -r genoma.fasta -q gen.fasta -o aln_out.sam -s utrs.fasta
  ```

  