CELL1
!pip install -q condacolab
import condacolab
condacolab.install()

CELL2
!mamba create -n prokka_env -c conda-forge -c bioconda prokka -y
!conda run -n prokka_env prokka --version

CELL3
from google.colab import files
uploaded = files.upload()

CELL4
!unzip -o "genome.zip" -d genomes_fasta/
!ls genomes_fasta | head

CELL5
%%bash
mkdir -p prokka_out

for f in "genomes_fasta/genome/"*.fna "genomes_fasta/genome/"*.fasta; do
  [ -f "$f" ] || continue
  base=$(basename "$f" | cut -d. -f1)
  echo "Annotating $base..."
  conda run -n prokka_env prokka --outdir prokka_out/${base} --prefix ${base} --cpus 2 --compliant "$f"
done

ls prokka_out/*/*.gff | head

CELL6
%%bash
zip -r prokka_results.zip prokka_out/

CELL6
from google.colab import files
files.download("prokka_results.zip")