# =============================================================================
# KaMcEs - configuracao unica
#
# Edite SO este arquivo. Todos os estagios o carregam automaticamente; voce
# nunca precisa abrir os scripts para mudar um caminho.
# =============================================================================

# --- modulos do cluster ------------------------------------------------------
MODULE_PYTHON="Python/3.13.5"
MODULE_R="R/4.5.1"
MODULE_SOURMASH="sourmash"

# Ordem importa: os modulos sao carregados na ordem em que aparecem no job, e o
# ULTIMO a mexer no PATH e quem define qual python3 roda. Se o modulo do
# sourmash trouxer o proprio interpretador e ele for de outra versao, o pyarrow
# instalado em py_pkgs nao vai importar. Rode  ./kamces doctor  para ver qual
# python3 cada estagio realmente usa e o que ele consegue importar.
#
# Deixe vazio ("") para nao carregar um modulo.

# --- raizes ------------------------------------------------------------------
BASE="/Storage/data1/yuri.dantas"
KAMCES_DIR="$BASE/scripts/KaMcEs"        # onde este pacote foi descompactado
DATA="$BASE/data"
PY_PKGS="$BASE/scripts/py_pkgs"          # modo antigo: pip install --target

# Ambiente virtual. Se existir, tem precedencia sobre PY_PKGS e o py_pkgs sai
# do PYTHONPATH. E o jeito confiavel de ter scipy e sourmash junto com pyarrow:
# 'pip install --target' nao lida bem com extensoes compiladas.
#   ./kamces setup --venv
KAMCES_VENV="$BASE/scripts/kamces_venv"
R_ENV_DIR="$BASE/scripts/R_env"          # ggplot2, arrow, scattermore...
LOGS="$BASE/scripts/logs_kamces"

# --- parametros biologicos ---------------------------------------------------
# Os cinco valores abaixo aceitam ser sobrepostos pelo AMBIENTE:
#     KS=50 ./kamces catalogs
# roda so o k novo, sem editar arquivo nenhum. Sem o ${VAR:-...}, o config
# sobrescrevia qualquer valor exportado e a variavel de ambiente era ignorada
# em silencio - o pior tipo de "nao funcionou e nao disse por que".
KS="${KS:-19 31 50}"                     # tamanhos de k
SCALEDS="${SCALEDS:-50 100 200 500}"     # valores de scaled
SPECIES="${SPECIES:-officinarum robustum spontaneum barberi}"
WINDOW_SIZE="${WINDOW_SIZE:-10000}"      # janela deslizante (pb)
STEP="${STEP:-100}"                      # passo entre janelas (pb)

# --- estagio 1: catalogos ----------------------------------------------------
CATALOGS_DIR="$DATA/samples/catalogs/final_catalogs"
CATALOGS_PARCIAIS="$DATA/samples/catalogs/parciais"  # fatias antes do merge
CATALOGS_MS_ZIPS="$CATALOGS_DIR/zips_por_acesso"    # um zip por (acesso, k)
CATALOGS_MS_TMP="$CATALOGS_DIR/tmp"                 # FASTQ em transito (1 par)
SRA_BASE_URL="${SRA_BASE_URL:-https://labbces.cena.usp.br/sra_datasets}"
GRUPOS="${GRUPOS:-1}"                    # fatias por especie (1 = sem fatiar)
LOTE_READS="${LOTE_READS:-1000}"         # reads por chamada de add_sequence
# O sliding procura a referencia como {especie}_{k}_scaled{N}.sig (aceita
# tambem _s{N}). Se os nomes ai dentro vierem do merge.sh antigo, eles NAO
# batem - rode  ./kamces refs  para conferir antes de submeter.
# Os FASTQ sao lidos DIRETO DA REDE, descomprimidos em memoria e descartados
# read a read: com FASTQ_CACHE vazio nenhum .fastq.gz toca o disco. Preencha
# so se voce ja tiver os arquivos baixados e quiser reaproveita-los.
FASTQ_CACHE=""                           # opcional: pasta com FASTQ ja baixados
TENTATIVAS_FASTQ="${TENTATIVAS_FASTQ:-3}"  # reabre um FASTQ que caiu no meio
ESPERA_FASTQ="${ESPERA_FASTQ:-15}"         # segundos entre tentativas
# Credenciais do servidor de FASTQ. NUNCA deixe a senha aqui num arquivo que
# vai para o git ou para um zip. Exporte no seu shell:
#     export SRA_USER=srruser
#     export SRA_PASSWORD='...'
SRA_USER="${SRA_USER:-srruser}"

# --- estagio 2: sliding window -----------------------------------------------
GENOME_DIR="$DATA/samples/chr_SP"         # 105 FASTA, um por cromossomo/haplotipo
# Arvore de PARQUET (pos-conversao). Nao aponte para a de CSVs: os CSVs
# convertidos foram apagados e o refazer acharia que tudo esta ausente.
INDEXES_DIR="$BASE/Chr_SP_indexes_parquet"
PART_WINDOWS=20000                       # janelas por parte antes de rotacionar
                                         # (com step=100 sao 10x mais janelas:
                                         #  5000 gerava ~230 partes por arquivo)

# --- estagio 4: graficos -----------------------------------------------------
# Uma pasta por fonte.
PLOTS_SIM_DIR="$DATA/grafs_sliding"

# O espectro de k-mers tem pasta propria: a tabela e as tres figuras nao se
# misturam nem com os catalogos nem com os graficos por cromossomo.
HIST_DIR="$DATA/hist"

# Os tres conjuntos de referencia escrevem em pastas separadas: sao
# experimentos diferentes e os numeros nao se comparam entre si.
#   normal   -> INDEXES_DIR        e PLOTS_SIM_DIR
#   uniform  -> INDEXES_UNIF_DIR   e PLOTS_UNIF_DIR
#   valid    -> INDEXES_VALID_DIR  e PLOTS_VALID_DIR
INDEXES_UNIF_DIR="${INDEXES_DIR}_uniformizado"
INDEXES_VALID_DIR="${INDEXES_DIR}_valid"
INDEXES_VALID_UNIF_DIR="${INDEXES_DIR}_valid_uniformizado"
PLOTS_UNIF_DIR="${PLOTS_SIM_DIR}_uniformizado"
PLOTS_VALID_DIR="${PLOTS_SIM_DIR}_valid"
PLOTS_VALID_UNIF_DIR="${PLOTS_SIM_DIR}_valid_uniformizado"

# CROMOSSOMOS DE VALIDACAO: pares rotulo:caminho, separados por espaco.
#
# Sao cromossomos de identidade CONHECIDA. Rodar o sliding neles e a
# conferencia do metodo: se uma janela de officinarum nao for pintada como
# officinarum, o problema nao esta no cromossomo hibrido.
#
# O rotulo vira PREFIXO do nome na saida (of_Chr01, sp_Chr01) - os dois
# conjuntos tem Chr01..Chr10 com o mesmo nome, e sem o prefixo uma tarefa
# sobrescreveria a outra.
GENOME_VALID_DIRS="of:/Storage/data1/gabriely.santos/DotPlotly/dados/parentais/fasta/of/Chr01-10 sp:/Storage/data1/gabriely.santos/DotPlotly/dados/parentais/fasta/sp/NpX/Chr01-10"

# QUAIS ARQUIVOS DE CADA PASTA. Padroes de glob separados por espaco.
#
# A pasta do spontaneum tem, alem dos 10 cromossomos em .fa, um .fasta do
# Chr01 e dois arquivos de outra nomenclatura (Cr01AD-sp.fa, Cr01-10A-sp.fasta).
# Sem filtro, os quatro entram - e o Chr01 aparece duas vezes (.fa e .fasta),
# o que faria duas tarefas escreverem o MESMO Parquet.
GENOME_PADRAO="*.fa *.fasta *.fa.gz *.fasta.gz"    # pasta principal
GENOME_VALID_PADRAO="Chr*.fa"                      # so os Chr01-10 em .fa

# --- execucao ----------------------------------------------------------------
QUEUE="all.q"
MAX_JOBS=12                              # -tc dos arrays (tarefas simultaneas)
# 12 e o -tc PEDIDO; quem decide quantas de fato rodam juntas e o MEM_FREE
# abaixo. Uma tarefa do sliding carrega os catalogos inteiros na memoria (no
# scaled 200 foram ~20 GB de RSS para 299 milhoes de hashes), entao 12 x 24G e
# o que o escalonador tem de achar livre - ele vai espacar as tarefas conforme
# os nos aguentarem, em vez de deixar 12 brigarem por RAM.
                                         # Com nada consumable na fila, este e
                                         # o UNICO controle real de quantas
                                         # tarefas dividem um no. Dimensione
                                         # por `qhost` (MEMTOT do no) dividido
                                         # pelo pico de uma tarefa.
# --- recursos pedidos ao qsub ------------------------------------------------
# H_VMEM e RLIMIT_AS: espaco de enderecamento VIRTUAL, nao memoria residente.
# numpy/Arrow/sourmash reservam mapeamentos virtuais grandes que nunca sao
# tocados, entao da para bater neste limite usando pouca RAM de verdade.
# Deixe VAZIO ("") para nao pedir limite nenhum.
#
# Nesta fila (all.q), conferido com `qconf -sc` e `qconf -sq all.q`:
#   h_vmem  requestable SIM, consumable NAO, limite da fila INFINITY
# Ou seja: a fila nao impoe teto nenhum e o SGE nao usa este valor para
# decidir quantos jobs cabem por no. O unico teto que existia era o que o
# proprio KaMcEs pedia aqui - era ele que produzia o
# "ArrowMemoryError: malloc of size N failed". Por isso o padrao e vazio.
H_VMEM=""

# MEM_FREE nao limita o processo: diz ao escalonador quanta RAM o job precisa
# encontrar livre no no. E o que espaca os jobs quando H_VMEM esta vazio; sem
# um dos dois, o SGE empilha tarefas ate o OOM killer do Linux agir (SIGKILL
# seco, exit 137, sem traceback). Vazio = nao pede.
#
# ATENCAO: nesta fila mem_free tambem NAO e consumable. Isso quer dizer que
# ele filtra os nos pela RAM livre no instante do agendamento, mas nao
# RESERVA nada - varias tarefas podem ser despachadas para o mesmo no e so
# depois crescerem juntas. Sem consumable, quem de fato controla quantas
# tarefas dividem um no e o MAX_JOBS (-tc) aqui embaixo.
MEM_FREE="24G"

# H_RT: tempo de parede maximo. Com step=100 os jobs ficam ~10x mais longos -
# se a fila tiver um h_rt padrao curto, e ele que mata. Vazio = nao pede.
# Nesta fila h_rt e INFINITY: nao ha o que ajustar.
H_RT=""
ARROW_THREADS=1                          # ver README: aborto 134
FORMATO="parquet"                        # parquet | csv
COMPRESSION="zstd"
