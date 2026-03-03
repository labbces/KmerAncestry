import os
import subprocess
import requests
import sys
import time
import threading
from tqdm import tqdm


# Amostras
SRR_SAMPLES = {
    'S. officinarum':[
    "SRR32744622_1", "SRR32744622_2", "SRR32744623_1", "SRR32744623_2",
    "SRR32744568_1", "SRR32744568_2", "SRR32744569_1", "SRR32744569_2",
    "SRR32744502_1", "SRR32744502_2", "SRR32744503_1", "SRR32744503_2",
    "SRR32744402_1", "SRR32744402_2", "SRR32744403_1", "SRR32744403_2",
    "SRR32743801_1", "SRR32743801_2", "SRR32740712_1", "SRR32740712_2",
    "SRR32740713_1", "SRR32740713_2", "SRR32740695_1", "SRR32740695_2",
    "SRR32740696_1", "SRR32740696_2", "SRR32673054_1", "SRR32673054_2",
    "SRR32636969_1", "SRR32636969_2", "SRR32634002_1", "SRR32634002_2",
    "SRR32633743_1", "SRR32633743_2", "SRR32632985_1", "SRR32632985_2",
    "SRR32632781_1", "SRR32632781_2", "SRR32632778_1", "SRR32632778_2",
    "SRR32632685_1", "SRR32632685_2", "SRR13165725_1", "SRR13165725_2",
    "SRR13165728_1", "SRR13165728_2"],

    'S. spontaneum':["SRR32744254_1", "SRR32744254_2", "SRR32744255_1", "SRR32744255_2",
    "SRR32744256_1", "SRR32744256_2", "SRR32744241_1", "SRR32744241_2",
    "SRR32744242_1", "SRR32744242_2", "SRR32743804_1", "SRR32743804_2",
    "SRR32743805_1", "SRR32743805_2", "SRR12935194_1", "SRR12935194_2"],

    'S. barberi':["SRR12929232_1", "SRR12929232_2", "SRR32637094_1", "SRR32637094_2",
    "SRR32745130_1", "SRR32745130_2", "SRR32745131_1", "SRR32745131_2",
    "SRR32745132_1", "SRR32745132_2", "SRR32745133_1", "SRR32745133_2",
    "SRR32745134_1", "SRR32745134_2", "SRR32745135_1", "SRR32745135_2",
    "SRR32744643_1", "SRR32744643_2", "SRR32744644_1", "SRR32744644_2",
    "SRR32744645_1", "SRR32744645_2", "SRR32744570_1", "SRR32744570_2",
    "SRR32744571_1", "SRR32744571_2", "SRR32744572_1", "SRR32744572_2",
    "SRR32744526_1", "SRR32744526_2", "SRR32744527_1", "SRR32744527_2",
    "SRR32744528_1", "SRR32744528_2", "SRR32744504_1", "SRR32744504_2",
    "SRR32744505_1", "SRR32744505_2", "SRR32744506_1", "SRR32744506_2",
    "SRR32744468_1", "SRR32744468_2", "SRR32744469_1", "SRR32744469_2",
    "SRR32744470_1", "SRR32744470_2", "SRR32673125_1", "SRR32673125_2",
    "SRR32637548_1", "SRR32637548_2", "SRR32637126_1", "SRR32637126_2"],

    'S. robustum':["SRR12934169_1", "SRR12934169_2", "SRR32763979_1", "SRR32763979_2",
    "SRR32763980_1", "SRR32763980_2", "SRR32744547_1", "SRR32744547_2",
    "SRR32744548_1", "SRR32744548_2", "SRR32744463_1", "SRR32744463_2",
    "SRR32744464_1", "SRR32744464_2", "SRR32744347_1", "SRR32744347_2",
    "SRR32744348_1", "SRR32744348_2", "SRR32744339_1", "SRR32744339_2",
    "SRR32744340_1", "SRR32744340_2", "SRR32744250_1", "SRR32744250_2",
    "SRR32744251_1", "SRR32744251_2", "SRR32744013_1", "SRR32744013_2",
    "SRR32744014_1", "SRR32744014_2", "SRR32744016_1", "SRR32744016_2",
    "SRR32744015_1", "SRR32744015_2", "SRR32743592_1", "SRR32743592_2",
    "SRR32743593_1", "SRR32743593_2", "SRR32743551_1", "SRR32743551_2",
    "SRR32743552_1", "SRR32743552_2"]
}

class ToolSpinner:
    def __init__(self, message):
        self.message = message
        self.stop_event = threading.Event()
        self.thread = threading.Thread(target=self._spin)

    def _spin(self):
        chars = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        i = 0
        while not self.stop_event.is_set():
            sys.stdout.write(f"\r{chars[i % len(chars)]} {self.message}...")
            sys.stdout.flush()
            time.sleep(0.1)
            i += 1

    def start(self): self.thread.start()
    def stop(self):
        self.stop_event.set()
        self.thread.join()
        sys.stdout.write("\r")

def download_file(url, auth, filename):
    response = requests.get(url, auth=auth, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    with open(filename, 'wb') as f, tqdm(
        desc=f"Baixando {filename}",
        total=total_size,
        unit='iB',
        unit_scale=True,
        colour='green'
    ) as bar:
        for data in response.iter_content(chunk_size=8192):
            bar.update(f.write(data))

def run_step(cmd, desc):
    spinner = ToolSpinner(f"{desc}")
    spinner.start()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    spinner.stop()
    
    if result.returncode != 0:
        print(f"Erro em {desc}: {result.stderr[:200]}")
        return False
    print(f"{desc} concluído.")
    return True

def main():
    # Barra de progresso )
    pbar_geral = tqdm(SRR_SAMPLES.items(), desc="Progresso", colour='blue')
    
    for species, srr in pbar_geral:
        fastq_in = f"{srr}.fastq.gz"
        fastq_out = f"{srr}_clean.fastq.gz"
        
        # 1. Download
        if not os.path.exists(fastq_in):
            download_file(f"{BASE_URL}{fastq_in}", (USER, PASS), fastq_in)
        
        # 2. BBDuk
        bb_cmd = (
            f"bbduk.sh -Xmx2g in={fastq_in} out={fastq_out} "
            f"ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen=50 t=2"
        )
        if not run_step(bb_cmd, f"Limpando {species} (BBDuk)"):
            continue

        # 3. KMC Count (Parâmetros aleatórios para testar)
        tmp_dir = f"tmp_{srr}"
        if not os.path.exists(tmp_dir): os.makedirs(tmp_dir)
        kmc_cmd = f"kmc -k31 -m2 -ci1 -cs10000 {fastq_out} {srr}_kmers {tmp_dir}"
        if not run_step(kmc_cmd, f"Contando K-mers {species} (KMC)"):
            continue

        # 4. KMC Dump
        dump_cmd = f"kmc_dump {srr}_kmers {srr}_kmer_catalog.txt"
        run_step(dump_cmd, f"Exportando Catálogo {species}")

    print("\n Finalizado.")

if __name__ == "__main__":
    main()
