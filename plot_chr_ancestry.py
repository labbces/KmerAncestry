#!/usr/bin/env python3
"""Barras por Chr_short/k, cores por Espécie e painéis por Catálogo.

Dependências: pandas, matplotlib, seaborn>=0.12.
Exemplo: python plot_contencao.py dados.tsv --genoma AP85 -o contencao
"""

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def natural_key(text):
    return tuple(int(p) if p.isdigit() else p.casefold()
                 for p in re.split(r"(\d+)", str(text)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arquivo", help="Tabela TSV, ou CSV com --sep ';'.")
    parser.add_argument("--sep", default="\t", help="Separador; padrão: tabulação.")
    parser.add_argument("--genoma", help="Selecionar um valor da coluna Genoma.")
    parser.add_argument("-o", "--saida", default="contencao",
                        help="Prefixo de saída, sem extensão (padrão: contencao).")
    parser.add_argument("--col-wrap", type=int, default=1,
                        help="Painéis por linha (padrão: 1).")
    args = parser.parse_args()
    if args.col_wrap < 1:
        parser.error("--col-wrap deve ser pelo menos 1.")

    sep = "\t" if args.sep == r"\t" else args.sep
    df = pd.read_csv(args.arquivo, sep=sep, dtype=str, keep_default_na=False)
    df.columns = df.columns.str.strip()
    required = ["Genoma", "Chr_short", "Catálogo", "k", "Espécie", "Contenção"]
    missing = set(required) - set(df.columns)
    if missing:
        parser.error(f"Colunas ausentes: {', '.join(sorted(missing))}.")
    for column in required:
        df[column] = df[column].str.strip()
    if df.empty or (df[required] == "").any().any():
        parser.error("A tabela está vazia ou tem campos obrigatórios vazios.")

    genomes = list(df["Genoma"].unique())
    if args.genoma is not None:
        df = df.loc[df["Genoma"] == args.genoma].copy()
        if df.empty:
            parser.error(f"Genoma não encontrado. Disponíveis: {', '.join(genomes)}")
    elif len(genomes) > 1:
        parser.error("Há vários genomas. Selecione com --genoma. Disponíveis: "
                     + ", ".join(genomes))

    try:
        df["Contenção"] = pd.to_numeric(
            df["Contenção"].str.replace(",", ".", regex=False), errors="raise")
        k = pd.to_numeric(df["k"], errors="raise")
        if k.isna().any() or ((k <= 0) | (k % 1 != 0)).any():
            parser.error("k deve conter números inteiros positivos.")
        df["k"] = k.astype(int)
    except ValueError as exc:
        parser.error(f"Valor numérico inválido: {exc}")
    if not df["Contenção"].map(lambda x: float('-inf') < x < float('inf')).all():
        parser.error("Contenção deve conter apenas valores numéricos finitos.")

    keys = ["Catálogo", "Chr_short", "k", "Espécie"]
    if df.duplicated(keys).any():
        parser.error("Há linhas repetidas para Catálogo/Chr_short/k/Espécie. "
                     "Defina como tratar as repetições antes de plotar.")

    # Ordem natural dos cromossomos e numérica dos tamanhos de k-mer.
    pairs = sorted(set(zip(df["Chr_short"], df["k"])),
                   key=lambda p: (natural_key(p[0]), p[1]))
    labels = {p: f"{p[0]} | k={p[1]}" for p in pairs}
    df["Grupo"] = [labels[p] for p in zip(df["Chr_short"], df["k"])]
    species = sorted(df["Espécie"].unique(), key=natural_key)
    catalogs = sorted(df["Catálogo"].unique(), key=natural_key)

    sns.set_theme(style="whitegrid", context="notebook")
    panel_width = max(7, 0.65 * len(pairs))
    g = sns.catplot(
        data=df, kind="bar", x="Grupo", y="Contenção", hue="Espécie",
        col="Catálogo", col_wrap=args.col_wrap,
        order=[labels[p] for p in pairs], hue_order=species, col_order=catalogs,
        errorbar=None, dodge=True, palette="colorblind",
        height=4.5, aspect=panel_width / 4.5, sharex=True, sharey=False,
        legend_out=True,
    )
    g.set_axis_labels("Cromossomo | tamanho do k-mer", "Contenção")
    g.set_titles("Catálogo: {col_name}")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", labelrotation=90, labelbottom=True)
        ax.set_ylim(bottom=0)
        ax.grid(axis="x", visible=False)
    g.figure.suptitle(f"Genoma: {df['Genoma'].iloc[0]}", y=1.02)
    g.figure.subplots_adjust(hspace=0.9, wspace=0.35)

    prefix = Path(args.saida)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        output = Path(f"{prefix}.{extension}")
        g.figure.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Salvo: {output}")
    plt.close(g.figure)


if __name__ == "__main__":
    main()

