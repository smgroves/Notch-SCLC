import enrichr
import string_db
import networkx as nx
import time

def prune(G_orig, prune_sources = True, prune_sinks = True):
    G = G_orig.copy()
    n = len(G.nodes())
    nold = n + 1

    while (n != nold):
        nold = n
        for tf in list(G.nodes()):
            if prune_sources == True:
                if G.in_degree(tf) == 0:  G.remove_node(tf)
            if prune_sinks == True:
                if G.out_degree(tf) == 0: G.remove_node(tf)
            else:
                if G.in_degree(tf) == 0 and G.out_degree(tf) == 0:G.remove_node(tf)
        n = len(G.nodes())
    return G


def prune_info(G_orig, prune_self_loops=True):
    G = G_orig.copy()
    for tf in list(G.nodes()):
        edges = G.adj[tf]
        for target in list(edges.keys()):
            if tf == target and prune_self_loops:
                G.remove_edge(tf, target)
                continue
            if 'db' not in edges[target]:
                G.remove_edge(tf, target)
            elif len(edges[target]['db']) < 2:
                G.remove_edge(tf, target)
    return prune(G)


def prune_to_chea(G_orig, prune_self_loops=True):
    G = G_orig.copy()
    for tf in list(G.nodes()):
        edges = G.adj[tf]
        for target in list(edges.keys()):
            if tf == target and prune_self_loops:
                G.remove_edge(tf, target)
                continue
            if 'db' in edges[target]:
                if not True in ['ChEA' in i for i in edges[target]['db']]: G.remove_edge(tf, target)
    #                if len(edges[target]['db']) < 2: G.remove_edge(tf, target)
    return prune(G)

def make_network(tfs, outdir):
    tfs = ['ASCL1', 'ASCL2', 'ATF2', 'ATF3', 'AVIL', 'BACH2', 'BCL3',
           'BHLHE22', 'BRCA1', 'BRD4', 'CEBPA', 'CEBPB', 'CEBPD', 'CEBPE',
           'CERS4', 'CLOCK', 'CNOT3', 'CTCF', 'CXXC1', 'DLX5', 'DLX6', 'E2F1',
           'E2F4', 'E2F5', 'E2F6', 'EBF1', 'EGR1', 'ELF3', 'ELK3', 'EOMES',
           'EP300', 'EPAS1', 'ERCC6', 'ETS2', 'ETV4', 'ETV6', 'FLI1', 'FOSL1',
           'FOXA1', 'FOXA2', 'FOXJ3', 'FOXM1', 'FOXN4', 'FOXP3', 'GATA2',
           'GATA3', 'GATA4', 'GATA6', 'GLI1', 'GRIP1', 'HDAC2', 'HES1',
           'HES6', 'HEY1', 'HIVEP3', 'HMG20B', 'HOXA1', 'HOXB5', 'HSF1',
           'INSM1', 'ISL1', 'ISL2', 'JUN', 'JUND', 'KAT2A', 'KAT8', 'KDM2B',
           'KDM6A', 'KLF1', 'KLF12', 'KLF2', 'KLF3', 'KLF4', 'KLF5', 'KLF6',
           'KMT2C', 'LHX2', 'LYAR', 'MAF', 'MAFG', 'MAX', 'MBD1', 'MBD2',
           'MECOM', 'MECP2', 'MGA', 'MITF', 'MLX', 'MLXIP', 'MLXIPL', 'MNT',
           'MTA1', 'MXD1', 'MXD3', 'MXD4', 'MXI1', 'MYB', 'MYBL2', 'MYC',
           'MYCL', 'MYCN', 'MYT1', 'NANOG', 'NCOR1', 'NELFE', 'NFATC4',
           'NFKB2', 'NFKBIZ', 'NFYA', 'NFYB', 'NHLH1', 'NHLH2', 'NKX2-1',
           'NONO', 'NR0B2', 'NR1H3', 'NR2F2', 'ONECUT2', 'OVOL2', 'PATZ1',
           'PAX5', 'PAX9', 'PBX3', 'PITX2', 'PKNOX2', 'PML', 'POU2F1',
           'POU4F1', 'POU4F3', 'POU6F2', 'PPARG', 'PROX1', 'PURG', 'RARG',
           'RBP1', 'RBPJ', 'RCOR2', 'RELA', 'RELB', 'REST', 'RFX7', 'RORC',
           'RUNX1', 'RUNX2', 'SALL2', 'SCRT1', 'SCRT2', 'SIN3A', 'SIX5',
           'SMAD2', 'SMAD3', 'SMAD4', 'SMAD9', 'SMARCA4', 'SOX11', 'SOX18',
           'SOX2', 'SOX3', 'SP100', 'SP110', 'SP5', 'SP6', 'SPI1', 'ST18',
           'STAG1', 'STAT4', 'STAT5A', 'STAT6', 'TAF1', 'TAF7', 'TAF7L',
           'TBP', 'TBPL1', 'TBX10', 'TBX6', 'TCF12', 'TCF15', 'TCF21', 'TCF3',
           'TCF4', 'TEAD2', 'TEAD4', 'TFCP2L1', 'TFE3', 'TGIF2', 'TOX',
           'TOX3', 'TRP63', 'TSHZ2', 'TTF2', 'USF1', 'VEZF1', 'XBP1', 'XRN2',
           'YAP1', 'ZBTB16', 'ZBTB18', 'ZBTB21', 'ZBTB33', 'ZBTB7C', 'ZBTB8B',
           'ZFP217', 'ZFP281', 'ZFP3', 'ZFP41', 'ZFP64', 'ZFX', 'ZKSCAN2',
           'ZMIZ1', 'ZSCAN12']

    G = nx.DiGraph()
    # prelim_G = nx.DiGraph()
    # with open("/Users/sarahmaddox/Dropbox (Vanderbilt)/Quaranta_Lab/SCLC/Network/mothers_network.csv") as infile:
    #     for line in infile:
    #         line = line.strip().split(',')
    #         prelim_G.add_edge(line[0], line[1])

    for tf in tfs: G.add_node(tf)

    for tf in tfs:
        enrichr.build_tf_network(G, tf, tfs)
        time.sleep(1)

    # for edge in prelim_G.edges():
    #     if edge[0] in tfs and edge[1] in tfs:
    #         G.add_edge(edge[0], edge[1])

    outfile = open("_0_network.csv", "w")
    for edge in G.edges(): outfile.write("%s,%s\n" % (edge[0], edge[1]))
    outfile.close()

    Gp = prune(G, prune_sinks=False, prune_sources=False)
    Gp.add_edge('NEUROD1',"MYC",db=["ChEA_2013","ChEA_2015"]) #artificially add NEUROD1 --> MYC connection based on Borromeo et al. 2016

    outfile = open("/Users/smgroves/Dropbox (VU Basic Sciences)/pycharm_workspace/NetworkTools copy/archetype_networking/RPM_thesis/_1_network.csv", "w")
    for edge in G.edges(): outfile.write("%s,%s\n" % (edge[0], edge[1]))
    outfile.close()
    Gpp = prune_info(Gp)

    # add code to keep INS and GCG even though they don't have out-going edges
    Gpp = prune_to_chea(Gp)

    outfile = open("/Users/smgroves/Dropbox (VU Basic Sciences)/pycharm_workspace/NetworkTools copy/archetype_networking/RPM_thesis/_2_network.csv", "w")
    for edge in Gpp.edges(): outfile.write("%s,%s\n" % (edge[0], edge[1]))
    outfile.close()
