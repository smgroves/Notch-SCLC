import requests
import numpy as np
import pandas as pd

# request.text tsv returns: stringId_A	stringId_B	preferredName_A	preferredName_B	ncbiTaxonId score	nscore	fscore	pscore	ascore	escore	dscore	tscore
# score: combined score - computed by combining the probabilities from the different evidence channels and corrected for the probability of randomly observing an interaction
#n: gene neighborhood - genes occur repeatedly in close neighborhood in genomes
#f: gene fusion: individual gene fusion events per species
#p: phylogenetic profile: presence or absence of linked proteins across species
#a: coexpression: co-expresed in the same or other species
#e: experimental: significant protein interaction datasets from PPI databases
#d: database:: protein interaction groups from curated databases
#t: textmining: extracted from abstracts of literature
def get_string_network_interactions(my_genes, score_threshold = 0.4, species = 9606, caller_id = 'Quaranta network list', verbose = True, scores = ['e','d','t']):
    """
    :param my_genes: list of genes to evaluate
    :param score_threshold: threshold for experimental_score (confidence of interaction)
    :param species: species NCBI taxon identifier (human is 9606, mouse is 10090)
    :param caller_id: your identifier
    :return: list of network connections between inputted genes
    """
    ##################################################################
    ## For the given list of proteins print out only the interactions
    # between these protein which have medium or higher confidence experimental score
    ##################################################################
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    ## Set parameters
    # my_genes = ["CDC42","CDK1","KIF23","PLK1",
    #             "RAC2","RACGAP1","RHOA","RHOB"]
    params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : species, # species NCBI identifier
        "caller_identity" : caller_id # your app name
    }

    ## Call STRING
    response = requests.post(request_url, data=params)

    network = []
    for line in response.text.strip().split("\n"):

        l = line.strip().split("\t")
        p1, p2 = l[2], l[3]

        ## filter the interaction according to experimental score
        experimental_score = float(l[10])
        database_score = float(l[11])
        text_score = float(l[12])
        check = []
        if 'e' in scores: check.append(experimental_score)
        if 'd' in scores: check.append(database_score)
        if 't' in scores: check.append(text_score)

        if np.max(check) > score_threshold:
            ## print
            if verbose:
                print("\t".join([p1, p2, "experimental score: %.3f" % experimental_score,"textmining score: %.3f" % text_score,"database score: %.3f" % database_score]))
            network.append([p1,p2])
    return network

def get_string_network_interactions_full(my_genes, score_threshold = 0.4, species = 9606, caller_id = 'Quaranta network list', verbose = True):
    """
    :param my_genes: list of genes to evaluate
    :param score_threshold: threshold for experimental_score (confidence of interaction)
    :param species: species NCBI taxon identifier (human is 9606, mouse is 10090)
    :param caller_id: your identifier
    :return: list of network connections between inputted genes
    """
    ##################################################################
    ## For the given list of proteins print out only the interactions
    # between these protein which have medium or higher confidence experimental score
    ##################################################################
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv"
    method = "network"

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    ## Set parameters
    # my_genes = ["CDC42","CDK1","KIF23","PLK1",
    #             "RAC2","RACGAP1","RHOA","RHOB"]
    params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : species, # species NCBI identifier
        "caller_identity" : caller_id # your app name
    }

    ## Call STRING
    response = requests.post(request_url, data=params)
    return response.text


def pp_interaction_enrichment(my_genes,species = 9606, caller_id = 'Quaranta network list'):
    ##############################################################
    ## The script prints out the p-value of STRING protein-protein
    ## interaction enrichment method for the given set of proteins
    ##############################################################

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "ppi_enrichment"

    ## Construct the request
    request_url = "/".join([string_api_url, output_format, method])

    ## Set parameters
    # my_genes = ['7227.FBpp0074373', '7227.FBpp0077451', '7227.FBpp0077788',
    #             '7227.FBpp0078993', '7227.FBpp0079060', '7227.FBpp0079448']

    params = {
        "identifiers": "%0d".join(my_genes),  # your proteins
        "species": species,  # species NCBI identifier
        "caller_identity": caller_id # your app name
    }

    ## Call STRING
    response = requests.post(request_url, data=params)

    ## Parse and print the respons Parse and print the responsee
    for line in response.text.strip().split("\n"):
        pvalue = line.split("\t")[5]
        print("P-value:", pvalue)



network_full_tsv = get_string_network_interactions_full(['ASCL1','BCL6','BEND6','MYC','NOTCH1','NOTCH2','NOTCH3','NOTCH4','DLL1','DLL3','DLL4','HES1',
                                           'HES6','JAG1','JAG2','LSD1','MAML','MYCL','MYCN','RBPJ','REST','RIN1','SIRT1'],  score_threshold=.4,
                                          verbose = False)
network_file = open('../data/full_network_table.csv','w')
network_file.write(network_full_tsv)
network_file.close()

network = get_string_network_interactions(['ASCL1','BCL6','BEND6','MYC','NOTCH1','NOTCH2','NOTCH3','NOTCH4','DLL1','DLL3','DLL4','HES1',
                                           'HES6','JAG1','JAG2','LSD1','MAML','MYCL','MYCN','RBPJ','REST','RIN1','SIRT1'],  score_threshold=.4,
                                          scores = ['e','d','t'], verbose = False)

#printing even items because it is doubling list for some reason...
[print(i) for x,i in enumerate(network) if x % 2 == 0]

