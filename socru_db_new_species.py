#!/usr/bin/env python3

from Bio import Entrez
from plumbum import local
import ssl,os,sys,time,hashlib,argparse,logging



ssl._create_default_https_context = ssl._create_unverified_context #it was giving me SSL errors

#Task1 Download complete genome from NCBI for a species

parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="Name of the species in the form Genus_species, for example: 'Salmonella_enterica'")
parser.add_argument("-v", "--verbose",help = "Verbose", action="store_true")
args = parser.parse_args()
if args.verbose:
    logging.basicConfig(level=logging.INFO)
print("Output directory:" + args.species_name)
bacteriaName = args.species_name

try:
    os.mkdir(bacteriaName)
except OSError:
    logging.info("Output dir exists!")
    sys.exit() 
else:
    logging.info("Output dir created!")

retmax=99999
Entrez.email = 'purnima.pachori@quadram.ac.uk'
search_term = bacteriaName + '[orgn] AND complete genome[title]'
handle = Entrez.esearch(db='nucleotide', term=search_term,retmax=retmax)
genome_ids = Entrez.read(handle)['IdList']
#genome_ids = genome_ids[0:3] #Uncomment this for testing purposes only 
#print(len(genome_ids))

for i in genome_ids:
    record = Entrez.efetch(db="nucleotide", property="complete genome",  id=i, rettype="fasta", retmode="text")
    completename = os.path.join(bacteriaName, bacteriaName + "_"+ i + ".fasta")
    logging.info('Writing:{}'.format(completename))
    with open(completename, 'w') as f:
        f.write(record.read())


logging.info ("Sorting the downloaded sequences by their ids")
#Task2 the sequence with lowest accession number to create Socru db
sorted_list=sorted(genome_ids,reverse=False)
#print(sorted_list[0])
socru_db_seq = bacteriaName + "_" +sorted_list[0] + ".fasta"
#print(socru_db_seq)


#Task3 create Socru database
#usage: socru_create [options] output_directory assembly.fasta
socru_db_create = local['socru_create']
(socru_db_create[bacteriaName + "/" + bacteriaName , bacteriaName + "/" + socru_db_seq])()
logging.info ("Socru database created!")

#Task4 run Socru on all the complete genomes against the DB created above
#usage: socru [options] species assembly.fasta

logging.info("Running Socru...")
socru=local['socru']
catBlastTopHits = os.path.join (bacteriaName , "combined_top_blast_hits.txt")
ofh = open(catBlastTopHits, 'w')
for i in genome_ids:
    completename2 = os.path.join(bacteriaName, bacteriaName + "_"+ i + ".fasta")
    blast_hit_out = bacteriaName + "/" + "top_blast_hits_" + i + ".txt"
    (socru['--db_dir',bacteriaName ,bacteriaName, completename2, "--top_blast_hits" , blast_hit_out , "-o",  bacteriaName + "/" + "socru_out_" + i + ".txt" ])()
    logging.info('Writing:{}'.format(catBlastTopHits))
    ifh = open(blast_hit_out,"r")
    ofh.write(ifh.read())
    ifh.close()
ofh.close()

logging.info ("Socru run finished!")


#Task5 run Socru_shrink_db on the combined blast_top_hits output
#usage: socru_shrink_database [options] blast_results input_db output_db
socru_shrink_db = local['socru_shrink_database']
(socru_shrink_db[bacteriaName + "/" + "combined_top_blast_hits.txt" ,  bacteriaName + "/" + bacteriaName, bacteriaName + "/" + "shrunk_" + bacteriaName])()
logging.info ("Socru shrink db created!")

#Task6 Run socru against the shrunk db

logging.info("Running Socru on shrunk db...")
socru=local['socru']
catBlastTopHitsshrunk = os.path.join (bacteriaName , "shrunk_combined_top_blast_hits.txt")
ofhs = open(catBlastTopHitsshrunk, 'w')
for i in genome_ids:
    completename3 = os.path.join(bacteriaName, bacteriaName + "_"+ i + ".fasta")
    blast_hit_out_shrunk = bacteriaName + "/" + "shrunk_top_blast_hits_" + i + ".txt"
    (socru['--db_dir',bacteriaName ,bacteriaName, completename3, "--top_blast_hits" , bacteriaName + "/" + "shrunk_top_blast_hits_" + i + ".txt" , "-o",  bacteriaName + "/" + "shrunk_socru_out_" + i + ".txt" ])()
    logging.info('Writing:{}'.format(catBlastTopHitsshrunk))
    ifhs = open(blast_hit_out_shrunk,"r")
    ofhs.write(ifhs.read())
    ifhs.close()
ofhs.close()
logging.info ("Socru shrink run finished!")

#Task7 Compare md5 of the output (socru and shrunk socru)
for i in genome_ids:
    fnsocru = os.path.join (bacteriaName , "socru_out_" + i + ".txt" )
    fnsocrushrunk = os.path.join (bacteriaName  , "shrunk_socru_out_" + i + ".txt" )
    socrumd5=hashlib.md5(open(fnsocru,'rb').read()).hexdigest()
    shrunksocrumd5=hashlib.md5(open(fnsocrushrunk,'rb').read()).hexdigest()
    if socrumd5 == shrunksocrumd5:
        print ("socru md5: " + socrumd5, "socru_shrunk md5: " + shrunksocrumd5, "*Hurray, no change in result!")
    else:
        print ("socru md5: " + socrumd5, "socru_shrunk md5: " + shrunksocrumd5, "*WARNING: Result is different!")
