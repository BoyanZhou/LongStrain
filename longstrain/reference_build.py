"""
Build a bowtie reference for a given species from NCBI database or self-made reference
Software dependency
python/3.6.5
bowtie2/2.3.5.1
samtools/1.9

Update log:
Date: 2021/06/30
Author: Boyan Zhou
1. Complete the whole function
"""


import pandas as pd
import os


class SpeciesRef:
    def __init__(self, species_name):
        """
        self.representative_genome: likes {'assembly_accession': 'GCF_009759685.1', 'refseq_category':
        'representative genome', 'taxid': '470', 'species_taxid': '470', 'organism_name': 'Acinetobacter baumannii',
        'infraspecific_name': 'strain=ATCC 19606',
        'ftp_path': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/759/685/GCF_009759685.1_ASM975968v1'}

        self.other_genome: likes {'GCF_000018445.1': {'assembly_accession': 'GCF_000018445.1', 'refseq_category': 'na',
         'taxid': '405416', 'species_taxid': '470', 'organism_name': 'Acinetobacter baumannii ACICU',
         'infraspecific_name': 'strain=ACICU',
         'ftp_path': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/445/GCF_000018445.1_ASM1844v1'}}

        :param species_name: like "Helicobacter pylori"
        """
        self.species_name = species_name
        self.species_name_joint = "_".join(species_name.split(" "))     # like "Helicobacter_pylori"
        self.complete_genome = False
        self.reference_genome = {}
        self.representative_genome = {}
        self.other_genome = {}

    @staticmethod
    def _genome_info(an_assembly_record):
        return {"assembly_accession": list(an_assembly_record["# assembly_accession"])[0],
                "refseq_category": list(an_assembly_record["refseq_category"])[0],
                "taxid": str(int(an_assembly_record["taxid"])),
                "species_taxid": str(int(an_assembly_record["species_taxid"])),
                "organism_name": list(an_assembly_record["organism_name"])[0],
                "infraspecific_name": list(an_assembly_record["infraspecific_name"])[0],
                "ftp_path": list(an_assembly_record["ftp_path"])[0]}

    def get_ref(self, assembly_table, logger):
        """
        :param assembly_table: assembly_complete_genome, a pandas table only includes complete genome
        which row is this species
        :param logger: a logging object
        :return: stores genome information in the Class SpeciesRef
        """
        print(assembly_table)
        print(self.species_name)
        species_index = [i.startswith(self.species_name) for i in assembly_table["organism_name"]]  # row index, bool
        assembly_species = assembly_table[species_index]
        print(assembly_species)
        if assembly_species.shape[0] == 0:
            logger.info(f"There is no {self.species_name} in the database or {self.species_name} "
                        f"doesn't have a complete genome.\n")
        else:
            self.complete_genome = True

            # get reference_genome
            reference_genome_species = assembly_species[assembly_species["refseq_category"] == "reference genome"]
            if reference_genome_species.shape[0] == 0:
                logger.info(f"{self.species_name} doesn't have a reference genome (not mean no complete genome).\n")
            else:
                self.reference_genome.update(SpeciesRef._genome_info(reference_genome_species[0:1]))
                logger.info(f"{self.species_name} has a reference genome.\n")

            # get representative_genome
            representative_genome_species = assembly_species[assembly_species["refseq_category"] ==
                                                             "representative genome"]
            if representative_genome_species.shape[0] == 0:
                logger.info(f"{self.species_name} doesn't have a representative genome.\n")
            else:
                self.representative_genome.update(SpeciesRef._genome_info(representative_genome_species[0:1]))
                logger.info(f"{self.species_name} has a representative genome.\n")

            # get common_genome
            common_genome_species = assembly_species[assembly_species["refseq_category"] == "na"]
            if common_genome_species.shape[0] == 0:
                logger.info(f"{self.species_name} doesn't have a common genome (not mean no complete genome).\n")
            else:
                for index, row in common_genome_species.iterrows():
                    self.other_genome.update({row["# assembly_accession"]:
                                             {"assembly_accession": row["# assembly_accession"],
                                              "refseq_category": row["refseq_category"],
                                              "taxid": str(row["taxid"]),
                                              "species_taxid": str(row["species_taxid"]),
                                              "organism_name": row["organism_name"],
                                              "infraspecific_name": row["infraspecific_name"],
                                              "ftp_path": row["ftp_path"]}})
                logger.info(f"{self.species_name} has common genome.\n")

    @staticmethod
    def _check_and_build(species_name, logger):
        """ check bowtie build and unzipped fna """
        fna_gz_lens = {}
        fna_gz = ""
        for file in os.listdir("./"):
            if file.split("_")[-1] == "genomic.fna.gz":
                fna_gz_lens.update({file: len(file.split("_"))})
        # get zipped genomic file
        if fna_gz_lens != {}:
            fna_gz = min(fna_gz_lens.items(), key=lambda x: x[1])[0]
        else:
            logger.info(f"Reference file is not found for {species_name}! Can't create new bowtie build!\n")

        if ("_".join(species_name.split(" ")) + ".1.bt2") in os.listdir("./"):
            logger.info("Bowtie build exists. Not create new bowtie build.\n")
        else:
            if fna_gz != "":
                logger.info(f"Reference data exists. Building bowtie reference using {fna_gz} ... ...\n")
                os.system(f"bowtie2-build {fna_gz} {'_'.join(species_name.split(' '))}")
        # unzip the genomic fna file
        fas_full_path = "_".join(species_name.split(" ")) + ".fas"
        if fas_full_path in os.listdir("./"):
            logger.info("Fasta reference exists. Not create new fasta reference.\n\n")
        else:
            logger.info(f"Creating new fasta reference using {fna_gz} ... ...\n\n")
            os.system(f"gunzip -c {fna_gz} > {fas_full_path}")

        """ build samtools faidx """
        # but don't know why it does not work
        if not os.path.exists(fas_full_path + ".fai"):
            os.system(f"samtools faidx {fas_full_path}")

    def bowtie_build(self, target_path, logger):
        """
        Build the basic reference for the species
        :param target_path:
        :param logger: a logging object
        :return:
        """
        # check the existence of this species' directory
        if os.path.isdir(os.path.join(target_path, "_".join(self.species_name.split(" ")))):
            logger.info(f"{self.species_name} exists in the path. Not create new directory.\n")
        else:
            # create new directory to store reference
            os.chdir(target_path)
            os.system("mkdir " + "_".join(self.species_name.split(" ")))
        os.chdir(os.path.join(target_path, "_".join(self.species_name.split(" "))))

        # check the existence of this species' fna ref
        if True in [g.endswith("genomic.fna.gz") for g in os.listdir("./")]:
            logger.info("Reference data exists. Don't need download from NCBI.\n")
        else:
            # download from NCBI
            if self.complete_genome:
                if self.reference_genome != {}:
                    logger.info(f"Building reference for {self.species_name} using reference genome ... ...\n")
                    os.system(f"wget -r -np -nd {self.reference_genome['ftp_path']}")

                elif self.representative_genome != {}:
                    logger.info(f"Building reference for {self.species_name} using representative genome ... ...\n")
                    os.system(f"wget -r -np -nd {self.representative_genome['ftp_path']}")

                else:
                    logger.info(f"Building reference for {self.species_name} using common genome ... ...\n")
                    os.system(f"wget -r -np -nd {self.other_genome[sorted(self.other_genome.keys())[0]]['ftp_path']}")
            else:
                logger.info(f"Warning! Did not build reference for {self.species_name} "
                            f"because of no complete genome!\n\n")

        # real build process
        SpeciesRef._check_and_build(self.species_name, logger)

    @staticmethod
    def _strain_build(species_name_joint, strain_info, target_path, logger):
        species_name_split = " ".join(species_name_joint.split("_"))

        strain_name = str(strain_info["infraspecific_name"])
        if strain_name == "nan":
            strain_name = strain_info["assembly_accession"]
        strain_name = strain_name.replace(" ", "_")
        strain_name = strain_name.replace("/", "_")
        strain_name = strain_name.replace(";", "_")
        strain_name = strain_name.replace("(", "_")
        strain_name = strain_name.replace(")", "_")
        if strain_name.startswith("strain="):
            strain_name = strain_name[7:]
        strain_name = species_name_joint + "_" + strain_name

        print(f"Strain name is {strain_name}")
        # make directory for this strain
        os.chdir(target_path)
        os.system(f"mkdir {strain_name}")
        os.chdir(os.path.join(target_path, strain_name))

        # get data from NCBI
        logger.info(f"Building reference for strain {strain_name} ... ...\n")
        os.system(f"wget -r -np -nd {strain_info['ftp_path']}")

        # check and build bowtie reference
        SpeciesRef._check_and_build(species_name_split, logger)

    def strains_build(self, target_path, logger, strains_num=3):
        """
        Build bowtie database for multiple strains
        :param target_path:
        :param logger: a logging object
        :param strains_num: number of strain database that we need to create
        :return:
        """
        # make a directory named by species under target path
        os.chdir(target_path)
        if not os.path.exists(self.species_name_joint):
            os.system(f"mkdir {self.species_name_joint}")
        species_target_path = os.path.join(target_path, self.species_name_joint)
        os.chdir(species_target_path)
        if self.representative_genome != {}:
            # representative genome exists
            logger.info(f"The representative genome of {self.species_name} exists.\n")
            if self.reference_genome != {}:
                # reference genome exists
                logger.info(f"Build database using reference genome.\n")
                SpeciesRef._strain_build(self.species_name_joint, self.reference_genome, species_target_path, logger)
                strains_num -= 1
            # strains_num is the number of strain that we still need
            if len(self.other_genome.keys()) >= strains_num:
                logger.info(f"There enough strains in the database for {self.species_name}.\n")
                for i in range(strains_num):
                    SpeciesRef._strain_build(self.species_name_joint,
                                             self.other_genome[sorted(self.other_genome.keys())[i]],
                                             species_target_path, logger)
            else:
                logger.info(f"Lack {strains_num - len(self.other_genome.keys())} strains in the database for "
                            f"{self.species_name}.\n")
                for i in sorted(self.other_genome.keys()):
                    SpeciesRef._strain_build(self.species_name_joint, self.other_genome[i], species_target_path, logger)
        elif self.reference_genome != {}:
            # basic reference use "reference genome", we can only use other genomes
            if len(self.other_genome.keys()) >= strains_num:
                logger.info(f"There enough strains in the database for {self.species_name}.\n")
                for i in range(strains_num):
                    SpeciesRef._strain_build(self.species_name_joint,
                                             self.other_genome[sorted(self.other_genome.keys())[i]],
                                             species_target_path, logger)
            else:
                logger.info(f"Lack {strains_num - len(self.other_genome.keys())} strains in the database for "
                            f"{self.species_name}.\n")
                for i in sorted(self.other_genome.keys()):
                    SpeciesRef._strain_build(self.species_name_joint, self.other_genome[i], species_target_path, logger)
        else:
            # basic reference use "other genome", we can only use from the second other genomes
            if len(self.other_genome.keys()) - 1 > strains_num:
                logger.info(f"There enough strains in the database for {self.species_name}.\n")
                for i in range(strains_num):
                    SpeciesRef._strain_build(self.species_name_joint,
                                             self.other_genome[sorted(self.other_genome.keys())[i + 1]],
                                             species_target_path, logger)
            elif len(self.other_genome.keys()) - 1 > 0:
                logger.info(f"Lack {strains_num - len(self.other_genome.keys())} strains in the database for "
                            f"{self.species_name}.\n")
                for i in sorted(self.other_genome.keys())[1:]:
                    SpeciesRef._strain_build(self.species_name_joint, self.other_genome[i], species_target_path, logger)
            else:
                logger.info(f"There is no other strains in the database for {self.species_name}. "
                            f"Not build bowtie database.\n")


# get assembly info of NCBI database
# assembly_summary =
# "/gpfs/data/lilab/home/zhoub03/software/kraken2/NCBI_standard/library/bacteria/assembly_summary.txt"
# assembly = pd.read_table(assembly_summary, skiprows=[0])  # skip the first row of annotation
# assembly_complete_genome = assembly[assembly["assembly_level"] == "Complete Genome"]

def ref_build(out_path, target_species_list, assembly_summary_path, logger):
    """
    :param out_path: path containing species folder, "out_path/species_name"
    :param target_species_list: list of target species, ["Bifidobacterium longum", "Bacteroides uniformis"]
    :param assembly_summary_path: Example "/gpfs/data/lilab/home/zhoub03/software/kraken2/NCBI_standard/
    library/bacteria/assembly_summary.txt"
    :param logger: a logging object
    :return: build a database from a list of target species
    """
    # get the table containing "Complete Genome"
    assembly = pd.read_table(assembly_summary_path, skiprows=[0])  # skip the first row of annotation
    assembly_complete_genome = assembly[assembly["assembly_level"] == "Complete Genome"]

    # add bowtie2 to environment variable
    # for each species in the list
    for species_name in target_species_list:
        # which row is this species; row index, bool
        species_index = [i.startswith(species_name) for i in assembly_complete_genome["organism_name"]]
        assembly_species = assembly_complete_genome[species_index]      # assembly records of that species
        species_ref = SpeciesRef(species_name)
        species_ref.get_ref(assembly_species, logger)
        species_ref.bowtie_build(out_path, logger)
