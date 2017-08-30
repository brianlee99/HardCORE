#!/usr/bin/python
import itertools, pprint, json, sys

def main():
    pairs = {}  # number of times a pair occurs in the same family
    size = 0  # number of pan-genome families
    max_count = 0  # count of the genome with most # of gene families
    genome_family_count = {}
    all_genomes = set()  # sorted list of all genomes
    to_json = {}

    with open("PanGenome.usearch") as f:
        for line in f:
            size += 1
            line = line.strip().split('\t')
            genomes = []
            for genome in line:
                # split to the left of underline
                genome = genome.split('_')[0][1:]
                genomes.append(genome)
                all_genomes.add(genome)

                if genome in genome_family_count.keys():
                    genome_family_count[genome] += 1
                else:
                    genome_family_count[genome] = 1

            genomes.sort()
            comb_gen = itertools.combinations(genomes, 2)
            for pair in comb_gen:
                # e.g. pair: ('Indikan.2539.BMH', 'Infantis.2887.BMH')
                pair = pair[0] + " " + pair[1]
                if pair in pairs:
                    pairs[pair] += 1
                else:
                    pairs.setdefault(pair, 1)

    # iterate over genome_family_count to find highest value
    highestGenome = ''
    for genome, count in genome_family_count.items():
        if count > max_count:
            max_count = count
            highestGenome = genome
    
    all_genomes = sorted(list(all_genomes))
    to_json["size"] = size
    to_json["pairs"] = pairs
    to_json["genomes"] = all_genomes
    to_json["max_count"] = max_count
    
    with open("data_pre.json", "w") as out:
        json.dump(to_json, out, indent=4)

if __name__ == "__main__":
    # folder_name = sys.argv[1]
    main()
