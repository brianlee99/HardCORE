#!/usr/bin/python
import itertools, pprint, json, sys, os

def main():
    tree = []
    all_genomes = []
    count = 0
    # get a complete list of genomes for this run
    with open("STRAINS") as f:
        for line in f:
            all_genomes.append(line.strip())
            count += 1

    # initialize empty nodes
    for i in range(count):
        node = {
            "name": i+1,
            "children": []
        }
        tree.append(node)

    os.chdir("Pan_Genome")
    with open("PanGenome.usearch") as f:
        for line in f:
            line = line.strip()

            line_arr = line.strip().split('\t')
            genomes = set()
            for genome in line_arr:
                # get genome name
                genome = genome.split('_')[0][1:]
                genomes.add(genome)
            
            included_genomes = sorted(list(genomes))
            num_genomes = len(genomes)
            # look inside the appropriate node
            node = tree[num_genomes-1]
            children = node['children']

            found = False
            # update existing node
            for child in children:
                if child["name"] == included_genomes:
                    found = True
                    child["size"] += 1
                    break
            # add new node
            if not found:
                # calculate list of missing genomes
                missing_genomes = []

                for genome in all_genomes:
                    if genome not in included_genomes:
                        missing_genomes.append(genome)

                children.append({
                    "name": included_genomes,
                    "missing": missing_genomes,
                    "size": 1
                })

    tree_final = {
        "name": "pan_genome",
        "children": tree
    }
    
    with open("flare.json", "w") as out:
        json.dump(tree_final, out, indent=4)

if __name__ == "__main__":
    folder_name = sys.argv[1]
    main()
