function createOut(data, cutoff) {
    const pairs = data.pairs;
    const size = data.max_count;
    const genomes = data.genomes;

    let output = {};
    let links = [];
    
    for (let pair in pairs) {
        const splitPair = pair.split(' ');
        const source = splitPair[0];
        const target = splitPair[1];
        const count = pairs[pair];
        if (count/size >= cutoff) {
            const linkEntry = {
                'source': source,
                'target': target
            };
            links.push(linkEntry);
        }
    }

    // for every genome in genomes
    const nodes = genomes.map(function(genome, i) {
        return {
            'id': genome,
            'group': i+1
        };
    });

    output.links = links;
    output.nodes = nodes;

    return output;
}