import yaml


def parse_histogram(infile):
    data = []

    with open(infile) as inputdata:
        for line in inputdata:
            if line.startswith('#'):
                line = line.strip().split(':')
                if line[0] == "#numreads":
                    numreads = int(line[1])
                elif line[0] == "#mean":
                    mean = float(line[1])
                elif line[0] == "#stdev":
                    stdev = float(line[1])
                else:
                    raise Exception()
                continue

            line = line.strip().split(',')
            i = int(line[0])
            val = float(line[1])
            data.append((i, val))

    return data, mean, stdev, numreads


def merge_histo(indata, merged_data, numreads):
    for (i, val) in indata:
        if not i in merged_data:
            merged_data[i] = 0
        merged_data[i] += (val * numreads)
    return merged_data


def normalize_histo(merged_data, total_reads):
    data = []
    indices = sorted(merged_data.keys())
    for idx in indices:
        value = merged_data[idx]
        value = value / total_reads
        data.append((idx, value))
    return data


def prune_histogram(histogram):
    # towards the tail end, most cells will be 0
    # dividing by total reads will make most of these almost 0
    # remove these
    if not histogram:
        return histogram
    for idx in range(len(histogram) - 1, -1, -1):
        if float(histogram[idx][1]) >= 0.0001:
            break

    histogram = histogram[:idx]

    return histogram


def write_histo_file(data, outfile):
    with open(outfile, 'w') as histo_file:
        for i, val in data:
            histo_file.write("{}\t{}\n".format(i, val))


def write_metadata(mean, stdev, outfile):
    with open(outfile, 'w') as fileoutput:
        yaml.safe_dump({'mean': mean, 'stdev': stdev}, fileoutput)


def merge_histograms(infiles, outfile, metadata):
    merged_data = {}
    total_reads = 0

    means = 0
    stdevs = 0

    if isinstance(infiles, dict):
        infiles = infiles.values()

    # if input is a single file
    if isinstance(infiles, str):
        infiles = [infiles]

    for infile in infiles:
        data, mean, stdev, numreads = parse_histogram(infile)

        merged_data = merge_histo(data, merged_data, numreads)

        total_reads += numreads

        means += (mean * numreads)
        stdevs += (stdev * numreads)

    final_histo = normalize_histo(merged_data, total_reads)
    final_histo = prune_histogram(final_histo)

    mean = means / total_reads
    stdev = stdevs / total_reads

    write_histo_file(final_histo, outfile)

    write_metadata(mean, stdev, metadata)
