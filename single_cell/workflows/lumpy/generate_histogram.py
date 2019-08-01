import pysam
import sys
import logging

def get_flag_mask():
    required = 97
    restricted = 3484
    flag_mask = required | restricted

    return flag_mask


def mean_std(L):
    s = sum(L)
    mean = s / float(len(L))
    sq_sum = 0.0
    for v in L:
        sq_sum += (v - mean) ** 2.0
    var = sq_sum / float(len(L))
    return mean, var ** 0.5


def median(L):
    if len(L) % 2 == 1:
        return L[int(len(L) / 2)]  # cast to int since divisions always return floats in python3
    mid = int(len(L) / 2) - 1
    return (L[mid] + L[mid + 1]) / 2.0


def unscaled_upper_mad(xs):
    """Return a tuple consisting of the median of xs followed by the
    unscaled median absolute deviation of the values in xs that lie
    above the median.
    """
    xs.sort()
    med = median(xs)
    umad = median([x - med for x in xs if x > med])
    return med, umad


def get_read_groups(infile):
    with pysam.AlignmentFile(infile) as samfile:
        header = samfile.header
        readgroups = [val['ID'] for val in header['RG']]
    return readgroups


def read_bam_file(infile, readgroups, N, skip):

    flag_mask = get_flag_mask()
    required = 97

    skip_count = {rg: 0 for rg in readgroups}
    with pysam.AlignmentFile(infile) as samfile:
        data = {v: [] for v in readgroups}
        counts = {v: 0 for v in readgroups}
        reads_per_rg = {v:0 for v in readgroups}
        for i, read in enumerate(samfile.fetch()):
            readgroup = read.get_tag('RG')
            reads_per_rg[readgroup] += 1

            if skip_count[readgroup] < (skip - 1):
                skip_count[readgroup] += 1
                continue

            if all([val >= N for val in counts.values()]):
                continue

            flag = read.flag

            refname = read.reference_id

            mate_refname = read.next_reference_id

            isize = read.template_length

            want = mate_refname == refname and flag & flag_mask == required and isize >= 0
            if want:
                data[readgroup].append(isize)
                counts[readgroup] += 1

    return reads_per_rg, data, counts


def calculate_histogram(readgroups, data, counts, min_elements, mads, read_length, X, reads_per_rg):
    means=[]
    stdevs=[]

    finaldata = {v: None for v in readgroups}
    for readgroup, rgdata in data.items():
        if len(rgdata) < min_elements:
            warn_str = "cannot generate histogram for readgroup {}".format(readgroup)
            logging.getLogger("lumpy.histogram").warn(warn_str)
            continue

        c = counts[readgroup]
        med, umad = unscaled_upper_mad(rgdata)
        upper_cutoff = med + mads * umad

        L = [v for v in rgdata if v < upper_cutoff]
        new_len = len(L)
        removed = c - new_len
        sys.stderr.write("Removed %d outliers with isize >= %d\n" %
                         (removed, upper_cutoff))
        c = new_len

        mean, stdev = mean_std(L)

        start = read_length
        end = int(mean + X * stdev)

        H = [0] * (end - start + 1)
        s = 0

        for x in L:
            if (x >= start) and (x <= end):
                j = int(x - start)
                H[j] = H[int(x - start)] + 1
                s += 1

        result = {i: (H[i], s) for i in range(end - start)}

        finaldata[readgroup] = result

        means.append(mean * reads_per_rg[readgroup])
        stdevs.append(stdev * reads_per_rg[readgroup])

    if means:
        mean = sum(means)/sum(reads_per_rg.values())
    else:
        mean = 0

    if stdevs:
        stdev = sum(stdevs)/sum(reads_per_rg.values())
    else:
        stdev = 0

    return finaldata, mean, stdev


def merge_readgroups(readgroups, histogram_data):
    keys = set()
    for readgroup, rgdata in histogram_data.items():
        if not rgdata:
            continue
        for isize in rgdata:
            keys.add(isize)
    if not keys:
        return
    keys = sorted(keys)

    finalresult = []
    for i in keys:
        j = 0
        k = 0
        for rg in readgroups:
            if not histogram_data[rg]:
                continue
            values = histogram_data[rg].get(i, None)
            if not values:
                continue
            j += values[0]
            k += values[1]

            finalresult.append((i, (float(j) / k)))
    return finalresult

def write_output(data, mean , stdev, outfile, numreads):

    with open(outfile, 'w') as output:
        output.write('#mean:{}\n'.format(mean))
        output.write('#stdev:{}\n'.format(stdev))
        output.write("#numreads:{}\n".format(numreads))
        if data:
            for idx, value in data:
                output.write("{},{}\n".format(idx,value))


def gen_histogram(
        infile, outfile, N=10000, skip=100000,
        min_elements=1000, mads=10, X=4, read_length=101,
):

    readgroups = get_read_groups(infile)

    reads_per_rg, isizes, counts = read_bam_file(infile, readgroups, N, skip)

    histodata, mean, stdev = calculate_histogram(readgroups, isizes, counts, min_elements, mads, read_length, X, reads_per_rg)

    numreads = sum(reads_per_rg.values())
    merged_histo = merge_readgroups(readgroups, histodata)

    write_output(merged_histo, mean, stdev, outfile, numreads)