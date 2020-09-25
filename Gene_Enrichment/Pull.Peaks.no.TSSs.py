#need to pull peaks that are not near TSSs (farther away than 500bp)
#can read in TSSs at peaks, all peaks and then remove peaks if in TSS file
#to run script: python3 Pull.Peaks.no.TSSs.py <tss observed overlap file> <all peaks xls file> <output>

import sys


#read in tss observed file
#returns dictionary with key == isoform id and value == [chr number, peak start, peak end]
def read_tss_with_peaks():
    tss_file = sys.argv[1]
    peaks_at_tss = []
    with open(tss_file, 'r') as tss_pos:
        for line in tss_pos:
            if line.startswith("Chr"):
                continue
            else:
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[1])
                peak_end = int(new_line[2])
                list_value = [chr_num, peak_start, peak_end]
                if list_value in peaks_at_tss:
                    continue
                elif list_value not in peaks_at_tss:
                    peaks_at_tss.append(list_value)
    return peaks_at_tss


#read in all peaks
#Returns list of peaks with [chr num, peak start, peak end]
def read_all_peaks():
    peaks_file = sys.argv[2]
    peaks_list = []
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[1])
                peak_end = int(new_line[2])
                list_value = [chr_num, peak_start, peak_end]
                peaks_list.append(list_value)
    return peaks_list


#sorting peaks
#returns list of peaks not found within 500bp of TSSs
def sort():
    all_peaks = read_all_peaks()
    peaks_at_tss = read_tss_with_peaks()
    peaks_not_at_tss = []
    for peak in all_peaks:
        if peak in peaks_at_tss:
            continue
        elif peak not in peaks_at_tss:
            peaks_not_at_tss.append(peak)
    return peaks_not_at_tss

#write output
def write():
    peaks_not_at_tss = sort()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Chr.Num\tPeak.Start\tPeak.End\n"
        out.write(header)
        for peak in peaks_not_at_tss:
            final = "%s\t%s\t%s\n" % (peak[0], str(peak[1]), str(peak[2]))
            out.write(final)



write()
