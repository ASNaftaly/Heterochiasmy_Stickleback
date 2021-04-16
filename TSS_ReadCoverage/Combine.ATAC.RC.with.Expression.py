#combining read coverage data and expression data for ATAC/expression correlations
#will read in output from Condense.Isoforms.by.TSS.py for ATAC read coverage and tsv output from kallisto (will use TPM values here)
#to run script: python3 Combine.ATAC.RC.with.Expression.py <atac read coverage file (averaged)> <tsv expression file (pulling TPM values)> <output file: isoform.id \t average read coverage \t tpm value for isoform
#Author: Alice Naftaly, March 2021


import sys

#read in read coverage file:
#returns dictionary with key == isoform and value == average read coverage in 1kb around TSS
def read_coverage():
    rc_file = sys.argv[1]
    rc_dict = {}
    with open(rc_file, 'r') as rc:
        for line in rc:
            new_line = line.split()
            isoform =  new_line[1]
            ave_rc = float(new_line[2])
            rc_dict.update({isoform:ave_rc})
    return rc_dict


#read in expression file
#returns dictionary with key == isoform and value == tpm
def read_expression():
    exp_file = sys.argv[2]
    tpm_dict = {}
    with open(exp_file, 'r') as exp:
        for line in exp:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                tpm = float(new_line[4])
                tpm_dict.update({isoform:tpm})
    return tpm_dict


#combine tpm and ATAC read coverage
def combine():
    tpm_values = read_expression()
    read_cov = read_coverage()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for iso in read_cov:
            if iso in tpm_values:
                single_read_cov = read_cov[iso]
                single_tpm = tpm_values[iso]
                final = "%s\t%s\t%s\n" % (str(iso), str(single_read_cov), str(single_tpm))
                out.write(final)

combine()
