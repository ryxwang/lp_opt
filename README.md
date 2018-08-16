# LP_OPT
**Input files**
* Hi-C interaction file form  

15200000  15200000  78.7486823612  
15200000  15210000  34.7780532997  
15210000  15210000  122.873216928   
……  

1st col is the left edge of bin1 at a given resolution.
2nd col is the left edge of bin2 at a given resolution.
3rd col is the interaction between bin1 and bin2, raw or normalized.

* CTCF file form  

30  
37  
43  
44  
48  
…  

A list of CTCF indices which is a subset of 1:n, where n is the size of the Hi-C input matrix

**Output files**
* TAD coordinates  

1.5755e+07	1.6225e+07  
1.6435e+07	1.6895e+07  
1.7405e+07	1.8585e+07  
……  

Start and end coordinates of TADs.

**Usage**
* Calling TADs on a segment of chr21 (Hi-C data from Rao et al. 2014):  

`hic_file = './input/chr21_GM12878_example.txt';`  
`ctcf_file = './input/chr21_GM12878_ctcf_example.txt';`   
`bin_size = 10000;`  
`q=[0.9,0.5,0.5]; %thresholds for each level`  
`[bd_start, bd_end, pval] = tad_call(hic_file, ctcf_file, bin_size, q);`  

* Calling TADs on a segment of chr21 jointly for two cell types:

`hic_file1 = './input/chr21_GM12878_example.txt';`  
`hic_file2 = './input/chr21_HUVEC_example.txt';`  
`hic_files={hic_file1, hic_file2};`  
`ctcf_file1 = './input/chr21_GM12878_ctcf_example.txt';`  
`ctcf_file2 = './input/chr21_HUVEC_ctcf_example.txt';`  
`ctcf_files={ctcf_file1, ctcf_file2};`  
`bin_size = 10000;`  
`q=[0.9,0.5,0.5]; %thresholds for each level`  
`[bd_start, bd_end, pval] = tad_call_combined(hic_files, ctcf_files, bin_size, q);`  



