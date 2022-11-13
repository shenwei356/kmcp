
Species information of the ground-truth. The third column is the taxid of the species rank.

    $ csvtk cut -t -f taxid names2.tsv \
        | sed 1d \
        | taxonkit reformat --data-dir taxdump/ -I 1 -f '{s}' -t \
        | tee gs.tsv
    1357705 Pseudoalteromonas phage HM1     1357705
    1357706 Pseudoalteromonas virus HP1     2845950
    1357707 Pseudoalteromonas phage HS1     1357707
    1357708 Pseudoalteromonas phage HS2     1357708
    1357710 Pseudoalteromonas phage HS6     1357710
    1327977 Cellulophaga phage phi38:1      1327977
    1327983 Cellulophaga phage phi18:3      1327983
    1327999 Cellulophaga phage phi38:2      1327999
    1327992 Cellulophaga virus ST   1918720
    1327982 Cellulophaga virus Cba181       1918192
    2886930 Escherichia virus phiX174       10847
    10849   Escherichia virus alpha3        1945589

Viral species in the MetaPhlAn3 database:
    
    $ cat mpa_v30_CHOCOPhlAn_201901_marker_info.txt \
        | grep k__Viruses \
        | rush -k "echo {@^(\d+)}" > mpa3-viral.taxid
    
    $ csvtk uniq -Ht mpa3-viral.taxid > mpa3-viral.taxid.uniq
    
    $ cat mpa3-viral.taxid.uniq \
        | taxonkit reformat --data-dir taxdump/ -I 1 -f '{s}' -t \
        > mpa3-viral.tsv
        
    # how many species are found in the database?
    $ cat mpa3-viral.tsv \
        | csvtk uniq -Ht -f 3 \
        | csvtk grep -Ht -f 3 -P <(cut -f 3 gs.tsv)
    1327983 Cellulophaga phage phi18:3      1327983
    1327977 Cellulophaga phage phi38:1      1327977
    1918720 Cellulophaga virus ST   1918720
    1918192 Cellulophaga virus Cba181       1918192

Viral species in the KMCP database:

    
    $ cut -f 2 kmcp-taxid.map \
        | csvtk uniq -Ht \
        | taxonkit reformat --data-dir taxdump/ -I 1 -f '{s}' -t \
        > kmcp-viral.tsv
    
    # how many species are found in the database?
    $ cat kmcp-viral.tsv \
        | csvtk uniq -Ht -f 3 \
        | csvtk grep -Ht -f 3 -P <(cut -f 3 gs.tsv)
    10847   Escherichia virus phiX174       10847
    10849   Escherichia virus alpha3        1945589
    756282  Cellulophaga virus ST   1918720
    1327977 Cellulophaga phage phi38:1      1327977
    1327982 Cellulophaga virus Cba181       1918192
    1327983 Cellulophaga phage phi18:3      1327983
    1327999 Cellulophaga phage phi38:2      1327999
    1357705 Pseudoalteromonas phage HM1     1357705
    1357706 Pseudoalteromonas virus HP1     2845950
    1357707 Pseudoalteromonas phage HS1     1357707
    1357708 Pseudoalteromonas phage HS2     1357708
    1357710 Pseudoalteromonas phage HS6     1357710
