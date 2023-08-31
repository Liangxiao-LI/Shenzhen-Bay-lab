
/* comment header */

title 'my little tittle'

proc means data= ;
run;

proc ds2;
data;

    method init();
        
    end;

enddata;
run;
quit;


proc print data=;
run;









LIBNAME mylib XLSX '/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants.xls';

PROC IMPORT DATAFILE=mylib.'6vxx_variants$'
            OUT='6vxx'
            DBMS=XLSX REPLACE;
RUN;