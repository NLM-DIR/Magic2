function cleanSpace(_s) {
    _s2 = "" ;
    _iMax = length(_s) ;
    for (_i = 1 ; _i <= _iMax ; _i++)
    {
	_z = substr(_s,_i) ;
	_j = index(_z,"  ") ;
	if (_j == 0)
	{
	    _s2 = _s2 _z ;
	    break ;
	}
	_s2 = _s2 substr(_z,1,_j) ;
	_j++ ;
	while (substr(_s2,_j,1) == " ")
	    _j++ ;
	_i = _i + _j - 2 ;
    } 
    while (substr (_s2,1,1) == " ") _s2 = substr (_s2, 2) ;
    while (substr(_s2,1,1) == "<") # jump the decorations 
    {  
	i = index (_s2, ">") ;  
	if (i < 1) i = 1 ; # avoid loops in case of misconfigured html
	_s2 = substr (_s2, i + 1) ;
	while (substr (_s2,1,1) == " ") _s2 = substr (_s2, 2) ;
    }

    return _s2  ;
}
{ gsub (/\"/,"\\\"",$0); }
{ gsub (/http:\//,"http:\\/",$0) ;}

/pubmed/{
    z = $0 ;
    i = index($0,"/pubmed/") ;
    while (i > 0)
    {
	z = substr (z,i+8);
	j = index(z,"\"") ;
	pm= pm "Reference \"pm" substr(z,1,j-2) "\"\n";  ;
	i = index(z,"/pubmed/") ;
    }
}

/Accession:/{
    i = index($0,"Accession:") ;
    if (i > 0)
    {
	z = substr ($0,i+11);
	i = index(z,"<") ;
	prj= substr(z,1,i-2) ;
    }

    i = index($0,"<div class=\\\"Title\\\"><h2>") ;
    if (i > 0)
    {
	z = substr ($0,i+25);
	i = index(z,"</div>") ;
	title = cleanSpace(substr(z,1,i-1)) ;
	gsub ("</h2>","",title) ;
	gsub ("<h3>",". ",title) ;
	gsub ("</h3>","",title) ;
    }
	
    i = index($0,"<div id=\\\"DescrAll\\\"") ; k = 10 ;
    if (i<1) {i = index($0,"<div class=\\\"Description\\\"") ; k=16 ;}
    if (i > 0)
    {
	z = substr ($0,i+k);
	i = index(z,">") ;
	z = substr (z,i+1);
	i = index(z,"</div>") ;
	j = index(z,"<a class") ; if (i > 1 && j > i)j=i ;
	descr = cleanSpace(substr(z,1,j-1)) ;
    }

    i=index($0," GEO: ") ;
    if (i > 0)
    {
	z = substr ($0,i+6);
	i=index(z,"<") ;
	z = substr (z,1,i-1);
	geo = z ;
    }
}
END {
    if (ab && dcr && index(ab,dcr) > 0)
	dcr  =  ""  ;
    printf("SRP %s\n-D Title\n-D Abstract\n-D Description\n-D Identifier\n-D Reference\n-D GEO\n", prj) ;
    
    if (length(title) > 1)
	printf("Title \"%s\"\n",title) ;
    if (length(descr) > 1)
	printf("Abstract \"%s\"\n",descr) ;
    if (length(id) > 1)
	printf("Identifier \"%s\"\n",id) ;
    if (length(geo) > 1)
	printf("GEO \"%s\"\n",geo) ;
    if (length(pm) > 1)
	print pm ;
    printf("\n") ;
    if (length(ab) > 1)
    { 
	printf("LongText \"%s\"\n",srp) ;
	print ab ;
        printf("***LongTextEnd***\n\n") ;
    }
}
