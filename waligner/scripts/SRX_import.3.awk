BEGIN { papers = "Reference" ; TTT = "" ; samples = ""; }
{ gsub (/\"/,"\\\"",$0); }
{ gsub (/http:\//,"http:\\/",$0) ;}

/>Contributor\(s\)</{ inAuth = 1 ;  next ; }  # Contributors are presented as a table
/>Organization name</ { inOrg = 1 ; next ; }
/>Country</ { inCountry = 1 ; next ; }
/>City</ { inCity = 1 ; next ; }

{ if (inAuth == 1)
    {
	nn = split ($0,aa,"</a>") ;
	for (ii = 1 ; ii <= nn ; ii++)
	{
	    z = aa[ii] ;
            while (substr (z,1,1) == " " || substr (z,1,1) == "," )
		z = substr(z, 2) ;
            while (substr (z,1,1) == "<")
	    {
		j = index (z, ">") ;
		if (j < 1) j = 1 ;
		z = substr(z, j+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0) 
		z = substr (z, 1, j-1) ;
            j = length (z) ;
	    if (j > 1) 
		authors = authors ", " z ;
	}
    }
    if (inOrg == 1)
    {
	nn = split ($0,aa,"</a>") ;
	for (ii = 1 ; ii <= nn ; ii++)
	{
	    z = aa[ii] ;
            while (substr (z,1,1) == " " || substr (z,1,1) == "," )
		z = substr(z, 2) ;
            while (substr (z,1,1) == "<")
	    {
		j = index (z, ">") ;
		if (j < 1) j = 1 ;
		z = substr(z, j+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0) 
		z = substr (z, 1, j-1) ;
            j = length (z) ;
	    if (j > 1) 
		org = org  ", " z ;
	}
    }
    if (inCountry == 1)
    {
	nn = split ($0,aa,"</a>") ;
	for (ii = 1 ; ii <= nn ; ii++)
	{
	    z = aa[ii] ;
            while (substr (z,1,1) == " " || substr (z,1,1) == "," )
		z = substr(z, 2) ;
            while (substr (z,1,1) == "<")
	    {
		j = index (z, ">") ;
		if (j < 1) j = 1 ;
		z = substr(z, j+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0) 
		z = substr (z, 1, j-1) ;
            j = length (z) ;
	    if (j > 1) 
		country = country  ", " z ;
	}
    }
    if (inCity == 1)
    {
	nn = split ($0,aa,"</a>") ;
	for (ii = 1 ; ii <= nn ; ii++)
	{
	    z = aa[ii] ;
            while (substr (z,1,1) == " " || substr (z,1,1) == "," )
		z = substr(z, 2) ;
            while (substr (z,1,1) == "<")
	    {
		j = index (z, ">") ;
		if (j < 1) j = 1 ;
		z = substr(z, j+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0) 
		z = substr (z, 1, j-1) ;
            j = length (z) ;
	    if (j > 1) 
		city = city ", " z ;
	}
    }
}

/>Title</{inExpType=1;next;}
/>Summary</{inExpType=1;next;}
/>Overall design</{inExpType=1;next;}
/>Experiment type</{inExpType=1;next;}
{ if(inExpType==1)
    {
	i = index($0,">") ;
	if (i>0)
	{
	    z = substr($0,i+1) ;
	    j = index(z,"<") ;
	    if (j > 1) { TTT = TTT "T \"" substr(z,1,j-1) "\"\n" ;}
	}
	inExpType = 0 ;
    }
}

/acc=GSM/{i=index($0,"acc=GSM");z=substr($0,i+4,50);j=index(z,"\"");samples=samples "Sample " substr(z,1,j-2) "\n";}

/<\/tr>/{ inAuth = 0 ; inOrg = 0 ;  inCountry = 0 ; inCity = 0 ; }

/pubmed_id/ {
	nn = split ($0,aa,"</a>") ;
	for (ii = 1 ; ii <= nn ; ii++)
	{
	    z = aa[ii] ;
            while (substr (z,1,1) == " " || substr (z,1,1) == "," )
		z = substr(z, 2) ;
            while (substr (z,1,1) == "<")
	    {
		j = index (z, ">") ;
		if (j < 1) j = 1 ;
		z = substr(z, j+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0) 
		z = substr (z, 1, j-1) ;
            j = length (z) ;
	    if (j > 1) 
	    papers = papers "\nReference pm" z ;
    }
}

END {
    printf ("GEO %s\n-D Author\n-D Reference\n", geo) ;
    z = "" ; dx = 3 ;
    if (length (authors) > 1)
    { z = z substr (authors,dx) ; dx = 1 ; }
    if (dx == 1) { dx = 2 ; z = z ";" ; } 
    if (length (org) > 1)
    { z = z substr (org,dx) ; dx = 1 ; }
    if (length (city) > 1)
    { z = z substr (city,dx) ; dx = 1 ; }
    if (length (city) > 1)
    { z = z substr (country,dx) ; dx = 1 ; }
    if (length (z) > 1) 
	printf("Author \"%s.\"\n", z) ;
   if (length (papers) > 10) 
       print papers ;
   if (length (samples) > 10) 
       printf("%s", samples) ;
   if (length (TTT) > 10) 
       printf("%s", TTT) ;
    printf ("\n") ;
}
